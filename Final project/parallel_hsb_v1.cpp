#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <mpi.h>

int mpi_size;
int mpi_rank;
int nproc_x = 4, nproc_y = 4,nproc_z=4;
enum {ghost_cell_request, ghost_cell_answer};

bool verbose = false;

class spin_system {
    public:
    int flips = 100; // Number of flips the system will simulate.
    int n_spins = 64; // Number of spins in The system.
    int n_dims = 3; // Number of dimensions the spins are placed in.
    float N_spins_row; // Number of rows in the square/cube

    int x_offsets[6] = {-1,0,0,1,0,0}; // For calculating neighbor indices in energy calculation
    int y_offsets[6] = {0,-1,1,0,0,0};
    int z_offsets[6] = {0,0,0,0,-1,1};

    int x_dist = 1;
    int y_dist = 1;
    int z_dist = 1;

    int nearest_neighbours = 1; // Number of nearest neighbour interactions to be calculated
    int xlen, ylen, zlen;
    double J = 1; // Magnetization parameter, used in calculating energy.
    double H; // Total energy of the system;
    double B = 0; // Magnetic field in z direction
    double Temperature = 1; // Temperature of the system.
    std::string filename = "parallel_out.txt"; // Output file.
    int write = 1;
    std::vector<std::vector<double>> position; // Three by n_spins matrix, defining the spin's 3d position.
    std::vector<std::vector<double>> spin; // Three by n_spins matrix, defining the spin vector for each spin.
    std::vector<std::vector<int>>    neighbours; // 2*n_dims by n_spins matrix, defining the neighbour indices of each cell, so they need only be calculated once.
    std::vector<double>              energy; // Energy in each cell, derivaed at the end.
    spin_system(std::vector <std::string> argument){
    for (long unsigned int i = 1; i<argument.size() ; i += 2){
            std::string arg = argument.at(i);
            if(arg=="-h"){ // Write help
                std::cout << "Heisenberg_simulation\n --flips <number of flips performed>\n --nspins <number of spins simulated>\n --ndims <number of dimensions to simulate> \n"
                          << " --ofile <filename>\n --magnet <strength of external magnetic field in z direction>\n"
                          << " --temp <temperature> \n" << " --writeout <write to data file (1 for true, 0 for false)>\n";

                exit(0);
                break;
            } else if(arg=="--flips"){
                flips = std::stoi(argument[i+1]);
            } else if(arg=="--nspins"){
                n_spins = std::stoi(argument[i+1]);
            } else if(arg=="--ndims"){
                n_dims = std::stoi(argument[i+1]);
            } else if(arg=="--ofile"){
                filename = argument[i+1];
            } else if(arg=="--magnet"){
                B = std::stoi(argument[i+1]);
            } else if(arg=="--temp"){
                Temperature = std::stod(argument[i+1]);
            } else if(arg=="--writeout"){
                write = std::stoi(argument[i+1]);
            } else{
                std::cout << "---> error: the argument type is not recognized \n";
            }
        }
    N_spins_row = cbrt(double(n_spins)); //Equal size in all dimensions
    int n_spins_row = round(N_spins_row);
    xlen = n_spins_row;
    ylen = n_spins_row;
    zlen = n_spins_row;
    if(verbose) std::cout << "Nspins row " << n_spins_row << std::endl;
    }
};

class local_spins{
    public:
    int x_offsets[6] = {-1,0,0,1,0,0}; // For calculating neighbor indices in energy calculation
    int y_offsets[6] = {0,-1,1,0,0,0};
    int z_offsets[6] = {0,0,0,0,-1,1};

    int x_dist,y_dist,z_dist;

    int nearest_neighbours; // Number of nearest neighbour interactions to be calculated
    int xlen, ylen, zlen;
    int pad_xlen, pad_ylen, pad_zlen;
    int offset_x, offset_y, offset_z;
    int n_spins;
    double J; // Magnetization parameter, used in calculating energy.
    double H; // Total energy of the system;
    double B; // Magnetic field in z direction
    double Temperature; // Temperature of the system.
    std::string filename; // Output file.

    int no_in_padded_layer, no_in_layer;

    local_spins(spin_system &sys,
            int local_xlen, int local_ylen, int local_zlen,
            int offx, int offy, int offz){
        x_dist = sys.x_dist;
        y_dist = sys.y_dist;
        z_dist = sys.z_dist;
        
        nearest_neighbours = sys.nearest_neighbours;
        xlen = local_xlen;
        ylen = local_ylen;
        zlen = local_zlen;
        
        pad_xlen = xlen+2;
        pad_ylen = ylen+2;
        pad_zlen = zlen+2;
        
        n_spins = xlen*ylen*zlen;
        no_in_layer = xlen*ylen;
        no_in_padded_layer = (xlen+2)*(ylen+2);

        offset_x = offx;
        offset_y = offy;
        offset_z = offz;
        J = sys.J;
        H = sys.H;
        B = sys.B;
        Temperature = sys.Temperature;
        filename = sys.filename;
    };
    std::vector<std::vector<double>> position; // Three by n_spins matrix, defining the spin's 3d position.
    std::vector<std::vector<double>> spin; // Three by n_spins matrix, defining the spin vector for each spin.
    std::vector<std::vector<int>>    neighbours; // 2*n_dims by n_spins matrix, defining the neighbour indices of each cell, so they need only be calculated once.


    int index_to_padded_index(int index){
        int x = index%(no_in_layer)%xlen + 1;
        int y = (index%(no_in_layer))/xlen + 1;
        int z = index/(no_in_layer) + 1;

        return z*no_in_padded_layer + y*pad_xlen + x; 
    }

    int padded_index_to_index(int index){
        int x,y,z;
        padded_index_to_padded_coordinates(index,x,y,z);
        x -= 1;
        y -= 1;
        z -= 1;
        return z * xlen * ylen + y * xlen + x;
    }
    void padded_index_to_padded_coordinates(int index, int& x, int& y, int& z){
        x = index % pad_xlen; // Which row the spin is in
        y = (index/pad_ylen)%pad_ylen; // Which column the spin is in
        z = index / no_in_padded_layer;
    }

    void index_to_coordinates(int index, int& x, int& y, int& z){
        x = index % xlen;
        y = (index/ylen)%ylen;
        z = index / (ylen * xlen);
    }

    void padded_coordinates_to_padded_index(int &index,int x, int y, int z){
        index = x%pad_xlen + (y%pad_ylen) * pad_xlen + (z%pad_zlen) * no_in_padded_layer;
    }
};


// Function that generates rectangular positions for alle the spins in the system, 
void generate_positions_box(local_spins &sys){
        for (double k=0; k<sys.pad_zlen; k++)
        for (double j=0; j<sys.pad_ylen; j++)
        for (double i=0; i<sys.pad_xlen; i++)
            sys.position.push_back({double(i*sys.x_dist), double(j*sys.y_dist), double(k*sys.z_dist)});
};

// Function that generates random directions for all the spins in the system
void generate_spin_directions(local_spins &sys){
    
    for (int i = 0; i<sys.pad_zlen*sys.no_in_padded_layer; i++){
        srand(i); // Seed is here to make it perform nicely when comparing to parallel
        double spin_azimuthal = (double) rand()/RAND_MAX * M_PI;
        srand(i*rand()); // Seed is here to make it perform nicely when comparing to parallel
        double spin_polar = (double) rand()/RAND_MAX * 2. * M_PI;
        
        sys.spin.push_back({sin(spin_azimuthal)*cos(spin_polar), 
                            sin(spin_azimuthal)*sin(spin_polar),
                            cos(spin_azimuthal)});
    }   
};

void generate_neighbours(local_spins &sys){
    
    for (int spin = 0; spin< sys.pad_zlen*sys.no_in_padded_layer; spin++){
        // Find position in square / cube
        int spin_x, spin_y, spin_z;
        sys.padded_index_to_padded_coordinates(spin, spin_x, spin_y, spin_z);

        // Find indices of neighbours
        std::vector<int> spin_interactions;
        for(int i = 0; i < 6; i++){
            spin_interactions.push_back((spin_x + sys.x_offsets[i])%sys.pad_xlen + 
                                        (spin_y + sys.y_offsets[i])%sys.pad_ylen * sys.pad_xlen + 
                                        (spin_z + sys.z_offsets[i])%sys.pad_zlen * sys.no_in_padded_layer);
            //std::cout<< "RANK: " << mpi_rank << ". neighbour " << i << " Of padded local index " << spin << " is " << spin_interactions[i] << std::endl;
        }
        sys.neighbours.push_back(spin_interactions);
    }
}
// Function that calculates the energy of a single spin in 2d
double energy_calculation_nd(local_spins &sys, int spin, MPI_Comm& cart_comm){
    double energy = 0;
    double dot_product;
    for (int i=0; i<6; i++){
        // Calculate the energy with the nearest neighbour with no corners
        dot_product = sys.spin[spin][0]*sys.spin[sys.neighbours[spin][i]][0] 
                        + sys.spin[spin][1]* sys.spin[sys.neighbours[spin][i]][1]
                        + sys.spin[spin][2]* sys.spin[sys.neighbours[spin][i]][2];
        energy -= sys.J/2*dot_product;
    }
    energy += sys.B*sys.spin[spin][2];
    return energy;
};

// Calculate the total energy of the system
void Calculate_h(local_spins& sys, MPI_Comm cart_comm){
    sys.H = 0; // Set H to zero before recalculating it
    double mag_energy = 0;
    for (int i=0; i<sys.n_spins; i++){
        int pad_i = sys.index_to_padded_index(i);
        sys.H += energy_calculation_nd(sys, pad_i, cart_comm)*0.5; // Half the energy, because we calculate on all the spins
        mag_energy += sys.spin[pad_i][2]; 
    }
    sys.H += sys.B*mag_energy * 0.5; // Half of the magnetization energy is removed above
};

// Write the spin configurations in the output file.
void Writeoutput(spin_system& sys, std::ofstream& file, MPI_Comm cart_comm){  
    // Loop over all spins, and write out position and spin direction<
    file << "Position_x " << "Position_y " << "Position_z " << "Spin_x " <<  "Spin_y " <<  "Spin_z " << "Energy" << std::endl;
    for (int i = 0; i<sys.n_spins; i++){
        file << sys.position[i][0] << " " << sys.position[i][1] << " "  << sys.position[i][2] << " "
            << sys.spin[i][0] << " " << sys.spin[i][1] << " "  << sys.spin[i][2] << " " << sys.energy[i]
            << std::endl;
    }
};


void exchange_ghost_cells(local_spins &local_sys,
                        MPI_Aint &sdispls, MPI_Aint &rdispls, 
                        MPI_Datatype &sendtypes, MPI_Datatype &recvtypes,
                        MPI_Comm cart_comm){
    int counts[6] = {1,1,1,1,1,1};
    // Make Proof of concept work
    std::vector<double> sx;
    std::vector<double> sy;
    std::vector<double> sz;
    
    for (uint64_t i=0; i<local_sys.spin.size(); i++){
        sx.push_back(local_sys.spin[i][0]);
        sy.push_back(local_sys.spin[i][1]);
        sz.push_back(local_sys.spin[i][2]);
    }
    // Send ghostcells
    MPI_Neighbor_alltoallw (sx.data(), counts,  &sdispls, &sendtypes,
                            sx.data(), counts,  &rdispls, &recvtypes, cart_comm); 
    MPI_Neighbor_alltoallw (sy.data(), counts,  &sdispls, &sendtypes,
                            sy.data(), counts,  &rdispls, &recvtypes, cart_comm); 
    MPI_Neighbor_alltoallw (sz.data(), counts,  &sdispls, &sendtypes,
                            sz.data(), counts,  &rdispls, &recvtypes, cart_comm); 
    for(uint64_t i=0; i<local_sys.spin.size();i++){
        local_sys.spin[i] = {sx[i],sy[i],sz[i]};
    }

};

void Simulate(spin_system& sys, local_spins& localsys,MPI_Aint &sdispls, MPI_Aint &rdispls, 
                        MPI_Datatype &sendtypes, MPI_Datatype &recvtypes,
                        MPI_Comm cart_comm,
                        int neighbors[6],spin_system& global_sys){
    double old_energy, new_energy, spin_azimuthal, spin_polar, probability_of_change;
    std::vector<double> old_state(3);
    int not_flipped = 0;
    int flipped = 0;
    if(verbose) std::cout << old_state[0] << " " << old_state[1] << " " << old_state[2] << std::endl;
    
    
    if(verbose)std::cout << "Temp: " <<sys.Temperature
                            << std::endl;

    int local_iterations = sys.flips/mpi_size;
    int request_ghost_index;
    exchange_ghost_cells(localsys,sdispls, rdispls, 
                        sendtypes, recvtypes,
                        cart_comm);
    
    
    for (int iteration=0; iteration<local_iterations; iteration++){
        
        bool flip = false;
        
        // Choose a random spin site
        srand(iteration+mpi_rank);
        int rand_site = rand()%(localsys.n_spins);
        rand_site = localsys.index_to_padded_index(rand_site);
        
        // Calculate its old energy
        old_energy = energy_calculation_nd(localsys, rand_site, cart_comm);

        // Store its old state. 
        old_state[0] = localsys.spin[rand_site][0];
        old_state[1] = localsys.spin[rand_site][1];
        old_state[2] = localsys.spin[rand_site][2];
      
        // Generate new state
        spin_azimuthal = (double) rand()/RAND_MAX * M_PI;
        srand(mpi_rank*iteration + iteration);
        spin_polar = (double) rand()/RAND_MAX * 2. * M_PI;
        localsys.spin[rand_site] = {sin(spin_azimuthal)*cos(spin_polar), 
                            sin(spin_azimuthal)*sin(spin_polar),
                            cos(spin_azimuthal)};
        flipped++; 
        flip = true;
        // Calculate if it lowers energy
        new_energy = energy_calculation_nd(localsys, rand_site, cart_comm);
        if (new_energy > old_energy){
            // If not, see if it should be randomised in direction
            if (verbose) std::cout << "New Energy: " << new_energy << " Old energy: " << old_energy << std::endl;
            probability_of_change = exp(-(new_energy-old_energy)/localsys.Temperature); // FIgure out probability of change

                        srand((mpi_rank+1)*(iteration+1)*2);
            if (probability_of_change < (double) rand()/RAND_MAX){
                // If not, revert to old state
                localsys.spin[rand_site] = {
                    old_state[0], old_state[1], old_state[2]
                };
                not_flipped++;
                flipped--;
                flip = false;
                new_energy = old_energy;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        exchange_ghost_cells(localsys,sdispls, rdispls, 
                        sendtypes, recvtypes,
                        cart_comm);
    }
    if(verbose) std::cout <<"Finished my jobs /"<<mpi_rank<<"\n";
    MPI_Barrier(cart_comm);
    if(verbose){
        std::cout << "Not flipped no. is " << not_flipped << std::endl;
        std::cout << "Flipped no. is " << flipped << std::endl;
        std::cout << "Total energy: " << localsys.H << std::endl;
    }
    Calculate_h(localsys, cart_comm);
}

//=============================================================================================
//=========================   MAIN FUNCTION   =================================================
//=============================================================================================
int main(int argc, char* argv[]){
    if (verbose) std::cout << "Hello Heisenberg!" << std::endl;

    //MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (verbose) {
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    if(verbose) std::cout << "Heisenberg running on " << processor_name
              << ", rank " << mpi_rank << " out of " << mpi_size << std::endl;
    }

    //Initialise and load config
    spin_system global_sys({argv, argv+argc});

    //Setup MPI
    int dims[3] = {nproc_z, nproc_y, nproc_x};
    int periods[3] = {1,1,1};
    int coords[3];
    MPI_Dims_create(mpi_size, 3, dims);
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods,
            0, &cart_comm);
    MPI_Cart_coords(cart_comm, mpi_rank, 3, coords);

    int nleft, nright, nbottom, ntop, nfront, nback;
    MPI_Cart_shift(cart_comm, 2,1,&nleft,&nright);
    MPI_Cart_shift(cart_comm, 1,1,&nbottom,&ntop);
    MPI_Cart_shift(cart_comm, 0,1,&nfront,&nback);
    int neigbors[6] = {nleft, nright, nbottom, ntop, nfront, nback};

    const long int offset_x = global_sys.xlen * coords[2] / nproc_x -1;
    const long int offset_y = global_sys.ylen * coords[1] / nproc_y -1;
    const long int offset_z = global_sys.zlen * coords[0] / nproc_z -1;

    const long int end_x = global_sys.xlen * (coords[2]+1) / nproc_x +1; 
    const long int end_y = global_sys.ylen * (coords[1]+1) / nproc_y +1; 
    const long int end_z = global_sys.zlen * (coords[0]+1) / nproc_z +1;

    if(verbose) std::cout << mpi_rank << " " << end_x << " " << offset_x << " "<< end_y << " " << offset_y <<" "<< end_z<<" " << offset_z<<std::endl;
    local_spins local_sys(global_sys,
            end_x-offset_x-2, end_y-offset_y-2, end_z-offset_z-2,
            offset_x, offset_y, offset_z);
    //========================================================================================
    //========================= START OF GHOST CELL COMMUNICATION SETUP ======================
    //========================================================================================

    /* The following send blocks are defined as follows:
     * h_type sends to the nearest MPI block in the horizontal direction
     * v_type sends to the nearest MPI block in the vertical direction
     * d_type sends to the nearest MPI block in the depth direction
     * 
     * The blocks define the parts of data that will be sent in 
     * MPI_Neigbor_alltoallw.
     * 
     * HEAVILY INSPIRED BY https://github.com/essentialsofparallelcomputing/Chapter8/blob/master/GhostExchange/CartExchange3D_Neighbor/CartExchange.cc
     * LINE 110 AND FORWARD.
     */
    
    // Define subarray types for ghost cell exchanges
    const int array_sizes[] = {local_sys.pad_xlen,local_sys.pad_ylen, local_sys.pad_zlen};
    int subarray_sizes_h[] = { local_sys.zlen,local_sys.ylen,1};
    int subarray_h_start[] = {1,1,0};
    MPI_Datatype h_type;
    MPI_Type_create_subarray (3, array_sizes, subarray_sizes_h, subarray_h_start,
                            MPI_ORDER_C, MPI_DOUBLE, &h_type);
    MPI_Type_commit(&h_type);

    int subarray_sizes_v[] = {local_sys.zlen, 1, local_sys.xlen};
    int subarray_v_start[] = {1,0,1};
    MPI_Datatype v_type;
    MPI_Type_create_subarray (3, array_sizes, subarray_sizes_v, subarray_v_start,
                            MPI_ORDER_C, MPI_DOUBLE, &v_type);
    MPI_Type_commit(&v_type);

    int subarray_sizes_d[] = { 1, local_sys.zlen,local_sys.xlen};
    int subarray_d_start[] = {0,1,1};
    MPI_Datatype d_type;
    MPI_Type_create_subarray (3, array_sizes, subarray_sizes_d, subarray_d_start,
                            MPI_ORDER_C, MPI_DOUBLE, &d_type);
    MPI_Type_commit(&d_type);
    int element_size = sizeof(double);
    int nhalo = 1;
    int xyplane_mult = local_sys.pad_ylen*local_sys.pad_xlen*element_size; //8 because datatype is 3 doubles, 
    int xstride_mult = local_sys.pad_xlen*element_size;
    // Define displacements of send and receive in bottom top left right.
    MPI_Aint sdispls[6] = { nhalo          * xyplane_mult,
                            local_sys.zlen * xyplane_mult,
                            nhalo          * xstride_mult,
                            local_sys.ylen * xstride_mult,
                            nhalo          * element_size,
                            local_sys.xlen * element_size                          
                            };
    MPI_Aint rdispls[6] = {0,
                           (local_sys.zlen+1) * xyplane_mult,
                           0,
                           (local_sys.ylen+1) * xstride_mult,
                           0,
                           (local_sys.zlen+1) * element_size,
                           };
    
    MPI_Datatype sendtypes[6] = { d_type, d_type, v_type, v_type, h_type, h_type};
    MPI_Datatype recvtypes[6] = { d_type, d_type, v_type, v_type, h_type, h_type};
    
    //========================================================================================
    //========================= END OF GHOST CELL COMMUNICATION SETUP ========================
    //========================================================================================
    //Generate system 
    generate_positions_box(local_sys);
    generate_spin_directions(local_sys);
    generate_neighbours(local_sys);
    
    Calculate_h(local_sys, cart_comm);

    
    auto begin = std::chrono::steady_clock::now();
    // Run Simulation
    Simulate(global_sys,local_sys,*sdispls, *rdispls, 
                        *sendtypes, *recvtypes,
                        cart_comm,
                        neigbors, global_sys);
    auto end = std::chrono::steady_clock::now();

    
    MPI_Barrier(cart_comm);
    MPI_Reduce(&local_sys.H, &global_sys.H, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);
    if (mpi_rank == 0){ std::cout << "Final_energy: " << "Elapsed_time" << "Temperature " << "B_field " << "System_size " << "No_of_ranks " << "No_of_flips " << "Version " <<std::endl;
                        std::cout << global_sys.H << " " << (end-begin).count() / 1000000000.0 << " " << global_sys.Temperature << " " << global_sys.B <<
                              " " << global_sys.n_spins << " " << mpi_size << " " << global_sys.flips << " " << 1 << std::endl;
                          }

    // The rest is mainly for output file
    std::vector<double> px, py, pz;
    std::vector<double> sx, sy, sz;
    std::vector<double> energy;
    std::vector<double> globpx, globpy, globpz;
    std::vector<double> globsx, globsy, globsz;
    std::vector<double> globenergy;
    if (mpi_rank ==0) {
        globpx.reserve(global_sys.n_spins);
        globpy.reserve(global_sys.n_spins);
        globpz.reserve(global_sys.n_spins);
        globsx.reserve(global_sys.n_spins);
        globsy.reserve(global_sys.n_spins);
        globsz.reserve(global_sys.n_spins);
        globenergy.reserve(global_sys.n_spins);
        
    }

    int temp;
    for (int i = 0; i<local_sys.n_spins; i++){
                temp = local_sys.index_to_padded_index(i);
                px.push_back(local_sys.position[temp][0]);
                py.push_back(local_sys.position[temp][1]);
                pz.push_back(local_sys.position[temp][2]);
                sx.push_back(local_sys.spin[temp][0]);
                sy.push_back(local_sys.spin[temp][1]);
                sz.push_back(local_sys.spin[temp][2]);
                energy.push_back(energy_calculation_nd(local_sys,temp,cart_comm));
        }
    MPI_Barrier(cart_comm);
    
    MPI_Gather(px.data(), px.size(), MPI_DOUBLE, 
                globpx.data(), local_sys.n_spins, MPI_DOUBLE, 
                0, cart_comm);
    
    MPI_Gather(py.data(), py.size(), MPI_DOUBLE, 
                globpy.data(), local_sys.n_spins, MPI_DOUBLE, 
                0, cart_comm);
    MPI_Gather(pz.data(), pz.size(), MPI_DOUBLE, 
                globpz.data(), local_sys.n_spins, MPI_DOUBLE, 
                0, cart_comm);
    MPI_Gather(sx.data(), sx.size(), MPI_DOUBLE, 
                globsx.data(), local_sys.n_spins, MPI_DOUBLE, 
                0, cart_comm);
    MPI_Gather(sy.data(), sy.size(), MPI_DOUBLE, 
                globsy.data(), local_sys.n_spins, MPI_DOUBLE, 
                0, cart_comm);
    MPI_Gather(sz.data(), sz.size(), MPI_DOUBLE, 
                globsz.data(), local_sys.n_spins, MPI_DOUBLE, 
                0, cart_comm); 
    MPI_Gather(energy.data(), energy.size(), MPI_DOUBLE, 
                globenergy.data(), local_sys.n_spins, MPI_DOUBLE, 
                0, cart_comm); 
    
    if (mpi_rank == 0){
        for (int i=0; i<global_sys.n_spins; i++){
            global_sys.spin.push_back({globsx[i], globsy[i], globsz[i]});
            global_sys.position.push_back({globpx[i], globpy[i], globpz[i]});
            global_sys.energy.push_back(globenergy[i]);
        }
        std::ofstream file(global_sys.filename); // open file
        Writeoutput(global_sys, file, cart_comm);
    }
    MPI_Type_free(&h_type);
    MPI_Type_free(&v_type);
    MPI_Type_free(&d_type);
    MPI_Barrier(cart_comm);

    MPI_Finalize();
    return 0;
}
