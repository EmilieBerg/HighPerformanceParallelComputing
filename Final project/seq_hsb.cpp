#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <cstdlib>

bool verbose = false;
class spin_system {
    public:
    int flips = 100; // Number of flips the system will simulate.
    int n_spins = 216; // Number of spins in The system.
    int n_dims = 3; // Number of dimensions the spins are placed in.
    int n_spins_row; // Number of rows in the square/cube

    int row_offsets[6] = {-1,0,0,1,0,0}; // For calculating neighbor indices in energy calculation
    int col_offsets[6] = {0,-1,1,0,0,0};
    int dep_offsets[6] = {0,0,0,0,-1,1};

    int x_dist = 10;
    int y_dist = 10;
    int z_dist = 10;
    int nearest_neighbours = 1; // Number of nearest neighbour interactions to be calculated
    double J = 1; // Magnetization parameter, used in calculating energy.
    double H; // Total energy of the system;
    double B = 0; // Magnetic field in z direction
    double Temperature = 1; // Temperature of the system.
    std::string filename = "seq_out.txt"; // Output file.
    int write = 1;
    std::vector<std::vector<double>> position; // Three by n_spins matrix, defining the spin's 3d position.
    std::vector<std::vector<double>> spin; // Three by n_spins matrix, defining the spin vector for each spin.
    std::vector<std::vector<int>>    neighbours; // 2*n_dims by n_spins matrix, defining the neighbour indices of each cell, so they need only be calculated once.
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
    n_spins_row = pow(float(n_spins),1./float(n_dims)); //Equal size in all dimensions
    }
};


// Function that generates rectangular positions for alle the spins in the system, 
void generate_positions_box(spin_system &sys){
    switch (sys.n_dims)
    {
    case 1:
        for (double i=0; i<sys.n_spins; i++)
            sys.position.push_back({i*sys.x_dist, 0, 0});
        break;
    case 2:
        for (double i=0; i<sqrt(sys.n_spins)+2; i++)
        for (double j=0; j<sqrt(sys.n_spins)+2; j++)
            sys.position.push_back({i*sys.x_dist, j*sys.y_dist, 0});
        break;
    case 3:
        for (double i=0; i<cbrt(sys.n_spins); i++)
        for (double j=0; j<cbrt(sys.n_spins); j++)
        for (double k=0; k<cbrt(sys.n_spins); k++)
            sys.position.push_back({i*sys.x_dist, j*sys.y_dist, k*sys.z_dist});
        break;
    default:
        break;
    }

};

// Function that generates random directions for all the spins in the system
void generate_spin_directions(spin_system &sys){
    
    for (int i = 0; i<sys.n_spins; i++){
        srand(i); // Seed is here to make it perform nicely when comparing to parallel
        double spin_azimuthal = (double) rand()/RAND_MAX * M_PI;
        srand(i*rand()); // Seed is here to make it perform nicely when comparing to parallel
        double spin_polar = (double) rand()/RAND_MAX * 2. * M_PI;
        
        sys.spin.push_back({sin(spin_azimuthal)*cos(spin_polar), 
                            sin(spin_azimuthal)*sin(spin_polar),
                            cos(spin_azimuthal)});
    }   
};

void generate_neighbours(spin_system &sys){
    for (int spin = 0; spin<sys.n_spins; spin++){
        // Find position in square / cube
        int spin_row = spin/sys.n_spins_row; // Which row the spin is in
        if(sys.n_dims>2) spin_row = spin_row % sys.n_spins_row;
        int spin_col = spin%sys.n_spins_row; // Which column the spin is in
        int spin_dep = spin / (sys.n_spins_row*sys.n_spins_row);
        
        // Find indices of neighbours
        std::vector<int> spin_interactions;
        for(int i = 0; i < 2* sys.n_dims; i++){
            spin_interactions.push_back((spin_row + sys.row_offsets[i])%sys.n_spins_row*sys.n_spins_row + 
                                        (spin_col + sys.col_offsets[i])%sys.n_spins_row + 
                                        (spin_dep + sys.dep_offsets[i])%sys.n_spins_row * (sys.n_spins_row * sys.n_spins_row));
            if (spin_interactions[i]<0) spin_interactions[i] *= -1; // Should'nt be necessary, but modulo is not good with - 
        }
        sys.neighbours.push_back(spin_interactions);
    }
}
// Function that calculates the energy of a single spin in 2d
double energy_calculation_nd(spin_system &sys, int spin){
    
    double energy = 0;
    double dot_product;
    for (int i=0; i<2*sys.n_dims; i++){
        // Calculate the energy with the nearest neighbour with no corners
        dot_product = sys.spin[spin][0]*sys.spin[sys.neighbours[spin][i]][0] 
                        + sys.spin[spin][1]* sys.spin[sys.neighbours[spin][i]][1]
                        + sys.spin[spin][2]* sys.spin[sys.neighbours[spin][i]][2];
        energy += - sys.J/2*dot_product;
    }
    return energy;
};

// Calculate the total energy of the system
void Calculate_h(spin_system& sys){
    sys.H = 0; // Set H to zero before recalculating it
    for (int i=0; i<sys.n_spins; i++){
        sys.H += energy_calculation_nd(sys, i)*0.5; // Half the energy, because we calculate on all the spins
    }
};

// Write the spin configurations in the output file.
void Writeoutput(spin_system& sys, std::ofstream& file){  
    // Loop over all spins, and write out position and spin direction
    file << "Position_x " << "Position_y " << "Position_z " << "Spin_x " <<  "Spin_y " <<  "Spin_z " <<  "Spin_energy " << "Temperature " << "n_spins" << std::endl;
    for (int i = 0; i<sys.n_spins; i++){
        if (i == 0) {
            file << sys.position[i][0] << " " << sys.position[i][1] << " "  << sys.position[i][2] << " "
                << sys.spin[i][0] << " " << sys.spin[i][1] << " "  << sys.spin[i][2] << " "
                << energy_calculation_nd(sys,i) << " " << sys.Temperature << " " << sys.n_spins
                << std::endl;
        } else {
            file << sys.position[i][0] << " " << sys.position[i][1] << " "  << sys.position[i][2] << " "
            << sys.spin[i][0] << " " << sys.spin[i][1] << " "  << sys.spin[i][2] << " "
            << energy_calculation_nd(sys,i)
            << std::endl;
        }
    }
};

void Simulate(spin_system& sys){
    double old_energy, new_energy, spin_azimuthal, spin_polar, probability_of_change;
    std::vector<double> old_state (3);
    int not_flipped = 0;
    int flipped = 0;
    if(verbose) std::cout << old_state[0] << " " << old_state[1] << " " << old_state[2] << std::endl;
    
    auto begin = std::chrono::steady_clock::now();
    //std::cout << "Temp: " <<sys.Temperature//<< " Boltzmann: " << Boltzmann 
    //                        << std::endl;
    for (int iteration=0; iteration<sys.flips; iteration++){
        // Choose a random spin site
        srand(iteration);
        int rand_site = rand()%(sys.n_spins);
        
        // Calculate it's old energy
        old_energy = energy_calculation_nd(sys, rand_site);

        // Store it's old state. 
        old_state[0] = sys.spin[rand_site][0];
        old_state[1] = sys.spin[rand_site][1];
        old_state[2] = sys.spin[rand_site][2];
        
        // Generate new state
        spin_azimuthal = (double) rand()/RAND_MAX * M_PI;
        srand(iteration + 1);
        spin_polar = (double) rand()/RAND_MAX * 2. * M_PI;
        sys.spin[rand_site] = {sin(spin_azimuthal)*cos(spin_polar), 
                            sin(spin_azimuthal)*sin(spin_polar),
                            cos(spin_azimuthal)};
        flipped++; 
        // Calculate if it lowers energy
        new_energy = energy_calculation_nd(sys, rand_site);
        if (new_energy > old_energy){
            // Commenting out this bit for now, since p is always zero anyway. Then it can go a bit faster
            
            // If not, see if it should be randomised in direction
            if (verbose) std::cout << "New Energy: " << new_energy << " Old energy: " << old_energy << std::endl;
            probability_of_change = exp(-(new_energy-old_energy)/sys.Temperature); // FIgure out probability of change
            //std::cout << "Change prob: " << probability_of_change << " Exp factor :" << -(new_energy-old_energy)/(Boltzmann*sys.Temperature) << std::endl;
            srand(iteration*2);
            if (probability_of_change < (double) rand()/RAND_MAX){
                // If not, revert to old state
                sys.spin[rand_site] = {
                    old_state[0], old_state[1], old_state[2]
                };
                not_flipped += 1;
                flipped--;
                new_energy = old_energy;
            }

            
            // Remove next few lines when commenting the above back in
            /*sys.spin[rand_site] = {
                    old_state[0], old_state[1], old_state[2]
                };
            not_flipped += 1;
            new_energy = old_energy;*/
        }
        
        // Change H to represent the total energy of the system. Gave wrong results. Unclear why. Currently commented out. H is just calculated at the end, as it is not used anywhere in the loop anyway.
        //sys.H = sys.H - old_energy + new_energy; 

    }
    Calculate_h(sys);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Final_energy: " << "Elapsed_time" << "Temperature " << "B_field " << "System_size " << "No_of_ranks " << "No_of_flips " << "Version " <<std::endl;
    std::cout << sys.H << " " << (end-begin).count() / 1000000000.0 << " " << sys.Temperature << " " << sys.B <<
                              " " << sys.n_spins << " " << 1 << " " << sys.flips << " " << 0 << std::endl;
}
//=============================================================================================
//=========================   MAIN FUNCTION   =================================================
//=============================================================================================
int main(int argc, char* argv[]){
    //std::cout << "Hello Heisenberg!" << std::endl;

    //Load config
    spin_system sys({argv, argv+argc});

    //Generate system
    generate_positions_box(sys);
    generate_spin_directions(sys);
    generate_neighbours(sys);

    //Magic
    Calculate_h(sys);
    Simulate(sys);

    if (sys.write) {
        std::ofstream file(sys.filename); // open file
        Writeoutput(sys, file);
    }
}
