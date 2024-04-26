#include <iostream>
#include <random>
#include <chrono>
#include <thread>
#include <array>

// To run an MPI program we always need to include the MPI headers
#include <mpi.h>

const int NTASKS=5000;  // number of tasks
const int RANDOM_SEED=1234;

void master (int nworker) {
    std::array<int, NTASKS> task, result;

    // set up a random number generator
    std::random_device rd;
    //std::default_random_engine engine(rd());
    std::default_random_engine engine;
    engine.seed(RANDOM_SEED);
    // make a distribution of random integers in the interval [0:30]
    std::uniform_int_distribution<int> distribution(0, 30);

    for (int& t : task) {
        t = distribution(engine);   // set up some "tasks"
    }

    // initialize: send a task to each worker
    for (int i = 0; i < nworker; i++) {
            MPI_Send(&task[i], 1, MPI_INT, i+1, 0, MPI_COMM_WORLD); }

    int sent_tasks = nworker;  // number of tasks sent to workers
    int recieved = 0;  // number of recieved results

    while (sent_tasks < NTASKS) {  // loop until all tasks have been sent to a worker
        MPI_Status status;
        MPI_Recv(&result[recieved], 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);  // recieve task from any worker

        int source = status.MPI_SOURCE;  // get worker source

        MPI_Send(&task[sent_tasks], 1, MPI_INT, source, 0, MPI_COMM_WORLD);  // send new task to the worker
        
        sent_tasks++;  // one more task has been sent
        recieved++;  // one more task has been recieved
        
        }
    int stop_signal = -1;
    for (int i = 0; i < nworker; i++) {
        MPI_Status status;
        MPI_Recv(&result[recieved], 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);  // recieve last results from worker
        int source = status.MPI_SOURCE;
        MPI_Send(&stop_signal, 1, MPI_INT, source, 0, MPI_COMM_WORLD);  // terminate current worker
        recieved++;
    }

    // Print out a status on how many tasks were completed by each worker
    for (int worker=1; worker<=nworker; worker++) {
        int tasksdone = 0; int workdone = 0;
        for (int itask=0; itask<NTASKS; itask++)
        if (result[itask]==worker) {
            tasksdone++;
            workdone += task[itask];
        }
        std::cout << "Master: Worker " << worker << " solved " << tasksdone << 
                    " tasks" << "corresponding to the work " << workdone << "\n";    
    }
}

// call this function to complete the task. It sleeps for task milliseconds
void task_function(int task) {
    std::this_thread::sleep_for(std::chrono::milliseconds(task));
}

void worker (int rank) {
  
    while (true) {
    int received_task;
    MPI_Recv(&received_task, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // recieve task from master
    if (received_task == -1) break;  // stop if master send signal to stop

    task_function(received_task);  // perform task

    MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);  // sent result to master
    }
    
}

int main(int argc, char *argv[]) {
    int nrank, rank;

    MPI_Init(&argc, &argv);                // set up MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nrank); // get the total number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rank of this process

    if (rank == 0)       // rank 0 is the master
        master(nrank-1); // there is nrank-1 worker processes
    else                 // ranks in [1:nrank] are workers
        worker(rank);

    MPI_Finalize();      // shutdown MPI
}
