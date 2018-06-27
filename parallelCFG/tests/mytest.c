#include "mpi.h"

int main(int argc, char* argv[])
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    int msg, total;
    
    if(rank == 0) {
        total = 0;
            #pragma pcfg_match (alph, 1) (beta, 1)
            MPI_Send(&rank, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    else if(rank == 1) { 
        // total = rank*100;
        #pragma pcfg_match (beta, 1) (alph, 1)
        MPI_Recv(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        // total += msg;
    }

    else if(rank == 2) { 
        // total = rank*100;
        #pragma pcfg_match (ga, 1) (alph, 1)
        MPI_Send(&msg, 1, MPI_INT, 3, 0, MPI_COMM_WORLD);
        // total += msg;
    }

    else if(rank == 3) { 
        // total = rank*100;
        #pragma pcfg_match (alph, 1) (ga, 1)
        MPI_Recv(&msg, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
        // total += msg;
    }   
    
    #pragma pcfg_merge
    MPI_Finalize();

    return 0;
}

