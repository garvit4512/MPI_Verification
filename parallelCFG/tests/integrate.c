// #include <stdio.h>
#include "mpi.h"

int main () {
  
  int argc; 
  char **argv;
  
  int rank, size;
  int sbuff, rbuff;
  MPI_Status status;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
    if(rank == 0) {
#pragma pcfg_match (alpha, 1) (beta, 1)
        MPI_Send(&sbuff, 1, MPI_INT, 2, 1, MPI_COMM_WORLD);
#pragma pcfg_match (alpha, 0) (star, 1)
        MPI_Recv(&rbuff, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
#pragma pcfg_match (alpha, 2) (gamma, 1)
        MPI_Send(&sbuff, 1, MPI_INT, 2, 1, MPI_COMM_WORLD);
#pragma pcfg_match (alpha, 0) (star, 2)
        MPI_Recv(&rbuff, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
    }
    else if (rank == 1) {
#pragma pcfg_match (beta, 2) (alpha, 3)
        MPI_Send(&sbuff, 1, MPI_INT, 2, 1, MPI_COMM_WORLD);
#pragma pcfg_match (beta, 1) (alpha, 1)
        MPI_Recv(&sbuff, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
     }
    else if (rank == 2) {
#pragma pcfg_match (gamma, 2) (alpha, 3)
       MPI_Send(&rbuff, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
#pragma pcfg_match (gamma, 1) (alpha, 2)
       MPI_Recv(&rbuff, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
    }
    MPI_Finalize();
    return 0;
}

