#include "mpi.h"

#include<iostream>
using namespace std;

#define tag 100

main(int argc, char *argv[])  {
  int rank, size;



  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);



  int outbuf = rank;
  int dest[4] = {0,1,2,3};
  int recv[10]{};
  MPI_Request reqs[8];
  MPI_Status stats[8];

  for (int i=0; i<4; i++) {
     MPI_Isend(&outbuf, 1, MPI_INT, i, tag,
               MPI_COMM_WORLD, &reqs[i]);
     MPI_Irecv(&recv[i], 1, MPI_INT, i, tag,
               MPI_COMM_WORLD, &reqs[i+4]);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //MPI_Waitall(8, reqs, stats);
  if (rank==0) {for(int i = 0; i<10; i++) cout<<recv[i]<<endl;}
  if (rank==0) {for(int i = 0; i<8; i++) cout<<reqs[i]<<endl;}

	

  MPI_Finalize();
}
