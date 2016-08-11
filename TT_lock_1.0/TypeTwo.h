#pragma once
#include "data_defs.h"
#include "Basis.h"
#include <set>
#include <map>
#include <mpi.h>
class TypeTwo{
public:
  double omega;         // Relaxation factor
  int SIZE;            // Number of threads
  int RANK;            // Rank of thread

  int addL;             // length of valAddress and upAddress
  int sumL;            // length of twoN, sum and aIndex;
  int mapL;           // Length of mapOne and mapTwo;
  int buffL;          // Length of recvBuff
  int sendL;          // Number of BValues that must be sent
  int Ncomm;          // Length of sendRanks.

  double** valAddress = NULL;  // addresses to values
  double** upAddress = NULL;   // addresses to update

  int* aIndex = NULL;          // Length: twoN
  int* twoN = NULL;            // Length: twoN
  double* sum = NULL;          // Length: twoN

  int* mapOne = NULL;          // Length: mapL;
  int* mapTwo = NULL;          // Length: mapL;

  double* recvBuff = NULL;    // Length: buffL;
  int* recvCounts = NULL;    // Needed for sending, length: SIZE
  int* displs = NULL;        // Needed for sending, length: SIZE  (displacement)

  int* basisDistr;      // Length: SIZE+1

  int* sendRanks;          // Array containing the ranks where sharing will be carried out.

  map<int, set<int> > H;  // H map

  MPI_Status* statSend;
  MPI_Status* statRecv;
  MPI_Request* send_request;
  MPI_Request* recv_request;


  TypeTwo(Grid& grid, Basis* basisArray, int M, int S, int R, int* bD);
  ~TypeTwo();

  void localSum();        // Sum local type 2 cells
  void TTupdate();        // Update type two cells
  void setupSending();    // Create mapOne, mapTwo and recvBuff.
  void sendAndRecieve();  // Share recvBuff

  void print();
};
