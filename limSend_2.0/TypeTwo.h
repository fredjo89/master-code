#pragma once
#include "data_defs.h"
#include "Basis.h"
#include <set>
#include <map>

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

  double** valAddress;  // addresses to values
  double** upAddress;   // addresses to update

  int* aIndex;          // Length: twoN
  int* twoN;            // Length: twoN
  double* sum;          // Length: twoN

  int* mapOne;          // Length: mapL;
  int* mapTwo;          // Length: mapL;

  double* recvBuff;    // Length: buffL;
  int* recvCounts;    // Needed for sending, length: SIZE
  int* displs;        // Needed for sending, length: SIZE  (displacement)

  int* basisDistr;      // Length: SIZE+1

  int* sendRanks;          // Array containing the ranks where sharing will be carried out.
  int Ncomm;                // Number of threads the thread will send/recv to.

  map<int, set<int> > H;  // H map




  TypeTwo();
  TypeTwo(Grid& grid, Basis* basisArray, int M, int S, int R, int* bD);
  ~TypeTwo();

  void localSum();        // Sum local type 2 cells
  void TTupdate();        // Update type two cells
  void setupSending();    // Create mapOne, mapTwo and recvBuff.
  void sendAndRecieve();  // Share recvBuff

  void print();
};
