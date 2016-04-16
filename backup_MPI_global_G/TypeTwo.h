#pragma once
#include "data_defs.h"


class TypeTwo{
public:
  double** valAddress;  // addresses to values
  double** upAddress;   // addresses to update

  int* aIndex;          // Length: twoN
  int* twoN;            // Length: twoN
  double* sum;          // Length: twoN

  int addL;             // length of valAddress and upAddress
  int sumL;            // length of twoN, sum and aIndex;

  int SIZE;            // Number of threads
  int RANK;            // Rank of thread

  int* basisDistr;

  int mapL;           // Length of mapOne and mapTwo;
  int* mapOne;
  int* mapTwo;

  int buffL;          // Length of recvBuff
  double* recvBuff;
  int* recvCounts;    // Needed for sending, length: SIZE
  int* displs;        // Needed for sending, length: SIZE  (displacement)



  TypeTwo(){valAddress=NULL; upAddress=NULL; aIndex=NULL; twoN = NULL; sum = NULL;
    mapOne = NULL; mapTwo = NULL; recvBuff = NULL;  addL = 0; sumL = 0; buffL = 0; }
  TypeTwo(Basis* basisArray, int M, int S, int R, int* bD);
  ~TypeTwo(){delete[] valAddress; delete[] upAddress,  delete[] aIndex;
            delete[] twoN; delete[] sum; delete[] mapOne; delete[] mapTwo; delete[] recvBuff;
            delete[] recvCounts; delete[] displs;}

  void print();

  void localSum();        // Sum local type 2 cells
  void globalSum();       // Sum global type 2 cells
  void TTupdate();        // Update type two cells

  void setupSending();    // Create mapOne, mapTwo and recvBuff.
  void sendAndRecieve();  // Share recvBuff


};
