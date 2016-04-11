#pragma once
#include "data_defs.h"
#include <iostream>
using namespace std;

class TypeTwo{
public:
  double** valAddress;  // addresses to values
  double** upAddress;   // addresses to update

  int* aIndex;          // Length: twoN
  int* twoN;            // Length: twoN
  double* sum;          // Length: twoN

  int addL;             // length of valAddress and upAddress
  int sumL;            // length of twoN, sum and aIndex;

  double** mapOne;
  double** mapTwo;
  double* recvBuff;

  TypeTwo(){valAddress=NULL; upAddress=NULL; aIndex=NULL; twoN = NULL;
   sum = NULL; mapOne = NULL; mapTwo = NULL; recvBuff = NULL;  addL = 0; sumL = 0; }
  TypeTwo(Basis* basisArray, int M);
  ~TypeTwo(){delete[] valAddress; delete[] upAddress,  delete[] aIndex;
                                            delete[] twoN; delete[] sum; }

  void print();

  void localSum();    // Sum local type 2 cells
  void globalSum();   // Sum global type 2 cells
  void TTupdate();      // Update type two cells

  void setupSending();
  void sendAndRecieve();


};
