#pragma once
#include "data_defs.h"
#include <iostream>
#include <map>
using namespace std;

class TypeTwo{
public:
  double** valAddress;  // addresses to values
  double** upAddress;   // addresses to update

  int* aIndex;          // Length: twoN
  int* twoN;            // Length: twoN
  double* sum;          // Length: twoN
  
  int addL;             // length of address;
  int sumL;            // length of twoN, sum and aIndex;

  TypeTwo(){valAddress=NULL; upAddress=NULL; aIndex=NULL; twoN = NULL;
    sum = NULL; addL = 0; sumL = 0; }
  TypeTwo(Basis* basisArray, int M);
  ~TypeTwo(){delete[] valAddress; delete[] upAddress,  delete[] aIndex;
                                            delete[] twoN; delete[] sum; }

  void print();
  void sumAndModify();
};
