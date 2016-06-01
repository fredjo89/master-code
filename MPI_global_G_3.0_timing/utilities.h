#pragma once

#include "Basis.h"

#include <map>
#include <set>
#include <vector>
#include "data_defs.h"

int blockToThread(int block, int SIZE, int M);

int threadStart(int RANK, int SIZE, int M );

void printArray(int* a, int l);

void printArray(double* a, int l);

void makeH(Grid& grid, map<int, set<int> >& H, int RANK, int*basisDistr);

void printH(map<int, set<int> >& H);


void gatherBasis(Grid& grid, int* basisDistr, Basis* basisArray, int RANK, int SIZE,
double* global_GSolution);

void computeDiscrepancy(Grid& grid, double* global_GSolution,
  double* ompSolution  );

void printSet(set<int> mySet); 
