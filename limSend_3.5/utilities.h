#pragma once

#include "Basis.h"
#include "data_defs.h"

#include <iostream>
#include <map>
#include <set>
using namespace std;

void renormalize(Basis* basisArray, int nBasis, int RANK, int SIZE, int N);
int blockToThread(int block, int SIZE, int M);
int threadStart(int RANK, int SIZE, int M );
void makeH(Grid& grid, map<int, set<int> >& H, map<int,int>& typeTwoInfo,
                                                    int RANK, int*basisDistr);
void gatherBasis(Grid& grid, int* basisDistr, Basis* basisArray, int RANK,
                                        int SIZE, double* global_GSolution);
void computeDiscrepancy(Grid& grid, double* global_GSolution, double* ompSolution  );



void connStrenght(Grid& grid);

void checkNormalization(Basis* basisArray, int nBasis, int RANK, int SIZE, int N);
void printH(map<int, set<int> >& H);
void printArray(int* a, int l);
void printArray(double* a, int l);
void printSet(set<int> mySet);
void printMap(map<int,int> mymap);
void printMap(map<int, set<double*> >& myMap);
