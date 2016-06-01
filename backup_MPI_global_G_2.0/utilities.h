#pragma once

#include <map>
#include <set>
#include <vector>
#include "data_defs.h"

int blockToThread(int block, int SIZE, int M);

int threadStart(int RANK, int SIZE, int M );

void printArray(int* a, int l);

void printArray(double* a, int l);

void makeH(Grid& grid, map<int, set<int> >& H, int RANK, int*basisDistr);
