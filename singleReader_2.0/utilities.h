#pragma once
#include "data_defs.h"

void setupSupportMatrix(Grid& grid, Matrix& mat);

void makeMapping(Grid& grid, Matrix& mat);

void gatherBasis(Grid& localGrid, Grid& globalGrid, Graph& graph, int RANK, int SIZE);

void computeDiscrepancy(Grid& grid, double* ompSol);

void destroyMatrix(Matrix& mat);

void printDistribution(Graph& graph, Node& node, int RANK, int SIZE);

void reOrderGrid(Grid& grid, int* basisDistr, int SIZE);

void reOrderGrid(Grid& grid);
