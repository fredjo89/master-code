#pragma once
#include "data_defs.h"

void print(Grid& grid);

void print(Matrix& mat);

void print(BInfo& B);

void print(H_data& H);

void print(Graph& grap);

void print(Node& node);

void print(Options& options);

void computeDiscrepancy(Grid& grid, double* ompSol);

void printDistribution(Graph& graph, Node& node, int RANK, int SIZE);

void writeBasisDistr(Graph& graph);

void writeBasis(Grid& g_grid);


void computeError(Grid& grid, double* solution, double* infNorm, double* twoNorm, int it);

void writeErrors(double* infNorm, double* twoNorm, int it);
