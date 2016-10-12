// Functions used in the development of the program. Not used in the finished program.

#pragma once
#include "data_defs.h"

// Print to screen functions.
void print(Grid& grid);
void print(Matrix& mat);
void print(BInfo& B);
void print(H_data& H);
void print(Graph& grap);
void print(Node& node);
void print(Options& options);

// Compute Discrepancy between two prolongation operators.
void computeDiscrepancy(Grid& grid, double* ompSol);

// Print / write basis function distribtuion among processes.
void printDistribution(Graph& graph, Node& node, int RANK, int SIZE);
void writeBasisDistr(Graph& graph);
