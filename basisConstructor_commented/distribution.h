// Functions used to distribute the problem.

#pragma once
#include "data_defs.h"
// makeCoarseGraph(): Creates the distribution of basis functions among processes, and finds the
// data that much be exchanged.
// Only used by process 0.
// Uses the following functions:
// makeH(), makeFineGraph(), METIS_partition(), findSendCells().
void makeCoarseGraph(Grid& grid, Matrix& mat, int& SIZE, Graph& graph);

// Creates the data structure H_data.
void makeH(Grid& grid, H_data& H);

// Creates the graph describing basis function computations and dependence.
// fineGraph is the graph to be partitioned to distribute the computations.
// makeFineGraph2() is alternative to makeFineGraph.
void makeFineGraph(Grid& grid, H_data& H, Graph& fineGraph);
void makeFineGraph2(Grid& grid, H_data& H, Graph& fineGraph);

// Partitions fineGraph by the use of METIS.
void METIS_partition(Graph& graph, Graph& fineGraph);

// findSendCells(): After partition has been completed, the function finds what cell-sums
// must be sent for each process and to where.
void findSendCells(H_data& H, Graph& graph);

// setupNodes(): Creates the Node data structure on each process, with key information
// of data to be distributed.
void setupNodes(Graph& graph, Node& node, int& RANK);

// distributeProblem(): Distributes the problem data form process 0 to all others.
// Uses sortNDremvDup().
void distributeProblem(Grid& g_grid, Grid& l_grid, Matrix& g_mat, Matrix& l_mat,
Graph& graph, Node& node);

// sortNDremvDup(): Sort int-array and removes duplicates.
void sortNDremvDup(int* input, int& n_input, int* output, int& n_output);

// Gather completed basis function from all processes to process 0.
void gatherBasis(Grid& l_grid, Grid& g_grid, Graph& graph, int RANK, int SIZE);
