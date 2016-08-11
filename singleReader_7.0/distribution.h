#pragma once
#include "data_defs.h"

void makeCoarseGraph(Grid& grid, Matrix& mat, int& SIZE, Graph& graph);

void makeH(Grid& grid, H_data& H);

void makeFineGraph(Grid& grid, H_data& H, Graph& fineGraph);

void makeFineGraph2(Grid& grid, H_data& H, Graph& fineGraph);

void METIS_partition(Graph& graph, Graph& fineGraph);

void findSendCells(H_data& H, Graph& graph);

void setupNodes(Graph& graph, Node& node, int& RANK);

void distributeProblem(Grid& g_grid, Grid& l_grid, Matrix& g_mat, Matrix& l_mat,
Graph& graph, Node& node);

void sortNDremvDup(int* input, int& n_input, int* output, int& n_output);

void gatherBasis(Grid& l_grid, Grid& g_grid, Graph& graph, int RANK, int SIZE);
