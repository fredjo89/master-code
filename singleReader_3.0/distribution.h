#pragma once
#include "data_defs.h"

void makeCoarseGraph(Grid& grid, Matrix& mat, int& SIZE, Graph& graph);

void makeH(Grid& grid, H_data& H);

void makeFineGraph(Grid& grid, H_data& H, Graph& fineGraph);

void makeFineGraph2(Grid& grid, H_data& H, Graph& fineGraph);

void makeFineGraph3(Grid& grid, H_data& H, Graph& fineGraph);

void METIS_partition(Graph& graph, Graph& fineGraph);

void findSendCells(H_data& H, Graph& graph);

void setupNodes(Graph& graph, Node& node, int& RANK);

void distributeProblem(Grid& globalGrid, Grid& localGrid, Matrix& globalMat, Matrix& localMat,
Graph& graph, Node& node);
