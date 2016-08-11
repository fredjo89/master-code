#pragma once
#include <iostream>
#include "data_defs.h"
using namespace std;

void distributeProblem(Grid& globalGrid, Grid& localGrid, Matrix& globalMat, Matrix& localMat, Graph& graph, Node& node, int RANK, int SIZE);

void setupNodes(Graph& graph, Node& node, int RANK);

void makeCoarseGraph(Grid& grid, Matrix& mat, int SIZE, Graph& graph);

void makeFineGraph(H_data& H, Grid& grid, Graph& fineGraph);

void METIS_partition(int M, int* nodeWeights, int n_edges, int* edgeOffsets, int* edges, int* edgeWeights, int* basisDistr, int SIZE );

void makeH(Grid& grid, H_data& H);

void findSendCells(H_data& H, Graph& graph);



void reOrderGrid(Grid& grid, int* basisDistr, int SIZE);

void reOrderGrid(Grid& grid);

void destroyMatrix(Matrix& mat);

void printDistribution(Graph& graph, Node& node, int RANK, int SIZE);
