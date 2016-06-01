#pragma once
#include "data_defs.h"


void orderGS_Basis(Grid& grid, Matrix& mat, double* basis, Options opt);
void getOrders(Grid& grid, Matrix& mat, int* orders);

void makeOrderStructures(Grid& grid, Matrix& mat, Grid& orderGrid, Matrix& orderMat, int* orders);

void checkUnity(double* basis, int* support, int n_basis, int N, int M);
