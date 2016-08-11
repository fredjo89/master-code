#pragma once
#include "data_defs.h"

void setupSupportMatrix(Grid& grid, Matrix& mat);

void setupSupportMatrix2(Grid& grid, Matrix& mat);

void makeMapping(Grid& grid, Matrix& mat);

void makeTypeTwo_local(Grid& grid, Node& node, TypeTwo& TypeTwo);

void makeTypeTwo_sending(Grid& grid, TypeTwo& tt, int RANK);

void makeTypeTwo_serial(Grid& grid, TypeTwo& tt);

void jacobi(Grid& grid, Matrix& mat, double& omega);

void localSum(Grid& grid, TypeTwo& tt, int RANK);

void sendAndRecieve(Grid& grid, TypeTwo& tt, int RANK);

void TTnormalize(Grid& grid, TypeTwo& tt, int RANK);
