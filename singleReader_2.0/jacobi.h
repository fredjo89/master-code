#pragma once
#include "data_defs.h"

void jacobi(Grid& grid, Matrix& mat, double& omega);

void makeTypeTwo(Grid& grid, Node& node, int& SIZE, TypeTwo& TypeTwo);

void makeTypeTwo_local(Grid& grid, Node& node, TypeTwo& TypeTwo);

void makeTypeTwo_sending(Grid& grid, TypeTwo& typeTwo);

void makeTypeTwo_serial(Grid& grid, TypeTwo& typeTwo);

void localSum(Grid& grid, TypeTwo& typeTwo, int RANK);

void sendAndRecieve(Grid& grid, TypeTwo& typeTwo, int RANK);

void TTnormalize(Grid& grid, TypeTwo& typeTwo, int RANK);
