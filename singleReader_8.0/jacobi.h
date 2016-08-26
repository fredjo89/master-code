#pragma once
#include "data_defs.h"

void setupSupportMatrix(Grid& grid, Matrix& mat);

void setupSupportMatrix2(Grid& grid, Matrix& mat);

void makeMapping(Grid& grid, Matrix& mat);

void makeBInfo_local(Grid& grid, Node& node, BInfo& B);

void makeBInfo_sending(Grid& grid, BInfo& B, int RANK);

void makeBInfo_serial(Grid& grid, BInfo& B);

void jacobi(Grid& grid, Matrix& mat, Options& opt);

void localSum(Grid& grid, BInfo& B);

void sendAndRecieve(BInfo& B);

void BNormalize(Grid& grid, BInfo& B);
