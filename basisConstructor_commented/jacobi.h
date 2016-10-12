// Functions performing iteration procedure and setup of BInfo.

#pragma once
#include "data_defs.h"

// Creates matrix structure for each basis function.
void setupSupportMatrix(Grid& grid, Matrix& mat);

// Changes sup_index in Matrix to map to locations by index.
void makeMapping(Grid& grid, Matrix& mat);

// Creates local part of BInfo.
void makeBInfo_local(Grid& grid, Node& node, BInfo& B);

// Creates message-passing-info part of BInfo. Performs message-passing.
void makeBInfo_sending(Grid& grid, BInfo& B, int RANK);

// Makes BInfo for a single process run.
void makeBInfo_serial(Grid& grid, BInfo& B);

// Jacobi sweep.
void jacobi(Grid& grid, Matrix& mat, Options& opt);

// Computes local sum of boundary cells.
void localSum(Grid& grid, BInfo& B);

// message-passing to complete sum of boundary cells.
void sendAndRecieve(BInfo& B);

// Normalize boundary cells.
void BNormalize(Grid& grid, BInfo& B, Options& opt);
