#include "data_defs.h"



void ompBasis(Grid& grid, Matrix& mat, double* basis, Options opt, double& setup_time, double& jacobi_time);
void setupCoarseMapping(Grid& grid, Matrix& mat);
void renormalize(Grid& grid,Matrix& mat, double* basis);
