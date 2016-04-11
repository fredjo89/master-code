#include "data_defs.h"



void ompBasis(Grid& grid, Matrix& mat, double* basis, Options opt, double& OMP_totalTime, double& OMP_JacobiTime);
void setupCoarseMapping(Grid& grid, Matrix& mat);
void renormalize(Grid& grid,Matrix& mat, double* basis);
