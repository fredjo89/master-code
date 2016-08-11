#include "data_defs.h"

void ompBasis(Grid& grid, Matrix& mat, double* basis, Options opt, double& setupTime, double& iterTime);
void setupCoarseMapping(Grid& grid, Matrix& mat);
void OMPrenormalize(Grid& grid,Matrix& mat, double* basis);

void ompBasis_new(Grid& grid, Matrix& mat, double* basis, Options opt, double& setupTime, double& iterTime);
