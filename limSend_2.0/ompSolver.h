#include "data_defs.h"

void ompBasis(Grid& grid, Matrix& mat, double* basis, Options opt);
void setupCoarseMapping(Grid& grid, Matrix& mat);
void OMPrenormalize(Grid& grid,Matrix& mat, double* basis);
