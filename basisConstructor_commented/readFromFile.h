// Read write data from file.

#pragma once

#include "data_defs.h"
#include <fstream>
#include <cstring>

int readInfo(Grid * grid, std::string file_path);
int readMatrix(Grid * grid, Matrix * mat, std::string file_path);
int readBasisOperator(Grid * grid, double * basis, std::string file_path);
int writeBasis(Grid& g_grid, std::string file_path);
int openFlatFile(std::string file_path, std::ifstream * file);
int openFlatFile(std::string file_path, std::ofstream * file);
