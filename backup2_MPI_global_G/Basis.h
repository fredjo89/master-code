#pragma once

#include "data_defs.h"
#include <iostream>
using namespace std;

class Basis{
public:
	int number;					// number of the basis function
	int size; 					// number of cells in the support region
	int maxRow;					// max nonzero coefficients in one row

	double* values;			// pressure values
	int* support; 			// global numbering of cells in the support region
	int* celltypes;			// correspinding celltypes

	double* update; 		// update values
	int* mat_indices; 	// indices of matrix to each cell
	double* mat_coef;		// corresponding matrix coefficients

	Basis(Grid&, Matrix&, double* basis, int number);
	Basis(){values = NULL; support = NULL; celltypes = NULL; mat_indices = NULL;
															mat_coef = NULL; update = NULL; }
	~Basis(){delete[] mat_indices; delete[] mat_coef;  delete[] update; }


	// Shallow copying and swapping. rhs will be deleted after.
	Basis& operator=(Basis& rhs);

	void makeLocalMapping();
	void jacobiProduct();

	void print();
	void printPressure();
};
