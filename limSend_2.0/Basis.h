#pragma once

#include "data_defs.h"

class Basis{
public:
	double omega; 			// Relaxation factor
	double tolerance; 	// tolerance;

	int number;					// number of the basis function
	int size; 					// number of cells in the support region
	int maxRow;					// max nonzero coefficients in one row

	double* values;			// pressure values
	int* support; 			// global numbering of cells in the support region
	int* celltypes;			// correspinding celltypes

	double* update; 		// update values
	int* mat_indices; 	// indices of matrix to each cell
	double* mat_coef;		// corresponding matrix coefficients

	Basis(Grid&, Matrix&, double* basis, int number, Options& opt);
	Basis(){values = NULL; support = NULL; celltypes = NULL;
					update = NULL; mat_indices = NULL; mat_coef = NULL;  }
	~Basis(){delete[] mat_indices; delete[] mat_coef;  delete[] update; }


	// Shallow copying and swapping. rhs will be deleted after.
	Basis& operator=(Basis& rhs);

	void makeBasis(Grid&, Matrix&, double* basis, int number, Options& opt);
	void makeLocalMapping();
	bool jacobiProduct(bool check_convergence);


	void print();
	void printPressure();
};
