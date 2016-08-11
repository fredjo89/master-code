#pragma once

#include "data_defs.h"

class Basis{
public:
	double omega; 			// Relaxation factor
	double tolerance; 	// tolerance;
	bool*  underTol;		// check if update is under tolerance.

	int number;					// number of the basis function
	int size; 					// number of cells in the support region
	int maxRow;					// max nonzero coefficients in one row

	double* values;			// pressure values
	int* support; 			// global numbering of cells in the support region
	int* celltypes;			// correspinding celltypes

	double* update; 		// update values
	int* mat_indices; 	// indices of matrix to each cell
	double* mat_coef;		// corresponding matrix coefficients

	Basis();
	~Basis();

	void makeBasis(Grid&, Matrix&, double* basis, int number, Options& opt, bool& underTol);

	void jacobiProduct();
	void jacobiProduct_lockG(); 	// lock global boundary cells.
	void jacobiProduct_ConvCheck();


	void print();
	void printPressure();
};
