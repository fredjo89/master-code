#pragma once
#include <iostream>
#include <vector>
#include <map>
using namespace std;

/* ************************************************************************** */
struct Grid
{
	int N;				// number of fine cells
	int M;				// number of coarse blocks
	int n_basis;

	int* support; 		// n_basis elements
	int* celltypes; 	// n_basis elements
	int* offsets; 		// M + 1 elements

	Grid() {support = NULL; celltypes = NULL; offsets = NULL;}
	~Grid(){ delete[] support; delete[] celltypes; delete[] offsets;  }
	void print();
};
/* ************************************************************************** */
struct Matrix
{
	int 		N;
	int 		maxRow;
	int 		n_basis;

	double* conn;
	int* 		j_index;

	int *  loc_index;
 double *  loc_conn;

	Matrix(){conn = NULL; j_index = NULL; loc_index = NULL; loc_conn = NULL; }
	~Matrix() { delete[] conn; delete[] j_index; delete[] loc_index;
																								delete[] loc_conn;}
	void print();
};

/* ************************************************************************** */
struct Options
{
	double tolerance;
	int maxIter;
	double omega;
	int reNorm;
	int checkTol;

	Options(){tolerance = 0.001; maxIter = 500; omega =  (double) 2/3;
																				reNorm = 100; checkTol = 1; }
	void print();
};





/* ************************************************************************** */
