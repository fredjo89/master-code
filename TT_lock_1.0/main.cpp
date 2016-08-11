//The limSend_3.0 code


#include "data_defs.h"
#include "readFromFile.h"
#include "ompSolver.h"
#include "utilities.h"
#include "TypeTwo.h"
#include "Basis.h"

#include <iostream>
#include <mpi.h>
#include <cmath>
#include <iomanip>
using namespace std;

// MAIN********************************************************************** //
int main(int argc, char* argv[])
{
	int RANK, SIZE, TAG = 50;
	double start, start2, start3; 									// timing variable
	double T1 = 0, T2 = 0, T3 = 0, T4 = 0;					// timing variable

	int check_convergence, localDone, globalDone;		// check convergence variables
	bool underTol, underTol_G;
	Options opt, opt2;
	opt.maxIter = 2000;
	opt.tolerance = -1;
	opt.checkTol = 1;


	opt2.maxIter = 10000;
	opt2.tolerance = 0.0000000001;
	opt2.checkTol = 1;


	Grid grid;
	Matrix mat;
	double* initBasis = NULL;
	double* ompSolution = NULL;
	double* ompSolution_Converged = NULL;
	double* limSendSolution = NULL;

	const string infile = argv[1];
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);


// Read from file************************************************************ //
	if(readInfo(&grid, infile)) return 1; if(readMatrix(&grid, &mat, infile)) return 1;
	initBasis = new double[grid.n_basis]; if(readBasisOperator(&grid, initBasis, infile)) return 1;

// get converged OMP solution************************************************ //
	if (RANK==0 ){
		ompSolution_Converged = new double[grid.n_basis];
		for ( int i=0; i<grid.n_basis; i++){ ompSolution_Converged[i] = initBasis[i]; }
		ompBasis(grid, mat, ompSolution_Converged, opt2, T1, T2);
	}




// get OMP solution********************************************************** //
	if (RANK==0 ){
		ompSolution = new double[grid.n_basis];
		for ( int i=0; i<grid.n_basis; i++){ ompSolution[i] = initBasis[i]; }
		ompBasis(grid, mat, ompSolution, opt, T1, T2);
	}

// MPI_global_G_2.0 starts*************************************************** //

	if (RANK==0) limSendSolution = new double[grid.n_basis];

	int basisDistr[SIZE+1];
	for ( int i=0; i< SIZE; i++) basisDistr[i] = threadStart(i, SIZE, grid.M);
	basisDistr[SIZE] = grid.M;
	int firstBasis = basisDistr[RANK];
	int nBasis = basisDistr[RANK+1]-firstBasis;
	Basis  basisArray[nBasis];
	for ( int i=0 ; i<nBasis; i++) basisArray[i].makeBasis(grid, mat, initBasis, firstBasis+i, opt, underTol);
	TypeTwo typeTwo(grid, basisArray, nBasis, SIZE, RANK, basisDistr);
	typeTwo.setupSending();


	gatherBasis(grid, basisDistr, basisArray,  RANK,  SIZE, limSendSolution);
	computeDiscrepancy( grid, limSendSolution, ompSolution_Converged);
	for ( int i=1; i<=opt.maxIter; i++){
		if ( i%10==0 ){
			for (int j=0; j<nBasis; j++){ basisArray[j].jacobiProduct(); }
			typeTwo.localSum(); typeTwo.TTupdate();
			if (i%opt.reNorm == 0) renormalize(basisArray,nBasis,RANK,SIZE,grid.N);
		}
		else{
			for (int j=0; j<nBasis; j++){ basisArray[j].jacobiProduct_lockG(); }
		}
		if (i%5==0){
			gatherBasis(grid, basisDistr, basisArray,  RANK,  SIZE, limSendSolution);
			computeDiscrepancy( grid, limSendSolution, ompSolution_Converged);
			//cout<<i<<endl;
		}
	}









	delete [] initBasis;
	delete [] limSendSolution;
	delete [] ompSolution;
	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		MPI_Barrier(MPI_COMM_WORLD);
