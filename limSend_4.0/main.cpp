//	The limSend_4.0 code
// 	Two changes from limSend.3.0
// 	1. 	The normalization is done in the new and improved way, also for the ompSolver
// 	2. 	The ompSolver has been slightly modified to handle grids where type 1
// 			cells has been removed.
//	-	Also, the help function printPoperator() has been created to help visualize
//		the prolongation operator.
//	-	No convergence check
//	-	No renormalization (No longer required)

#include "data_defs.h"
#include "readFromFile.h"
#include "ompSolver.h"
#include "utilities.h"
#include "TypeTwo.h"
#include "Basis.h"

#include <metis.h>			///Metis package

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
	Options opt;
	opt.maxIter = 100;
	opt.tolerance = -1;
	opt.checkTol = -1;

	Grid grid;
	Matrix mat;
	double* initBasis = NULL;
	double* ompSolution = NULL;
	double* limSendSolution = NULL;

	const string infile = argv[1];
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);


// Read from file************************************************************ //
	if(readInfo(&grid, infile)) return 1;
	if(readMatrix(&grid, &mat, infile)) return 1;
	initBasis = new double[grid.n_basis];
	if(readBasisOperator(&grid, initBasis, infile)) return 1;


// get OMP solution********************************************************** //
	if (RANK==0 ){
		ompSolution = new double[grid.n_basis];
		for ( int i=0; i<grid.n_basis; i++){ ompSolution[i] = initBasis[i]; }
		ompBasis_new(grid, mat, ompSolution, opt, T1, T2);
	}

// MPI_global_G_2.0 starts*************************************************** //

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	int basisDistr[SIZE+1];
	for ( int i=0; i< SIZE; i++) basisDistr[i] = threadStart(i, SIZE, grid.M);
	basisDistr[SIZE] = grid.M;
	int firstBasis = basisDistr[RANK];
	int nBasis = basisDistr[RANK+1]-firstBasis;

	Basis  basisArray[nBasis];
	for ( int i=0 ; i<nBasis; i++){
		basisArray[i].makeBasis(grid, mat, initBasis, firstBasis+i, opt, underTol);
	}

	TypeTwo typeTwo(grid, basisArray, nBasis, SIZE, RANK, basisDistr);
	typeTwo.setupSending();

	for ( int i=1; i<=opt.maxIter; i++){
		for (int j=0; j<nBasis; j++){ basisArray[j].jacobiProduct(); }
		typeTwo.localSum_two();
		typeTwo.sendAndRecieve();
		typeTwo.TTupdate_two();
	}

// Gather basis functions on thread 0 and check discrepancy****************** //
	if (RANK==0) limSendSolution = new double[grid.n_basis];
	gatherBasis(grid, basisDistr, basisArray,  RANK,  SIZE, limSendSolution);
	if (RANK==0) computeDiscrepancy( grid, limSendSolution, ompSolution);

	MPI_Barrier(MPI_COMM_WORLD);
	if (RANK==0) cout<<MPI_Wtime()-start<<endl;

	delete [] initBasis;
	delete [] limSendSolution;
	delete [] ompSolution;
	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		MPI_Barrier(MPI_COMM_WORLD);
// 		//cout<<setprecision(6);
