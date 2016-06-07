//The limSend_1.0 code


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

int RANK, SIZE, TAG = 100;

// *****************************MAIN***************************************** //
int main(int argc, char* argv[])
{
	const string infile = argv[1];
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

	// ***********************Read from file*********************************** //
	Grid grid;
	Matrix mat;
	if(readInfo(&grid, infile)) return 1;
	if(readMatrix(&grid, &mat, infile)) return 1;
	double globalBasis[grid.n_basis];
	if(readBasisOperator(&grid, globalBasis, infile)) return 1;

	Options opt;
	opt.maxIter = 10000;
	opt.tolerance = -1;
	double start; // timing variable

	double* ompSolution = NULL;
	double* limSendSolution = NULL;

	// *************************get OMP solution******************************* //
	if (RANK==0){
		start = MPI_Wtime();
		ompSolution = new double[grid.n_basis];
		for (int i=0; i<grid.n_basis; i++){ ompSolution[i] = globalBasis[i]; }
		ompBasis(grid, mat, ompSolution, opt);
		//cout<<"OMP time:\t"<<MPI_Wtime()-start<<endl;
	}


	// *******************MPI_global_G_2.0 starts****************************** //

	// Create the basisDistr vector.
	int basisDistr[SIZE+1];
	for (int i = 0; i< SIZE; i++){ basisDistr[i] =  threadStart(i, SIZE, grid.M); }
	basisDistr[SIZE] = grid.M;

	int firstBasis = basisDistr[RANK];					// Number of the first basis function
	int nBasis = basisDistr[RANK+1]-firstBasis;	// Number of basis functions on the current thread

	// Create Basis datastructures.
	Basis  basisArray[nBasis];
	for (int i = 0 ; i<nBasis; i++){
		basisArray[i].makeBasis(grid, mat, globalBasis, firstBasis+i, opt.omega);
		basisArray[i].makeLocalMapping();
	}


	// Create the TypeTwo structure.
	TypeTwo typeTwo(grid, basisArray, nBasis, SIZE, RANK, basisDistr );
	typeTwo.setupSending();

	double TIME = 0;

	// Iteration.
	for (int i =0; i<opt.maxIter; i++){
		start = MPI_Wtime();
		for (int j = 0; j<nBasis; j++){ basisArray[j].jacobiProduct(); }
		TIME += MPI_Wtime() - start;

		typeTwo.localSum();
		typeTwo.sendAndRecieve();
		typeTwo.TTupdate();
		}

		for (int j = 0; j<nBasis; j++){ basisArray[j].print(); }

	// ********Gather basis functions on thread 0 and check discrepancy******** //
	//if (RANK==0) limSendSolution = new double[grid.n_basis];
	//gatherBasis(grid, basisDistr, basisArray,  RANK,  SIZE, limSendSolution);
	//if (RANK==0) computeDiscrepancy( grid, limSendSolution, ompSolution  );

	cout<<TIME<<endl;


	delete [] limSendSolution;
	delete [] ompSolution;
	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		double start = MPI_Wtime();
//		MPI_Barrier(MPI_COMM_WORLD);
