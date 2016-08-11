//The limSend_3.5 code.
//The ompSolver has been changed so that it handles grids where type 1 cells has
// been removed. The change results in a non-optimal implementation.


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
	Options opt;
	opt.maxIter = 100;
	opt.tolerance = 0.001;
	opt.checkTol = 10;
	opt.reNorm = 100;

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
		ompBasis(grid, mat, ompSolution, opt, T1, T2);
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

	MPI_Barrier(MPI_COMM_WORLD);
	T3=MPI_Wtime()-start;

	start = MPI_Wtime();
	for ( int i=1; i<=opt.maxIter; i++){
		check_convergence = (i % opt.checkTol) == 0 && opt.tolerance > 0;
		if (!check_convergence){
			for (int j=0; j<nBasis; j++){ basisArray[j].jacobiProduct(); }
			typeTwo.localSum();
			if (SIZE!=1) typeTwo.sendAndRecieve();
			typeTwo.TTupdate();
			if (i%opt.reNorm == 0) renormalize(basisArray,nBasis,RANK,SIZE,grid.N);
		}
		else {
			underTol = true;
			for (int j=0; j<nBasis; j++){ basisArray[j].jacobiProduct_ConvCheck(); }
			typeTwo.localSum();
			if (SIZE!=1) typeTwo.sendAndRecieve();
			typeTwo.TTupdate();
			if (i%opt.reNorm == 0) renormalize(basisArray,nBasis,RANK,SIZE,grid.N);
			if (SIZE!=1){
				MPI_Allreduce(&underTol, &underTol_G, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
				if (underTol_G){
					if (RANK==0) cout<<"limSend_3.0 converged after "<<i<<" iterations."<<endl;
					break;
				}
			}
			else if (underTol){
				cout<<"limSend_3.0 converged after "<<i<<" iterations."<<endl;
				break;
			}
		}
		if (i == opt.maxIter && RANK==0){
			//cout<<"limSend_3.0 did not converge after "<<i<<" iterations."<<endl;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	T4=MPI_Wtime()-start;


// Gather basis functions on thread 0 and check discrepancy****************** //
	if (RANK==0) limSendSolution = new double[grid.n_basis];
	gatherBasis(grid, basisDistr, basisArray,  RANK,  SIZE, limSendSolution);
	if (RANK==0) computeDiscrepancy( grid, limSendSolution, ompSolution);

	if (RANK==0){
		//cout<<setprecision(6)<<T2<<"\t"<<T4<<endl;
	}


	for (int i = 0; i<grid.n_basis; i++){
		//cout<<limSendSolution[i]<<"\t"<<ompSolution[i]<<endl;
	}

	//grid.print();
	//mat.print();


	delete [] initBasis;
	delete [] limSendSolution;
	delete [] ompSolution;
	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		MPI_Barrier(MPI_COMM_WORLD);