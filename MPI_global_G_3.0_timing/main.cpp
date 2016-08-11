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
#include <stdexcept>
using namespace std;

int RANK, SIZE, TAG = 100;
const string infile = "/home/shomeb/f/fredjoha/Desktop/master-code/MRST/mrst-core/examples/data/tmp/basis_custom/input/";

/* *****************************MAIN***************************************** */
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

	// ***********************Read from file*********************************** //
	Grid grid;
	Matrix mat;
	if(readInfo(&grid, infile)){return 1;};
	if(readMatrix(&grid, &mat, infile)){return 1;};
	double globalBasis[grid.n_basis];
	if(readBasisOperator(&grid, globalBasis, infile)){return 1;};

	Options opt;
	opt.maxIter = 1000;
	opt.tolerance = -1;
	double start, start2;  // timing variable

	double temp1 = 0; double temp2 = 0;
	// *************************get OMP solution******************************* //
	double* ompSolution = NULL;
	if (RANK==-1){
		start = MPI_Wtime();
		ompSolution = new double[grid.n_basis];
		for (int i=0; i<grid.n_basis; i++){ ompSolution[i] = globalBasis[i]; }
		ompBasis(grid, mat, ompSolution, opt, temp1, temp2);
		//cout<<"OMP time:\t"<<MPI_Wtime()-start<<endl;
	}




	// *******************MPI_global_G_3.0 starts****************************** //
		double T1 = 0, T2 = 0;
		start = MPI_Wtime();

		int basisDistr[SIZE+1];
		for (int i = 0; i< SIZE; i++){ basisDistr[i] =  threadStart(i, SIZE, grid.M); }
		basisDistr[SIZE] = grid.M;

		int firstBasis = basisDistr[RANK];				// Number of the first basis function
		int lastBasis = basisDistr[RANK+1]-1;			// Number of the last basis function
		int nBasis = lastBasis-firstBasis+1;		  // Number of basis functions on the current thread

		Basis  basisArray[nBasis];
		for (int i = 0 ; i<nBasis; i++){
			Basis temp(grid, mat, globalBasis, firstBasis+i);
			basisArray[i] = temp;
			basisArray[i].makeLocalMapping_NEW();
		}

		// Making  H for type 2 cells
		map<int, set<int> > H;
		makeH(grid,H, RANK, basisDistr);


		// Making the TypeTwo structure
		TypeTwo typeTwo(basisArray, nBasis, SIZE, RANK, basisDistr, H);
		typeTwo.setupSending();

		MPI_Barrier(MPI_COMM_WORLD);
		T1+= MPI_Wtime() - start;
		start = MPI_Wtime();

		for (int i =0; i<opt.maxIter; i++){
			for (int j = 0; j<nBasis; j++){ basisArray[j].jacobiProduct(); }
			typeTwo.localSum();
			typeTwo.sendAndRecieve();
			typeTwo.TTupdate();
		}

		MPI_Barrier(MPI_COMM_WORLD);
		T2+= MPI_Wtime() - start;
		if (RANK==0)cout<<T1<<endl<<T2/opt.maxIter<<endl<<endl;

		// ********Gather basis functions on thread 0 and check discrepancy******** //
		//double global_GSolution[grid.n_basis];
		//gatherBasis(grid, basisDistr, basisArray,  RANK,  SIZE, global_GSolution);
		//if (RANK==0) computeDiscrepancy( grid, global_GSolution, ompSolution  );





	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		double start = MPI_Wtime();
//		MPI_Barrier(MPI_COMM_WORLD);
