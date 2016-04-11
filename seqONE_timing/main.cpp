/* A new sequential basis constructor that seems to work and is compared with
	 the OMP one. It uses the data structures Basis and TypeTwo. The Jacobi
	 iterations seems to be slightly faster than  the original OMP code.
	 Main idea behind the program:
	 	- Each basis function performs the Jacobi iteration in the normal way, and
		is represented by its own data structure.
		- After all basis functions has performed its Jacobi computation, a single
		data strcuture named TypeTwo performs the normalization of the type 2 cells.
*/

#include "data_defs.h"
#include "readFromFile.h"
#include "ompSolver.h"
#include "utilities.h"
#include "TypeTwo.h"

#include <iostream>
#include <string>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

int RANK, SIZE, TAG = 100;
const string infile = "/home/shomeb/f/fredjoha/Desktop/master-coding/MRST/mrst-core/examples/data/tmp/basis_custom/input/";




/* *****************************MAIN***************************************** */
/* *****************************MAIN***************************************** */
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

	/* Timing variables */
	const double NRUNS = 100;
	double OMP_totalTime = 0;
	double OMP_JacobiTime = 0;
	double seqONE_totalTime = 0;
	double seqONE_JacobiTime = 0;
	double makeBasis_Time = 0;
	double makeTypeTwo_Time = 0;



/* *************************Read from file*********************************** */
	Grid grid;
	Matrix mat;
	if(readInfo(&grid, infile)){return 1;};
	if(readMatrix(&grid, &mat, infile)){return 1;};
	double* globalBasis = new double[grid.n_basis]();
	if(readBasisOperator(&grid, globalBasis, infile)){return 1;};


	/* *************************get OMP solution******************************* */
	Options opt;
	opt.maxIter = 1;
	opt.tolerance = -1;
	double* ompSolution = new double[grid.n_basis];

	for (int k = 0; k<1; k++){
		for (int i =0; i<grid.n_basis; i++){ ompSolution[i] = globalBasis[i]; }
		ompBasis(grid, mat, ompSolution, opt, OMP_totalTime, OMP_JacobiTime);
	}

	/* *************************My new method********************************** */
	double* seqONESolution;
	for (int k = 0; k<NRUNS; k++){
			seqONESolution = new double[grid.n_basis];
			for (int i =0; i<grid.n_basis; i++){ seqONESolution[i] = globalBasis[i]; }

			double start = MPI_Wtime();

			Basis*  basisArray = new Basis[grid.M];
			for (int i = 0 ; i<grid.M; i++){
				Basis temp(grid, mat, seqONESolution, i);
				basisArray[i] = temp;
				basisArray[i].makeLocalMapping();
			}

			makeBasis_Time+=MPI_Wtime()-start;

			double startTT = MPI_Wtime();
			TypeTwo typeTwo(basisArray, grid.M);

			makeTypeTwo_Time+=MPI_Wtime()-startTT;

			double startJacobi = MPI_Wtime();
			for (int i =0; i<opt.maxIter; i++){
				for (int i = 0; i<grid.M; i++) basisArray[i].jacobiProduct();
				typeTwo.sumAndModify();
			}
			seqONE_JacobiTime+=MPI_Wtime()-startJacobi;

			seqONE_totalTime+=MPI_Wtime()-start;

			delete[] basisArray;
	}
	/*
	cout<<endl;
	cout<<"seqONE Total Time:\t\t"<<endl<<seqONE_totalTime/NRUNS<<endl;
	cout<<"Make Basis Time:\t\t"<<endl<<makeBasis_Time/NRUNS<<endl;
	cout<<"Make TypeTwo Time:\t\t"<<endl<<makeTypeTwo_Time/NRUNS<<endl;
	cout<<"Jacobi Iteration Time:\t\t"<<endl<<seqONE_JacobiTime/NRUNS<<endl<<endl;
	*/
	cout<<endl<<endl;
	cout<<seqONE_totalTime/NRUNS<<endl;
	cout<<makeBasis_Time/NRUNS<<endl;
	cout<<makeTypeTwo_Time/NRUNS<<endl;
	cout<<seqONE_JacobiTime/NRUNS<<endl<<endl;



	/* *************************Compute discrepancy**************************** */
	double error = 0;
	for (int i = 0; i<grid.n_basis; i++){
		double temp = abs(seqONESolution[i] - ompSolution[i]);
		if (temp>error) error = abs(temp);
	}
	cout<<"Discrepancy: "<<error<<endl;

	delete[] seqONESolution;
	delete[] ompSolution;
	delete[] globalBasis;
	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		double start = MPI_Wtime();
//		MPI_Barrier(MPI_COMM_WORLD);
