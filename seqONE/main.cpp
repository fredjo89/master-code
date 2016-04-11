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
	for (int i =0; i<grid.n_basis; i++){ ompSolution[i] = globalBasis[i]; }

	double start = MPI_Wtime();
	ompBasis(grid, mat, ompSolution, opt);
	cout<<"OMP total time:\t\t\t"<<MPI_Wtime()-start<<endl;

	/* *************************My new method********************************** */
	start = MPI_Wtime();

	Basis*  basisArray = new Basis[grid.M];
	for (int i = 0 ; i<grid.M; i++){
		Basis temp(grid, mat, globalBasis, i);
		basisArray[i] = temp;
		basisArray[i].makeLocalMapping();
	}

	TypeTwo typeTwo(basisArray, grid.M);

	for (int i =0; i<opt.maxIter; i++){
		for (int i = 0; i<grid.M; i++) basisArray[i].jacobiProduct();
		typeTwo.sumAndModify();
	}

	cout<<"seqOne Total Time:\t\t"<<MPI_Wtime()-start<<endl;


	/* *************************Compute discrepancy**************************** */
	double error = 0;
	for (int i = 0; i<grid.n_basis; i++){
		double temp = abs(globalBasis[i] - ompSolution[i]);
		if (temp>error) error = abs(temp);
	}
	cout<<"Discrepancy: "<<error<<endl;

	delete[] ompSolution;
	delete[] globalBasis;
	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		double start = MPI_Wtime();
//		MPI_Barrier(MPI_COMM_WORLD);
