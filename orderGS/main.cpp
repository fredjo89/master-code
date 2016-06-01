/* Ordered Jacobi iteration, NOT GS!!
Step towards ordered GS iteration.
Dones NOT work for large problems! Memory trouble!
*/

#include "data_defs.h"
#include "readFromFile.h"
#include "ompSolver.h"
#include "utilities.h"
#include "orderGS_Solver.h"


#include <iostream>
#include <string>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <set>
using namespace std;

int RANK, SIZE, TAG = 100;
const string infile = "/home/shomeb/f/fredjoha/Desktop/master-code/MRST/mrst-core/examples/data/tmp/basis_custom/input/";




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
	double globalBasis[grid.n_basis];
	if(readBasisOperator(&grid, globalBasis, infile)){return 1;};

	int M = grid.M;
	int N = grid.N;
	int n_basis = grid.n_basis;
	int maxRow = mat.maxRow;
	int* offsets = grid.offsets;

	Options opt;
	opt.maxIter = 1000;
	opt.tolerance = 0.001;

	//****************************get OMP solution**********************************

	double ompSolution[n_basis];
	for (int i =0; i<n_basis; i++){ ompSolution[i] = globalBasis[i]; }

	ompBasis(grid, mat, ompSolution, opt);



	//****************************Ordered Jacobi************************************

	double GS_Solution[n_basis];
	for (int i =0; i<n_basis; i++){ GS_Solution[i] = globalBasis[i]; }

	orderGS_Basis(grid, mat, GS_Solution, opt);

	//*******************************Compute discrepancy****************************



	double error = 0;
	for (int i = 0; i<grid.n_basis; i++){
		bool mError = false;
		double temp = abs(GS_Solution[i] - ompSolution[i]);
		if (temp>error) error = temp;
		if (error>1 && !mError){
			cout<<"MAJOR ERROR"<<endl;
			mError = true;
		}
		if (temp!=temp && !mError){
			cout<<"NaN OCCURED!"<<endl;
			mError = true;
		}
	}
	cout<<"Discrepancy: "<<error<<endl;







	cout<<endl;
	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		double start = MPI_Wtime();
//		MPI_Barrier(MPI_COMM_WORLD);
