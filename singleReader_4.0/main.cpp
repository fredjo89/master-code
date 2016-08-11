//	The singleReader_1.0 code
//	Goals for the code:
//	-	Only once processor will read from file and distribute to the other processors
//	-	METIS will be used to divide basis functions between processors in a good way
//
//	A requirement for the code to work will probably be that each support region is sortet from
//	cellnumber with lowest to highest number. If this is not the case, then a sorting can be applied
//	initially.
//	Assumes no cells of type 1.

#include "data_defs.h"
#include "readFromFile.h"
#include "ompSolver.h"
#include "distribution.h"
#include "jacobi.h"
#include "utilities.h"

#include <iostream>
#include <mpi.h>
#include <cmath>
#include <iomanip>
using namespace std;

// MAIN****************************************************************************************** //
int main(int argc, char* argv[])
{
	// Declear variables************************************************************************* //
	int RANK, SIZE, TAG = 50;
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

	double start1, start2, start3;
	double T0 = 0, T1 = 0, T2 = 0, T3 = 0, T4 = 0, T5 = 0, T6 = 0, T7 = 0, T8 = 0;
	Options opt;
	opt.maxIter = 1;
	opt.tolerance = -1;
	opt.checkTol = 1;
	Grid globalGrid, localGrid;
	Matrix globalMat, localMat;
	Graph graph;
	Node node;
	TypeTwo typeTwo;
	double* ompSolution = NULL;

	if (RANK==0 ){
		// Read from file************************************************************************ //
		const string infile = argv[1];
		if(readInfo(&globalGrid, infile)) return 1;
		if(readMatrix(&globalGrid, &globalMat, infile)) return 1;
		globalMat.n_basis = globalGrid.n_basis;


		// Get OMP solution********************************************************************** //
		ompSolution = new double[globalGrid.n_basis];
		for ( int i=0; i<globalGrid.n_basis; i++){ ompSolution[i] = globalGrid.basis[i]; }
		ompBasis_new(globalGrid, globalMat, ompSolution, opt, T1, T2);
		T1= 0; T2 = 0;
	}

	if (SIZE==1){
		start1 = MPI_Wtime();

		globalGrid.updates = new double[globalGrid.n_basis];
		setupSupportMatrix(globalGrid, globalMat);
		makeMapping(globalGrid, globalMat);
		makeTypeTwo_serial(globalGrid, typeTwo);

		for (int i = 0; i<opt.maxIter; i++){
			jacobi(globalGrid, globalMat,  opt.omega);
			localSum(globalGrid, typeTwo, RANK);
			TTnormalize(globalGrid, typeTwo, RANK);

		}

		computeDiscrepancy(globalGrid, ompSolution);

	}
	else{
		MPI_Barrier(MPI_COMM_WORLD);
		start1 = MPI_Wtime();
		if (RANK==0){
			makeCoarseGraph(globalGrid, globalMat, SIZE, graph);
			setupNodes(graph,  node,  RANK);
			reOrderGrid(globalGrid, graph.basisDistr, SIZE);
			distributeProblem(globalGrid, localGrid, globalMat, localMat, graph, node);

		}
		else{
			setupNodes(graph,  node,  RANK);
			distributeProblem(globalGrid, localGrid, globalMat, localMat, graph, node);
		}

		makeMapping(localGrid, localMat);
		makeTypeTwo_local(localGrid, node, typeTwo);
		makeTypeTwo_sending(localGrid, typeTwo);

		for (int i = 0; i<opt.maxIter; i++){
			jacobi(localGrid, localMat,  opt.omega);
			localSum(localGrid, typeTwo, RANK);
			sendAndRecieve(localGrid, typeTwo, RANK);
			TTnormalize(localGrid, typeTwo, RANK);

		}

		if (RANK==0){
			localGrid.basis = NULL;
			localGrid.offsets = NULL;
		}
		else{
			localGrid.celltypes = NULL;
			localGrid.offsets = NULL;
		}

		start1 = MPI_Wtime();

		gatherBasis(localGrid, globalGrid, graph, RANK, SIZE);
		if (RANK==0){
			reOrderBasis(globalGrid);
			computeDiscrepancy(globalGrid, ompSolution);
		}
	}





	// Wrapup************************************************************************************ //
	MPI_Barrier(MPI_COMM_WORLD);
	delete [] ompSolution;
	if (RANK==0) cout<<endl;
	MPI_Finalize();
};






/* ****** Useful commands ******** */
//		MPI_Barrier(MPI_COMM_WORLD);
// 		//cout<<setprecision(6);
//		Initiallizing zero array: int array[length] = {0};
//		Initiallizing allocated memory to zero: int* array = new int[length]();
