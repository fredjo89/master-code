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
	int RANK, SIZE;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

	double start1, start2, start3;
	double t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0, t8 = 0;
	Options opt;
	opt.maxIter = 10;
	opt.tolerance = -1;
	opt.checkTol = 1;
	Grid g_grid, l_grid;
	Matrix g_mat, l_mat;
	Graph graph;
	Node node;
	TypeTwo tt;
	double* ompSol = NULL;

	if (RANK==0 ){
		// Read from file************************************************************************ //
		const string infile = argv[1];
		if(readInfo(&g_grid, infile)) return 1;
		if(readMatrix(&g_grid, &g_mat, infile)) return 1;

		// Get OMP solution********************************************************************** //
		ompSol = new double[g_grid.n_sup];
		for ( int i=0; i<g_grid.n_sup; i++){ ompSol[i] = g_grid.basis[i]; }
		ompBasis_new(g_grid, g_mat, ompSol, opt, t1, t2);
		t1 = 0; t2 = 0;
	}

	if (SIZE==1){
		g_grid.updates = new double[g_grid.n_sup];
		setupSupportMatrix(g_grid, g_mat);
		makeMapping(g_grid, g_mat);
		makeTypeTwo_serial(g_grid, tt);

		for (int i = 0; i<opt.maxIter; i++){
			start1 = MPI_Wtime();
			jacobi(g_grid, g_mat,  opt.omega);
			t1 += MPI_Wtime() - start1;

			localSum(g_grid, tt, RANK);
			TTnormalize(g_grid, tt, RANK);
		}

		computeDiscrepancy(g_grid, ompSol);
	}
	else{
		if (RANK==0){
			makeCoarseGraph(g_grid, g_mat, SIZE, graph);
			setupNodes(graph,  node,  RANK);
			distributeProblem(g_grid, l_grid, g_mat, l_mat, graph, node);
		}
		else{
			setupNodes(graph,  node,  RANK);
			distributeProblem(g_grid, l_grid, g_mat, l_mat, graph, node);
		}

		makeMapping(l_grid, l_mat);
		makeTypeTwo_local(l_grid, node, tt);
		makeTypeTwo_sending(l_grid, tt, RANK);

		for (int i = 0; i<opt.maxIter; i++){
			start1 = MPI_Wtime();
			jacobi(l_grid, l_mat,  opt.omega);
			t1 += MPI_Wtime() - start1;
			localSum(l_grid, tt, RANK);
			sendAndRecieve(l_grid, tt, RANK);
			TTnormalize(l_grid, tt, RANK);

		}


		gatherBasis(l_grid, g_grid, graph, RANK, SIZE);
		if (RANK==0){
			computeDiscrepancy(g_grid, ompSol);
		}

	}

	//if (RANK==0 && SIZE!=1 ) writeBasisDistr(graph);
	//if (RANK==0 ) writeBasis(g_grid);

	if (RANK==0) cout<<t1<<endl;


	// Wrapup************************************************************************************ //
	MPI_Barrier(MPI_COMM_WORLD);
	delete [] ompSol;
	if (RANK==0) cout<<endl;


	MPI_Finalize();
	return (0);

};



/* ****** Useful commands ******** */
//		MPI_Barrier(MPI_COMM_WORLD);
// 		//cout<<setprecision(6);
//		Initiallizing zero array: int array[length] = {0};
