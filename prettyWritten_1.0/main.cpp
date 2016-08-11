// singleReader_6.0 rewritten to be more representable / easier to read.

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
	int RANK, SIZE;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

	Options opt;
	opt.maxIter = 100;
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
		double t1, t2;
		ompSol = new double[g_grid.n_sup];
		for ( int i=0; i<g_grid.n_sup; i++){ ompSol[i] = g_grid.basis[i]; }
		ompBasis_new(g_grid, g_mat, ompSol, opt, t1, t2);
	}

	if (SIZE==1){
		g_grid.updates = new double[g_grid.n_sup];
		setupSupportMatrix(g_grid, g_mat);
		makeMapping(g_grid, g_mat);
		makeTypeTwo_serial(g_grid, tt);

		for (int i = 0; i<opt.maxIter; i++){
			jacobi(g_grid, g_mat,  opt.omega);
			localSum(g_grid, tt);
			TTnormalize(g_grid, tt);
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
			jacobi(l_grid, l_mat,  opt.omega);
			localSum(l_grid, tt);
			sendAndRecieve(l_grid, tt, RANK);
			TTnormalize(l_grid, tt);

		}

		gatherBasis(l_grid, g_grid, graph, RANK, SIZE);
		if (RANK==0){
			computeDiscrepancy(g_grid, ompSol);
		}

	}



	// Wrapup************************************************************************************ //
	MPI_Barrier(MPI_COMM_WORLD);
	delete [] ompSol;
	if (RANK==0) cout<<endl;
	MPI_Finalize();
	return (0);
};
