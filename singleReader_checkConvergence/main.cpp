// - This program is created to compute the convergence rate of the basis functions, by computing
//   the distance from the converged basis function in each step.
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
	opt.maxIter = 2000;
	opt.tolerance = -1;
	opt.checkTol = 1;
	Grid g_grid, l_grid;
	Matrix g_mat, l_mat;
	Graph graph;
	Node node;
	BInfo B;
	double* ompSol = NULL;

	Options opt2;
	opt2.maxIter = 20000;
	opt2.tolerance = 0.0001;
	opt2.checkTol = 1;


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

		g_grid.boundary = new bool[g_grid.n_sup]();
		for (int i = 0; i < g_grid.n_sup; i++ ){
			if ( g_grid.celltypes[i]==2 ){
				g_grid.boundary[i] = true;
			}
		}
		delete[] g_grid.celltypes;
		g_grid.celltypes = NULL;
	}


	if (SIZE==1){
		g_grid.updates = new double[g_grid.n_sup];
		setupSupportMatrix(g_grid, g_mat);
		makeMapping(g_grid, g_mat);
		makeBInfo_serial(g_grid, B);

		double* infNorm = new double[opt2.maxIter]();
		double* twoNorm= new double[opt2.maxIter]();


		for (int i = 0; i<opt2.maxIter; i++){
			computeError(g_grid, ompSol, infNorm, twoNorm, i);

			jacobi(g_grid, g_mat, opt2);
			localSum(g_grid, B);
			BNormalize(g_grid, B);

			if (opt2.underTol){
				cout<<"MPI-code converged after "<< i + 1 <<" iterations."<<endl;
				break;
			}
		}

		for (int i = 1; i < opt2.maxIter; i++){
			infNorm[i] /=infNorm[0];
			twoNorm[i]/=twoNorm[0];
		}
		infNorm[0] = 1;
		twoNorm[0] = 1;

		writeErrors(infNorm, twoNorm, opt2.maxIter);

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
		makeBInfo_local(l_grid, node, B);
		makeBInfo_sending(l_grid, B, RANK);

		for (int i = 0; i<opt.maxIter; i++){
			jacobi(l_grid, l_mat,  opt);
			localSum(l_grid, B);
			sendAndRecieve(B);
			BNormalize(l_grid, B);

			if (opt.tolerance > 0 && opt.counter ==0 ){
				MPI_Allreduce(&opt.underTol, &opt.underTol_G, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
				if (opt.underTol_G){
					if (RANK==0){
						cout<<"MPI-code converged after "<< i + 1 <<" iterations."<<endl;
					}
					break;
				}
			}

		}

		gatherBasis(l_grid, g_grid, graph, RANK, SIZE);

		if (RANK==0){
			//computeDiscrepancy(g_grid, ompSol);
		}
	}




	//if (RANK==0 && SIZE!=1 ) writeBasisDistr(graph);
	//if (RANK==0 ) writeBasis(g_grid);
	//if (RANK==0) print(graph);

	writeBasis(g_grid);


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
