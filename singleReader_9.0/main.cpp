// - Convergence test fixed
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
	double start1, start2, start3;

	double t_read =0, t_total =0, t_setup =0, t_iter =0;
	double t_jacobi =0, t_localSum =0, t_sendRcv =0, t_TTnorm =0;
	double t_makeCG =0, t_setupN =0, t_distr =0, t_makeMap =0, t_makeTT_local =0, t_makeTT_send =0;
	double t_gather =0;
	double t_convTest = 0;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

	Options opt;
	opt.maxIter = 400;
	opt.tolerance = 0.001;
	opt.checkTol = 10;
	Grid g_grid, l_grid;
	Matrix g_mat, l_mat;
	Graph graph;
	Node node;
	BInfo B;

	if (RANK==0 ){
		// Read from file************************************************************************ //
		start1 = MPI_Wtime();

		const string infile = argv[1];
		if(readInfo(&g_grid, infile)) return 1;
		if(readMatrix(&g_grid, &g_mat, infile)) return 1;

		g_grid.boundary = new bool[g_grid.n_sup]();
		for (int i = 0; i < g_grid.n_sup; i++ ){
			if ( g_grid.celltypes[i] == 2 ){
				g_grid.boundary[i] = true;
			}
		}
		delete[] g_grid.celltypes;
		g_grid.celltypes = NULL;

		t_read = MPI_Wtime() - start1;
	}


	if (SIZE==1){
		start1 = MPI_Wtime();

		g_grid.updates = new double[g_grid.n_sup];
		setupSupportMatrix(g_grid, g_mat);
		makeMapping(g_grid, g_mat);
		makeBInfo_serial(g_grid, B);

		t_setup = MPI_Wtime() - start1;

		for (int i = 0; i<opt.maxIter; i++){
			jacobi(g_grid, g_mat, opt);
			localSum(g_grid, B);
			BNormalize(g_grid, B, opt);

			if (opt.underTol){
				cout<<"MPI-code converged after "<< i + 1 <<" iterations."<<endl;
				break;
			}
		}

		t_iter = MPI_Wtime() - start1 - t_setup;
		t_total = MPI_Wtime() - start1;
	}
	else{
		start1 = MPI_Wtime();

		if (RANK==0){
			start2 = MPI_Wtime();
			makeCoarseGraph(g_grid, g_mat, SIZE, graph);
			t_makeCG = (MPI_Wtime() - start2);

			start2 = MPI_Wtime();
			setupNodes(graph,  node,  RANK);
			t_setupN = (MPI_Wtime() - start2);

			start2 = MPI_Wtime();
			distributeProblem(g_grid, l_grid, g_mat, l_mat, graph, node);
			t_distr = (MPI_Wtime() - start2);
		}
		else{
			setupNodes(graph,  node,  RANK);
			distributeProblem(g_grid, l_grid, g_mat, l_mat, graph, node);
		}

		start2 = MPI_Wtime();
		makeMapping(l_grid, l_mat);
		t_makeMap = (MPI_Wtime() - start2);

		start2 = MPI_Wtime();
		makeBInfo_local(l_grid, node, B);
		t_makeTT_local = (MPI_Wtime() - start2);

		start2 = MPI_Wtime();
		makeBInfo_sending(l_grid, B, RANK);
		t_makeTT_send = (MPI_Wtime() - start2);

		MPI_Barrier(MPI_COMM_WORLD);
		t_setup = MPI_Wtime() - start1;


		for (int i = 0; i<opt.maxIter; i++){

			if (RANK == 0) cout<<i<<endl;

			start2 = MPI_Wtime();
			jacobi(l_grid, l_mat,  opt);
			t_jacobi += (MPI_Wtime() - start2);

			start2 = MPI_Wtime();
			localSum(l_grid, B);
			t_localSum += (MPI_Wtime() - start2);

			start2 = MPI_Wtime();
			sendAndRecieve(B);
			t_sendRcv += (MPI_Wtime() - start2);

			start2 = MPI_Wtime();
			BNormalize(l_grid, B, opt);
			t_TTnorm += (MPI_Wtime() - start2);

			start2 = MPI_Wtime();
			if (opt.tolerance > 0 && opt.counter ==0 ){
				MPI_Allreduce(&opt.underTol, &opt.underTol_G, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
				if (opt.underTol_G){
					if (RANK==0){
						cout<<"MPI-code converged after "<< i + 1 <<" iterations."<<endl;
					}
					break;
				}
			}
			t_convTest+= (MPI_Wtime() - start2);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		t_iter = MPI_Wtime() - start1 - t_setup;

		start2 = MPI_Wtime();
		gatherBasis(l_grid, g_grid, graph, RANK, SIZE);
		t_gather = MPI_Wtime() - start2;

		t_total = MPI_Wtime() - start1;

	}

	if (RANK==0){
		cout.precision(4);
		cout<<"______________________________"<<endl;
		cout<<"SIZE:\t"<<SIZE<<endl;
		cout<<"N:\t"<<g_grid.N<<endl;
		cout<<"M:\t"<<g_grid.M<<endl;
		cout<<"n_sup:\t"<<g_grid.n_sup<<endl;
		cout<<"n_sendCells:\t"<<graph.n_sendCells<<endl<<endl;


		cout<<"t_read:\t"<<t_read<<endl;
		cout<<"t_total:\t"<<t_total<<endl<<endl;

		cout<<"t_setup:\t"<<t_setup<<endl<<endl;
		cout<<"t_makeCG:\t"<<t_makeCG<<endl;
		cout<<"t_setupN:\t"<<t_setupN<<endl;
		cout<<"t_distr:\t"<<t_distr<<endl;
		cout<<"t_makeMap:\t"<<t_makeMap<<endl;
		cout<<"t_makeTT_local:\t"<<t_makeTT_local<<endl;
		cout<<"t_makeTT_send:\t"<<t_makeTT_send<<endl<<endl;

		cout<<"t_iter:\t"<<t_iter<<endl;

		cout<<"t_jacobi:\t"<<t_jacobi<<endl;
		cout<<"t_localSum:\t"<<t_localSum<<endl;
		cout<<"t_sendRcv:\t"<<t_sendRcv<<endl;
		cout<<"t_TTnorm:\t"<<t_TTnorm<<endl<<endl;

		cout<<"t_gather:\t"<<t_gather<<endl<<endl;

		cout<<"t_convTest:\t"<<t_convTest<<endl;

		cout<<"______________________________"<<endl;
	}

	writeBasis(g_grid);


	// Wrapup************************************************************************************ //
	MPI_Barrier(MPI_COMM_WORLD);
	if (RANK==0) cout<<endl;
	MPI_Finalize();
	return (0);
};



/* ****** Useful commands ******** */
//		MPI_Barrier(MPI_COMM_WORLD);
// 		//cout<<setprecision(6);
//		Initiallizing zero array: int array[length] = {0};
