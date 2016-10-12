// A program that creates prolongation operator used by the MsRSB method, and can use combinations
// of MPI processes and OpenMP threads.
// Input: file path to where the problem data is located.
// Output: prolongation operator / basis functions.
// Iteration options can be changed by modified the struct Options.
// The input data is the same as Olav Møyners prolongation-construction operator which is contained
// in MRST, except the cells of type 1 have been removed from support regions.

#include "data_defs.h"
#include "readFromFile.h"
#include "distribution.h"
#include "jacobi.h"
#include "utilities.h"

#include <mpi.h>
#include <iomanip>
using namespace std;

// ********************************************************************************************** //
int main(int argc, char* argv[]){
	// Timing variables
	double t_total = MPI_Wtime();
	double t_read = 0, t_setup = 0, t_iter = 0;

	// Initiallizing MPI
	int RANK; 	// Rank of process.
	int SIZE;	// Number of processes.
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
	MPI_Comm_rank (MPI_COMM_WORLD, &RANK);

	// Decleare variables (defined in data_defs.h)
	Options opt;	// Holds options like tolerance and max iterations.
	// g_grid / g_mat (g -> "global"): only used by process 0.
	// l_grid / l_mat (l -> "local"): local problem used by each process (not used if SIZE==1).
	Grid g_grid, l_grid;
	Matrix g_mat, l_mat;
	// BInfo: holds info about union boundary. used by each process.
	BInfo B;
	// Graph: holds info about data distribution. only used by process 0. (not used if SIZE==1).
	Graph graph;
	// Node: holds info about distributed data. used by each process. (not used if SIZE==1).
	Node node;

	// Read from file**************************************************************************** //
	if (RANK==0){
		t_read = MPI_Wtime();
		// Read data
		const string fn = argv[1];
		if(readInfo(&g_grid, fn)) return 1;
		if(readMatrix(&g_grid, &g_mat, fn)) return 1;
		t_read = MPI_Wtime() - t_read;
	}

	// If #MPIprocesses==1 (serial or pure OpenMP application)
	if (SIZE==1){
		t_setup = MPI_Wtime();

		g_grid.updates = new double[g_grid.n_sup];	// Buffer storing updated basis.
		// setupSupportMatrix() + makeMapping() -> creates sup_index and sup_coef.
		setupSupportMatrix(g_grid, g_mat);
		makeMapping(g_grid, g_mat);
		makeBInfo_serial(g_grid, B); // Creates BInfo.

		t_setup = MPI_Wtime() - t_setup;
		t_iter = MPI_Wtime();

		// Iteration
		for (int i = 0; i < opt.maxIter; i++){
			jacobi(g_grid, g_mat, opt);	// Update operator
			localSum(g_grid, B);		// Sum Union Boundary
			BNormalize(g_grid, B, opt);	// Normalize Union Boundary

			// Convergence check
			if (opt.underTol){
				cout<<"Operator converged after "<< i + 1 <<" iterations."<<endl;
				break;
			}
			else if (i==opt.maxIter-1){
				cout<<"Operator did not converge after "<< i + 1 <<" iterations."<<endl;
			}
		}
		t_iter = MPI_Wtime() - t_iter;
	}
	// If #MPIprocesses!=1 (MPI application)
	else{
		t_setup = MPI_Wtime();

		if (RANK==0){
			makeCoarseGraph(g_grid, g_mat, SIZE, graph); // Make graph determening basis partition.
			setupNodes(graph,  node,  RANK);	// Create Node (Includes first message-passing)
			distributeProblem(g_grid, l_grid, g_mat, l_mat, graph, node);	// Distribute problem.
		}
		else{
			setupNodes(graph,  node,  RANK);	// Create Node (includes first message-passing)
			distributeProblem(g_grid, l_grid, g_mat, l_mat, graph, node); // Receive problem.
		}
		makeMapping(l_grid, l_mat);
		// makeBInfo_local() + makeBInfo_sending() -> creates BInfo
		makeBInfo_local(l_grid, node, B);
		makeBInfo_sending(l_grid, B, RANK);

		t_setup = MPI_Wtime() - t_setup;
		t_iter = MPI_Wtime();

		// Iteration.
		for (int i = 0; i<opt.maxIter; i++){
			jacobi(l_grid, l_mat, opt); // Update Operator.
			localSum(l_grid, B);	// Sum Unon boundary found locally.
			sendAndRecieve(B);		// Send / recieve Union boundary info.
			BNormalize(l_grid, B, opt); // Normalize Union Boundary.

			// Convergence check
			if (opt.tolerance > 0 && opt.counter == 0 ){
				MPI_Allreduce(&opt.underTol, &opt.underTol_G, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
				if (opt.underTol_G){
					if (RANK==0){
						cout<<"Operator converged after "<< i + 1 <<" iterations."<<endl;
					}
					break;
				}
			}
			else if (i==opt.maxIter-1){
				if (RANK==0){
					cout<<"Operator did not converge after "<< i + 1 <<" iterations."<<endl;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		t_iter = MPI_Wtime() - t_iter;
		// Gather operator on process 0.
		gatherBasis(l_grid, g_grid, graph, RANK, SIZE);

	}

	// Write operator to file.
	if (RANK==0){
		const string fn = argv[1];
		writeBasis(g_grid, fn);
	}

	// Wrapup************************************************************************************ //
	t_total = MPI_Wtime() - t_total;
	if (RANK==0){
		cout.precision(4);
		cout<<"______________________________"<<endl
		<<"Total time:\t"<<t_total<<endl
		<<"Read time:\t"<<t_read<<endl
		<<"Setup time:\t"<<t_setup<<endl
		<<"Iteration time:\t"<<t_iter<<endl
		<<"______________________________"<<endl;
	}

	MPI_Finalize();
	return (0);
};
