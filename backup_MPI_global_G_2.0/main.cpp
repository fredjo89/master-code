/* Taking seqONE a step further.
*/

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

/* *************************Read from file*********************************** */
	Grid grid;
	Matrix mat;
	if(readInfo(&grid, infile)){return 1;};
	if(readMatrix(&grid, &mat, infile)){return 1;};
	double globalBasis[grid.n_basis];
	if(readBasisOperator(&grid, globalBasis, infile)){return 1;};

	Options opt;
	opt.maxIter = 1;
	opt.tolerance = -1;
	double start; // timing variable

	/* *************************get OMP solution******************************* */
	double  ompSolution[grid.n_basis];
	if (RANK==0){
		start = MPI_Wtime();
		for (int i=0; i<grid.n_basis; i++){ ompSolution[i] = globalBasis[i]; }
		ompBasis(grid, mat, ompSolution, opt);
		cout<<"OMP time:\t"<<MPI_Wtime()-start<<endl;
	}

	/* *******************MPI_global_G starts********************************** */
	MPI_Barrier(MPI_COMM_WORLD);
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
			basisArray[i].makeLocalMapping();
		}


		// Making  H for type 2 cells
		map<int, set<int> > H;
		makeH(grid,H, RANK, basisDistr);



		TypeTwo typeTwo(basisArray, nBasis, SIZE, RANK, basisDistr, H);
		typeTwo.setupSending();


		for (int i =0; i<opt.maxIter; i++){
			for (int i = 0; i<nBasis; i++){ basisArray[i].jacobiProduct(); }
			typeTwo.localSum();
			typeTwo.sendAndRecieve();
			typeTwo.TTupdate();
		}

	/* ********Gather basis functions on thread 0 and check discrepancy******** */
	int recvCounts [SIZE]; int displs [SIZE];
	int temp = 0; int k = 0;
	for (int i = 0; i < grid.M; i++){
		if (i==basisDistr[k+1]){
			recvCounts[k] = temp;
			temp = (grid.offsets[i+1]-grid.offsets[i]);
			k++;
		}
		else{ temp+=(grid.offsets[i+1]-grid.offsets[i]); }
	}
	recvCounts[SIZE-1] = temp;
	displs[0] = 0;
	for (int i =1; i<SIZE; i++){ displs[i] = displs[i-1] + recvCounts[i-1]; }
	double global_GSolution[grid.n_basis];
	MPI_Gatherv(basisArray[0].values, recvCounts[RANK], MPI_DOUBLE,
					global_GSolution, recvCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (RANK==0) cout<<"TIME:\t\t"<<MPI_Wtime()- start<<endl;

	/* *************************Compute discrepancy**************************** */
	double error = 0;
	if (RANK==0){
		bool mError = false;
		for (int i = 0; i<grid.n_basis; i++){
			double temp = abs(global_GSolution[i] - ompSolution[i]);
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
		cout<<"Discrepancy:\t"<<error<<endl;
	}

	if (RANK==0) cout<<endl<<endl;

	MPI_Comm graph_comm;
	int index[4] = {1, 3, 5, 6};
	int edges[6] = {1, 0, 2, 1, 3, 2};
	int reorder = 1;
	MPI_Graph_create(MPI_COMM_WORLD, SIZE, index, edges, reorder, &graph_comm);

	if (RANK==0){
		int count = -1;

		int currentRank = 0;

		MPI_Graph_neighbors_count(graph_comm, currentRank, &count);
		cout<<count<<endl;

		int neighbors[count];
		MPI_Graph_neighbors(graph_comm, currentRank,count, neighbors);

		for (int i = 0; i<count; i++){
			cout<<"Neighbor "<<i<<":\t"<<neighbors[i]<<endl;
		}

		int INDEX = 4;
		int EDGES = 6;
		int index2[4];
		int edges2[6];

		MPI_Graph_get(graph_comm,INDEX, EDGES,index2,edges2 );

		for (int i = 0; i<INDEX; i++)cout<<index2[i]<<"\t"; cout<<endl;
		for (int i = 0; i<EDGES; i++)cout<<edges2[i]<<"\t"; cout<<endl;

		int topo_type;
		MPI_Topo_test(graph_comm, &topo_type);
		cout<<topo_type<<endl;
	}

	int sendNumber[4]; for (int i = 0; i<4; i++){ sendNumber[i]= RANK*RANK+1;}
	int recvNumber[10] {};

	MPI_Alltoall(&sendNumber, 1, MPI_INT, &recvNumber, 1, MPI_INT, graph_comm );







	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		double start = MPI_Wtime();
//		MPI_Barrier(MPI_COMM_WORLD);