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
using namespace std;

int RANK, SIZE, TAG = 100;
const string infile = "/home/shomeb/f/fredjoha/Desktop/master-code/MRST/mrst-core/examples/data/tmp/basis_custom/input/";

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

	double global_G_totalTime = 0;
	double global_G_JacobiTime = 0;
	double global_G_setupTime = 0;

	double mainProductTime = 0;
	double localSumTime = 0;
	double sendAndRecieveTime = 0;
	double TTupdateTime = 0;
	double gatherTime = 0;

/* *************************Read from file*********************************** */
	Grid grid;
	Matrix mat;
	if(readInfo(&grid, infile)){return 1;};
	if(readMatrix(&grid, &mat, infile)){return 1;};
	double globalBasis[grid.n_basis];
	if(readBasisOperator(&grid, globalBasis, infile)){return 1;};

	Options opt;
	opt.maxIter = 500;
	opt.tolerance = -1;
	double start; // timing variable


	/* *************************get OMP solution******************************* */
	double  ompSolution[grid.n_basis];
	if (RANK==0){
		for (int k =0; k<1; k++){
			for (int i=0; i<grid.n_basis; i++){ ompSolution[i] = globalBasis[i]; }
			start = MPI_Wtime();
			ompBasis(grid, mat, ompSolution, opt);
			OMP_totalTime+=MPI_Wtime()-start;
		}
	}

	/* *******************MPI_global_G starts********************************** */
		double global_GSolution[grid.n_basis];

		for (int iTT = 0; iTT<NRUNS; iTT++){

				double solutionCopy[grid.n_basis];

				for (int i=0; i<grid.n_basis; i++){
					solutionCopy[i] = globalBasis[i];
					global_GSolution[i] = 0;
				}

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
					Basis temp(grid, mat, solutionCopy, firstBasis+i);
					basisArray[i] = temp;
					basisArray[i].makeLocalMapping();
				}


				TypeTwo typeTwo(basisArray, nBasis, SIZE, RANK, basisDistr);
				typeTwo.setupSending();

				global_G_setupTime+=MPI_Wtime()-start;

				double start2 = MPI_Wtime();

				double t1;
				for (int i =0; i<opt.maxIter; i++){
					t1 = MPI_Wtime();
					for (int i = 0; i<nBasis; i++){ basisArray[i].jacobiProduct(); }
					mainProductTime+=MPI_Wtime() - t1;

					t1 = MPI_Wtime();
					typeTwo.localSum();
					localSumTime+=MPI_Wtime() - t1;

					t1 = MPI_Wtime();
					typeTwo.sendAndRecieve();
					sendAndRecieveTime+=MPI_Wtime() - t1;

					t1 = MPI_Wtime();
					typeTwo.TTupdate();
					TTupdateTime+=MPI_Wtime() - t1;
				}


				global_G_JacobiTime+= MPI_Wtime() - start2;

				t1 = MPI_Wtime();
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

				MPI_Gatherv(basisArray[0].values, recvCounts[RANK], MPI_DOUBLE,
								global_GSolution, recvCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			gatherTime+=MPI_Wtime()-t1;
			 global_G_totalTime+=MPI_Wtime()-start;

		}

	if (RANK==0){

		/*

		cout<<"Total time:\t\t"<<global_G_totalTime/NRUNS<<endl;
		cout<<"Total Jacobi time:\t"<<global_G_JacobiTime/NRUNS<<endl<<endl;

		cout<<"Setup time:\t\t"<<global_G_setupTime/NRUNS<<endl<<endl;

		cout<<"mainProductTime:\t"<<mainProductTime/NRUNS<<endl;

		cout<<"localSumTime:\t\t"<<localSumTime/NRUNS<<endl;

		cout<<"sendAndRecieveTime:\t"<<sendAndRecieveTime/NRUNS<<endl;

		cout<<"TTupdateTime:\t\t"<<TTupdateTime/NRUNS<<endl<<endl;

		cout<<"gatherTime:\t\t"<<gatherTime<<endl<<endl;

		cout<<"Time Discrepancy in %:\t"<<(global_G_totalTime-global_G_setupTime-mainProductTime-localSumTime-sendAndRecieveTime-TTupdateTime-gatherTime)
		/global_G_totalTime*100<<endl;

		*/
		cout<<global_G_totalTime/NRUNS<<endl;
		cout<<global_G_JacobiTime/NRUNS<<endl;
		cout<<global_G_setupTime/NRUNS<<endl<<endl;
		cout<<mainProductTime/NRUNS<<endl;
		cout<<localSumTime/NRUNS<<endl;
		cout<<sendAndRecieveTime/NRUNS<<endl;
		cout<<TTupdateTime/NRUNS<<endl<<endl;
		cout<<gatherTime<<endl<<endl;

	}
	/* *************************Compute discrepancy**************************** */
	double error = 0;
	if (RANK==0){
		for (int i = 0; i<grid.n_basis; i++){
			double temp = abs(global_GSolution[i] - ompSolution[i]);
			if (temp>error) error = temp;
			if (error>1) cout<<"MAJOR ERROR!"<<endl;
		}
		cout<<"Discrepancy: "<<error<<endl;
	}


	MPI_Finalize();
};


/* ****** Useful commands ******** */
//		double start = MPI_Wtime();
//		MPI_Barrier(MPI_COMM_WORLD);
