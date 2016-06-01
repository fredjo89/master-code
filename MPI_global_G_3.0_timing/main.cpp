


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

	// ***********************Read from file*********************************** //
	Grid grid;
	Matrix mat;
	if(readInfo(&grid, infile)){return 1;};
	if(readMatrix(&grid, &mat, infile)){return 1;};
	double globalBasis[grid.n_basis];
	if(readBasisOperator(&grid, globalBasis, infile)){return 1;};

	Options opt;
	opt.maxIter = 99;
	opt.tolerance = -1;
	double start, start2;  // timing variable

	double OMP_setup_time  = 0;
	double OMP_jacobi_time = 0;

	double MPI_setup_time = 0;
	double MPI_jacobi_time = 0;

	int IT = 1000;

	// *************************get OMP solution******************************* //
	double  ompSolution[grid.n_basis];
	for (int k = 0; k<1; k++){
		double temp1, temp2;
		if (RANK==0){
			for (int i=0; i<grid.n_basis; i++){ ompSolution[i] = globalBasis[i]; }
			ompBasis(grid, mat, ompSolution, opt, temp1, temp2);

			OMP_setup_time+=temp1;
			OMP_jacobi_time+=temp2;
		}
	}


	double Basis_time = 0;
	double H_time = 0;
	double TT_time = 0;

	double product_time = 0;
	double localSum_time = 0;
	double send_time = 0;
	double TTupdate_time = 0;

	// *******************MPI_global_G_2.0 starts****************************** //
	for (int k = 0; k<IT; k++){
		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();

		int basisDistr[SIZE+1];
		for (int i = 0; i< SIZE; i++){ basisDistr[i] =  threadStart(i, SIZE, grid.M); }
		basisDistr[SIZE] = grid.M;

		int firstBasis = basisDistr[RANK];				// Number of the first basis function
		int lastBasis = basisDistr[RANK+1]-1;			// Number of the last basis function
		int nBasis = lastBasis-firstBasis+1;		  // Number of basis functions on the current thread


		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();
		Basis  basisArray[nBasis];
		for (int i = 0 ; i<nBasis; i++){
			Basis temp(grid, mat, globalBasis, firstBasis+i);
			basisArray[i] = temp;
			basisArray[i].makeLocalMapping();
		}
		MPI_Barrier(MPI_COMM_WORLD);
		Basis_time+=MPI_Wtime()-start;

		// Making  H for type 2 cells
		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();
		map<int, set<int> > H;
		makeH(grid,H, RANK, basisDistr);
		MPI_Barrier(MPI_COMM_WORLD);
		H_time+=MPI_Wtime()-start;


		// Making the TypeTwo structure
		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();
		TypeTwo typeTwo(basisArray, nBasis, SIZE, RANK, basisDistr, H);
		typeTwo.setupSending();
		MPI_Barrier(MPI_COMM_WORLD);
		TT_time+=MPI_Wtime()-start;



		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();
		for (int i =0; i<opt.maxIter; i++){

			MPI_Barrier(MPI_COMM_WORLD);
			start2 = MPI_Wtime();
			MPI_Barrier(MPI_COMM_WORLD);
			for (int j = 0; j<nBasis; j++){ basisArray[j].jacobiProduct(); }
			MPI_Barrier(MPI_COMM_WORLD);
			product_time+=MPI_Wtime()-start2;
			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			start2 = MPI_Wtime();
			MPI_Barrier(MPI_COMM_WORLD);
			typeTwo.localSum();
			MPI_Barrier(MPI_COMM_WORLD);
			localSum_time+=MPI_Wtime()-start2;
			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			start2 = MPI_Wtime();
			MPI_Barrier(MPI_COMM_WORLD);
			typeTwo.sendAndRecieve();
			MPI_Barrier(MPI_COMM_WORLD);
			send_time+=MPI_Wtime()-start2;
			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			start2 = MPI_Wtime();
			MPI_Barrier(MPI_COMM_WORLD);
			typeTwo.TTupdate();
			MPI_Barrier(MPI_COMM_WORLD);
			TTupdate_time+=MPI_Wtime()-start2;
			MPI_Barrier(MPI_COMM_WORLD);

		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_jacobi_time+=MPI_Wtime()-start;


		// ********Gather basis functions on thread 0 and check discrepancy******** //
		double global_GSolution[grid.n_basis];
		gatherBasis(grid, basisDistr, basisArray,  RANK,  SIZE, global_GSolution);

		//if (RANK==0) cout<<"TIME:\t\t"<<MPI_Wtime()- start1<<endl;
		//if (RANK==0) computeDiscrepancy( grid, global_GSolution, ompSolution  );

	}

	//cout<<TT_time<<endl;


	 Basis_time /=IT;
	 H_time /=IT;
	 TT_time /=IT;

	 MPI_jacobi_time /=IT;
	 product_time /=IT;
	 localSum_time /=IT;
	 send_time /=IT;
	 TTupdate_time /=IT;

	 MPI_jacobi_time /=opt.maxIter;
	 product_time /=opt.maxIter;
	 localSum_time /=opt.maxIter;
	 send_time /=opt.maxIter;
	 TTupdate_time /=opt.maxIter;


	 double diff = (MPI_jacobi_time - product_time - localSum_time - send_time - TTupdate_time);

	 /*
	 double Basis_time = 0;
	 double H_time = 0;
	 double TT_time = 0;

	 double product_time = 0;
	 double localSum_time = 0;
	 double send_time = 0;
	 double TTupdate_time = 0;
	 */

	if (RANK==0){

		cout << Basis_time <<endl;
		cout <<  H_time<<endl;
		cout <<  TT_time <<endl<<endl;


		cout << product_time <<endl;
		cout <<  localSum_time <<endl;
		cout <<  send_time <<endl;
		cout <<  TTupdate_time <<endl<<endl;


		/*
		cout << "Basis_time:\t\t"<<Basis_time <<endl;
		cout << "H_time:\t\t\t"<< H_time<<endl;
		cout << "TT_time:\t\t"<< TT_time <<endl<<endl;

		//cout << "MPI_jacobi_time:\t"<< MPI_jacobi_time  <<endl;
		//cout << "diff:\t\t\t"<< diff <<endl<<endl;
		cout << "product_time:\t\t"<< product_time <<endl;
		cout << "localSum_time:\t\t"<< localSum_time <<endl;
		cout << "send_time:\t\t"<< send_time <<endl;
		cout << "TTupdate_time:\t\t"<< TTupdate_time <<endl<<endl;
		*/

	}


// cout << ":\t"<<  <<endl;








	MPI_Finalize();
	cout<<endl<<endl;
};


/* ****** Useful commands ******** */
//		double start = MPI_Wtime();
//		MPI_Barrier(MPI_COMM_WORLD);
