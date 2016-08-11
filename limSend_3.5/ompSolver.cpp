#include "ompSolver.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
using namespace std;


/* ************************************************************************** */
void ompBasis(Grid& grid, Matrix& mat, double* basis, Options opt, double& setupTime, double& iterTime){

  double start = MPI_Wtime();

  const int NF = grid.N;
  const int NC = grid.M;
  const int NB = grid.n_basis;
  const int MR = mat.maxRow;

  mat.loc_index = new int[NB*MR];
  mat.loc_conn = new double[NB*MR];

	setupCoarseMapping(grid, mat);


  // Cache is used for intermediate storage of basis
  double* cache = new double[NB];
  // Used for sum of updates for renormalization
  double* update_sum = new double[NF]{};
  // Storage for the update values
  double* basis_update = new double[NB]{};

  int maxIter = opt.maxIter;
  int checkTol = opt.checkTol;
  int reNorm = opt.reNorm;
  double tolerance = opt.tolerance;
  double omega = opt.omega;

  setupTime += MPI_Wtime() - start;
  start = MPI_Wtime();

/* ******************JACOBI ITERATION *************************************** */
  for (int iter = 1; iter <=maxIter; iter++){
      bool check_convergence = (iter % checkTol) == 0 && tolerance > 0;
      bool done = check_convergence;

      /* **************STEP 1: Main Iteration********************** */
      for (int coarseIx = 0; coarseIx < NC; coarseIx++) {
        int start = grid.offsets[coarseIx];
        int stop  = grid.offsets[coarseIx+1];
        for (int j = start; j < stop; j++) {
          cache[j] = basis[j];
          for (int k = 0; k < MR; k++) {
            int jx = mat.loc_index[j*MR + k];


            if (jx != -1) cache[j] +=  mat.loc_conn[j*MR + k]* basis[jx];
          }
        }
      }

  		/* **************STEP 2: Modification ONE ********************* */
  		for (int i = 0; i < NB; i++) {
  			int type = grid.celltypes[i];
  			if(type == 0) basis_update[i] += omega*cache[i];
  			else if (type == 2) {
  				basis_update[i] += omega*cache[i];
  				update_sum[grid.support[i]] += omega*cache[i];
  			}
  			else basis_update[i] = 0;
  		}

  	/* **************STEP 3: Modification TWO ********************* */
  		for (int i = 0; i < NB; i++) {
  			if (grid.celltypes[i] == 2) {
  				double us = update_sum[grid.support[i]];
  				basis_update[i] = (basis_update[i] - basis[i] * us) / (1.0 + us);
  			}
  		}

  		/* ************STEP 4: Update Basis Functions**************** */
  		for (int i = 0; i < NB; i++) {
  			basis[i] -= basis_update[i];
  			if (check_convergence && fabs(basis_update[i]/omega) > tolerance
                                      && grid.celltypes[i] == 0) done = false;
  			basis_update[i] = 0;
  		}

  		/* ***************STEP 5: Renormalize************************ */
  		if(iter % reNorm == 0 ) OMPrenormalize(grid, mat, basis);

  		/* **********STEP 6: set update sum to zero****************** */
  		for (int i = 0; i < NF; i++) update_sum[i] = 0;

  		/* **************STEP 7: Check if done*********************** */
  		if (done) {
        cout<<"OMP-code converged after "<<iter<<" iterations."<<endl;
        break;
      }
      else if (iter == maxIter){
        //cout<<"OMP-code did not converge after "<<iter<<" iterations."<<endl;
      }
  }

  iterTime += MPI_Wtime() - start;


/* ******************END JACOBI ITERATION *********************************** */
  delete[] basis_update;
  delete[] update_sum;
  delete[] cache;
  delete[] mat.loc_index; mat.loc_index = NULL;
  delete[] mat.loc_conn; mat.loc_conn = NULL;
}



/* ************************************************************************** */
void setupCoarseMapping(Grid& grid, Matrix& mat){
  mat.n_basis = grid.n_basis;

  int l;
  bool done;


  const int NF = grid.N;
  const int NC = grid.M;
  const int NB = grid.n_basis;
  const int MR = mat.maxRow;

  for (int i = 0; i < NB*MR; i++) {
    mat.loc_index[i] = -1;
    mat.loc_conn[i] = 0;
  }

  for (int coarseIx = 0; coarseIx < NC; coarseIx++) {
    int start = grid.offsets[coarseIx];
    int stop  = grid.offsets[coarseIx + 1];

    for (int j = start; j < stop; j++) {
        int c = grid.support[j];
        for (int k = 0; k < MR; k++) {
            int c_n = mat.j_index[c*MR + k];

            if(c_n == -1) break;

            done = false;


            // This branch assumes support is sorted.
            if (c_n > c) {
                for (l = j; l < stop; l++) {
                    if (grid.support[l] == c_n) {
                        done = true;
                        break;
                    }
                }
            }
            else {
                for (l = j - 1; l >= start; l--) {
                    if (grid.support[l] == c_n) {
                        done = true;
                        break;
                    }
                }
            }



            if(done){
                mat.loc_index[j*MR + k] = l;
                mat.loc_conn[j*MR + k] = mat.conn[c*MR + k];
            }
        }
    }
  }

}

/* ************************************************************************** */
void OMPrenormalize(Grid& grid,Matrix& mat, double* basis){
	int NB = grid.n_basis;

  double* sumval = new double[grid.N]{};
	for(int i = 0; i < NB; i++){ sumval[grid.support[i]] += basis[i]; }
	for(int i = 0; i < NB; i++){ basis[i] = basis[i]/sumval[grid.support[i]]; }
  delete[] sumval;
}
