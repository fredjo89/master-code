#include "ompSolver.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
using namespace std;

/* ************************************************************************** */
void ompBasis_new(Grid& grid, Matrix& mat, double* basis, Options opt, double& setupTime, double& iterTime){

  double start = MPI_Wtime();

  const int NF = grid.N;
  const int NC = grid.M;
  const int NB = grid.n_sup;
  const int MR = mat.maxRow;

  mat.sup_index = new int[NB*MR];
  mat.sup_coef = new double[NB*MR];

  setupCoarseMapping(grid, mat);

  // Cache is used for intermediate storage of basis
  double* cache = new double[NB];
  // Used for sum of updates for renormalization
  double* update_sum = new double[NF]();
  // Storage for the update values
  double* basis_update = new double[NB]();

  int maxIter = opt.maxIter;
  int checkTol = opt.checkTol;
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
            int jx = mat.sup_index[j*MR + k];
            if (jx != -1) {
              cache[j] +=  mat.sup_coef[j*MR + k]* basis[jx];
            }
          }
        }
      }

      for (int i = 0; i<NB; i++){
        int type = grid.celltypes[i];
        if (type!=1){
          basis_update[i] = basis[i]-omega*cache[i];
          if (check_convergence && fabs(cache[i]) > tolerance
                                        && grid.celltypes[i] == 0) done = false;
        }
        if (type==2){
          update_sum[grid.support[i]] += basis_update[i];
        }
      }

      for (int i = 0; i<NB; i++){
        if (grid.celltypes[i] == 2) {
          basis_update[i]/=update_sum[grid.support[i]];
        }
      }


      for (int i = 0; i < NB; i++) {
        basis[i] = basis_update[i];
  		}


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
  delete[] mat.sup_index; mat.sup_index = NULL;
  delete[] mat.sup_coef; mat.sup_coef = NULL;
}

/* ************************************************************************** */
void setupCoarseMapping(Grid& grid, Matrix& mat){
  mat.n_sup = grid.n_sup;

  const int NF = grid.N;
  const int NC = grid.M;
  const int NB = grid.n_sup;
  const int MR = mat.maxRow;

  for (int i = 0; i < NB*MR; i++) {
    mat.sup_index[i] = -1;
    mat.sup_coef[i] = 0;
  }

  for (int coarseIx = 0; coarseIx < NC; coarseIx++) {
    int start = grid.offsets[coarseIx];
    int stop  = grid.offsets[coarseIx + 1];

    for (int j = start; j < stop; j++) {
      int c = grid.support[j];
      for (int k = 0; k < MR; k++) {
        int c_n = mat.mat_index[c*MR + k];
        int l;
        bool done = false;
        if(c_n == -1) break;
        if(1==1){
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
        }else{
          for (l = start; l < stop; l++) {
            if (grid.support[l] == c_n) {
              done = true;
              break;
            }
          }
        }


        if(done){
          mat.sup_index[j*MR + k] = l;
          mat.sup_coef[j*MR + k] = mat.mat_coef[c*MR + k];
        }
      }
    }
  }
}

/* ************************************************************************** */
void OMPrenormalize(Grid& grid,Matrix& mat, double* basis){
	int NB = grid.n_sup;

  double* sumval = new double[grid.N]();
	for(int i = 0; i < NB; i++){ sumval[grid.support[i]] += basis[i]; }
	for(int i = 0; i < NB; i++){ basis[i] = basis[i]/sumval[grid.support[i]]; }
  delete[] sumval;
}
