#include "orderGS_Solver.h"
#include "ompSolver.h"
#include "TypeTwo.h"
#include "Basis.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <set>
using namespace std;


//********************************************************************************
void orderGS_Basis(Grid& grid, Matrix& mat, double* basis, Options opt){

  int N = grid.N;               // Number of cells
  int M = grid.M;               // Number of blocks
  int n_basis = grid.n_basis;   // Number of elements in basis array
  int maxRow = mat.maxRow;      // Max element in a row
  int* offsets = grid.offsets;  // offsets array

  int* orders = new int [n_basis];         // Order array

  getOrders( grid, mat, orders);

  int* celltypes = new int[n_basis];

  for (int i = 0; i<M; i++){
    for (int j = offsets[i]; j<offsets[i+1]; j++){
      int currentCell = orders[j];
      int ctype = -5;
      for (int k = offsets[i]; k<offsets[i+1]; k++){
        if (currentCell==grid.support[k]) {
          ctype = grid.celltypes[k];
          break;
        }
      }
      celltypes[j] = ctype;
    }
  }

  int mapOne[n_basis]; // Map from initial order to new order
  int mapTwo[n_basis]; // Map from  new order to initial order
  for (int i = 0; i<M; i++){
    for (int j =offsets[i]; j<offsets[i+1]; j++){
      int k;
      for ( k = offsets[i]; k<offsets[i+1]; k++){
        if (orders[k]==grid.support[j])
        break;
      }
      mapOne[k] = j;
      mapTwo[j] = k;
    }
  }

  double orderBasis [n_basis];
  for (int i = 0; i<n_basis; i++){
    orderBasis[i] = basis[mapOne[i]];
  }


  delete [] grid.support;
  delete[] grid.celltypes;
  grid.support = orders;
  grid.celltypes = celltypes;

  // Create basis structures
  Basis*  basisArray = new Basis[grid.M];
  for (int i = 0 ; i<grid.M; i++){
    Basis temp(grid, mat, orderBasis, i);
    basisArray[i] = temp;
    basisArray[i].makeLocalMapping();
  }

  TypeTwo typeTwo(basisArray, grid.M);

  // Jacobi iteration
  for (int i =0; i<opt.maxIter; i++){
    bool done = true;

    for (int i = 0; i<grid.M; i++){
      if ( !basisArray[i].jacobiProduct(opt.tolerance) ) done = false;
    }
    typeTwo.sumAndModify();

    if((i + 1) % opt.reNorm == 0 && i > 1) renormalize(grid, mat, orderBasis);


    if (done){
      cout<<"Ordered Jacobi converged after "<< i+1 <<" iterations."<<endl;
      break;
    }
    if (i+1 == opt.maxIter){
      cout<<"Ordered Jacobi did not converge after "<<i+1<<" iterations."<<endl;
    }

  }

  // Reorder for output
  for (int i = 0; i<n_basis; i++){
    basis[i] = orderBasis[mapTwo[i]];
  }

}
//********************************************************************************

void getOrders(Grid& grid, Matrix& mat, int* orders){

  int M = grid.M;
  int N = grid.N;
  int n_basis = grid.n_basis;

  int myVec[N*M];

	for (int it = 0; it<M; it++){
		set<int> mySet;
		myVec[it*N] = grid.centers[it];
		mySet.insert(grid.centers[it]);

		int currentSize = 1;
		int doneCells = 0;
		int i = 0;
		while (mySet.size()<N ){
			for (int j = doneCells; j<currentSize; j++){
				int index = myVec[it*N + doneCells]*mat.maxRow;
				for ( int k = index; k<index + mat.maxRow; k++){
					int insert = mat.j_index[k];
					if (insert!= -1 && mySet.find(insert)==mySet.end()  ){
						myVec[ it*N + currentSize] = insert;
						mySet.insert(insert);
						currentSize++;
					}
				}
				doneCells++;
			}
		}
	}

	int m = 0; // current Basis number
	for (int i = 0; i<N*M; i++){
		if (i % N == 0 && i!=0) m++;

		int currentCell = myVec[i];
		bool inSupport = false;
		for (int j = grid.offsets[m]; j<grid.offsets[m+1]; j++ ){
			if (currentCell==grid.support[j]){
			inSupport = true;
			break;
			}
		}
		if (!inSupport) myVec[i] = -5;
	}

	int k = 0;
	for (int i = 0; i<N*M; i++){
		if (myVec[i]!=-5){
			orders[k]= myVec[i];
			k++;
		}
	}

}
//******************************************************************************** 


void makeOrderStructures(Grid& grid, Matrix& mat, Grid& orderGrid, Matrix& orderMat, int* orders){

  // Make Grid
  int M = grid.M;
  int N = grid.N;
  int n_basis = grid.n_basis;
  int maxRow = mat.maxRow;

  orderGrid.N = N;
  orderGrid.M = M;
  orderGrid.n_basis = n_basis;

  int* offsets = new int[M+1];
  for (int i = 0; i<M+1; i++) offsets[i]=grid.offsets[i];

  int* centers = new int[M];
  for (int i = 0; i<M; i++) centers[i]=grid.centers[i];

  orderGrid.offsets = offsets;
  orderGrid.centers = centers;

  int* celltypes = new int[n_basis]{};

  for (int i = 0; i<M; i++){
    for (int j = offsets[i]; j<offsets[i+1]; j++){
      int currentCell = orders[j];
      int ctype = -5;
      for (int k = offsets[i]; k<offsets[i+1]; k++){
        if (currentCell==grid.support[k]) {
          ctype = grid.celltypes[k];
          break;
        }
      }
      celltypes[j] = ctype;
    }
  }

  orderGrid.support = orders;
  orderGrid.celltypes = celltypes;

  // Make Matrix
  orderMat.N = N;
  orderMat.maxRow = maxRow;
  orderMat.n_basis = n_basis;

  double* newConn = new double [N*maxRow];
  int* newj_index = new int[N*maxRow];
  for (int i = 0; i<N*maxRow; i++){
    newConn[i] = mat.conn[i];
    newj_index[i] = mat.j_index[i];
  }

  orderMat.conn = newConn;
  orderMat.j_index = newj_index;

  int* loc_index = new int[n_basis*maxRow];
  double* loc_conn = new double [n_basis*maxRow];


  int m = 0;
  for (int i = 0; i<n_basis; i++){
    if (i==offsets[m+1]){
      m++;
    }

    int cellNum = orderGrid.support[i];
    int* tempOne = &orderMat.j_index[maxRow*cellNum];
    double* tempTwo = &orderMat.conn[maxRow*cellNum];
    int k = 0 ;
    for (int j = maxRow*i; j<maxRow*(i+1); j++){
      if (tempOne[k]==-1) loc_index[j] = -1;
      else loc_index[j] = tempOne[k]; // + offsets[m] ;
      loc_conn[j] = tempTwo[k];
      k++;
    }
  }

  orderMat.loc_index = loc_index;
  orderMat.loc_conn = loc_conn;



}
