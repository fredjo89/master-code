#include "utilities.h"
#include "data_defs.h"
#include "Basis.h"

#include <mpi.h>
#include <map>
#include <vector>
#include <iostream>
#include <set>
#include <cmath>
using namespace std;


int blockToThread(int block, int SIZE, int M){
  if (block>=M){return -1;}

  int blocks_per_thread = M/SIZE;

  if (M%SIZE==0){ return block/blocks_per_thread; }

  int rest = M - blocks_per_thread*SIZE;

  int blockDivider = (blocks_per_thread+1)*rest -1;

  if (block<=blockDivider) return block/(blocks_per_thread+1);

  int blocksLeft = block - blockDivider;

  return blockDivider/(blocks_per_thread+1)+1 + (blocksLeft-1)/blocks_per_thread;


}

void makeH(Grid& grid, map<int, set<int> >& H, int RANK, int*basisDistr){
  int basisNumber = 0;
  for (int i = 0; i<grid.n_basis; i++){
    if (i==grid.offsets[basisNumber+1]) basisNumber++;
    if (grid.celltypes[i]==2){
      H[grid.support[i]].insert(basisNumber);
    }
  }

  set<int> basis;
  for (int i = basisDistr[RANK]; i<basisDistr[RANK+1]; i++){
    basis.insert(i);
  }

  map<int, set<int> > newH;

  map<int, set<int> >::iterator it;
  set<int>::iterator it2;
  set<int>::iterator it3;

  for (it=H.begin(); it!=H.end(); it++){
    set<int> temp = it->second;

    for (it2=basis.begin(); it2!=basis.end(); it2++){
      if (temp.find(*it2)!=temp.end() ){
        for (it3=temp.begin(); it3!=temp.end(); it3++){
          if (basis.find(*it3)==basis.end()){
            newH[it->first]=it->second;
            break;
          }
        }
      }
    }
  }
  H=newH;
}

int threadStart(int RANK, int SIZE, int M ){
  int basisEach = M / SIZE;
  int rest = M - basisEach * SIZE;
  if (RANK < rest ) return (basisEach+1)*RANK;
  else return (basisEach+1)*rest + basisEach*(RANK-rest);
}

void printArray(int* a, int l){
  for (int i =0; i<l; i++){
    cout<<a[i]<<endl;
  }
}

void printArray(double* a, int l){
  for (int i =0; i<l; i++){
    cout<<a[i]<<endl;
  }
}

void printH(map<int, set<int> >& H){
  cout<<"______________________________"<<endl<<"Printing H"<<endl;
  cout<<"SIZE:\t"<<H.size()<<endl;
  map<int, set<int> >::iterator it;
  set<int>::iterator it2;
  for (it = H.begin(); it!=H.end(); it++){
    cout<<it->first<<":\t";

    for (it2 = it->second.begin(); it2!= it->second.end(); it2++){
      cout<<*it2<<"\t";
    }
    cout<<endl;
  }
    cout<<"______________________________"<<endl;
}

void gatherBasis(Grid& grid, int* basisDistr, Basis* basisArray, int RANK, int SIZE,
double* global_GSolution){
  // ********Gather basis functions on thread 0 and check discrepancy******** //
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
}



void computeDiscrepancy(Grid& grid, double* GSol,  double* ompSol  ){
  double error = 0;
  bool mError = false;
  for (int i = 0; i<grid.n_basis; i++){
    double temp = abs(GSol[i] - ompSol[i]);
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
