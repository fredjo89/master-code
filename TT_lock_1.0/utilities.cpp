#include "utilities.h"
#include "data_defs.h"
#include "Basis.h"

#include <mpi.h>
#include <map>
#include <vector>
#include <iostream>
#include <set>
#include <cmath>
#include <algorithm>
using namespace std;

/* ************************************************************************** */
void renormalize(Basis* basisArray, int nBasis, int RANK, int SIZE, int N){
  double* sendBuff = new double[N]{};
  double* recvBuff  = new double[N]{};
  for (int i = 0; i<nBasis; i++){
    for (int j = 0; j<basisArray[i].size; j++){
      if (basisArray[i].celltypes[j]!=1) {
        sendBuff[basisArray[i].support[j]]+=basisArray[i].values[j];
      }
    }
  }

  if (SIZE!=0) MPI_Allreduce(sendBuff, recvBuff, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int i = 0; i<nBasis; i++){
    for (int j = 0; j<basisArray[i].size; j++){
      if (basisArray[i].celltypes[j]!=1) {
        basisArray[i].values[j]/=recvBuff[basisArray[i].support[j]];
      }
    }
  }
  delete[] sendBuff;
  delete[] recvBuff;
}


/* ************************************************************************** */
void makeH(Grid& grid, map<int, set<int> >& H, map<int,int>& typeTwoInfo, int RANK, int* basisDistr){
  int N = grid.N;
  int M = grid.M;
  int n_basis = grid.n_basis;
  int* support = grid.support;
  int* celltypes = grid.celltypes;
  int* offsets = grid.offsets;

  int maxBasis = basisDistr[RANK+1];
  int minBasis = basisDistr[RANK];
  int TTcells[n_basis];  // cell numbers of type 2 cells will be stored here;
  int bNumbers[n_basis];  // corresponding basis numbers will be stored here;
  int N_TT = 0;         // The total number of type 2 cells will be stored here.

  int power1 = 100;
  int power2 = 100;
  while (M>power1) power1*=10;
  while (N>power2) power2*=10;
  int base = power1*power2;

  int k = 0; // basis number counter;
  for (int i = 0; i<n_basis; i++){
    if (i==offsets[k+1]) k++;
    if (celltypes[i]==2){
      TTcells[N_TT] = support[i]*power1 +k;
      N_TT++;
    }
  }

  std::sort(TTcells, TTcells+N_TT);
  for (int i = 0; i<N_TT; i++){
    double temp = TTcells[i];
    TTcells[i] = temp/power1;
    bNumbers[i] = temp - TTcells[i]*power1;
  }

  // Creates H and typeTwo Info. Loop through all type 2 cell instances.
  int BN[M];
  bool inThread, mustSend;
  int i = 0;
  while (i<N_TT){
    int cell = TTcells[i];
    inThread = false;
    mustSend = false;
    BN[0] = bNumbers[i];
    if ( minBasis<=BN[0] && BN[0]<maxBasis ) inThread = true;
    else mustSend = true;

    if (inThread) typeTwoInfo[cell]+=1;

    int k = 1;
    while(cell==TTcells[i+k]){
      BN[k] = bNumbers[i+k];
      if ( minBasis<=BN[k] && BN[k]<maxBasis ){
        inThread = true;
        typeTwoInfo[cell]+=1;
      }
      else mustSend = true;
      k++;
    }

    if (inThread && mustSend){
      for (int j = 0; j<k; j++){
        if ( BN[j]<minBasis || maxBasis<=BN[j] ) H[cell].insert(BN[j]);
      }
    }
    i = i + k;
  }


  if (RANK==2);
}


void connStrenght(Grid& grid){
  int N = grid.N;
  int M = grid.M;
  int n_basis = grid.n_basis;
  int* support = grid.support;
  int* celltypes = grid.celltypes;
  int* offsets = grid.offsets;

  int TTcells[n_basis];  // cell numbers of type 2 cells will be stored here;
  int bNumbers[n_basis];  // corresponding basis numbers will be stored here;
  int N_TT = 0;         // The total number of type 2 cells will be stored here.

  int power1 = 100;
  int power2 = 100;
  while (M>power1) power1*=10;
  while (N>power2) power2*=10;
  int base = power1*power2;

  int k = 0; // basis number counter;
  for (int i = 0; i<n_basis; i++){
    if (i==offsets[k+1]) k++;
    if (celltypes[i]==2){
      TTcells[N_TT] = support[i]*power1 +k;
      N_TT++;
    }
  }

  std::sort(TTcells, TTcells+N_TT);
  for (int i = 0; i<N_TT; i++){
    double temp = TTcells[i];
    TTcells[i] = temp/power1;
    bNumbers[i] = temp - TTcells[i]*power1;
  }

  // ConnStrenght Matrix
  int strength[M*M]; for (int i = 0; i<M*M; i++) strength[i] = 0;

  int i= 0; k = 0;

  while (i<N_TT){
    k = 1;
    while (TTcells[i+k]==TTcells[i]){k++; }
    if (k==2){
      int b1 = bNumbers[i];
      int b2 = bNumbers[i+1];
      strength[M*b1+b2]+=1;
      strength[M*b2+b1]+=1;
    }
    i+=k;
  }


  for (int i = 0; i<M*M; i++){
    cout<<strength[i]<<"\t";
    if ((i+1)%M==0)cout<<endl;
  }

}



/* ************************************************************************** */
void gatherBasis(Grid& grid, int* basisDistr, Basis* basisArray, int RANK,
                                          int SIZE, double* global_GSolution){
  int recvCounts [SIZE]; int displs [SIZE];
  int temp = 0; int k = 0;
  for (int i = 0; i < grid.M; i++){
    if (i==basisDistr[k+1]){
      recvCounts[k] = temp;
      temp = (grid.offsets[i+1]-grid.offsets[i]);
      k++;
    }
    else{ temp+=(grid.offsets[i+1]-grid.offsets[i]); }
  } recvCounts[SIZE-1] = temp;
  displs[0] = 0;
  for (int i =1; i<SIZE; i++){ displs[i] = displs[i-1] + recvCounts[i-1]; }
  MPI_Gatherv(basisArray[0].values, recvCounts[RANK], MPI_DOUBLE,
          global_GSolution, recvCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/* ************************************************************************** */
void computeDiscrepancy(Grid& grid, double* GSol,  double* ompSol  ){
  bool mError = false;
  double LINF = 0; double L2 = 0;
  double ompSum = 0;
  for (int i = 0; i<grid.n_basis; i++){
    double temp = abs(GSol[i] - ompSol[i]);
    L2+=temp*temp;
    ompSum+=ompSol[i]*ompSol[i];

    if (temp>LINF) LINF = temp;
    if (LINF>1 && !mError){
      cout<<"MAJOR ERROR"<<endl;
      mError = true;
    }
    if (temp!=temp && !mError){
      cout<<"NaN OCCURED!"<<endl;
      mError = true;
    }
  }

  L2 = sqrt(L2/ompSum);



  cout<<LINF<<endl;
}

/* ************************************************************************** */
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

/* ************************************************************************** */
int threadStart(int RANK, int SIZE, int M ){
  int basisEach = M / SIZE;
  int rest = M - basisEach * SIZE;
  if (RANK < rest ) return (basisEach+1)*RANK;
  else return (basisEach+1)*rest + basisEach*(RANK-rest);
}














void printArray(int* a, int l){ for (int i =0; i<l; i++) cout<<a[i]<<"\t"; cout<<endl; }
void printArray(double* a, int l){ for (int i =0; i<l; i++) cout<<a[i]<<"\t"; cout<<endl;}

/* ************************************************************************** */
void printMap(map<int, set<double*> >& myMap){
  cout<<"______________________________"<<endl<<"Printing Map"<<endl
  <<"SIZE:\t"<<myMap.size()<<endl;
  map<int, set<double*> >::iterator it;
  set<double*>::iterator it2;
  for (it = myMap.begin(); it!=myMap.end(); it++){
    cout<<it->first<<":\t";
    for (it2 = it->second.begin(); it2!= it->second.end(); it2++)cout<<*it2<<"\t";
    cout<<endl;
  } cout<<"______________________________"<<endl;
}

/* ************************************************************************** */
void printH(map<int, set<int> >& H){
  cout<<"______________________________"<<endl<<"Printing H"<<endl
  <<"SIZE:\t"<<H.size()<<endl;
  map<int, set<int> >::iterator it;
  set<int>::iterator it2;
  for (it = H.begin(); it!=H.end(); it++){
    cout<<it->first<<":\t";
    for (it2 = it->second.begin(); it2!= it->second.end(); it2++)cout<<*it2<<"\t";
    cout<<endl;
  } cout<<"______________________________"<<endl;
}

/* ************************************************************************** */
void printMap(map<int,int> mymap){
  for (map<int,int>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
    cout << it->first << ":\t" << it->second << '\n';
}

/* ************************************************************************** */
void printSet(set<int> mySet){
  cout<<"Printing set"<<endl;
  set<int>::iterator it;
  for (it = mySet.begin(); it!=mySet.end(); it++) cout<<*it<<"\t";
  cout<<endl;
}
