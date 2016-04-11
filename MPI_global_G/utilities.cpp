#include "utilities.h"

#include <iostream>
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



int threadStart(int RANK, int SIZE, int M ){
  int basisEach = M / SIZE;

  int rest = M - basisEach * SIZE;

  if (RANK < rest ){
    return (basisEach+1)*RANK;
  }
  else{
    return (basisEach+1)*rest + basisEach*(RANK-rest);
  }

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
