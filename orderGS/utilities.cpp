#include "utilities.h"
#include "data_defs.h"
//#include "data_defs.cpp"  // Why do i need this??
#include <iostream>
#include <set>
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
