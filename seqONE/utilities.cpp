#include "utilities.h"


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
