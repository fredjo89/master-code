#include "TypeTwo.h"

#include "Basis.h"
#include "data_defs.h"
#include "utilities.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <mpi.h>
#include <stdexcept>
#include <cmath>
using namespace std;


TypeTwo::~TypeTwo(){
  delete[] valAddress; delete[] upAddress,
  delete[] aIndex; delete[] twoN; delete[] sum;
  delete[] mapOne; delete[] mapTwo;
  delete[] recvBuff; delete[] recvCounts; delete[] displs;
  delete[] sendRanks;
  delete[] statSend; delete[] statRecv; delete[] send_request; delete[] recv_request;
}

// ************************************************************************** //
TypeTwo::TypeTwo(Grid& grid, Basis* basisArray, int M, int S, int R, int* bD) {
  addL=0; sumL=0; buffL=0; Ncomm = 0; mapL = 0;
  int  k;
  int index, cell;
  map<int,int>::iterator it;
  int* support;
  int* celltypes;
  double* values;
  double* update;

  omega = basisArray[0].omega;
  SIZE = S;
  RANK = R;
  basisDistr = bD;

  map<int,int> typeTwoInfo;
  makeH(grid,H, typeTwoInfo, RANK, basisDistr);

  sendL = H.size();
  sumL = typeTwoInfo.size();
  aIndex = new int[sumL];
  twoN = new int[sumL];
  sum = new double [sumL]{};

  k = 0;
  // First placing the cells that will be sent
  for (it = typeTwoInfo.begin(); it != typeTwoInfo.end(); it++){
    if ( H.find(it->first) != H.end() ) {
      twoN[k] = it->first;
      aIndex[k] = addL;
      addL+=it->second;
      k++;
    }
  }
  // Then placing the rest.
  for (it = typeTwoInfo.begin(); it != typeTwoInfo.end(); it++){
    if ( H.find(it->first) == H.end() ) {
      twoN[k] = it->first;
      aIndex[k] = addL;
      addL+=it->second;
      k++;
    }
  }

  valAddress = new double*[addL]{};
  upAddress = new double*[addL]{};

  for ( int i=0; i<M; i++){
    support = basisArray[i].support;
    celltypes = basisArray[i].celltypes;
    values = basisArray[i].values;
    update = basisArray[i].update;

    for ( int j=0; j<basisArray[i].size; j++){
      if (celltypes[j]==2){
        cell = support[j];
        k = 0;
        for (k=0; k<sumL; k++){
          if (cell==twoN[k]) break;
        }
        index = aIndex[k];
        while (valAddress[index]!=NULL) index++;
        valAddress[index] = &values[j];
        upAddress[index] = &update[j];
      }
    }
  }

}

// setupSending()************************************************************ //
void TypeTwo::setupSending(){
  int k, l, cellNumber;
  int* temp1, *temp2, *temp3;
  int M = basisDistr[SIZE];       // Total number of basis functions.
  map<int, set<int> >::iterator it1;
  set<int>::iterator it2;

  temp1 = new int[M]{};
  temp2 = new int[M]{};

  // Finds all threads that requires communication
  k = 0;
  for (it1=H.begin(); it1!=H.end(); it1++){
    for (it2=it1->second.begin(); it2!=it1->second.end(); it2++){
      if (temp1[*it2] ==0){ temp1[*it2] = 1; k++; }
    }
  }
  l = 0;
  for (int i=0; i<M; i++) {
    if (temp1[i]==1) {temp2[l] = i; l++;}
    temp1[i] = 0;
  }
  for (int j = 0; j<k; j++){
    for (int i = 0; i<SIZE+1; i++){
      if (temp2[j]>=basisDistr[i] && temp2[j]<basisDistr[i+1]){
        if (temp1[i]==0){ temp1[i] = 1; Ncomm++; }
      }
    }
  }
  sendRanks = new int[Ncomm]; k = 0;
  for ( int i=0; i<SIZE; i++){
    if (temp1[i]==1) { sendRanks[k] = i; k++; }
  }delete[] temp1; delete[] temp2;



  // Preforms a sending for each thread to figure out how much data it will
  // recieve from each process.
  recvCounts = new int[Ncomm];

  statSend = new MPI_Status[Ncomm]; statRecv = new MPI_Status[Ncomm];
  send_request = new MPI_Request[Ncomm]; recv_request = new MPI_Request[Ncomm];


  for ( int i=0; i<Ncomm; i++){
    MPI_Isend(&sendL, 1, MPI_INT, sendRanks[i], 100,MPI_COMM_WORLD, &send_request[i]);
    MPI_Irecv(&recvCounts[i], 1, MPI_INT, sendRanks[i], 100, MPI_COMM_WORLD, &recv_request[i]);
  }
  MPI_Waitall(Ncomm, send_request, statSend);
  MPI_Waitall(Ncomm, recv_request, statRecv);


  //if (RANK!=-1) cout<<RANK<<"\t"<<Ncomm<<endl;
  // By now we have:
  // sendRanks: All threads that it must send and recieve to
  // recvCounts: The number of elements it will recieve from each threads

  for ( int i=0; i<Ncomm; i++) buffL+=recvCounts[i];

  recvBuff = new double [buffL];
  temp1 = new int [buffL];
  temp2 = new int[buffL];
  temp3 = new int [buffL];
  displs = new int[Ncomm]; displs[0] = 0;
  for ( int i=1; i<Ncomm; i++) displs[i] = displs[i-1]+recvCounts[i-1];



  for ( int i=0; i<Ncomm; i++){
    MPI_Isend(twoN, sendL, MPI_INT, sendRanks[i], 100,MPI_COMM_WORLD, &send_request[i]);
    MPI_Irecv(&temp1[displs[i]], recvCounts[i], MPI_INT, sendRanks[i], 100, MPI_COMM_WORLD, &recv_request[i]);
  } MPI_Waitall(Ncomm, send_request, statSend); MPI_Waitall(Ncomm, recv_request, statRecv);

  for ( int i=0; i<buffL; i++){
    for ( int j=0; j<sumL; j++){
      if (twoN[j]==temp1[i]){
        temp2[mapL] = i;
        temp3[mapL] = j;
        mapL++;
        break;
      }
    }
  }

  mapOne = new int[mapL];
  mapTwo = new int[mapL];

  for (int i=0; i<mapL; i++){
    mapOne[i]=temp3[i];
    mapTwo[i]=temp2[i];
  }

  delete[] temp1; delete[] temp2; delete[] temp3;
  H.clear();
}

// ************************************************************************** //
void TypeTwo::sendAndRecieve(){
  for ( int i=0; i<Ncomm; i++){
    MPI_Isend(sum, sendL, MPI_DOUBLE, sendRanks[i], 100,MPI_COMM_WORLD, &send_request[i]);
    MPI_Irecv(&recvBuff[displs[i]], recvCounts[i], MPI_DOUBLE, sendRanks[i], 100, MPI_COMM_WORLD, &recv_request[i]);
  } MPI_Waitall(Ncomm, send_request, statSend); MPI_Waitall(Ncomm, recv_request, statRecv);

  for ( int i=0; i<mapL; i++) sum[mapOne[i]]+=recvBuff[mapTwo[i]];
}

// ************************************************************************** //
void TypeTwo::localSum(){
  double oneSum = 0; int k = 0;
  for ( int i=0; i<addL; i++){
    if (i!=aIndex[k+1]) oneSum+=*upAddress[i];
    else{
      sum[k] = oneSum;
      oneSum = *upAddress[i];
      k++;
    }
  }
  sum[k] = oneSum;
}

// ************************************************************************** //
void TypeTwo::TTupdate(){
  int k = 1;
  for ( int i=0; i<addL; i++){
    if (i==aIndex[k]) k++;
    double temp = omega*sum[k-1];
    *valAddress[i]-= (*upAddress[i]*omega-*valAddress[i]*temp) / (1+temp);
  }
}


// ************************************************************************** //
void TypeTwo::print(){
  int i, j, k, l;
  cout<<"================================================================="<<endl
  <<"PRINTING TypeTwo"<<endl<<"RANK: "<<RANK<<endl<<"SIZE: "<<SIZE<<endl<<"basisDistr: "<<endl;
  for ( i = 0; i<SIZE+1; i++) cout<<basisDistr[i]<<endl
  <<"Number OF TYPE 2 CELLS: "<<sumL<<endl<<"TOTAL NUMBER OF OCCURENCESS: "<<addL<<endl
  <<"PRINTING aInxed:"<<endl; for ( i = 0; i<sumL; i++){ cout<<aIndex[i]<<endl; }
  cout<<"PRINTING twoN:"<<endl; for (  i = 0; i<sumL; i++){ cout<<twoN[i]<<endl; }
  cout<<"PRINTING sum:"<<endl;
  for ( i = 0; i<sumL; i++){ cout<<sum[i]<<endl; }
  cout<<"PRINTING valAddress and upAddress:"<<endl
  <<"cell\tval\t\tvalue\tup\t\t value"<<endl
  <<"================================================================="<<endl;
  k =0;
  for ( i = 0; i<addL; i++){
    if (i==aIndex[k+1]) k++;
    cout<<twoN[k]<<":\t"<<valAddress[i]<<":\t"<<*valAddress[i]<<"\t"
    <<upAddress[i]<<":\t"<<*upAddress[i]<<endl;
  }
}
