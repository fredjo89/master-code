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

TypeTwo::TypeTwo(){
  addL=0; sumL=0; buffL=0; Ncomm = 0;
  valAddress=NULL; upAddress=NULL;
  aIndex=NULL; twoN=NULL; sum=NULL;
  mapOne=NULL; mapTwo=NULL;
  recvBuff=NULL; recvCounts=NULL; displs=NULL;
  sendRanks = NULL;
}

TypeTwo::~TypeTwo(){
  delete[] valAddress; delete[] upAddress,
  delete[] aIndex; delete[] twoN; delete[] sum;
  delete[] mapOne; delete[] mapTwo;
  delete[] recvBuff; delete[] recvCounts; delete[] displs;
  delete[] sendRanks;
}

// ************************************************************************** //
TypeTwo::TypeTwo(Grid& grid, Basis* basisArray, int M, int S, int R, int* bD) : TypeTwo(){
  int i, j, k, l;
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
  for (it = typeTwoInfo.begin(); it!=typeTwoInfo.end(); it++){
    if ( H.find(it->first) != H.end() ) {
      twoN[k] = it->first;
      aIndex[k] = addL;
      addL+=it->second;
      k++;
    }
  }
  // Then placing the rest.
  for (it = typeTwoInfo.begin(); it!=typeTwoInfo.end(); it++){
    if ( H.find(it->first) == H.end() ) {
      twoN[k] = it->first;
      aIndex[k] = addL;
      addL+=it->second;
      k++;
    }
  }

  valAddress = new double*[addL]{};
  upAddress = new double*[addL]{};

  for ( i = 0; i<M; i++){
    support = basisArray[i].support;
    celltypes = basisArray[i].celltypes;
    values = basisArray[i].values;
    update = basisArray[i].update;

    for ( j = 0; j<basisArray[i].size; j++){
      if (celltypes[j]==2){
        cell = support[j];
        k = 0;
        for (k = 0; k<sumL; k++){
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

// *******************setupSending()***************************************** //
void TypeTwo::setupSending(){
  int i, j, k, l;
  int cellNumber;

  // Finds all basis function that requires communication
  set<int> commBasis; // All basis functions which requires communication

  map<int, set<int> >::iterator it1;
  set<int>::iterator it2;
  for (it1=H.begin(); it1!=H.end(); it1++){
    for (it2=it1->second.begin(); it2!=it1->second.end(); it2++) commBasis.insert(*it2);
  }

  int temp[SIZE]; for ( i = 0; i<SIZE; i++) temp[i] = 0;

  // Finds all threads that requires communication
  for (it2=commBasis.begin(); it2!=commBasis.end(); it2++){
    for ( i = 0; i<SIZE+1; i++){
      if (*it2>=basisDistr[i] && *it2<basisDistr[i+1]){
        if (temp[i]==0){
          temp[i] = 1;
          Ncomm++;
        }
      }
    }
  }

  sendRanks = new int[Ncomm];
  k = 0;
  for ( i = 0; i<SIZE; i++){
    if (temp[i]==1) {
      sendRanks[k] = i;
      k++;
    }
  }

  // Preforms a sending for each thread to figure out how much data it will
  // recieve from each process.
  recvCounts = new int[Ncomm];

  MPI_Status statSend[Ncomm], statRecv[Ncomm];
  MPI_Request	send_request[Ncomm],recv_request[Ncomm];
  for ( i = 0; i<Ncomm; i++){
    MPI_Isend(&sendL, 1, MPI_INT, sendRanks[i], 100,MPI_COMM_WORLD, &send_request[i]);
    MPI_Irecv(&recvCounts[i], 1, MPI_INT, sendRanks[i], 100, MPI_COMM_WORLD, &recv_request[i]);
  }
  MPI_Waitall(Ncomm, send_request, statSend);
  MPI_Waitall(Ncomm, recv_request, statRecv);

  // By now we have:
  // sendRanks: All threads that it must send and recieve to
  // recvCounts: The number of elements it will recieve from each threads

  for ( i = 0; i<Ncomm; i++) buffL+=recvCounts[i];

  recvBuff = new double [buffL];

  double twoNdouble[sendL];
  for ( i = 0; i<sendL; i++) twoNdouble[i] = twoN[i];

  displs = new int[Ncomm]; displs[0] = 0;
  for ( i = 1; i<Ncomm; i++) displs[i] = displs[i-1]+recvCounts[i-1];


  for ( i = 0; i<Ncomm; i++){
    MPI_Isend(twoNdouble, sendL, MPI_DOUBLE, sendRanks[i], 100,MPI_COMM_WORLD, &send_request[i]);
    MPI_Irecv(&recvBuff[displs[i]], recvCounts[i], MPI_DOUBLE, sendRanks[i], 100, MPI_COMM_WORLD, &recv_request[i]);
  }
  MPI_Waitall(Ncomm, send_request, statSend);
  MPI_Waitall(Ncomm, recv_request, statRecv);

  map<int,int> recvIndices;
  for ( i = 0; i <buffL; i++){
    cellNumber = recvBuff[i];
    for ( j = 0; j<sumL; j++){
      if (twoN[j]==cellNumber){
        recvIndices[i] = j;
        break;
      }
    }
  }

  mapL = recvIndices.size();
  mapOne = new int[mapL];
  mapTwo = new int[mapL];

  map<int,int>::iterator it3;
  k = 0;
  for (it3 = recvIndices.begin(); it3!=recvIndices.end(); it3++){
    mapOne[k] = it3->second;
    mapTwo[k] = it3->first;
    k++;
  }

  H.clear();
}

// ************************************************************************** //
void TypeTwo::sendAndRecieve(){
  int i;
  MPI_Status statSend[Ncomm], statRecv[Ncomm];
  MPI_Request	send_request[Ncomm],recv_request[Ncomm];
  for ( i = 0; i<Ncomm; i++){
    MPI_Isend(sum, sendL, MPI_DOUBLE, sendRanks[i], 100,MPI_COMM_WORLD, &send_request[i]);
    MPI_Irecv(&recvBuff[displs[i]], recvCounts[i], MPI_DOUBLE, sendRanks[i], 100, MPI_COMM_WORLD, &recv_request[i]);
  }
  MPI_Waitall(Ncomm, send_request, statSend);
  MPI_Waitall(Ncomm, recv_request, statRecv);

  for ( i = 0; i<mapL; i++) sum[mapOne[i]]+=recvBuff[mapTwo[i]];
}

// ************************************************************************** //
void TypeTwo::localSum(){
  int i, k;
  double oneSum = 0;
  k = 0;
  for ( i = 0; i<addL; i++){
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
  int i, k;
  k = 1;
  for ( i = 0; i<addL; i++){
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
