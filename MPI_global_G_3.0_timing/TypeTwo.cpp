
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


TypeTwo::TypeTwo(Basis* basisArray, int M, int S, int R, int* bD,
                                      map<int,set<int> >& h) : TypeTwo(){
  H = h;
  sendL = H.size();

  omega = basisArray[0].omega;

  SIZE = S;
  RANK = R;
  basisDistr = bD;
  map<int,int> typeTwoInfo;

  for (int i = 0; i<M; i++){
    int* support = basisArray[i].support;
    int* celltypes = basisArray[i].celltypes;
    for (int j = 0; j<basisArray[i].size; j++){
      if (celltypes[j]==2){
        typeTwoInfo[support[j]]+=1;
      }
    }
  }

  sumL = typeTwoInfo.size();
  aIndex = new int[sumL];
  twoN = new int[sumL];
  sum = new double [sumL]{};


  map<int,int>::iterator it;
  int k=0;
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

  for (int i = 0; i<M; i++){
    int* support = basisArray[i].support;
    int* celltypes = basisArray[i].celltypes;
    double* values = basisArray[i].values;
    double* update = basisArray[i].update;

    for (int j = 0; j<basisArray[i].size; j++){
      if (celltypes[j]==2){
        int cell = support[j];
        int counter = 0;
        for (counter = 0; counter<sumL; counter++){
          if (cell==twoN[counter]) break;
        }
        int index = aIndex[counter];
        while (valAddress[index]!=NULL) index++;
        valAddress[index] = &values[j];
        upAddress[index] = &update[j];
      }
    }
  }
}

// *******************setupSending()***************************************** //
void TypeTwo::setupSending(){

  // Finds all basis function that requires communication
  set<int> commBasis; // All basis functions which requires communication
  map<int, set<int> >::iterator it1;
  set<int>::iterator it2;
  for (it1=H.begin(); it1!=H.end(); it1++){
    for (it2=it1->second.begin(); it2!=it1->second.end(); it2++) commBasis.insert(*it2);
  }
  for (int i = basisDistr[RANK]; i<basisDistr[RANK+1]; i++) commBasis.erase(i);

  // Finds all threads that requires communication
  for (it2=commBasis.begin(); it2!=commBasis.end(); it2++){
    for (int i = 0; i<SIZE+1; i++){
      if (*it2>=basisDistr[i] && *it2<basisDistr[i+1]){
        sendRanks.insert(i);
      }
    }
  }

  // Preforms a sending for each thread to figure out how much data it will
  // recieve from each process.
  Ncomm = sendRanks.size();
  recvCounts = new int[Ncomm];

  MPI_Status statSend[Ncomm], statRecv[Ncomm];
  MPI_Request	send_request[Ncomm],recv_request[Ncomm];
  int k = 0;
  for (it2 = sendRanks.begin(); it2!=sendRanks.end(); it2++){
    MPI_Isend(&sendL, 1, MPI_INT, *it2, 100,MPI_COMM_WORLD, &send_request[k]);
    MPI_Irecv(&recvCounts[k], 1, MPI_INT, *it2, 100, MPI_COMM_WORLD, &recv_request[k]);
    k++;
  }
  MPI_Waitall(Ncomm, send_request, statSend);
  MPI_Waitall(Ncomm, recv_request, statRecv);


  // By now we have:
  // sendRanks: All threads that it must send and recieve to
  // recvCounts: The number of elements it will recieve from each threads

  for (int i = 0; i<sendRanks.size(); i++) buffL+=recvCounts[i];

  recvBuff = new double [buffL];

  double twoNdouble[sendL];
  for (int i =0; i<sendL; i++) twoNdouble[i] = twoN[i];

  displs = new int[sendRanks.size()]; displs[0] = 0;
  for (int i = 1; i<sendRanks.size(); i++) displs[i] = displs[i-1]+recvCounts[i-1];

  k = 0;
  for (it2 = sendRanks.begin(); it2!=sendRanks.end(); it2++){
    MPI_Isend(twoNdouble, sendL, MPI_DOUBLE, *it2, 100,MPI_COMM_WORLD, &send_request[k]);
    MPI_Irecv(&recvBuff[displs[k]], recvCounts[k], MPI_DOUBLE, *it2, 100, MPI_COMM_WORLD, &recv_request[k]);
    k++;
  }
  MPI_Waitall(Ncomm, send_request, statSend);
  MPI_Waitall(Ncomm, recv_request, statRecv);

  map<int,int> recvIndices;
  for (int i = 0; i <buffL; i++){
    int cellNumber = recvBuff[i];
    for (int j = 0; j<sumL; j++){
      if (twoN[j]==cellNumber){
        recvIndices[i] = j;
        // break??
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
}
// ************************************************************************** //

void TypeTwo::sendAndRecieve(){

  MPI_Status statSend[Ncomm], statRecv[Ncomm];
  MPI_Request	send_request[Ncomm],recv_request[Ncomm];
  set<int>::iterator it;
  int k = 0;
  for (it = sendRanks.begin(); it!=sendRanks.end(); it++){
    MPI_Isend(sum, sendL, MPI_DOUBLE, *it, 100,MPI_COMM_WORLD, &send_request[k]);
    MPI_Irecv(&recvBuff[displs[k]], recvCounts[k], MPI_DOUBLE, *it, 100, MPI_COMM_WORLD, &recv_request[k]);
    k++;
  }
  MPI_Waitall(Ncomm, send_request, statSend);
  MPI_Waitall(Ncomm, recv_request, statRecv);

  for (int i = 0; i<mapL; i++) sum[mapOne[i]]+=recvBuff[mapTwo[i]];
}

void TypeTwo::localSum(){
  double oneSum = 0;
  int k = 0;
  for (int i =0; i<addL; i++){
    if (i!=aIndex[k+1]) oneSum+=*upAddress[i];
    else{
      sum[k] = oneSum;
      oneSum = *upAddress[i];
      k++;
    }
  }
  sum[k] = oneSum;



}



void TypeTwo::TTupdate(){
  int k = 1;
  for (int i = 0; i<addL; i++){
    if (i==aIndex[k]) k++;
    double temp = omega*sum[k-1];
    *valAddress[i]-= (*upAddress[i]*omega-*valAddress[i]*temp) / (1+temp);
  }
}

void TypeTwo::print(){
  cout<<"================================================================="<<endl
  <<"PRINTING TypeTwo"<<endl;
  cout<<"RANK: "<<RANK<<endl<<"SIZE: "<<SIZE<<endl;
  cout<<"basisDistr: "<<endl;
  for (int i = 0; i<SIZE+1; i++) cout<<basisDistr[i]<<endl;
  cout<<"Number OF TYPE 2 CELLS: "<<sumL<<endl
  <<"TOTAL NUMBER OF OCCURENCESS: "<<addL<<endl;
  cout<<"PRINTING aInxed:"<<endl; for (int i = 0; i<sumL; i++){ cout<<aIndex[i]<<endl; }
  cout<<"PRINTING twoN:"<<endl; for (int i = 0; i<sumL; i++){ cout<<twoN[i]<<endl; }
  cout<<"PRINTING sum:"<<endl;
  for (int i = 0; i<sumL; i++){ cout<<sum[i]<<endl; }
  cout<<"PRINTING valAddress and upAddress:"<<endl;
  cout<<"cell\tval\t\tvalue\tup\t\t value"<<endl
  <<"================================================================="<<endl;
  int k =0;
  for (int i = 0; i<addL; i++){
    if (i==aIndex[k+1]) k++;
    cout<<twoN[k]<<":\t";
    cout<<valAddress[i]<<":\t"<<*valAddress[i]<<"\t";
    cout<<upAddress[i]<<":\t"<<*upAddress[i]<<endl;
  }
}
