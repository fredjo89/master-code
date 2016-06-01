
#include "TypeTwo.h"
#include "Basis.h"
#include "data_defs.h"
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <mpi.h>
#include <stdexcept>
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


void TypeTwo::localSum(){
  double oneSum = 0;
  int k = 0;
  for (int i =0; i<addL; i++){
    if (i!=aIndex[k+1]){
      oneSum+=*upAddress[i];
    }
    else{
      sum[k] = oneSum;
      k++;
      oneSum = *upAddress[i];
    }
  }
  sum[k] = oneSum;
}


void TypeTwo::setupSending(){
  recvCounts = new int[SIZE];

  // replace sumL with sendL
  MPI_Allgather(&sendL, 1, MPI_INT, recvCounts, 1, MPI_INT,  MPI_COMM_WORLD);

  for (int i =0; i<SIZE; i++) buffL+=recvCounts[i];

  recvBuff = new double [buffL];

  // relace sumL with sendL
  double twoNdouble[sendL];
  for (int i =0; i<sendL; i++)twoNdouble[i] = twoN[i];

  displs = new int[SIZE]{};
  for (int i = 1; i<SIZE; i++){
    displs[i]= displs[i-1]+recvCounts[i-1];
  }

  // replace sumL with sendL
  MPI_Allgatherv(twoNdouble, sendL, MPI_DOUBLE,recvBuff, recvCounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

  map<int,int> recvIndices;
  for (int i = 0; i <buffL; i++){
    if (i<displs[RANK] || ( RANK!=(SIZE-1) && displs[RANK+1]<=i ) ){
      int cellNumber = recvBuff[i];
      for (int j = 0; j<sumL; j++){
        if (twoN[j]==cellNumber){
          recvIndices[i] = j;
        }
      }
    }
  }

  mapL = recvIndices.size();
  mapOne = new int[mapL];
  mapTwo = new int[mapL];

  map<int,int>::iterator it;
  int k = 0;
  for (it = recvIndices.begin(); it!=recvIndices.end(); it++){
    mapOne[k] = it->second;
    mapTwo[k] = it->first;
    k++;
  }

}

void TypeTwo::sendAndRecieve(){
  // replace sumL with sendL
  MPI_Allgatherv(sum, sendL, MPI_DOUBLE,recvBuff, recvCounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

  for (int i = 0; i<mapL; i++) sum[mapOne[i]]+=recvBuff[mapTwo[i]];
}

void TypeTwo::TTupdate(){
  int k =1;
  for (int i =0; i<addL; i++){
    if (i==aIndex[k]) k++;
    double temp = omega*sum[k-1];
    *valAddress[i]-= (*upAddress[i]*omega-*valAddress[i]*temp) / (1+temp);
  }
}

void TypeTwo::print(){
  cout<<"PRINTING TypeTwo"<<endl
  <<"================================================================="<<endl;
  cout<<"RANK: "<<RANK<<endl<<"SIZE: "<<SIZE<<endl;
  cout<<"basisDistr: "<<endl;
  for (int i = 0; i<SIZE+1; i++) cout<<basisDistr[i]<<endl;
  cout<<"Number OF TYPE 2 CELLS: "<<sumL<<endl
  <<"TOTAL NUMBER OF OCCURENCESS: "<<addL<<endl;
  cout<<"PRINTING twoN:"<<endl;
  for (int i = 0; i<sumL; i++){ cout<<twoN[i]<<endl; }
  cout<<"PRINTING aInxed:"<<endl;
  for (int i = 0; i<sumL; i++){ cout<<aIndex[i]<<endl; }
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
