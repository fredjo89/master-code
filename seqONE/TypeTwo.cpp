
#include "TypeTwo.h"
#include "data_defs.h"
#include <iostream>
#include <iomanip>
#include <map>
using namespace std;


TypeTwo::TypeTwo(Basis* basisArray, int M) : TypeTwo(){
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
  int k =0;
  for (it = typeTwoInfo.begin(); it!=typeTwoInfo.end(); it++){
    twoN[k] = it->first;
    aIndex[k] = addL;
    addL+=it->second;
    k++;
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


void TypeTwo::sumAndModify(){
  double omega = (double) 2/3;
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

  k =1;
  for (int i =0; i<addL; i++){
    if (i==aIndex[k]) k++;
    double temp = omega*sum[k-1];
    *valAddress[i]-= (*upAddress[i]*omega-*valAddress[i]*temp) / (1+temp);
  }
}




void TypeTwo::print(){
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
