
#include "data_defs.h"
#include <iostream>
#include <iomanip>
#include <map>
using namespace std;


void Grid::print(){
  cout<<"PRINTING GRID"<<endl<<"Number of cells: "<<N<<endl
  <<"Number of blocks: "<<M<<endl<<"Number of elements in support: "
  <<n_basis<<endl<<"Offsets:"<<endl;
  for (int i = 0; i<M+1; i++){ cout<<offsets[i]<<"  "; } cout<<endl;

  cout<<"Centers:"<<endl;
  for (int i = 0; i<M; i++){ cout<<centers[i]<<"  "; } cout<<endl;

  cout<<"Support:"<<endl;
  int it = 1;
  for (int i = 0; i<n_basis; i++){
    if (i==offsets[it]){ cout<<endl<<endl; it++; }
    cout<<support[i]<<"  ";
  } cout<<endl;
  cout<<"Celltypes:"<<endl;
  it = 1;
  for (int i = 0; i<n_basis; i++){
    if (i==offsets[it]){ cout<<endl<<endl; it++; }
    cout<<celltypes[i]<<"  ";
  } cout<<endl;
}

void Matrix::print(){
    cout.precision(2);
    cout<<"PRINTING MATRIX"<<endl<<"Number of cells " <<N<<endl
    <<"MaxRow: "<<maxRow<<endl
    <<"j_index:"<<endl;
    for (int i = 0; i<maxRow*N; i++){
      cout<<j_index[i]<<"\t";
      if ((i+1)%maxRow==0){cout<<endl;}
    }
    cout<<"Coefficients:"<<endl;
    for (int i = 0; i<maxRow*N; i++){
      cout<<conn[i]<<"\t";
      if ((i+1)%maxRow==0){cout<<endl;}
    }
    if (loc_index!=NULL && loc_conn!=NULL){
      cout<<"loc_index:"<<endl;
      for (int i = 0; i<maxRow*n_basis; i++){
        cout<<loc_index[i]<<"\t";
        if ((i+1)%maxRow==0){cout<<endl;}
      }
      cout<<"loc_conn:"<<endl;
      for (int i = 0; i<maxRow*n_basis; i++){
        cout<<loc_conn[i]<<"\t";
        if ((i+1)%maxRow==0){cout<<endl;}
      }
    }
}

void Options::print(){
  cout<<"PRINTING OPTIONS"<<endl<<"tolerance: "<<tolerance<<endl
  <<"maxIter: "<<maxIter<<endl<<"omega: "<<omega<<endl
  <<"reNorm: "<<reNorm<<endl<<"checkTol: "<<checkTol<<endl;
}
