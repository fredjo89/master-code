
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

/* ************************************************************************** */
Basis::Basis(Grid& grid, Matrix& mat, double* globalBasis, int n){
  number = n;
  size = grid.offsets[n+1]-grid.offsets[n];
  maxRow = mat.maxRow;

  support = &grid.support[grid.offsets[n]];
  celltypes = &grid.celltypes[grid.offsets[n]];
  values = &globalBasis[grid.offsets[n]];

  update = new double [size];
  mat_indices = new int[size * maxRow];
  mat_coef = new double[size*maxRow];

  int* it1 = mat_indices;
  double* it2 = mat_coef;
  for (int i =0; i<size; i++){
    int currentIndex = support[i]*maxRow;
    for (int j =0; j<maxRow; j++){
      *it1 = mat.j_index[currentIndex+j];
      it1++;
      *it2 = mat.conn[currentIndex+j];
      it2++;
    }
  }
}

/* ************************************************************************** */
Basis& Basis::operator=(Basis& rhs){
  number = rhs.number;
  size = rhs.size;
  maxRow = rhs.maxRow;

  swap(values,rhs.values);
  swap(support,rhs.support);
  swap(celltypes, rhs.celltypes);
  swap(mat_indices, rhs.mat_indices);
  swap(mat_coef, rhs.mat_coef);
  swap(update, rhs.update);

  return *this;
}

/* ************************************************************************** */
void Basis::makeLocalMapping(){
  map<int,int> GtoL;
  for (int i =0; i<size; i++){
    GtoL[support[i]] = i;
  }
  GtoL[-1]=-1;

  map<int,int>::iterator mapPtr;
  for (int i = 0; i<size*maxRow; i++){
    mapPtr = GtoL.find(mat_indices[i]);
    if (mapPtr != GtoL.end()){
      mat_indices[i] = mapPtr->second;
    }
    else{
      mat_indices[i] = -2;
    }
  }
}

/* ************************************************************************** */
void Basis::jacobiProduct(){
  double omega = (double) 2/3;
  for (int i = 0; i<size; i++){
    update[i] = values[i];
    if(celltypes[i]!=1){
      int it = i*maxRow;
      for (int j =0; j<maxRow; j++){
        int index = mat_indices[it+j];
        if (index==-1) break;
        if (index!=-2) update[i]+=mat_coef[it+j]*values[index];
      }
    }
  }
  for (int i =0; i<size; i++){
    if (celltypes[i]==0) values[i]-=omega*update[i];
  }
}



/* ************************************************************************** */
void Basis::print(){
  if (values==NULL){cout<<"INVALID BASIS"<<endl; return;}
  cout << setprecision(2) << fixed
  <<"PRINTING BASIS"<<endl<<"Number:\t\t"<<number<<endl
  <<"Size:\t\t"<<size<<endl<<"maxRow:\t\t"<<maxRow<<endl;
  cout<<"values:\t\t";
  for (int i = 0; i<size; i++){
    cout<<values[i]<<"\t";
  }
  cout<<endl<<"update:\t\t";
  for (int i = 0; i<size; i++){
    cout<<update[i]<<"\t";
  }
  cout<<endl<<"Cells:\t\t";
  for (int i = 0; i<size; i++){
    cout<<support[i]<<"\t";
  }
  cout<<endl<<"celltypes:\t";
  for (int i = 0; i<size; i++){
    cout<<celltypes[i]<<"\t";
  }
  cout<<endl<<endl<<"mat_indices:\t";
  for (int i = 0; i<size*maxRow; i++){
    cout<<mat_indices[i]<<"\t";
    if ((i+1)%maxRow==0) cout<<endl<<"\t\t";
  }
  cout<<endl<<"mat_coef:\t";
  for (int i = 0; i<size*maxRow; i++){
    cout<<mat_coef[i]<<"\t";
    if ((i+1)%maxRow==0) cout<<endl<<"\t\t";
  }
  cout<<endl;
}
/* ************************************************************************** */
void Basis::printPressure(){
  cout << setprecision(3) << fixed<<"PRINTING BASIS"<<endl<<"values:\t\t";
  for (int i = 0; i<size; i++){
    cout<<values[i]<<"\t";
  }
  cout<<endl<<"Cells:\t\t";
  for (int i = 0; i<size; i++){ cout<<support[i]<<"\t"; }
  cout<<endl<<"celltypes:\t";
  for (int i = 0; i<size; i++){ cout<<celltypes[i]<<"\t"; }
  cout<<endl;
}
