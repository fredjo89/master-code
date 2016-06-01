#include "Basis.h"
#include "data_defs.h"

#include <iomanip>
#include <map>
#include <iostream>
#include <math.h>
using namespace std;

/* ************************************************************************** */
Basis::Basis(Grid& grid, Matrix& mat, double* globalBasis, int n){
  omega = (double) 2/3;
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

Basis::Basis(Grid& grid, Matrix& mat, double* globalBasis, int n, double omg)
: Basis(grid, mat,globalBasis,n){omega = omg; }

/* ************************************************************************** */
Basis& Basis::operator=(Basis& rhs){
  number = rhs.number;
  size = rhs.size;
  maxRow = rhs.maxRow;
  omega = rhs.omega;

  values = rhs.values;
  support = rhs.support;
  celltypes = rhs.celltypes;

  swap(mat_indices, rhs.mat_indices);
  swap(mat_coef, rhs.mat_coef);
  swap(update, rhs.update);

  return *this;
}

/* ************************************************************************** */
void Basis::makeLocalMapping(){
  map<int,int> GtoL;
  for (int i =0; i<size; i++){ GtoL[support[i]] = i; }
  GtoL[-1]=-1;
  map<int,int>::iterator mapPtr;
  for (int i = 0; i<size*maxRow; i++){
    mapPtr = GtoL.find(mat_indices[i]);   // Probably costly
    if (mapPtr != GtoL.end()){
      mat_indices[i] = mapPtr->second;
    }
    else{
      mat_indices[i] = -2;
    }
  }
}

/* ************************************************************************** */
bool Basis::jacobiProduct(double tol){
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

  bool done = true;

  for (int i =0; i<size; i++){
    if (celltypes[i]==0){
      update[i] = omega*update[i];
      if ( done == true && fabs(update[i]) > tol){
        done = false;
      }
      values[i]-=update[i];
    }
  }


  return done;
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
