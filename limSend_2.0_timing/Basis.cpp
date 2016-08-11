#include "Basis.h"
#include "data_defs.h"

#include <iomanip>
#include <map>
#include <iostream>
#include <math.h>
using namespace std;

/* ************************************************************************** */
Basis::Basis(Grid& grid, Matrix& mat, double* globalBasis, int n, Options& opt){
  omega = opt.omega;
  tolerance = opt.tolerance;

  number = n;

  int start = grid.offsets[n];
  int end = grid.offsets[n+1];

  size = end - start;
  maxRow = mat.maxRow;

  support = &grid.support[start];
  celltypes = &grid.celltypes[start];
  values = &globalBasis[start];           // Optimize memory usage here?

  update = new double [size];
  mat_indices = new int[size * maxRow];
  mat_coef = new double[size*maxRow];

  int k = 0;
  for (int i =0; i<size; i++){
    int currentIndex = support[i]*maxRow;
    for (int j =0; j<maxRow; j++){
      mat_indices[k] = mat.j_index[currentIndex+j];
      mat_coef[k]    = mat.conn[currentIndex+j];
      k++;
    }
  }
}

/* ************************************************************************** */
Basis& Basis::operator=(Basis& rhs){
  omega = rhs.omega; number = rhs.number; size = rhs.size; maxRow = rhs.maxRow;
  values = rhs.values; support = rhs.support; celltypes = rhs.celltypes;
  swap(mat_indices, rhs.mat_indices);
  swap(mat_coef, rhs.mat_coef);
  swap(update, rhs.update);
  return *this;
}

/* ************************************************************************** */
void Basis::makeBasis(Grid& grid, Matrix& mat, double* initB,int n, Options& opt, bool& under){
  omega = opt.omega;
  tolerance = opt.tolerance;
  underTol = &under;

  number = n;

  int start = grid.offsets[n];
  int end = grid.offsets[n+1];

  size = end - start;
  maxRow = mat.maxRow;

  support = &grid.support[start];
  celltypes = &grid.celltypes[start];
  values = &initB[start];           // Optimize memory usage here?

  update = new double [size];
  mat_indices = new int[size * maxRow];
  mat_coef = new double[size*maxRow];

  int currentIndex, k = 0;
  for (int i =0; i<size; i++){
    currentIndex = support[i]*maxRow;
    for (int j = 0; j<maxRow; j++){
      mat_indices[k] = mat.j_index[currentIndex+j];
      mat_coef[k]    = mat.conn[currentIndex+j];
      k++;
    }
  }
}

/* ************************************************************************** */
void Basis::makeLocalMapping(){
  int max = support[0]; int min = support[0];
  for (int i = 0; i<size; i++){
    if (support[i]>max) max = support[i];
    else if (support[i]<min) min = support[i];
  }

  int myMap[max - min + 1];
  for (int i = 0; i<max-min+1; i++) myMap[i] = -2;
  for (int i = 0; i<size; i++) myMap[support[i]-min] = i;

  int index;
  for (int i = 0; i<size*maxRow; i++){
    index = mat_indices[i];
      if (index!=-1){
        if (index>=min && index<=max) mat_indices[i] = myMap[index - min];
        else mat_indices[i] = -2;
      }
  }
}





/* ************************************************************************** */
void Basis::jacobiProduct(){
  for (int i = 0; i<size; i++){
    update[i] = values[i];
    if(celltypes[i]!=1){
      int it = i*maxRow;
      for (int j =0; j<maxRow; j++){
        int index = mat_indices[it+j];
        if (index==-1) break;
        if (index!=-2){ update[i]+=mat_coef[it+j]*values[index]; }
      }
    }
  }

  for (int i =0; i<size; i++){
    if (celltypes[i]==0) values[i]-=omega*update[i];
  }
}

/* ************************************************************************** */
void Basis::jacobiProduct_ConvCheck(){
  for (int i = 0; i<size; i++){
    update[i] = values[i];
    if(celltypes[i]!=1){
      int it = i*maxRow;
      for (int j =0; j<maxRow; j++){
        int index = mat_indices[it+j];
        if (index==-1) break;
        if (index!=-2){ update[i]+=mat_coef[it+j]*values[index]; }
      }
    }
  }

  for (int i =0; i<size; i++){
    if (celltypes[i]==0) {
      values[i]-=omega*update[i];
      if (*underTol && fabs(update[i])>tolerance) *underTol = false;
    }
  }
}



/* ************************************************************************** */
void Basis::print(){
  if (values==NULL){cout<<"INVALID BASIS"<<endl; return;}
  cout << setprecision(2) << fixed
  <<"PRINTING BASIS"<<endl<<"Number:\t\t"<<number<<endl
  <<"Size:\t\t"<<size<<endl<<"maxRow:\t\t"<<maxRow<<endl;
  cout<<"values:\t\t";
  for (int i = 0; i<size; i++) cout<<values[i]<<"\t";
  cout<<endl<<"update:\t\t";
  for (int i = 0; i<size; i++) cout<<update[i]<<"\t";
  cout<<endl<<"Cells:\t\t";
  for (int i = 0; i<size; i++) cout<<support[i]<<"\t";
  cout<<endl<<"celltypes:\t";
  for (int i = 0; i<size; i++) cout<<celltypes[i]<<"\t";
  cout<<endl<<endl<<"mat_indices:\t";
  for (int i = 0; i<size*maxRow; i++){
    cout<<mat_indices[i]<<"\t";
    if ((i+1)%maxRow==0) cout<<endl<<"\t\t";
  }
  cout<<endl<<"mat_coef:\t";
  for (int i = 0; i<size*maxRow; i++){
    cout<<mat_coef[i]<<"\t";
    if ((i+1)%maxRow==0) cout<<endl<<"\t\t";
  } cout<<endl;
}
/* ************************************************************************** */
void Basis::printPressure(){
  cout << setprecision(3) << fixed<<"PRINTING BASIS"<<endl<<"values:\t\t";
  for (int i = 0; i<size; i++)cout<<values[i]<<"\t";
  cout<<endl<<"Cells:\t\t";
  for (int i = 0; i<size; i++){ cout<<support[i]<<"\t"; }
  cout<<endl<<"celltypes:\t";
  for (int i = 0; i<size; i++){ cout<<celltypes[i]<<"\t"; } cout<<endl;
}
