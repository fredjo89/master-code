#include "Basis.h"
#include "data_defs.h"
#include "utilities.h"

#include <iomanip>
#include <map>
#include <iostream>
#include <math.h>
using namespace std;

Basis::Basis(){
  values = NULL; support = NULL; celltypes = NULL;
  update = NULL; mat_indices = NULL; mat_coef = NULL;
}
Basis::~Basis(){
  delete[] mat_indices;
  delete[] mat_coef;
  delete[] update;
}


// makeBasis***************************************************************** //
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

  update = new double[size];
  mat_indices = new int[size*maxRow];
  mat_coef = new double[size*maxRow];

  int currentIndex, k = 0;
  for (int i=0; i<size; i++){
    currentIndex = support[i]*maxRow;
    for (int j = 0; j<maxRow; j++){
      mat_indices[k] = mat.j_index[currentIndex+j];
      mat_coef[k]    = mat.conn[currentIndex+j];
      k++;
    }
  }

  // Make Local Mapping //
  int max = support[0]; int min = support[0];
  for (int i=0; i<size; i++){
    if (support[i]>max) max = support[i];
    else if (support[i]<min) min = support[i];
  }

  int* myMap = new int[max-min+1];
  for (int i=0; i<max-min+1; i++) myMap[i] = -2;
  for (int i=0; i<size; i++) myMap[support[i]-min] = i;

  int index;
  for (int i=0; i<size*maxRow; i++){
    index = mat_indices[i];
      if (index!=-1){
        if (index>=min && index<=max) mat_indices[i] = myMap[index - min];
        else mat_indices[i] = -2;
      }
  }
  delete[] myMap;
}


// jacobiProduct************************************************************* //
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
  for (int i=0; i<size; i++){
    if (celltypes[i]==0) values[i]-=omega*update[i];
  }
}

// jacobiProduct_lockG******************************************************* //
void Basis::jacobiProduct_lockG(){
  for (int i=0; i<size; i++){
    if(celltypes[i]==0){
      update[i] = values[i];
      int it = i*maxRow;
      for (int j =0; j<maxRow; j++){
        int index = mat_indices[it+j];
        if (index==-1) break;
        if (index!=-2){ update[i]+=mat_coef[it+j]*values[index]; }
      }
    }
  }
  for (int i=0; i<size; i++){
    if (celltypes[i]==0) values[i]-=omega*update[i];
  }
}

// jacobiProduct_ConvCheck*************************************************** //
void Basis::jacobiProduct_ConvCheck(){
  for (int i=0; i<size; i++){
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
  for (int i=0; i<size; i++){
    if (celltypes[i]==0) {
      values[i]-=omega*update[i];
      if (*underTol && fabs(update[i])>tolerance) *underTol = false;
    }
  }
}



// ************************************************************************** //
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
// ************************************************************************** //
void Basis::printPressure(){
  cout << setprecision(3) << fixed<<"PRINTING BASIS"<<endl<<"values:\t\t";
  for (int i = 0; i<size; i++)cout<<values[i]<<"\t";
  cout<<endl<<"Cells:\t\t";
  for (int i = 0; i<size; i++){ cout<<support[i]<<"\t"; }
  cout<<endl<<"celltypes:\t";
  for (int i = 0; i<size; i++){ cout<<celltypes[i]<<"\t"; } cout<<endl;
}
