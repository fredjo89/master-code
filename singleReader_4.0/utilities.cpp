#include "utilities.h"
#include "data_defs.h"

#include <mpi.h>
#include <cmath>
using namespace std;

// ********************************************************************************************** //
void setupSupportMatrix(Grid& grid, Matrix& mat){
    mat.loc_index = new int[mat.maxRow * mat.n_basis];
    mat.loc_conn = new double[mat.maxRow * mat.n_basis];
    int i = 0;  int j = 0;  int k = grid.support[0]*mat.maxRow;
    for (int it = 0; it < mat.n_basis*mat.maxRow; it++){
        if (j==mat.maxRow){
            j = 0; i++;
            k = grid.support[i]*mat.maxRow;
        }
        mat.loc_index[it] = mat.j_index[k + j];
        mat.loc_conn[it] = mat.conn[k + j];
        j++;
    }
}

// ********************************************************************************************** //
void makeMapping(Grid& grid, Matrix& mat){
    int max, min, index, start, end;
    int* myMap;

    for (int i = 0; i<grid.M; i++){
        start = grid.offsets[i];    end = grid.offsets[i+1];
        max = grid.support[start];  min = grid.support[start];

        for (int j = start; j<end; j++ ){
            if (grid.support[j]>max)        max = grid.support[j];
            else if (grid.support[j]<min)   min = grid.support[j];
        }

        myMap = new int[max-min+1];
        for (int j=0; j<max-min+1; j++)     myMap[j] = -2;
        for (int j=start; j<end; j++)       myMap[grid.support[j]-min] = j;

        for (int j=start*mat.maxRow; j<end*mat.maxRow; j++){
            index = mat.loc_index[j];
            if (index!=-1){
                if (index>=min && index<=max)   mat.loc_index[j] = myMap[index - min];
                else                            mat.loc_index[j] = -2;
            }
        }
        delete[] myMap;
    }
}

// ********************************************************************************************** //
void gatherBasis(Grid& localGrid, Grid& globalGrid, Graph& graph, int RANK, int SIZE){

    if(RANK==0){
        MPI_Status status;
        int k = graph.n_basis_Nodes[0];
        for (int i = 1; i<SIZE; i++){
            int offsets = globalGrid.offsets[k];
            MPI_Recv(&globalGrid.basis[offsets], graph.nodeWeights[i], MPI_DOUBLE, i, 666, MPI_COMM_WORLD, &status);
            k += graph.n_basis_Nodes[i];
        }
    }
    else{
        MPI_Send(localGrid.basis, localGrid.n_basis, MPI_DOUBLE, 0, 666, MPI_COMM_WORLD);
    }

}

// ********************************************************************************************** //
void computeDiscrepancy(Grid& grid, double* ompSol){
  double error = 0; bool mError = false;
  for (int i = 0; i<grid.n_basis; i++){
    double temp = abs(grid.basis[i] - ompSol[i]);
    if (temp>error) error = temp;
    if (error>1 && !mError){
      cout<<"MAJOR ERROR"<<endl;
      mError = true;
    }
    if (temp!=temp && !mError){
      cout<<"NaN OCCURED!"<<endl;
      cout<<grid.basis[i]<<"\t"<<ompSol[i]<<endl;
      mError = true;
    }
  } cout<<"Discrepancy:\t"<<error<<endl;
}

// ********************************************************************************************** //
void destroyMatrix(Matrix& mat){
    mat.N = -1;         mat.maxRow = -1;        mat.n_basis = -1;
    delete[] mat.conn;  delete[] mat.j_index;   delete[] mat.loc_index;     delete[] mat.loc_conn;
    mat.conn = NULL;    mat.j_index = NULL;     mat.loc_index = NULL;       mat.loc_conn = NULL;
}

// ********************************************************************************************** //
void printDistribution(Graph& graph, Node& node, int RANK, int SIZE){
    MPI_Barrier(MPI_COMM_WORLD);
    if (RANK==0) graph.print();
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i<SIZE; i++){
        if (RANK==i) node.print();
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// ********************************************************************************************** //
void reOrderGrid(Grid& grid, int* basisDistr, int SIZE){
    if (grid.order!=NULL){ cout<<"Invalid use of reOrder"<<endl; return; }

    int N = grid.N;
    int M = grid.M;
    int n_basis = grid.n_basis;
    int* order = grid.order;
    int* support = grid.support;
    int* celltypes = grid.celltypes;
    double* basis = grid.basis;
    int* offsets = grid.offsets;

    int* newSupport = new int[n_basis];
    int* newCelltypes = new int[n_basis];
    double* newBasis = new double[n_basis];
    int* newOffsets = new int[M+1]; newOffsets[0] = 0;

    order = new int[M];
    int k = 0;
    for (int i = 0; i<SIZE; i++){
        for (int j = 0; j<M; j++){
            if (i == basisDistr[j]){
                order[k] = j;
                k++;
            }
        }
    }

    for (int i=0; i<M; i++){
        int oldStart = offsets[order[i]];
        int oldEnd = offsets[order[i]+1];
        newOffsets[i+1] = newOffsets[i] + oldEnd - oldStart;
        k = 0;
        for (int j = newOffsets[i]; j<newOffsets[i+1]; j++){
            newSupport[j] = support[oldStart + k];
            newCelltypes[j] = celltypes[oldStart + k];
            newBasis[j] = basis[oldStart + k];
            k++;
        }
    }

    grid.order = order;

    delete[] support;   grid.support = newSupport;
    delete[] celltypes; grid.celltypes = newCelltypes;
    delete[] basis;     grid.basis = newBasis;
    delete[] offsets;   grid.offsets = newOffsets;
}

// ********************************************************************************************** //
void reOrderGrid(Grid& grid){
    if (grid.order==NULL){ cout<<"Invalid use of reOrder()"<<endl; return; }

    int N = grid.N;
    int M = grid.M;
    int n_basis = grid.n_basis;
    int* order = grid.order;
    int* support = grid.support;
    int* celltypes = grid.celltypes;
    double* basis = grid.basis;
    int* offsets = grid.offsets;

    int* newSupport = new int[n_basis];
    int* newCelltypes = new int[n_basis];
    double* newBasis = new double[n_basis];
    int* newOffsets = new int[M+1]; newOffsets[0] = 0;

    for (int i=0; i<M; i++) newOffsets[order[i]+1]= offsets[i+1] - offsets[i];
    for (int i=0; i<M; i++) newOffsets[i+1] = newOffsets[i+1] + newOffsets[i];

    for (int i=0; i<M; i++ ){
        int k = newOffsets[order[i]];
        for (int j = offsets[i]; j<offsets[i+1]; j++){
            newSupport[k] = support[j];
            newCelltypes[k] = celltypes[j];
            newBasis[k] = basis[j];
            k++;
        }
    }

    delete[] order;     grid.order = NULL;
    delete[] support;   grid.support = newSupport;
    delete[] celltypes; grid.celltypes = newCelltypes;
    delete[] basis;     grid.basis = newBasis;
    delete[] offsets;   grid.offsets = newOffsets;
}



// ********************************************************************************************** //
void reOrderBasis(Grid& grid){
    if (grid.order==NULL){ cout<<"Invalid use of reOrder()"<<endl; return; }

    double* newBasis = new double[grid.n_basis];
    int* newOffsets = new int[grid.M+1]; newOffsets[0] = 0;

    for (int i=0; i<grid.M; i++) newOffsets[grid.order[i]+1]= grid.offsets[i+1] - grid.offsets[i];
    for (int i=0; i<grid.M; i++) newOffsets[i+1] = newOffsets[i+1] + newOffsets[i];

    for (int i=0; i<grid.M; i++ ){
        int k = newOffsets[grid.order[i]];
        for (int j = grid.offsets[i]; j<grid.offsets[i+1]; j++){
            newBasis[k] = grid.basis[j];
            k++;
        }
    }

    delete[] grid.order;     grid.order = NULL;
    delete[] grid.basis;     grid.basis = newBasis;
    delete[] grid.offsets;   grid.offsets = newOffsets;
}
