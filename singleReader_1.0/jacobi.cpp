#include "jacobi.h"
#include "data_defs.h"

#include <mpi.h>
#include <iostream>
#include <cmath>
using namespace std;

// ********************************************************************************************** //
void setupSupportMatrix(Grid& grid, Matrix& mat){
    int maxRow = mat.maxRow;
    int n_basis = mat.n_basis;
    int* support = grid.support;

    mat.loc_index = new int[maxRow*n_basis];
    mat.loc_conn = new double[maxRow*n_basis];

    for (int i = 0; i<n_basis; i++){
        int k = support[i];
        for (int j = 0; j<maxRow; j++){
            mat.loc_index[i*maxRow + j] = mat.j_index[k*maxRow+j];
            mat.loc_conn[i*maxRow + j] = mat.conn[k*maxRow+j];
        }
    }
}

// ********************************************************************************************** //
void makeMapping(Grid& grid, Matrix& mat){
    int M = grid.M;
    int maxRow = mat.maxRow;

    int* offsets = grid.offsets;
    int* support = grid.support;
    int* loc_index = mat.loc_index;

    int max, min, index;
    int* myMap;

    for (int i = 0; i<M; i++){
        int start = offsets[i];
        int end = offsets[i+1];
        int length = end - start;

        max = support[start];   min = support[start];
        for (int j = start; j<end; j++ ){
            if (support[j]>max) max = support[j];
            else if (support[j]<min) min = support[j];
        }

        myMap = new int[max-min+1];
        for (int j=0; j<max-min+1; j++) myMap[j] = -2;
        for (int j=start; j<end; j++) myMap[support[j]-min] = j;

        for (int j=start*maxRow; j<end*maxRow; j++){
          index = loc_index[j];
            if (index!=-1){
              if (index>=min && index<=max) loc_index[j] = myMap[index - min];
              else loc_index[j] = -2;
            }
        }

        delete[] myMap;
    }
}

// ********************************************************************************************** //
void jacobi(Grid& grid, Matrix& mat, double omega){
    int maxRow = mat.maxRow;
    int n_basis = mat.n_basis;
    double* basis = grid.basis;
    int* loc_index = mat.loc_index;
    double* loc_conn = mat.loc_conn;

    double* update = new double[n_basis];

    for (int i = 0; i<n_basis; i++){
        update[i] = basis[i];
        int it = i*maxRow;
        for (int j = 0; j<maxRow; j++){
            int index = loc_index[it+j];
            if (index==-1) break;
            if (index!=-2){ update[i]+=loc_conn[it+j]*basis[index]; }
        }
    }

    for (int i=0; i<n_basis; i++){
      basis[i]-=omega*update[i];
    }


    delete[] update;
}

void makeTT(Grid& grid, Node& node, TT& tt){
    // unpacking
    int n_sendCells = node.n_sendCells;
    int* sendCells = node.sendCells;
    int n_edges = node.n_edges;
    int* edges = node.edges;

    // ****************************************************************************************** //
    // MAKING THE VARIABLES:
    int* cNumbers;              // The cell-numbers of all type-two cells on the node
    int* cIndices;              // The index location of all type-two cells on the node
    int* cIndices_Offsets;      // Offsets array to the above
    double* sums;                  // where type 2 sums are stored.
    int n_cNumbers = 0;         // Length of cNumbers
    int n_cIndices = 0;         // Length of cIndices

    // some temp-variables
    int* temp1 = new int[grid.n_basis];
    int* temp2 = new int[grid.n_basis];
    int* temp3;
    int* temp4;

    // finding all type-two cells and where they are placed
    for (int i = 0; i<grid.n_basis; i++){
        if (grid.celltypes[i]==2){
            temp1[n_cIndices] = grid.support[i];
            temp2[n_cIndices] = i;
            n_cIndices++;
        }
    }

    // allocate memory
    cNumbers = new int[n_cIndices];
    cIndices_Offsets = new int[n_cIndices+1]();
    cIndices = new int[n_cIndices];

    // copying sendCells
    for (int i = 0; i<n_sendCells; i++){
        cNumbers[i] = sendCells[i];
        n_cNumbers++;
    }

    temp3 = new int[n_cIndices];
    for (int i=0; i<n_cIndices; i++){
        bool isIn = false;
        for (int j=0; j<n_cNumbers; j++){
            if (cNumbers[j]==temp1[i]){
                temp3[i] = j;
                cIndices_Offsets[j+1]++;
                isIn = true;    break;
            }
        }
        if (!isIn){
            temp3[i] = n_cNumbers;
            cNumbers[n_cNumbers]=temp1[i];
            cIndices_Offsets[n_cNumbers+1]++;
            n_cNumbers++;
        }
    }

    for (int i = 0; i<n_cNumbers; i++) cIndices_Offsets[i+1]+=cIndices_Offsets[i];

    temp4 = new int[n_cNumbers]();
    for (int i = 0; i<n_cIndices; i++){
        int cIndex = temp3[i];
        cIndices[ cIndices_Offsets[cIndex] + temp4[cIndex] ] = temp2[i];
        temp4[cIndex]++;
    }

    sums = new double[n_cNumbers];
    delete[] temp1; delete[] temp2; delete[] temp3; delete[] temp4;
    // ****************************************************************************************** //
    // SETTING UP SENDING

    tt.statSend = new MPI_Status[n_edges];;
    tt.statRecv = new MPI_Status[n_edges];;
    tt.send_request = new MPI_Request[n_edges];
    tt.recv_request = new MPI_Request[n_edges];

    int* recvBuff_Offsets = new int[n_edges + 1]; recvBuff_Offsets[0] = 0;

    // First sending
    for (int i = 0; i<n_edges; i++){
        MPI_Isend(&n_sendCells, 1, MPI_INT, edges[i], 666,MPI_COMM_WORLD, &tt.send_request[i]);
        MPI_Irecv(&recvBuff_Offsets[i+1], 1, MPI_INT, edges[i], 666, MPI_COMM_WORLD, &tt.recv_request[i]);
    }
    MPI_Waitall(n_edges, tt.send_request, tt.statSend );
    MPI_Waitall(n_edges, tt.recv_request, tt.statRecv);


    for (int i = 0; i<n_edges; i++) recvBuff_Offsets[i+1]+=recvBuff_Offsets[i];
    int n_recvBuff = recvBuff_Offsets[n_edges];
    int* recvBuff = new int [n_recvBuff];

    // Second sending
    for (int i = 0; i<n_edges; i++){
        MPI_Isend(cNumbers, n_sendCells, MPI_INT, edges[i], 666,MPI_COMM_WORLD, &tt.send_request[i]);
        MPI_Irecv(&recvBuff[recvBuff_Offsets[i]], recvBuff_Offsets[i+1]-recvBuff_Offsets[i], MPI_INT, edges[i], 666, MPI_COMM_WORLD, &tt.recv_request[i]);
    }
    MPI_Waitall(n_edges, tt.send_request, tt.statSend);
    MPI_Waitall(n_edges, tt.recv_request, tt.statRecv);

    int* temp_Offsets = new int[n_sendCells+1]; temp_Offsets[0] = 0;
    int* tempIndices = new int[n_recvBuff];
    int n_tempIndices = 0;

    for (int i = 0; i<n_sendCells; i++){
        temp_Offsets[i+1] = temp_Offsets[i];
        for (int j = 0; j<n_recvBuff; j++){
            if (recvBuff[j] == cNumbers[i]){
                tempIndices[n_tempIndices] = j;
                n_tempIndices++;
                temp_Offsets[i+1]++;
            }
        }
    }

    // ****************************************************************************************** //
    // CREATING THE TT STRUCTURE

    tt.n_cNumbers = n_cNumbers;
    tt.cNumbers = cNumbers;
    tt.sums = sums;

    tt.n_cIndices = n_cIndices;
    tt.cIndices = cIndices;
    tt.cIndices_Offsets = cIndices_Offsets;

    tt.n_sendCells = n_sendCells;
    tt.n_edges = node.n_edges;
    tt.edges =  node.edges;

    tt.n_recvBuff = n_recvBuff;
    tt.recvBuff = new double[n_recvBuff];
    tt.recvBuff_Offsets = recvBuff_Offsets;

    tt.n_recvIndices = n_tempIndices;
    tt.recvIndices = tempIndices;
    tt.recvIndices_Offsets = temp_Offsets;

}

// ********************************************************************************************** //
void localSum(Grid& grid, TT& tt, int RANK){
    for (int i = 0; i<tt.n_cNumbers; i++) tt.sums[i] = 0;

    int* cIndices = tt.cIndices;
    int* cIndices_Offsets = tt.cIndices_Offsets;

    for (int i = 0; i<tt.n_cNumbers; i++){
        for (int j =cIndices_Offsets[i]; j<cIndices_Offsets[i+1]; j++ ){
            tt.sums[i]+=grid.basis[cIndices[j]];
        }

    }


}

// ********************************************************************************************** //
void sendAndRecieve(Grid& grid, TT& tt, int RANK){

    int n_edges = tt.n_edges;
    int* edges = tt.edges;

    int n_sendCells = tt.n_sendCells;
    double* sums = tt.sums;

    int n_recvBuff = tt.n_recvBuff;
    double* recvBuff = tt.recvBuff;
    int* recvBuff_Offsets = tt.recvBuff_Offsets;


    for (int i = 0; i<n_edges; i++){
        MPI_Isend(sums, n_sendCells, MPI_DOUBLE, edges[i], 666,MPI_COMM_WORLD, &tt.send_request[i]);
        MPI_Irecv(&recvBuff[recvBuff_Offsets[i]], recvBuff_Offsets[i+1]-recvBuff_Offsets[i], MPI_DOUBLE, edges[i], 666, MPI_COMM_WORLD, &tt.recv_request[i]);
    }
    MPI_Waitall(n_edges, tt.send_request, tt.statSend);
    MPI_Waitall(n_edges, tt.recv_request, tt.statRecv);

    int n_recvIndices = tt.n_recvIndices;
    int* recvIndices = tt.recvIndices;
    int* recvIndices_Offsets = tt.recvIndices_Offsets;

    for (int i = 0; i<n_sendCells; i++){
        for (int j=recvIndices_Offsets[i]; j<recvIndices_Offsets[i+1]; j++ ){
            sums[i] += recvBuff[recvIndices[j]];

        }
    }

}

// ********************************************************************************************** //
void TTnormalize(Grid& grid, TT& tt, int RANK){
    int* cIndices = tt.cIndices;
    int* cIndices_Offsets = tt.cIndices_Offsets;

    for (int i = 0; i<tt.n_cNumbers; i++){
        for (int j =cIndices_Offsets[i]; j<cIndices_Offsets[i+1]; j++ ){
            grid.basis[cIndices[j]]/=tt.sums[i];
        }

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

        for (int i = 0; i<localGrid.n_basis; i++){
            globalGrid.basis[i] = localGrid.basis[i];
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




void TT::print(){
    cout.precision(2);
    cout<<"______________________________"<<endl<<"Printing TT"<<endl<<endl;

    if (n_cNumbers!=-1){
        cout<<"n_cNumbers:\t"<<n_cNumbers<<endl;
        cout<<"cNumbers:\t";
        for (int i = 0; i<n_cNumbers; i++) cout<<cNumbers[i]<<"\t";
        cout<<endl;
        cout<<"sums:\t\t";
        for (int i = 0; i<n_cNumbers; i++) cout<<sums[i]<<"\t";
        cout<<endl<<endl;
    }

    if (n_cIndices!=-1){
        cout<<"n_cIndices:\t"<<n_cIndices<<endl;
        cout<<"cIndices:"<<endl<<"\t\t";
        for (int i = 0; i<n_cNumbers; i++){
            for (int j=cIndices_Offsets[i]; j<cIndices_Offsets[i+1]; j++){
                cout<<cIndices[j]<<"\t";
            }
            cout<<endl<<"\t\t";
        }
        cout<<endl;
    }

    if (n_sendCells!=-1){
        cout<<"n_sendCells:\t"<<n_sendCells<<endl;
    }

    if (n_edges!=-1){
        cout<<"n_edges:\t"<<n_edges<<endl;
        cout<<"edges:\t\t";
        for (int i = 0; i<n_edges; i++) cout<<edges[i]<<"\t";
        cout<<endl<<endl;
    }

    if (n_recvBuff!=-1){
        cout<<"n_recvBuff:\t"<<n_recvBuff<<endl;
        cout<<"recvBuff:"<<endl<<"\t\t";
        for (int i = 0; i<n_edges; i++){
            for (int j = recvBuff_Offsets[i]; j<recvBuff_Offsets[i+1]; j++){
                cout<<recvBuff[j]<<"\t";
            }
            cout<<endl<<"\t\t";
        }
        cout<<endl<<endl;
    }

    if (n_recvIndices!=-1){
        cout<<"n_recvIndices:\t"<<n_recvIndices<<endl;
        cout<<"recvIndices:"<<endl<<"\t\t";
        for (int i = 0; i<n_recvIndices; i++){
            for (int j = recvIndices_Offsets[i]; j<recvIndices_Offsets[i+1]; j++){
                cout<<recvIndices[j]<<"\t";
            }
            cout<<endl<<"\t\t";
        }
        cout<<endl;
    }

    cout<<"______________________________"<<endl;
}












//yolo
