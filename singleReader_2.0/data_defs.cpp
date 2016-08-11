#include "data_defs.h"

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// ********************************************************************************************** //
Grid::Grid(){
    N = -1;             M = -1;             n_basis = -1;
    order = NULL;       support = NULL;     celltypes = NULL;
    basis = NULL;       updates = NULL;     offsets = NULL;
}
Grid::~Grid(){
    delete[] order;     delete[] support;   delete[] celltypes;
    delete[] basis;     delete[] updates;   delete[] offsets;
}

// ********************************************************************************************** //
Matrix::Matrix(){
    N = -1; maxRow = -1; n_basis = -1;
    conn = NULL; j_index = NULL;  loc_index = NULL; loc_conn = NULL;
}
Matrix::~Matrix() {
    delete[] conn;  delete[] j_index;    delete[] loc_index;    delete[] loc_conn;
}

// ********************************************************************************************** //
TypeTwo::TypeTwo(){
    n_cNumbers = 0;
    n_cIndices = 0;
    n_sendCells = 0;    n_edges = 0;
    n_recvBuff = 0;
    n_recvIndices = 0;

    cNumbers = NULL;        sums = NULL;
    cIndices = NULL;        cIndices_Offsets = NULL;
    edges = NULL;
    recvBuff = NULL;        recvBuff_Offsets = NULL;
    recvIndices = NULL;     recvIndices_Offsets = NULL;
    statSend = NULL;        statRecv = NULL;    send_request = NULL;    recv_request = NULL;
}
TypeTwo::~TypeTwo(){
    delete[] cNumbers;      delete[] sums;
    delete[] cIndices;      delete[] cIndices_Offsets;
    delete[] edges;
    delete[] recvBuff;      delete[] recvBuff_Offsets;
    delete[] recvIndices;   delete[] recvIndices_Offsets;
    delete[] statSend; 	    delete[] statRecv;    delete[] send_request;    delete[]recv_request;
}

// ********************************************************************************************** //
H_data::H_data(){
    n_TT = 0;
    TTcells = NULL;     offsets = NULL;     bNumbers = NULL;
}
H_data::~H_data(){
    delete[] TTcells;   delete[] offsets;   delete[] bNumbers;
}

// ********************************************************************************************** //
Graph::Graph(){
    n_nodes = 0;        n_edges = 0;
    maxRow = -1;        M = -1;
    n_sendCells = 0;
    nodeWeights = NULL;     edgeOffsets = NULL;     edges = NULL;   edgeWeights = NULL;
    basisDistr = NULL;      n_basis_Nodes = NULL;
    sendOffsets = NULL;     sendCells = NULL;
    package = NULL;
}
Graph::~Graph(){
    delete[] nodeWeights;   delete[] edgeOffsets;   delete[] edges;     delete[] edgeWeights;
    delete[] basisDistr;    delete[] n_basis_Nodes;
    delete[] sendOffsets;   delete[] sendCells;
    delete[] package;
}

// ********************************************************************************************** //
Node::Node(){
    number = -1;    weight = -1;    n_edges = 0;
    maxRow = -1;    n_basis = 0;    n_sendCells = 0;
    edges = NULL;   sendCells = NULL;
}
Node::~Node(){
    delete[] edges;     delete[] sendCells;
}

// ********************************************************************************************** //
Options::Options(){
    tolerance = 0.001;  maxIter = 500;  omega =  (double) 2/3;
    reNorm = 100;       checkTol = 1;
}

// ********************************************************************************************** //
void Grid::print(){
    cout.precision(4);
    cout<<"______________________________"<<endl<<"PRINTING GRID"<<endl;
    if (N!=-1) cout<<"Number of cells: "<<N<<endl;
    if (M!=-1) cout<<"Number of blocks: "<<M<<endl;
    if (n_basis!=-1) cout<<"Number of elements in support: "<<n_basis<<endl;
    if (order!=NULL){
        cout<<"Order:"<<endl;
        for (int i = 0; i<M; i++){cout<<order[i]<<" ";} cout<<endl;
    }
    if (offsets!=NULL){
        cout<<"Offsets:"<<endl;
        for (int i = 0; i<M+1; i++){ cout<<offsets[i]<<"  "; } cout<<endl;
    }
    if (support!=NULL){
        cout<<"Support:"<<endl;
        int it = 1;
        for (int i = 0; i<n_basis; i++){
            if (i==offsets[it]){ cout<<endl<<endl; it++; }
            cout<<support[i]<<"  ";
        } cout<<endl;
    }
    if (celltypes!=NULL){
        cout<<"Celltypes:"<<endl;
        int it = 1;
        for (int i = 0; i<n_basis; i++){
            if (i==offsets[it]){ cout<<endl<<endl; it++; }
            cout<<celltypes[i]<<"  ";
        } cout<<endl;
    }
    if (basis!=NULL){
        cout<<"Basis:"<<endl;
        int it = 1;
        for (int i = 0; i<n_basis; i++){
            if (i==offsets[it]){ cout<<endl<<endl; it++; }
            cout<<basis[i]<<"  ";
        } cout<<endl;
    }
    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void Grid::printOrder(){
    if (order==NULL){ cout<<"Default order"<<endl; return; }
    for (int i = 0; i<M; i++) cout<<order[i]<<"\t"; cout<<endl;
}

// ********************************************************************************************** //
void Matrix::print(){
    cout.precision(2);
    cout<<"PRINTING MATRIX"<<endl;
    if (N!=-1) cout<<"Number of cells " <<N<<endl;
    if (N!=-1) cout<<"MaxRow: "<<maxRow<<endl;
    if (j_index!=NULL){
        cout<<"j_index:"<<endl;
        for (int i = 0; i<maxRow*N; i++){
          cout<<j_index[i]<<"\t";
          if ((i+1)%maxRow==0){cout<<endl;}
        }
    }
    if (conn!=NULL){
        cout<<"Coefficients:"<<endl;
        for (int i = 0; i<maxRow*N; i++){
          cout<<conn[i]<<"\t";
          if ((i+1)%maxRow==0){cout<<endl;}
        }
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

// ********************************************************************************************** //
void TypeTwo::print(){
    cout.precision(2);
    cout<<"______________________________"<<endl<<"Printing TT"<<endl<<endl;

    if (n_cNumbers!=0){
        cout<<"n_cNumbers:\t"<<n_cNumbers<<endl;
        cout<<"cNumbers:\t";
        for (int i = 0; i<n_cNumbers; i++) cout<<cNumbers[i]<<"\t";
        cout<<endl;
        if (sums!=NULL){
            cout<<"sums:\t\t";
            for (int i = 0; i<n_cNumbers; i++) cout<<sums[i]<<"\t";
            cout<<endl<<endl;
        }
    }

    if (n_cIndices!=0){
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

    if (n_sendCells!=0){
        cout<<"n_sendCells:\t"<<n_sendCells<<endl;
    }

    if (n_edges!=0){
        cout<<"n_edges:\t"<<n_edges<<endl;
        cout<<"edges:\t\t";
        for (int i = 0; i<n_edges; i++) cout<<edges[i]<<"\t";
        cout<<endl<<endl;
    }

    if (n_recvBuff!=0){
        cout<<"n_recvBuff:\t"<<n_recvBuff<<endl;
        cout<<"recvBuff:"<<endl<<"\t\t";
        for (int i = 0; i<n_edges; i++){
            for (int j = recvBuff_Offsets[i]; j<recvBuff_Offsets[i+1]; j++){
                cout<<recvBuff[j]<<"\t";
            }
            cout<<endl<<"\t\t";
        }cout<<endl<<endl;
    }

    if (n_recvIndices!=0){
        cout<<"n_recvIndices:\t"<<n_recvIndices<<endl;
        cout<<"recvIndices:"<<endl<<"\t\t";
        for (int i = 0; i<n_recvIndices; i++){
            for (int j = recvIndices_Offsets[i]; j<recvIndices_Offsets[i+1]; j++){
                cout<<recvIndices[j]<<"\t";
            }
            cout<<endl<<"\t\t";
        }cout<<endl;
    }
    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void H_data::print(){
    cout<<"______________________________"<<endl<<"Printing H"<<endl;
    cout<<"Number of TT-cells:\t"<<n_TT<<endl;
    for (int i = 0; i<n_TT; i++){
        cout<<TTcells[i]<<":\t";
        for (int j = offsets[i]; j<offsets[i+1]; j++){
            cout<<bNumbers[j]<<"\t";
        }cout<<endl;
    }
    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void Graph::print(){
    cout<<"______________________________"<<endl<<"Printing Graph"<<endl;
    if (n_nodes!=0) cout<<"Number of nodes:\t"<<n_nodes<<endl;
    if (n_edges!=0) cout<<"Number of edges:\t"<<n_edges<<endl;
    if (nodeWeights!=NULL){
        cout<<"Node Weights:"<<endl;
        for (int i = 0; i<n_nodes; i++) cout<<nodeWeights[i]<<"\t";
    }
    if (edges!=NULL){
        cout<<endl<<"Edges:"<<endl;
        int k = 1;
        for (int i = 0; i<n_edges; i++) {
            if (i==edgeOffsets[k]){
                k++; cout<<endl;
            }
            cout<<edges[i]<<"\t";
        }
    }
    if (edgeWeights!=NULL){
        cout<<endl<<"Edge weights:"<<endl;
        int k = 1;
        for (int i = 0; i<n_edges; i++) {
            if (i==edgeOffsets[k]){
                k++; cout<<endl;
            }
            cout<<edgeWeights[i]<<"\t";
        }
        cout<<endl<<endl;
    }
    if (maxRow!=-1) cout<<"maxRow:\t"<<maxRow<<endl;
    if (M!=-1) cout<<"Number of basis functions:\t"<<M<<endl;
    if (basisDistr!=NULL){
        cout<<"Basis distribtion:"<<endl;
        for (int i=0; i<M; i++) cout<<basisDistr[i]<<"\t"; cout<<endl;
        cout<<"Basis distribtion, other format:"<<endl;
        int sqrtM = sqrt(M);
        for (int i = 0; i<sqrtM; i++){
            for (int j = M-(i+1)*sqrtM; j<M-i*sqrtM; j++) cout<<basisDistr[j]<<"\t";
            cout<<endl;
        }
    }
    if (n_basis_Nodes!=NULL){
        cout<<"Number of basis functions on each node:"<<endl;
        for (int i = 0; i<n_nodes; i++) cout<<n_basis_Nodes[i]<<"\t"; cout<<endl;
    }
    if (n_sendCells!=0){
        cout<<"n_sendCells:\t"<<n_sendCells<<endl;
        cout<<"Cells to be sent:"<<endl;
        for (int i = 0; i<n_nodes; i++){
            cout<<"RANK:\t"<<i<<":\t";
            for (int j = sendOffsets[i]; j<sendOffsets[i+1]; j++){
                cout<<sendCells[j]<<"\t";
            }
            cout<<endl;
        }
    }
    if (package!=NULL){
        cout<<"Package:"<<endl;
        cout<<"#edges\t#support\tmaxRow\t#M\t#sendCells"<<endl;
        for (int i = 0; i<n_nodes*5; i++){
            cout<<package[i]<<"\t";
            if ((i+1)%5==0) cout<<endl;
            if ((i+1)%5==2) cout<<"\t";
        }
        cout<<endl;
    }

    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void Node::print(){
    cout<<"______________________________"<<endl<<"Printing Node"<<endl;
    if (number!=-1){
        cout<<"Node Number (RANK):\t"<<number<<endl;
    }
    if (weight!=-1){
        cout<<"weight:\t"<<weight<<endl;
    }
    if (n_edges!=0){
        cout<<"n_edges:\t"<<n_edges<<endl;
        cout<<"edges:\t";
        for (int i = 0; i<n_edges; i++) cout<<edges[i]<<"\t"; cout<<endl;
    }
    if (maxRow!=-1){
        cout<<"maxRow:\t"<<maxRow<<endl;
    }
    if (n_basis!=0){
        cout<<"n_basis:\t"<<n_basis<<endl;
    }
    if (n_sendCells!=0){
        cout<<"n_sendCells:\t"<<n_sendCells<<endl;
        if (sendCells!=NULL){
            cout<<"sendCells:\t";
            for (int i = 0; i<n_sendCells; i++) cout<<sendCells[i]<<"\t"; cout<<endl;
        }
    }
    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void Options::print(){
  cout<<"PRINTING OPTIONS"<<endl<<"tolerance: "<<tolerance<<endl
  <<"maxIter: "<<maxIter<<endl<<"omega: "<<omega<<endl
  <<"reNorm: "<<reNorm<<endl<<"checkTol: "<<checkTol<<endl;
}
