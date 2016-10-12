#include "utilities.h"
#include "data_defs.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

// ********************************************************************************************** //
void print(Grid& grid){
    cout.precision(4);
    cout<<"______________________________"<<endl<<"PRINTING GRID"<<endl;
    if (grid.N!=-1) cout<<"Number of cells: "<<grid.N<<endl;
    if (grid.M!=-1) cout<<"Number of blocks: "<<grid.M<<endl;
    if (grid.n_sup!=-1) cout<<"Number of elements in support: "<<grid.n_sup<<endl;
    if (grid.offsets!=NULL){
        cout<<"Offsets:"<<endl;
        for (int i = 0; i<grid.M+1; i++){ cout<<grid.offsets[i]<<"  "; } cout<<endl;
    }
    if (grid.support!=NULL){
        cout<<"Support:"<<endl;
        int it = 1;
        for (int i = 0; i<grid.n_sup; i++){
            if (i==grid.offsets[it]){ cout<<endl<<endl; it++; }
            cout<<grid.support[i]<<"  ";
        } cout<<endl;
    }
    if (grid.boundary!=NULL){
        cout<<"boundary:"<<endl;
        int it = 1;
        for (int i = 0; i<grid.n_sup; i++){
            if (i==grid.offsets[it]){ cout<<endl<<endl; it++; }
            cout<<grid.boundary[i]<<"  ";
        } cout<<endl;
    }
    if (grid.basis!=NULL){
        cout<<"Basis:"<<endl;
        int it = 1;
        for (int i = 0; i<grid.n_sup; i++){
            if (i==grid.offsets[it]){ cout<<endl<<endl; it++; }
            cout<<grid.basis[i]<<"  ";
        } cout<<endl;
    }
    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void print(Matrix& mat){
    cout.precision(2);
    cout<<"PRINTING MATRIX"<<endl;
    if (mat.N!=-1) cout<<"Number of cells " <<mat.N<<endl;
    if (mat.maxRow!=-1) cout<<"MaxRow: "<<mat.maxRow<<endl;
    if (mat.mat_index!=NULL){
        cout<<"mat_index:"<<endl;
        for (int i = 0; i<mat.maxRow*mat.N; i++){
          cout<<mat.mat_index[i]<<"\t";
          if ((i+1)%mat.maxRow==0){cout<<endl;}
        }
    }
    if (mat.mat_coef!=NULL){
        cout<<"Coefficients:"<<endl;
        for (int i = 0; i<mat.maxRow*mat.N; i++){
          cout<<mat.mat_coef[i]<<"\t";
          if ((i+1)%mat.maxRow==0){cout<<endl;}
        }
    }
    if (mat.sup_index!=NULL && mat.sup_coef!=NULL){
      cout<<"sup_index:"<<endl;
      for (int i = 0; i<mat.maxRow*mat.n_sup; i++){
        cout<<mat.sup_index[i]<<"\t";
        if ((i+1)%mat.maxRow==0){cout<<endl;}
      }
      cout<<"sup_coef:"<<endl;
      for (int i = 0; i<mat.maxRow*mat.n_sup; i++){
        cout<<mat.sup_coef[i]<<"\t";
        if ((i+1)%mat.maxRow==0){cout<<endl;}
      }
    }
}

// ********************************************************************************************** //
void print(BInfo& B){
    cout.precision(2);
    cout<<"______________________________"<<endl<<"Printing TT"<<endl<<endl;

    if (B.n_cNumbers!=0){
        cout<<"n_cNumbers:\t"<<B.n_cNumbers<<endl;
        cout<<"cNumbers:\t";
        for (int i = 0; i<B.n_cNumbers; i++) cout<<B.cNumbers[i]<<"\t";
        cout<<endl;
        if (B.sums!=NULL){
            cout<<"sums:\t\t";
            for (int i = 0; i<B.n_cNumbers; i++) cout<<B.sums[i]<<"\t";
            cout<<endl<<endl;
        }
    }

    if (B.n_cIndices!=0){
        cout<<"n_cIndices:\t"<<B.n_cIndices<<endl;
        cout<<"cIndices:"<<endl<<"\t\t";
        for (int i = 0; i<B.n_cNumbers; i++){
            for (int j=B.cIndices_Offsets[i]; j<B.cIndices_Offsets[i+1]; j++){
                cout<<B.cIndices[j]<<"\t";
            }
            cout<<endl<<"\t\t";
        }
        cout<<endl;
    }

    if (B.n_sendCells!=0){
        cout<<"n_sendCells:\t"<<B.n_sendCells<<endl;
    }

    if (B.n_edges!=0){
        cout<<"n_edges:\t"<<B.n_edges<<endl;
        cout<<"edges:\t\t";
        for (int i = 0; i<B.n_edges; i++) cout<<B.edges[i]<<"\t";
        cout<<endl<<endl;
    }

    if (B.n_recvBuff!=0){
        cout<<"n_recvBuff:\t"<<B.n_recvBuff<<endl;
        cout<<"recvBuff:"<<endl<<"\t\t";
        for (int i = 0; i<B.n_edges; i++){
            for (int j = B.recvBuff_Offsets[i]; j<B.recvBuff_Offsets[i+1]; j++){
                cout<<B.recvBuff[j]<<"\t";
            }
            cout<<endl<<"\t\t";
        }cout<<endl<<endl;
    }

    if (B.n_recvIndices!=0){
        cout<<"n_recvIndices:\t"<<B.n_recvIndices<<endl;
        cout<<"recvIndices:"<<endl<<"\t\t";
        for (int i = 0; i<B.n_recvIndices; i++){
            for (int j = B.recvIndices_Offsets[i]; j<B.recvIndices_Offsets[i+1]; j++){
                cout<<B.recvIndices[j]<<"\t";
            }
            cout<<endl<<"\t\t";
        }cout<<endl;
    }
    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void print(H_data& H){
    cout<<"______________________________"<<endl<<"Printing H"<<endl;
    cout<<"Number of TT-cells:\t"<<H.n_TT<<endl;
    for (int i = 0; i<H.n_TT; i++){
        cout<<H.TTcells[i]<<":\t";
        for (int j = H.offsets[i]; j<H.offsets[i+1]; j++){
            cout<<H.bNumbers[j]<<"\t";
        }cout<<endl;
    }
    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void print(Graph& grap){
    cout<<"______________________________"<<endl<<"Printing Graph"<<endl;
    if (grap.n_nodes!=0) cout<<"Number of nodes:\t"<<grap.n_nodes<<endl;
    if (grap.n_edges!=0) cout<<"Number of edges:\t"<<grap.n_edges<<endl;
    if (grap.nodeWeights!=NULL){
        cout<<"Node Weights:"<<endl;
        for (int i = 0; i<grap.n_nodes; i++) cout<<grap.nodeWeights[i]<<"\t";
    }
    if (grap.edges!=NULL){
        cout<<endl<<"Edges:"<<endl;
        int k = 1;
        for (int i = 0; i<grap.n_edges; i++) {
            if (i==grap.edgeOffsets[k] ){
                k++; cout<<endl;
                if (grap.edgeOffsets[k-1] == grap.edgeOffsets[k]){
                    k++; cout<<endl;
                }
            }
            cout<<grap.edges[i]<<"\t";
        }
    }
    if (grap.edgeWeights!=NULL){
        cout<<endl<<"Edge weights:"<<endl;
        int k = 1;
        for (int i = 0; i<grap.n_edges; i++) {
            if (i==grap.edgeOffsets[k]){
                k++; cout<<endl;
            }
            cout<<grap.edgeWeights[i]<<"\t";
        }
        cout<<endl<<endl;
    }
    if (grap.maxRow!=-1) cout<<"maxRow:\t"<<grap.maxRow<<endl;
    if (grap.M!=-1) cout<<"Number of basis functions:\t"<<grap.M<<endl;
    if (grap.basisDistr!=NULL){
        cout<<"Basis distribtion:"<<endl;
        for (int i=0; i<grap.M; i++) cout<<grap.basisDistr[i]<<"\t"; cout<<endl;
        cout<<"Basis distribtion, other format:"<<endl;
        int sqrtM = sqrt(grap.M);
        for (int i = 0; i<sqrtM; i++){
            for (int j = grap.M-(i+1)*sqrtM; j<grap.M-i*sqrtM; j++) cout<<grap.basisDistr[j]<<"\t";
            cout<<endl;
        }
    }
    if (grap.order!=NULL){
        cout<<"order:"<<endl;
        for (int i=0; i < grap.M; i++) cout<<grap.order[i]<<"\t"; cout<<endl;
    }
    if (grap.n_basis_Nodes!=NULL){
        cout<<"Number of basis functions on each node:"<<endl;
        for (int i = 0; i<grap.n_nodes; i++) cout<<grap.n_basis_Nodes[i]<<"\t"; cout<<endl;
    }
    if (grap.n_sendCells!=0){
        cout<<"n_sendCells:\t"<<grap.n_sendCells<<endl;
        cout<<"Cells to be sent:"<<endl;
        for (int i = 0; i<grap.n_nodes; i++){
            cout<<"RANK:\t"<<i<<":\t";
            for (int j = grap.sendOffsets[i]; j<grap.sendOffsets[i+1]; j++){
                cout<<grap.sendCells[j]<<"\t";
            }
            cout<<endl;
        }
    }
    if (grap.package!=NULL){
        cout<<"Package:"<<endl;
        cout<<"#edges\t#support\tmaxRow\t#M\t#sendCells"<<endl;
        for (int i = 0; i<grap.n_nodes*5; i++){
            cout<<grap.package[i]<<"\t";
            if ((i+1)%5==0) cout<<endl;
            if ((i+1)%5==2) cout<<"\t";
        }
        cout<<endl;
    }

    cout<<endl<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void print(Node& node){
    cout<<"______________________________"<<endl<<"Printing Node"<<endl;
    if (node.number!=-1){
        cout<<"Node Number (RANK):\t"<<node.number<<endl;
    }
    if (node.weight!=-1){
        cout<<"weight:\t"<<node.weight<<endl;
    }
    if (node.n_edges!=0){
        cout<<"n_edges:\t"<<node.n_edges<<endl;
        cout<<"edges:\t";
        for (int i = 0; i<node.n_edges; i++) cout<<node.edges[i]<<"\t"; cout<<endl;
    }
    if (node.maxRow!=-1){
        cout<<"maxRow:\t"<<node.maxRow<<endl;
    }
    if (node.n_basis!=0){
        cout<<"n_basis:\t"<<node.n_basis<<endl;
    }
    if (node.n_sendCells!=0){
        cout<<"n_sendCells:\t"<<node.n_sendCells<<endl;
        if (node.sendCells!=NULL){
            cout<<"sendCells:\t";
            for (int i = 0; i<node.n_sendCells; i++) cout<<node.sendCells[i]<<"\t"; cout<<endl;
        }
    }
    cout<<"______________________________"<<endl;
}

// ********************************************************************************************** //
void print(Options& options){
  cout<<"PRINTING OPTIONS"<<endl<<"tolerance: "<<options.tolerance<<endl
  <<"maxIter: "<<options.maxIter<<endl<<"omega: "<<options.omega<<endl
  <<"checkTol: "<<options.checkTol<<endl;
}

// ********************************************************************************************** //
void computeDiscrepancy(Grid& grid, double* ompSol){
  double error = 0; bool mError = false;
  for (int i = 0; i<grid.n_sup; i++){
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
void printDistribution(Graph& graph, Node& node, int RANK, int SIZE){
    MPI_Barrier(MPI_COMM_WORLD);
    if (RANK==0) print(graph);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i<SIZE; i++){
        if (RANK==i) print(node);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// ********************************************************************************************** //
void writeBasisDistr(Graph& graph){

    int maxWeight = graph.nodeWeights[0];
    int minWeight = graph.nodeWeights[0];

    for (int i = 0; i < graph.n_nodes; i++){
        if (graph.nodeWeights[i]>maxWeight){
            maxWeight = graph.nodeWeights[i];
        }
        else if (graph.nodeWeights[i]<minWeight){
            minWeight = graph.nodeWeights[i];
        }
    }

    int maxEdgeWeight = graph.edgeWeights[0];
    int minEdgeWeight = graph.edgeWeights[0];

    for (int i = 0; i < graph.n_edges; i++){
        if (graph.edgeWeights[i]>maxEdgeWeight){
            maxEdgeWeight = graph.edgeWeights[i];
        }
        else if (graph.edgeWeights[i]<minEdgeWeight){
            minEdgeWeight = graph.edgeWeights[i];
        }
    }

    int* n_edgeArray = new int[graph.n_edges];

    for (int i = 0; i < graph.n_nodes; i++){
        n_edgeArray[i] = graph.edgeOffsets[i+1]-graph.edgeOffsets[i];
    }

    int totalSend = 0;

    for (int i = 0 ; i < graph.n_edges; i ++){
        totalSend += graph.edgeWeights[i];
    }

    ofstream outfile("/home/shomeb/f/fredjoha/Desktop/master-code/matlab/plots/basisDistr.txt");
    if (!outfile.is_open()) {cout << "Error opening file(s)\n";return;}

    int M = graph.M;
    int* basisDistr = graph.basisDistr;

    outfile<<M<<endl;
    for (int i = 0; i < M; i++){
        outfile<<basisDistr[i]<<" ";
    }
    outfile.close();
}
