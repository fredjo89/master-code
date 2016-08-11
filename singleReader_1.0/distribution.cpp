#include "distribution.h"
#include "data_defs.h"

#include <iostream>
#include <algorithm>
#include <metis.h>
#include <mpi.h>
#include <cmath>
using namespace std;

// ********************************************************************************************** //
void distributeProblem(Grid& globalGrid, Grid& localGrid, Matrix& globalMat, Matrix& localMat, Graph& graph, Node& node, int RANK, int SIZE){
    if(RANK==0){
        int k = graph.n_basis_Nodes[0];
        for (int i = 1; i<SIZE; i++){
            int offsets = globalGrid.offsets[k];
            MPI_Send(&globalGrid.offsets[k], graph.n_basis_Nodes[i]+1, MPI_INT, i, 666, MPI_COMM_WORLD );
            MPI_Send(&globalGrid.support[offsets], graph.nodeWeights[i], MPI_INT, i, 666, MPI_COMM_WORLD );
            MPI_Send(&globalGrid.celltypes[offsets], graph.nodeWeights[i], MPI_INT, i, 666, MPI_COMM_WORLD );
            MPI_Send(&globalGrid.basis[offsets], graph.nodeWeights[i], MPI_DOUBLE, i, 666, MPI_COMM_WORLD );

            int maxRow = graph.maxRow;
            MPI_Send(&globalMat.loc_index[offsets*maxRow], graph.nodeWeights[i]*maxRow, MPI_INT, i, 666, MPI_COMM_WORLD );
            MPI_Send(&globalMat.loc_conn[offsets*maxRow], graph.nodeWeights[i]*maxRow, MPI_DOUBLE, i, 666, MPI_COMM_WORLD );

            k += graph.n_basis_Nodes[i];
        }

        localGrid.M = node.n_basis;
        localGrid.n_basis = node.weight;
        localGrid.support = new int[node.weight];
        localGrid.celltypes = new int[node.weight];
        localGrid.basis = new double[node.weight];
        localGrid.offsets = new int[node.n_basis+1];

        for (int i = 0; i<node.n_basis+1; i++){
            localGrid.offsets[i] = globalGrid.offsets[i];
        }

        for (int i = 0; i<localGrid.n_basis; i++){
            localGrid.support[i] = globalGrid.support[i];
            localGrid.celltypes[i] = globalGrid.celltypes[i];
            localGrid.basis[i] = globalGrid.basis[i];
        }

        localMat.maxRow = node.maxRow;
        localMat.n_basis = node.weight;
        localMat.loc_index = new int[localMat.maxRow*localMat.n_basis];
        localMat.loc_conn = new double[localMat.maxRow*localMat.n_basis];

        for (int i = 0; i<localGrid.n_basis*localMat.maxRow; i++){
            localMat.loc_index[i] = globalMat.loc_index[i];
            localMat.loc_conn[i] = globalMat.loc_conn[i];
        }

        // Clearing up memory which is not needed anymore
        destroyMatrix(globalMat);
    }
    else{
        localGrid.M = node.n_basis;
        localGrid.n_basis = node.weight;
        localGrid.support = new int[node.weight];
        localGrid.celltypes = new int[node.weight];
        localGrid.basis = new double[node.weight];
        localGrid.offsets = new int[node.n_basis+1];

        localMat.maxRow = node.maxRow;
        localMat.n_basis = node.weight;
        localMat.loc_index = new int[localMat.maxRow*localMat.n_basis];
        localMat.loc_conn = new double[localMat.maxRow*localMat.n_basis];

        MPI_Status* status;
        MPI_Recv(localGrid.offsets, node.n_basis+1, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
        MPI_Recv(localGrid.support, localGrid.n_basis, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
        MPI_Recv(localGrid.celltypes, localGrid.n_basis, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
        MPI_Recv(localGrid.basis, localGrid.n_basis, MPI_DOUBLE, 0, 666, MPI_COMM_WORLD, status );
        for (int i = 1; i<localGrid.M+1; i++){
            localGrid.offsets[i]=localGrid.offsets[i]-localGrid.offsets[0];
        }
        localGrid.offsets[0] = 0;

        MPI_Recv(localMat.loc_index, localMat.maxRow*localMat.n_basis, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
        MPI_Recv(localMat.loc_conn, localMat.maxRow*localMat.n_basis, MPI_DOUBLE, 0, 666, MPI_COMM_WORLD, status );
    }
}

// ********************************************************************************************** //
void setupNodes(Graph& graph, Node& node, int RANK){
    // 3 sendings. This can be done better.
    int SIZE = graph.n_nodes;
    if (RANK==0){
        for (int i = 1; i<SIZE; i++){
            MPI_Send(&graph.package[i*5], 5, MPI_INT, i, 666, MPI_COMM_WORLD );
        }
        for (int i = 1; i<SIZE; i++){
            MPI_Send(&graph.edges[graph.edgeOffsets[i]], graph.edgeOffsets[i+1]-graph.edgeOffsets[i], MPI_INT, i, 666, MPI_COMM_WORLD );
        }
        for (int i = 1; i<SIZE; i++){
            MPI_Send(&graph.sendCells[graph.sendOffsets[i]], graph.sendOffsets[i+1]-graph.sendOffsets[i], MPI_INT, i, 666, MPI_COMM_WORLD );
        }
        node.number = RANK;
        node.n_edges = graph.package[0];
        node.weight = graph.package[1];
        node.maxRow = graph.package[2];
        node.n_basis = graph.package[3];
        node.n_sendCells = graph.package[4];
        node.edges = new int[node.n_edges];
        node.sendCells = new int[node.n_sendCells];
        for (int i = 0; i<node.n_edges; i++) node.edges[i] = graph.edges[i];
        for (int i = 0; i<node.n_sendCells; i++) node.sendCells[i] = graph.sendCells[i];
    }
    else{
        int buff[5];
        MPI_Status* status;
        MPI_Recv(buff, 5, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
        node.number = RANK;
        node.n_edges = buff[0];
        node.weight = buff[1];
        node.maxRow = buff[2];
        node.n_basis = buff[3];
        node.n_sendCells = buff[4];
        node.edges = new int[node.n_edges];
        node.sendCells = new int[node.n_sendCells];
        MPI_Recv(node.edges, node.n_edges, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
        MPI_Recv(node.sendCells, node.n_sendCells, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
    }
}

// ********************************************************************************************** //
void makeCoarseGraph(Grid& grid, Matrix& mat, int SIZE, Graph& graph){
    int M = grid.M;         // Number of fine nodes (Number of blocks)

    // Create H
    H_data H;
    makeH(grid, H);

    // Create fine graph
    Graph fineGraph;
    makeFineGraph(H, grid, fineGraph);

    // Create basisDistr
    int* basisDistr = new int[M];
    METIS_partition(M,fineGraph.nodeWeights, fineGraph.n_edges, fineGraph.edgeOffsets, fineGraph.edges, fineGraph.edgeWeights, basisDistr, SIZE );

    // Create graph
    int* coarseStrength = new int[SIZE*SIZE]();
    int fineNode = 0;
    int coarseNode, coarseEdge;
    for (int i = 0; i<fineGraph.n_edges; i++){
        if (i==fineGraph.edgeOffsets[fineNode+1]) fineNode++;
        coarseNode = basisDistr[fineNode];
        coarseEdge = basisDistr[fineGraph.edges[i]];
        if (coarseNode!=coarseEdge){
            if (coarseStrength[coarseNode*SIZE + coarseEdge]==0) graph.n_edges++;
            coarseStrength[coarseNode*SIZE + coarseEdge] += fineGraph.edgeWeights[i];
        }
    }

    graph.n_nodes = SIZE;
    graph.nodeWeights = new int[SIZE]();
    graph.edgeOffsets = new int[SIZE+1];
    graph.edgeOffsets[0] = 0; graph.edgeOffsets[1] = 0;
    graph.edges = new int[graph.n_edges];
    graph.edgeWeights = new int[graph.n_edges];

    for (int i = 0; i<M; i++) graph.nodeWeights[basisDistr[i]] += fineGraph.nodeWeights[i];

    int n_edges = 0;   // edge number
    int k = 0;         // row Number
    int l = 0;         // column number
    for (int i = 0; i<SIZE*SIZE; i++){
        if (l==SIZE){ l = 0; k++; graph.edgeOffsets[k+1]=graph.edgeOffsets[k]; }
        if (coarseStrength[i]!=0){
            graph.edges[n_edges] = l;
            graph.edgeWeights[n_edges] = coarseStrength[i];
            graph.edgeOffsets[k+1]++;
            n_edges++;
        }
        l++;
    }
    delete[] coarseStrength;

    // CREATE THE REST OF THE STUFF NEEDED IN GRAPH
    graph.maxRow = mat.maxRow;
    graph.M = M;
    graph.basisDistr = basisDistr;
    graph.n_basis_Nodes = new int[SIZE]();
    for (int i = 0; i<M; i++) graph.n_basis_Nodes[basisDistr[i]]++;

    findSendCells(H,graph);

    graph.package = new int[SIZE*5];
    k = 0;
    for (int i = 0; i<SIZE; i++){
        // number of edges for each node
        graph.package[k]=graph.edgeOffsets[i+1]-graph.edgeOffsets[i];   k++;
        // number of elements in support for each node
        graph.package[k]=graph.nodeWeights[i];                          k++;
        // maxRow (each node must know this before recieving the matrix data)
        graph.package[k]=graph.maxRow;                                  k++;
        // number of basis functions given to the node
        graph.package[k]=graph.n_basis_Nodes[i];                        k++;
        // number of cells that must be sent
        graph.package[k]=graph.sendOffsets[i+1]-graph.sendOffsets[i];   k++;
     }

}

// ********************************************************************************************** //
void findSendCells(H_data& H, Graph& graph){
    int SIZE = graph.n_nodes;


    // Making a new temperarily H-structure
    H_data sendH;
    sendH.TTcells = new int[H.n_TT];
    sendH.offsets = new int[H.n_TT+1];   sendH.offsets[0] = 0;
    sendH.bNumbers = new int[H.offsets[H.n_TT]];

    // Allocating memory to graph.sendOffsets
    graph.sendOffsets = new int[SIZE+1]();

    // A temp-variables
    int temp[SIZE];


    // Change H by mapping to ranks numbers instead of basis numbers
    for (int i=0; i<H.offsets[H.n_TT]+1; i++){
        H.bNumbers[i] = graph.basisDistr[H.bNumbers[i]];
    }

    // Creating the new H-structure.
    // Basically strips the original H from unnececary information related to sending
    for (int i = 0; i< H.n_TT; i++){
        int counter = 0;
        int start = H.offsets[i];
        int end = H.offsets[i+1];

        for (int j = start; j<end; j++){
            bool isIn = false;
            for (int k = 0; k<counter; k++){
                if (H.bNumbers[j]==temp[k]) {
                    isIn = true;    break;
                }
            }
            if (!isIn){
                temp[counter] = H.bNumbers[j];  counter++;
            }
        }
        if (counter>1){
            sendH.TTcells[sendH.n_TT] = H.TTcells[i];
            start = sendH.offsets[sendH.n_TT];
            sendH.offsets[sendH.n_TT + 1] = start + counter;

            for (int j = 0; j<counter; j++){
                sendH.bNumbers[start+j] = temp[j];
                graph.sendOffsets[temp[j]+1]++;
            }
            sendH.n_TT++;
        }

    }



    // Creating n_sendCells, sendOffsets and sendCells belonging to graph.
    for (int i = 0; i<SIZE; i++){
        graph.sendOffsets[i+1] = graph.sendOffsets[i+1]+graph.sendOffsets[i];
    }

    graph.n_sendCells = sendH.offsets[sendH.n_TT];
    graph.sendCells = new int[graph.n_sendCells];

    for (int i = 0; i<SIZE; i++) temp[i] = 0;
    for (int i = 0; i<sendH.n_TT; i++){
        for (int j = sendH.offsets[i]; j<sendH.offsets[i+1]; j++){
            int rank = sendH.bNumbers[j];
            graph.sendCells[ graph.sendOffsets[rank] + temp[rank] ] = sendH.TTcells[i];
            temp[rank]++;
        }
    }

}


// ********************************************************************************************** //
void makeFineGraph(H_data& H, Grid& grid, Graph& fineGraph){
    int M = grid.M;
    int* offsets = grid.offsets;

    fineGraph.nodeWeights = new int[M];
    fineGraph.edgeOffsets = new int[M+1]();
    fineGraph.edges = new int[M*M];
    fineGraph.edgeWeights = new int[M*M]();

    for (int i = 1; i<M+1; i++) fineGraph.nodeWeights[i-1] = offsets[i]-offsets[i-1];

    int* strength = new int[M*M]();

    for (int i=0; i<H.n_TT; i++){
        int start = H.offsets[i];
        int end = H.offsets[i+1];
        for (int j=start; j<end; j++){
            for (int k = start; k<end; k++){
                if (j!=k) strength[M*H.bNumbers[j] + H.bNumbers[k]]+=1;
            }
        }
    }

    int n_edges = 0;
    int k = 0;      // row Number
    int l = 0;      // column number;
    for (int i = 0; i<M*M; i++){
        if (l==M){ l = 0; k++; fineGraph.edgeOffsets[k+1]=fineGraph.edgeOffsets[k]; }
        if (strength[i]!=0){
            fineGraph.edges[n_edges] = l;
            fineGraph.edgeWeights[n_edges] = strength[i];
            fineGraph.edgeOffsets[k+1]++;
            n_edges++;
        }
        l++;
    }

    fineGraph.n_nodes = M;
    fineGraph.n_edges = n_edges;

    delete[] strength;
}


// ********************************************************************************************** //
void METIS_partition(int M, int* nodeWeights, int n_edges, int* edgeOffsets, int* edges, int* edgeWeights, int* basisDistr, int SIZE ){

    // Number of vertices.
    idx_t nvtxs = M;
    // Number of balancing constraints. Should be at least 1.
    idx_t ncon = 1;
    // Edge Offsets (offsets of adjncy)
    idx_t xadj[M+1];        for (int i = 0; i<M+1; i++) xadj[i] = edgeOffsets[i];
    // Edges
    idx_t adjncy[n_edges];  for (int i = 0; i<n_edges; i++) adjncy[i] = edges[i];
    // Veights of vertices
    idx_t vwgt[M];          for (int i = 0; i<M; i++) vwgt[i] = nodeWeights[i];
    // See in manual
    idx_t* vsize = NULL;
    // Weight of edges
    idx_t adjwgt[n_edges];  for (int i = 0; i<n_edges; i++) adjwgt[i] = edgeWeights[i];
    // Number of parts to partion the graph.
    idx_t nparts = SIZE;
    // See in manual
    real_t* tpwgts = NULL;
    // See in manual
    real_t* ubvec = NULL;
    // Return value, total communication volume
    idx_t objval;
    // Return array, the partition.
    idx_t* part = new idx_t[M];

    // options, see in manual
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_CONTIG] = 1;
    options[METIS_OPTION_MINCONN] = 1;

    //METIS_PartGraphRecursive( &nvtxs, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, &nparts, tpwgts, ubvec, options, &objval, part);

    METIS_PartGraphKway( &nvtxs, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, &nparts, tpwgts, ubvec, options, &objval, part);

    for (int i = 0; i<M; i++ ) basisDistr[i] = part[i];
}

// ********************************************************************************************** //
void makeH(Grid& grid, H_data& H){
  int N = grid.N;
  int M = grid.M;
  int n_basis = grid.n_basis;
  int* support = grid.support;
  int* celltypes = grid.celltypes;
  int* offsets = grid.offsets;

  int TTcells[n_basis];                // cell numbers of type 2 cells will be stored here;
  int* bNumbers = new int[n_basis];     // corresponding basis numbers will be stored here;
  int N_TT = 0;                         // The total number of type 2 cells will be stored here.

  int power1 = 100;
  int power2 = 100;
  while (M>power1) power1*=10;
  while (N>power2) power2*=10;
  int base = power1*power2;

  int k = 0; // basis number counter;
  for (int i = 0; i<n_basis; i++){
    if (i==offsets[k+1]) k++;
    if (celltypes[i]==2){
      TTcells[N_TT] = support[i]*power1 +k;
      N_TT++;
    }
  }

  sort(TTcells, TTcells+N_TT);
  for (int i = 0; i<N_TT; i++){
    double temp = TTcells[i];
    TTcells[i] = temp/power1;
    bNumbers[i] = temp - TTcells[i]*power1;
  }

  H.TTcells = new int[N_TT]; H.TTcells[0] = TTcells[0];
  H.offsets = new int[N_TT]; H.offsets[0] = 0;
  H.bNumbers = bNumbers;

  for (int i = 1; i<N_TT; i++){
      if (H.TTcells[H.n_TT]!=TTcells[i]){
          H.n_TT++;
          H.TTcells[H.n_TT]=TTcells[i];
          H.offsets[H.n_TT] = i;
      }
  }
  H.n_TT++;
  H.offsets[H.n_TT] = N_TT;
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
void destroyMatrix(Matrix& mat){
  delete[] mat.conn; delete[] mat.j_index;  delete[] mat.loc_index; delete[] mat.loc_conn;
  mat.conn = NULL; mat.j_index = NULL;  mat.loc_index = NULL; mat.loc_conn = NULL;
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














// yolo
