#include "distribution.h"
#include "data_defs.h"
#include "utilities.h"

#include <iostream>
#include <algorithm>
#include <metis.h>
#include <mpi.h>
#include <cmath>
using namespace std;

// ********************************************************************************************** //
void makeCoarseGraph(Grid& grid, Matrix& mat, int& SIZE, Graph& graph){
    double start1, start2, start3;
    double T0 = 0, T1 = 0, T2 = 0, T3 = 0, T4 = 0, T5 = 0, T6 = 0, T7 = 0, T8 = 0;

    start1 = MPI_Wtime();

    H_data H;
    Graph fineGraph;
    int* coarseStrength = new int[SIZE*SIZE]();
    graph.n_nodes = SIZE;
    graph.nodeWeights = new int[SIZE]();
    graph.edgeOffsets = new int[SIZE+1];    graph.edgeOffsets[0] = 0; graph.edgeOffsets[1] = 0;
    graph.maxRow = mat.maxRow;
    graph.M = grid.M;
    graph.basisDistr = new int[graph.M];
    graph.n_basis_Nodes = new int[SIZE]();
    graph.package = new int[SIZE*5];

    start2 = MPI_Wtime();
    makeH(grid, H);
    T1 = MPI_Wtime() - start2;

    start2 = MPI_Wtime();
    //makeFineGraph(grid, H, fineGraph);
    makeFineGraph2(grid, H, fineGraph);
    //makeFineGraph3(grid, H, fineGraph);
    T2 = MPI_Wtime() - start2;

    start2 = MPI_Wtime();
    METIS_partition(graph, fineGraph);
    T3 = MPI_Wtime() - start2;

    // Create coarseStrength
    int fineNode = 0;
    int coarseNode, coarseEdge;
    for (int i = 0; i<fineGraph.n_edges; i++){
        if (i==fineGraph.edgeOffsets[fineNode+1]) fineNode++;
        coarseNode = graph.basisDistr[fineNode];
        coarseEdge = graph.basisDistr[fineGraph.edges[i]];
        if (coarseNode!=coarseEdge){
            if (coarseStrength[coarseNode*SIZE + coarseEdge]==0) graph.n_edges++;
            coarseStrength[coarseNode*SIZE + coarseEdge] += fineGraph.edgeWeights[i];
        }
    }
    graph.edges = new int[graph.n_edges];
    graph.edgeWeights = new int[graph.n_edges];
    for (int i=0; i<graph.M; i++) graph.nodeWeights[graph.basisDistr[i]]+=fineGraph.nodeWeights[i];
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
    for (int i = 0; i<graph.M; i++) graph.n_basis_Nodes[graph.basisDistr[i]]++;

    findSendCells(H,graph);

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
     delete[] coarseStrength;

     T8 = MPI_Wtime() - start1;

     /*
     cout
     <<"T1:\t"<<T1/T8*100<<endl
     <<"T2:\t"<<T2/T8*100<<endl
     <<"T3:\t"<<T3/T8*100<<endl;

     cout<<"Time not traced:\t"<<(T8-T1-T2-T3)/T8*100<<endl;
     */


}

// ********************************************************************************************** //
void makeH(Grid& grid, H_data& H){
    int N_TT = 0;
    int TT_max = 0;
    int* temp1 = new int[grid.N]();
    int* temp2 = new int[grid.n_basis];
    int* temp2_offsets = new int[grid.M + 1]; temp2_offsets[0] = 0;
    int k = 0;
    for (int i = 0; i < grid.M; i++){
        for (int j =  grid.offsets[i]; j < grid.offsets[i+1]; j++){
            if (grid.celltypes[j] == 2){
                if (temp1[ grid.support[j] ]==0) H.n_TT++;
                temp1[ grid.support[j] ]++;
                temp2[N_TT] = grid.support[j];
                N_TT++;
            }
        }
        temp2_offsets[i+1] = N_TT;
    }
    H.TTcells = new int[H.n_TT];
    H.offsets = new int[H.n_TT+1]; H.offsets[0] = 0;
    H.bNumbers = new int[N_TT];
    k = 0;
    for (int i = 0; i < grid.N; i++){
        if (temp1[i]!=0){
            H.TTcells[k] = i;
            H.offsets[k+1] = H.offsets[k] + temp1[i];
            if (temp1[i] > TT_max) TT_max = temp1[i];
            k++;
        }
    }
    int* temp3 = new int[grid.N*TT_max];
    int* temp4 = new int[grid.N]();
    for (int i = 0; i < grid.M; i++){
        for (int j = temp2_offsets[i]; j < temp2_offsets[i+1]; j++){
            temp3[ temp2[j]*TT_max + temp4[temp2[j]] ] = i;
            temp4[temp2[j]]++;
        }
    }
    k = 0;
    for (int i = 0; i < grid.N; i++){
        for (int j = TT_max*i; j < TT_max*i +temp4[i]; j++){
            H.bNumbers[k] = temp3[j];
            k++;
        }
    }

    H.TT_max = TT_max;
    delete[] temp1; delete temp2; delete[] temp3; delete[] temp4; delete[] temp2_offsets;
}

// ********************************************************************************************** //
void makeFineGraph(Grid& grid, H_data& H, Graph& fineGraph){
    int* strength = new int[grid.M*grid.M]();
    for (int i=0; i<H.n_TT; i++){
        int start = H.offsets[i];
        int end = H.offsets[i+1];
        for (int j=start; j<end; j++){
            for (int k = start; k<end; k++){
                if (j!=k) {
                    if (strength[grid.M*H.bNumbers[j] + H.bNumbers[k]]==0) fineGraph.n_edges++;
                    strength[grid.M*H.bNumbers[j] + H.bNumbers[k]]+=1;
                }
            }
        }
    }
    fineGraph.n_nodes = grid.M;
    fineGraph.nodeWeights = new int[grid.M];
    fineGraph.edgeOffsets = new int[grid.M+1]();
    fineGraph.edges = new int[fineGraph.n_edges];
    fineGraph.edgeWeights = new int[fineGraph.n_edges];
    for (int i = 1; i<grid.M+1; i++) fineGraph.nodeWeights[i-1] = grid.offsets[i]-grid.offsets[i-1];
    int edgeCounter = 0;
    int k = 0;      // row Number
    int l = 0;      // column number;
    for (int i = 0; i<grid.M*grid.M; i++){
        if (l==grid.M){ l = 0; k++; fineGraph.edgeOffsets[k+1]=fineGraph.edgeOffsets[k]; }
        if (strength[i]!=0){
            fineGraph.edges[edgeCounter] = l;
            fineGraph.edgeWeights[edgeCounter] = strength[i];
            fineGraph.edgeOffsets[k+1]++;
            edgeCounter++;
        }
        l++;
    }
    delete[] strength;
}

// ********************************************************************************************** //
void makeFineGraph2(Grid& grid, H_data& H, Graph& fineGraph){
    fineGraph.n_nodes = grid.M;
    fineGraph.nodeWeights = new int[grid.M];
    fineGraph.edgeOffsets = new int[grid.M+1]();
    for (int i = 1; i<grid.M+1; i++) fineGraph.nodeWeights[i-1] = grid.offsets[i]-grid.offsets[i-1];
    int maxEdges = H.TT_max*(H.TT_max-1)*H.n_TT;
    int* temp1 = new int[maxEdges];
    int power = 100;
    while(grid.M>power) power*=10;

    int c1 = 0;
    for (int i=0; i<H.n_TT; i++){
        int start = H.offsets[i];
        int end = H.offsets[i+1];
        for (int j=start; j<end; j++){
            for (int k = start; k<end; k++){
                if (j!=k) {
                    temp1[c1] = H.bNumbers[j]*power + H.bNumbers[k];
                    c1++;
                }
            }
        }
    }
    sort(temp1, temp1 + c1);

    int* temp2 = new int[maxEdges];
    temp2[0] = temp1[0];
    fineGraph.n_edges = 1;
    for (int i = 1; i < c1; i++){
        if (temp2[fineGraph.n_edges-1]!=temp1[i]){
            temp2[fineGraph.n_edges]=temp1[i];
            fineGraph.n_edges++;
        }
    }
    fineGraph.edges = new int[fineGraph.n_edges];
    fineGraph.edgeWeights = new int[fineGraph.n_edges];
    for (int i = 0; i < fineGraph.n_edges; i++ ) fineGraph.edgeWeights[i] = 1;

    fineGraph.edgeOffsets[0] = 0;
    int k = 0;
    for (int i = 0; i < fineGraph.n_edges; i++){
        int s = temp2[i]/power;
        int t = temp2[i] - power*s;
        if (s!=k){
            fineGraph.edgeOffsets[k+1] = i;
            k++;
        }
        fineGraph.edges[i] = t;
    }
    fineGraph.edgeOffsets[k+1] = fineGraph.n_edges;
}

// ********************************************************************************************** //
void makeFineGraph3(Grid& grid, H_data& H, Graph& fineGraph){
    fineGraph.n_nodes = grid.M;
    fineGraph.nodeWeights = new int[grid.M];
    fineGraph.edgeOffsets = new int[grid.M+1]();
    for (int i = 1; i<grid.M+1; i++) fineGraph.nodeWeights[i-1] = grid.offsets[i]-grid.offsets[i-1];

    int maxEdges = H.TT_max*(H.TT_max-1)*H.n_TT;
    int* temp1 = new int[maxEdges];
    int power = 10;
    while(grid.M>power) power*=10;

    bool* temp2 = new bool[grid.M*power]();

    int c1 = 0;
    for (int i=0; i<H.n_TT; i++){
        int start = H.offsets[i];
        int end = H.offsets[i+1];
        for (int j=start; j<end; j++){
            for (int k = start; k<end; k++){
                if (j!=k) {
                    int theNumber = H.bNumbers[j]*power + H.bNumbers[k];
                    if (!temp2[theNumber]){
                        temp2[theNumber] = true;
                        temp1[c1] = theNumber;
                        c1++;
                    }
                }
            }
        }
    }
    sort(temp1, temp1 + c1);

    fineGraph.n_edges = c1;
    fineGraph.edges = new int[fineGraph.n_edges];
    fineGraph.edgeWeights = new int[fineGraph.n_edges];
    for (int i = 0; i < fineGraph.n_edges; i++ ) fineGraph.edgeWeights[i] = 1;

    fineGraph.edgeOffsets[0] = 0;


    int k = 0;
    for (int i = 0; i < fineGraph.n_edges; i++){
        int s = temp1[i]/power;
        int t = temp1[i] - power*s;

        fineGraph.edges[i] = t;

        if (s!=k){
            fineGraph.edgeOffsets[k+1] = i;
            k++;
        }

    }
    fineGraph.edgeOffsets[k+1] = fineGraph.n_edges;
    delete[] temp1; delete[] temp2;
}










// ********************************************************************************************** //
void METIS_partition(Graph& graph, Graph& fineGraph){
    // See METIS manual for info about the functions.
    idx_t ncon = 1;
    idx_t objval;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_CONTIG] = 1;
    options[METIS_OPTION_MINCONN] = 1;

    METIS_PartGraphRecursive(
    //METIS_PartGraphKway(
        &graph.M,
        &ncon,
        fineGraph.edgeOffsets,
        fineGraph.edges,
        fineGraph.nodeWeights,
        NULL,
        fineGraph.edgeWeights,
        &graph.n_nodes,
        NULL,
        NULL,
        options,
        &objval,
        graph.basisDistr
    );
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
    int* temp = new int[SIZE];

    // Change H by mapping to ranks numbers instead of basis numbers
    for (int i=0; i<H.offsets[H.n_TT]; i++) {
        H.bNumbers[i] = graph.basisDistr[H.bNumbers[i]];
    }
    // Creating the new H-structure.
    // Basically strips the original H from unnececary information related to sending
    for (int i = 0; i < H.n_TT; i++){
        int start = H.offsets[i];
        int end = H.offsets[i+1];
        int counter = 0;
        for (int j = start; j < end; j++){
            bool isIn = false;
            for (int k = 0; k<counter; k++){
                if (H.bNumbers[j] == temp[k]) {
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
        graph.sendOffsets[i+1] = graph.sendOffsets[i+1] + graph.sendOffsets[i];
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
    delete[] temp;
}

// ********************************************************************************************** //
void setupNodes(Graph& graph, Node& node, int& RANK){
    // 3 sendings. This can be done better.
    if (RANK==0){
        for (int i = 1; i<graph.n_nodes; i++){
            MPI_Send(&graph.package[i*5], 5, MPI_INT, i, 666, MPI_COMM_WORLD );
        }
        for (int i = 1; i<graph.n_nodes; i++){
            MPI_Send(&graph.edges[graph.edgeOffsets[i]],
            graph.edgeOffsets[i+1]-graph.edgeOffsets[i], MPI_INT, i, 666, MPI_COMM_WORLD );
        }
        for (int i = 1; i<graph.n_nodes; i++){
            MPI_Send(&graph.sendCells[graph.sendOffsets[i]],
            graph.sendOffsets[i+1]-graph.sendOffsets[i], MPI_INT, i, 666, MPI_COMM_WORLD );
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
void distributeProblem(Grid& globalGrid, Grid& localGrid, Matrix& globalMat, Matrix& localMat,
Graph& graph, Node& node){
    // Distributes the problem in the following way:
    // Processor 0 preforms three sendings to each other processor.
    // First sending: Sends support, celltypes and offsets in Grid
    // Second sending: Sends j_index in Matrix
    // Third sending: Sends basis in Grid and conn in matrix.
    // After the sending is complete, processor 0 creates its own local data.
    // All processors also assembles loc_index and loc_conn after the sending is complete.
    if(node.number==0){
        int bNumber_start = 0;
        int start, n_basis, k, c1, c2;
        int n_sendBuff1, n_sendBuff2, n_sendBuff3;
        int* sendBuff1;
        int* sendBuff2;
        double* sendBuff3;
        int n_cUnique;
        int* cUnique;
        for (int i = 1; i < graph.n_nodes; i++){
            bNumber_start += graph.n_basis_Nodes[i-1];
            start = globalGrid.offsets[bNumber_start];
            n_basis = globalGrid.offsets[bNumber_start + graph.n_basis_Nodes[i]] - start;
            // First sending
            n_sendBuff1 = 2*n_basis + graph.n_basis_Nodes[i] + 1;
            sendBuff1 = new int[n_sendBuff1];
            k = start;
            for (int j = 0; j < n_basis; j++){
                sendBuff1[j] = globalGrid.support[k];
                k++;
            }
            k = start;
            for (int j = n_basis; j < 2*n_basis; j++){
                sendBuff1[j] = globalGrid.celltypes[k];
                k++;
            }
            k = bNumber_start;
            for (int j = 2*n_basis; j < n_sendBuff1; j++){
                sendBuff1[j] = globalGrid.offsets[k];
                k++;
            }
            MPI_Send(sendBuff1, n_sendBuff1, MPI_INT, i, 666, MPI_COMM_WORLD);
            // Find unique cellnumbers
            n_cUnique = 0;
            cUnique = new int[n_basis];
            sortNDremvDup(sendBuff1, n_basis, cUnique, n_cUnique);
            // Second sending
            n_sendBuff2 = globalMat.maxRow*n_cUnique;
            sendBuff2 = new int[n_sendBuff2];
            k = 0;                                          // index counter in cUnique
            c1 = 0;                                         // column counter
            c2 = cUnique[k];                                // cellNumber
            for (int j = 0; j < n_sendBuff2; j++){
                if (c1==globalMat.maxRow){
                    c1 = 0;
                    k++;
                    c2 = cUnique[k];
                }
                sendBuff2[j] = globalMat.j_index[ c2*globalMat.maxRow +c1 ];
                c1++;
            }
            MPI_Send(sendBuff2, n_sendBuff2, MPI_INT, i, 666, MPI_COMM_WORLD);
            // Third sending
            n_sendBuff3 = n_basis + globalMat.maxRow*n_cUnique;
            sendBuff3 = new double[n_sendBuff3];
            k = start;
            for (int j = 0; j < n_basis; j++){
                sendBuff3[j] = globalGrid.basis[k];
                k++;
            }
            k = 0;
            c1 = 0;
            c2 = cUnique[k];
            for (int j = n_basis; j < n_sendBuff3; j++){
                if (c1==globalMat.maxRow){
                    c1 = 0;
                    k++;
                    c2 = cUnique[k];
                }
                sendBuff3[j] = globalMat.conn[ c2*globalMat.maxRow +c1 ];
                c1++;
            }
            MPI_Send(sendBuff3, n_sendBuff3, MPI_DOUBLE, i, 666, MPI_COMM_WORLD);
            delete[] sendBuff1; delete[] sendBuff2; delete[] sendBuff3; delete[] cUnique;
        }
        localGrid.basis = globalGrid.basis;
        localGrid.offsets = globalGrid.offsets;
        localGrid.M = node.n_basis;
        localGrid.n_basis = node.weight;
        localGrid.support = new int[node.weight];
        localGrid.celltypes = new int[node.weight];
        localGrid.updates = new double[node.weight];
        for (int i = 0; i<localGrid.n_basis; i++){
            localGrid.support[i] = globalGrid.support[i];
            localGrid.celltypes[i] = globalGrid.celltypes[i];
        }
        delete[] globalGrid.support;    globalGrid.support = NULL;
        delete[] globalGrid.celltypes;  globalGrid.celltypes = NULL;
        localMat.maxRow = node.maxRow;
        localMat.n_basis = node.weight;
        localMat.loc_index = new int[localMat.maxRow*localMat.n_basis];
        localMat.loc_conn = new double[localMat.maxRow*localMat.n_basis];
        k = 0;
        c1 = 0;
        c2 = localGrid.support[k]*globalMat.maxRow;
        for (int j = 0; j < localGrid.n_basis*localMat.maxRow; j++){
            if (c1==localMat.maxRow){
                c1 = 0;
                k++;
                c2 = localGrid.support[k]*globalMat.maxRow;
            }
            localMat.loc_index[j] = globalMat.j_index[c2 + c1];
            localMat.loc_conn[j] = globalMat.conn[c2 + c1];
            c1++;
        }
        destroyMatrix(globalMat);
    }
    else{
        MPI_Status* status;
        int n_recvBuff1, n_recvBuff2, n_recvBuff3;
        int* recvBuff1;
        int* recvBuff2;
        double* recvBuff3;
        int n_cUnique;
        int* cUnique;
        int min, max, L, k, c1, c2;
        int* map;
        localGrid.M = node.n_basis;
        localGrid.n_basis = node.weight;
        localMat.maxRow = node.maxRow;
        localMat.n_basis = node.weight;
        localGrid.updates = new double[localGrid.n_basis];
        localGrid.basis = new double[localGrid.n_basis];
        localMat.loc_index = new int[localMat.maxRow*localGrid.n_basis];
        localMat.loc_conn = new double[localMat.maxRow*localGrid.n_basis];
        // First recieving  (support, cellnumbers and offsets in localGrid)
        n_recvBuff1 =  2*localGrid.n_basis + localGrid.M + 1;
        recvBuff1 = new int[n_recvBuff1];
        MPI_Recv(recvBuff1, n_recvBuff1, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
        localGrid.support = &recvBuff1[0];
        localGrid.celltypes = &recvBuff1[ localGrid.n_basis ];
        localGrid.offsets = &recvBuff1[ 2*localGrid.n_basis ];
        for (int i = 1; i < localGrid.M+1; i++){
            localGrid.offsets[i]=localGrid.offsets[i]-localGrid.offsets[0];
        }
        localGrid.offsets[0] = 0;
        // Find unique cellnumbers
        cUnique = new int[localGrid.n_basis];
        n_cUnique = 0;
        sortNDremvDup(localGrid.support, localGrid.n_basis, cUnique, n_cUnique);
        localGrid.N = n_cUnique;
        localMat.N = n_cUnique;
        // Second recieving (j_index in localMat)
        n_recvBuff2 = localMat.maxRow*n_cUnique;
        recvBuff2 = new int[n_recvBuff2];
        MPI_Recv(recvBuff2, n_recvBuff2, MPI_INT, 0, 666, MPI_COMM_WORLD, status );
        localMat.j_index = recvBuff2;
        // Third recieving (basis in localGrid and conn in localMat)
        n_recvBuff3 = localGrid.n_basis + localMat.maxRow*n_cUnique;
        recvBuff3 = new double[n_recvBuff3];
        MPI_Recv(recvBuff3, n_recvBuff3, MPI_DOUBLE, 0, 666, MPI_COMM_WORLD, status );
        for (int i = 0; i < localGrid.n_basis; i++) {
            localGrid.basis[i] = recvBuff3[i];
        }
        localMat.conn = &recvBuff3[localGrid.n_basis];
        // Creating loc_index and loc_conn in localMat
        min = cUnique[0];
        max = cUnique[n_cUnique-1];
        L = max - min + 1;
        map = new int [L];
        for (int i = 0; i < n_cUnique; i++) {
            map[ cUnique[i] - min ] = i;
        }
        k = 0;
        c1 = 0;
        c2 = map[ localGrid.support[k] - min ]*localMat.maxRow;
        for (int i = 0; i < localMat.maxRow*localGrid.n_basis; i++){
            if (c1==localMat.maxRow){
                c1 = 0; k++;
                c2 = map[ localGrid.support[k] - min ]*localMat.maxRow;
            }
            localMat.loc_index[i] = localMat.j_index[c2 + c1];
            localMat.loc_conn[i] = localMat.conn[c2 +c1];
            c1++;
        }
        delete[] localMat.j_index;
        localMat.j_index = NULL;
        delete[] recvBuff3;
        localMat.conn = NULL;
        delete[] map;
        delete[] cUnique;
    }
}

void sortNDremvDup(int* input, int& n_input, int* output, int& n_output){
    int max = input[0];
    int min = input[0];
    for (int i = 1; i < n_input; i ++){
        if (input[i] > max) max = input[i];
        else if (input[i] < min) min = input[i];
    }
    int L = max - min + 1;
    bool* temp1 = new bool[L]();
    n_output = 0;
    for (int i = 0; i < n_input; i++){
        if (! temp1[ input[i] - min  ]){
            temp1[input[i]- min] = true;
            output[n_output] = input[i];
            n_output++;
        }
    }
    sort(output, output + n_output );
    delete[] temp1;
}













// yolo
