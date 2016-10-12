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
    graph.n_nodes = SIZE;
    graph.maxRow = mat.maxRow;
    graph.M = grid.M;
    graph.basisDistr = new int[graph.M];
    graph.order = new int[grid.M];
    graph.n_basis_Nodes = new int[SIZE]();
    graph.package = new int[SIZE*5];

    // Create H.
    H_data H;
    makeH(grid, H);

    // Create f_graph ("fine graph").
    Graph f_graph;
    makeFineGraph(grid, H, f_graph);

    // Partition f_graph.
    METIS_partition(graph, f_graph);

    // Create graph.order.
    int k = 0;
    for (int i = 0; i<graph.n_nodes; i++){
        for (int j = 0; j<grid.M; j++){
            if (i == graph.basisDistr[j]){
                 graph.order[k] = j;
                k++;
            }
        }
    }

    // Create connection matrix between nodes (process graph in matrix form).
    int* connMat = new int[SIZE*SIZE]();
    int basis_n = 0;
    int node_n = graph.basisDistr[basis_n];
    for (int i = 0; i < f_graph.n_edges; i++){
        if ( i == f_graph.edgeOffsets[ basis_n+1 ] ){
            basis_n++;
            node_n = graph.basisDistr[basis_n];
        }
        int edge_n = graph.basisDistr[ f_graph.edges[i] ];
        if ( node_n != edge_n ){
            if (connMat[node_n*SIZE + edge_n] == 0){
                graph.n_edges++;
            }
            connMat[ node_n*SIZE + edge_n ] += f_graph.edgeWeights[i];
        }
    }

    graph.nodeWeights = new int[SIZE]();
    graph.edgeOffsets = new int[SIZE+1];
    graph.edgeOffsets[0] = 0;
    graph.edgeOffsets[1] = 0;
    graph.edges = new int[graph.n_edges];
    graph.edgeWeights = new int[graph.n_edges];

    for (int i = 0; i < graph.M; i++) {
        graph.nodeWeights[graph.basisDistr[i]] += f_graph.nodeWeights[i];
    }

    int edge_n = 0;
    int row = 0;
    int col = 0;
    for ( int i = 0; i < SIZE*SIZE; i++ ){
        if ( col == SIZE ){
            col = 0;
            row++;
            graph.edgeOffsets[row+1] = graph.edgeOffsets[row];
        }
        if ( connMat[i]!=0 ){
            graph.edges[edge_n] = col;
            graph.edgeWeights[edge_n] = connMat[i];
            graph.edgeOffsets[row+1]++;
            edge_n++;
        }
        col++;
    }
    for (int i = 0; i < graph.M; i++) {
        graph.n_basis_Nodes[ graph.basisDistr[i] ]++;
    }

    findSendCells(H,graph);

    k = 0;
    for (int i = 0; i < SIZE; i++){
        // number of edges for each node
        graph.package[k] = graph.edgeOffsets[i+1]-graph.edgeOffsets[i];
        k++;
        // number of elements in support for each node
        graph.package[k] = graph.nodeWeights[i];
        k++;
        // maxRow (each node must know this before recieving the matrix data)
        graph.package[k] = graph.maxRow;
        k++;
        // number of basis functions given to the node
        graph.package[k] = graph.n_basis_Nodes[i];
        k++;
        // number of cells that must be sent
        graph.package[k] = graph.sendOffsets[i+1]-graph.sendOffsets[i];
        k++;
     }

     delete[] connMat;
}

// ********************************************************************************************** //
void makeH(Grid& grid, H_data& H){
    // temp1[i]: number of basis functions having boundary-cell i in their support.
    // [ temp2[temp2_offsets[j]], temp2[temp2_offsets[j+1]-1] ]: boundary cell-numbers in
    // basis function j.
    int N_TT = 0;   // sum of boundary cells in all support regions.
    int TT_max = 0;
    int* temp1 = new int[grid.N]();
    int* temp2 = new int[grid.n_sup];
    int* temp2_offsets = new int[grid.M + 1];
    temp2_offsets[0] = 0;
    int k = 0;
    for ( int i = 0; i < grid.M; i++ ){
        for (int j =  grid.offsets[i]; j < grid.offsets[i+1]; j++){
            if ( grid.boundary[j] ){
                if (temp1[ grid.support[j] ]==0){
                    H.n_TT++;
                }
                temp1[ grid.support[j] ]++;
                temp2[N_TT] = grid.support[j];
                N_TT++;
            }
        }
        temp2_offsets[i+1] = N_TT;
    }
    H.TTcells = new int[H.n_TT];
    H.offsets = new int[H.n_TT+1];
    H.offsets[0] = 0;
    H.bNumbers = new int[N_TT];
    k = 0;
    for ( int i = 0; i < grid.N; i++ ){
        if (temp1[i] != 0){
            H.TTcells[k] = i;
            H.offsets[k+1] = H.offsets[k] + temp1[i];
            if (temp1[i] > TT_max) {
                TT_max = temp1[i];
            }
            k++;
        }
    }
    int* temp3 = new int[ grid.N*TT_max ];
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
    delete[] temp1;
    delete[] temp2;
    delete[] temp3;
    delete[] temp4;
    delete[] temp2_offsets;
}

// ********************************************************************************************** //
void makeFineGraph(Grid& grid, H_data& H, Graph& f_graph){
    // maxEdge: maximum number of edges possible.
    int maxEdges = H.TT_max*(H.TT_max-1)*H.n_TT;

    // edgeMap: all edges will be stored here, each mapped to a unique number.
    int* edgeMap = new int[maxEdges];
    int power = 100;
    while(grid.M>power) {
        power*=10;
    }

    int c1 = 0;
    for (int i=0; i<H.n_TT; i++){
        int start = H.offsets[i];
        int end = H.offsets[i+1];
        for (int j=start; j<end; j++){
            for (int k = start; k<end; k++){
                if (j!=k) {
                    edgeMap[c1] = H.bNumbers[j]*power + H.bNumbers[k];
                    c1++;
                }
            }
        }
    }

    sort(edgeMap, edgeMap + c1);

    // edgeMap_2: same as edgeMap, but duplicated edges removed
    int* edgeMap_2 = new int[maxEdges];
    // edgeWghts_temp: temp buffer storing edge weights
    int* edgeWghts_temp = new int[maxEdges];
    edgeMap_2[0] = edgeMap[0];
    edgeWghts_temp[0] = 1;
    int n_edges = 0;
    for (int i = 1; i < c1; i++){
        if (edgeMap_2[n_edges]!=edgeMap[i]){
            n_edges++;
            edgeMap_2[n_edges]=edgeMap[i];
            edgeWghts_temp[n_edges]=1;
        }
        else{
            edgeWghts_temp[n_edges]++;
        }
    }
    n_edges++;

    f_graph.n_nodes = grid.M;
    f_graph.n_edges = n_edges;

    f_graph.nodeWeights = new int[ f_graph.n_nodes ];
    f_graph.edgeOffsets = new int[ f_graph.n_nodes+1 ]();
    f_graph.edges = new int[f_graph.n_edges];
    f_graph.edgeWeights = new int[f_graph.n_edges];

    // Creating nodeWeights
    for (int i = 1; i<grid.M+1; i++) f_graph.nodeWeights[i-1] = grid.offsets[i]-grid.offsets[i-1];

    // Creating edges and edge offsets
    f_graph.edgeOffsets[0] = 0;
    f_graph.edgeOffsets[f_graph.n_nodes] = f_graph.n_edges;
    int k = 0;
    for (int i = 0; i < f_graph.n_edges; i++){
        int s = edgeMap_2[i]/power;
        int t = edgeMap_2[i] - power*s;
        if (s!=k){
            f_graph.edgeOffsets[k+1] = i;
            k++;
        }
        f_graph.edges[i] = t;
    }

    // Creating edgeWeights
    for (int i = 0; i < f_graph.n_edges; i++ ){
        f_graph.edgeWeights[i] = edgeWghts_temp[i];
    }

    delete[] edgeMap;
    delete[] edgeMap_2;
    delete[] edgeWghts_temp;
}

// ********************************************************************************************** //
void makeFineGraph2(Grid& grid, H_data& H, Graph& f_graph){
    // upper bound of the number of edges between blocks
    int maxEdges = H.TT_max*(H.TT_max-1)*H.n_TT;
    int* temp1 = new int[maxEdges];
    int power = 10;
    while( grid.M  > power ){
        power*=10;
    }
    bool* temp2 = new bool[ grid.M*power ]();
    for ( int i = 0; i < H.n_TT; i++ ){
        for ( int j = H.offsets[i]; j < H.offsets[i+1]; j++ ){
            for ( int k = H.offsets[i]; k < H.offsets[i+1]; k++ ){
                if (j != k) {
                    int c = H.bNumbers[j]*power + H.bNumbers[k];
                    if (!temp2[c]){
                        temp2[c] = true;
                        temp1[ f_graph.n_edges ] = c;
                        f_graph.n_edges++;
                    }
                }
            }
        }
    }
    sort(temp1, temp1 + f_graph.n_edges);
    f_graph.n_nodes = grid.M;
    f_graph.nodeWeights = new int[grid.M];
    f_graph.edgeOffsets = new int[grid.M+1]();
    f_graph.edges = new int[f_graph.n_edges];
    f_graph.edgeWeights = new int[f_graph.n_edges];
    for (int i = 1; i<grid.M+1; i++) {
        f_graph.nodeWeights[i-1] = grid.offsets[i]-grid.offsets[i-1];
    }
    for (int i = 0; i < f_graph.n_edges; i++ ){
            f_graph.edgeWeights[i] = 1;
    }
    f_graph.edgeOffsets[0] = 0;
    f_graph.edgeOffsets[grid.M] = f_graph.n_edges;
    int k = 0;
    for (int i = 0; i < f_graph.n_edges; i++){
        int s = temp1[i]/power;
        int t = temp1[i] - power*s;
        f_graph.edges[i] = t;
        if (s != k){
            f_graph.edgeOffsets[k+1] = i;
            k++;
        }
    }
    delete[] temp1;
    delete[] temp2;
}

// ********************************************************************************************** //
void METIS_partition(Graph& graph, Graph& f_graph){
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
        f_graph.edgeOffsets,
        f_graph.edges,
        f_graph.nodeWeights,
        NULL,
        f_graph.edgeWeights,
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
    sendH.offsets = new int[H.n_TT+1];
    sendH.offsets[0] = 0;
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
                    isIn = true;
                    break;
                }
            }
            if (!isIn){
                temp[counter] = H.bNumbers[j];
                counter++;
            }
        }
        if (counter > 1){
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
    for (int i = 0; i<SIZE; i++) {
        temp[i] = 0;
    }
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
    int tag = 100; // not used for anything, but required as parameter in MPI_Send / MPI_Recv.
    // Process 0 sends initial problem info to other processes.
    if (RANK==0){
        // Send initial package with key info.
        for (int i = 1; i<graph.n_nodes; i++){
            MPI_Send(&graph.package[i*5], 5, MPI_INT, i, tag, MPI_COMM_WORLD );
        }
        // Send dependent processes.
        for (int i = 1; i<graph.n_nodes; i++){
            MPI_Send(&graph.edges[graph.edgeOffsets[i]],
            graph.edgeOffsets[i+1]-graph.edgeOffsets[i], MPI_INT, i, tag, MPI_COMM_WORLD );
        }
        // Send cell-numbers that will be exhanged.
        for (int i = 1; i<graph.n_nodes; i++){
            MPI_Send(&graph.sendCells[graph.sendOffsets[i]],
            graph.sendOffsets[i+1]-graph.sendOffsets[i], MPI_INT, i, tag, MPI_COMM_WORLD );
        }
        // Create node.
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

        // graph.package no longer needed.
        delete[] graph.package;
        graph.package = NULL;
    }
    else{
        // Recieve initial data from process 0.
        int* buffer = new int[5];
        MPI_Status* status = new MPI_Status[2]();
        MPI_Recv(buffer, 5, MPI_INT, 0, tag, MPI_COMM_WORLD, &status[0] );
        node.number = RANK;
        node.n_edges = buffer[0];
        node.weight = buffer[1];
        node.maxRow = buffer[2];
        node.n_basis = buffer[3];
        node.n_sendCells = buffer[4];
        node.edges = new int[ node.n_edges ];
        node.sendCells = new int[ node.n_sendCells ];
        MPI_Recv(node.edges, node.n_edges, MPI_INT, 0, tag, MPI_COMM_WORLD, &status[1] );
        MPI_Recv(node.sendCells, node.n_sendCells, MPI_INT, 0, tag, MPI_COMM_WORLD, status );
        delete[] buffer;
    }
}

// ********************************************************************************************** //
// A horribly long function that distributes problem data.
void distributeProblem(Grid& g_grid, Grid& l_grid, Matrix& g_mat, Matrix& l_mat,
Graph& graph, Node& node){
    // Send from process 0.
    if(node.number==0){
        int* order = graph.order;
        int* offsets = g_grid.offsets;
        int bNumber_start = 0;
        // Looping through each process.
        for (int i = 1; i < graph.n_nodes; i++){
            bNumber_start += graph.n_basis_Nodes[i-1];
            int bNumber_end = bNumber_start +  graph.n_basis_Nodes[i];
            int start = offsets[bNumber_start];
            int n_sup = graph.nodeWeights[i];

            // First sending (support, cellnumbers and offsets in l_grid)
            int n_sendBuff1 = 2*n_sup + graph.n_basis_Nodes[i] + 1;
            int* sendBuff1 = new int[n_sendBuff1];
            int k = 0;
            for (int j = bNumber_start; j < bNumber_end; j++){
                for (int l=offsets[order[j]]; l<offsets[order[j]+1]; l ++){
                    sendBuff1[k] = g_grid.support[l];
                    sendBuff1[n_sup + k] = g_grid.boundary[l];
                    k++;
                }
            }
            k = 2*n_sup;
            sendBuff1[k] = 0;
            for (int j = bNumber_start; j < bNumber_end; j++){
                sendBuff1[k+1] = sendBuff1[k] + offsets[ order[j]+1 ] - offsets[ order[j] ];
                k++;
            }
            MPI_Send(sendBuff1, n_sendBuff1, MPI_INT, i, 666, MPI_COMM_WORLD);

            // Find unique cellnumbers
            int n_cUnique = 0;
            int* cUnique = new int[n_sup];
            sortNDremvDup(sendBuff1, n_sup, cUnique, n_cUnique);

            // Second sending (mat_index in l_mat).
            int n_sendBuff2 = g_mat.maxRow*n_cUnique;
            int* sendBuff2 = new int[n_sendBuff2];
            k = 0;
            int col = 0;
            int cell_n = cUnique[k];                              
            for (int j = 0; j < n_sendBuff2; j++){
                if (col==g_mat.maxRow){
                    col = 0;
                    k++;
                    cell_n = cUnique[k];
                }
                sendBuff2[j] = g_mat.mat_index[ cell_n*g_mat.maxRow +col ];
                col++;
            }
            MPI_Send(sendBuff2, n_sendBuff2, MPI_INT, i, 666, MPI_COMM_WORLD);

            // Third sending (basis in l_grid and mat_coef in l_mat).
            int n_sendBuff3 = n_sup + g_mat.maxRow*n_cUnique;
            double* sendBuff3 = new double[n_sendBuff3];
            k = 0;
            for (int j = bNumber_start; j < bNumber_end; j++){
                for (int l = offsets[ order[j] ]; l < offsets[ order[j]+1 ]; l ++){
                    sendBuff3[k] = g_grid.basis[l];
                    k++;
                }
            }
            k = 0;
            col = 0;
            cell_n = cUnique[k];
            for (int j = n_sup; j < n_sendBuff3; j++){
                if (col == g_mat.maxRow){
                    col = 0;
                    k++;
                    cell_n = cUnique[k];
                }
                sendBuff3[j] = g_mat.mat_coef[ cell_n*g_mat.maxRow +col ];
                col++;
            }
            MPI_Send(sendBuff3, n_sendBuff3, MPI_DOUBLE, i, 666, MPI_COMM_WORLD);
            delete[] sendBuff1; delete[] sendBuff2; delete[] sendBuff3; delete[] cUnique;

        }
        // Making l_grid and l_mat on process 0.
        l_grid.M = node.n_basis;
        l_grid.n_sup = node.weight;
        l_grid.support = new int[ l_grid.n_sup ];
        l_grid.boundary = new bool[ l_grid.n_sup ];
        l_grid.basis = new double[ l_grid.n_sup ];
        l_grid.updates = new double[ l_grid.n_sup ];
        l_grid.offsets = new int[l_grid.M + 1];
        l_grid.offsets[0] = 0;
        bNumber_start = 0;
        int bNumber_end = graph.n_basis_Nodes[0];
        int k = 0;
        for (int j = bNumber_start; j < bNumber_end; j++){
            l_grid.offsets[j+1] = l_grid.offsets[j] + offsets[order[j]+1] - offsets[order[j]];
            for (int l = offsets[ order[j] ]; l < offsets[ order[j]+1 ]; l ++){
                l_grid.support[k] = g_grid.support[l];
                l_grid.boundary[k] = g_grid.boundary[l];
                l_grid.basis[k] = g_grid.basis[l];
                k++;
            }
        }
        l_mat.maxRow = g_mat.maxRow;
        l_mat.n_sup = node.weight;
        l_mat.sup_index = new int[l_mat.maxRow*l_mat.n_sup];
        l_mat.sup_coef = new double[l_mat.maxRow*l_mat.n_sup];
        k = 0;
        int col = 0;
        int cell_n = l_grid.support[k]*g_mat.maxRow;
        for (int j = 0; j < l_grid.n_sup*l_mat.maxRow; j++){
            if (col == l_mat.maxRow){
                col = 0;
                k++;
                cell_n = l_grid.support[k]*g_mat.maxRow;
            }
            l_mat.sup_index[j] = g_mat.mat_index[cell_n + col];
            l_mat.sup_coef[j] = g_mat.mat_coef[cell_n + col];
            col++;
        }

        // Delete variables in g_grid no longer needed.
        delete[] g_grid.support;
        delete[] g_grid.boundary;
        delete[] g_grid.basis;
        delete[] g_mat.mat_index;
        delete[] g_mat.mat_coef;
        g_grid.support = NULL;
        g_grid.boundary = NULL;
        g_grid.basis = NULL;
        g_mat.mat_index = NULL;
        g_mat.mat_coef = NULL;
    }
    // Recv from process 0.
    else if (node.weight != 0){
        MPI_Status* status = new MPI_Status[3]();
        l_grid.M = node.n_basis;
        l_grid.n_sup = node.weight;

        l_grid.support = new int[l_grid.n_sup];
        l_grid.boundary = new bool[l_grid.n_sup]();
        l_grid.basis = new double[l_grid.n_sup];
        l_grid.updates = new double[l_grid.n_sup];
        l_grid.offsets = new int[l_grid.M + 1];

        l_mat.maxRow = node.maxRow;
        l_mat.n_sup = node.weight;
        l_mat.sup_index = new int[l_mat.maxRow*l_grid.n_sup];
        l_mat.sup_coef = new double[l_mat.maxRow*l_grid.n_sup];

        // First recieving  (support, cellnumbers and offsets in l_grid)
        int n_recvBuff1 =  2*l_grid.n_sup + l_grid.M + 1;
        int* recvBuff1 = new int[n_recvBuff1];
        MPI_Recv(recvBuff1, n_recvBuff1, MPI_INT, 0, 666, MPI_COMM_WORLD, &status[0] );
        for (int i = 0; i < l_grid.n_sup; i++){
            l_grid.support[i] = recvBuff1[i];
        }
        int k = 0;
        for (int i = l_grid.n_sup; i < 2*l_grid.n_sup; i++){
            if (recvBuff1[i]==1){
                l_grid.boundary[k] = true;
            }
            k++;
        }
        k = 0;
        for (int i = 2*l_grid.n_sup; i < n_recvBuff1; i++){
            l_grid.offsets[k] = recvBuff1[i];
            k++;
        }

        // Find unique cellnumbers
        int n_cUnique = 0;
        int* cUnique = new int[l_grid.n_sup];
        sortNDremvDup(l_grid.support, l_grid.n_sup, cUnique, n_cUnique);
        l_grid.N = n_cUnique;
        l_mat.N = n_cUnique;

        // Second receiving (mat_index in l_mat)
        int n_recvBuff2 = l_mat.maxRow*n_cUnique;
        l_mat.mat_index = new int[n_recvBuff2];
        MPI_Recv(l_mat.mat_index, n_recvBuff2, MPI_INT, 0, 666, MPI_COMM_WORLD, &status[1] );

        // Third receiving (basis in l_grid and mat_coef in l_mat)
        int n_recvBuff3 = l_grid.n_sup + l_mat.maxRow*n_cUnique;
        double* recvBuff3 = new double[n_recvBuff3];
        MPI_Recv(recvBuff3, n_recvBuff3, MPI_DOUBLE, 0, 666, MPI_COMM_WORLD, &status[2] );
        for (int i = 0; i < l_grid.n_sup; i++) {
            l_grid.basis[i] = recvBuff3[i];
        }
        l_mat.mat_coef = &recvBuff3[l_grid.n_sup];

        // Creating sup_index and sup_coef in l_mat
        int min = cUnique[0];
        int max = cUnique[n_cUnique-1];
        int L = max - min + 1;
        int* map = new int [L];
        for (int i = 0; i < n_cUnique; i++) {
            map[ cUnique[i] - min ] = i;
        }
        k = 0;
        int col = 0;
        int index = map[ l_grid.support[k] - min ]*l_mat.maxRow;
        for (int i = 0; i < l_mat.maxRow*l_grid.n_sup; i++){
            if (col == l_mat.maxRow){
                col = 0;
                k++;
                index = map[ l_grid.support[k] - min ]*l_mat.maxRow;
            }
            l_mat.sup_index[i] = l_mat.mat_index[index + col];
            l_mat.sup_coef[i] = l_mat.mat_coef[index +col];
            col++;
        }

        delete[] recvBuff1;
        delete[] l_mat.mat_index;
        delete[] recvBuff3;
        delete[] map;
        delete[] cUnique;

        l_mat.mat_index = NULL;
        l_mat.mat_coef = NULL;
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

// ********************************************************************************************** //
void gatherBasis(Grid& l_grid, Grid& g_grid, Graph& graph, int RANK, int SIZE){
    int tag = 100;

    if(RANK==0){
        g_grid.basis = new double[g_grid.n_sup];
        double* recvBuff = new double[g_grid.n_sup];
        for (int i = 0; i < graph.nodeWeights[0]; i++){
            recvBuff[i] = l_grid.basis[i];
        }
        MPI_Status status;
        int offsets = 0;
        for (int i = 1; i<SIZE; i++){
            offsets+=graph.nodeWeights[i-1];
            if (graph.nodeWeights[i] != 0){
                MPI_Recv(&recvBuff[offsets], graph.nodeWeights[i], MPI_DOUBLE, i, tag,
                MPI_COMM_WORLD, &status);
            }
        }

        int* ordOffsets = new int[ g_grid.M + 1 ]; ordOffsets[0] = 0;
        for (int i = 0; i < g_grid.M; i++){
            ordOffsets[i+1] = ordOffsets[i] +
            g_grid.offsets[graph.order[i]+1] - g_grid.offsets[graph.order[i]];
        }

        for (int i = 0; i < g_grid.M; i++){
            int start = g_grid.offsets[ graph.order[i] ];
            int k = 0;
            for (int j = ordOffsets[i]; j < ordOffsets[i+1]; j++){
                g_grid.basis[start + k] = recvBuff[j];
                k++;
            }
        }

        delete[] recvBuff;
        delete[] ordOffsets;
    }
    else if (l_grid.n_sup != -1){
        MPI_Send(l_grid.basis, l_grid.n_sup, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
}
