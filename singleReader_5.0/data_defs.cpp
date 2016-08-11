#include "data_defs.h"

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// ********************************************************************************************** //
Grid::Grid(){
    N = -1;
    M = -1;
    n_sup = -1;

    support = NULL;
    celltypes = NULL;
    basis = NULL;
    updates = NULL;
    offsets = NULL;
}
Grid::~Grid(){
    delete[] support;
    delete[] celltypes;
    delete[] basis;
    delete[] updates;
    delete[] offsets;
}

// ********************************************************************************************** //
Matrix::Matrix(){
    N = -1;
    maxRow = -1;
    n_sup = -1;

    conn = NULL;
    j_index = NULL;

    loc_index = NULL;
    loc_conn = NULL;
}
Matrix::~Matrix() {
    delete[] conn;
    delete[] j_index;
    delete[] loc_index;
    delete[] loc_conn;
}

// ********************************************************************************************** //
TypeTwo::TypeTwo(){
    n_cNumbers = 0;
    cNumbers = NULL;
    sums = NULL;

    n_cIndices = 0;
    cIndices = NULL;
    cIndices_Offsets = NULL;

    n_sendCells = 0;
    n_edges = 0;
    edges = NULL;

    n_recvBuff = 0;
    recvBuff = NULL;
    recvBuff_Offsets = NULL;

    n_recvIndices = 0;
    recvIndices = NULL;
    recvIndices_Offsets = NULL;

    statSend = NULL;
    statRecv = NULL;
    send_request = NULL;
    recv_request = NULL;
}

TypeTwo::~TypeTwo(){
    delete[] cNumbers;
    delete[] sums;

    delete[] cIndices;
    delete[] cIndices_Offsets;

    delete[] edges;

    delete[] recvBuff;
    delete[] recvBuff_Offsets;

    delete[] recvIndices;
    delete[] recvIndices_Offsets;

    delete[] statSend;
    delete[] statRecv;
    delete[] send_request;
    delete[]recv_request;
}

// ********************************************************************************************** //
H_data::H_data(){
    n_TT = 0;
    TT_max = 0;
    TTcells = NULL;

    offsets = NULL;
    bNumbers = NULL;
}
H_data::~H_data(){
    delete[] TTcells;
    delete[] offsets;
    delete[] bNumbers;
}

// ********************************************************************************************** //
Graph::Graph(){
    n_nodes = 0;
    n_edges = 0;
    nodeWeights = NULL;
    edgeOffsets = NULL;
    edges = NULL;
    edgeWeights = NULL;

    maxRow = -1;
    M = -1;
    basisDistr = NULL;
    order = NULL;
    n_basis_Nodes = NULL;

    n_sendCells = 0;
    sendOffsets = NULL;
    sendCells = NULL;

    package = NULL;
}
Graph::~Graph(){
    delete[] nodeWeights;
    delete[] edgeOffsets;
    delete[] edges;
    delete[] edgeWeights;

    delete[] basisDistr;
    delete[] order;
    delete[] n_basis_Nodes;

    delete[] sendOffsets;
    delete[] sendCells;

    delete[] package;
}

// ********************************************************************************************** //
Node::Node(){
    number = -1;
    weight = -1;
    n_edges = 0;
    edges = NULL;

    maxRow = -1;
    n_basis = 0;
    n_sendCells = 0;
    sendCells = NULL;
}
Node::~Node(){
    delete[] edges;
    delete[] sendCells;
}

// ********************************************************************************************** //
Options::Options(){
    tolerance = 0.001;
    maxIter = 500;
    omega =  (double) 2/3;
    checkTol = 1;
}
