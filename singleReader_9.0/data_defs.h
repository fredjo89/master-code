#pragma once
#include <mpi.h>

// ********************************************************************************************** //
struct Grid{
	int N;				// number of fine cells
	int M;				// number of basis functions
	int n_sup;

	int* support; 		// n_sup elements
	bool* boundary;
	double* basis;		// n_sup elements
	double* updates; 	// n_sup elements
	int* offsets; 		// M + 1 elements


	int* celltypes; 		// used only by the read from file and the OMP solver.

	Grid();
	~Grid();
};

// ********************************************************************************************** //
struct Matrix{
	int N;
	int n_sup;
	int maxRow;

	int*	mat_index;		//	j_index
	double*	mat_coef;		//	conn

	int*  	sup_index;		// loc_index
 	double*	sup_coef;		// loc_conn

	Matrix();
	~Matrix();
};

// ********************************************************************************************** //
struct BInfo{
	// LOCAL STUFF
	int n_cNumbers;         	// Length of cNumbers and sums
	int* cNumbers; 				// The cell-numbers of all type-two cells on the node
	double* sums;               // where type 2 sums are stored.

	int n_cIndices;         	// Length of cIndices
	int* cIndices;              // The index location of all type-two cells on the node
	int* cIndices_Offsets;      // Offsets array to the above

	// SENDING AND RECIEVING STUFF
	int n_sendCells;			// Number of cells that must be sent
	int n_edges;
	int	*edges;

	int n_recvBuff;
	double* recvBuff;
	int* recvBuff_Offsets;

	int n_recvIndices;
	int* recvIndices;
	int* recvIndices_Offsets;

	MPI_Status* statSend;
	MPI_Status* statRecv;
	MPI_Request* send_request;
	MPI_Request* recv_request;

	BInfo();
	~BInfo();
};

// ********************************************************************************************** //
struct H_data{
	int n_TT; 			// Number of type 2 cells
	int TT_max; 		// Maximum number of basis functions that share a singel type 2 cell.

	int* TTcells;		// Array of type 2 cell numbers

	int* offsets; 		// offsets array to bNumbers
	int* bNumbers; 		// array of basisNumbers

	H_data();
	~H_data();
};

// ********************************************************************************************** //
struct Graph{
	int n_nodes; 				// Number of nodes
	int n_edges;				// Number of edges

	int* nodeWeights;			// Node weights
	int* edgeOffsets;			// Edge offsets
	int* edges; 				// Edges
	int* edgeWeights; 			// Edge Weights

    int maxRow;
    int M; 						// Number of basis functions
    int* basisDistr; 			// basis distribution
	int* order;
    int* n_basis_Nodes;         // Number of basis functions on each node.

	// info about sendcells
	int n_sendCells;
	int* sendOffsets;
	int* sendCells;

	int* package;               // the package that will be sent in the first round

	Graph();
	~Graph();
};

// ********************************************************************************************** //
struct Node{
	int number; 			// Node number (RANK)
	int weight;             // Weight of the node
    int n_edges;            // Number of edges
    int* edges;             // Edges array.

    int maxRow;             // maxRow as in Matrix
    int n_basis;            // number of basis functions held by the node
	int n_sendCells; 		// Number of cells that must be sent.
	int* sendCells; 		// The numbers of the cells that must be sent.

    Node();
    ~Node();
};

// ********************************************************************************************** //
struct Options{
	double tolerance;
	int maxIter;
	double omega;
	int checkTol;

	int counter;
	// Used as boolean values, but are ints because MPI does not have a bool-type.
	int underTol;
	int underTol_G;

	Options();
};
















// yolo
