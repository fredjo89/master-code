// Data structures.

#pragma once
#include <mpi.h>

// ********************************************************************************************** //
// Holds info of grid
struct Grid{
	int N;				// # fine cells.
	int M;				// # basis functions.
	int n_sup;			// Sum of elements in any support region.

	int* support; 		// Cell-numbers in support regions (size: n_sup)
	bool* boundary;		// Type of cell; true if in union boundary, false if not (size: n_sup)
	double* basis;		// basis functions (size: n_sup)
	double* updates; 	// update buffer (size: n_sup)
	int* offsets; 		// Offsets of above arrays (size: M+1)

	Grid();
	~Grid();
};

// ********************************************************************************************** //
// Holds info of matrix
struct Matrix{
	int N;		// # fine cells.
	int n_sup;	// Sum of elements in any support region.
	int maxRow;	// Max number of nonzero entries in any row of the matrix, excluding diagonal.

	int*	mat_index;		//	Colum-indices of matrix (size: maxRow*N).
	double*	mat_coef;		//	Entries of matrix (size: maxRow*N).

	int*  	sup_index;		// Support Colum-indices / location in memory.
 	double*	sup_coef;		// Support matrix entries.

	Matrix();
	~Matrix();
};

// ********************************************************************************************** //
// Holds info of boundary cells, and of sending and receiving of data.
struct BInfo{
	// Local variables.
	int n_cNumbers;         	// # cNumbers / sums.
	int* cNumbers; 				// Cell-numbers of all boundary cells on the process (uniquely).
	double* sums;               // Sums of all boundary cells on the process.

	int n_cIndices;         	// # cIndices
	int* cIndices;              // Location of all boundary cells on the process (size: n_cIndices).
	int* cIndices_Offsets;      // Offsets of cIndices / where one cell starts and another begins.

	// Sending and receiving variables (only used for more than one MPI process).
	int n_sendCells;			// # cells that must be sent.
	int n_edges;				// # processes to exchange data with.
	int	*edges;					// Processes to exhange data with.

	int n_recvBuff;				// # recvBuff.
	double* recvBuff;			// Buffer to store recieved data.
	int* recvBuff_Offsets;		// Offsets of recvBuff (where receiveda data is stored).

	int n_recvIndices;			// # recvIndices.
	int* recvIndices;			// Location of of where in recvBuff array required data is stored.
	int* recvIndices_Offsets;	// Offsets of recvIndices(where one cell starts and another begins).

	// Variables used in MPI sending / receiving.
	MPI_Status* statSend;
	MPI_Status* statRecv;
	MPI_Request* send_request;
	MPI_Request* recv_request;

	BInfo();
	~BInfo();
};

// ********************************************************************************************** //
// Holds global info about boundary cells. Only used by process 0. Only used in setup procedure.
// TTcells[i] -> a single boundary cell-number.
// [ bNumbers[offsets[i]], bNumbers[offsets[i+1]]-1  ] -> Basis numbers having TTcells[i]
// in their support.
struct H_data{
	int n_TT; 			// Total number of boundary cells in the grid.
	int TT_max; 		// Maximum number of basis functions sharing a boundary cells.

	int* TTcells;		// All cell-numbers belonging to the boundary.
	int* bNumbers; 		// Sorted basis-numbers holding each boundary cell.
	int* offsets; 		// offsets array to bNumbers (where one cell stats and another begins).

	H_data();
	~H_data();
};

// ********************************************************************************************** //
// Graph is used to store a graph in two situations:
// 1: Graph representing basis function computation and dependency. This graph will be
// partitioned to find the distribution of the problem.
// 2: Graph representing processes and their dependence. What basis function belong to what
// process, which of them will exchange data, and also what data will be exchanged.
// Below we expalin the situation 2 use.
// Only used by process 0.
struct Graph{
	int n_nodes; 			// # nodes (same as number of MPI-processes (SIZE)).
	int n_edges;			// # Edges between processes (which processes will exhange data).

	int* nodeWeights;		// sum of basis function elements in each process.
	int* edges; 			// Edges between dependent processes.
	int* edgeWeights; 		// Edge-weights corresponding to edges (# elements that will be sent ).
	int* edgeOffsets;		// Offsets to edge and edgeWeights.

    int maxRow;				// Same as in Matrix.
    int M; 					// Total number of basis functions.
    int* basisDistr; // Basis distribution array #1. Basis function i goes to process basisDistr[i].
	int* order;			// Basis distribution array #2. Basis functions [order[i],order[i+1]-1]
						// goes to process i. Used to reorder basis at the end (gatherBasis()).
    int* n_basis_Nodes; 	// # basis function in each process. (size: # processes).

	// info about cell-sums to be sent.
	int n_sendCells;	// Total number of cell-sums that must be sent.
	int* sendCells;		// Cell-numbers of sums that must be sent.
	int* sendOffsets;	// Offsets of sendCells.
	// [ sendCells[sendOffsets[i]], sendCells[sendOffsets[i+1]-1] ] are cell-numbers-sums that
	// process i will send.

	// package (size: number of processes): First problem-data each recieving process receives.
	// Contains 5 numbers to each process:
	// 1. Number of processes to exchange data with.
	// 2. Total number of elements in support regions of its distributed basis functions.
	// 3. maxRow.
	// 4. Number of basis functions it will be distributed.
	// 5. Number of cell-sums it will send.
	int* package;

	Graph();
	~Graph();
};

// ********************************************************************************************** //
// Node is created by each prorcess in the setup procedure, and holds key data about what
// info it will receive.
struct Node{
	int number; 			// Node number (RANK)
	int weight;             // Weight of process (sum of elements in its assigned basis functions).
    int n_edges;            // Number of processes it will send /recieve data to.
    int* edges;             // Process-numbers of the processes it will send / receive data to.

    int maxRow;             // Same as in Matrix.
    int n_basis;            // number of basis functions assigned to the process.
	int n_sendCells; 		// Number of cell-sums that must be sent.
	int* sendCells; 		// The numbers of the cell-sums that must be sent.

    Node();
    ~Node();
};

// ********************************************************************************************** //
// Holds options regarding the itration procedure.
struct Options{
	double tolerance;	// Stop criteria
	int maxIter;		// Maximum number of iterations.
	double omega;		// Relaxation factor.
	int checkTol;		// How often to check for convergence

	int s; 				// # of extra smoothings to non-boundary cells.

	// Variables below are used under convergence check.
	int counter;	// How many iterations since last convergence check.
	int underTol;	// True if current process has reached tolerance.
	int underTol_G;	// True if all processes has reached tolerance.
	// underTol and underTol_G are used as boolean variables, but are ints because MPI does not
	// have a bool-type.

	Options();
};
