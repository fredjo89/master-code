#pragma once
#include "data_defs.h"

#include <iostream>
#include <mpi.h>
using namespace std;

// ********************************************************************************************** //
struct TT{

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


	TT(){
		n_cNumbers = -1; 	n_cIndices = -1; 	n_sendCells = -1; 	n_edges = -1;
		n_recvBuff = -1;	n_recvIndices = -1;
		cNumbers = NULL; 	sums = NULL; 		cIndices = NULL;		cIndices_Offsets = NULL;
		edges = NULL; 		recvBuff = NULL;
		recvIndices = NULL; 	recvIndices_Offsets = NULL;
		statSend = NULL; 	statRecv = NULL; 	send_request = NULL;	recv_request = NULL;
	}

	~TT(){
		delete[] cNumbers;			delete[] sums; 					delete[] cIndices;
		delete[] cIndices_Offsets;	delete[] edges; 				delete[] recvBuff;
		delete[] recvIndices;		delete[] recvIndices_Offsets;
		delete[] statSend; 			delete[] statRecv;
		delete[] send_request; 		delete[] recv_request;
	}

	void print();

};

void makeTT(Grid& grid, Node& node, TT& tt);

void localSum(Grid& grid, TT& tt, int RANK);

void sendAndRecieve(Grid& grid, TT& tt, int RANK);

void TTnormalize(Grid& grid, TT& tt, int RANK);

void gatherBasis(Grid& localGrid, Grid& globalGrid, Graph& graph, int RANK, int SIZE);


void setupSupportMatrix(Grid& grid, Matrix& mat);

void makeMapping(Grid& grid, Matrix& mat);

void jacobi(Grid& grid, Matrix& mat, double omega);

void computeDiscrepancy(Grid& grid, double* ompSol);
