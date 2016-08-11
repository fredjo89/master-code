#pragma once
#include <iostream>
using namespace std;

// ********************************************************************************************** //
struct Grid{
	int N;				// number of fine cells
	int M;				// number of coarse blocks
	int n_basis;

	int* order; 		// Order of basis functions

	int* support; 		// n_basis elements
	int* celltypes; 	// n_basis elements
	double* basis;		// n_basis elements

	int* offsets; 		// M + 1 elements

	Grid();
	~Grid();
	void print();
	void printOrder();
};

// ********************************************************************************************** //
struct Matrix{
	int N;
	int maxRow;
	int n_basis;

	double*	conn;
	int*	j_index;

	int*  	loc_index;
 	double*	loc_conn;

	Matrix();
	~Matrix();
	void print();
};

// ********************************************************************************************** //
struct H_data{
	int n_TT; 			// Number of type 2 cells

	int* TTcells;		// Array of type 2 cell numbers

	int* offsets; 		// offsets array to bNumbers
	int* bNumbers; 		// array of basisNumbers

	H_data();
	~H_data();
	void print();
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
    int* n_basis_Nodes;         // Number of basis functions on each node.

	// info about send cells
	int n_sendCells;
	int* sendOffsets;
	int* sendCells;

	int* package;               // the package that will be sent in the first round

	Graph();
	~Graph();
	void print();
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
    void print();
};

// ********************************************************************************************** //
struct Options{
	double tolerance;
	int maxIter;
	double omega;
	int reNorm;
	int checkTol;

	Options();
	void print();
};
















// yolo
