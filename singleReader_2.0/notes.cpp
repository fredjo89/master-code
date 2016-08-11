INITIAL MEMORY USAGE:

Grid{
    int     N:          Number of cells
    int     M:          Number of blocks
    int     n_basis:    Number of cellnumbers in support

    int*    order:      NULL
    int*    support:    Length: n_basis
    int*    celltypes:  Length: n_basis
    double* basis:      Length: n_basis
    double* updates:    NULL
    int*    offsets:    Length: M+1

    TOTAL MEMORY USAGE:
    CONSTANTS:          3
    int:                2*n_basis + M+1
    double:             n_basis
}

MATRIX{
    int     N:          Number of cells
    int     maxRow:     Max number of elements in row
    int     n_basis:    Number of cellnumbers in support

    int*    j_index:    Length: maxRow*N
    double* conn:       Length: maxRow*N

    int*    loc_index:  NULL
    double* loc_conn:   NULL

    TOTAL MEMORY USAGE:
    CONSTANTS:          3
    int:                maxRow*N
    double:             maxRow*N
}

// ********************************************************************************************** //
AFTER makeCoarseGraph HAS BEEN EXECUTED:

GRAPH{
	int    n_nodes:        SIZE (Number of threads)
	int    n_edges:        Number of edges between threads 

	int*   nodeWeights:	   Length: SIZE
	int*   edgeOffsets:    Length: SIZE+1
	int*   edges:          Length: n_edges (<= SIZE*(SIZE-1))
    int*   edgeWeights:    Length: n_edges


    int     maxRow:         Max number of elements in row
    int     M:              Number of blocks
    int*    basisDistr:     Length: M
    int*    n_basis_Nodes:  Length: SIZE

	int    n_sendCells:    Length of sendCells: (<= SIZE*N)
	int*   sendOffsets:    Length: SIZE+1
	int*   sendCells:      Length: n_sendCells

	int*   package:        Length: 5*SIZE

    TOTAL MEMORY USAGE:
    CONSTANTS:          5
    int:                n_sendCells +  M + 2*SIZE(SIZE-1) + 9*SIZE + 2
};
































// yolo
