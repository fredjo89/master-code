#include "jacobi.h"
#include "data_defs.h"
#include "utilities.h"

#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <algorithm>
using namespace std;

// ********************************************************************************************** //
void setupSupportMatrix(Grid& grid, Matrix& mat){
    int n_sup = grid.n_sup;
    int maxRow = mat.maxRow;
    mat.sup_index = new int[maxRow * n_sup];
    mat.sup_coef = new double[maxRow * n_sup];
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_sup; i++){
        int index = grid.support[i];
        for (int j = 0; j < mat.maxRow; j++){
            mat.sup_index[i*maxRow + j] = mat.mat_index[index*maxRow + j];
            mat.sup_coef[i*maxRow + j] = mat.mat_coef[index*maxRow + j];
        }
    }
}

// ********************************************************************************************** //
void makeMapping(Grid& grid, Matrix& mat){
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < grid.M; i++){
        int start = grid.offsets[i];
        int end = grid.offsets[i+1];
        int max = grid.support[start];
        int min = grid.support[start];

        for (int j = start; j < end; j++ ){
            if (grid.support[j] > max){
                max = grid.support[j];
            }
            else if (grid.support[j] < min){
                min = grid.support[j];
            }
        }

        int* myMap = new int[max-min+1];
        for (int j = 0; j < max-min+1; j++){
            myMap[j] = -2;
        }

        for (int j = start; j < end; j++) {
            myMap[grid.support[j]-min] = j;
        }

        for (int j = start*mat.maxRow; j < end*mat.maxRow; j++){
            int index = mat.sup_index[j];
            if (index!=-1){
                if (index >= min && index <= max){
                    mat.sup_index[j] = myMap[index - min];
                }
                else {
                    mat.sup_index[j] = -2;
                }
            }
        }
        delete[] myMap;
    }
}

// ********************************************************************************************** //
void makeBInfo_local(Grid& grid, Node& node, BInfo& B){
    if (node.weight == 0) return;   // If process is not assinged any work, return.

    int max = grid.support[0];
    int min = grid.support[0];
    for (int i = 1; i < grid.n_sup; i++ ){
        if (grid.support[i] > max)      max = grid.support[i];
        else if (grid.support[i] < min) min = grid.support[i];
    }
    int L = max - min + 1;

    int* TTcells = new int[L]();
    int* temp1 = new int[2*grid.n_sup];
    int maxTT = 0;

    for (int i = 0; i < grid.n_sup; i++){
        if ( grid.boundary[i] ){
            temp1[2 * B.n_cIndices] = grid.support[i];
            temp1[2 * B.n_cIndices+1] = i;
            if (TTcells[grid.support[i] - min] == 0) B.n_cNumbers++;
            TTcells[grid.support[i] - min]++;
            if (TTcells[grid.support[i] - min] > maxTT) maxTT = TTcells[grid.support[i] - min];
            B.n_cIndices++;
        }
    }

    B.cNumbers = new int[B.n_cNumbers];
    B.cIndices = new int[B.n_cIndices];
    B.cIndices_Offsets = new int[B.n_cNumbers + 1];   B.cIndices_Offsets[0] = 0;

    for (int i = 0; i < node.n_sendCells; i++) B.cNumbers[i] = node.sendCells[i];
    sort(B.cNumbers, B.cNumbers + node.n_sendCells );

    int c1 = 0;
    int c2 = node.n_sendCells;
    for (int i = 0; i < L; i++){
        if (TTcells[i] != 0){
            if (i + min == B.cNumbers[c1] && c1 < node.n_sendCells){
                B.cIndices_Offsets[c1 + 1]= TTcells[i];
                c1++;
            }
            else{
                B.cNumbers[c2]= i + min;
                B.cIndices_Offsets[c2 + 1]= TTcells[i];
                c2++;
            }
        }
    }
    for (int i = 0; i < B.n_cNumbers; i++){
        B.cIndices_Offsets[i + 1] += B.cIndices_Offsets[i];
    }

    int* temp2 = new int[L*maxTT]; for (int i = 0; i < L*maxTT; i++) temp2[i] = -1;
    int* temp3 = new int[L]();

    for (int i = 0; i < B.n_cIndices; i++){
        int index = temp1[2*i]-min;
        temp2[ index*maxTT + temp3[index] ] = temp1[2*i + 1];
        temp3[index]++;
    }
    c1 = 0; c2 = 0;
    for (int i = 0; i < L; i++){
        int cellNumber = i + min;
        if (cellNumber == B.cNumbers[c1]){
            for (int j = i*maxTT; j < (i+1)*maxTT; j++){
                if (temp2[j] != -1){
                    B.cIndices[c2] = temp2[j];
                    c2++;
                }
                temp2[j] = -1;
            }
            c1++;
        }
    }
    for (int i = 0; i < L*maxTT; i++){
        if (temp2[i] != -1){
            B.cIndices[c2] = temp2[i];
            c2++;
        }
    }

    B.sums = new double[B.n_cNumbers];
    delete[] TTcells; delete[] temp1; delete[] temp2; delete[] temp3;
    B.n_sendCells = node.n_sendCells;
    B.n_edges = node.n_edges;
    B.edges = node.edges; node.edges = NULL;
    delete[] node.sendCells;    node.sendCells = NULL;
}

// ********************************************************************************************** //
void makeBInfo_sending(Grid& grid, BInfo& B, int RANK){
    if (grid.n_sup == 0) return;

    // 1. Initial declarations
    B.statSend = new MPI_Status[B.n_edges]();
    B.statRecv = new MPI_Status[B.n_edges]();
    B.send_request = new MPI_Request[B.n_edges]();
    B.recv_request = new MPI_Request[B.n_edges]();

    // 2. Create n_recvBuff, recvBuff and recvBuff_Offsets in BInfo.
    B.recvBuff_Offsets = new int[B.n_edges + 1];
    B.recvBuff_Offsets[0] = 0;
    for (int i = 0; i<B.n_edges; i++){
        MPI_Isend(&B.n_sendCells, 1, MPI_INT, B.edges[i], 666,
        MPI_COMM_WORLD, &B.send_request[i]);
        MPI_Irecv(&B.recvBuff_Offsets[i+1], 1, MPI_INT, B.edges[i], 666,
        MPI_COMM_WORLD, &B.recv_request[i]);
    }
    MPI_Waitall(B.n_edges, B.send_request, B.statSend);
    MPI_Waitall(B.n_edges, B.recv_request, B.statRecv);
    for (int i = 0; i<B.n_edges; i++){
        B.recvBuff_Offsets[i+1]+=B.recvBuff_Offsets[i];
    }

    B.n_recvBuff = B.recvBuff_Offsets[B.n_edges];
    B.recvBuff = new double[B.n_recvBuff];

    // 3. Create n_recvIndices, recvIndices and recvIndices_Offsets in BInfo.
    B.n_recvIndices = 0;
    B.recvIndices = new int[B.n_recvBuff];
    B.recvIndices_Offsets = new int[B.n_sendCells+1];
    B.recvIndices_Offsets[0] = 0;
    int* temp_recvBuff = new int[B.n_recvBuff];
    for (int i = 0; i<B.n_edges; i++){
        MPI_Isend(B.cNumbers, B.n_sendCells, MPI_INT, B.edges[i],
        666,MPI_COMM_WORLD, &B.send_request[i]);
        MPI_Irecv(&temp_recvBuff[B.recvBuff_Offsets[i]],
        B.recvBuff_Offsets[i+1]-B.recvBuff_Offsets[i],
        MPI_INT, B.edges[i], 666, MPI_COMM_WORLD, &B.recv_request[i]);
    }
    MPI_Waitall(B.n_edges, B.send_request, B.statSend);
    MPI_Waitall(B.n_edges, B.recv_request, B.statRecv);
    for (int i = 0; i<B.n_sendCells; i++){
        B.recvIndices_Offsets[i+1] = B.recvIndices_Offsets[i];
        for (int j = 0; j<B.n_recvBuff; j++){
            if (temp_recvBuff[j] == B.cNumbers[i]){
                B.recvIndices[B.n_recvIndices] = j;
                B.n_recvIndices++;
                B.recvIndices_Offsets[i+1]++;
            }
        }
    }
    delete[] temp_recvBuff;
}

// ********************************************************************************************** //
void makeBInfo_serial(Grid& grid, BInfo& B){
    int* TTcells = new int[grid.N]();
    int* temp1 = new int[2*grid.n_sup];
    int maxTT = 0;

    for (int i = 0; i<grid.n_sup; i++){
        if ( grid.boundary[i] ){
            temp1[2*B.n_cIndices] = grid.support[i];
            temp1[2*B.n_cIndices+1] = i;
            if (TTcells[grid.support[i]]==0) B.n_cNumbers++;
            TTcells[grid.support[i]]++;
            if (TTcells[grid.support[i]]>maxTT) maxTT = TTcells[grid.support[i]];
            B.n_cIndices++;
        }
    }

    B.cNumbers = new int[B.n_cNumbers];
    B.cIndices = new int[B.n_cIndices];
    B.cIndices_Offsets = new int[B.n_cNumbers+1];   B.cIndices_Offsets[0] = 0;

    int k = 0;
    for (int i = 0; i<grid.N; i++){
        if (TTcells[i]!=0){
            B.cNumbers[k] = i;
            B.cIndices_Offsets[k+1] = B.cIndices_Offsets[k] + TTcells[i];
            k++;
        }
    }

    int* temp2 = new int[grid.N*maxTT];
    int* temp3 = new int[grid.N]();
    for (int i=0; i<grid.N*maxTT; i++) temp2[i] = -1;
    for (int i=0; i<B.n_cIndices; i++){
        int index = temp1[2*i];
        temp2[maxTT*index+temp3[index]] = temp1[2*i+1];
        temp3[index]++;
    }

    k = 0;
    for (int i = 0; i<grid.N*maxTT; i++){
        if (temp2[i]!=-1){
            B.cIndices[k] = temp2[i];
            k++;
        }
    }

    B.sums = new double[B.n_cNumbers];
    delete[] TTcells; delete[] temp1; delete[] temp2; delete[] temp3;
}

// ********************************************************************************************** //
void jacobi(Grid& grid, Matrix& mat, Options& opt){
    opt.counter++;

    // Extra smoothing to non-boundary cells.
    for (int k = 0; k < opt.s; k++){
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < grid.n_sup; i++){
            if (!grid.boundary[i]){
                grid.updates[i] = grid.basis[i];
                int it = i*mat.maxRow;
                for (int j = 0; j < mat.maxRow; j++){
                    int index = mat.sup_index[it+j];
                    if (index==-1) break;
                    if (index!=-2){
                        grid.updates[i] += mat.sup_coef[it+j]*grid.basis[index];
                    }
                }
            }
        }
        #pragma omp parallel for schedule(static)
        for (int i=0; i<grid.n_sup; i++){
        if (!grid.boundary[i])
            grid.basis[i]-=opt.omega*grid.updates[i];
        }
    }

    // Jacobi sweep
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < grid.n_sup; i++){
        grid.updates[i] = grid.basis[i];
        int it = i*mat.maxRow;
        for (int j = 0; j < mat.maxRow; j++){
            int index = mat.sup_index[it+j];
            if (index==-1) break;
            if (index!=-2){
                grid.updates[i] += mat.sup_coef[it+j]*grid.basis[index];
            }
        }
    }

    // Convergence check.
    if (opt.tolerance < 0 || opt.counter < opt.checkTol){
        #pragma omp parallel for schedule(static)
        for (int i=0; i<grid.n_sup; i++){
            grid.basis[i]-=opt.omega*grid.updates[i];
        }
    }
    else{
        opt.counter = 0;
        opt.underTol = true;
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < grid.n_sup; i++){
            if ( opt.underTol && !grid.boundary[i] && opt.tolerance < abs(grid.updates[i]) ){
                opt.underTol = false;
            }
            grid.basis[i] -= opt.omega*grid.updates[i];
        }
    }
}

// ********************************************************************************************** //
void localSum(Grid& grid, BInfo& B){
    #pragma omp parallel for schedule(static)
    for (int i = 0; i<B.n_cNumbers; i++){
        B.sums[i] = 0;
        for (int j =B.cIndices_Offsets[i]; j<B.cIndices_Offsets[i+1]; j++ ){
            B.sums[i]+=grid.basis[B.cIndices[j]];
        }
    }
}

// ********************************************************************************************** //
void sendAndRecieve(BInfo& B){
    for (int i = 0; i < B.n_edges; i++){
        MPI_Isend(
            B.sums, B.n_sendCells, MPI_DOUBLE,
            B.edges[i], 100, MPI_COMM_WORLD, &B.send_request[i]
        );
        MPI_Irecv(
            &B.recvBuff[ B.recvBuff_Offsets[i] ],
            B.recvBuff_Offsets[i+1] - B.recvBuff_Offsets[i],
            MPI_DOUBLE, B.edges[i], 100, MPI_COMM_WORLD, &B.recv_request[i]
        );
    }
    MPI_Waitall(B.n_edges, B.send_request, B.statSend);
    MPI_Waitall(B.n_edges, B.recv_request, B.statRecv);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < B.n_sendCells; i++){
        for (int j = B.recvIndices_Offsets[i]; j < B.recvIndices_Offsets[i+1]; j++ ){
            B.sums[i] += B.recvBuff[ B.recvIndices[j] ];
        }
    }
}

// ********************************************************************************************** //
void BNormalize(Grid& grid, BInfo& B, Options& opt){
    if (!opt.underTol || opt.counter != 0  ){
        #pragma omp parallel for schedule(static)
        for (int i = 0; i<B.n_cNumbers; i++){
            for (int j =B.cIndices_Offsets[i]; j<B.cIndices_Offsets[i+1]; j++ ){
                grid.basis[B.cIndices[j]]/=B.sums[i];
            }
        }
    }
    else {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i<B.n_cNumbers; i++){
            double prev;
            for (int j =B.cIndices_Offsets[i]; j<B.cIndices_Offsets[i+1]; j++ ){

                if (opt.underTol){
                    //prev = grid.basis[B.cIndices[j]] + opt.omega*grid.updates[B.cIndices[j]];
                }

                grid.basis[B.cIndices[j]]/=B.sums[i];

                if ( opt.underTol && abs(prev - grid.basis[B.cIndices[j]] ) > opt.tolerance ){
                    //opt.underTol = false;
                }
            }
        }
    }
}
