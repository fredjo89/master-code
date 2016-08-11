#include "jacobi.h"
#include "data_defs.h"
#include "utilities.h"

#include <mpi.h>
#include <cmath>
#include <algorithm>
using namespace std;

// ********************************************************************************************** //
void setupSupportMatrix(Grid& grid, Matrix& mat){
    mat.sup_index = new int[mat.maxRow * grid.n_sup];
    mat.sup_coef = new double[mat.maxRow * grid.n_sup];
    int i = 0;
    int j = 0;
    int k = grid.support[0]*mat.maxRow;
    for (int it = 0; it < grid.n_sup*mat.maxRow; it++){
        if (j==mat.maxRow){
            j = 0;
            i++;
            k = grid.support[i]*mat.maxRow;
        }
        mat.sup_index[it] = mat.mat_index[k + j];
        mat.sup_coef[it] = mat.mat_coef[k + j];
        j++;
    }
}


// ********************************************************************************************** //
// (Can be deleted)
void setupSupportMatrix2(Grid& grid, Matrix& mat){
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
    for (int i = 0; i<grid.M; i++){
        int start = grid.offsets[i];
        int end = grid.offsets[i+1];
        int max = grid.support[start];
        int min = grid.support[start];

        for (int j = start; j<end; j++ ){
            if (grid.support[j]>max){
                max = grid.support[j];
            }
            else if (grid.support[j]<min){
                min = grid.support[j];
            }
        }

        int* myMap = new int[max-min+1];
        for (int j=0; j<max-min+1; j++){
            myMap[j] = -2;
        }
        for (int j=start; j<end; j++){
            myMap[grid.support[j]-min] = j;
        }
        for (int j=start*mat.maxRow; j<end*mat.maxRow; j++){
            int index = mat.sup_index[j];
            if (index!=-1){
                if (index>=min && index<=max){
                    mat.sup_index[j] = myMap[index - min];
                }
                else{
                    mat.sup_index[j] = -2;
                }
            }
        }
        delete[] myMap;
    }
}

// ********************************************************************************************** //
void makeTypeTwo_local(Grid& grid, Node& node, TypeTwo& tt){
    int max = grid.support[0];
    int min = grid.support[0];
    for (int i = 1; i < grid.n_sup; i++ ){
        if (grid.support[i] > max){
            max = grid.support[i];
        }
        else if (grid.support[i] < min){
            min = grid.support[i];
        }
    }
    int L = max - min + 1;

    int* TTcells = new int[L]();
    int* temp1 = new int[2*grid.n_sup];
    int maxTT = 0;

    for (int i = 0; i < grid.n_sup; i++){
        if (grid.celltypes[i] == 2){
            temp1[2 * tt.n_cIndices] = grid.support[i];
            temp1[2 * tt.n_cIndices+1] = i;
            if (TTcells[grid.support[i] - min] == 0) tt.n_cNumbers++;
            TTcells[grid.support[i] - min]++;
            if (TTcells[grid.support[i] - min] > maxTT) maxTT = TTcells[grid.support[i] - min];
            tt.n_cIndices++;
        }
    }

    tt.cNumbers = new int[tt.n_cNumbers];
    tt.cIndices = new int[tt.n_cIndices];
    tt.cIndices_Offsets = new int[tt.n_cNumbers + 1];   tt.cIndices_Offsets[0] = 0;

    for (int i = 0; i < node.n_sendCells; i++) tt.cNumbers[i] = node.sendCells[i];
    sort(tt.cNumbers, tt.cNumbers + node.n_sendCells );

    int c1 = 0;
    int c2 = node.n_sendCells;
    for (int i = 0; i < L; i++){
        if (TTcells[i] != 0){
            if (i + min == tt.cNumbers[c1] && c1 < node.n_sendCells){
                tt.cIndices_Offsets[c1 + 1]= TTcells[i];
                c1++;
            }
            else{
                tt.cNumbers[c2]= i + min;
                tt.cIndices_Offsets[c2 + 1]= TTcells[i];
                c2++;
            }
        }
    }
    for (int i = 0; i < tt.n_cNumbers; i++){
        tt.cIndices_Offsets[i + 1] += tt.cIndices_Offsets[i];
    }

    int* temp2 = new int[L*maxTT]; for (int i = 0; i < L*maxTT; i++) temp2[i] = -1;
    int* temp3 = new int[L]();

    for (int i = 0; i < tt.n_cIndices; i++){
        int index = temp1[2*i]-min;
        temp2[ index*maxTT + temp3[index] ] = temp1[2*i + 1];
        temp3[index]++;
    }
    c1 = 0; c2 = 0;
    for (int i = 0; i < L; i++){
        int cellNumber = i + min;
        if (cellNumber == tt.cNumbers[c1]){
            for (int j = i*maxTT; j < (i+1)*maxTT; j++){
                if (temp2[j] != -1){
                    tt.cIndices[c2] = temp2[j];
                    c2++;
                }
                temp2[j] = -1;
            }
            c1++;
        }
    }
    for (int i = 0; i < L*maxTT; i++){
        if (temp2[i] != -1){
            tt.cIndices[c2] = temp2[i];
            c2++;
        }
    }

    tt.sums = new double[tt.n_cNumbers];
    delete[] TTcells; delete[] temp1; delete[] temp2; delete[] temp3;
    tt.n_sendCells = node.n_sendCells;
    tt.n_edges = node.n_edges;
    tt.edges = node.edges; node.edges = NULL;
    delete[] node.sendCells;    node.sendCells = NULL;
}

// ********************************************************************************************** //
void makeTypeTwo_sending(Grid& grid, TypeTwo& tt, int RANK){
    // 1. Initial declarations
    tt.statSend = new MPI_Status[tt.n_edges]();
    tt.statRecv = new MPI_Status[tt.n_edges]();
    tt.send_request = new MPI_Request[tt.n_edges]();
    tt.recv_request = new MPI_Request[tt.n_edges]();

    // 2. Create n_recvBuff, recvBuff and recvBuff_Offsets in TypeTwo.
    tt.recvBuff_Offsets = new int[tt.n_edges + 1];
    tt.recvBuff_Offsets[0] = 0;
    for (int i = 0; i<tt.n_edges; i++){
        MPI_Isend(&tt.n_sendCells, 1, MPI_INT, tt.edges[i], 666,
        MPI_COMM_WORLD, &tt.send_request[i]);
        MPI_Irecv(&tt.recvBuff_Offsets[i+1], 1, MPI_INT, tt.edges[i], 666,
        MPI_COMM_WORLD, &tt.recv_request[i]);
    }
    MPI_Waitall(tt.n_edges, tt.send_request, tt.statSend);
    MPI_Waitall(tt.n_edges, tt.recv_request, tt.statRecv);
    for (int i = 0; i<tt.n_edges; i++){
        tt.recvBuff_Offsets[i+1]+=tt.recvBuff_Offsets[i];
    }


    tt.n_recvBuff = tt.recvBuff_Offsets[tt.n_edges];
    tt.recvBuff = new double[tt.n_recvBuff];
    // 3. Create n_recvIndices, recvIndices and recvIndices_Offsets in TypeTwo.
    tt.n_recvIndices = 0;
    tt.recvIndices = new int[tt.n_recvBuff];
    tt.recvIndices_Offsets = new int[tt.n_sendCells+1];
    tt.recvIndices_Offsets[0] = 0;
    int* temp_recvBuff = new int[tt.n_recvBuff];
    for (int i = 0; i<tt.n_edges; i++){
        MPI_Isend(tt.cNumbers, tt.n_sendCells, MPI_INT, tt.edges[i],
        666,MPI_COMM_WORLD, &tt.send_request[i]);
        MPI_Irecv(&temp_recvBuff[tt.recvBuff_Offsets[i]],
        tt.recvBuff_Offsets[i+1]-tt.recvBuff_Offsets[i],
        MPI_INT, tt.edges[i], 666, MPI_COMM_WORLD, &tt.recv_request[i]);
    }
    MPI_Waitall(tt.n_edges, tt.send_request, tt.statSend);
    MPI_Waitall(tt.n_edges, tt.recv_request, tt.statRecv);
    for (int i = 0; i<tt.n_sendCells; i++){
        tt.recvIndices_Offsets[i+1] = tt.recvIndices_Offsets[i];
        for (int j = 0; j<tt.n_recvBuff; j++){
            if (temp_recvBuff[j] == tt.cNumbers[i]){
                tt.recvIndices[tt.n_recvIndices] = j;
                tt.n_recvIndices++;
                tt.recvIndices_Offsets[i+1]++;
            }
        }
    }
    delete[] temp_recvBuff;
}

// ********************************************************************************************** //
void makeTypeTwo_serial(Grid& grid, TypeTwo& tt){
    int* TTcells = new int[grid.N]();
    int* temp1 = new int[2*grid.n_sup];
    int maxTT = 0;

    for (int i = 0; i<grid.n_sup; i++){
        if (grid.celltypes[i]==2){
            temp1[2*tt.n_cIndices] = grid.support[i];
            temp1[2*tt.n_cIndices+1] = i;
            if (TTcells[grid.support[i]]==0) tt.n_cNumbers++;
            TTcells[grid.support[i]]++;
            if (TTcells[grid.support[i]]>maxTT) maxTT = TTcells[grid.support[i]];
            tt.n_cIndices++;
        }
    }

    tt.cNumbers = new int[tt.n_cNumbers];
    tt.cIndices = new int[tt.n_cIndices];
    tt.cIndices_Offsets = new int[tt.n_cNumbers+1];   tt.cIndices_Offsets[0] = 0;

    int k = 0;
    for (int i = 0; i<grid.N; i++){
        if (TTcells[i]!=0){
            tt.cNumbers[k] = i;
            tt.cIndices_Offsets[k+1] = tt.cIndices_Offsets[k] + TTcells[i];
            k++;
        }
    }

    int* temp2 = new int[grid.N*maxTT];
    int* temp3 = new int[grid.N]();
    for (int i=0; i<grid.N*maxTT; i++) temp2[i] = -1;
    for (int i=0; i<tt.n_cIndices; i++){
        int index = temp1[2*i];
        temp2[maxTT*index+temp3[index]] = temp1[2*i+1];
        temp3[index]++;
    }

    k = 0;
    for (int i = 0; i<grid.N*maxTT; i++){
        if (temp2[i]!=-1){
            tt.cIndices[k] = temp2[i];
            k++;
        }
    }

    tt.sums = new double[tt.n_cNumbers];
    delete[] TTcells; delete[] temp1; delete[] temp2; delete[] temp3;
}

// ********************************************************************************************** //
void jacobi(Grid& grid, Matrix& mat, double& omega){

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

    #pragma omp parallel for schedule(static)
    for (int i=0; i<grid.n_sup; i++){
        grid.basis[i]-=omega*grid.updates[i];
    }
}

// ********************************************************************************************** //
void localSum(Grid& grid, TypeTwo& tt, int RANK){

    #pragma omp parallel for schedule(static)
    for (int i = 0; i<tt.n_cNumbers; i++) tt.sums[i] = 0;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i<tt.n_cNumbers; i++){
        for (int j =tt.cIndices_Offsets[i]; j<tt.cIndices_Offsets[i+1]; j++ ){
            tt.sums[i]+=grid.basis[tt.cIndices[j]];
        }
    }
}

// ********************************************************************************************** //
void sendAndRecieve(Grid& grid, TypeTwo& tt, int RANK){
    for (int i = 0; i<tt.n_edges; i++){
        MPI_Isend(tt.sums, tt.n_sendCells, MPI_DOUBLE,
        tt.edges[i], 666,MPI_COMM_WORLD, &tt.send_request[i]);
        MPI_Irecv(&tt.recvBuff[tt.recvBuff_Offsets[i]],
        tt.recvBuff_Offsets[i+1]-tt.recvBuff_Offsets[i], MPI_DOUBLE,
        tt.edges[i], 666, MPI_COMM_WORLD, &tt.recv_request[i]);
    }
    MPI_Waitall(tt.n_edges, tt.send_request, tt.statSend);
    MPI_Waitall(tt.n_edges, tt.recv_request, tt.statRecv);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i<tt.n_sendCells; i++){
        for (int j=tt.recvIndices_Offsets[i]; j<tt.recvIndices_Offsets[i+1]; j++ ){
            tt.sums[i] += tt.recvBuff[tt.recvIndices[j]];
        }
    }
}

// ********************************************************************************************** //
void TTnormalize(Grid& grid, TypeTwo& tt, int RANK){

    #pragma omp parallel for schedule(static)
    for (int i = 0; i<tt.n_cNumbers; i++){
        for (int j =tt.cIndices_Offsets[i]; j<tt.cIndices_Offsets[i+1]; j++ ){
            grid.basis[tt.cIndices[j]]/=tt.sums[i];
        }
    }
}
















//yolo
