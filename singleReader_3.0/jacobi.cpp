#include "jacobi.h"
#include "data_defs.h"

#include <mpi.h>
#include <cmath>
#include <algorithm>
using namespace std;

// ********************************************************************************************** //
void jacobi(Grid& grid, Matrix& mat, double& omega){
    for (int i = 0; i < mat.n_basis; i++){
        grid.updates[i] = grid.basis[i];
        int it = i*mat.maxRow;
        for (int j = 0; j < mat.maxRow; j++){
            int index = mat.loc_index[it+j];
            if (index==-1) break;
            if (index!=-2){
                grid.updates[i] += mat.loc_conn[it+j]*grid.basis[index];
            }
        }
    }

    for (int i=0; i<mat.n_basis; i++){
        grid.basis[i]-=omega*grid.updates[i];
    }
}

// ********************************************************************************************** //
void makeTypeTwo(Grid& grid, Node& node, int& SIZE, TypeTwo& typeTwo){
    // CREATING LOCAL STUFF ********************************************************************* //
    // 1. Finding all TypeTwo cells and where they are placed and create n_cIndices in TypeTwo
    int* temp1 = new int[grid.n_basis];
    int* temp2 = new int[grid.n_basis];
    for (int i = 0; i<grid.n_basis; i++){
        if (grid.celltypes[i]==2){
            temp1[typeTwo.n_cIndices] = grid.support[i];
            temp2[typeTwo.n_cIndices] = i;
            typeTwo.n_cIndices++;
        }
    }
    typeTwo.cNumbers = new int[typeTwo.n_cIndices];
    typeTwo.cIndices_Offsets = new int[typeTwo.n_cIndices+1]();
    typeTwo.cIndices = new int[typeTwo.n_cIndices];
    // 2. Create cNumbers, sums cIndices_Offsets and n_cNumbers in TypeTwo
    for (int i = 0; i<node.n_sendCells; i++){
        typeTwo.cNumbers[i] = node.sendCells[i];
        typeTwo.n_cNumbers++;
    }
    int* temp3 = new int[typeTwo.n_cIndices];
    for (int i=0; i<typeTwo.n_cIndices; i++){
        bool isIn = false;
        for (int j=0; j<typeTwo.n_cNumbers; j++){
            if (typeTwo.cNumbers[j]==temp1[i]){
                temp3[i] = j;
                typeTwo.cIndices_Offsets[j+1]++;
                isIn = true;
                break;
            }
        }
        if (!isIn){
            temp3[i] = typeTwo.n_cNumbers;
            typeTwo.cNumbers[typeTwo.n_cNumbers]=temp1[i];
            typeTwo.cIndices_Offsets[typeTwo.n_cNumbers+1]++;
            typeTwo.n_cNumbers++;
        }
    }
    typeTwo.sums = new double[typeTwo.n_cNumbers];
    for (int i = 0; i<typeTwo.n_cNumbers; i++){
        typeTwo.cIndices_Offsets[i+1]+=typeTwo.cIndices_Offsets[i];
    }
    // 3. Create cIndices in TypeTwo
    int* temp4 = new int[typeTwo.n_cNumbers]();
    for (int i = 0; i<typeTwo.n_cIndices; i++){
        int cIndex = temp3[i];
        typeTwo.cIndices[ typeTwo.cIndices_Offsets[cIndex] + temp4[cIndex] ] = temp2[i];
        temp4[cIndex]++;
    }
    delete[] temp1; delete[] temp2; delete[] temp3; delete[] temp4;
    if (SIZE==1) return;
    // CREATING SENDING AND RECIEVING STUFF ***************************************************** //
    // 1. Initial declarations
    typeTwo.statSend = new MPI_Status[node.n_edges];
    typeTwo.statRecv = new MPI_Status[node.n_edges];
    typeTwo.send_request = new MPI_Request[node.n_edges];
    typeTwo.recv_request = new MPI_Request[node.n_edges];
    // 2. Create n_sendCells, n_edges and edges in TypeTwo. Delete arrays in node.
    typeTwo.n_sendCells = node.n_sendCells;
    typeTwo.n_edges = node.n_edges;
    typeTwo.edges = node.edges;
    node.edges = NULL;
    delete[] node.sendCells;
    node.sendCells = NULL;
    // 3. Create n_recvBuff, recvBuff and recvBuff_Offsets in TypeTwo.
    typeTwo.recvBuff_Offsets = new int[typeTwo.n_edges + 1];
    typeTwo.recvBuff_Offsets[0] = 0;
    for (int i = 0; i<typeTwo.n_edges; i++){
        MPI_Isend(&typeTwo.n_sendCells, 1, MPI_INT, typeTwo.edges[i], 666,
        MPI_COMM_WORLD, &typeTwo.send_request[i]);
        MPI_Irecv(&typeTwo.recvBuff_Offsets[i+1], 1, MPI_INT, typeTwo.edges[i], 666,
        MPI_COMM_WORLD, &typeTwo.recv_request[i]);
    }
    MPI_Waitall(typeTwo.n_edges, typeTwo.send_request, typeTwo.statSend);
    MPI_Waitall(typeTwo.n_edges, typeTwo.recv_request, typeTwo.statRecv);
    for (int i = 0; i<typeTwo.n_edges; i++){
        typeTwo.recvBuff_Offsets[i+1]+=typeTwo.recvBuff_Offsets[i];
    }
    typeTwo.n_recvBuff = typeTwo.recvBuff_Offsets[typeTwo.n_edges];
    typeTwo.recvBuff = new double[typeTwo.n_recvBuff];
    // 4. Create n_recvIndices, recvIndices and recvIndices_Offsets in TypeTwo.
    typeTwo.n_recvIndices = 0;
    typeTwo.recvIndices = new int[typeTwo.n_recvBuff];
    typeTwo.recvIndices_Offsets = new int[typeTwo.n_sendCells+1];
    typeTwo.recvIndices_Offsets[0] = 0;
    int* temp_recvBuff = new int[typeTwo.n_recvBuff];
    for (int i = 0; i<typeTwo.n_edges; i++){
        MPI_Isend(typeTwo.cNumbers, typeTwo.n_sendCells, MPI_INT, typeTwo.edges[i],
        666,MPI_COMM_WORLD, &typeTwo.send_request[i]);
        MPI_Irecv(&temp_recvBuff[typeTwo.recvBuff_Offsets[i]],
        typeTwo.recvBuff_Offsets[i+1]-typeTwo.recvBuff_Offsets[i],
        MPI_INT, typeTwo.edges[i], 666, MPI_COMM_WORLD, &typeTwo.recv_request[i]);
    }
    MPI_Waitall(typeTwo.n_edges, typeTwo.send_request, typeTwo.statSend);
    MPI_Waitall(typeTwo.n_edges, typeTwo.recv_request, typeTwo.statRecv);
    for (int i = 0; i<typeTwo.n_sendCells; i++){
        typeTwo.recvIndices_Offsets[i+1] = typeTwo.recvIndices_Offsets[i];
        for (int j = 0; j<typeTwo.n_recvBuff; j++){
            if (temp_recvBuff[j] == typeTwo.cNumbers[i]){
                typeTwo.recvIndices[typeTwo.n_recvIndices] = j;
                typeTwo.n_recvIndices++;
                typeTwo.recvIndices_Offsets[i+1]++;
            }
        }
    }
    delete[] temp_recvBuff;
}




// ********************************************************************************************** //
void makeTypeTwo_local(Grid& grid, Node& node, TypeTwo& typeTwo){
    int max = grid.support[0];
    int min = grid.support[0];
    for (int i = 1; i < grid.n_basis; i++ ){
        if (grid.support[i] > max)      max = grid.support[i];
        else if (grid.support[i] < min) min = grid.support[i];
    }
    int L = max - min + 1;

    int* TTcells = new int[L]();
    int* temp1 = new int[2*grid.n_basis];
    int maxTT = 0;

    for (int i = 0; i < grid.n_basis; i++){
        if (grid.celltypes[i] == 2){
            temp1[2 * typeTwo.n_cIndices] = grid.support[i];
            temp1[2 * typeTwo.n_cIndices+1] = i;
            if (TTcells[grid.support[i] - min] == 0) typeTwo.n_cNumbers++;
            TTcells[grid.support[i] - min]++;
            if (TTcells[grid.support[i] - min] > maxTT) maxTT = TTcells[grid.support[i] - min];
            typeTwo.n_cIndices++;
        }
    }

    typeTwo.cNumbers = new int[typeTwo.n_cNumbers];
    typeTwo.cIndices = new int[typeTwo.n_cIndices];
    typeTwo.cIndices_Offsets = new int[typeTwo.n_cNumbers + 1];   typeTwo.cIndices_Offsets[0] = 0;

    for (int i = 0; i < node.n_sendCells; i++) typeTwo.cNumbers[i] = node.sendCells[i];
    sort(typeTwo.cNumbers, typeTwo.cNumbers + node.n_sendCells );

    int c1 = 0;
    int c2 = node.n_sendCells;
    for (int i = 0; i < L; i++){
        if (TTcells[i] != 0){
            if (i + min == typeTwo.cNumbers[c1] && c1 < node.n_sendCells){
                typeTwo.cIndices_Offsets[c1 + 1]= TTcells[i];
                c1++;
            }
            else{
                typeTwo.cNumbers[c2]= i + min;
                typeTwo.cIndices_Offsets[c2 + 1]= TTcells[i];
                c2++;
            }
        }
    }
    for (int i = 0; i < typeTwo.n_cNumbers; i++){
        typeTwo.cIndices_Offsets[i + 1] += typeTwo.cIndices_Offsets[i];
    }

    int* temp2 = new int[L*maxTT]; for (int i = 0; i < L*maxTT; i++) temp2[i] = -1;
    int* temp3 = new int[L]();

    for (int i = 0; i < typeTwo.n_cIndices; i++){
        int index = temp1[2*i]-min;
        temp2[ index*maxTT + temp3[index] ] = temp1[2*i + 1];
        temp3[index]++;
    }
    c1 = 0; c2 = 0;
    for (int i = 0; i < L; i++){
        int cellNumber = i + min;
        if (cellNumber == typeTwo.cNumbers[c1]){
            for (int j = i*maxTT; j < (i+1)*maxTT; j++){
                if (temp2[j] != -1){
                    typeTwo.cIndices[c2] = temp2[j];
                    c2++;
                }
                temp2[j] = -1;
            }
            c1++;
        }
    }
    for (int i = 0; i < L*maxTT; i++){
        if (temp2[i] != -1){
            typeTwo.cIndices[c2] = temp2[i];
            c2++;
        }
    }

    typeTwo.sums = new double[typeTwo.n_cNumbers];
    delete[] TTcells; delete[] temp1; delete[] temp2; delete[] temp3;
    typeTwo.n_sendCells = node.n_sendCells;
    typeTwo.n_edges = node.n_edges;
    typeTwo.edges = node.edges; node.edges = NULL;
    delete[] node.sendCells;    node.sendCells = NULL;
}

// ********************************************************************************************** //
void makeTypeTwo_sending(Grid& grid, TypeTwo& typeTwo){
    // 1. Initial declarations
    typeTwo.statSend = new MPI_Status[typeTwo.n_edges];
    typeTwo.statRecv = new MPI_Status[typeTwo.n_edges];
    typeTwo.send_request = new MPI_Request[typeTwo.n_edges];
    typeTwo.recv_request = new MPI_Request[typeTwo.n_edges];
    // 2. Create n_recvBuff, recvBuff and recvBuff_Offsets in TypeTwo.
    typeTwo.recvBuff_Offsets = new int[typeTwo.n_edges + 1];
    typeTwo.recvBuff_Offsets[0] = 0;
    for (int i = 0; i<typeTwo.n_edges; i++){
        MPI_Isend(&typeTwo.n_sendCells, 1, MPI_INT, typeTwo.edges[i], 666,
        MPI_COMM_WORLD, &typeTwo.send_request[i]);
        MPI_Irecv(&typeTwo.recvBuff_Offsets[i+1], 1, MPI_INT, typeTwo.edges[i], 666,
        MPI_COMM_WORLD, &typeTwo.recv_request[i]);
    }
    MPI_Waitall(typeTwo.n_edges, typeTwo.send_request, typeTwo.statSend);
    MPI_Waitall(typeTwo.n_edges, typeTwo.recv_request, typeTwo.statRecv);
    for (int i = 0; i<typeTwo.n_edges; i++){
        typeTwo.recvBuff_Offsets[i+1]+=typeTwo.recvBuff_Offsets[i];
    }
    typeTwo.n_recvBuff = typeTwo.recvBuff_Offsets[typeTwo.n_edges];
    typeTwo.recvBuff = new double[typeTwo.n_recvBuff];
    // 3. Create n_recvIndices, recvIndices and recvIndices_Offsets in TypeTwo.
    typeTwo.n_recvIndices = 0;
    typeTwo.recvIndices = new int[typeTwo.n_recvBuff];
    typeTwo.recvIndices_Offsets = new int[typeTwo.n_sendCells+1];
    typeTwo.recvIndices_Offsets[0] = 0;
    int* temp_recvBuff = new int[typeTwo.n_recvBuff];
    for (int i = 0; i<typeTwo.n_edges; i++){
        MPI_Isend(typeTwo.cNumbers, typeTwo.n_sendCells, MPI_INT, typeTwo.edges[i],
        666,MPI_COMM_WORLD, &typeTwo.send_request[i]);
        MPI_Irecv(&temp_recvBuff[typeTwo.recvBuff_Offsets[i]],
        typeTwo.recvBuff_Offsets[i+1]-typeTwo.recvBuff_Offsets[i],
        MPI_INT, typeTwo.edges[i], 666, MPI_COMM_WORLD, &typeTwo.recv_request[i]);
    }
    MPI_Waitall(typeTwo.n_edges, typeTwo.send_request, typeTwo.statSend);
    MPI_Waitall(typeTwo.n_edges, typeTwo.recv_request, typeTwo.statRecv);
    for (int i = 0; i<typeTwo.n_sendCells; i++){
        typeTwo.recvIndices_Offsets[i+1] = typeTwo.recvIndices_Offsets[i];
        for (int j = 0; j<typeTwo.n_recvBuff; j++){
            if (temp_recvBuff[j] == typeTwo.cNumbers[i]){
                typeTwo.recvIndices[typeTwo.n_recvIndices] = j;
                typeTwo.n_recvIndices++;
                typeTwo.recvIndices_Offsets[i+1]++;
            }
        }
    }
    delete[] temp_recvBuff;
}


// ********************************************************************************************** //
void makeTypeTwo_serial(Grid& grid, TypeTwo& typeTwo){
    int* TTcells = new int[grid.N]();
    int* temp1 = new int[2*grid.n_basis];
    int maxTT = 0;

    for (int i = 0; i<grid.n_basis; i++){
        if (grid.celltypes[i]==2){
            temp1[2*typeTwo.n_cIndices] = grid.support[i];
            temp1[2*typeTwo.n_cIndices+1] = i;
            if (TTcells[grid.support[i]]==0) typeTwo.n_cNumbers++;
            TTcells[grid.support[i]]++;
            if (TTcells[grid.support[i]]>maxTT) maxTT = TTcells[grid.support[i]];
            typeTwo.n_cIndices++;
        }
    }

    typeTwo.cNumbers = new int[typeTwo.n_cNumbers];
    typeTwo.cIndices = new int[typeTwo.n_cIndices];
    typeTwo.cIndices_Offsets = new int[typeTwo.n_cNumbers+1];   typeTwo.cIndices_Offsets[0] = 0;

    int k = 0;
    for (int i = 0; i<grid.N; i++){
        if (TTcells[i]!=0){
            typeTwo.cNumbers[k] = i;
            typeTwo.cIndices_Offsets[k+1] = typeTwo.cIndices_Offsets[k] + TTcells[i];
            k++;
        }
    }

    int* temp2 = new int[grid.N*maxTT];
    int* temp3 = new int[grid.N]();
    for (int i=0; i<grid.N*maxTT; i++) temp2[i] = -1;
    for (int i=0; i<typeTwo.n_cIndices; i++){
        int index = temp1[2*i];
        temp2[maxTT*index+temp3[index]] = temp1[2*i+1];
        temp3[index]++;
    }

    k = 0;
    for (int i = 0; i<grid.N*maxTT; i++){
        if (temp2[i]!=-1){
            typeTwo.cIndices[k] = temp2[i];
            k++;
        }
    }

    typeTwo.sums = new double[typeTwo.n_cNumbers];
    delete[] TTcells; delete[] temp1; delete[] temp2; delete[] temp3;
}

// ********************************************************************************************** //
void localSum(Grid& grid, TypeTwo& typeTwo, int RANK){
    for (int i = 0; i<typeTwo.n_cNumbers; i++) typeTwo.sums[i] = 0;
    for (int i = 0; i<typeTwo.n_cNumbers; i++){
        for (int j =typeTwo.cIndices_Offsets[i]; j<typeTwo.cIndices_Offsets[i+1]; j++ ){
            typeTwo.sums[i]+=grid.basis[typeTwo.cIndices[j]];
        }
    }
}

// ********************************************************************************************** //
void sendAndRecieve(Grid& grid, TypeTwo& typeTwo, int RANK){
    for (int i = 0; i<typeTwo.n_edges; i++){
        MPI_Isend(typeTwo.sums, typeTwo.n_sendCells, MPI_DOUBLE,
        typeTwo.edges[i], 666,MPI_COMM_WORLD, &typeTwo.send_request[i]);
        MPI_Irecv(&typeTwo.recvBuff[typeTwo.recvBuff_Offsets[i]],
        typeTwo.recvBuff_Offsets[i+1]-typeTwo.recvBuff_Offsets[i], MPI_DOUBLE,
        typeTwo.edges[i], 666, MPI_COMM_WORLD, &typeTwo.recv_request[i]);
    }
    MPI_Waitall(typeTwo.n_edges, typeTwo.send_request, typeTwo.statSend);
    MPI_Waitall(typeTwo.n_edges, typeTwo.recv_request, typeTwo.statRecv);

    for (int i = 0; i<typeTwo.n_sendCells; i++){
        for (int j=typeTwo.recvIndices_Offsets[i]; j<typeTwo.recvIndices_Offsets[i+1]; j++ ){
            typeTwo.sums[i] += typeTwo.recvBuff[typeTwo.recvIndices[j]];
        }
    }
}

// ********************************************************************************************** //
void TTnormalize(Grid& grid, TypeTwo& typeTwo, int RANK){
    for (int i = 0; i<typeTwo.n_cNumbers; i++){
        for (int j =typeTwo.cIndices_Offsets[i]; j<typeTwo.cIndices_Offsets[i+1]; j++ ){
            grid.basis[typeTwo.cIndices[j]]/=typeTwo.sums[i];
        }
    }
}
















//yolo
