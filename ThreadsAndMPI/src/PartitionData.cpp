//
// Created by Jobs on 4/24/2015.
//

#include "PartitionData.h"
#include "MPIHolder.h"
#include <iostream>
#include <stdio.h>
#include <mpi.h>
#include <vector>
using namespace std;
void PartitionData::broadcastMatrix(MatrixXd &data) {
    //MPI_Init(NULL, NULL);
    #ifdef __DISTRIBUTED
    int rank;
    int size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    int headerPacket[1];
    if ( 0 == rank )
    {
        if ( data.rows() != 0 )
        {
            headerPacket[0] = data.rows();
            headerPacket[1] = data.cols();
        }
        else
        {
            headerPacket[0] = 0;
            headerPacket[1] = 0;
        }
    }
    MPI_Bcast( headerPacket, 2, MPI_INTEGER, 0, MPI_COMM_WORLD );

    int rows = headerPacket[0];
    int columns = headerPacket[1];

    if ( rows != 0 )
    {
         if ( data.rows()==0 ) {
             data.resize( rows, columns );

        }
//        cout<<"=========================before b cast================================"<<endl;
//        cout<<data<<"Rank "<< rank<<endl;
        MPI_Bcast( data.data(), rows * columns, MPI_DOUBLE, 0, MPI_COMM_WORLD );
//        cout<<"=========================after b cast================================"<<endl;
//        cout<<data<<"Rank "<< rank<<endl;
    }
    #endif
   // MPI_Finalize();
}
void PartitionData::broadcastNumbersOfThreads(int& maxNumberOfThreads, vector<int> threadsDistribution) {

    int rank;
    int size;
    //int* maxThreads = new int[5];
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    int* maxThreads = new int[1];

    int* sendCounts = new int[ size ];
    int* sendBase = new int[ size ];
    sendBase[ 0 ] = 0;
    int numOfThreadsPerProc = threadsDistribution.size() / size;
    int extraThreads = threadsDistribution.size() % size;
    int minDataCount = numOfThreadsPerProc;
    for ( int i = 0; i < size; ++i )
    {
        sendCounts[ i ] = minDataCount;
        if ( i < extraThreads ) sendCounts[ i ] += 1;
        if ( 0 != i ) sendBase[ i ] = sendBase[ i - 1 ] + sendCounts[ i - 1 ];
    }

    int* dataSet= new int[threadsDistribution.size()];

    for(int k=0;k<threadsDistribution.size();k++) {
                dataSet[k]=threadsDistribution[k];

    }

    MPI_Scatterv(
            dataSet, sendCounts, sendBase, MPI_INT,
            maxThreads, sendCounts[ rank ], MPI_INT,
            0, MPI_COMM_WORLD );

    cout<<"Threds split method :"<<maxThreads[0] <<endl;
     maxNumberOfThreads=maxThreads[0];

}

void PartitionData::broadcastNumbersOfThreads(int& maxNumberOfThreads, int totalNumberOfThreads) {

    int rank;
    int size;
    //int* maxThreads = new int[5];
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    int* maxThreads = new int[1];

    int* sendCounts = new int[ size ];
    int* sendBase = new int[ size ];
    sendBase[ 0 ] = 0;
    int numOfThreadsPerProc = totalNumberOfThreads / size;
    int extraThreads = totalNumberOfThreads % size;
    int minDataCount = numOfThreadsPerProc;
    for ( int i = 0; i < size; ++i )
    {
        sendCounts[ i ] = minDataCount;
        if ( i < extraThreads ) sendCounts[ i ] += 1;
        if ( 0 != i ) sendBase[ i ] = i;
    }

    int* dataSet= new int[size];

    for(int k=0;k<size;k++) {
        dataSet[k]=1;

    }

   /* if(rank==0){
        for(int k=0;k<size;k++) {
            //cout<<"================================"<<endl;
            //cout<<sendCounts[k]<<endl;

        }
    }*/

    MPI_Scatterv(
            sendCounts, dataSet, sendBase, MPI_INT,
            maxThreads, dataSet[ rank ], MPI_INT,
            0, MPI_COMM_WORLD );

    //cout<<"  :"<<maxThreads[0] <<endl;
    maxNumberOfThreads=maxThreads[0];

}


void PartitionData::distributeMatrix(MatrixXd& mat) {


    int rank;
    int size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    int headerPacket[2];

    if ( 0 == rank )
    {
        headerPacket[0] = mat.rows();
        headerPacket[1] = mat.cols();
    }
    MPI_Bcast( headerPacket, 2, MPI_INTEGER, 0, MPI_COMM_WORLD );

    int numRows = headerPacket[0];
    int numColumns = headerPacket[1] / size;
    if ( rank < (headerPacket[1] % size) ) ++numColumns;

    MPIHolder recvMatrix(numRows, numColumns, headerPacket[0], headerPacket[1], rank, size , NULL);

    // prepare for the scatter communication
    int* sendCounts = new int[ size ];
    int* sendBase = new int[ size ];
    sendBase[ 0 ] = 0;
    int columnsPerProc = headerPacket[ 1 ] / size;
    int extraColumns = headerPacket[ 1 ] % size;
    int minDataCount = columnsPerProc * headerPacket[ 0 ];
    for ( int i = 0; i < size; ++i )
    {
        sendCounts[ i ] = minDataCount;
        if ( i < extraColumns ) sendCounts[ i ] += headerPacket[ 0 ];
        if ( 0 != i ) sendBase[ i ] = sendBase[ i - 1 ] + sendCounts[ i - 1 ];
    }
    double* metaData= mat.data();
    // send the data
    MPI_Scatterv(
            metaData, sendCounts, sendBase, MPI_DOUBLE,
            recvMatrix.data, sendCounts[ rank ], MPI_DOUBLE,
            0, MPI_COMM_WORLD );

//    cout<< "." << sendCounts[ rank ]<<endl;
    MatrixXd received =  Map<MatrixXd>( recvMatrix.data, numRows, numColumns);
    mat.resize(received.rows(),received.cols());
    mat=received;

    delete[] sendCounts;
    delete[] sendBase;
   // MPI_Finalize();

}

void PartitionData::gatherMatrix(MatrixXd &data,int totalRows, int totalCols) {

    #ifdef __DISTRIBUTED
    int rank;
    int size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    double* d = data.data();

    int columns = (int)data.cols();
    int rows = (int)data.rows();
    int dataLength = columns * rows;

    MatrixXd recvData;
    int totalDataLength = totalRows*totalCols;

    if ( rank == 0 )
    {
       // recvData = new T[ totalDataLength ];
        recvData.resize(totalRows, totalCols);

    }

    // prepare to scatter
    int* recvCounts = new int[ size ];
    int* recvBase = new int[ size ];
    recvBase[ 0 ] = 0;
    int columnsPerProc = totalCols / size;
    int extra = totalCols % size;
    int minDataCount = columnsPerProc * rows;
    for ( int i = 0; i < size; ++i )
    {
        recvCounts[ i ] = minDataCount;
        if ( i < extra ) recvCounts[ i ] += rows;
        if ( 0 != i ) recvBase[ i ] = recvBase[ i - 1 ] + recvCounts[ i - 1 ];
    }

    MPI_Gatherv( d, dataLength, MPI_DOUBLE,  // send information
                 recvData.data(), recvCounts, recvBase, MPI_DOUBLE, // recv info
                 0, MPI_COMM_WORLD );

    if ( rank == 0 )
    {
        data.resize(totalRows, totalCols);
        data= recvData;
        delete[] d;
    }
    #endif

}
