#ifndef MPIWRAPPER_CPP_UNIVERSITY_OF_OREGON_NIC
#define MPIWRAPPER_CPP_UNIVERSITY_OF_OREGON_NIC

#include "MPIWrapper.h"
#include <stdio.h>
#include "TypeIdentifier.h"
using namespace std;
#ifdef __DISTRIBUTED

int getMPIRank()
{
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    return rank;
}

int getMPISize()
{
    int size = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    return size;
}

template <class T>
void broadcastDataMPI( T* data, int length, int root )
{
    MPI_Bcast( data, length, TypeIdentifier<T>::Type(), root, MPI_COMM_WORLD );
}


void allreduceDataMPIM( VectorXd input, VectorXd& output, int length )
{
    // printf("inside allreduce finction ..............................");
     MPI_Allreduce( input.data(), output.data(), length,
        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
}
void allreduceDataMPIMNew( MatrixXd input, MatrixXd& output, int length )
{
     printf("allreduceDataMPIMNew ..............................");
     MPI_Allreduce( input.data(), output.data(), length*length,
        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
}
void initializeMPI( int argc, char** argv, int& rank, int& size )
{
    //MPI_Init( &argc, &argv );
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
}

void finalizeMPI()
{
    MPI_Finalize();
}

#else

int getMPIRank()
{
    return 0;
}

int getMPISize()
{
    return 1;
}

template <class T>
void broadcastDataMPI( T* data, int length, int root )
{
    return;
}


void allreduceDataMPIM( VectorXd input, VectorXd& output, int length )
{
    for ( int i = 0; i < length; ++i )
    {
        output[i] = input[i];
    }
}
void allreduceDataMPIMNew( MatrixXd input, MatrixXd& output, int length )
{
     printf("inside allreduce finction ..............................");
    output=input;
}
void initializeMPI( int argc, char** argv, int& rank, int& size )
{
    rank = 0;
    size = 1;
}

void finalizeMPI()
{
}

#endif

#ifdef INSTANTIATE_TEMPLATES
template void broadcastDataMPI<double>( double* data, int length, int root );


template void broadcastDataMPI<float>( float* data, int length, int root );

#endif

#endif
// MPIWRAPPER_CPP_UNIVERSITY_OF_OREGON_NIC
