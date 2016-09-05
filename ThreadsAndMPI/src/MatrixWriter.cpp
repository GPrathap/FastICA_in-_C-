#ifndef MATRIXWRITER_CPP_UNIVERSITY_OF_OREGON_NIC
#define MATRIXWRITER_CPP_UNIVERSITY_OF_OREGON_NIC

#include <fstream>
#include "DataWriter.h"
#include "MatrixWriter.h"


using namespace std;

template <class T>
bool MatrixWriter<T>::writeMatrix(
    MatrixXd& matrix,
    DataFormat::DataFormat format,
    const std::string& fileName,
    bool transpose )
{
    // set the mode for file io
    ios::openmode mode = ios::out;
    if ( format != DataFormat::TEXT ) mode |= ios::binary;

    typename DataWriter<T>::Endian endian;
    typename DataWriter<T>::Format writerFormat;

    switch ( format )
    {
        case DataFormat::BIG:
            endian = DataWriter<T>::BigEndian;
            writerFormat = DataWriter<T>::Binary;
            break;
        case DataFormat::LITTLE:
            endian = DataWriter<T>::LittleEndian;
            writerFormat = DataWriter<T>::Binary;
            break;
        case DataFormat::NATIVE:
            endian = DataWriter<T>::NativeEndian;
            writerFormat = DataWriter<T>::Binary;
            break;
        case DataFormat::TEXT:
            endian = DataWriter<T>::NativeEndian;
            writerFormat = DataWriter<T>::Ascii;
            break;
        default:
            endian = DataWriter<T>::NativeEndian;
            writerFormat = DataWriter<T>::Ascii;
    }

    ofstream outputFile( fileName.c_str(), mode );
    if ( !outputFile.good() )
    {
        return false;
    }

    DataWriter<T> writer( outputFile, writerFormat, endian );

    int count = 0;

    if ( !transpose )
    {
         writer.write(matrix);
    }
    else
    {
         writer.write(matrix);
    }
        
    outputFile.close();
    return ( count == 0 );
}

#ifdef INSTANTIATE_TEMPLATES
template class MatrixWriter<double>;
template class MatrixWriter<float>;
#endif

#endif
// MATRIXWRITER_CPP_UNIVERSITY_OF_OREGON_NIC
