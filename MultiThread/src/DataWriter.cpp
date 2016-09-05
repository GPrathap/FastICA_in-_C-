#ifndef DATAWRITER_CPP
#define DATAWRITER_CPP

#include "DataWriter.h"
#ifdef MINGW
#include "WinHeader.h"
#endif

template < class T >
DataWriter<T>::DataWriter( 
    std::ostream& os, 
    typename DataWriter<T>::Format format,
    typename DataWriter<T>::Endian endian )
:   m_sizeOfT( sizeof( T ) ),
    m_outputStream( os ),
    m_format( format ),
    m_endian( endian ),
    m_nextValueCopy( new char[ m_sizeOfT ] )
{
    if ( BYTE_ORDER == BIG_ENDIAN )
    {
        m_native = DataWriter<T>::BigEndian;
    }
    else
    {
        m_native = DataWriter<T>::LittleEndian;
    }
    if ( m_endian == DataWriter<T>::NativeEndian )
    {
        m_endian = m_native;
    }
}

template <class T>
DataWriter<T>::~DataWriter()
{
    delete[] m_nextValueCopy;
}

template <class T>
bool
DataWriter<T>::write( T& value )
{
    //return write( &value );
    return 0;
}

template <class T>
bool
DataWriter<T>::write( MatrixXd value )
{
    if ( m_format == DataWriter<T>::Ascii )
    {

        m_outputStream << value << " ";
    }
    else
    {
        m_outputStream << value << " ";
    }
    return true;
}

template <class T>
char*
DataWriter<T>::writeBuffer( char* buffer, T& value )
{
    char* v = reinterpret_cast< char* >( &value );
    for ( int i = 0; i < sizeof( T ); ++i )
    {
        buffer[i] = v[i];
    }
    return ( buffer + sizeof( T ) );
}

template < class T >
void
DataWriter<T>::reverseBytes( char* from, char* to )
{
    for ( int i = 0; i < m_sizeOfT; ++i )
    {
        to[ i ] = from[ m_sizeOfT - i - 1 ];
    }
}

template <class T>
bool
DataWriter<T>::linebreak()
{
    if ( m_format == DataWriter<T>::Ascii )
    {
        m_outputStream << std::endl;
    }
    return true;
}

#ifdef INSTANTIATE_TEMPLATES
template class DataWriter<double>;
template class DataWriter<float>;
template class DataWriter<int>;
template class DataWriter<short>;
#endif

#endif
// DATAWRITER_CPP_UNIVERSITY_OF_OREGON_NIC
