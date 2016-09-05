//
// Created by Jobs on 4/24/2015.
//

#ifndef HIPERSAT1_1_PARTITIONDATA_H
#define HIPERSAT1_1_PARTITIONDATA_H
#include <Eigen/Dense>
#include <vector>
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::Map;


class PartitionData {
public:
// forget this stuff. It's from the old signal cleaner days
    typedef char            FileType;
    enum
    {
        c_fileType_binaryBigEndian = 'B',
        c_fileType_binaryLittleEndian = 'L',
        c_fileType_binaryNativeEndian = 'P'
    };

    typedef char            FileOrientation;
    enum
    {
        c_fileOrientationColumnMajor = 'C',
        c_fileOrientationRowMajor = 'O'
    };


    void loadFile(
            int                 rank,
            const std::string&  dataSrc,
            MatrixXd data,
            int                 ch,
            int                 columns,
            int                 rem,
            int                 obs,
            int                 n_cpus,
            FileType            fileType,
            FileOrientation     fileOrientation );

    void storeFile( int rank, MatrixXd& mat, const std::string& m_whitenedFilename , int n_cpus=1, int* shuffle=NULL, int* invVector=NULL) ;
    static void distributeMatrix(MatrixXd& data );
    static void broadcastMatrix(MatrixXd& data );
    static void gatherMatrix(MatrixXd& ,int totalRows, int totalCols);
    static void broadcastNumbersOfThreads(int& maxNumberOfThreads, std::vector<int> threadsDistribution );
    static void broadcastNumbersOfThreads(int& maxNumberOfThreads, int totalNumberOfThreads );
};


#endif //HIPERSAT1_1_PARTITIONDATA_H
