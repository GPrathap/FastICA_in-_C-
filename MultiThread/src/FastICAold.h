//
// Created by Jobs on 4/19/2015.
//

#ifndef NEW_FASTICAOLD_H
#define NEW_FASTICAOLD_H


#include "../src/FastICASettings.h"
#include "ContrastFunctions.h"
#include <Eigen/Dense>
#include "sharedData.h"
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;
template <class T>
class FastICAold
{


public:
    FastICAold(
            const FastICASettings<T>& settings,
            MatrixXd& data,
            unsigned long seed = 123456UL,
            bool verbose = false , int totalRows =0, int totalCols=0, MatrixXd& w_init=NULL, int start=0,int end=0, int threadRank=0);

    ~FastICAold();

    // get the computed mixing matrix
    void mixingMatrix( MatrixXd& matrix, bool spheringApplied = true );

    // set the whitening (aka sphering) matrix
    void setWhiteningMatrix( MatrixXd& matrix );

    void orthogonalize(VectorXd& w, MatrixXd A);
    void signChange(MatrixXd& m);
    void getMaxSign(VectorXd& v);
    //void allreduceDataMPI( VectorXd input, VectorXd& output, int length );
    void normalize( VectorXd& w);
    void Gaussian( MatrixXd x_, VectorXd& weight, VectorXd& weight_plus);
    void HyperbolicTan( MatrixXd x_, VectorXd& weight, VectorXd& weight_plus);
    void Tger(VectorXd w, VectorXd w1, MatrixXd& m);
    MatrixXd getRandomMatrix(int n);
    // compute the mixing matrix
    void runFastICA();
    MatrixXd  sym_decorrelation(MatrixXd& m);
    T estimateConvergence(MatrixXd w0, MatrixXd w1);
    int m_channels; // rows of data
    int m_samples;  // columns of data
    int totalCols;
    int totalRows;
    MatrixXd m_whiteningMatrix;
    MatrixXd m_invWhiteningMatrix;
    ContrastFunction<T>* m_contrastFunction;
    // MatrixXd m_data1(4,4);
private:
    // private methods
    void initialize( VectorXd& w, int channel, int retries );
    void loadGuess();



    //private members
    // settings
    FastICASettings<T> m_settings;

    // data
    MatrixXd m_data;
    MatrixXd n;
    MatrixXd m_mixingMatrix;

    MatrixXd m_initialGuess;

    //ContrastFunction<T>* m_contrastFunction;

    int m_threadRank;
    int m_commRank;
    int m_commSize;
    int m_start;
    int m_end;
    bool m_verbose;

};



#endif //NEW_FASTICAOLD_H
