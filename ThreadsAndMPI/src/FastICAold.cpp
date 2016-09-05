//
// Created by Jobs on 4/19/2015.
//

#include "FastICAold.h"



#include "MPIWrapper.h"

#include <math.h>
#include <cstring>
#include <stdio.h>
#include<iostream>
#include "Util.h"
using Eigen::ComputeFullV;
using Eigen::ComputeThinV;
using Eigen::ComputeThinU;
using Eigen::ArrayXXd;
using Eigen::EigenSolver;
using Eigen::MatrixXcd;
using namespace std;
template <class T> FastICAold<T>::FastICAold(
        const FastICASettings<T>& settings,
        MatrixXd& data,
        unsigned long seed,
        bool verbose, int tRows, int tCols, MatrixXd& w_init, int start, int end, int threadRank )


{
    m_settings=settings;
    m_contrastFunction=ContrastFunction<T>::getContrastFunction(m_settings.contrastFunction ) ;
    m_data.resize(data.rows(), data.cols());
    m_data= data;
    m_mixingMatrix.resize( data.rows(), data.rows() );
    m_channels=m_data.rows();
    m_samples=m_data.cols();
    totalCols=tCols;
    totalRows=tRows;
    m_commSize=getMPISize();
    m_verbose=verbose;
    m_whiteningMatrix.resize(data.rows(), data.rows());
    m_whiteningMatrix=w_init;
    m_start=start;
    m_end=end;
    m_threadRank=threadRank;
}

template <class T>
FastICAold<T>::~FastICAold()
{
    // if ( m_whiteningMatrix != 0 ) delete m_whiteningMatrix;
    // if ( m_invWhiteningMatrix != 0 ) delete m_invWhiteningMatrix;
    // if ( m_initialGuess != 0 ) delete m_initialGuess;
}

template <class T>
void
FastICAold<T>::mixingMatrix( MatrixXd& matrix, bool spheringApplied )
{
    if ( spheringApplied)
    {
        matrix=m_whiteningMatrix*m_mixingMatrix;
    }
    else
    {
        matrix = m_mixingMatrix;
    }
}

template <class T>
void
FastICAold<T>::setWhiteningMatrix( MatrixXd& matrix )
{
    m_whiteningMatrix.resize(matrix.rows(),matrix.cols());
    m_whiteningMatrix = matrix;
}


template <class T>
void FastICAold<T>::runFastICA()
{

       /* MatrixXd X1=m_data;
        MatrixXd Tw(m_channels,m_channels);
        MatrixXd Tw_d(m_channels,m_channels);
        MatrixXd w_init=m_whiteningMatrix;
        ArrayXXd gwtx_into_x1transpose_p;
        ArrayXXd gwtx_into_W;
        MatrixXd input_to_symmetric;
        MatrixXd W1;
        double lim;	//limit to check with tolerance
        bool success = false;
        MatrixXd W = Util::sym_decorrelation(w_init);
        for(int i=0;i<m_settings.maximumIterations;i++){

            MatrixXd dotprod;
            MatrixXd gwtx(m_channels,m_samples);
            MatrixXd g_wtx(m_channels,1);
            dotprod = W*X1;
            (*m_contrastFunction)( dotprod, gwtx, g_wtx);

            gwtx_into_x1transpose_p = (gwtx * X1.transpose()).array()/m_samples;
            allreduceDataMPIMNew(gwtx_into_x1transpose_p,Tw,m_channels);
            gwtx_into_W = Util::multiplyColumnWise(g_wtx,W);
            allreduceDataMPIMNew(gwtx_into_W,Tw_d,m_channels);
            input_to_symmetric=Tw-Tw_d;
            W1 = Util::sym_decorrelation(input_to_symmetric);
            lim = Util::estimateConvergence(W1,W);
            W = W1;
            if(lim<m_settings.convergenceTolerance){
                success = true;
                break;
            }

        }
        if(!success){
            cout<<"!!!!! did not converged, increase the max_iter count!!!!!"<<endl;
        }
        return W;*/

    MatrixXd X1=m_data;

    int maxIter=0;
    MatrixXd Tw(m_channels,m_channels);
    MatrixXd Tw_d(m_channels,m_channels);

    MatrixXd w_init=m_whiteningMatrix;
    //std::cerr<<w_init<<std::endl;
    ArrayXXd gwtx_into_x1transpose_p;
    ArrayXXd gwtx_into_W;
    MatrixXd input_to_symmetric;
    MatrixXd W1;
    double lim;
    bool success = false;
    MatrixXd W = Util::sym_decorrelation(w_init);

    //for(int i=0;i<m_settings.maximumIterations;i++){
    while(maxIter<m_settings.maximumIterations)
    {
        maxIter++;
        //cout<<"inside thread "<<m_threadRank<<endl;
        std::unique_lock<std::mutex> lkf(lockWeighMat);
        data_sends.wait(lkf,[this]{return checkeWeight[m_threadRank]!=1;});
        if(isTerminate==true){
            cerr<<"Thread:\t"<<this->m_threadRank<<" has been terminated"<<endl;
            break;
        }
        checkeWeight[m_threadRank]=1;
        W=weight_mat;
        lkf.unlock();

        MatrixXd dotprod;
        MatrixXd gwtx(m_channels,m_samples);
        MatrixXd g_wtx(m_channels,1);
        //cout<<"=====================================0================================"<<W.rows()<<W.cols()<<X1.rows()<<X1.cols()<<endl;
        dotprod = W*X1;
        //cout<<"=====================================1================================"<<endl;
        //cout<<"=====================================2================================"<<dotprod.rows()<<dotprod.cols()<<endl;
        (*m_contrastFunction)( dotprod, gwtx, g_wtx);
        //cout<<"=====================================3================================"<<endl;
        gwtx_into_x1transpose_p = (gwtx * X1.transpose()).array()/m_samples;
        //cout<<"=====================================4================================"<<endl;
        gwtx_into_W = Util::multiplyColumnWise(g_wtx,W);
        //cout<<"=====================================5================================"<<endl;
        std::unique_lock<std::mutex> lk(lockWorker);
        data_prepare.wait(lk,[this]{return checkedEvent[m_threadRank]!=1;});
        //cout<<"=====================================6================================"<<endl;
        MatrixXd temp1 =  gwtx_into_x1transpose_p.block(0,0,this->m_channels,this->m_channels);
        //cout<<"=====================================7================================"<<endl;
        gFunction+= temp1;
        //cout<<"=====================================8================================"<<endl;
        MatrixXd temp2 =gwtx_into_W.block(0,0,this->m_channels,this->m_channels);
        derivativeOfGFunction+=temp2;
        checkedEvent[m_threadRank]=1;
        IsPrepare=IsPrepare+1;
        std::cerr <<"Process is done: "<<IsPrepare<<" <<"<<std::endl;
        if(m_threadRank==0){
            numberOfIteratios++;
        }
        lk.unlock();
    }

}



template <class T>
void FastICAold<T>::initialize(VectorXd& w, int channel, int retries )
{
    //printf("step n-1 ------------------------------------------------------------------------------------>\n");
    if ( m_commRank != 0 )
    {
        if ( m_settings.initializationType == FastICASettings<T>::IDENTITY )
        {
            //printf("step n+2------------------------------------------------------------------------------------>\n");
            for ( int i = 0; i < m_channels; ++i )
            {
                if ( i == channel ){
                    w[i] = 1;
                }
                else{
                    w[i] = 0;
                }
            }
        }

    }
    //broadcastDataMPI( w, m_channels, 0 );
}

#ifdef INSTANTIATE_TEMPLATES
template class FastICAold<double>;
template class FastICAold<float>;
#endif



