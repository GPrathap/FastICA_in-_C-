//
// Created by Jobs on 4/19/2015.
//

#include "FastICAold.h"



//#include "MPIWrapper.h"

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
        bool verbose, int tRows, int tCols, MatrixXd& w_init, int start, int end, int threadRank)
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
    m_commSize=1;
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
void FastICAold<T>::runFastICA()
{

  


        MatrixXd X1=m_data;
        //std::cerr<<X1<<st d::endl;
        int maxIter=0;
        MatrixXd Tw(m_channels,m_channels);
        MatrixXd Tw_d(m_channels,m_channels);
        MatrixXd w_init=m_whiteningMatrix;
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

            std::unique_lock<std::mutex> lkf(lockWeighMat);
            data_sends.wait(lkf,[this]{return checkeWeight[m_threadRank]!=1;});
            if(isTerminate==true){
                //cout<<"Thread:\t"<<this->m_threadRank<<" has been terminated"<<endl;
                break;
            }
            checkeWeight[m_threadRank]=1;
            W=weight_mat;
            lkf.unlock();

            MatrixXd dotprod;
            MatrixXd gwtx(m_channels,m_samples);
            MatrixXd g_wtx(m_channels,1);
            dotprod = W*X1;
            (*m_contrastFunction)( dotprod, gwtx, g_wtx);
            gwtx_into_x1transpose_p = (gwtx * X1.transpose()).array()/m_samples;
            gwtx_into_W = Util::multiplyColumnWise(g_wtx,W);
            std::unique_lock<std::mutex> lk(lockWorker);
            data_prepare.wait(lk,[this]{return checkedEvent[m_threadRank]!=1;});
            MatrixXd temp1 =  gwtx_into_x1transpose_p.block(0,0,this->m_channels,this->m_channels);
            gFunction+= temp1;
            MatrixXd temp2 =gwtx_into_W.block(0,0,this->m_channels,this->m_channels);
            derivativeOfGFunction+=temp2;
            checkedEvent[m_threadRank]=1;
            IsPrepare=IsPrepare+1;
//            std::cerr <<"Process is done: "<<IsPrepare<<" <<"<<std::endl;
            if(m_threadRank==0){
                numberOfIteratios++;
            }
            lk.unlock();
        }
}


#ifdef INSTANTIATE_TEMPLATES
template class FastICAold<double>;
template class FastICAold<float>;
#endif



