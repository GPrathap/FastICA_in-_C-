#ifndef CONTRASTFUNCTIONS_CPP
#define CONTRASTFUNCTIONS_CPP

#include "ContrastFunctions.h"
#include <math.h>
#include <stdio.h>
#include<iostream>
#include <Eigen/Dense>

using namespace std;
using Eigen::ArrayXXd;
using Eigen::MatrixXcd;

template <class T>
ContrastFunction<T>*
ContrastFunction<T>::getContrastFunction( Type type )
{
    ContrastFunction* returnValue;
    switch ( type )
    {
        case HYPERBOLIC_TAN:
            returnValue = new HyperbolicTan<T>;
            break;
        case CUBIC:
            returnValue = new Cubic<T>;
            break;
        case GAUSSIAN:
            returnValue = new Gaussian<T>;
            break;
        default:
            returnValue = 0;
    }
    return returnValue;
}

// Purpose: compute the covariance matrix and its eigenvalues
template <class T>
void HyperbolicTan<T>::operator()
    ( MatrixXd x, MatrixXd& weight, MatrixXd& weight_plus)
{

   /* int n_channels = x->rows;
    int n_samples = x->columns;
    int total_samples = x->total_columns;

    T* whitened_x = x->data;
    T* w = weight->data;
    T* w_plus = weight_plus->data;

    int i = 0, j = 0;
    T* x_i = NULL;
    T xTw = 0;
    T ht = 0;

    // make sure w+ is cleared

    for (i = 0; i < n_samples; i++)
    {
        xTw = 0;
        x_i = whitened_x + i * n_channels;

        for (j = 0; j < n_channels; j++)
        {
            xTw += x_i[j] * w[j];
        }

        //    ht = tanh(xTw);
        ht = 1 - ( 2 / ( exp(2*xTw) + 1 ) );

        for (j = 0; j < n_channels; j++)
        {
            w_plus[j] += x_i[j] * ht - (1-ht*ht) * w[j];
        }
    }


    for (j = 0; j < n_channels; j++)
    {   // E(x g(wTx))
        w_plus[j] /= total_samples;
    }*/

}

// Purpose: compute the covariance matrix and its eigenvalues
template <class T>
void Cubic<T>::operator()
    ( MatrixXd x, MatrixXd& weight, MatrixXd& weight_plus)
{


    int n =x.rows();
    ArrayXXd y = x.array();
    y*=(y*y);
    weight= y.matrix();

    ArrayXXd y_d = x.array();
    y_d*=y_d;
    MatrixXd xx =(3*y_d).matrix();
    MatrixXd  means(n,1);

    int i;
    for(i=0;i<n;i++){
        means(i,0) = xx.row(i).mean();
    }

    weight_plus=means;
}

// Purpose: compute the covariance matrix and its eigenvalues
template <class T>
void Gaussian<T>::operator()
    ( MatrixXd x_, MatrixXd& weight, MatrixXd& weight_plus)
{
   /* T* whitened_x = x_->data;
    T* w = weight->data;
    T* w_plus = weight_plus->data;
    int n_samples = x_->columns;
    int n_channels = x_->rows;

    T* x = NULL;
    T wTx = 0;
    T wTx_2 = 0;
    T exp_wTx_2 = 0;

    int i = 0, j = 0;
    printf("============inside contention free ==================\n");
    for (i = 0; i < n_samples; i++)
    {
        wTx = 0;
        x = whitened_x + i * n_channels;

        // compute wTx
        for (j = 0; j < n_channels; j++)
        {
            wTx += w[j] * x[j];
            //printf("===%f=====%f===\n",w[j],x[j]);
        }
      // printf("\n");
       // printf("\n");
        printf(" wTx %f ",wTx);
        wTx_2 = wTx * wTx;
        exp_wTx_2 = exp(-wTx_2/2);
        printf("wTx_2 %f  ",wTx_2);
        printf("exp_wTx_2 %f\n",exp_wTx_2);
        for (j = 0; j < n_channels; j++)
        {   // E(x g(wTx))
            w_plus[j] +=  ( wTx * exp_wTx_2 * x[j] ) -
                          ((1 - wTx_2) * exp_wTx_2 * w[j]);
            printf(" %f ",w_plus[j]);
        }
        printf("=========================================================\n");
    }*/
}

#ifdef INSTANTIATE_TEMPLATES
template class ContrastFunction<double>;
template class ContrastFunction<float>;

template class HyperbolicTan<double>;
template class HyperbolicTan<float>;

template class Cubic<double>;
template class Cubic<float>;

template class Gaussian<double>;
template class Gaussian<float>;
#endif

#endif
// CONTRASTFUNCTIONS_CPP
