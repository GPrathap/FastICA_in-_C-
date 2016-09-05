#ifndef CONTRASTFUNCTIONS_H
#define CONTRASTFUNCTIONS_H


#include <Eigen/Dense>
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;
// The virtual wrapper for the contrast functions
// also includes a factory method to create a
// specific subclass.
template <class T>
class ContrastFunction
{
public:
    virtual void operator()( MatrixXd, MatrixXd&, MatrixXd& ) = 0;

    enum Type
    {
        HYPERBOLIC_TAN,
        CUBIC,
        GAUSSIAN
    };

    static ContrastFunction<T>* getContrastFunction( Type type );
};

template <class T>
class HyperbolicTan : public ContrastFunction<T>
{
    void operator()( MatrixXd, MatrixXd&, MatrixXd& );
};


template <class T>
class Cubic : public ContrastFunction<T>
{

    void operator()( MatrixXd, MatrixXd&, MatrixXd& );
};

template <class T>
class Gaussian : public ContrastFunction<T>
{
    void operator()( MatrixXd, MatrixXd&, MatrixXd& );
};


#endif

