
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


#include <fcntl.h>   
  
#include <stdio.h>
#include <stdlib.h>
#include <iostream>


#include <fstream>
#include <sys/time.h>
#include <time.h>
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;

using namespace std;
using namespace Eigen;


using Eigen::ArrayXXd ;


#define ArgCount 3
#define PRECISION 10
#define MAX_ITER 200
#define TOL 0.00000001 


class reader {
public:
    void  readInputData(int row,int column,MatrixXd& X,string fileName);
    reader(int f); 
};

