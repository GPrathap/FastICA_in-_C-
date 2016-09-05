//
// Created by Jobs on 4/24/2015.
//

#ifndef HIPERSAT1_1_MPIHOLDER_H
#define HIPERSAT1_1_MPIHOLDER_H
#include <Eigen/Dense>
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;

class MPIHolder {
     public:
        int rows;
        int columns;
        int size;
        int total_size;
        int total_rows;
        int total_columns;
        int n_cpus;
        int rank;
        int referenceCount;
        bool isAlias;
        double* data;
        MatrixXd received;

         MPIHolder(int r = 0, int c = 0, int t_r = 0, int t_c = 0,
             int rk = 0, int n_c = 1, double* m= NULL);
         ~MPIHolder();

         double** allocateMemory(int rows, int cols);
};


#endif //HIPERSAT1_1_MPIHOLDER_H
