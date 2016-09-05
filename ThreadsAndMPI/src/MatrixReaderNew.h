//
// Created by Jobs on 4/23/2015.
//

#ifndef NEW_MATRIXREADERNEW_H
#define NEW_MATRIXREADERNEW_H


#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <Eigen/Dense>
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;

using namespace std;


class MatrixReaderNew {
public:
    MatrixXd getInputData(string fileName);
    int ReadNumbers( const string & s, vector <double> & v );
    MatrixXd import_matrix_from_txt_file(string filename_X, vector <double>& v, int& rows, int& cols);
    MatrixReaderNew(int num);
    ~MatrixReaderNew();
};


#endif //NEW_MATRIXREADERNEW_H
