//
// Created by Jobs on 6/24/2015.
//

#ifndef NEW_UTIL_H
#define NEW_UTIL_H

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <iostream>
#include <map>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ComputeFullV;
using Eigen::ComputeThinV;
using Eigen::ComputeThinU;
using Eigen::JacobiSVD;
using Eigen::SelfAdjointEigenSolver;
using Eigen::ArrayXXd;

class Util
{
public:
    Util(int args);
    static void addMatrix(MatrixXd mat,int rows, MatrixXd& resultMat);
    //static void splitData(int rows, int cols, int numberOfThreads, std::map<int,int>& result);
    static void maxSign(VectorXd& v);
    static MatrixXd getRandomMatrix(int n);
    static MatrixXd getMean1(MatrixXd X,int n);
    static MatrixXd normalize(MatrixXd X,MatrixXd means,int rows);
    static void getMaxSign(VectorXd& v);
    static void signChange(MatrixXd& mat);
    static MatrixXd sym_decorrelation(MatrixXd& m);
    static void data_processing();
    static MatrixXd devide(MatrixXd u,MatrixXd d,int cols);
    static void splitData(int rows, int cols, int numberOfThreads, std::map<int,int>& result);
    static void changeSign(MatrixXd& mat);
    static double estimateConvergence(MatrixXd W1, MatrixXd W);
    static MatrixXd generateRandomMatrix(int n);
    static ArrayXXd multiplyColumnWise(MatrixXd g_wtx,MatrixXd W);
};


#endif //NEW_UTIL_H
