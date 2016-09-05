//
// Created by Jobs on 6/14/2015.
//

#ifndef NEW_SHAREDDATA_H
#define NEW_SHAREDDATA_H

#include <Eigen/Dense>
#include <condition_variable>
#include <mutex>
#include <vector>
#include <iostream>
#include <thread>

extern  int numberOfThreadsFinal;
extern  int IsPrepare;
extern bool isTerminate;
extern std::condition_variable data_prepare;
extern std::condition_variable data_sends;
extern std::vector<int> checkedEvent;
extern std::vector<int> checkeWeight;
extern std::mutex lockWorker;
extern std::mutex lockWeighMat;
extern Eigen::MatrixXd gFunction;
extern Eigen::MatrixXd derivativeOfGFunction;
extern Eigen::MatrixXd weight_mat;
extern int numberOfIteratios;
extern int weightIsReady;
#endif
