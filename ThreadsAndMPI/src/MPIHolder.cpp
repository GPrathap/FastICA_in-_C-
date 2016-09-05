//
// Created by Jobs on 4/24/2015.
//

#include "MPIHolder.h"
#include <iostream>
#include <stdio.h>

using namespace std;

MPIHolder::MPIHolder(int r, int c, int t_r, int t_c, int rk, int n_c, double* m) {
    rows=r;
    total_rows=t_r;
    columns=c;
    total_columns=t_c;
    rank=rk;
    size=r*c;
    total_size=t_c*t_r;
    n_cpus=n_c;
    if(m==NULL){
        data=new double[size];
        //cout<< "Rank ....."<< rank << endl;
        //cout<< "Rows ....."<< rows << endl;
       // cout<< "Columns ....."<< columns << endl;
    }else{
        data=m;
    }
    received.resize(0,0);
}

MPIHolder::~MPIHolder(){

};

double** MPIHolder::allocateMemory(int rows, int cols) {
    double** temp = 0;
    int d=0,j=0;
    temp = new double*[rows];

    for(d=0;d<rows;d++){
        temp[d]= new double[cols];
        for(j=0;j<cols;j++){
            temp[d][j]=0.0;
        }
    }
    return temp;
}