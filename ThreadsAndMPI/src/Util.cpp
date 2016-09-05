//
// Created by Jobs on 6/24/2015.
//

#include "Util.h"


Util::Util (int args)
{
    cout<<"init file"<<endl;
}

ArrayXXd Util::multiplyColumnWise(MatrixXd g_wtx,MatrixXd W){
    ArrayXXd W_in = W;
    ArrayXXd g_wtx_in =  g_wtx;
    int n = W_in.cols();
    int i;
    for(i=0;i<n;i++){
        W_in.col(i)*=g_wtx_in;
    }
    return W_in;
}
MatrixXd Util::devide(MatrixXd u,MatrixXd d,int cols){
    //each column of u should devide by each row of d
    int i;
    for(i=0;i<cols;i++){
        u.col(i) /= d(i,0);
    }
    return u;
}

void Util::addMatrix(MatrixXd mat,int rows, MatrixXd& resultMat)
{
    MatrixXd finalR(rows,rows);
    for(int i=0;i<mat.cols();i+=rows)
    {
        if(i==0){
            finalR=mat.block(0,i, rows,rows);
        }
        finalR+=mat.block(0,i, rows,rows);
    }
    resultMat=finalR;
}


void Util::splitData(int rows, int cols, int numberOfThreads, std::map<int,int>& result)
{
    int headerPacket[2];
    headerPacket[0]=rows;
    headerPacket[1]=cols;
    int* sendCounts = new int[ numberOfThreads ];
    int sendBase[numberOfThreads];
    sendBase[ 0 ] = 0;
    int columnsPerProc = headerPacket[ 1 ] / numberOfThreads;
    int extraColumns = headerPacket[ 1 ] % numberOfThreads;
    int minDataCount = columnsPerProc * headerPacket[ 0 ];
    for ( int i = 0; i < numberOfThreads; ++i )
    {
        sendCounts[ i ] = minDataCount;
        if ( i < extraColumns ) sendCounts[ i ] += headerPacket[ 0 ];
        if ( 0 != i ) sendBase[ i ] = sendBase[ i - 1 ] + sendCounts[ i - 1 ];
    }

    for(int i=0; i<numberOfThreads;i++)
    {
        if(i+1<numberOfThreads)
        {
            result[sendBase[i]/headerPacket[0]]=(sendBase[i+1]/headerPacket[0])-1;
        }
        else
        {
            result[sendBase[i]/headerPacket[0]]=cols-1;
        }
    }

    return ;
}


void Util::maxSign(VectorXd& v)
{
    double max=abs(v[0]);
    double realValue=v[0];
    int k;
    for(k=0;k<v.size();k++){
        if(max<=abs(v[k])){max=abs(v[k]);realValue=v[k];}
    }
    if(realValue>0){

        for(k=0;k<v.size();k++){
            v[k]=-1*v[k];
        }
    }
    return;
}

void Util::changeSign(MatrixXd& mat){
    for(int i=0;i<mat.cols();i++){
        VectorXd v1 = mat.col(i);

        maxSign(v1);
        mat.col(i).array()=v1;
    }
    return ;
}

MatrixXd Util::getRandomMatrix(int n) {
    return MatrixXd::Random(n,n);
}

MatrixXd Util::getMean1(MatrixXd X,int n){
    MatrixXd  means(n,1);
    int i;
    for(i=0;i<n;i++){
        means(i,0) = X.row(i).mean();
    }
    return means;

}

MatrixXd Util::normalize(MatrixXd X,MatrixXd means,int rows)
{
    int i;
    for(i=0;i<rows;i++){
        X.row(i) = X.row(i).array() - means(i,0);
    }
    return X;

}

MatrixXd Util::generateRandomMatrix(int n){
    //MatrixXd m(n,n);

    return MatrixXd::Random(n,n);
}

void Util::getMaxSign(VectorXd& v)
{

    double max=abs(v[0]);
    double realValue=v[0];
    int count=1;
    int k;
    for(k=0;k<v.size();k++){
        if(max<=abs(v[k])){
            max=abs(v[k]);
            realValue=v[k];
            if(v[k]>0)count++;
        }
    }
    if(realValue<0){

        for(k=0;k<v.size();k++){
            v[k]=-1*v[k];
        }
    }
    //cout<<count << v.size()<<endl;
    if(count==v.size()){
        for(k=0;k<v.size();k++){
            v[k]=-1*v[k];
        }
    }
    return;

}

void Util::signChange(MatrixXd& mat) {
    for(int i=0;i<mat.cols();i++){
        VectorXd v1 = mat.col(i);

        getMaxSign(v1);
        mat.col(i).array()=v1;
    }
    return ;
}


MatrixXd Util::sym_decorrelation(MatrixXd& m) {

    SelfAdjointEigenSolver<MatrixXd> es(m*m.transpose());
    VectorXd eigenvalues = es.eigenvalues();
    MatrixXd eieVec = es.eigenvectors();
    signChange(eieVec);
    MatrixXd u =eieVec;
    MatrixXd v =eieVec;
    for(int i=0;i<eigenvalues.size();i++){
        eigenvalues[i]=1/sqrt(eigenvalues[i]);
    }
    for(int i=0;i<u.rows();i++){
        for(int j=0;j<u.cols();j++){
            u(i,j)=u(i,j)*eigenvalues[j];
        }
    }
    MatrixXd m1 = v*u.transpose();
    MatrixXd m2 = m1*m;
    m.resize(eigenvalues.size(),eigenvalues.size());
    m=m2;
    return m;
}

double Util::estimateConvergence(MatrixXd W1, MatrixXd W)
{
    return ((((((W1*W.transpose()).diagonal()).array()).abs()) - 1).abs()).maxCoeff();
}
