//
// Created by Jobs on 4/23/2015.
//

#include "MatrixReaderNew.h"

MatrixXd MatrixReaderNew::getInputData(string fileName) {
    vector <double> v;
    int rows=0;
    int cols=0;

    return import_matrix_from_txt_file(fileName,v,rows,cols);
}
MatrixReaderNew::MatrixReaderNew(int num) {
    num=num*2;
}
MatrixReaderNew::~MatrixReaderNew() {

}
int MatrixReaderNew::ReadNumbers(const string &s, vector<double> &v) {
    istringstream is( s );
    double n;
    while( is >> n ) {
        v.push_back( n );
    }
    return v.size();
}
MatrixXd MatrixReaderNew::import_matrix_from_txt_file(string filename_X, vector<double> &v, int &rows,
                                                      int &cols) {
    ifstream file_X;
    string line;
    file_X.open(filename_X.c_str());
    if (file_X.is_open())
    {
        int i=0;
        getline(file_X, line);


        cols =ReadNumbers( line, v );
        cout << "cols:" << cols << endl;


        for ( i=1;i<32767;i++){
            if ( getline(file_X, line) == 0 ) break;
            ReadNumbers( line, v );

        }

        rows=i;
        cout << "rows :" << rows << endl;
        if(rows >32766) cout<< "N must be smaller than MAX_INT";

        file_X.close();
    }
    else{
        cout << "file open failed";
    }

    cout << "v:" << endl;
    MatrixXd inputData(rows,cols);
    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++) {
            inputData(i, j) = v[i * cols + j];
        }
    }
    return inputData;

}

/*

int main(){
    MatrixReaderNew matrixReaderNew(5);

    cout<<matrixReaderNew.getInputData("x.txt")<<endl;

}

 */