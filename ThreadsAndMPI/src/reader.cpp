
#include "reader.h"

reader::reader(int f) {
}
void reader::readInputData(int row,int column,MatrixXd& X,string filename_X ) {
  
	
	FILE * fp = fopen(filename_X.c_str(),"r");
	int  i,j;
	double temp;
  
    for(i=0;i<row;i++){
	  for(j=0;j<column;j++){
		  int d=fscanf(fp,"%lf",&temp);
		  d=0;
		  X(i,j) = temp;
	  }
	  while(fgetc(fp)!='\n');	  
   }
  
}

	

