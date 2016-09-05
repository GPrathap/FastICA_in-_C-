

#include <stdio.h>
#include "hFastICA.h"


#include "../src/ProgramTools.h"
#include "../src/FastICAold.h"
#include "../src/reader.h"
#include "../src/PartitionData.h"
#include "../src/DataFormat.h"
#include "../src/FastICAProgram.h"
#include "../src/Util.h"
#include <iostream>
#include <Eigen/Dense>
#include <sys/time.h>
#include <fstream> 
#include <time.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ComputeFullV;
using Eigen::ComputeThinV;
using Eigen::ComputeThinU;
using Eigen::JacobiSVD;
//function to measure time
typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp (){
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
      }



int numberOfThreadsFinal;
int IsPrepare=0;
bool isTerminate=false;

std::condition_variable data_prepare;
std::condition_variable data_sends;
std::vector<int> checkedEvent;
std::vector<int> checkeWeight;
std::mutex lockWorker;
std::mutex lockWeighMat;
Eigen::MatrixXd gFunction;
Eigen::MatrixXd derivativeOfGFunction;
Eigen::MatrixXd weight_mat;
int numberOfIteratios;
int maxNumberOfIterations;
int numberOfChannels;
int numberOfSamples;
double convergenceTolerance;
std::vector<std::thread> thredArray;
int mpi_proces_rank;
int weightIsReady = true;

void data_processing()
{

    double lim;
    MatrixXd Tw(numberOfChannels, numberOfChannels);
    MatrixXd Tw_d(numberOfChannels, numberOfChannels);
    while(true)
    {
        std::cout<<"";
        while(IsPrepare==numberOfThreadsFinal)
        {

            std::lock_guard<std::mutex> lk(lockWorker);

            MatrixXd input_to_symmetric=gFunction;
            MatrixXd input_to_symmetricDerivative=derivativeOfGFunction;
            MPI_Barrier(MPI_COMM_WORLD);
            allreduceDataMPIMNew(input_to_symmetric,Tw,numberOfChannels);
            //MPI_Barrier(MPI_COMM_WORLD);
            //allreduceDataMPIMNew(input_to_symmetricDerivative,Tw_d,numberOfChannels);
            //input_to_symmetric=Tw-Tw_d;
            //input_to_symmetric=input_to_symmetric-input_to_symmetricDerivative;
            MatrixXd W1 = Util::sym_decorrelation(input_to_symmetric);
            lim = Util::estimateConvergence(input_to_symmetric,W1);
            std::cerr<<"Lim:\t"<<lim<<"\tNumber of Iter: "<<convergenceTolerance<<std::endl;
            if((lim<convergenceTolerance)){
                isTerminate=true;
                std::fill(checkeWeight.begin(), checkeWeight.end(), 0);
                data_sends.notify_all();
                data_prepare.notify_all();
                if(mpi_proces_rank==0){
                    cout<<"Solution is converged ..............................................."<<endl;
                }

                goto stop;
            }else if(numberOfIteratios>=maxNumberOfIterations){
                if(mpi_proces_rank==0){
                    cout<<"Solution is not converged ..............................................."<<endl;
                }
                goto stop;
            }
            std::fill(checkedEvent.begin(), checkedEvent.end(), 0);
            IsPrepare=0;
            gFunction<< MatrixXd::Zero(numberOfChannels,numberOfChannels);
            derivativeOfGFunction<< MatrixXd::Zero(numberOfChannels,numberOfChannels);
            std::lock_guard<std::mutex> lkf(lockWeighMat);
            weight_mat = W1;
            std::fill(checkeWeight.begin(), checkeWeight.end(), 0);
            data_sends.notify_all();
            data_prepare.notify_all();

        }
    }

    stop:

    cout<<"End of processing thread"<<endl;

}


template <class T> int mainWrapper( SimpleCLParser& parser, int rank, int size )
{
    int rows = parser.readIntOption( "-c" );
    int columns = parser.readIntOption( "-s" );
    numberOfChannels= parser.readIntOption( "-c" );
    numberOfSamples=  parser.readIntOption( "-s" );
    int numOfProcesses = parser.readIntOption( "-pc" );
    //numberOfThreadsFinal = parser.readIntOption( "-th" );
    maxNumberOfIterations= parser.readIntOption( "-I" );
    convergenceTolerance = parser.readDoubleOption( "-t" );
    mpi_proces_rank=rank;
    cout<<"Final number of threads ...."<<rank<<endl;
    int n=rows;
    int p=columns;


    DataFormat::DataFormat outputDataFormat = getOutputFormat( parser );

    MatrixXd inputData(rows,columns);
    MatrixXd spheringMatrix(rows,rows);
    int maxNumberOfThreads;
    MatrixXd avg ;
    int totalRows;
    int totalCols;
    MatrixXd centered;
    MatrixXd w_init(n,n);	//random matrix
    MatrixXd W(n,n);	//result
    MatrixXd unmixedSignal;	//unmixing X using W
    MatrixXd S;
	timestamp_t t0 = get_timestamp();
	timestamp_t t1;
	timestamp_t t2;
	timestamp_t t3;
	timestamp_t t4;
	timestamp_t t5;
    std::vector<int> threadsDistribution;
    int l=0;

    for(l=0;l<numOfProcesses;l++){
        threadsDistribution.push_back(parser.readIntOption( "-th" ));
    }
	std::ofstream out("testResult.txt", std::ios::app);

    if ( rank == 0 ) {

		
        reader readf(5);
        readf.readInputData(rows,columns,inputData,parser.readStringOption("-i"));
        t5 = get_timestamp();
        
	    out<<(t5 - t0) / 1000000.0L<<endl;
	    cout<<"read input data" <<(t5 - t0) / 1000000.0L<<endl;
        MatrixXd mat = inputData.transpose();
       
        totalRows = inputData.rows();
        totalCols = inputData.cols();
        avg.resize(1,rows);
       
        avg =mat.colwise().mean();
        
        centered= mat.rowwise() - mat.colwise().mean();
        
        Eigen::JacobiSVD<MatrixXd> svd(centered, ComputeThinU | ComputeThinV);
        VectorXd d = svd.singularValues();
        MatrixXd u = svd.matrixV();
        Util::changeSign(u);
      
        spheringMatrix.resize(rows,rows);
       
        int g=0;
        for(g=0;g<d.size();g++){
            VectorXd v1 = u.col(g);
            v1= v1/d[g];
            u.col(g).array()=v1;
        }
       
        spheringMatrix=u.transpose();
     
        MatrixXd X2 = (centered*u).transpose();
     
        w_init=Util::getRandomMatrix(rows);
       
        inputData=X2*sqrt(columns);
        cout<<"whiteting mat ------------------------------------------------------------>"<<endl;
        cout << w_init << endl;

		t1 = get_timestamp();
		out<<(t1 - t5) / 1000000.0L<<endl;
	    //cout<<"pre data" <<(t1 - t5) / 1000000.0L<<endl;

    }
	
    PartitionData::distributeMatrix( inputData );
    PartitionData::broadcastMatrix( w_init );
    //PartitionData::broadcastNumbersOfThreads( maxNumberOfThreads, threadsDistribution );
    PartitionData::broadcastNumbersOfThreads( maxNumberOfThreads, parser.readIntOption( "-th" ) );

    //maxNumberOfThreads=numberOfThreadsFinal;
    MPI_Barrier(MPI_COMM_WORLD);
    numberOfThreadsFinal=maxNumberOfThreads;

    //weight_mat=w_init;
    //cerr<<"-----------------------------------------"<<endl;
    //cout<<w_init<<endl;
    gFunction.resize(numberOfChannels, numberOfChannels);
    gFunction<< MatrixXd::Zero(numberOfChannels,numberOfChannels);
    derivativeOfGFunction.resize(numberOfChannels, numberOfChannels);
    derivativeOfGFunction<< MatrixXd::Zero(numberOfChannels ,numberOfChannels);
    weight_mat.resize(numberOfChannels, numberOfChannels);
    cerr<<"---------number of threads -------------"<<maxNumberOfThreads<<"------------------"<<endl;

    weight_mat=w_init;
    FastICASettings<T> settings;
    loadSettings<T>( parser, settings );

    bool verbose = false;
    if ( rank == 0 )verbose = true;

    int rankOfThread=0;
    std::map<int,int> result;
    std::vector<int> resultPair;
    Util::splitData(numberOfChannels,inputData.cols(),maxNumberOfThreads,result);
    for(auto elem : result)
    {
        resultPair.push_back(elem.first);
        resultPair.push_back(elem.second);
        cerr<<elem.first<<"----"<<elem.second<<endl;

    }


    std::thread process_thread(data_processing);
    int getPart=0;
    int startP=0;
    //process_thread.join();

    for (int f=0;f<maxNumberOfThreads;f++ )
    {
        std::cerr<<rankOfThread<<endl;
        MatrixXd partOfInput=inputData.block(0, resultPair[getPart],numberOfChannels,resultPair[getPart+1]-resultPair[getPart]+1);
        FastICAold<T>* fastIca = new FastICAold<T>( settings,partOfInput, 123456UL, verbose, numberOfSamples,
                                                    numberOfChannels ,w_init, startP, 0,rankOfThread);
        getPart+=2;
        thredArray.push_back(std::thread (&FastICAold<T>::runFastICA,fastIca));
        checkedEvent.push_back(0);
        checkeWeight.push_back(0);
        startP=startP+numberOfChannels;
        rankOfThread++;
    }

    int joinedThreds=0;
   // MPI_Barrier(MPI_COMM_WORLD);
    for(auto &t:thredArray)
    {

        if(t.joinable()){
            t.join();
            cerr<<"Threads "<<joinedThreds<<"joined"<<endl;
            joinedThreds++;
        }

    }
    process_thread.join();

   cout<<"All threads have been joined"<<endl;

    //fastIca.setWhiteningMatrix( w_init );
    
   // W=
            //fastIca.runFastICA();

    if ( parser.readFlagOption( "-o" ) )
         {
             PartitionData::gatherMatrix(inputData, totalRows, totalCols);
             if(rank==0){
                  cout<<"Gather matrix...."<<rank<<endl;
                  cout<<inputData<<endl;
             }
         }


    if ( rank == 0 ) {

            cout<<"Finally rank zero is going to write"<<endl;
           /* unmixedSignal = (W*spheringMatrix)*centered.transpose();
            S = unmixedSignal.transpose();
            MatrixXd unmixingMatrix;
            fastIca.mixingMatrix(unmixingMatrix, false);
            t2 = get_timestamp();
			out<<(t2 - t1) / 1000000.0L<<endl;
			cout<<"fastica " <<(t2 - t1) / 1000000.0L <<endl;
            if (parser.readFlagOption("-om")) {
				cout<<"inside om"<<endl;
                MatrixXd mixingMatrix(rows, columns);
                mixingMatrix = unmixingMatrix.inverse();
                MatrixWriter<T>::writeMatrix(
                        mixingMatrix, outputDataFormat, parser.readStringOption("-om"));

            }



              if ( parser.readFlagOption( "-ow") )
             {
                 MatrixWriter<T>::writeMatrix(S, outputDataFormat, parser.readStringOption( "-ow" ) );
                  t3 = get_timestamp();
				  out<<(t3 - t2) / 1000000.0L<<endl;
				  cout<<"write data" <<(t3 - t2) / 1000000.0L<<endl;
				  t4= get_timestamp();
				  out<<(t4 - t0) / 1000000.0L<<endl;
				   cout<<"total time" <<(t4 - t0) / 1000000.0L <<endl;
             }

             // sphering matrix
             if ( parser.readFlagOption( "-os" ))
             {
                 MatrixWriter<T>::writeMatrix(spheringMatrix, outputDataFormat, parser.readStringOption( "-os" ) );
             }

             if ( parser.readFlagOption( "-og" ) )
             {
                 MatrixXd weightMatrix(rows,columns);
                 fastIca.mixingMatrix( weightMatrix, false );
                 MatrixWriter<T>::writeMatrix(
                    weightMatrix, outputDataFormat, parser.readStringOption( "-og" ) );
             }

             // the separated data
             if ( parser.readFlagOption( "-o" ) )
             {
                 MatrixXd mixingMatrix(rows,rows);
                 VectorXd unmixedAverages( rows );
                 fastIca.mixingMatrix( mixingMatrix, false );
                 inputData=mixingMatrix*inputData;
                 fastIca.mixingMatrix( mixingMatrix );
                 MatrixXd Outputata=inputData.transpose();
                 MatrixWriter<T>::writeMatrix(
                         Outputata, outputDataFormat, parser.readStringOption( "-o" ) );
             }
             */
    }
    
   // }
    return 0;
}

int main( int argc, char** argv )
{
    FastICAProgram program( argc, argv );
    program.setupParser();
    cout.precision(10);
    int rank, size;

    int returnValue;
    

//        printf("This is double ");
        returnValue = mainWrapper<double>(
            program.parser(), program.MPIrank(), program.MPIsize() );


    program.shutdown();
    return returnValue;
}

void setupParser( SimpleCLParser& parser, int rank )
{
    bool isRequired = true;
    bool notRequired = false;

    
    string usage = "\nFastICA v. 0.0.1";

    parser.addUsageHeader( usage );


    // input file options
    parser.addOption( SimpleCLParser::STRING, "-i", "", isRequired );
    parser.addOption( SimpleCLParser::STRING, "-if", "", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-c", "number of input channels", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-s", "number of samples per channel", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-pc", "number of samples per channel", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-th", "number of samples per channel", isRequired );
    parser.addOption( SimpleCLParser::FLAG, "-single", "use single precision computations", notRequired );

    // output options
    parser.addOption( SimpleCLParser::STRING, "-o", "unmixed data file name", notRequired );
    parser.addOption( SimpleCLParser::STRING, "-om", "mixing matrix file name", notRequired );
    parser.addOption( SimpleCLParser::STRING, "-ow", "unmixing matrix file name", notRequired );
    parser.addOption( SimpleCLParser::STRING, "-os", "sphering matrix file name", notRequired );
    parser.addOption( SimpleCLParser::STRING, "-og", "weight matrix file name", notRequired );
    parser.addOption( SimpleCLParser::STRING, "-of", "output data format ( big, little, native, text )", notRequired );

    // sphering options
    parser.addOption( SimpleCLParser::FLAG, "-sphering", "compute sphering matrix from input data", notRequired );
    parser.addOption( SimpleCLParser::STRING, "-is", "load precomputed sphering matrix from file", notRequired );

    // FastICA options
    parser.addOption( SimpleCLParser::DOUBLE, "-t", "convergence tolerance", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-I", "maximum numer of iterations per channel", isRequired );
    parser.addOption( SimpleCLParser::STRING, "-C", "contrast function ( cubic, hyptan, gaussian )", isRequired );
    parser.addOption( SimpleCLParser::STRING, "-g", "weight matrix initialization type ( identity, random, user )", isRequired );
    parser.addOption( SimpleCLParser::STRING, "-ig", "weight matrix initialization file name", notRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-r", "maximum number of retries (using random restarts)", isRequired );

    // usage
    parser.addOption( SimpleCLParser::FLAG, "-h", "usage information", notRequired );

    parser.parse();
    if ( parser.readFlagOption( "-h" ) )
    {
        if ( rank == 0 )
        {
            cout << parser.usage() << endl;
        }
        exit( 0 );
    }

    parser.checkRequired();
}

DataFormat::DataFormat getOutputFormat( SimpleCLParser& parser )
{
    DataFormat::DataFormat format = DataFormat::TEXT;
    cout << "ICdata format type output  " << format<< " ";
    if ( parser.readFlagOption( "-of" ) )
    {
        string formatName = parser.readStringOption( "-of" );
        if ( formatName == "big" )
        {
            format = DataFormat::BIG;
        }
        else if ( formatName == "little" )
        {
            format = DataFormat::LITTLE;
        }
        else if ( formatName == "native" )
        {
            format = DataFormat::NATIVE;
        }
        else if ( formatName == "text" )
        {
            format = DataFormat::TEXT;
        }
        else 
        {
            cerr << "Invalid output type. Defaulting to text" << endl;
        }
    }
    return format;
}

template <class T>
void loadSettings( SimpleCLParser& parser, FastICASettings<T>& settings)
{
    settings.convergenceTolerance = parser.readDoubleOption( "-t" );
    settings.maximumIterations = parser.readIntOption( "-I" );
    string cf = parser.readStringOption( "-C" );
    if ( cf == "cubic" )
    {
        settings.contrastFunction = ContrastFunction<T>::CUBIC;
    }
    else if ( cf == "hyptan" )
    {
        settings.contrastFunction = ContrastFunction<T>::HYPERBOLIC_TAN;
    }
    else if ( cf == "gaussian" )
    {
        settings.contrastFunction = ContrastFunction<T>::GAUSSIAN;
    }
    else
    {
        cerr << "unknown constrast function: " << cf << endl;
        cerr << "exiting" << endl;
        exit( 1 );
    }
    string initType = parser.readStringOption( "-g" );
    if ( initType == "identity" )
    {
        settings.initializationType = FastICASettings<T>::IDENTITY;
    }
    else if ( initType == "random" )
    {
        settings.initializationType = FastICASettings<T>::RANDOM;
    }
    else if ( initType == "user" )
    {
        settings.initializationType = FastICASettings<T>::USER_SPECIFIED;
    }
    else
    {
        cerr << "unknown initialization type " << initType << endl;
        cerr << "exiting" << endl;
        exit( 1 );
    }
    settings.userInitializationFile = parser.readStringOption( "-ig" );
    // TODO make this smarter
    settings.userInitializationFileFormat = DataFormat::BIG;
    settings.maximumRetries = parser.readIntOption( "-r" );
}

