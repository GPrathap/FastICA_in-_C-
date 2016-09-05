

#include <stdio.h>
#include "hFastICA.h"


#include "../src/ProgramTools.h"
#include "../src/FastICAold.h"
#include "../src/reader.h"
#include "../src/DataFormat.h"
#include "../src/FastICAProgram.h"
#include <iostream>
#include <Eigen/Dense>
#include <sys/time.h>
#include <fstream> 
#include <time.h>
#include <vector>
#include <iostream>
#include "../src/Util.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ComputeFullV;
using Eigen::ComputeThinV;
using Eigen::ComputeThinU;
using Eigen::JacobiSVD;

typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp ()
{
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int numberOfThreadsFinal;
int IsPrepare=0;
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
bool isTerminate =false;
int weightIsReady = true;

void data_processing()
{

    double lim;
    while(true)
    {

        std::cout<<"";
        while(IsPrepare==numberOfThreadsFinal)
        {
            std::lock_guard<std::mutex> lk(lockWorker);

            MatrixXd input_to_symmetric=gFunction;
            MatrixXd input_to_symmetricDerivative=derivativeOfGFunction;
            input_to_symmetric=input_to_symmetric-input_to_symmetricDerivative;
            MatrixXd W1 = Util::sym_decorrelation(input_to_symmetric);
            lim = Util::estimateConvergence(W1,weight_mat );
            std::cerr<<"Lim:\t"<<lim<<"\tNumber of Iter: "<<numberOfIteratios<<std::endl;
            if((lim<convergenceTolerance)){
                isTerminate=true;
                std::fill(checkeWeight.begin(), checkeWeight.end(), 0);
                data_sends.notify_all();
                data_prepare.notify_all();
                cout<<"Solution is converged ..............................................."<<endl;
                goto stop;
            }else if(numberOfIteratios>=maxNumberOfIterations){
                cout<<"Solution is not converged ..............................................."<<endl;
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

template <class T> int mainWrapper( SimpleCLParser& parser)
{
    timestamp_t t0 = get_timestamp();
    timestamp_t t1,t2,t3,t4,t5;

    numberOfChannels = parser.readIntOption( "-c" );
    numberOfSamples  = parser.readIntOption( "-s" );
    int numberOfThreads= parser.readIntOption( "-th" );
    maxNumberOfIterations= parser.readIntOption( "-I" );
    convergenceTolerance= parser.readDoubleOption( "-t" );
    numberOfThreadsFinal=numberOfThreads;
    MatrixXd inputDataShared(numberOfChannels,numberOfSamples);

    DataFormat::DataFormat outputDataFormat = getOutputFormat( parser );

    MatrixXd inputData(numberOfChannels,numberOfSamples);
    MatrixXd spheringMatrix(numberOfChannels, numberOfChannels);

    gFunction.resize(numberOfChannels, numberOfChannels);
    gFunction<< MatrixXd::Zero(numberOfChannels,numberOfChannels);
    derivativeOfGFunction.resize(numberOfChannels, numberOfChannels);
    derivativeOfGFunction<< MatrixXd::Zero(numberOfChannels ,numberOfChannels);
    weight_mat.resize(numberOfChannels, numberOfChannels);

    MatrixXd avg ;
    MatrixXd centered;
    MatrixXd w_init(numberOfChannels, numberOfChannels);	//random matrix
    MatrixXd W(numberOfChannels, numberOfChannels);	//result
    MatrixXd unmixedSignal;	//unmixing X using W
    MatrixXd S;


	std::ofstream out("testResult.txt", std::ios::app);  

    reader::readInputData(numberOfChannels ,numberOfSamples,inputData,parser.readStringOption("-i"));

    t1 = get_timestamp();
    out<<(t1 - t0) / 1000000.0L<<endl;
    //cout<<"read input data" <<(t5 - t0) / 1000000.0L<<endl;

    MatrixXd mat = inputData.transpose();

    avg.resize(1,numberOfChannels);
    avg =mat.colwise().mean();

    centered= mat.rowwise() - mat.colwise().mean();

    Eigen::JacobiSVD<MatrixXd> svd(centered, ComputeThinU | ComputeThinV);
    VectorXd d = svd.singularValues();
    MatrixXd u = svd.matrixV();
    Util::changeSign(u);

    spheringMatrix.resize(numberOfChannels, numberOfChannels);

    int g=0;
    for(g=0;g<d.size();g++){
        VectorXd v1 = u.col(g);
        v1= v1/d[g];
        u.col(g).array()=v1;
    }

    spheringMatrix=u.transpose();
    //todo
    MatrixXd X2 = (centered*u).transpose();
    w_init=Util::getRandomMatrix(numberOfChannels);

    inputData=X2*sqrt(numberOfSamples);
   // cout<<inputData.transpose()<<endl;
    weight_mat=w_init;
    t2 = get_timestamp();
    out<<(t2 - t1) / 1000000.0L<<endl;
    cout<<"pre-processing time :" <<(t5 - t0) / 1000000.0L<<endl;
    std::map<int,int> result;
    std::vector<int> resultPair;
    Util::splitData(numberOfChannels,numberOfSamples,numberOfThreads,result);
    for(auto elem : result)
    {
        resultPair.push_back(elem.first);
        resultPair.push_back(elem.second);
    }
    FastICASettings<T> settings;
    loadSettings<T>( parser, settings );
    bool verbose = false;
    int startP=0;

    int rankOfThread=0;

    std::thread process_thread(data_processing);
    int getPart=0;
    for (int f=0;f<numberOfThreads;f++ )
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

    for(auto &t:thredArray)
    {
        t.join();
    }
    process_thread.join();
    unmixedSignal = (weight_mat*spheringMatrix)*centered.transpose();
    S = unmixedSignal.transpose();
    t3 = get_timestamp();
    out<<(t3 - t2) / 1000000.0L<<endl;
    cout<<"fastica " <<(t3 - t2) / 1000000.0L <<endl;
     if ( parser.readFlagOption( "-ow") )
     {
           //  MatrixWriter<T>::writeMatrix(S, outputDataFormat, parser.readStringOption( "-ow" ) );
             t4 = get_timestamp();
             out<<45.789<<endl;
             cout<<"write data" <<(t4 - t3) / 1000000.0L<<endl;
             t5= get_timestamp();
             out<<456.78<<endl;
             cout<<"total time" <<(t5 - t0) / 1000000.0L <<endl;
     }

    return 0;
}

int main( int argc, char** argv )
{
    FastICAProgram program( argc, argv );
    program.setupParser();
    cout.precision(10);
    int returnValue;
    returnValue = mainWrapper<double>(program.parser());
    return returnValue;
}

void setupParser( SimpleCLParser& parser, int rank )
{
    bool isRequired = true;
    bool notRequired = false;
    string usage = "\nFastICA v. 0.0.1";
    parser.addUsageHeader( usage );
    parser.addOption( SimpleCLParser::STRING, "-i", "", isRequired );
    parser.addOption( SimpleCLParser::STRING, "-if", "", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-c", "number of input channels", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-s", "number of samples per channel", isRequired );
    parser.addOption( SimpleCLParser::STRING, "-ow", "unmixing matrix file name", notRequired );
    parser.addOption( SimpleCLParser::FLAG, "-sphering", "compute sphering matrix from input data", notRequired );
    parser.addOption( SimpleCLParser::STRING, "-is", "load precomputed sphering matrix from file", notRequired );
    parser.addOption( SimpleCLParser::DOUBLE, "-t", "convergence tolerance", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-I", "maximum numer of iterations per channel", isRequired );
    parser.addOption( SimpleCLParser::STRING, "-C", "contrast function ( cubic, hyptan, gaussian )", isRequired );
    parser.addOption( SimpleCLParser::STRING, "-g", "weight matrix initialization type ( identity, random, user )", isRequired );
    parser.addOption( SimpleCLParser::INTEGER, "-r", "maximum number of retries (using random restarts)", isRequired );
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

