//
// Created by Jobs on 6/28/2015.
//

#ifndef NEW_FASTICAPROGRAM_H
#define NEW_FASTICAPROGRAM_H


class FastICAProgram : public ProgramTools
{
public:
    FastICAProgram( int argc, char** argv )
            : ProgramTools( argc, argv )
    {
    }

    virtual void setupParser()
    {
        bool isRequired = true;
        bool notRequired = false;
        m_parser.addUsageHeader( "FASTICA" );


        // input file options
        m_parser.addOption( SimpleCLParser::STRING,
                            "-i", "input data file name",
                            isRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-if", "input data format ( big, little, native, text )",
                            isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-c", "number of input channels",
                            isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-s", "number of samples per channel",
                            isRequired );

        m_parser.addOption( SimpleCLParser::FLAG,
                            "-single", "use single precision computations",
                            notRequired );

        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-th", "Number of threads are required",
                            isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-pc", "Number of processes are required",
                            isRequired );
        // output options
        m_parser.addOption( SimpleCLParser::STRING,
                            "-o", "unmixed data file name", notRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-om", "mixing matrix file name", notRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-ow", "unmixing matrix file name", notRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-os", "sphering matrix file name", notRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-og", "weight matrix file name", notRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-of", "output data format ( big, little, native, text )",
                            notRequired );

        // sphering options
        m_parser.addOption( SimpleCLParser::FLAG,
                            "-sphering", "compute sphering matrix from input data", notRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-is", "load precomputed sphering matrix from file", notRequired );

        // FastICA options
        m_parser.addOption( SimpleCLParser::DOUBLE,
                            "-t", "convergence tolerance", isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-I", "maximum numer of iterations per channel", isRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-C", "contrast function ( cubic, hyptan, gaussian )", isRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-g", "weight matrix initialization type ( identity, random, user )", isRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-ig", "weight matrix initialization file name", notRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-r", "maximum number of retries (using random restarts)", isRequired );

        setup();

    }
};

#endif //NEW_FASTICAPROGRAM_H
