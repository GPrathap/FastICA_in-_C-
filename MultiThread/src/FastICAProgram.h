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
        m_parser.addOption( SimpleCLParser::STRING,
                            "-i", "input data file name", isRequired );
        m_parser.addOption( SimpleCLParser::STRING, "-if", "input data format ( big, little, native, text )",
                            isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-c", "number of input channels", isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-th", "number of threads", isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-s", "number of samples per channel", isRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-ow", "unmixing matrix file name", notRequired );
        m_parser.addOption( SimpleCLParser::FLAG,
                            "-sphering", "compute sphering matrix from input data", notRequired );
        m_parser.addOption( SimpleCLParser::DOUBLE,
                            "-t", "convergence tolerance", isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-I", "maximum numer of iterations per channel", isRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-C", "contrast function ( cubic, hyptan, gaussian )", isRequired );
        m_parser.addOption( SimpleCLParser::STRING,
                            "-g", "weight matrix initialization type ( identity, random, user )", isRequired );
        m_parser.addOption( SimpleCLParser::INTEGER,
                            "-r", "maximum number of retries (using random restarts)", isRequired );
        setup();
    }
};


#endif //NEW_FASTICAPROGRAM_H
