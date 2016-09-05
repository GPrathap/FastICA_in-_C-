#ifndef NICFASTICA_H_UNIVERSITY_OF_OREGON_NIC
#define NICFASTICA_H_UNIVERSITY_OF_OREGON_NIC

#include "../src/SimpleCLParser.h"
#include "../src/MatrixWriter.h"
#include "../src/FastICASettings.h"
#include "../src/DataFormat.h"

void setupParser( SimpleCLParser& parser, int rank );
DataFormat::DataFormat getInputFormat( SimpleCLParser& parser );
DataFormat::DataFormat getOutputFormat( SimpleCLParser& parser );

template <class T>
void loadSettings( SimpleCLParser& parser, FastICASettings<T>& settings );

#endif
// NICFASTICA_H_UNIVERSITY_OF_OREGON_NIC
