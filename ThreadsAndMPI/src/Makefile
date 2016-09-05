include ../system/make.common
include ../system/make.$(SYSTEM)

# Files
SOURCES= \
    TypeIdentifier.cpp      reader.cpp \
    SimpleCLParser.cpp      MPIHolder.cpp \
    DataWriter.cpp          ContrastFunctions.cpp\
    PartitionData.cpp       Util.cpp\
    MatrixWriter.cpp        FastICAold.cpp \
    FastICASettings.cpp \
    MPIWrapper.cpp \
    ProgramTools.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# directories
OBJDIR = ../obj
SYSOBJDIR = $(OBJDIR)/$(SYSTEM)
BINDIR = ../bin

# target objects
SYSOBJECTS = $(addprefix $(SYSOBJDIR)/,$(OBJECTS))

# Defines
DEFINES =  $(TAU_DEFINE) $(TAU_DEFS) $(MPI_DEFINE) \
    -DINSTANTIATE_TEMPLATES $(OMP_DEFINE) $(OMP_FLAGS) $(MATRIX_DEFINE)

# Includes
INCLUDE = $(MATRIX_INCLUDE) $(TAU_INCLUDE) $(MPI_INCLUDE)

# Link flags
LINK = $(MATRIX_LIBS) $(MATRIX_LIB_PATH) $(TAU_LIBS) $(MPI_LIBS) \
    $(OMP_LIBS)

# The main target. No other targets should appear before this one
all: $(OBJDIR) $(SYSOBJDIR) $(SYSOBJECTS)

$(OBJDIR):
	mkdir $(OBJDIR)

$(SYSOBJDIR):
	mkdir $(SYSOBJDIR)

$(SYSOBJDIR)/%.o: %.cpp %.h
	$(CXX) $(DEFINES) $(INCLUDE) $(CXX_FLAGS) -std=c++11 -c $< -o $@

tools:
	cd ../tools; make SYSTEM=$(SYSTEM) all

selftests:
	cd ../test; make SYSTEM=$(SYSTEM) selftests

clean:
	$(RM) $(SYSOBJDIR) $(SYSTEM_TOOLS) $(TOOLS)
