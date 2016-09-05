include ../system/make.common
include ../system/make.$(SYSTEM)

# Files
SOURCES= \
    TypeIdentifier.cpp      reader.cpp \
    SimpleCLParser.cpp      Util.cpp\
    DataWriter.cpp          ContrastFunctions.cpp\
    MatrixWriter.cpp        FastICAold.cpp \
    FastICASettings.cpp \
    ProgramTools.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# directories
OBJDIR = ../obj
SYSOBJDIR = $(OBJDIR)/$(SYSTEM)
BINDIR = ../bin

# target objects
SYSOBJECTS = $(addprefix $(SYSOBJDIR)/,$(OBJECTS))

# Defines
DEFINES = $(LAPACK_DEFINE) $(TAU_DEFINE) $(TAU_DEFS) $(MPI_DEFINE) \
    -DINSTANTIATE_TEMPLATES $(OMP_DEFINE) $(OMP_FLAGS) $(MATRIX_DEFINE)

# Includes
INCLUDE = $(MATRIX_INCLUDE) $(TAU_INCLUDE) $(MPI_INCLUDE)

# Link flags
LINK = $(MATRIX_LIBS) $(MATRIX_LIB_PATH) $(TAU_LIBS) $(MPI_LIBS) \
    $(OMP_LIBS) $(CLIBS)

# The main target. No other targets should appear before this one
all: $(OBJDIR) $(SYSOBJDIR) $(SYSOBJECTS)

$(OBJDIR):
	mkdir $(OBJDIR)

$(SYSOBJDIR):
	mkdir $(SYSOBJDIR)

$(SYSOBJDIR)/%.o: %.cpp %.h
	$(CXX) $(DEFINES) $(INCLUDE) $(CXX_FLAGS) -c $< -o $@

tools:
	cd ../tools; make SYSTEM=$(SYSTEM) all

selftests:
	cd ../test; make SYSTEM=$(SYSTEM) selftests

clean:
	$(RM) $(SYSOBJDIR) $(SYSTEM_TOOLS) $(TOOLS)
