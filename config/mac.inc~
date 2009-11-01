############################################
# Specify Directories and output name
############################################

SRCDIR	= libs

############################################
#  Paths to Include Files
############################################

# GNU Scientific Library
GSLINC	= /opt/local/include/gsl
 
# ATLAS Include
ATLASINC = /usr/include/atlas

############################################
#  Paths to Libraries
############################################

# GNU Scientific Library
GSLLIB	= /opt/local/lib
 
# ATLAS
ATLASLIB = /usr/lib64/atlas

############################################
# Libraries
############################################

# GNU Scientific Library
GSL = -lgsl

# Lapack
LAPACK = -llapack

# Blas
BLAS   = -lptf77blas -lptcblas -latlas -lpthread -lgfortran

# OpenMP
OMP= -lgomp


############################################
# Compiler and Flags
############################################
CPP	= g++
CFLAGS  = -I$(GSLINC) -I$(ATLASINC) -I/opt/local/include -fopenmp
DEBUGFLAGS = -g
OPTFLAGS = -O3

# For open source atlas

# LDFLAGS	= -L$(GSLLIB) -L$(ATLASLIB) $(GSL) $(LAPACK) $(BLAS)

# Accelerate Framework

  LDFLAGS = -framework Accelerate -L/opt/local/lib -L$(GSLLIB) $(GSL) $(OMP)

# For MKL

# LDFLAGS = -L$(GSLLIB) $(GSL) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
