############################################
# Specify Directories and output name
############################################

SRCDIR	= lib

############################################
#  Paths to Include Files
############################################

# GNU Scientific Library
GSLINC	= /usr/include/gsl
 
# ATLAS Include
ATLASINC = /usr/local/atlas/include

############################################
#  Paths to Libraries
############################################

# GNU Scientific Library
GSLLIB	= /usr/local/lib
 
# ATLAS
ATLASLIB = /usr/local/atlas/lib

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
CFLAGS  = -I$(GSL) -I$(ATLASINC)  -fopenmp
DEBUGFLAGS = -g
OPTFLAGS = -O3

# For open source atlas

LDFLAGS	= -L$(GSLLIB) -L$(ATLASLIB) $(GSL) $(LAPACK) $(BLAS) $(OMP)

# For MKL

# LDFLAGS = -L$(GSLLIB) $(GSL) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
