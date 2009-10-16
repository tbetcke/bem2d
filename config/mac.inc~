############################################
# Specify Directories and output name
############################################

SRCDIR	= libs

############################################
#  Paths to Include Files
############################################

# GNU Scientific Library
GSLINC	= /usr/include/gsl
 
# ATLAS Include
ATLASINC = /usr/include/atlas

############################################
#  Paths to Libraries
############################################

# GNU Scientific Library
GSLLIB	= /usr/lib64
 
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



############################################
# Compiler and Flags
############################################
CPP	= g++
CFLAGS  = -I$(ATLASINC)
DEBUGFLAGS = -g
OPTFLAGS = -O3

# For open source atlas

LDFLAGS	= -L$(GSLLIB) -L$(ATLASLIB) $(GSL) $(LAPACK) $(BLAS)

# For MKL

# LDFLAGS = -L$(GSLLIB) $(GSL) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
