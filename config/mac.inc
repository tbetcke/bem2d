############################################
# Specify Directories and output name
############################################

SRCDIR	= lib

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

# SCLAPACK
SCLAPACKLIB= /Users/tbetcke/libs/scalapack


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

# ScaLapack
SCLAPACK= -lscalapack -lblacsC -lblacs -lblacsC -lgfortran -lmpi_f77 -lgfortran

############################################
# Compiler and Flags
############################################
CPP	= mpic++
CFLAGS  = -I$(GSLINC) -I$(ATLASINC) -I/opt/local/include
DEBUGFLAGS = -g
OPTFLAGS = -O3

# For open source atlas

# LDFLAGS	= -L$(GSLLIB) -L$(ATLASLIB) $(GSL) $(LAPACK) $(BLAS)

# Accelerate Framework

  LDFLAGS = -framework Accelerate -L$(SCLAPACKLIB) -L/opt/local/lib -L$(GSLLIB) $(GSL) $(SCLAPACK)

# For MKL

# LDFLAGS = -L$(GSLLIB) $(GSL) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
