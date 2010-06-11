############################################
#  Paths to Include Files
############################################

# GNU Scientific Library
# files expected in $(GSLINC)/gsl
GSLINC	= /opt/local/include


############################################
#  Paths to Libraries
############################################

# GNU Scientific Library
GSLLIB	= /opt/local/lib
 
# BLAS
BLASLIB =/Users/tbetcke/lib

# LAPACK
LAPACKLIB =/Users/tbetcke/lib

# SCLAPACK
SCLAPACKLIB= /Users/tbetcke/lib


############################################
# Libraries
############################################

# GNU Scientific Library
GSL = -lgsl

# Lapack
LAPACK = -llapack

# Blas
BLAS   = -lgoto2

# OpenMP
OMP= -lgomp

# ScaLapack
SCLAPACK= -lscalapack -lblacsC -lblacs -lblacsC -lmpi_f77 -lgfortran

############################################
# Compiler and Flags
############################################
CPP	= g++
F77     = gfortran
FFLAGS  = -m64 -O
INCLUDES = -I$(GSLINC)
CFLAGS  += $(INCLUDES)
DEBUGFLAGS = -g -save-temps
OPTFLAGS = -O3
LDFLAGS = -L$(GSLLIB) -L$(LAPACKLIB) -L$(BLASLIB)
LDFLAGS += $(GSL) $(LAPACK) $(BLAS) $(OMP) -lgfortran


