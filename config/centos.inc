############################################
#  Paths to Include Files
############################################

# GNU Scientific Library
# files expected in $(GSLINC)/gsl
GSLINC	= /usr/local/include


############################################
#  Paths to Libraries
############################################

# GNU Scientific Library
GSLLIB	= /usr/local/lib
 
# BLAS
BLASLIB =/usr/local/atlas/lib

# LAPACK
LAPACKLIB = /usr/local/atlas/lib

# SCLAPACK
SCLAPACKLIB= /usr/local/lib


############################################
# Libraries
############################################

# GNU Scientific Library
GSL = -lgsl

# Lapack
LAPACK = /usr/local/atlas/lib/liblapack.a

# Blas
BLAS   = -lptcblas -lptf77blas -latlas -lpthread

# OpenMP
OMP= -lgomp

# ScaLapack
SCLAPACK= -lscalapack -lblacsC -lblacs -lblacsC -lmpi_f77 -lgfortran

############################################
# Compiler and Flags
############################################
CPP	= mpic++
F77     = gfortran
INCLUDES = -I$(GSLINC)
CFLAGS  += $(INCLUDES)
DEBUGFLAGS = -g -save-temps
OPTFLAGS = -O3
LDFLAGS = -L$(GSLLIB) -L$(SCLAPACKLIB) -L$(LAPACKLIB) -L$(BLAS)
LDFLAGS += $(GSL) $(SCLAPACK) $(LAPACK) $(BLAS) $(OMP)


