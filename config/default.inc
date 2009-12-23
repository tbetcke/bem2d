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
BLASLIB =

# LAPACK
LAPACKLIB =

# SCLAPACK
SCLAPACKLIB= /Users/tbetcke/libs/scalapack


############################################
# Libraries
############################################

# GNU Scientific Library
GSL = -lgsl

# Lapack
LAPACK = -framework Accelerate

# Blas
BLAS   = 

# OpenMP
OMP= -lgomp

# ScaLapack
SCLAPACK= -lscalapack -lblacsC -lblacs -lblacsC -lmpi_f77 -lgfortran

############################################
# Compiler and Flags
############################################
CPP	= mpic++
INCLUDES = -I$(GSLINC)
CFLAGS  += $(INCLUDES)
DEBUGFLAGS = -g -save-temps
OPTFLAGS = -O3
LDFLAGS = -L$(GSLLIB) -L$(SCLAPACKLIB) -L$(LAPACKLIB) -L$(BLAS)
LDFLAGS += $(GSL) $(SCLAPACK) $(LAPACK) $(BLAS) $(OMP)


