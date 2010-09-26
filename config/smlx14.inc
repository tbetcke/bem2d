############################################
#  Paths to Include Files
############################################

# GNU Scientific Library
# files expected in $(GSLINC)/gsl
GSLINC	= /usr/include


############################################
#  Paths to Libraries
############################################

# GNU Scientific Library
GSLLIB	= /usr/lib
 
# BLAS
BLASLIB =/apps/intel/Compiler/11.1/073/mkl/lib/em64t

# LAPACK
LAPACKLIB =/apps/intel/Compiler/11.1/073/mkl/lib/em64t

# SCLAPACK
SCLAPACKLIB=/apps/intel/Compiler/11.1/073/mkl/lib/em64t


############################################
# Libraries
############################################

# GNU Scientific Library
GSL = -lgsl

# Lapack
LAPACK = -lmkl_solver_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

# Blas
BLAS   = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

# OpenMP
OMP= -liomp5

# ScaLapack

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


