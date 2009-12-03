#include "bem2d_mpi.h"

namespace bem2d {
	
	BlacsSystem* BlacsSystem::pInstance_=0;
	
	BlacsSystem::BlacsSystem(int nprow, int npcol,int mb, int nb) :
	mb_(mb), nb_(nb)
	{
		
		// Initialize Scalapack
		sl_init_(&ictxt_,&nprow,&npcol);
		
		// Fill the context information
		MPI_Comm_size(MPI_COMM_WORLD,&mpiprocs_);
		MPI_Comm_rank(MPI_COMM_WORLD,&mpirank_);
		
		Cblacs_gridinfo(ictxt_,&nprow_,&npcol_,&myrow_,&mycol_);
		
	}
	
		
		
}