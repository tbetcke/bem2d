#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include "mpi.h"



int main(int argc, char** argv){
	// File to test different things
	
	int myrank_mpi, nprocs_mpi;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank_mpi);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs_mpi);
	
	
	int nprow=2; int npcol=2;
	int myrow, mycol;
	int ictxt;
	char order='R';
	
	sl_init_(&ictxt, &nprow, &npcol);
	
	std::cout << ictxt << std::endl;
	
	Cblacs_gridinfo(ictxt,&nprow,&npcol,&myrow,&mycol);

	std::cout << "Rank: " << myrank_mpi << " Size: " << nprocs_mpi << " Row: " << myrow << " Col: "<< mycol <<std::endl;

	
	Cblacs_gridexit(0);	
	MPI_Finalize(); exit(0);
	
	int n=500;
	
	
	
	
	bem2d::Circle cobj;
	bem2d::AnalyticCurve<bem2d::Circle> circle(n,cobj);
	bem2d::pGeometry pgeom=circle.GetGeometry();
	
	
	
	/*
	bem2d::DiskShapePiecewiseConst circle(n,1.0);
	*/
	
	
	bem2d::freqtype k=10;
	bem2d::PolBasis::AddBasis(0,pgeom);
	std::cout << pgeom->size() << std::endl; 
	
	
	bem2d::Point direction=bem2d::normalize(bem2d::Point(0,-1));
	bem2d::PlaneWave pw(direction,k);
	bem2d::CombinedPlaneWave cpw(direction,k);
	
	
	std::cout << pw.k() << std::endl;

	
	bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
	soundsoft.SetQuadOption(5,5,0.15);
	soundsoft.set_polygons(pgeom);
	soundsoft.set_plotInterior();

	soundsoft.Discretize();
	
	std::cout << "Solve system" << std::endl;

	soundsoft.Solve();

	std::cout << "System solved" << std::endl;
	
	std::cout << "Condition Number: " << soundsoft.L2Condition() << std::endl;
	
	int xpts=100; int ypts=100;
	bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,3,-2,3,"trapping"));
	soundsoft.SetOutput(pout);
	//soundsoft.WriteAll();
	
	bem2d::WriteMatrix("/Users/tbetcke/svn/numerical_coercivity/matlab/diskmatrix10",soundsoft.GetMatrix());	
	bem2d::WriteMatrix("/Users/tbetcke/svn/numerical_coercivity/matlab/iddiskmatrix10",soundsoft.GetIdent());
		 
	 
	return 0;
}
