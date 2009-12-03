#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>



int main(int argc, char** argv){
	// File to test different things
	
	MPI_Init(&argc, &argv);
	
	int nprow=2; int npcol=2; int mb=3; int nb=10;
	bem2d::BlacsSystem* b=bem2d::BlacsSystem::Initialize(nprow,npcol,mb,nb);
	
	int myrow=b->get_myrow(); int mycol=b->get_mycol();
	 
	std::cout << "Rank: " << b->get_mpirank() << " Size: " << b->get_mpiprocs() << " Row: " << b->get_myrow() << " Col: "<< b->get_mycol() <<std::endl;
	bem2d::BlacsSystem::Release();
	
	int n=104; int np0=0;
	std::cout << bem2d::numroc_(&n,&mb,&myrow,&np0,&nprow) << std::endl;
	
	MPI_Finalize(); return 0; 
	/*
	int n=500;
	
	
	
	
	bem2d::Circle cobj;
	bem2d::AnalyticCurve<bem2d::Circle> circle(n,cobj);
	bem2d::pGeometry pgeom=circle.GetGeometry();
	
	
	
	
	bem2d::DiskShapePiecewiseConst circle(n,1.0);
	
	
	
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
	*/
}
