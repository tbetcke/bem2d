#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>



int main(int argc, char** argv){
	// File to test different things
	
	MPI_Init(&argc, &argv);
	
		
	int nprow=1; int npcol=2; int mb=33; int nb=33;
	bem2d::BlacsSystem* b=bem2d::BlacsSystem::Initialize(nprow,npcol,mb,nb);

	int myrow=b->get_myrow(); int mycol=b->get_mycol();

	/*	
	 
	std::cout << "Rank: " << b->get_mpirank() << " Size: " << b->get_mpiprocs() << " Row: " << b->get_myrow() << " Col: "<< b->get_mycol() <<std::endl;
	
	bem2d::Matrix A(5);
	
	// Fill the matrix
	
	if ((myrow==0)&(mycol==0)){
		(*A.data)[0]=1; (*A.data)[3]=2; (*A.data)[6]=5;
		(*A.data)[1]=6; (*A.data)[4]=9; (*A.data)[7]=7;
		(*A.data)[2]=3; (*A.data)[5]=2; (*A.data)[8]=8;
	}
	
	if ((myrow==0)&(mycol==1)){
		(*A.data)[0]=3; (*A.data)[3]=4;
		(*A.data)[1]=2; (*A.data)[4]=1;
		(*A.data)[2]=1; (*A.data)[5]=0;
	}
	
	if ((myrow==1)&(mycol==0)){
		(*A.data)[0]=3; (*A.data)[2]=2; (*A.data)[4]=5;
		(*A.data)[1]=9; (*A.data)[3]=0; (*A.data)[5]=4;
	}
		
	if ((myrow==1)&(mycol==1)){
		(*A.data)[0]=8; (*A.data)[2]=7;
		(*A.data)[1]=2; (*A.data)[3]=1;
	}

	bem2d::Matrix B(5,2);
	if ((myrow==0)&(mycol==0)){
		(*B.data)[0]=1; (*B.data)[3]=2;
	}
	
	bem2d::Matrix C(5,2);
	C=bem2d::SolveSystem(A,B);
	
	if ((myrow==0)&(mycol==0)){
		std::cout << "Result: " << (*C.data)[4] << std::endl;
		std::cout << "Global index for (3,1): " << b->L2gc(3) << std::endl;
	}
	
	//for (int i=0;i<9;i++) std::cout << A.desc[i] << " ";
	//std::cout << std::endl;
*/
	int n=500;
	
	bem2d::Circle cobj;
	bem2d::AnalyticCurve<bem2d::Circle> circle(n,cobj);
	bem2d::pGeometry pgeom=circle.GetGeometry();
	
	
	
	
	//bem2d::DiskShapePiecewiseConst circle(n,1.0);
	
	
	
	bem2d::freqtype k=5;
	bem2d::PolBasis::AddBasis(0,pgeom);
	std::cout << pgeom->size() << std::endl; 
	
	
	bem2d::Point direction=bem2d::normalize(bem2d::Point(1,0));
	bem2d::PlaneWave pw(direction,k);
	bem2d::CombinedPlaneWave cpw(direction,k);
	
	
	std::cout << pw.k() << std::endl;

	std::cout << "Initialize Soundsoft problem" << std::endl;
	
	
	bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
	soundsoft.SetQuadOption(5,5,0.15);
	soundsoft.set_polygons(pgeom);
	soundsoft.set_plotInterior();

	std::cout << "Discretize System" << std::endl;
	
	soundsoft.Discretize();
	
	std::cout << "Solve system" << std::endl;

	soundsoft.Solve();

	std::cout << "System solved" << std::endl;
	
	//std::cout << "Condition Number: " << soundsoft.L2Condition() << std::endl;
	
	// Test evaluation of solution
	
	//std::vector<bem2d::Point> p; p.push_back(bem2d::Point(1.5,0.0));
	//bem2d::pcvector sol=soundsoft.EvalSol(p);
	
	//bem2d::pMatrix A=soundsoft.GetMatrix();
	//std::cout << myrow << " " << mycol << " " << (*A->data)[250*250-1] << std::endl;

	
	
	
	int xpts=100; int ypts=100;
	bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,2,-2,2,"disk"));
	soundsoft.SetOutput(pout);
	
	soundsoft.WriteAll();
	
	//bem2d::WriteMatrix("/Users/tbetcke/svn/numerical_coercivity/matlab/diskmatrix10",soundsoft.GetMatrix());	
	//bem2d::WriteMatrix("/Users/tbetcke/svn/numerical_coercivity/matlab/iddiskmatrix10",soundsoft.GetIdent());
	
	bem2d::BlacsSystem::Release();	
	MPI_Finalize(); 
	
	 
}
