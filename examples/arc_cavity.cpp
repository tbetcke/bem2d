#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>


int main(int argc, char** argv)
{
	
#ifdef BEM2DMPI
	MPI_Init(&argc, &argv);
	
	
	int nprow=2;
	int npcol=1;
	int mb=100;
	int nb=100;
	bem2d::BlacsSystem* b=bem2d::BlacsSystem::Initialize(nprow,npcol,mb,nb);
	
	if (!b) {
		std::cout <<  "Could not create Blacs context" << std::endl;
		MPI_Finalize();
		exit(1);
	}
	if ((b->get_myrow()==-1)&&(b->get_mycol()==-1)) {
		MPI_Finalize();
		exit(0);
	}
	
	
	
	int myrow=b->get_myrow();
	int mycol=b->get_mycol();
#endif

	bem2d::freqtype k={20,0};
	int ppw=10;
	int L=0;
	double l=3./4;
	double ang=bem2d::PI/4;
	bem2d::Point p0(cos(ang),sin(ang));
	bem2d::Point p1(l*cos(ang),l*sin(ang));
	bem2d::Point p2(l*cos(ang),1);
	bem2d::Point p3(1.2,1);
	bem2d::Point p4(1.2,-0.2);
	bem2d::Point p5(l,-0.2);
	bem2d::Point p6(l,0);
	bem2d::Point p7(1,0);
	
	
	
	bem2d::pCurve Arc(new bem2d::EllipseArc(1,1,ang,0));
	bem2d::AnalyticCurve ellipseArc(ppw,k,Arc,0,L);
	bem2d::Line l0(p7,p6,ppw,k,L);
	bem2d::Line l1(p6,p5,ppw,k,L);
	bem2d::Line l2(p5,p4,ppw,k,L);
	bem2d::Line l3(p4,p3,ppw,k,L);
	bem2d::Line l4(p3,p2,ppw,k,L);
	bem2d::Line l5(p2,p1,ppw,k,L);
	bem2d::Line l6(p1,p0,ppw,k,L);
					  
					  
					  
	std::vector<bem2d::pGeometry> geoms;
	geoms.push_back(l0.GetGeometry());
	geoms.push_back(l1.GetGeometry());
	geoms.push_back(l2.GetGeometry());
	geoms.push_back(l3.GetGeometry());
	geoms.push_back(l4.GetGeometry());
	geoms.push_back(l5.GetGeometry());
	geoms.push_back(l6.GetGeometry());
	geoms.push_back(ellipseArc.GetGeometry());
					  
	
	
	
	bem2d::pGeometry pgeom(new bem2d::Geometry(geoms));
	/*
	 
	 std::vector<bem2d::pElement> elems =pgeom->elements();
	 for (int i=1;i<elems.size();i++){
	 std::cout << elems[i]->First() << " " << elems[i]->Last() << " " << elems[i]->First()-elems[i-1]->Last() << std::endl;
	 }
	 */
	bem2d::PolBasis::AddBasis(0,pgeom);
	
	bem2d::Point direction=bem2d::normalize(bem2d::Point(1,0));
	bem2d::PlaneWave pw(direction,k);
	bem2d::CombinedPlaneWave cpw(direction,k);
	
	bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
	soundsoft.SetQuadOption(5,5,0.15);
	soundsoft.set_polygons(pgeom);
	soundsoft.set_plotInterior();
	
	
	
	std::cout << "Discretizing with k=" << bem2d::complex(k.re,k.im) << " and n=" << pgeom->size() << std::endl;
	clock_t start, finish;
	double time;
	
	start=clock();
	soundsoft.Discretize();
	soundsoft.Solve();
	finish=clock();
	time=(double(finish)-double(start))/CLOCKS_PER_SEC;
	
#ifdef BEM2DMPI
	if (b->IsRoot()) {
#endif
		std::cout << "Computing time (minutes): " << time/60 << std::endl;
#ifdef BEM2DMPI
	}
#endif
	double cond;
	double norm;
	soundsoft.NormCond(norm,cond);
	
#ifdef BEM2DMPI
	if (b->IsRoot()) {
#endif
		std::cout << "Condition Number: " << cond << std::endl;
		std::cout << "Norm: " << norm << std::endl;
#ifdef BEM2DMPI
	}
#endif
	
	int xpts=200;
	int ypts=200;
	bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,0,2,-1,1.5,"arc_cavity"));
	soundsoft.SetOutput(pout);
	soundsoft.WriteAll();
	
	
	
#ifdef BEM2DMPI
	bem2d::BlacsSystem::Release();
	MPI_Finalize();
#endif
	
	
	
	return 0;
}






