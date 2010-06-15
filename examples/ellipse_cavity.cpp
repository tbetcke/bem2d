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

        int n=10;
		double a= 1;
        double ah=1.3;
        double b=0.5;
        double bh=0.6;
        double t0=0.01*bem2d::PI;
        double t1=acos(-a/ah*cos(t0));
        double alpha=bh*sin(t1)-b*sin(t0);
		bem2d::freqtype k={n*bem2d::PI/2/b,0};
	
        bem2d::Point p0(-a*cos(t0),b*sin(t0));
        bem2d::Point p1(-a*cos(t0),alpha+b*sin(t0));
        bem2d::Point p2(-a*cos(t0),-alpha-b*sin(t0));
        bem2d::Point p3(-a*cos(t0),-b*sin(t0));

        bem2d::pCurve cobj(new bem2d::Circle);
        bem2d::AnalyticCurve circle(n,cobj);
 
       
        bem2d::pCurve Arc(new bem2d::EllipseArc(a,b,bem2d::PI-t0,-bem2d::PI+t0));
        bem2d::AnalyticCurve ellipseArc(10,k,Arc,0,5);
        bem2d::pCurve Arc2(new bem2d::EllipseArc(ah,bh,-t1,t1));
        bem2d::AnalyticCurve ellipseArc2(10,k,Arc2,0,5);
        bem2d::Line l0(p1,p0,10,k,5);
        bem2d::Line l1(p3,p2,10,k,5);
	


        std::vector<bem2d::pGeometry> geoms;
		geoms.push_back(ellipseArc2.GetGeometry());
        geoms.push_back(l0.GetGeometry());
		geoms.push_back(ellipseArc.GetGeometry());
        geoms.push_back(l1.GetGeometry());
        bem2d::pGeometry pgeom(new bem2d::Geometry(geoms));


        std::vector<bem2d::pElement> elems =pgeom->elements();
        for (int i=1;i<elems.size();i++){
            std::cout << elems[i]->First() << " " << elems[i]->Last() << " " << elems[i]->First()-elems[i-1]->Last() << std::endl;
        }

        bem2d::PolBasis::AddBasis(2,pgeom);

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
        bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,3,-2,2,"ellipse_cavity"));
        soundsoft.SetOutput(pout);
        soundsoft.WriteAll();



#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif



        return 0;
}



