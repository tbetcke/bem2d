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


        int m=4;
        bem2d::freqtype k=5*m;
        int n=20*m;

        double a=0.31;
        double c=1;
        double l=c-a;



        std::vector<bem2d::Point> trapping;
        trapping.push_back(bem2d::Point(0,0));
        trapping.push_back(bem2d::Point(-c,0));
        trapping.push_back(bem2d::Point(-c,-l));
        trapping.push_back(bem2d::Point(l,-l));
        trapping.push_back(bem2d::Point(l,2*c-l));
        trapping.push_back(bem2d::Point(-c,2*c-l));
        trapping.push_back(bem2d::Point(-c,2*a));
        trapping.push_back(bem2d::Point(0,2*a));
        bem2d::Polygon poly(trapping,n);
        bem2d::pGeometry pgeom=poly.GetGeometry();



        bem2d::PolBasis::AddBasis(0,pgeom);


        bem2d::Point direction=bem2d::normalize(bem2d::Point(1,0));
        bem2d::PlaneWave pw(direction,k);
        bem2d::CombinedPlaneWave cpw(direction,k);

        bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
        soundsoft.SetQuadOption(5,5,0.15);
        soundsoft.set_polygons(pgeom);
        soundsoft.set_plotInterior();



        std::cout << "Discretizing with k=" << k << " and n=" << pgeom->size() << std::endl;
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
        bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,3,-2,2,"trapping"));
        soundsoft.SetOutput(pout);
        soundsoft.WriteAll();



#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif



        return 0;
}



