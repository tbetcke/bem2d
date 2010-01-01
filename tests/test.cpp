#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>



int main(int argc, char** argv)
{


#ifdef BEM2DMPI
        MPI_Init(&argc, &argv);


        int nprow=4;
        int npcol=2;
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
        int n=4000;


        bem2d::Circle cobj;
        bem2d::AnalyticCurve<bem2d::Circle> circle(n,cobj);
        bem2d::pGeometry pgeom=circle.GetGeometry();




        //bem2d::DiskShapePiecewiseConst circle(n,1.0);



        bem2d::freqtype k=50;
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

        soundsoft.DiscretizeMatrix();

        //std::cout << "Solve system" << std::endl;

        //soundsoft.Solve();

        //std::cout << "System solved" << std::endl;

        //std::cout << "Condition Number: " << soundsoft.L2Condition() << std::endl;

        // Test evaluation of solution

        //std::vector<bem2d::Point> p; p.push_back(bem2d::Point(1.5,0.0));
        //bem2d::pcvector sol=soundsoft.EvalSol(p);

        //bem2d::pMatrix A=soundsoft.GetMatrix();
        //std::cout << myrow << " " << mycol << " " << (*A->data)[250*250-1] << std::endl;




        //int xpts=100;
        //int ypts=100;
        //bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,2,-2,2,"disk"));
        //soundsoft.SetOutput(pout);

	std::cout << "Compute norm/cond" << std::endl;
        double norm,cond;
        soundsoft.NormCond(norm,cond);
        std::cout << norm << " " << cond << std::endl;

        //soundsoft.WriteAll();

        //bem2d::WriteMatrix("/Users/tbetcke/svn/numerical_coercivity/matlab/diskmatrix10",soundsoft.GetMatrix());
        //bem2d::WriteMatrix("/Users/tbetcke/svn/numerical_coercivity/matlab/iddiskmatrix10",soundsoft.GetIdent());


#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

