#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>


int main(int argc, char** argv)
{
 
  bem2d::freqtype k={160,0}; // Wavenumber
        int n=(int) 3000;     // Size of the linear system


        clock_t start, finish;
        double time;
        start=clock();
 

#ifdef BEM2DMPI
        MPI_Init(&argc, &argv);


        int nprow=4; // Number of rows in process grid
        int npcol=2; // Number of columns in process grid
        int mb=100;  // Row Block size
        int nb=100;  // Column Block size
        bem2d::BlacsSystem* b=bem2d::BlacsSystem::Initialize(nprow,npcol,mb,nb);

	// Exit if Context could not be created or process does not belong to context

        if (!b) {
                std::cout <<  "Could not create Blacs context" << std::endl;
                MPI_Finalize();
                exit(1);
        }
        if ((b->get_myrow()==-1)&&(b->get_mycol()==-1)) {
                MPI_Finalize();
                exit(0);
        }
#endif

        bem2d::pCurve cobj(new bem2d::Circle);
        bem2d::AnalyticCurve circle(n,cobj);
        bem2d::pGeometry pgeom=circle.GetGeometry();

        bem2d::PolBasis::AddBasis(0,pgeom); // Add constant basis functions


	// Set up direction of incoming wave.
        // Needed as input parameter for soundsoft scattering class

        bem2d::Point direction=bem2d::normalize(bem2d::Point(1,0));
        bem2d::PlaneWave pw(direction,k);
        bem2d::CombinedPlaneWave cpw(direction,k);

        bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
        soundsoft.SetQuadOption(5,5,0.15);
        soundsoft.set_polygons(pgeom);
        soundsoft.set_plotInterior();

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Discretize System" << std::endl;
	}
#else
	std::cout << "Discretize System" << std::endl;
#endif

        soundsoft.DiscretizeMatrix();

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Compute norm and condition number of operator" << std::endl;
	}
#else
	std::cout << "Compute norm and condition number of operator" << std::endl;
#endif

        double norm,cond;
        soundsoft.NormCond(norm,cond);

        finish=clock();
        time=(double(finish)-double(start))/CLOCKS_PER_SEC/60;


#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Norm: " << norm << " Norm of Inverse: " << cond/norm << " Condition Nr.: " << cond <<  std::endl;
	  std::cout << "Overalll time (minutes): " << time << std::endl;
	}
#else
	  std::cout << "Norm: " << norm << " Norm of Inverse: " << cond/norm << " Condition Nr.: " << cond <<  std::endl;
	  std::cout << "Overalll time (minutes): " << time << std::endl;

#endif

#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

