#include<iostream>
#include<fstream>
#include<iomanip>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>


int main(int argc, char** argv)
{
 
  bem2d::freqtype k={1,0};
 
  std::vector<int> ppwvec;
  ppwvec.push_back(5);
  ppwvec.push_back(10);
  ppwvec.push_back(50);

  ppwvec.push_back(100);
  ppwvec.push_back(500);
  ppwvec.push_back(1000);
  ppwvec.push_back(1500);
  ppwvec.push_back(2000);

        clock_t start, finish;
        double time;
        start=clock();
 
	bem2d::dvector result_real;
	bem2d::dvector result_imag;


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


	for (int j=0;j<ppwvec.size();j++){


	int ppw=ppwvec[j];
	double eta1=k.re; // Coupling between conj. double and single layer pot.
        bem2d::pCurve cobj(new bem2d::Circle);
	int n=(int)(cobj->Length()*k.re*ppw/2.0/bem2d::PI);
        bem2d::AnalyticCurve circle(n,cobj);
        bem2d::pGeometry pgeom=circle.GetGeometry();

        bem2d::LinBasis::AddBasis(pgeom); // Add constant basis functions



	// Set up direction of incoming wave.
        // Needed as input parameter for soundsoft scattering class

        bem2d::Point direction=bem2d::normalize(bem2d::Point(1,0));
        bem2d::PlaneWave pw(direction,k);
        bem2d::CombinedPlaneWave cpw(direction,k);

        bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
        soundsoft.SetQuadOption(7,30,0.15);
        soundsoft.set_polygons(pgeom);
        soundsoft.set_plotInterior();

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Discretize System" << std::endl;
	}
#else
	std::cout << "Discretize System" << std::endl;
#endif

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

	std::vector<bem2d::Point> pts;
	pts.push_back(bem2d::Point(-1.5,0));
	bem2d::pcvector output=soundsoft.EvalSol(pts);

	std::cout << (*output)[0] << std::endl;
	result_real.push_back(std::real((*output)[0]));
	result_imag.push_back(std::imag((*output)[0]));
	 

	}

	std::string fr="result_real";
        std::ofstream outr(fr.c_str());
	std::string fi="result_imag";
	std::ofstream outi(fi.c_str());

	for (int j=0;j<result_real.size();j++){
	  outr << std::setprecision(16) << result_real[j] << std::endl;
	  outi << std::setprecision(16) << result_imag[j] << std::endl;
	}

	outr.close();
	outi.close();


	/*
        int xpts=100;
        int ypts=100;
        bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,2,-2,2,"disk"));
        soundsoft.SetOutput(pout);
        soundsoft.WriteAll();
	*/



#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

