#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>
#include<vector>
#include<fstream>

int main(int argc, char** argv)
{

  int ppw=10;     // Point per wavelength
  std::string filename="cracknormcond10.txt";
  
 
  std::vector<bem2d::freqtype> freqs;
  freqs.push_back(5);
  freqs.push_back(10);
  freqs.push_back(20);
  freqs.push_back(40);
  freqs.push_back(80);
  freqs.push_back(160);
  freqs.push_back(320);
  freqs.push_back(640);
  freqs.push_back(1280);
  

  std::vector<double> norm_sl(freqs.size());
  std::vector<double> cond_sl(freqs.size());



        clock_t start, finish;
        double time;
        start=clock(); 

#ifdef BEM2DMPI
        MPI_Init(&argc, &argv);


        int nprow=4; // Number of rows in process grid
        int npcol=2; // Number of columns in process grid
        int mb=24;  // Row Block size
        int nb=24;  // Column Block size
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

	for (int j=0;j<freqs.size();j++){

	double k=(double)freqs[j];
	double eta1=k; // Coupling between conj. double and single layer pot.
	double eta2=cbrt(k*k);
	int n=(int)((int)k*ppw/2.0/bem2d::PI);

	// Create the geometry

	std::vector<bem2d::pElement> elements(n);
	for (int i=0;i<n;i++){
	  bem2d::Point p1=(1.0*i)/n*bem2d::Point(0,1);
	  bem2d::Point p2=(1.0*(i+1))/n*bem2d::Point(0,1);
	  elements[i]=bem2d::pElement(new bem2d::ConstElement(p1,p2,i));
	  elements[i]->set_next(i+1);
	  elements[i]->set_prev(i-1);
	}
	bem2d::pGeometry pgeom(new bem2d::Geometry(elements));
        bem2d::PolBasis::AddBasis(0,pgeom); // Add constant basis functions

	// Discretize the single and double layer potential

	bem2d::SingleLayer sl(k);

	bem2d::QuadOption quadopts;

	quadopts.L=3;
        quadopts.N=5;
        quadopts.sigma=0.15;

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Discretize Kernels with n=" << pgeom->size() << std::endl;
	}
#else
	std::cout << "Discretize Kernels with n=" << pgeom->size() << std::endl;
#endif



	bem2d::Matrix dsl=*(DiscreteKernel(*pgeom,quadopts,sl));
	bem2d::Matrix Id=*(EvalIdent(*pgeom, quadopts));
	
	dsl=2.0*bem2d::ChangeBasis(dsl,Id);

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Compute norms and condition numbers" << std::endl;
	}
#else
	std::cout << "Compute norms and condition numbers" << std::endl;
#endif
	

	bem2d::L2NormCond(dsl,norm_sl[j],cond_sl[j]);
	
	}
        finish=clock();
        time=(double(finish)-double(start))/CLOCKS_PER_SEC/60;


#ifdef BEM2DMPI
	if (b->IsRoot()){
#endif

	  std::ofstream out(filename.c_str());
 	  out << "Single Layer" << std::endl;

	  for (int j=0;j<freqs.size();j++){
	    out << "k=" << freqs[j] << " Norm: " << norm_sl[j] << " Norm of Inverse: " << cond_sl[j]/norm_sl[j] << " Condition Nr.: " << cond_sl[j] <<  std::endl;
	  }




	  out << "Overalll time (minutes): " << time << std::endl;
	  std::cout << "Overall time (minutes): " << time << std::endl;
	  out.close();
#ifdef BEM2DMPI	
}
#endif


#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

