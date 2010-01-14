#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>
#include<vector>
#include<fstream>
#include<sstream>

int main(int argc, char** argv)
{

  int ppw=10;     // Point per wavelength
  std::string file1="diskrange";
  std::string file2="diskeig";

  int numrange_n=20; // Number of discretization points for num. range.
  
 
  std::vector<bem2d::freqtype> freqs;
  freqs.push_back(10);
  freqs.push_back(20);

  /*
  freqs.push_back(50);
  freqs.push_back(100);
  freqs.push_back(200);
  freqs.push_back(500);
  */

        clock_t start, finish;
        double time;
        start=clock();
 

#ifdef BEM2DMPI
        MPI_Init(&argc, &argv);


        int nprow=2; // Number of rows in process grid
        int npcol=1; // Number of columns in process grid
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
        bem2d::pCurve cobj(new bem2d::Circle);
	int n=(int)(cobj->Length()*k*ppw/2.0/bem2d::PI);
        bem2d::AnalyticCurve circle(n,cobj);
        bem2d::pGeometry pgeom=circle.GetGeometry();

        bem2d::PolBasis::AddBasis(0,pgeom); // Add constant basis functions


	// Discretize the single and double layer potential

	bem2d::SingleLayer sl(k);
	bem2d::ConjDoubleLayer cdl(k);

	bem2d::QuadOption quadopts;

	quadopts.L=3;
        quadopts.N=5;
        quadopts.sigma=0.15;

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Discretize Kernels with n=" << n << std::endl;
	}
#else
	std::cout << "Discretize Kernels with n=" << n << std::endl;
#endif



	bem2d::Matrix dsl=*(DiscreteKernel(*pgeom,quadopts,sl));
	bem2d::Matrix dcdl=*(DiscreteKernel(*pgeom,quadopts,cdl));
	bem2d::Matrix Id=*(EvalIdent(*pgeom, quadopts));
	bem2d::Matrix combined1=Id+2.0*dcdl-bem2d::complex(0,2.0)*eta1*dsl;
	
	combined1=bem2d::ChangeBasis(combined1,Id);

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Compute Eigenvalues" << std::endl;
	}
#else
	std::cout << "Compute Eigenvalues" << std::endl;
#endif
	
	bem2d::pcvector eigvals;
	bem2d::Eigenvalues(combined1,eigvals);

	// Write out eigenvalues

#ifdef BEM2DMPI
	if (b->IsRoot()){
#endif

	  std::ostringstream os;
	  os << file2 << k; 
	  std::string s=os.str();
	  std::ofstream o1(s.c_str());
	  for (int i=0;i<eigvals->size();i++) o1 << std::real((*eigvals)[i])
			     << " " << std::imag((*eigvals)[i]) << std::endl;
	  o1.close();
#ifdef BEM2DMPI
	}
#endif

	

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Compute Numerical Range" << std::endl;
	}
#else
	std::cout << "Compute Numerical Range" << std::endl;
#endif


	std::ostringstream os2;
	os2 << file1 << k;
	NumRange(combined1, numrange_n, os2.str());

	
	}
        finish=clock();
        time=(double(finish)-double(start))/CLOCKS_PER_SEC/60;

#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

