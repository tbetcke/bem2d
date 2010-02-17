#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>
#include<vector>
#include<fstream>
#include<sstream>

int main(int argc, char** argv)
{

  std::string file="/home/tbetcke/svn/numerical_coercivity/data/squarehconv";

  int numrange_n=50; // Number of discretization points for num. range.
  int computenorm=0; // Set to 1 to compute norm and condition number
  bem2d::freqtype k=1;
 
  std::vector<int> ppwvec;
  ppwvec.push_back(10);
  ppwvec.push_back(50);
  ppwvec.push_back(100);
  ppwvec.push_back(500);
  ppwvec.push_back(1000);
  ppwvec.push_back(1500);
  ppwvec.push_back(2000);
  //  ppwvec.push_back(1280);





        clock_t start, finish;
        double time;
        start=clock();
 

        std::vector<bem2d::Point> square;
        square.push_back(bem2d::Point(0,0));
        square.push_back(bem2d::Point(1,0));
        square.push_back(bem2d::Point(1,1));
        square.push_back(bem2d::Point(0,1));

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

	for (int j=0;j<ppwvec.size();j++){

	int ppw=ppwvec[j];
	double eta1=k; // Coupling between conj. double and single layer pot.
	bem2d::Polygon poly(square,ppw,k,0,0.15);
        bem2d::pGeometry pgeom=poly.GetGeometry();
	


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
	  std::cout << "Discretize Kernels with n=" << pgeom->size() << std::endl;
	}
#else
	std::cout << "Discretize Kernels with n=" << pgeom->size() << std::endl;
#endif



	bem2d::Matrix dsl=*(DiscreteKernel(*pgeom,quadopts,sl));
	bem2d::Matrix dcdl=*(DiscreteKernel(*pgeom,quadopts,cdl));
	bem2d::Matrix Id=*(EvalIdent(*pgeom, quadopts));
	bem2d::Matrix combined1=Id+2.0*dcdl-bem2d::complex(0,2.0)*eta1*dsl;
	
	combined1=bem2d::ChangeBasis(combined1,Id);

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Compute Eigenvalues and norm" << std::endl;
	}
#else
	std::cout << "Compute Eigenvalues and norm" << std::endl;
#endif

	double norm; double cond;
	
	bem2d::pcvector eigvals;
	bem2d::Eigenvalues(combined1,eigvals);

	if (computenorm){
	bem2d::L2NormCond(combined1,norm,cond);


	// Write out norm and condition number

#ifdef BEM2DMPI
	if (b->IsRoot()){
#endif

	  std::ostringstream osnormcond;
	  osnormcond << file << "_normcond_" << ppw; 
	  std::string s0=osnormcond.str();
	  std::ofstream o(s0.c_str());
	  o << norm << std::endl << cond << std::endl;
	  o.close();
#ifdef BEM2DMPI
	}
#endif
	}

	// Write out eigenvalues

#ifdef BEM2DMPI
	if (b->IsRoot()){
#endif

	  std::ostringstream os;
	  os << file << "_eig_" << ppw; 
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
	os2 << file << "_range_" << ppw;
	NumRange(combined1, numrange_n, os2.str());

	
	}
        finish=clock();
        time=(double(finish)-double(start))/CLOCKS_PER_SEC/60;

#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

