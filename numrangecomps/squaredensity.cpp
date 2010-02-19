// Compute the eigenfunctions of the BEM operator on the square and write
// them to a file
// Does not work with MPI!!

#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>
#include<vector>
#include<fstream>
#include<sstream>

int main(int argc, char** argv)
{
#ifdef BEM2DMPI
  exit(0);
#endif

  int ppw=200;     // Point per wavelength
  std::string file="/Users/tbetcke/svn/numerical_coercivity/data/squaredensity";

  int computenorm=0; // Set to 1 to compute norm and condition number
  
 
  std::vector<bem2d::freqtype> freqs;
  //freqs.push_back(10);
  freqs.push_back(1);
  //freqs.push_back(100);
  //freqs.push_back(200);

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
	bem2d::Polygon poly(square,ppw,k,10,0.15);
        bem2d::pGeometry pgeom=poly.GetGeometry();

        bem2d::PolBasis::AddBasis(2,pgeom); // Add constant basis functions


	// Discretize the single and double layer potential

	bem2d::SingleLayer sl(k);
	bem2d::ConjDoubleLayer cdl(k);

	bem2d::QuadOption quadopts;

	quadopts.L=3;
        quadopts.N=10;
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
	  std::cout << "Compute Eigenvalues and vectors" << std::endl;
	}
#else
	std::cout << "Compute Eigenvalues and vectors" << std::endl;
#endif

	double norm; double cond;
	
	bem2d::pcvector eigvals;
	bem2d::pMatrix evectors;
	bem2d::Eigenvalues(combined1,eigvals,evectors);

	// Extract the first nvec eigenvectors
	int nvec=50;

	bem2d::Matrix ev2(evectors->dim[0],nvec);
	for (int i=0;i<evectors->dim[0]*nvec;i++)
	  (*ev2.data)[i]=(*evectors->data)[i];

	// Write out eigenvalues

#ifdef BEM2DMPI
	if (b->IsRoot()){
#endif

	  std::ostringstream os;
	  os << file << "_eig_" << k; 
	  std::string s=os.str();
	  std::ofstream o1(s.c_str());
	  for (int i=0;i<eigvals->size();i++) o1 << std::real((*eigvals)[i])
			     << " " << std::imag((*eigvals)[i]) << std::endl;
	  o1.close();

#ifdef BEM2DMPI
	}
#endif

	// Write out eigenfunctions evaluated on boundary

	std::ostringstream os2;
	os2 << file << "_density_" << k;
	std::string s2=os2.str();
	std::ofstream o2(s.c_str());
	bem2d::WriteDensity(s2,ev2,pgeom,100);	

	}
        finish=clock();
        time=(double(finish)-double(start))/CLOCKS_PER_SEC/60;

#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

