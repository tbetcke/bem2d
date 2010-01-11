#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>
#include<vector>

int main(int argc, char** argv)
{

  int ppw=40;     // Point per wavelength

 
  std::vector<bem2d::freqtype> freqs;
  freqs.push_back(5);
  freqs.push_back(10);
  freqs.push_back(20);

  std::vector<double> norm_sl(freqs.size());
  std::vector<double> norm_dl(freqs.size());
  std::vector<double> norm_combined(freqs.size());
  std::vector<double> cond_sl(freqs.size());
  std::vector<double> cond_dl(freqs.size());
  std::vector<double> cond_combined(freqs.size());



        clock_t start, finish;
        double time;
        start=clock();
 

#ifdef BEM2DMPI
        MPI_Init(&argc, &argv);


        int nprow=2; // Number of rows in process grid
        int npcol=1; // Number of columns in process grid
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

	for (int j=0;j<freqs.size();j++){

	double k=(double)freqs[j];
	double eta=k; // Coupling between conj. double and single layer pot.
        bem2d::pCurve cobj(new bem2d::Circle);
	int n=(int)(cobj->Length()*k*ppw/2/bem2d::PI);
        bem2d::AnalyticCurve circle(n,cobj);
        bem2d::pGeometry pgeom=circle.GetGeometry();

        bem2d::PolBasis::AddBasis(0,pgeom); // Add constant basis functions


	// Discretize the single and double layer potential

	bem2d::SingleLayer sl(k);
	bem2d::ConjDoubleLayer cdl(k);
	bem2d::DoubleLayer dl(k);

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
	bem2d::Matrix ddl=(*DiscreteKernel(*pgeom,quadopts,dl));
	bem2d::Matrix dcdl=*(DiscreteKernel(*pgeom,quadopts,cdl));
	bem2d::Matrix Id=*(EvalIdent(*pgeom, quadopts));
	bem2d::Matrix combined=Id+2.0*dcdl-bem2d::complex(0,2.0)*eta*dsl;
	
	dsl=2.0*bem2d::ChangeBasis(dsl,Id);
	ddl=2.0*bem2d::ChangeBasis(ddl,Id);
	dcdl=2.0*bem2d::ChangeBasis(dcdl,Id);
	combined=bem2d::ChangeBasis(combined,Id);

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Compute norms and condition numbers" << std::endl;
	}
#else
	std::cout << "Compute norms and condition numbers" << std::endl;
#endif
	

	bem2d::L2NormCond(dsl,norm_sl[j],cond_sl[j]);
	bem2d::L2NormCond(ddl,norm_dl[j],cond_dl[j]);
	bem2d::L2NormCond(combined,norm_combined[j],cond_combined[j]);
	
	}
        finish=clock();
        time=(double(finish)-double(start))/CLOCKS_PER_SEC/60;


#ifdef BEM2DMPI
	if (b->IsRoot()){
#endif

 	  std::cout << "Single Layer" << std::endl;

	  for (int j=0;j<freqs.size();j++){
	    std::cout << "k=" << freqs[j] << " Norm: " << norm_sl[j] << " Norm of Inverse: " << cond_sl[j]/norm_sl[j] << " Condition Nr.: " << cond_sl[j] <<  std::endl;
	  }

 	  std::cout << "Double Layer" << std::endl;

	  for (int j=0;j<freqs.size();j++){
	    std::cout << "k=" << freqs[j] << " Norm: " << norm_dl[j] << " Norm of Inverse: " << cond_dl[j]/norm_dl[j] << " Condition Nr.: " << cond_dl[j] <<  std::endl;
	  }

 	  std::cout << "Combined Layer" << std::endl;

	  for (int j=0;j<freqs.size();j++){
	    std::cout << "k=" << freqs[j] << " Norm: " << norm_combined[j] << " Norm of Inverse: " << cond_combined[j]/norm_combined[j] << " Condition Nr.: " << cond_combined[j] <<  std::endl;
	  }

	  std::cout << "Overalll time (minutes): " << time << std::endl;
#ifdef BEM2DMPI	
}
#endif


#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

