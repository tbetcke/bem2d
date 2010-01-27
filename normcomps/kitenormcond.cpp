#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>
#include<vector>
#include<fstream>

int main(int argc, char** argv)
{

  int ppw=10;     // Point per wavelength
  std::string filename="kitenormcond10.txt";
  
 
  std::vector<bem2d::freqtype> freqs;
  freqs.push_back(5);
  freqs.push_back(10);
  freqs.push_back(20);
  freqs.push_back(40);
  freqs.push_back(80);
  freqs.push_back(160);
  freqs.push_back(320);
  freqs.push_back(640);
  


  std::vector<double> norm_sl(freqs.size());
  std::vector<double> norm_dl(freqs.size());
  std::vector<double> norm_combined1(freqs.size());
  std::vector<double> norm_combined2(freqs.size());
  std::vector<double> cond_sl(freqs.size());
  std::vector<double> cond_dl(freqs.size());
  std::vector<double> cond_combined1(freqs.size());
  std::vector<double> cond_combined2(freqs.size());



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
        bem2d::pCurve kobj(new bem2d::Kite);
	int n=(int)(kobj->Length()*k*ppw/2.0/bem2d::PI);
        bem2d::AnalyticCurve kite(n,kobj);
        bem2d::pGeometry pgeom=kite.GetGeometry();

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
	bem2d::Matrix combined1=Id+2.0*dcdl-bem2d::complex(0,2.0)*eta1*dsl;
	bem2d::Matrix combined2=Id+2.0*dcdl-bem2d::complex(0,2.0)*eta2*dsl;
	
	dsl=2.0*bem2d::ChangeBasis(dsl,Id);
	ddl=2.0*bem2d::ChangeBasis(ddl,Id);
	dcdl=2.0*bem2d::ChangeBasis(dcdl,Id);
	combined1=bem2d::ChangeBasis(combined1,Id);
	combined2=bem2d::ChangeBasis(combined2,Id);

#ifdef BEM2DMPI
	if (b->IsRoot()){
	  std::cout << "Compute norms and condition numbers" << std::endl;
	}
#else
	std::cout << "Compute norms and condition numbers" << std::endl;
#endif
	

	bem2d::L2NormCond(dsl,norm_sl[j],cond_sl[j]);
	bem2d::L2NormCond(ddl,norm_dl[j],cond_dl[j]);
	bem2d::L2NormCond(combined1,norm_combined1[j],cond_combined1[j]);
	bem2d::L2NormCond(combined2,norm_combined2[j],cond_combined2[j]);
	
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

 	  out << "Double Layer" << std::endl;
	  

	  for (int j=0;j<freqs.size();j++){
	    out << "k=" << freqs[j] << " Norm: " << norm_dl[j] << " Norm of Inverse: " << cond_dl[j]/norm_dl[j] << " Condition Nr.: " << cond_dl[j] <<  std::endl;
	  }

 	  out << "Combined Layer eta=k" << std::endl;

	  for (int j=0;j<freqs.size();j++){
	    out << "k=" << freqs[j] << " Norm: " << norm_combined1[j] << " Norm of Inverse: " << cond_combined1[j]/norm_combined1[j] << " Condition Nr.: " << cond_combined1[j] <<  std::endl;
	  }

 	  out << "Combined Layer eta=k^(2/3)" << std::endl;

	  for (int j=0;j<freqs.size();j++){
	    out << "k=" << freqs[j] << " Norm: " << norm_combined2[j] << " Norm of Inverse: " << cond_combined2[j]/norm_combined2[j] << " Condition Nr.: " << cond_combined2[j] <<  std::endl;
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

