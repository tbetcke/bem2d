// Compute resonances of the trapping domain 
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
  std::string file="trapresonance";

  double krmin=4;
  double krmax=8;
  double kimin=-2;
  double kimax=0;
  int krsteps=10;
  int kisteps=10;
  std::vector<bem2d::pMatrix> matrices(krsteps*kisteps);
  std::vector<double> norminv(krsteps*kisteps);
 
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


	// Do the resonance computation

	std::cout << "Discretize Matrices" << std::endl;

	for (int i=0;i<krsteps;i++)
	  for (int j=0;j<kisteps;j++)
	    {
	      double kr=krmin+i*(krmax-krmin)/krsteps;
	      double ki=kimin+j*(kimax-kimin)/kisteps;

	      bem2d::freqtype k={kr,ki};
	      bem2d::freqtype kmax={krmax,kimax};
	      double eta1=std::abs(bem2d::complex(kr,ki)); // Coupling between conj. double and single layer pot.
	      bem2d::Polygon poly(trapping,ppw,kmax,10,0.15);
	      bem2d::pGeometry pgeom=poly.GetGeometry();
	      bem2d::PolBasis::AddBasis(2,pgeom);


	// Discretize the single and double layer potential

	      bem2d::SingleLayer sl(k);
	      bem2d::ConjDoubleLayer cdl(k);

	      bem2d::QuadOption quadopts;

	      quadopts.L=3;
	      quadopts.N=5;
	      quadopts.sigma=0.15;
	      bem2d::Matrix dsl=*(DiscreteKernel(*pgeom,quadopts,sl));
	      bem2d::Matrix dcdl=*(DiscreteKernel(*pgeom,quadopts,cdl));
	      bem2d::Matrix Id=*(EvalIdent(*pgeom, quadopts));
	      bem2d::Matrix combined1=Id+2.0*dcdl-bem2d::complex(0,2.0)*eta1*dsl;
	      matrices[i*kisteps+j]=bem2d::pMatrix(new bem2d::Matrix(bem2d::ChangeBasis(combined1,Id)));
	    }

	std::cout << "Compute inv-norms" << std::endl;

	// Compute the smallest singular values
#pragma omp parallel for
	for (int i=0;i<kisteps*krsteps;i++)
	  {
	    bem2d::pdvector singvals;
	    bem2d::SingularValues(*matrices[i],singvals);
	    int n=singvals->size();
	    norminv[i]=1/(*singvals)[n-1];
	  }

	std::cout << "Write out results" << std::endl;

	  std::ostringstream krname;
	  std::ostringstream kiname;
	  std::ostringstream valsname;
	  krname << file << "_kr";
	  kiname << file << "_ki";
	  valsname << file << "_vals";

	  std::string krs=krname.str();
	  std::string kis=kiname.str();
	  std::string valss=valsname.str();

	  std::ofstream okr(krs.c_str());
	  std::ofstream oki(kis.c_str());
	  std::ofstream ovals(valss.c_str());

	  for (int i=0;i<kisteps;i++)
	    {
	      for (int j=0;j<krsteps;j++)
		{
		  double kr=krmin+i*(krmax-krmin)/krsteps;
		  double ki=kimin+j*(kimax-kimin)/kisteps;
	      
		  okr << kr << " ";
		  oki << ki << " ";
		  ovals << norminv[i*kisteps+j] << " ";
		}
	      okr << std::endl;
	      oki << std::endl;
	      ovals << std::endl;
	    }
	  okr.close();
	  oki.close();
	  ovals.close();

	  return 0;
}

