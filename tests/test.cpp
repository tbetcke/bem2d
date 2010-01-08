#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>



int main(int argc, char** argv)
{


#ifdef BEM2DMPI
        MPI_Init(&argc, &argv);


        int nprow=4;
        int npcol=4;
        int mb=200;
        int nb=200;
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
        int n=100;


        bem2d::pCurve cobj(new bem2d::Kite);
        bem2d::AnalyticCurve circle(n,cobj);
        bem2d::pGeometry pgeom=circle.GetGeometry();
	std::cout << circle.Length() << std::endl;
	circle.ParameterizeArc(n);


        //bem2d::DiskShapePiecewiseConst circle(n,1.0);

	/*

        bem2d::freqtype k=5;
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


	std::cout << " Do some tests" << std::endl;
	bem2d::pMatrix pk=soundsoft.GetMatrix();
	bem2d::pMatrix pm=soundsoft.GetIdent();

	*pk=bem2d::ChangeBasis(*pk,*pm);
	bem2d::Matrix H=0.5*(*pk+bem2d::ConjTranspose(*pk));

	bem2d::pdvector pevalues;
	bem2d::pMatrix pevectors;
       
	std::cout <<" Compute HErmitian eigenvalues" << std::endl;

	
	bem2d::HermitianEigenvalues(H,pevalues,pevectors);

	
	std::cout << "Exctract a column and compute inner products" << std::endl;

	bem2d::Matrix z=bem2d::ExtractColumn(*pevectors,0);
	bem2d::Matrix Hz=H*z;
	bem2d::complex lambda=DotProduct(z,H*z);
	bem2d::complex lambda2=DotProduct(z,z);

	std::cout << lambda/lambda2 << " " << (*pevalues)[0] << std::endl;

	
	bem2d::pcvector eigvalues;
	std::cout << " Compute Eigenvalues" << std::endl;
	bem2d::Eigenvalues(*pk,*pm,eigvalues);
	
	
#ifdef BEM2DMPI
	if (b->IsRoot()){
	  for (int i=0;i<pevalues->size();i++) std::cout << (*pevalues)[i] << std::endl;
	}
#else
	  for (int i=0;i<pevalues->size();i++) std::cout << (*pevalues)[i] << std::endl;
#endif
	*/

        //int xpts=100;
        //int ypts=100;
        //bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,2,-2,2,"disk"));
        //soundsoft.SetOutput(pout);



        //bem2d::WriteMatrix("/Users/tbetcke/svn/numerical_coercivity/matlab/diskmatrix10",soundsoft.GetMatrix());
        //bem2d::WriteMatrix("/Users/tbetcke/svn/numerical_coercivity/matlab/iddiskmatrix10",soundsoft.GetIdent());


#ifdef BEM2DMPI
        bem2d::BlacsSystem::Release();
        MPI_Finalize();
#endif

}

