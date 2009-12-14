#include <tr1/array>
#include <cstdlib>
#include "boost/shared_ptr.hpp"
#include "bem2d_quadrature.h"

namespace bem2d
{


void Gauss(dvector& x, dvector& w, int N) throw (LapackError)
{
  // Points xa and weights w for N point Gaussian quadrature
  // in the interval [0,1]


  char c = 'V';
  int info;


  x.clear();
  w.clear(); // Make sure that x and w are empty
  x.resize(N);
  w.resize(N);

  if (N == 1)
    {
      x[0]=0.5;
      w[0]=1.0;
      return;
    }

  dvector alpha(N,0);
  dvector beta(N-1,0);
  dvector v(N*N,0);
  dvector work(2*N-2,0);


  /* Fill the Array */

  for (std::size_t i = 0; i < N - 1; i++)
    {
      beta[i] = .5 / sqrt(1. - 1. / (4. * (i + 1)*(i + 1)));
    }

  /* Call Lapack */

  dstev_(&c, &N, &alpha[0], &beta[0], &v[0], &N, &work[0], &info);

  /* Return results */

  for (std::size_t i = 0; i < N; i++)
    {
      x[i] = (1.0 + alpha[i]) / 2.0;
      w[i] = v[i * N] * v[i * N];
    }

  if (info) throw LapackError();

}

void MapPoints(dvector& x, dvector& w, double a, double b)
{

  for (std::size_t i = 0; i < x.size(); i++)
    {
      x[i] = a + (b - a) * x[i];
      w[i] = (b - a) * w[i];
    }

}

void MapPoints2d(dvector& x, dvector& y, dvector& w, double a, double b, double c, double d)
{

  for (std::size_t i=0; i< x.size(); i++)
    {
      x[i]=a+(b-a)*x[i];
      y[i]=c+(d-c)*y[i];
      w[i]=(b-a)*(d-c)*w[i];
    }
}


	pMatrix EvalIdent(const Geometry& geom, QuadOption opts){
		boost::shared_ptr<Geometry::flat_basis_map> bfuns=geom.FlatMap();
		std::size_t N=bfuns->size();
		Gauss1D g1d(opts.N);
		const dvector x=g1d.x();
		const dvector w=g1d.w();
		
#ifdef BEM2DMPI
		BlacsSystem* b=BlacsSystem::Instance();
		int nrow=b->MSize(N);
		int ncol=b->NSize(N);
#else
		int nrow=N;
		int ncol=N;
#endif
		pMatrix pId(new Matrix(N));
		
		// Now iterate through all the elements
#pragma omp parallel for		
		for (int i=0;i<nrow;i++){
#ifdef BEM2DMPI
			int gi=b->L2gr(N);
#else
			int gi=i;
#endif
			pElement pi=(*bfuns)[gi].first;
			pBasis fi=(*bfuns)[gi].second;
			for (int j=0;j<ncol;j++){
#ifdef BEM2DMPI
				int gj=b->L2gc(N);
#else
				int gj=j;
#endif
				pElement pj=(*bfuns)[gj].first;
				if (pi->index()==pj->index()){
					pBasis fj=(*bfuns)[gj].second;
					complex value=0;
					for (int t=0;t<x.size();t++){
						double s=length(pj->Deriv(x[t]));
						value+=s*std::conj((*fi)(x[t]))*(*fj)(x[t])*w[t];
					}
					(*pId->data)[N*j+i]=value;
				}
			}
		}
		return pId;
	}
	





}
