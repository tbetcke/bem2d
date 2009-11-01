#ifndef _SOUNDSOFTSCATTERING_H_
#define _SOUNDSOFTSCATTERING_H_

#include <vector>
#include "bem2d_defs.h"
#include "bem2d_geometry.h"
#include "bem2d_fun.h"
#include "bem2d_outputhandler.h"
#include "bem2d_quadrature.h"
#include "bem2d_mathroutines.h"
#include "bem2d_cblas.h"

namespace bem2d
{

template<typename T1, typename T2>
class SoundSoftScattering
{
public:
  SoundSoftScattering(pGeometry geom, freqtype kvalue, T1 inc, T2 rhs);
  virtual ~SoundSoftScattering();


  void SetQuadOption(int L, int N, double sigma);
  virtual void Discretize();

  void Solve();

  pcvector EvalIncident(const std::vector<Point>& points);
  pcvector EvalSol(const std::vector<Point>& points);

  void SetOutput(pOutputHandler output);
  void WriteIncident();
  void WriteScattered();
  void WriteFull();
  void WriteAll();

private:
  pGeometry pgeom_;
  freqtype k_;

  T1 incident_;
  T2 frhs_;

  pOutputHandler pout_;

  QuadOption quadopts_;
  Matrix A_;
  pcvector prhs_;
  pcvector psol_;


};

template<typename T1, typename T2>
SoundSoftScattering<T1,T2>::SoundSoftScattering(pGeometry pgeom, freqtype k, T1 inc, T2 rhs): pgeom_(pgeom),
    k_(k), A_(pgeom->size()), frhs_(rhs), incident_(inc)
{
  quadopts_.L=3;
  quadopts_.N=5;
  quadopts_.sigma=0.15;
}

template<typename T1, typename T2>
SoundSoftScattering<T1,T2>::~SoundSoftScattering() {}


template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::SetQuadOption(int L, int N, double sigma)
{
  quadopts_.L=L;
  quadopts_.N=N;
  quadopts_.sigma=sigma;
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::Discretize()
{


  bem2d::CombinedSingleConjDouble scdl(k_,frhs_.eta());
  A_.data=DiscreteKernel(*pgeom_,quadopts_,scdl);
  pcvector pident=DotProdBasFuns(*pgeom_,quadopts_,IdFun());
  cblas_zdscal(A_.dim*A_.dim,2.0,&(*A_.data)[0],1);

  for (int i=0; i<A_.dim; i++) (*A_.data)[A_.dim*i+i]+=(*pident)[i];

  prhs_=bem2d::DotProdBasFuns(*pgeom_,quadopts_,frhs_);
  cblas_zdscal(A_.dim,2.0,&(*prhs_)[0],1);

  /*
  std::cout << "A Matrix" << std::endl;
  for (int i=0;i<A.dim;i++) std::cout << (*A.data)[i] << std::endl;

  std::cout << "rhs";
  for (int i=0;i<A.dim;i++) std::cout << (*prhs)[i] << std::endl;
  std::cout << std::endl;
  */
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::Solve()
{
  psol_=pcvector(new cvector(*prhs_));
  SolveSystem(A_.data,psol_);

  /*
  std::cout << "Sol" << std::endl;
  for (int i=0;i<A.dim;i++) std::cout << (*psol)[i] << std::endl;
   */
}

template<typename T1,typename T2>
pcvector SoundSoftScattering<T1,T2>::EvalIncident(const std::vector<Point>& points)
{
  pcvector pz(new cvector(points.size()));
#pragma omp for
  for (int i=0; i<points.size(); i++) (*pz)[i]=incident_(points[i]);
  return pz;
}

template<typename T1, typename T2>
pcvector SoundSoftScattering<T1,T2>::EvalSol(const std::vector<Point>& points)
{
  boost::shared_ptr<bem2d::Geometry::flat_basis_map> bfuns=pgeom_->FlatMap();
  bem2d::SingleLayer sl(k_);
  pcvector kerneleval=pcvector(new cvector(points.size()));

  bem2d::Gauss1D g1d(quadopts_.N);

  /*
  for (int i=0;i<A.dim;i++) {
  	evalkernel((*bfuns)[i], points, *result, (*psol)[i], dl,g1d);
  }
   */
  EvalKernel(*bfuns,points,*kerneleval,*psol_,sl,g1d);

  /*
  std::cout << "Single Layer Kernel";
  for (int i=0;i<points.size();i++) std::cout << (*kerneleval)[i] << std::endl;
  */


  pcvector result=EvalIncident(points);
  complex alpha=-1.0;
  cblas_zaxpy(points.size(),&alpha, &(*kerneleval)[0], 1, &(*result)[0], 1);
  return result;


}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::SetOutput(pOutputHandler output)
{
  pout_=output;
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::WriteIncident()
{

  pcvector pinc=EvalIncident(pout_->mesh());
  pout_->WriteIncident(*pinc);
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::WriteScattered()
{
  pcvector psolution=EvalSol(pout_->mesh());
  pcvector pincident=EvalIncident(pout_->mesh());
  pcvector pscattered=pcvector(new cvector(*psolution));
  complex alpha=-1.0;
  cblas_zaxpy(psolution->size(),&alpha, &(*pincident)[0], 1, &(*pscattered)[0], 1);
  pout_->WriteScattered(*pscattered);

}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::WriteFull()
{
  pcvector psolution=EvalSol(pout_->mesh());
  pout_->WriteFull(*psolution);
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::WriteAll()
{
  pcvector pincident=EvalIncident(pout_->mesh());
  pcvector pfull=EvalSol(pout_->mesh());

  pcvector pscattered=pcvector(new cvector(*pfull));
  complex alpha=-1.0;
  cblas_zaxpy(pincident->size(),&alpha, &(*pincident)[0], 1, &(*pscattered)[0], 1);
  pout_->WriteAll(*pincident,*pscattered,*pfull);
}




}


#endif // _SOUNDSOFTSCATTERING_H_
