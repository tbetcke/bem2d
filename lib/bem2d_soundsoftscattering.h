#ifndef _SOUNDSOFTSCATTERING_H_
#define _SOUNDSOFTSCATTERING_H_

#include <vector>
#include <limits>
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
		inline Matrix GetMatrix() const{
			return Matrix(A_);
		}
		inline Matrix GetIdent() const{
			return Matrix(Id_);
		}
		
		pcvector EvalIncident(const std::vector<Point>& points);
		pcvector EvalSol(const std::vector<Point>& points);
		inline double L2Condition() const{
			return L2Cond(A_);
		}
		
		void SetOutput(pOutputHandler output);
		void WriteIncident();
		void WriteScattered();
		void WriteFull();
		void WriteAll();
		void set_polygons(std::vector<pGeometry>& polygons);
		void set_polygons(pGeometry polygon);
		void set_plotInterior();
		
	private:
		pGeometry pgeom_;
		freqtype k_;
		
		T1 incident_;
		T2 frhs_;
		
		pOutputHandler pout_;
		
		QuadOption quadopts_;
		Matrix A_;
		Matrix Id_;
		pcvector prhs_;
		pcvector psol_;
		std::vector<pGeometry> polygons_;
		bool plotInterior_;
		
		
	};
	
	template<typename T1, typename T2>
	SoundSoftScattering<T1,T2>::SoundSoftScattering(pGeometry pgeom, freqtype k, T1 inc, T2 rhs): pgeom_(pgeom),
    k_(k), A_(pgeom->size()), Id_(pgeom->size()), frhs_(rhs), incident_(inc), plotInterior_(false)
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
		Id_.dim=(pgeom_->size());
		Id_.data=EvalIdent(*pgeom_, quadopts_);
		A_=Id_+2.0*A_;
		
		prhs_=bem2d::DotProdBasFuns(*pgeom_,quadopts_,frhs_);
		
		
		cblas_zdscal(A_.dim,2.0,&(*prhs_)[0],1);
		
	}
	
	template<typename T1, typename T2>
	void SoundSoftScattering<T1,T2>::Solve()
	{
		psol_=pcvector(new cvector(*prhs_));
		Matrix tmp(A_);
		SolveSystem(tmp.data,psol_);
		
	}
	
	template<typename T1,typename T2>
	pcvector SoundSoftScattering<T1,T2>::EvalIncident(const std::vector<Point>& points)
	{
		pcvector pz(new cvector(points.size()));
		std::vector<int> inpoly(points.size());
		if (plotInterior_ & (polygons_.size()>0)){
			InPolygon(polygons_, points,inpoly);
		}		
#pragma omp parallel for
		for (int i=0; i<points.size(); i++){
			if (inpoly[i]==0){
				(*pz)[i]=incident_(points[i]);
			}
			else {
				(*pz)[i]=std::numeric_limits<double>::quiet_NaN();
			}
		}

		return pz;
	}
	
	template<typename T1, typename T2>
	pcvector SoundSoftScattering<T1,T2>::EvalSol(const std::vector<Point>& points)
	{
		boost::shared_ptr<bem2d::Geometry::flat_basis_map> bfuns=pgeom_->FlatMap();
		bem2d::SingleLayer sl(k_);
		pcvector kerneleval=pcvector(new cvector(points.size()));
		
		bem2d::Gauss1D g1d(quadopts_.N);
		
		EvalKernel(*bfuns,points,*kerneleval,*psol_,sl,g1d);
		
		
		pcvector result=EvalIncident(points);
		complex alpha=-1.0;
		cblas_zaxpy(points.size(),&alpha, &(*kerneleval)[0], 1, &(*result)[0], 1);
		
		
		// If 
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
	
	
	template<typename T1, typename T2>
	void SoundSoftScattering<T1,T2>::set_polygons(std::vector<pGeometry>& polygons){
		polygons_=polygons;
	}
	
	template<typename T1, typename T2> 
	void SoundSoftScattering<T1,T2>::set_plotInterior(){
		plotInterior_=true;
	}
	
	template<typename T1, typename T2>
	void SoundSoftScattering<T1,T2>::set_polygons(pGeometry polygon){
		polygons_.push_back(polygon);
	}
	
	
}


#endif // _SOUNDSOFTSCATTERING_H_
