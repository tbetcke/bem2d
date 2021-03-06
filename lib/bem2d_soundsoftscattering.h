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
        virtual void DiscretizeMatrix();
        virtual void DiscretizeRhs();
        void Discretize();

        void Solve();
        inline pMatrix GetMatrix() const {
                return pA_;
        }
        inline pMatrix GetIdent() const {
                return pId_;
        }

        pcvector EvalIncident(const std::vector<Point>& points);
        pcvector EvalSol(const std::vector<Point>& points);
        inline void NormCond(double& norm, double& cond)  const {
	  Matrix k=ChangeBasis(*pA_,*pId_);
	  L2NormCond(k,norm,cond);
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
        pMatrix pA_;
        pMatrix pId_;
        pMatrix prhs_;
        pMatrix psol_;
        std::vector<pGeometry> polygons_;
        bool plotInterior_;


};

template<typename T1, typename T2>
SoundSoftScattering<T1,T2>::SoundSoftScattering(pGeometry pgeom, freqtype k, T1 inc, T2 rhs): pgeom_(pgeom),
                k_(k), frhs_(rhs), incident_(inc), plotInterior_(false)
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
void SoundSoftScattering<T1,T2>::DiscretizeMatrix()
{
        bem2d::CombinedSingleConjDouble scdl(k_,frhs_.eta());

        pA_=DiscreteKernel(*pgeom_,quadopts_,scdl);
        pId_=EvalIdent(*pgeom_, quadopts_);
        (*pA_)=(*pId_)+2.0*(*pA_);
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::DiscretizeRhs()
{
        prhs_=bem2d::DotProdBasFuns(*pgeom_,quadopts_,frhs_);
        (*prhs_)=2.0*(*prhs_);
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::Discretize()
{


        DiscretizeMatrix();
        DiscretizeRhs();
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::Solve()
{
        psol_=SolveSystem(*pA_,*prhs_);

}

template<typename T1,typename T2>
pcvector SoundSoftScattering<T1,T2>::EvalIncident(const std::vector<Point>& points)
{

        pcvector pz(new cvector(points.size()));
        std::vector<int> inpoly(points.size());
        if (plotInterior_ & (polygons_.size()>0)) {
                InPolygon(polygons_, points,inpoly);
        }
#pragma omp parallel for
        for (int i=0; i<points.size(); i++) {
                if (inpoly[i]==0) {
                        (*pz)[i]=incident_(points[i]);
                } else {
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
void SoundSoftScattering<T1,T2>::set_polygons(std::vector<pGeometry>& polygons)
{
        polygons_=polygons;
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::set_plotInterior()
{
        plotInterior_=true;
}

template<typename T1, typename T2>
void SoundSoftScattering<T1,T2>::set_polygons(pGeometry polygon)
{
        polygons_.push_back(polygon);
}


}


#endif // _SOUNDSOFTSCATTERING_H_

