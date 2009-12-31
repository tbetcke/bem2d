#ifndef _QUADRATURE_H
#define	_QUADRATURE_H

#include <map>
#include <complex>
#include "boost/shared_ptr.hpp"
#include "bem2d_defs.h"
#include "bem2d_exceptions.h"
#include "bem2d_kernel.h"
#include "bem2d_element.h"
#include "bem2d_basis.h"
#include "bem2d_quadpoints.h"
#include "bem2d_point.h"
#include "bem2d_geometry.h"
#include "bem2d_fun.h"
#include "bem2d_mathroutines.h"

#ifdef BEM2DMPI
#include "bem2d_mpi.h"
#endif

namespace bem2d
{


struct QuadOption {
        int L;
        int N;
        double sigma;
};


void Gauss(dvector& x, dvector& w, int N) throw (LapackError);
void MapPoints(dvector& x, dvector& w, double a, double b);
void MapPoints2d(dvector& x, dvector& y, dvector& w, double a, double b, double c, double d);

template<typename kernel,typename T>
complex Integrate(std::pair<pElement,pBasis>& testfun, std::pair<pElement,pBasis>& basfun, kernel g, T& g2d)
{

        pElement elem1=testfun.first;
        pElement elem2=basfun.first;
        pBasis fun1=testfun.second;
        pBasis fun2=basfun.second;
        const dvector x=g2d.x();
        const dvector y=g2d.y();
        const dvector w=g2d.w();
        std::size_t t1=elem1->index();
        std::size_t t2=elem2->index();

        complex result=0;
        for (int i=0; i<g2d.Size(); i++) {
                g.SetNormal(elem1->Normal(x[i]),elem2->Normal(y[i]));
                complex f1=std::conj((*fun1)(x[i]));
                complex f2=(*fun2)(y[i]);
                Point xp=elem1->Map(x[i]);
                Point yp=elem2->Map(y[i]);
                complex gvalue=g(xp,yp);
                double s1=length(elem1->Deriv(x[i]));
                double s2=length(elem2->Deriv(y[i]));
                result+= gvalue*f1*f2*s1*s2*w[i];
        }
        return result;

}

template<typename kernel>
complex IntegrateIdent(std::pair<pElement,pBasis> testfun, std::pair<pElement,pBasis> basfun, kernel g, AdaptedGauss1& g2d)
{
        pElement elem1=testfun.first;
        pElement elem2=basfun.first;
        pBasis fun1=testfun.second;
        pBasis fun2=basfun.second;
        const dvector x=g2d.x();
        const dvector y=g2d.y();
        const dvector w=g2d.w();

        complex result=0;
        for (int i=0; i<g2d.Size(); i++) {
                g.SetNormal(elem1->Normal(x[i]),elem2->Normal(y[i]));

                complex f1=std::conj((*fun1)(x[i]));
                complex f2=(*fun2)((1-y[i])*x[i]);
                Point xp=elem1->Map(x[i]);
                Point yp=elem2->Map((1-y[i])*x[i]);
                complex gvalue=g(xp,yp);
                double s1=length(elem1->Deriv(x[i]));
                double s2=length(elem2->Deriv((1-y[i])*x[i]));
                result+= x[i]*gvalue*f1*f2*s1*s2*w[i];


                f1=std::conj((*fun1)(1-x[i]));
                f2=(*fun2)(1-x[i]+y[i]*x[i]);
                Point xpp=elem1->Map(1-x[i]);
                Point ypp=elem2->Map(1-x[i]+y[i]*x[i]);
                gvalue=g(xpp,ypp);
                s1=length(elem1->Deriv(1-x[i]));
                s2=length(elem2->Deriv(1-x[i]+y[i]*x[i]));
                result+= x[i]*gvalue*f1*f2*s1*s2*w[i];


        }
        return result;


}

template<typename Bem2dfun>
complex Integrate(std::pair<pElement,pBasis> testfun, Bem2dfun fun, Gauss1D& g1d)
{
        pElement elem=testfun.first;
        pBasis tf=testfun.second;
        const dvector x=g1d.x();
        const dvector w=g1d.w();

        complex result=0;
        for (int i=0; i<g1d.Size(); i++) {
                fun.SetNormal(elem->Normal(x[i]));
                complex tfval=std::conj((*tf)(x[i]));
                Point xp=elem->Map(x[i]);
                double s=length(elem->Deriv(x[i]));
                result+=fun(xp)*s*tfval*w[i];
        }
        return result;
}

template<typename Bem2dfun>
pMatrix DotProdBasFuns(const Geometry& geom, QuadOption ops, Bem2dfun fun)
{

        int N=geom.size();
        pMatrix prhs(new Matrix(N,1));
        boost::shared_ptr<Geometry::flat_basis_map> bfuns=geom.FlatMap();
        Gauss1D g1d(ops.N);

#ifdef BEM2DMPI
        BlacsSystem* b=BlacsSystem::Instance();
        int nrow=b->MSize(N);
#else
        int nrow=N;
#endif

#pragma omp parallel for
        for (int i=0; i<nrow; i++) {
#ifdef BEM2DMPI
                int gi=b->L2gr(i);
#else
                int gi=i;
#endif
                (*prhs->data)[i]=Integrate((*bfuns)[gi],fun, g1d);
        }
        return prhs;
}

template<typename kernel>
void EvalKernel(std::pair<pElement,pBasis>& basfun, const std::vector<Point>& points,
                cvector& vals, complex alpha, kernel g, Gauss1D& g1d) throw (SizeError)
{
        pElement elem=basfun.first;
        pBasis bf=basfun.second;
        const dvector x=g1d.x();
        const dvector w=g1d.w();
        if (points.size()!=vals.size()) throw SizeError();
        int N=points.size();


        for (int j=0; j<N; j++) {
                for (int i=0; i<g1d.Size(); i++) {
                        g.SetNormal(elem->Normal(x[i]),elem->Normal(x[i]));
                        complex basval=(*bf)(x[i]);
                        Point xp=elem->Map(x[i]);
                        double s=length(elem->Deriv(x[i]));
                        vals[j]+=alpha*g(points[j],xp)*s*basval*w[i];
                }
        }

}

pMatrix EvalIdent(const Geometry& geom, QuadOption opts);

template<typename kernel>
void EvalKernel(const Geometry::flat_basis_map& basfuns, const std::vector<Point>& points,
                cvector& vals, Matrix& alpha, kernel g, const Gauss1D& g1d) throw (SizeError)
{

#ifdef BEM2DMPI
        BlacsSystem* b=BlacsSystem::Instance();
        int nrow=b->MSize(basfuns.size());
#else
        int nrow=basfuns.size();
#endif

        const dvector x=g1d.x();
        const dvector w=g1d.w();
        if (points.size()!=vals.size()) throw SizeError();
        int N=points.size();

#pragma omp parallel for
        for (int j=0; j<N; j++) {
                for (int itbas=0; itbas<nrow; itbas++) {
#ifdef BEM2DMPI
                        int global_itbas=b->L2gr(itbas);
#else
                        int global_itbas=itbas;
#endif
                        pElement elem=basfuns[global_itbas].first;
                        pBasis bf=basfuns[global_itbas].second;
                        kernel bg(g);
                        for (int i=0; i<g1d.Size(); i++) {
                                bg.SetNormal(elem->Normal(x[i]),elem->Normal(x[i]));
                                complex basval=(*bf)(x[i]);
                                Point xp=elem->Map(x[i]);
                                double s=length(elem->Deriv(x[i]));
                                vals[j]+=(*alpha.data)[itbas]*bg(points[j],xp)*s*basval*w[i];
                        }
                }
        }
#ifdef BEM2DMPI
        // Now distribute all information in vals to all nodes
        char ALL[]="ALL";
        char EMPTY[]=" ";
        Czgsum2d(b->get_ictxt(),ALL,EMPTY,points.size(),1,&vals[0],points.size(),-1,-1);
#endif

}


template<typename kernel>
pMatrix DiscreteKernel(const Geometry& geom, QuadOption opts, kernel g)
{

        boost::shared_ptr<Geometry::flat_basis_map> bfuns=geom.FlatMap();
        std::size_t N=bfuns->size();
        pMatrix dkernel(new Matrix(N));

        Gauss2D g0(opts.N);
        AdaptedGauss1 g1(opts.N,opts.L,opts.sigma);
        AdaptedGauss2 g2(opts.N,opts.L,opts.sigma);
        AdaptedGauss3 g3(opts.N,opts.L,opts.sigma);

#ifdef BEM2DMPI
        BlacsSystem* b=BlacsSystem::Instance();
        int nrow=b->MSize(N);
        int ncol=b->NSize(N);
#else
        int nrow=N;
        int ncol=N;
#endif

#pragma omp parallel for
        for (int j=0; j<ncol; j++) {
                for (int i=0; i<nrow; i++) {
#ifdef BEM2DMPI
                        int gi=b->L2gr(i); // Global i index
                        int gj=b->L2gc(j); // Global j index
#else
                        int gi=i;
                        int gj=j;
#endif

                        std::pair<pElement,pBasis> pi((*bfuns)[gi]);
                        std::pair<pElement,pBasis> pj((*bfuns)[gj]);
                        std::size_t indi=pi.first->index();
                        std::size_t indj=pj.first->index();
                        //std::cout << indi << " " << indj << " " << i << " " << j << std::endl;


                        if (indi==indj) {
                                // Identical elements
                                //std::cout << "Integrate Ident" << std::endl;
                                (*dkernel->data)[i+nrow*j]=IntegrateIdent(pi,pj,g,g1);
                        } else if (pi.first->next()==indj) {
                                // (1,0) situation
                                //std::cout << "Integrate (1,0)" << std::endl;
                                (*dkernel->data)[i+nrow*j]=Integrate(pi,pj,g,g3);
                        } else if (pi.first->prev()==indj) {
                                // (0,1) situation
                                //std::cout << "Integrate (0,1)" << std::endl;
                                (*dkernel->data)[i+nrow*j]=Integrate(pi,pj,g,g2);
                        } else {
                                // Elements are not neighbors
                                //std::cout << "Integrate Remote" << std::endl;
                                (*dkernel->data)[i+nrow*j]=Integrate(pi,pj,g,g0);
                        }
                        //std::cout << (*pmatrix)[i+N*j] << std::endl;

                }
        }
        return dkernel;

}




}
#endif	/* _QUADRATURE_H */

