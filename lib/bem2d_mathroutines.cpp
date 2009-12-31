#include<algorithm>
#include<cmath>
#include "bem2d_mathroutines.h"
#include "gsl/gsl_sf_bessel.h"
#include "bem2d_cblas.h"



#ifdef BEM2DMPI
#include "bem2d_mpi.h"
#endif


namespace bem2d
{


complex BesselH0(double x)
{
        double fr, fi;

        // Evaluate Bessel function

        fr = gsl_sf_bessel_J0(x);
        fi = gsl_sf_bessel_Y0(x);
        return complex(fr, fi);

}

complex BesselH1(double x)
{
        double fr, fi;

        // Evaluate Bessel function

        fr = gsl_sf_bessel_J1(x);
        fi = gsl_sf_bessel_Y1(x);
        return complex(fr, fi);

}



Matrix::Matrix(int n)
{

        dim[0]=n;
        dim[1]=n;
#ifdef BEM2DMPI
        BlacsSystem* b=BlacsSystem::Instance();
        int mb=b->get_mb();
        int nb=b->get_nb();
        int irsrc=0;
        int icsrc=0;
        int ictxt=b->get_ictxt();
        int lld=std::max(b->MSize(n),1);
        int info;
        descinit_(desc,&n,&n,&mb,&nb,&irsrc,&icsrc,&ictxt,&lld,&info);
        if (info) throw ScaLapackError();
        msize=std::max(b->MSize(n),1);
        nsize=std::max(b->NSize(n),1);
#else
        msize=n;
        nsize=n;
#endif
        data=pcvector(new cvector(msize*nsize));
}


Matrix::Matrix(int m, int n)
{
        dim[0]=m;
        dim[1]=n;
#ifdef BEM2DMPI
        BlacsSystem* b=BlacsSystem::Instance();
        int mb=b->get_mb();
        int nb=b->get_nb();
        int irsrc=0;
        int icsrc=0;
        int ictxt=b->get_ictxt();
        int lld=std::max(b->MSize(m),1);
        int info;
        descinit_(desc,&m,&n,&mb,&nb,&irsrc,&icsrc,&ictxt,&lld,&info);
        if (info) throw ScaLapackError();
        msize=std::max(b->MSize(m),1);
        nsize=std::max(b->NSize(n),1);
#else
        msize=m;
        nsize=n;
#endif
        data=pcvector(new cvector(msize*nsize));
}


Matrix::Matrix(const Matrix& m)
{
        dim[0]=m.dim[0];
        dim[1]=m.dim[1];
        msize=m.msize;
        nsize=m.nsize;
        data=pcvector(new cvector(*m.data));
#ifdef BEM2DMPI
        for (int i=0; i<9; i++) desc[i]=m.desc[i];
#endif

}


Matrix operator+(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch)
{
        if (lhs.dim[0]!=rhs.dim[0]) throw ArrayMismatch();
        if (lhs.dim[1]!=rhs.dim[1]) throw ArrayMismatch();
        Matrix result(rhs);
        complex alpha=1.0;
        cblas_zaxpy(lhs.msize*lhs.nsize,&alpha,&(*lhs.data)[0],1,&(*result.data)[0],1);
        return result;

}

Matrix operator-(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch)
{
        if (lhs.dim[0]!=rhs.dim[0]) throw ArrayMismatch();
        if (lhs.dim[1]!=rhs.dim[1]) throw ArrayMismatch();
        Matrix result(lhs);
        complex alpha=-1.0;
        cblas_zaxpy(lhs.msize*lhs.nsize,&alpha,&(*rhs.data)[0],1,&(*result.data)[0],1);
        return result;

}


Matrix operator*(const Matrix& lhs, const complex& alpha)
{
        Matrix result(lhs);
        cblas_zscal(lhs.msize*lhs.nsize,&alpha,&(*result.data)[0],1);
        return result;
}

Matrix operator*(const complex& alpha, const Matrix& rhs)
{
        return rhs*alpha;
}

Matrix operator*(const Matrix& lhs, const double& alpha)
{
        Matrix result(lhs);
        cblas_zdscal(lhs.msize*lhs.nsize,alpha,&(*result.data)[0],1);
        return result;
}

Matrix operator*(const double& alpha, const Matrix& rhs)
{
        return rhs*alpha;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch)
{
        if (lhs.dim[1]!=rhs.dim[0]) throw ArrayMismatch();
        Matrix result(lhs.dim[0],rhs.dim[1]);
        complex alpha=1.0;
        complex beta=0.0;
#ifdef BEM2DMPI
        char trans='N';
        int ione=1;
        pzgemm_(&trans,&trans,&lhs.dim[0],&rhs.dim[1],&lhs.dim[1],&alpha,&(*lhs.data)[0],&ione,&ione,lhs.desc,&(*rhs.data)[0],&ione,&ione,rhs.desc,&beta,&(*result.data)[0],&ione,&ione,result.desc);
#else
        cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lhs.dim[0],rhs.dim[1],lhs.dim[1],&alpha,&(*lhs.data)[0],lhs.dim[0],&(*rhs.data)[0],
                    rhs.dim[0],&beta,&(*result.data)[0],rhs.dim[0]);
#endif
        return result;

}

#ifdef BEM2DMPI
void L2NormCond(const Matrix& stiff, const Matrix& mass, double& norm, double& cond) throw (ScaLapackError)
{
#else
void L2NormCond(const Matrix& stiff, const Matrix& mass, double& norm, double& cond) throw (LapackError)
{
#endif


        Matrix tstiff(stiff);
        Matrix tmass(mass);

        // Compute Cholesky decomposition of mass matrix
#ifdef BEM2DMPI
        // Compute Cholesky decomposition of mass matrix
        char uplo='U';
        int ia=1;
        int ja=1;
        int info;
        pzpotrf_(&uplo, &tmass.dim[0], &(*tmass.data)[0],&ia, &ja,tmass.desc,&info);
        if (info) throw ScaLapackError();

        // Factor Cholesky decomposition into original matrix

        char side='R';
        char transa='N';
        char diag='N';
        complex alpha=1.0;

        pztrsm_(&side,&uplo,&transa,&diag,&tstiff.dim[0],&tstiff.dim[1],&alpha,&(*tmass.data)[0],&ia,&ja,tmass.desc,&(*tstiff.data)[0],&ia,&ja,tstiff.desc);
        side='L';
        transa='C';
        pztrsm_(&side,&uplo,&transa,&diag,&tstiff.dim[0],&tstiff.dim[1],&alpha,&(*tmass.data)[0],&ia,&ja,tmass.desc,&(*tstiff.data)[0],&ia,&ja,tstiff.desc);

        // tstiff now includes the mass matrix

        // Compute singular values

        char jobu='N';
        char jobvt='N';
        int n=tstiff.dim[1];
        double s[n];
        complex wsize;
        int lwork=-1;
        double rwork[1+4*n];



        // Compute Workspace Size
        pzgesvd_(&jobu,&jobvt,&tstiff.dim[0],&tstiff.dim[1],&(*tstiff.data)[0],&ia,&ja,tstiff.desc,s,0,&ia,&ja,tstiff.desc,0,&ia,&ja,tstiff.desc,&wsize,&lwork,rwork,&info);
        // Now do the actual call

        lwork=(int) real(wsize);
        complex work[lwork];

        pzgesvd_(&jobu,&jobvt,&tstiff.dim[0],&tstiff.dim[1],&(*tstiff.data)[0],&ia,&ja,tstiff.desc,s,0,&ia,&ja,tstiff.desc,0,&ia,&ja,tstiff.desc,work,&lwork,rwork,&info);

        norm=s[0];
        cond=s[0]/s[n-1];
#else

        // Compute Cholesky of mass matrix

        char uplo='U';
        int info;
        complex alpha=1.0;
        zpotrf_(&uplo,&tmass.dim[0],&(*tmass.data)[0],&tmass.dim[0],&info);

        // Factor Cholesky factor into stiffness matrix

        cblas_ztrsm(CblasColMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,tstiff.dim[0],tstiff.dim[0],&alpha,&(*tmass.data)[0],tmass.dim[0],&(*tstiff.data)[0],tstiff.dim[0]);

        cblas_ztrsm(CblasColMajor,CblasLeft,CblasUpper,CblasConjTrans,CblasNonUnit,tstiff.dim[0],tstiff.dim[0],&alpha,&(*tmass.data)[0],tmass.dim[0],&(*tstiff.data)[0],tstiff.dim[0]);

        // Compute singular values of stiffness matrix

        char jobu='N';
        char jobvt='N';
        int M=tstiff.dim[0];
        dvector s(M);
        int lwork=5*M;
        cvector work(lwork);
        dvector rwork(lwork);

        zgesvd_(&jobu, &jobvt, &M, &M, &(*tstiff.data)[0], &M, &s[0], NULL, &M, NULL, &M, &work[0], &lwork, &rwork[0], &info);
        if (info!=0) throw LapackError();
        norm=s[0];
        cond=s[0]/s[M-1];
#endif
}


#ifdef BEM2DMPI
pMatrix SolveSystem(Matrix& m, Matrix& rhs) throw (ScaLapackError)
#else
pMatrix SolveSystem(Matrix& m, Matrix& rhs) throw (LapackError)
#endif
{
        char trans='N';
        Matrix tmp(m);
        pMatrix res(new Matrix(rhs));
        int info;
#ifdef BEM2DMPI
        BlacsSystem* b=BlacsSystem::Instance();
        int ipiv[m.msize+b->get_mb()];
        int ione=1;
        pzgetrf_(&tmp.dim[0],&tmp.dim[1],&(*tmp.data)[0],&ione,&ione,tmp.desc,ipiv,&info);
        if (info) throw ScaLapackError();
        pzgetrs_(&trans,&tmp.dim[0],&res->dim[1],&(*tmp.data)[0],&ione,&ione,tmp.desc,ipiv,&(*res->data)[0],&ione,&ione,res->desc,&info);
        if (info) throw ScaLapackError();
#else
        int nrhs=1;
        //int ipiv[std::min(tmp.dim[0],tmp.dim[1])];
        int ipiv[tmp.dim[0]];
        zgetrf_(&tmp.dim[0], &tmp.dim[1], &(*tmp.data)[0], &tmp.dim[0], ipiv,&info);
        if (info) throw LapackError();
        zgetrs_(&trans,&tmp.dim[0],&res->dim[1],&(*tmp.data)[0],&tmp.dim[0],ipiv,&(*res->data)[0],&res->dim[0],&info);
        if (info) throw LapackError();
#endif
        return res;
}

void InPolygon(const std::vector<pGeometry> polygons, const std::vector<Point>& testpoints,std::vector<int>& inpoly)
{
        int N=testpoints.size();
        if (N!=inpoly.size()) inpoly.resize(N);
        for (int i=0; i<N; i++) inpoly[i]=0;
        for (int i=0; i<polygons.size(); i++) {
                pGeometry pgeom=polygons[i];
                PointVector polypoints((*pgeom).points());
                dvector px(polypoints.size());
                dvector py(polypoints.size());
                for (int j=0; j<polypoints.size(); j++) {
                        px[j]=polypoints[j].x;
                        py[j]=polypoints[j].y;
                }
#pragma omp parallel for
                for (int j=0; j<N; j++) {
                        inpoly[j]+=pnpoly(polypoints.size(), &px[0], &py[0], testpoints[j].x, testpoints[j].y);
                }
        }
}

void InPolygon(const pGeometry polygon, const std::vector<Point>& testpoints,std::vector<int>& inpoly)
{
        std::vector<pGeometry> polygons;
        polygons.push_back(polygon);
        InPolygon(polygons,testpoints,inpoly);
}



}
