#include<algorithm>
#include<cmath>
#include<complex>
#include<fstream>
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
void L2NormCond(const Matrix& stiff, double& norm, double& cond) throw (ScaLapackError)
{
#else
void L2NormCond(const Matrix& stiff, double& norm, double& cond) throw (LapackError)
{
#endif


        Matrix tstiff(stiff);

#ifdef BEM2DMPI

        // Compute singular values

        char jobu='N';
        char jobvt='N';
        int n=tstiff.dim[1];
        double s[n];
        complex wsize;
        int lwork=-1;
        double rwork[1+4*n];
	int info;
	int ia=1;
	int ja=1;



        // Compute Workspace Size
        pzgesvd_(&jobu,&jobvt,&tstiff.dim[0],&tstiff.dim[1],&(*tstiff.data)[0],&ia,&ja,tstiff.desc,s,0,&ia,&ja,tstiff.desc,0,&ia,&ja,tstiff.desc,&wsize,&lwork,rwork,&info);
        // Now do the actual call

        lwork=(int) real(wsize);
        complex work[lwork];

        pzgesvd_(&jobu,&jobvt,&tstiff.dim[0],&tstiff.dim[1],&(*tstiff.data)[0],&ia,&ja,tstiff.desc,s,0,&ia,&ja,tstiff.desc,0,&ia,&ja,tstiff.desc,work,&lwork,rwork,&info);

        norm=s[0];
        cond=s[0]/s[n-1];
#else


        // Compute singular values of stiffness matrix

        char jobu='N';
        char jobvt='N';
        int M=tstiff.dim[0];
        dvector s(M);
        int lwork=5*M;
        cvector work(lwork);
        dvector rwork(lwork);
	int info;

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

#ifdef BEM2DMPI
 void HermitianEigenvalues(const Matrix& k, pdvector& evalues) throw (ScaLapackError){
#else
   void HermitianEigenvalues(const Matrix& k, pdvector& evalues) throw (LapackError){
#endif
 
     Matrix tk(k);
 
     evalues=pdvector(new dvector(tk.dim[0])); 
 
#ifdef BEM2DMPI
     BlacsSystem* b=BlacsSystem::Instance();
	char jobz='N';
	char uplo='U';
	char range='A';
	int ia=1;
	int ja=1;
	int info;
	complex wsize;
	int lwork=-1;
	double rwsize;
	int lrwork=-1;
	int n=tk.dim[0];
	char cmach='S';
	int ictxt=b->get_ictxt();
	double abstol=2*pdlamch_(&ictxt,&cmach);
	int m,nz;
	double orfac=0;
	int iwsize;
	int liwork=-1;
	int ifail[n];
	int icluster[2*b->get_nprow()*b->get_npcol()];
	double gap[b->get_nprow()*b->get_npcol()];
	int il=1;
	int iu=1;
	int rlwork=-1;
 
	// Workspace query for eigenvalue computation
 
	pzheev_(&jobz,&uplo,&tk.dim[0],&(*tk.data)[0],&ia,&ja,tk.desc,&(*evalues)[0],
		NULL,&ia,&ja,tk.desc,&wsize,&lwork,&rwsize,&rlwork,
		&info);
 
	lwork=(int)std::real(wsize);
	complex work[lwork];
	rlwork=(int)rwsize;
	double rwork[rlwork];
 
	
	pzheev_(&jobz,&uplo,&tk.dim[0],&(*tk.data)[0],&ia,&ja,tk.desc,&(*evalues)[0],NULL,
	        &ia,&ja,tk.desc,work,&lwork,rwork,&rlwork,
		&info);


	if (info) throw ScaLapackError();
	
		
 
#else
     int itype=1;
     char jobz='N';
     char uplo='U';
     complex wsize;
     int lwork=-1;
     double rwork[3*tk.dim[0]-2];
     int info;
     
     // Get work size
 
     zheev_(&jobz,&uplo,&tk.dim[0],&(*tk.data)[0],&tk.dim[0],&(*evalues)[0],
	    &wsize,&lwork,rwork,&info);
 
     lwork=(int)std::real(wsize);
     complex work[lwork];
 
     zheev_(&jobz,&uplo,&tk.dim[0],&(*tk.data)[0],&tk.dim[0],&(*evalues)[0],
	    work,&lwork,rwork,&info);
 
 
#endif
   }



#ifdef BEM2DMPI
   void Eigenvalues(const Matrix& k, pcvector& evalues) throw (ScaLapackError){
#else
   void Eigenvalues(const Matrix& k, pcvector& evalues) throw (LapackError){
#endif

        Matrix tstiff(k);
	evalues=pcvector(new cvector(tstiff.dim[0]));

        // Compute Cholesky decomposition of mass matrix
#ifdef BEM2DMPI

        // Compute Hessenberg form of tstiff

	BlacsSystem* b=BlacsSystem::Instance();

	int ilo=1;
	int ihi=tstiff.dim[0];
	complex tau[b->NSize(tstiff.dim[0]-1)];
	complex wsize;
	int lwork=-1;
	int info;
	int ia=1;
	int ja=1;

	pzgehrd_(&tstiff.dim[0],&ilo,&ihi,&(*tstiff.data)[0],&ia,&ja,tstiff.desc,tau,
		 &wsize,&lwork,&info);

	lwork=(int)std::real(wsize);
	complex work1[lwork];

	pzgehrd_(&tstiff.dim[0],&ilo,&ihi,&(*tstiff.data)[0],&ia,&ja,tstiff.desc,tau,
		 work1,&lwork,&info);

	// Compute eigenvalues of Hessenberg form

	int wantt=0;
	int wantz=0;
	lwork=-1;
	
	pzlahqr_(&wantt,&wantz,&tstiff.dim[0],&ilo,&ihi,&(*tstiff.data)[0],
		 tstiff.desc,&(*evalues)[0],
		 &ilo,&ihi,NULL,tstiff.desc,&wsize,&lwork,NULL,NULL,&info);
        
	lwork=(int)std::real(wsize);
	complex work2[lwork];

	pzlahqr_(&wantt,&wantz,&tstiff.dim[0],&ilo,&ihi,&(*tstiff.data)[0],
		 tstiff.desc,&(*evalues)[0],
		 &ilo,&ihi,NULL,tstiff.desc,work2,&lwork,NULL,NULL,&info);
	if (info) throw ScaLapackError();

#else
	char jobvl='N';
	char jobvr='N';
	int n=tstiff.dim[0];
	int info;
	int lwork=4*n;
	complex work[4*n];
	double rwork[2*n];

	zgeev_(&jobvl,&jobvr,&n,&(*tstiff.data)[0],&n,&(*evalues)[0],
	       NULL,&n,NULL,&n,work,&lwork,rwork,&info);

#endif


   }

#ifndef BEM2DMPI

 void Eigenvalues(const Matrix& k, pcvector& evalues, pMatrix& evectors)
   throw (LapackError){
 

       Matrix tstiff(k);
       evalues=pcvector(new cvector(tstiff.dim[0]));
       evectors=pMatrix(new Matrix(tstiff.dim[0]));

	char jobvl='N';
	char jobvr='V';
	int n=tstiff.dim[0];
	int info;
	int lwork=4*n;
	complex work[4*n];
	double rwork[2*n];

	zgeev_(&jobvl,&jobvr,&n,&(*tstiff.data)[0],&n,&(*evalues)[0],
	       NULL,&n,&(*evectors->data)[0],&n,work,&lwork,rwork,&info);
	if (info) throw LapackError();

 }


#endif


#ifdef BEM2DMPI
   Matrix ChangeBasis(const Matrix& k, const Matrix& m) throw (ScaLapackError){
#else
     Matrix ChangeBasis(const Matrix& k, const Matrix& m) throw (LapackError){
#endif

        Matrix tstiff(k);
        Matrix tmass(m);

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

	return tstiff;

#else

        // Compute Cholesky of mass matrix

        char uplo='U';
        int info;
        complex alpha=1.0;
        zpotrf_(&uplo,&tmass.dim[0],&(*tmass.data)[0],&tmass.dim[0],&info);

        // Factor Cholesky factor into stiffness matrix

        cblas_ztrsm(CblasColMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,tstiff.dim[0],tstiff.dim[0],&alpha,&(*tmass.data)[0],tmass.dim[0],&(*tstiff.data)[0],tstiff.dim[0]);

        cblas_ztrsm(CblasColMajor,CblasLeft,CblasUpper,CblasConjTrans,CblasNonUnit,tstiff.dim[0],tstiff.dim[0],&alpha,&(*tmass.data)[0],tmass.dim[0],&(*tstiff.data)[0],tstiff.dim[0]);

	return tstiff;

#endif

     }

     void NumRange(const Matrix& m, int n, std::string filename){
 // Compute the boundary of the numerical range of a matrix m and write the points to a
 // file given by 'filename'. n is the number of discretization points from 0 to pi.

       double step=PI/n;
       dvector eigs(2*n);
       dvector steps(2*n);
       complex im=complex(0,1);

#pragma omp parallel for
       for (int i=0;i<n;i++){
	 double theta=i*step;
	 steps[i]=theta; steps[n+i]=theta;
	 Matrix mtheta=std::exp(im*theta)*m;
	 // Now take the Hermitian Part
	 Matrix H=0.5*(mtheta+ConjTranspose(mtheta));
         pdvector pevalues;
	 HermitianEigenvalues(H,pevalues);
	 
	 double lmin=*(std::min_element(pevalues->begin(),pevalues->end()));
	 double lmax=*(std::max_element(pevalues->begin(),pevalues->end()));
	 eigs[i]=lmin; eigs[n+i]=lmax;
       }
	
       // Now iterate through to get the points of the Num. Range approximation

       
       steps.push_back(0);
       eigs.push_back(eigs[0]); // Append first element. Makes next for loop easier
       cvector numrange(2*n+1);
       for (int i=0;i<2*n;i++){
	 complex e1=std::exp(-im*steps[i]);
	 complex e2=std::exp(-im*steps[i+1]);
	 double c1=std::cos(steps[i]); double c2=std::cos(steps[i+1]);
	 double s1=std::sin(steps[i]); double s2=std::sin(steps[i+1]);
	 double lam1=eigs[i];
	 double lam2=eigs[i+1];
	 double t1=(lam2-lam1*(c1*c2+s1*s2))/(s1*c2-s2*c1);
	 numrange[i]=e1*(lam1+im*t1);
       }

	 // Add first point again for convenience
       numrange[2*n]=numrange[0];


       // Now give out results in file.

#ifdef BEM2DMPI
       BlacsSystem* b=BlacsSystem::Instance();
       if (b->IsRoot()){
#endif
	 std::ofstream out(filename.c_str());
	 for (int i=0;i<2*n+1;i++){ 
	   out << std::real(numrange[i]) << " " << std::imag(numrange[i]) << std::endl;
	 }
	 out.close();
	

#ifdef BEM2DMPI 
       }
#endif


     }

   Matrix ConjTranspose(const Matrix& m){


     int M=m.dim[0];
     int N=m.dim[1];
     Matrix result(N,M);

#ifdef BEM2DMPI
     char trans='C';
     int one=1;
     complex alpha=1.0;
     complex beta=1.0;

     pzgeadd_(&trans,&N,&M,&alpha,&(*m.data)[0],&one,&one,m.desc,&beta,&(*result.data)[0],&one,&one,
	      result.desc);
#else
     for (int i=0;i<M;i++){
       for (int j=0;j<N;j++){
	 (*result.data)[j+N*i]=std::conj((*m.data)[i+M*j]);
       }
     }

#endif

     return result;

   }

   Matrix ExtractColumn(const Matrix& m, int j){
     
     Matrix result(m.dim[0],1);
#ifdef BEM2DMPI
     BlacsSystem* b=BlacsSystem::Instance();
     int lj=b->G2lc(j);
     // Get Process of column and check if it is the first column of the process grid
     int p=(j/b->get_nb()) % b->get_npcol();
     if (p==0){
       for (int i=0;i<result.msize;i++) (*result.data)[i]=(*m.data)[lj*m.msize+i];
     }
     else{
       // Send data over to the first column of processors
       if (b->get_mycol()==0) Czgerv2d(b->get_ictxt(),result.msize,1,&(*result.data)[0],result.msize,b->get_myrow(),p);
	 if (b->get_mycol()==p) Czgesd2d(b->get_ictxt(),m.msize,1,&(*m.data)[lj*m.msize],m.msize,b->get_myrow(),0); 
     }
     return result;
#else
     for (int i=0;i<result.msize;i++) (*result.data)[i]=(*m.data)[j*m.msize+i];
     return result;
#endif
       }

   complex DotProduct(const Matrix& a, const Matrix& b) throw (ArrayMismatch){

     if ((a.nsize!=1)||(b.nsize!=1)) throw ArrayMismatch();
     if ((a.msize!=b.msize)) throw ArrayMismatch();
     complex dotc;
#ifdef BEM2DMPI
     int one=1;
     pzdotc_(&a.msize,&dotc,&(*a.data)[0],&one,&one,a.desc,&one,&(*b.data)[0],&one,&one,b.desc,&one);
     // Send result over to all other processors in the row
     BlacsSystem* blacs=BlacsSystem::Instance();
     char scope[]="All";
     char top[]=" ";
     
     if (blacs->IsRoot()){ 
       Czgebs2d(blacs->get_ictxt(),scope,top,1,1,&dotc,1);
     }
     else{
       Czgebr2d(blacs->get_ictxt(),scope,top,1,1,&dotc,1,0,0);
     }
#else
     cblas_zdotc_sub(a.msize,&(*a.data)[0],1,&(*b.data)[0],1,&dotc);
#endif
     return dotc;
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

