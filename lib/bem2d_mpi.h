#ifndef _BEM2D_MPI_H
#define _BEM2D_MPI_H

#include "mpi.h"
#include "bem2d_defs.h"

// All MPI/ScaLapack related functions

// External functions

namespace bem2d
{

extern "C" {


        // Blacs Functions

        extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
        extern void   Cblacs_get( int context, int request, int* value);
        extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
        extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
        extern void   Cblacs_gridexit( int context);
        extern void   Cblacs_exit( int error_code);
        extern void   Czgsum2d(int icontxt, char* scope, char* top, int m, int n, complex* A, int lda, int rdest, int cdest);
  extern void   Czgesd2d(int icontxt, int m, int n, complex* A, int lda, int rdest, int cdest);
  extern void   Czgerv2d(int icontxt, int m, int n, complex* A, int lda, int rsrc, int csrc);
  extern void   Czgebs2d(int icontxt, char *scope, char *top, int m, int n, complex* A, int lda);
  extern void   Czgebr2d(int icontxt, char *scope, char *top, int m, int n, complex* A, int lda, int rsrc, int csrc);

        // Scalapack functions

        void sl_init_(int* ictxt, int* nprow, int* npcol);
        int numroc_(int*, int*, int*, int*, int*);
        void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
        void pzgemm_(const char*,const char*,const int*,const int*,const int*,const complex*,const complex*,const int*,const int*,const int*,const complex*,const int*,const int*,const int*,const complex*, complex*,const int*,const int*,int*);
        void pzgetrf_(const int*, const int*, complex*, const int*, const int*, const int*, int*, int*);
        void pzgetrs_(const char*, const int*, const int*, const complex*, const int*, const int*, const int*, const int*, complex*, const int*, const int*, const int*, int*);
        int indxl2g_(const int*, const int*, const int*, const int*, const int*);
  int indxg2l_(const int*, const int*, const int*, const int*, const int*);
        void pzpotrf_(const char*, const int*, complex* ,const int*, const int*, const int*, int*);
        void pztrsm_(const char*, const char*, const char*, const char*, const int*, const int*, const complex*, const complex*, const int*, const int*, const int*, complex*, const int*, const int*, const int*);
        void pzgesvd_(const char*, const char*, const int*, const int*, complex*, const int*, const int*, const int*, double*, complex*, const int*, const int*, const int*, complex*, const int*, const int*, const int*, complex*, const int*, double*, int*);
  void pzhegvx_(const int* ibtype, const char* jobz, const char* range, const char* uplo, const int* N, 
		complex* A, const int* IA, const int* JA, const int* desca, complex* B, const int* ib,
		const int* jb, const int* descb, const double* vl, const double* vu, const int* il,
		const int* iu, const double* abstol, int* M, int* NZ, double* w, const double* orfac,
		complex* Z, const int* iz, const int* jz, const int* descz, complex* work, const int* lwork,
		double* rwork, const int* lrwork, int* iwork, const int* liwork, int* ifail, int* icluster,
		double* gap, int* info);
  double pdlamch_(const int* ictxt, const char* cmach);

  void pzgeadd_(const char* trans, const int* M, const int* N, const complex* alpha, const complex* A, 
		const int* ia, const int* ja, const int* desca, const complex* beta, complex* C, const int* ic, 
		const int* jc, const int* descc);
  void pzdotc_(const int* n, complex* dotc, const complex* x, const int* ix, const int* jx, 
	       const int* descx,
	       const int* incx, const complex* y, const int* iy, const int* jy, const int* descy,
	       const int* incy);
  void pzgehrd_(const int* n, const int* ilo, const int* ihi, complex* A, const int* ia,
		const int* ja, const int* desca, complex* tau, complex* work,
		const int* lwork, int* info);
  void pzlahqr_(const int* wantt, const int* wantz, const int* n, const int* ilo,
		const int* ihi, complex* A, const int* desca, complex* w, const int* iloz,
	        const int* ihiz, complex* Z, const int* descz, complex* work, const int* lwork,
		int* iwork, const int* ilwork, int* info);
  void pzheev_(const char* jobz, const char* uplo, const int* n, complex* A, const int* ia,
	       const int* ja, const int* desca, double* w, complex* Z, const int* iz,
	       const int* jz, const int* descz, complex* work, const int* lwork,
	       double* rwork, const int* lrwork, int* info);
  void pzheevx_(const char* jobz, const char* range, const char* uplo,
		const int* n, complex* A, const int* ia, const int* ja,
		const int* desca, const double* vl, const double* vu,
		const int* il, const int* iu, const double* abstol,
		int* m, int* nz, double* w, const double* orfac,
		complex* Z, const int* iz, const int* jz, const int* descz,
		complex* work, int* lwork, double* rwork,
		const int* lrwork, int* iwork, const int* liwork,
		int* ifail, int* icluster, double* gap, int* info);
		

}


// Define Singleton that initializes the MPI/BLACS system and stores the process grid information

class BlacsSystem
{
public:
        static BlacsSystem* Initialize(int nprow, int npcol, int mb, int nb) {
                if (!pInstance_) {
                        int nprocs;
                        MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
                        if (nprow*npcol>nprocs) return pInstance_;
                        pInstance_=new BlacsSystem(nprow, npcol, mb, nb);
                }
                return pInstance_;
        }

        static BlacsSystem* Instance() {
                return pInstance_;
        }


        static void Release() {
                if (!pInstance_) return;
                int ictxt=BlacsSystem::Instance()->get_ictxt();
                delete pInstance_;
                pInstance_=0;
                Cblacs_gridexit(ictxt);
        }

        int MSize(int m);
        // Return number of rows stored locally

        int NSize(int n);
        // Return number of columns stored locally

        inline int get_ictxt() {
                return ictxt_;
        }
        inline int get_nprow() {
                return nprow_;
        }
        inline int get_npcol() {
                return npcol_;
        }
        inline int get_mpiprocs() {
                return mpiprocs_;
        }

        inline int get_myrow() {
                return myrow_;
        }
        inline int get_mycol() {
                return mycol_;
        }
        inline int get_mpirank() {
                return mpirank_;
        }
        inline int get_mb() {
                return mb_;
        }
        inline int get_nb() {
                return nb_;
        }
        inline int L2gr(int m) {
                // Map local row index m to global row index
                // (counting from zero)
                int zero=0;
                int mp1=m+1;
                return indxl2g_(&mp1,&mb_,&myrow_,&zero,&nprow_)-1;
        }
        inline int L2gc(int n) {
                // Map local row index m to global row index
                // (counting from zero)
                int zero=0;
                int np1=n+1;
                return indxl2g_(&np1,&nb_,&mycol_,&zero,&npcol_)-1;
        }
	inline int G2lr(int m){
	  // Map global row index m to local row index
	  // return -1 if index not on local processor
	  int p= (m/mb_) % nprow_;
	  if (p!=myrow_) return -1;
	  int l=m/(nprow_*mb_);
	  int x=m % mb_;
	  return l*mb_+x;

	}
	inline int G2lc(int n){
	  // Map global row index m to local row index
	  // return -1 if index not on local processor
	  int p=(n/nb_) % npcol_;
	  if (p!=mycol_) return -1;
	  int l=n/(npcol_*nb_);
	  int x=n % nb_;
	  return l*nb_+x;
	}


        inline bool IsRoot() {
                return ((myrow_==0)&&(mycol_==0));
        }


private:
        BlacsSystem();
        BlacsSystem(int nprow, int npcol, int mb, int nb);
        BlacsSystem(const BlacsSystem&);
        BlacsSystem& operator=(const BlacsSystem&);
        ~BlacsSystem() {};

        static BlacsSystem* pInstance_;

        int ictxt_;

        int nprow_;
        int npcol_;
        int mpiprocs_;

        int myrow_;
        int mycol_;
        int mpirank_;

        int mb_;
        int nb_;

};
}

#endif

