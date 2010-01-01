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

        // Scalapack functions

        void sl_init_(int* ictxt, int* nprow, int* npcol);
        int numroc_(int*, int*, int*, int*, int*);
        void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
        void pzgemm_(const char*,const char*,const int*,const int*,const int*,const complex*,const complex*,const int*,const int*,const int*,const complex*,const int*,const int*,const int*,const complex*, complex*,const int*,const int*,int*);
        void pzgetrf_(const int*, const int*, complex*, const int*, const int*, const int*, int*, int*);
        void pzgetrs_(const char*, const int*, const int*, const complex*, const int*, const int*, const int*, const int*, complex*, const int*, const int*, const int*, int*);
        int indxl2g_(const int*, const int*, const int*, const int*, const int*);
        void pzpotrf_(const char*, const int*, complex* ,const int*, const int*, const int*, int*);
        void pztrsm_(const char*, const char*, const char*, const char*, const int*, const int*, const complex*, const complex*, const int*, const int*, const int*, complex*, const int*, const int*, const int*);
        void pzgesvd_(const char*, const char*, const int*, const int*, complex*, const int*, const int*, const int*, double*, complex*, const int*, const int*, const int*, complex*, const int*, const int*, const int*, complex*, const int*, double*, int*);
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

