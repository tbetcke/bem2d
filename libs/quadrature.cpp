#include "quadrature.h"

namespace bem2d {

void gauss(dvector& x, dvector& w, int N) throw (lapack_error) {
    // Points xa and weights w for N point Gaussian quadrature
    // in the interval [0,1]


    /*
        double* beta;
        double* alpha;
        double* work;
        double* v;
     */
        char c = 'V';
        int info;


        x.clear(); w.clear(); // Make sure that x and w are empty

        if (N == 1) {
            x.push_back(0.5);
            w.push_back(1.0);
            return;
        }

        /* Allocate memory */

        /*
        alpha = (double*) calloc(N, sizeof (double));
        beta = (double*) calloc(N - 1, sizeof (double));
        v = (double*) malloc(N * N * sizeof (double));
        work = (double*) malloc((2 * N - 2) * sizeof (double));

        */

        dvector alpha(N,0);
        dvector beta(N-1,0);
        dvector v(N*N,0);
        dvector work(2*N-2,0);

        /* Fill the Array */

        for (int i = 0; i < N - 1; i++) {
            beta[i] = .5 / sqrt(1. - 1. / (4. * (i + 1)*(i + 1)));
        }

        /* Call Lapack */

        dstev_(&c, &N, &alpha[0], &beta[0], &v[0], &N, &work[0], &info);

        /* Return results */

        for (int i = 0; i < N; i++) {
            x[i] = (1.0 + alpha[i]) / 2.0;
            w[i] = v[i * N] * v[i * N];
        }

        /*
        free(alpha);
        free(beta);
        free(work);
        free(v);
        */

        if (info) throw lapack_error();

    }

}