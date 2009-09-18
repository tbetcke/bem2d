#include "quadrature.h"
#include <tr1/array>
#include <cstdlib>

namespace bem2d {

void gauss(dvector& x, dvector& w, int N) throw (lapack_error) {
    // Points xa and weights w for N point Gaussian quadrature
    // in the interval [0,1]


        char c = 'V';
        int info;


        x.clear(); w.clear(); // Make sure that x and w are empty
        x.resize(N); w.resize(N);

        if (N == 1) {
            x[0]=0.5;
            w[0]=1.0;
            return;
        }

        dvector alpha(N,0);
        dvector beta(N-1,0);
        dvector v(N*N,0);
        dvector work(2*N-2,0);


        /* Fill the Array */

        for (std::size_t i = 0; i < N - 1; i++) {
            beta[i] = .5 / sqrt(1. - 1. / (4. * (i + 1)*(i + 1)));
        }

        /* Call Lapack */

        dstev_(&c, &N, &alpha[0], &beta[0], &v[0], &N, &work[0], &info);

        /* Return results */

        for (std::size_t i = 0; i < N; i++) {
            x[i] = (1.0 + alpha[i]) / 2.0;
            w[i] = v[i * N] * v[i * N];
        }

        if (info) throw lapack_error();

    }

void mappoints(dvector& x, dvector& w, double a, double b){

    for (std::size_t i = 0; i < x.size(); i++) {
        x[i] = a + (b - a) * x[i];
        w[i] = (b - a) * w[i];
    }

}

void mappoints2d(dvector& x, dvector& y, dvector& w, double a, double b, double c, double d){

	for (std::size_t i=0; i< x.size(); i++){
		x[i]=a+(b-a)*x[i];
		y[i]=c+(d-c)*y[i];
		w[i]=(b-a)*(d-c)*w[i];
	}
}


}
