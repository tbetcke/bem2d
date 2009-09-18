#ifndef _QUADRATURE_H
#define	_QUADRATURE_H

#include "bem2ddefs.h"
#include "exceptions.h"

namespace bem2d {

    void gauss(dvector& x, dvector& w, int N) throw (lapack_error);
    void mappoints(dvector& x, dvector& w, double a, double b);
    void mappoints2d(dvector& x, dvector& y, dvector& w, double a, double b, double c, double d);

}
#endif	/* _QUADRATURE_H */

