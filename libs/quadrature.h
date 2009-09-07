#ifndef _QUADRATURE_H
#define	_QUADRATURE_H

#include "bem2ddefs.h"
#include "exceptions.h"

namespace bem2d {

    void gauss(dvector& x, dvector& w, int N) throw (lapack_error);

}
#endif	/* _QUADRATURE_H */

