#ifndef _CONSTBASIS_H
#define	_CONSTBASIS_H

namespace bem2d{

class ConstBasis {
public:
    inline static double evaluate(double t, int n){
        return 1.0;
    }

};
}
#endif	/* _CONSTBASIS_H */

