#ifndef _CURVE_H_
#define _CURVE_H_

#include "bem2d_point.h"
#include "bem2d_defs.h"
#include<cmath>

namespace bem2d
{

  class Curve
  {
  public:
    virtual Point Map(double t) const=0;
    virtual Point Deriv(double t) const=0;
    virtual ~Curve(){};
  };

typedef boost::shared_ptr<Curve> pCurve;


}


#endif // _CURVE_H_

