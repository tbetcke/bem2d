#ifndef _CURVETYPES_H_
#define _CURVETYPES_H_

#include "bem2d_point.h"
#include "bem2d_defs.h"
#include "bem2d_curve.h"
#include<cmath>

namespace bem2d
{


  class Circle: public Curve
{
public:
        inline Point Map(double t) const {
                return Point(cos(2*PI*t),sin(2*PI*t));
        }
        inline Point Deriv(double t) const {
                return 2*PI*Point(-1.0*sin(2*PI*t),cos(2*PI*t));
        }
};


  class Trefoil: public Curve
{
public:
        inline Point Map(double t) const {
                double r=1+0.3*cos(3*2*PI*t);
                return Point(r*cos(2*PI*t),r*sin(2*PI*t));
        }
        inline Point Deriv(double t) const {
                double r=1+0.3*cos(3*2*PI*t);
                double rp=-0.3*3*2*PI*sin(3*2*PI*t);
                return Point(rp*cos(2*PI*t)-r*2*PI*sin(2*PI*t),rp*sin(2*PI*t)+r*2*PI*cos(2*PI*t));
        }
};

  class Kite: public Curve
{
 public:
  inline Point Map(double t) const {
    return Point(cos(2*PI*t)+0.65*cos(4*PI*t)-0.65,1.5*sin(2*PI*t));
  }
  inline Point Deriv(double t) const {
    return Point(-2*PI*sin(2*PI*t)-0.65*4*PI*sin(4*PI*t),1.5*2*PI*cos(2*PI*t));
  }
};

  class Ellipse: public Curve
  {
  public:
    Ellipse(double rho);
    inline Point Map(double t) const {
      return Point(rho_*cos(2*PI*t),1/rho_*sin(2*PI*t));
    }
    inline Point Deriv(double t) const {
      return Point(-rho_*2*PI*sin(2*PI*t),2*PI/rho_*cos(2*PI*t));
    }
  private:
    double rho_;
  };

}


#endif // _CURVETYPES_H_

