#ifndef _BASIS_H_
#define _BASIS_H_

#include "bem2d_defs.h"

namespace bem2d
{

class Basis
{
public:
  virtual complex operator()(double t) const=0;
  virtual ~Basis();
};

typedef boost::shared_ptr<Basis> pBasis;

}


#endif // _BASIS_H_
