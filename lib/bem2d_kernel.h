#ifndef _KERNEL_H
#define _KERNEL_H

#include <complex>
#include "boost/utility.hpp"
#include "bem2d_point.h"
#include "bem2d_defs.h"

namespace bem2d
{

class SingleLayer
{
public:
  SingleLayer(freqtype k);
  SingleLayer(const SingleLayer& s);
  complex operator ()(Point x, Point y) const;
  inline void SetNormal(Point normal1, Point normal2)
  {
    n1_=normal1;
    n2_=normal2;
  }
  inline freqtype k() const
  {
    return k_;
  }

private:
  freqtype k_;
  Point n1_;
  Point n2_;
};

class DoubleLayer
{
public:
  DoubleLayer(freqtype k);
  DoubleLayer(const DoubleLayer& d);
  complex operator ()(Point x, Point y) const;
  inline void SetNormal(Point normal1, Point normal2)
  {
    n1_=normal1;
    n2_=normal2;
  }
  inline freqtype k() const
  {
    return k_;
  }
private:
  freqtype k_;
  Point n1_;
  Point n2_;
};

class ConjDoubleLayer
{
public:
  ConjDoubleLayer(freqtype k);
  ConjDoubleLayer(const ConjDoubleLayer& d);
  complex operator ()(Point x, Point y) const;
  inline void SetNormal(Point normal1, Point normal2)
  {
    n1_=normal1;
    n2_=normal2;
  }
  inline freqtype k() const
  {
    return k_;
  }
private:
  freqtype k_;
  Point n1_;
  Point n2_;
};

class CombinedSingleConjDouble
{
public:
  CombinedSingleConjDouble(freqtype k, double eta);
  CombinedSingleConjDouble(freqtype k);
  CombinedSingleConjDouble(const CombinedSingleConjDouble& scdl);
  complex operator()(Point x, Point y) const;
  inline void SetNormal(Point normal1, Point normal2)
  {
    n1_=normal1;
    n2_=normal2;
  }
  inline freqtype getk() const
  {
    return k_;
  }
private:
  freqtype k_;
  Point n1_;
  Point n2_;
  double eta_;
};
	
	class IdKernel{
	public:
		inline void SetNormal(Point normal1, Point normal2){}
		complex operator()(Point x, Point y) const{
			return 1.0;
		}
	};


}


#endif // _KERNEL_H