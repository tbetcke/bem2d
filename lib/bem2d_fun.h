#ifndef _BEM2DFUN_H_
#define _BEM2DFUN_H_

#include <complex>
#include "boost/shared_ptr.hpp"
#include "bem2d_defs.h"
#include "bem2d_point.h"
#include "bem2d_kernel.h"

namespace bem2d
{


class PlaneWave
{
public:
  PlaneWave() {};
  PlaneWave(Point direction, freqtype k);
  PlaneWave(const PlaneWave& p);

  inline complex operator()(Point p) const
  {
    complex i(0,1);
    return std::exp(k_*i*(direction_.x*p.x+direction_.y*p.y));
  }
  void SetNormal(Point normal) {};
  inline freqtype k() const
  {
    return k_;
  }
  inline Point direction() const
  {
    return direction_;
  }
private:
  Point direction_;
  freqtype k_;
};

class NormalPlaneWave
{
public:
  NormalPlaneWave() {};
  NormalPlaneWave(Point direction, freqtype k);
  NormalPlaneWave(const NormalPlaneWave& np);

  inline complex operator()(Point p) const
  {
    complex i(0,1);
    complex f=std::exp(k_*i*(direction_.x*p.x+direction_.y*p.y));
    return i*k_*f*(n_.x*direction_.x+n_.y*direction_.y);
  }
  inline void SetNormal(Point normal)
  {
    n_=normal;
  }
  inline freqtype k() const
  {
    return k_;
  }
  inline Point direction() const
  {
    return direction_;
  }
private:
  Point direction_;
  freqtype k_;
  Point n_;
};

class CombinedPlaneWave
{
public:
  CombinedPlaneWave(Point direction, freqtype k, double eta);
  CombinedPlaneWave(Point direction, freqtype k);
  CombinedPlaneWave(const CombinedPlaneWave& np);

  inline complex operator()(Point p) const
  {
    complex i(0,1);
    complex f=std::exp(k_*i*(direction_.x*p.x+direction_.y*p.y));
    return i*k_*f*(n_.x*direction_.x+n_.y*direction_.y)+i*eta_*f;
  }
  inline void SetNormal(Point normal)
  {
    n_=normal;
  }
  inline freqtype k() const
  {
    return k_;
  }
  inline Point direction() const
  {
    return direction_;
  }
  inline double eta() const
  {
    return eta_;
  }
private:
  Point direction_;
  freqtype k_;
  Point n_;
  double eta_;
};




class OutWave
{
public:
  OutWave(freqtype k);
  OutWave(const OutWave& owave);
  complex operator()(Point x) const;
  void SetNormal(Point normal) {};
  inline freqtype k() const
  {
    return k_;
  }
private:
  freqtype k_;
  SingleLayer s_;
};

class IdFun
{
public:
  inline complex operator()(Point p) const
  {
    return 1.0;
  }
  void SetNormal(Point normal) {};
};

}

#endif // _INCIDENTFUN_H_
