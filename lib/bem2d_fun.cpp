#include "bem2d_fun.h"

namespace bem2d
{


PlaneWave::PlaneWave(Point direction, freqtype k): direction_(direction), k_(k) {}
PlaneWave::PlaneWave(const PlaneWave& p): k_(p.k_), direction_(p.direction_) {}

NormalPlaneWave::NormalPlaneWave(Point direction, freqtype k):
                direction_(direction), k_(k) {}

NormalPlaneWave::NormalPlaneWave(const NormalPlaneWave& np):
                direction_(np.direction_), k_(np.k_) {}

CombinedPlaneWave::CombinedPlaneWave(Point direction, freqtype k, double eta):
                direction_(direction), k_(k), eta_(eta) {}

CombinedPlaneWave::CombinedPlaneWave(Point direction, freqtype k):
                direction_(direction), k_(k), eta_(k) {}


CombinedPlaneWave::CombinedPlaneWave(const CombinedPlaneWave& np):
                direction_(np.direction_), k_(np.k_), eta_(np.eta_) {}


OutWave::OutWave(freqtype k): k_(k), s_(k) {}

OutWave::OutWave(const OutWave& owave): k_(owave.k_), s_(owave.k_) {};

complex OutWave::operator()(Point p) const
{
        Point zero(0,0);
        return s_(zero,p);
}

}