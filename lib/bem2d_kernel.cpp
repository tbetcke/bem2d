#include "bem2d_kernel.h"
#include "bem2d_mathroutines.h"

namespace bem2d
{

SingleLayer::SingleLayer(freqtype k): k_(k)
{
}

SingleLayer::SingleLayer(const SingleLayer& s): k_(s.k())
{
}

complex SingleLayer::operator()(Point x, Point y) const
{

        double d=length(x-y);
        complex i(0,1);
        return i/4.0*(BesselH0(k_*d));


}

DoubleLayer::DoubleLayer(freqtype k): k_(k)
{
}

DoubleLayer::DoubleLayer(const DoubleLayer& d): k_(d.k())
{
}


complex DoubleLayer::operator()(Point x, Point y) const
{
        Point w=y-x;
        double d=length(w);
        complex i(0,1);

        return -i*k_/4.0*(BesselH1(k_*d)/d*(n2_.x*w.x+n2_.y*w.y));
}



ConjDoubleLayer::ConjDoubleLayer(freqtype k): k_(k)
{
}

ConjDoubleLayer::ConjDoubleLayer(const ConjDoubleLayer& cd): k_(cd.k())
{
}


complex ConjDoubleLayer::operator()(Point x, Point y) const
{
        Point w=x-y;
        double d=length(w);
        complex i(0,1);

        return -i*k_/4.0*(BesselH1(k_*d)/d*(n1_.x*w.x+n1_.y*w.y));
}

CombinedSingleConjDouble::CombinedSingleConjDouble(freqtype k, double eta):
                k_(k),
                eta_(eta)
{}

CombinedSingleConjDouble::CombinedSingleConjDouble(freqtype k): k_(k), eta_(k) {}

CombinedSingleConjDouble::CombinedSingleConjDouble(const CombinedSingleConjDouble& scdl):
                k_(scdl.k_), eta_(scdl.eta_) {}

complex CombinedSingleConjDouble::operator()(Point x, Point y) const
{
        Point w=x-y;
        double d=length(w);
        complex i(0,1);
        complex slayer=i/4.0*(BesselH0(k_*d));
        complex cdlayer=-i*k_/4.0*(BesselH1(k_*d)/d*(n1_.x*w.x+n1_.y*w.y));
        return cdlayer-i*eta_*slayer;
}

}