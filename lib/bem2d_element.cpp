#include "bem2d_element.h"

namespace bem2d
{

Element::~Element()
{
};

ConstElement::ConstElement(const Point& p1, const Point& p2,std::size_t index)
                : p1_(p1),
                d_(p2 - p1)
{
        set_index(index);
}

ConstElement::~ConstElement()
{
};

AnalyticCurveElement::AnalyticCurveElement(double t1, double t2, pCurve curve, std::size_t index):
                tstart_(t1),
                tend_(t2),
                curve_(curve)
{
        set_index(index);
}


}

