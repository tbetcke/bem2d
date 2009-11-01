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

}
