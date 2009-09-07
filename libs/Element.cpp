#include "Element.h"

namespace bem2d {

    Element::~Element() {
    };

    ConstElement::ConstElement(const Point& pp1, const Point& pp2)
    : p1(pp1),
    d(pp2 - pp1) {
    }

    ConstElement::~ConstElement() {
    };

}
