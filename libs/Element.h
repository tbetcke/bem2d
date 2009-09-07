#ifndef _ELEMENT_H
#define	_ELEMENT_H

#include <boost/utility.hpp>
#include "Point.h"

namespace bem2d {

    class Element : private boost::noncopyable {
    public:
        virtual const Point map(double t) const = 0;
        virtual const Point deriv(double t) const = 0;

        inline const Point first() {
            return map(0);
        };

        inline const Point last() {
            return map(1);
        };

        inline const Point normal(double t) const{
            Point p(deriv(t));
            Point pn(p.y, -p.x);

            return normalize(pn);
        }


        virtual ~Element();
    private:

    };

    class ConstElement : public Element {
    public:
        ConstElement(const Point& pp1, const Point& pp2);

        inline const Point map(double t) const {
            return p1 + d*t;
        };

        inline const Point deriv(double t) const {
            return d;
        };
        ~ConstElement();
    private:
        const Point p1; // Start Point
        const Point d; // End point - Start Point

    };

}
#endif	/* _ELEMENT_H */

