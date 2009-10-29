#ifndef _ELEMENT_H
#define	_ELEMENT_H

#include <boost/utility.hpp>
#include "Point.h"
#include<boost/shared_ptr.hpp>
#include <cstdlib>


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
		
		inline void setIndex(std::size_t ind){
			index=ind;
		}
		
		inline std::size_t getIndex() const {
			return index;
		}

		inline void setNext(std::size_t ind){
			next=ind;
		}
		
		inline std::size_t getNext() const {
			return next;
		}
		
		inline void setPrev(std::size_t ind){
			prev=ind;
		}
		
		inline std::size_t getPrev() const {
			return prev;
		}
		

        virtual ~Element();
    private:
		std::size_t index;
		std::size_t next;
		std::size_t prev;

    };
	
	
	typedef boost::shared_ptr<Element> pElement;

    class ConstElement : public Element {
    public:
        ConstElement(const Point& pp1, const Point& pp2, std::size_t index);

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
	
	template<typename T>
	class AnalyticCurveElement: public Element {
	public:
		AnalyticCurveElement(double t1, double t2, const T& curveobj, std::size_t index);
		
		inline const Point map(double t) const {
			return curve.map(tstart+t*(tend-tstart));
		}
		inline const Point deriv(double t) const {
			return (tend-tstart)*curve.deriv(tstart+t*(tend-tstart));
		}
		
	private:
		T curve;
		double tstart;
		double tend;
	};
	
	template<typename T>
	AnalyticCurveElement<T>::AnalyticCurveElement(double t1, double t2, const T& curveobj, std::size_t index):
	tstart(t1),
	tend(t2),
	curve(curveobj){
		setIndex(index);
	}

}
#endif	/* _ELEMENT_H */

