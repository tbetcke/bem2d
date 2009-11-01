#ifndef _ELEMENT_H
#define	_ELEMENT_H

#include <cstdlib>
#include "boost/utility.hpp"
#include "bem2d_defs.h"
#include "bem2d_point.h"


namespace bem2d {

    class Element : private boost::noncopyable {
    public:
        virtual const Point Map(double t) const = 0;
        virtual const Point Deriv(double t) const = 0;

        inline const Point First() {
            return Map(0);
        };

        inline const Point Last() {
            return Map(1);
        };

        inline const Point Normal(double t) const{
            Point p(Deriv(t));
            Point pn(p.y, -p.x);

            return normalize(pn);
        }
		
		inline void set_index(std::size_t index){
			index_=index;
		}
		
		inline std::size_t index() const {
			return index_;
		}

		inline void set_next(std::size_t index){
			next_=index;
		}
		
		inline std::size_t next() const {
			return next_;
		}
		
		inline void set_prev(std::size_t index){
			prev_=index;
		}
		
		inline std::size_t prev() const {
			return prev_;
		}
		

        virtual ~Element();
    private:
		std::size_t index_;
		std::size_t next_;
		std::size_t prev_;

    };
	
	
	typedef boost::shared_ptr<Element> pElement;

    class ConstElement : public Element {
    public:
        ConstElement(const Point& p1, const Point& p2, std::size_t index);

        inline const Point Map(double t) const {
            return p1_ + d_*t;
        };

        inline const Point Deriv(double t) const {
            return d_;
        };
        ~ConstElement();
    private:
        const Point p1_; // Start Point
        const Point d_; // End point - Start Point

    };
		
	template<typename T>
	class AnalyticCurveElement: public Element {
	public:
		AnalyticCurveElement(double t1, double t2, const T& curve, std::size_t index);
		
		inline const Point Map(double t) const {
			return curve_.Map(tstart_+t*(tend_-tstart_));
		}
		inline const Point Deriv(double t) const {
			return (tend_-tstart_)*curve_.Deriv(tstart_+t*(tend_-tstart_));
		}
		
	private:
		T curve_;
		double tstart_;
		double tend_;
	};
	
	
	template<typename T>
	AnalyticCurveElement<T>::AnalyticCurveElement(double t1, double t2, const T& curve, std::size_t index):
	tstart_(t1),
	tend_(t2),
	curve_(curve){
		set_index(index);
	}

}
#endif	/* _ELEMENT_H */

