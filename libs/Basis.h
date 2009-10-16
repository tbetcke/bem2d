#ifndef _BASIS_H
#define	_BASIS_H

#include<cstdlib>

namespace bem2d{
	
	class Basis {
	public:
		virtual double operator()(double t)=0;
		virtual ~Basis();
		inline void setIndex(std::size_t n)
		{
			_index=n;
		}
		inline std::size_t getIndex(std::size_t n){
			return _index;
		}
	private:
		std::size_t _index;
	};
		
	
	
	class ConstBasis: public Basis {
	public:
		inline double operator()(double t){
			return 1.0;
		}
		
	};
}
#endif	/* _BASIS_H */

