#ifndef QUADPOINTS_H_
#define QUADPOINTS_H_

#include "bem2ddefs.h"
#include "exceptions.h"
#include<cstdlib>


namespace bem2d {

struct QuadPoints1D{
	dvector x;
	dvector w;
	double a;
	double b;

};

struct QuadPoints2D{
	// Tensor Quadratur points in rectangle [a,b]x[c,d]
	// with coordinates (x[i],y[i]) and weights w[i]
	dvector x;
	dvector y;
	dvector w;
	double a;
	double b;
	double c;
	double d;

};


class Gauss1D {
public:
	Gauss1D(int N) throw (parameter_exception);
	// Create Gauss points between 0 and 1

	QuadPoints1D rescale(double a, double b) const throw (interval_exception);

	inline dvector::const_iterator getpoints() const {
		return x.begin();
	}

	inline dvector::const_iterator pointend() const {
		return x.end();
	}

	inline dvector::const_iterator getweights() const {
		return w.begin();
	}

	inline dvector::const_iterator weightend() const {
		return w.end();
	}

	inline int size() const {
		return N;
	}

	~Gauss1D();

private:
	dvector x;
	dvector w;
	int N;
};

class Gauss2D {
public:
	Gauss2D(int N) throw (parameter_exception);
	// Create N*N Tensor-Gauss points in [0,1]^2

	QuadPoints2D rescale(double a, double b, double c, double d) const throw (interval_exception);
	// Create Gauss Points in the rectangle [a,b]x[c,d]

	inline dvector::const_iterator getpointsx() const {
		return x.begin();
	}

	inline dvector::const_iterator getpointsy() const {
		return y.begin();
	}

	inline dvector::const_iterator pointendx() const {
		return x.end();
	}

	inline dvector::const_iterator pointendy() const {
		return y.end();
	}


	inline dvector::const_iterator getweights() const {
		return w.begin();
	}

	inline dvector::const_iterator weightend() const {
		return w.end();
	}

	inline int size() const {
		return N; // Number of tensor Gauss points
	}

	~Gauss2D();

private:
	dvector x;
	dvector y;
	dvector w;
	int N;
};

class AdaptedGauss1 {
	// Gauss Quadrature Points for a function that has singularities at the left side and lower side
	// of the box [0,1]x[0,1]
public:
	AdaptedGauss1(int N, int L, double sigma);

	inline dvector::const_iterator getpointsx() const {
		return x.begin();
	}

	inline dvector::const_iterator getpointsy() const {
		return y.begin();
	}

	inline dvector::const_iterator pointendx() const {
		return x.end();
	}

	inline dvector::const_iterator pointendy() const {
		return y.end();
	}


	inline dvector::const_iterator getweights() const {
		return w.begin();
	}

	inline dvector::const_iterator weightend() const {
		return w.end();
	}

	inline int size() const {
		return N; // Number of tensor Gauss points
	}

	~AdaptedGauss1();

private:
	dvector x;
	dvector y;
	dvector w;
	int N;
	int L;
	double sigma;
};

class AdaptedGauss2 {
	// Gauss Quadrature Points for a function that has singularities at the point 0
public:
	AdaptedGauss2(int N, int L, double sigma);

	inline dvector::const_iterator getpointsx() const {
		return x.begin();
	}

	inline dvector::const_iterator getpointsy() const {
		return y.begin();
	}

	inline dvector::const_iterator pointendx() const {
		return x.end();
	}

	inline dvector::const_iterator pointendy() const {
		return y.end();
	}


	inline dvector::const_iterator getweights() const {
		return w.begin();
	}

	inline dvector::const_iterator weightend() const {
		return w.end();
	}

	inline int size() const {
		return N; // Number of tensor Gauss points
	}

	~AdaptedGauss2();

private:
	dvector x;
	dvector y;
	dvector w;
	int N;
	int L;
	double sigma;
};


	class AdaptedGauss3 {
		// Gauss Quadrature Points for a function that has singularities at the point 0
	public:
		AdaptedGauss3(int N, int L, double sigma);
		
		inline dvector::const_iterator getpointsx() const {
			return x.begin();
		}
		
		inline dvector::const_iterator getpointsy() const {
			return y.begin();
		}
		
		inline dvector::const_iterator pointendx() const {
			return x.end();
		}
		
		inline dvector::const_iterator pointendy() const {
			return y.end();
		}
		
		
		inline dvector::const_iterator getweights() const {
			return w.begin();
		}
		
		inline dvector::const_iterator weightend() const {
			return w.end();
		}
		
		inline int size() const {
			return N; // Number of tensor Gauss points
		}
		
		~AdaptedGauss3();
		
	private:
		dvector x;
		dvector y;
		dvector w;
		int N;
		int L;
		double sigma;
	};
	
	
}
#endif /* QUADPOINTS_H_ */
