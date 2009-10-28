#ifndef _OUTPUTHANDLER_H_
#define _OUTPUTHANDLER_H_

#include "bem2ddefs.h"
#include<vector>
#include "Point.h"
#include <boost/shared_ptr.hpp>
#include <string>
#include "Point.h"

namespace bem2d {
	
	class Outputhandler {
	public:
		virtual ~Outputhandler(){};
		virtual void writeIncident(const cvector& vals){};
		virtual void writeScattered(const cvector& vals){};
		virtual void writeFull(const cvector& vals){};
		virtual void writeAll(const cvector& valsincident, const cvector& valsscattered, const cvector& valsfull){};
		inline const std::vector<Point>& getmesh() const {
			return mesh;
		}
		inline void setreal(){
			real=true;
		}
		inline void setimag(){
			real=false;
		}
		
	protected:
		pdvector realdata(const cvector& vals);
		pdvector imagdata(const cvector& vals);
		pdvector turntorealimag(const cvector& vals);
		std::vector<Point> mesh;
		bool real;
	};

	typedef boost::shared_ptr<Outputhandler> pOutputhandler;
	
	class Gplotoutput: public Outputhandler {
	public:
		Gplotoutput(int xptsvalue, int yptsvalue, int ax, int bx, int ay, int by, std::string fname);
		void writeIncident(const cvector& vals);
		void writeScattered(const cvector& vals);
		void writeFull(const cvector& vals);
		void writeAll(const cvector& valsincident, const cvector& valsscattered, const cvector& valsfull);		
	private:
		std::string name;
		int xpts; int ypts;
	};
	
}

#endif // _OUTPUTHANDLER_H_