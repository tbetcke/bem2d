#include "outputroutines.h"
#include<fstream>
#include "outputroutines.h"

namespace bem2d {
	
	
	void abrange(dvector& xx, double a, double b, int n){
		xx.resize(n);
		for (int j = 0; j < n; j++) xx[j] = a + (1.0 / (n - 1)) * j * (b - a);
	}
	
	boost::shared_ptr<std::vector<Point > > meshgrid(double ax, double bx,
													 double ay, double by, int xpts, int ypts){
		
		dvector xx;
		dvector yy;
		
		abrange(xx,ax,bx,xpts);
		abrange(yy,ay,by,ypts);
		
		boost::shared_ptr<std::vector<Point > > pvec(new std::vector<Point>);
		pvec->reserve(xpts*ypts);
		
        // Create the grid
		
        for (int i = 0; i < xpts; i++) {
            for (int j = 0; j < ypts; j++) {
				pvec->push_back(Point(xx[i],yy[j]));
            }
        }
		
		return pvec;
    }
	
	
	void gplotout(std::string name, const std::vector<Point>& points, const dvector& z,
				  int xpts, int ypts){
		
		std::ofstream out(name.c_str());
		
		for (int i=0;i<xpts;i++){
			for (int j=0;j<ypts;j++){
				out << points[i*ypts+j].x << " " << points[i*ypts+j].y << " " << z[i*ypts+j] << "\n";
			}
			out << "\n";
		}
		out.close();
		
		// Now write out a corresponding gnuplot script
		
		std::string scriptname=name + "_gplot";
		
		std::ofstream out2(scriptname.c_str());
		
		out2 << "set pm3d map\n";
		out2 << "set hidden3d\n";
		out2 << "set palette rgbformulae 33,13,10\n";
		out2 << "set size ratio -1\n";
		out2 << "set size square\n";
		out2 << "set cbrange [-1:1]\n";
		out2 << "splot \'" +  name + "\'\n";
		out2 << "pause -1\n";
		out2.close();
		
	}
	
}