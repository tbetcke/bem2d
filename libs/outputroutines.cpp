#include "outputroutines.h"
#include<fstream>
#include "outputroutines.h"

namespace bem2d {
	
	void gplotout(std::string name, std::vector<Point>& points, dvector& z,
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
		out2 << "set size ratio -1\n";
		out2 << "set size square\n";
		out2 << "set palette gray\n";
		out2 << "splot \'" +  name + "\'\n";
		out2 << "pause -1\n";
		out2.close();
		
	}
	
}