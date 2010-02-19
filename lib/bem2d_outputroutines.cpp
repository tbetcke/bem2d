#include<fstream>
#include<complex>
#include<vector>
#include<map>
#include "bem2d_outputroutines.h"

namespace bem2d
{


void Abrange(dvector& xx, double a, double b, int n)
{
        xx.resize(n);
        for (int j = 0; j < n; j++) xx[j] = a + (1.0 / (n - 1)) * j * (b - a);
}

boost::shared_ptr<std::vector<Point > > MeshGrid(double ax, double bx,
                double ay, double by, int xpts, int ypts)
{

        dvector xx;
        dvector yy;

        Abrange(xx,ax,bx,xpts);
        Abrange(yy,ay,by,ypts);

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


void GplotOut(std::string name, const std::vector<Point>& points, const dvector& z,
              int xpts, int ypts)
{

        std::ofstream out(name.c_str());

        for (int i=0; i<xpts; i++) {
                for (int j=0; j<ypts; j++) {
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
        out2 << "set datafile missing \"nan\"\n";
        out2 << "set palette rgbformulae 33,13,10\n";
        out2 << "set size ratio -1\n";
        out2 << "set size square\n";
        out2 << "set cbrange [-1:1]\n";
        out2 << "splot \'" +  name + "\'\n";
        out2 << "pause -1\n";
        out2.close();

}

void WriteMatrix(std::string fname, const Matrix m)
{
        std::string fr=fname+"_real";
        std::string fi=fname+"_imag";
        std::ofstream outr(fr.c_str());
        std::ofstream outi(fi.c_str());

        for (int i=0; i<m.dim[0]; i++) {
                for (int j=0; j<m.dim[0]; j++) {
                        outr << real((*m.data)[j*m.dim[0]+i]) << " ";
                        outi << imag((*m.data)[j*m.dim[0]+i]) << " ";

                }
                outr << std::endl;
                outi << std::endl;
        }
        outr.close();
        outi.close();
}

  void WriteDensity(std::string name, const Matrix& m, pGeometry pgeom, int npoints)
{

  // Get the elements map from the geometry
  std::vector<pElement> elements=pgeom->elements();
  
  // Get the flat basis map from the geometry
  boost::shared_ptr<Geometry::flat_basis_map> pflatmap=pgeom->FlatMap();
  int nbases=pflatmap->size();
  int nelem=elements.size();

  // Create a keymap from element ids to entries in the elements map

  std::map<std::size_t,std::size_t> elemkeys;

  for (int i=0;i<elements.size();i++) 
    elemkeys[elements[i]->index()]=i;

  std::vector<std::vector<Point> > points(nelem); // Stores the evaluation
                                                  // points for each element

  std::vector<std::vector<cvector>  > vals(nelem); // Stores the associated values
  

  for (int i=0;i<pflatmap->size();i++)
    {

      pElement elem=(*pflatmap)[i].first;
      pBasis base=(*pflatmap)[i].second;
      int ind=elemkeys[elem->index()];
      // Initialize if element was not yet previously used
      if (points[ind].size()==0)
	{
	  for (int j=0;j<npoints;j++) points[ind].push_back(elem->Map(j/npoints));
	  vals[ind].resize(npoints);
	  for (int j=0;j<npoints;j++) vals[ind][j].resize(m.dim[1]);
	}
      // Now add up the values
      for (int j=0;j<npoints;j++)
	for (int t=0;t<m.dim[1];t++)
	  vals[ind][j][t]+=(*base)(j/npoints)*(*m.data)[t*m.dim[0]+i];
    }

  // Now write data into a file

  std::ofstream out(name.c_str());

  for (int i=0;i<nelem;i++)
    for (int j=0;j<npoints;j++)
      {
	out << points[i][j].x << " " << points[i][j].y << " ";
	for (int t=0;t<m.dim[1];t++)
	  out << real(vals[i][j][t]) << " ";
	for (int t=0;t<m.dim[1];t++)
	  out << imag(vals[i][j][t]) << " ";
	out << std::endl;
      }

  out.close();
	
}


}

