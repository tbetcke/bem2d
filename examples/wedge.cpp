#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>


int main(int argc, char** argv)
{

        bem2d::freqtype k=10;
        int n=20*k;




        std::vector<bem2d::Point> wedge;
        wedge.push_back(bem2d::Point(0,0));
        wedge.push_back(bem2d::Point(6,0));
        wedge.push_back(bem2d::Point(1,1));
        wedge.push_back(bem2d::Point(0,6));
        bem2d::Polygon poly(wedge,n);
        bem2d::pGeometry pgeom=poly.GetGeometry();



        bem2d::PolBasis::AddBasis(0,pgeom);


        bem2d::Point direction=bem2d::normalize(bem2d::Point(1,1));
        bem2d::PlaneWave pw(direction,k);
        bem2d::CombinedPlaneWave cpw(direction,k);

        bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
        soundsoft.SetQuadOption(5,5,0.15);
        soundsoft.set_polygons(pgeom);
        soundsoft.set_plotInterior();



        std::cout << "Discretizing with k=" << k << " and n=" << pgeom->size() << std::endl;
        clock_t start, finish;
        double time;

        start=clock();
        soundsoft.Discretize();
        finish=clock();
        time=(double(finish)-double(start))/CLOCKS_PER_SEC;
        std::cout << "Computing time (minutes): " << time/60 << std::endl;

        // Output the matrix


        int xpts=300;
        int ypts=300;
        bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-4,8,-4,8,"wedge"));
        soundsoft.SetOutput(pout);
        soundsoft.WriteAll();


}
