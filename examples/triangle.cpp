#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>


int main(int argc, char** argv)
{

  bem2d::freqtype k={10,0};
        int n=3*k.re;




        std::vector<bem2d::Point> triangle;
        triangle.push_back(bem2d::Point(0,0));
        triangle.push_back(bem2d::Point(1,0));
        triangle.push_back(bem2d::Point(.5,.5*sqrt(3)));
        bem2d::Polygon poly(triangle,n);
        bem2d::pGeometry pgeom=poly.GetGeometry();



        bem2d::PolBasis::AddBasis(0,pgeom);


        bem2d::Point direction=bem2d::normalize(bem2d::Point(1,1));
        bem2d::PlaneWave pw(direction,k);
        bem2d::CombinedPlaneWave cpw(direction,k);

        bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
        soundsoft.SetQuadOption(5,5,0.15);
        soundsoft.set_polygons(pgeom);
        soundsoft.set_plotInterior();



        std::cout << "Discretizing with k=" << bem2d::complex(k.re,k.im) << " and n=" << pgeom->size() << std::endl;
        clock_t start, finish;
        double time;

        start=clock();
        soundsoft.Discretize();
        finish=clock();
        time=(double(finish)-double(start))/CLOCKS_PER_SEC;
        std::cout << "Computing time (minutes): " << time/60 << std::endl;



        int xpts=300;
        int ypts=300;
        bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-4,8,-4,8,"wedge"));
        soundsoft.SetOutput(pout);
        soundsoft.WriteAll();


}



