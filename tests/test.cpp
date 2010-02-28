#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>



int main(int argc, char** argv)
{

  std::cout << bem2d::AmosBesselH0(bem2d::complex(1,2)) << std::endl;
  std::cout << bem2d::AmosBesselH1(bem2d::complex(2,3)) << std::endl;
  return 0;

}

