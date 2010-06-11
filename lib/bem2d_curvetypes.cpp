#include "bem2d_curvetypes.h"


namespace bem2d {

  Ellipse::Ellipse(double rho): rho_(rho){}

  EllipseArc::EllipseArc(double a, double b, double theta1, double theta2):
    a_(a), b_(b), theta1_(theta1), theta2_(theta2) {}

  InvEllipse::InvEllipse(double alpha): alpha_(alpha){};

}
