#include "gsl/gsl_integration.h"
#include "bem2d_curve.h"

namespace bem2d {

double AbsDerivative(double t, void* c){
  return length(((Curve*)c)->Deriv(t));
}

  double Curve::Length(){

  gsl_integration_workspace* w=gsl_integration_workspace_alloc(1000);

  double result, error;
  gsl_function F;
  F.function=&AbsDerivative;
  F.params=this; // (curve) is shared pointer but need standard pointer
  gsl_integration_qag(&F,0,1,0,1E-4,1000,GSL_INTEG_GAUSS61,w,&result,&error);
  gsl_integration_workspace_free(w);

  return result;

}

}
