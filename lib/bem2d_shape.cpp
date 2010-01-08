#include "bem2d_shape.h"
#include "bem2d_defs.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv.h"

namespace bem2d
{

DiskShapePiecewiseConst::DiskShapePiecewiseConst(int n,double radius): elements_(n), points_(n)
{

        double angle=2*PI/n;
        for (int i=0; i<n; i++) points_[i]=boost::shared_ptr<bem2d::Point>(new bem2d::Point(radius*cos(i*angle),radius*sin(i*angle)));
        for (int i=0; i<n-1; i++) {
                elements_[i]=bem2d::pElement(new bem2d::ConstElement(*(points_[i]),*(points_[i+1]),i));
        }
        elements_[n-1]=bem2d::pElement(new bem2d::ConstElement(*(points_[n-1]),*(points_[0]),n-1));
        for (int i=1; i<n-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[0]->set_prev(n-1);
        elements_[n-1]->set_next(0);
        elements_[n-1]->set_prev(n-2);

}

pGeometry DiskShapePiecewiseConst::GetGeometry()
{
        return bem2d::pGeometry(new bem2d::Geometry(elements_));
}

Polygon::Polygon(const std::vector<Point>& points, int n)
{

        int esize=n*points.size();
        elements_.resize(esize);
        std::vector<Point> p(points);
        p.push_back(p[0]); // Add last element again - makes next for-loop easier
        for (int i=0; i<points.size(); i++) {
                Point direction=p[i+1]-p[i];
                for (int j=0; j<n; j++) {
                        Point p1=p[i]+(1.0*j)/n*direction;
                        Point p2=p[i]+(1.0*(j+1))/n*direction;
                        elements_[n*i+j]=bem2d::pElement(new bem2d::ConstElement(p1,p2,n*i+j));
                }
        }
        for (int i=1; i<elements_.size()-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[0]->set_prev(esize-1);
        elements_[n-1]->set_next(0);
        elements_[esize-1]->set_prev(esize-2);

}


pGeometry Polygon::GetGeometry()
{
        return bem2d::pGeometry(new bem2d::Geometry(elements_));
}


int InvAbsDerivative(double t, const double y[], double f[], void* c){
  f[0]=1./AbsDerivative(y[0],c);
  return GSL_SUCCESS;
}

  AnalyticCurve::AnalyticCurve(int n, pCurve curve): elements_(n), curve_(curve)
{
        ParameterizeArc(n);
        double h=1.0/n;
        for (int i=0; i<n-1; i++) {
                elements_[i]=pElement(new AnalyticCurveElement(arclengthparam[i], 
							       arclengthparam[i+1],curve,i));
        }
        elements_[n-1]=pElement(new AnalyticCurveElement(arclengthparam[n-1], 1,curve,n-1));
        for (int i=1; i<n-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[0]->set_prev(n-1);
        elements_[n-1]->set_next(0);
        elements_[n-1]->set_prev(n-2);
}

pGeometry AnalyticCurve::GetGeometry()
{
        return pGeometry(new bem2d::Geometry(elements_));
}



void AnalyticCurve::ParameterizeArc(int n){

  arclengthparam.resize(n+1);
  arclengthparam[0]=0;
  double L=curve_->Length();
  const gsl_odeiv_step_type* T = gsl_odeiv_step_rk2;
  gsl_odeiv_step* s = gsl_odeiv_step_alloc(T,1);
  gsl_odeiv_control* c=gsl_odeiv_control_y_new(1E-10,0);
  gsl_odeiv_evolve* e=gsl_odeiv_evolve_alloc(1);
  gsl_odeiv_system sys={&InvAbsDerivative,NULL,1,&(*curve_)};
  double t0=0.0, t1=L;

  double h=L/n;
  double y=0.0;
  double t=0;

  for (int i=1; i<=n;i++)
    {
      double ti=i*t1/n;
      while (t<ti){
	gsl_odeiv_evolve_apply(e,c,s,
			       &sys,
			       &t, ti,&h,&y);
      }
      arclengthparam[i]=y;
    }
  arclengthparam[n]=1; // Set last value to theoretically correct one.
  
			       
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  
}

}
