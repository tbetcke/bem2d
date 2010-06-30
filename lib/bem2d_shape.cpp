#include<algorithm>
#include<cmath>
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

  Polygon::Polygon(const std::vector<Point>& points, int ppw, freqtype k,int L, double sigma){
        std::vector<Point> p(points);
        p.push_back(p[0]); // Add last element again - makes next for-loop easier
	int elemcount=0;
        for (int i=0; i<points.size(); i++) {
                Point direction=p[i+1]-p[i];
		int n=(int)ceil(ppw*(double)k.re*length(direction)/2.0/PI);
		n=std::max(10,n); // At least 10 elements
		// Add refined version of first element
		if (L==0){
		  Point p1=p[i];
		  Point p2=p[i]+(1.0/n)*direction;
		  elements_.push_back(bem2d::pElement(new bem2d::ConstElement(p1,p2,elemcount)));
		  elemcount++;
		 }
		else{
		  // Add smallest element
		    double tsigma=1;
		    for (int m=0;m<L;m++) tsigma=tsigma*sigma;
		  {
		    Point p1=p[i];
		    Point p2=p[i]+(tsigma/n)*direction;
		    elements_.push_back(bem2d::pElement(new bem2d::ConstElement(p1,p2,elemcount)));
		    elemcount++;
		  }
		  // Add the other elements
		  for (int j=0;j<L;j++){
		    Point p1=p[i]+tsigma/n*direction;
		    tsigma=tsigma/sigma;
		    Point p2=p[i]+tsigma/n*direction;
		    elements_.push_back(bem2d::pElement(new bem2d::ConstElement(p1,p2,elemcount)));
		    elemcount++;
		}
		}
		  // Add the non-refined elements
                for (int j=1; j<n-1; j++) {
                        Point p1=p[i]+(1.0*j)/n*direction;
                        Point p2=p[i]+(1.0*(j+1))/n*direction;
                        elements_.push_back(bem2d::pElement(new bem2d::ConstElement(p1,p2,elemcount)));
			elemcount++;
                }
		// Add refined version of last element
		if (L==0){
		  Point p1=p[i]+(1.0*(n-1))/n*direction;
		  Point p2=p[i+1];
		  elements_.push_back(bem2d::pElement(new bem2d::ConstElement(p1,p2,elemcount)));
		  elemcount++;
		}
		else{
		  double tsigma=1;
		  // Add exponentially weighted elements
		  for (int j=0;j<L;j++){
		    Point p1=p[i]+(1.0*(n-1))/n*direction+(1-tsigma)/n*direction;
		    tsigma*=sigma;
		    Point p2=p[i]+(1.0*(n-1))/n*direction+(1-tsigma)/n*direction;
		    elements_.push_back(bem2d::pElement(new bem2d::ConstElement(p1,p2,elemcount)));
		    elemcount++;
		}
		  // Add largest element
		  {
		    Point p1=p[i]+(1.0*(n-1))/n*direction+(1-tsigma)/n*direction;
		    Point p2=p[i+1];
		    elements_.push_back(bem2d::pElement(new bem2d::ConstElement(p1,p2,elemcount)));
		    elemcount++;
		  }
		}	       
        }
        for (int i=1; i<elements_.size()-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[0]->set_prev(elemcount-1);
        elements_[elemcount-1]->set_next(0);
        elements_[elemcount-1]->set_prev(elemcount-2);

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
        elements_[esize-1]->set_next(0);
        elements_[esize-1]->set_prev(esize-2);

}


pGeometry Polygon::GetGeometry()
{
        return bem2d::pGeometry(new bem2d::Geometry(elements_));
}


Line::Line(Point p0, Point p1, int ppw, freqtype k, int L, double sigma)
{
    int elemcount=0;
            Point direction=p1-p0;
            int n=(int)ceil(ppw*(double)k.re*length(direction)/2.0/PI);
            n=std::max(10,n); // At least 10 elements
            // Add refined version of first element
            if (L==0){
              elements_.push_back(bem2d::pElement(new bem2d::ConstElement(p0,p0+1./n*direction,elemcount)));
              elemcount++;
            }
            else{
              // Add smallest element
                double tsigma=1;
                for (int m=0;m<L;m++) tsigma=tsigma*sigma;
              {
                Point pp1=p0;
                Point pp2=p0+(tsigma/n)*direction;
                elements_.push_back(bem2d::pElement(new bem2d::ConstElement(pp1,pp2,elemcount)));
                elemcount++;
              }
              // Add the other elements
              for (int j=0;j<L;j++){
                Point pp1=p0+tsigma/n*direction;
                tsigma=tsigma/sigma;
                Point pp2=p0+tsigma/n*direction;
                elements_.push_back(bem2d::pElement(new bem2d::ConstElement(pp1,pp2,elemcount)));
                elemcount++;
            }
            }
              // Add the non-refined elements
            for (int j=1; j<n-1; j++) {
                    Point pp1=p0+(1.0*j)/n*direction;
                    Point pp2=p0+(1.0*(j+1))/n*direction;
                    elements_.push_back(bem2d::pElement(new bem2d::ConstElement(pp1,pp2,elemcount)));
                    elemcount++;
            }
            // Add refined version of last element
            if (L==0){
              Point pp1=p0+(1.0*(n-1))/n*direction;
              Point pp2=p1;
              elements_.push_back(bem2d::pElement(new bem2d::ConstElement(pp1,pp2,elemcount)));
              elemcount++;
            }
            else{
              double tsigma=1;
              // Add exponentially weighted elements
              for (int j=0;j<L;j++){
                Point pp1=p0+(1.0*(n-1))/n*direction+(1-tsigma)/n*direction;
                tsigma*=sigma;
                Point pp2=p0+(1.0*(n-1))/n*direction+(1-tsigma)/n*direction;
                elements_.push_back(bem2d::pElement(new bem2d::ConstElement(pp1,pp2,elemcount)));
                elemcount++;
            }
              // Add smallest element
              {
                Point pp1=p0+(1.0*(n-1))/n*direction+(1-tsigma)/n*direction;
                Point pp2=p1;
                elements_.push_back(bem2d::pElement(new bem2d::ConstElement(pp1,pp2,elemcount)));
                elemcount++;
              }
            }
        for (int i=1; i<elements_.size()-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[0]->set_prev(elemcount-1);
        elements_[elemcount-1]->set_next(0);
        elements_[elemcount-1]->set_prev(elemcount-2);


}
 

Line::Line(Point p0, Point p1, int n){
            Point direction=p1-p0;
            for (int j=0; j<n; j++) {
                        Point pp1=p0+(1.0*j)/n*direction;
                        Point pp2=p0+(1.0*(j+1))/n*direction;
                        elements_.push_back(bem2d::pElement(new bem2d::ConstElement(pp1,pp2,j)));
                }
        
        for (int i=1; i<elements_.size()-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[n-1]->set_prev(n-2);
}

pGeometry Line::GetGeometry(){

        return bem2d::pGeometry(new bem2d::Geometry(elements_));
}

int InvAbsDerivative(double t, const double y[], double f[], void* c){
  f[0]=1./AbsDerivative(y[0],c);
  return GSL_SUCCESS;
}

  AnalyticCurve::AnalyticCurve(int n, pCurve curve,int closed): elements_(n), curve_(curve)
{

		double curveLength=curve_->Length();
		dvector partition(n+1);
	for (int i=0;i<=n;i++) partition[i]=i*curveLength/n;
        ParameterizeArc(partition);
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
        elements_[n-1]->set_prev(n-2);
        if (closed) {
            elements_[0]->set_prev(n-1);
            elements_[n-1]->set_next(0);
        }
}

	AnalyticCurve::AnalyticCurve(int ppw, freqtype k, pCurve curve, int closed, int L, double sigma): curve_(curve)
	{

		double curveLength=curve_->Length();
		dvector partition;
		partition.push_back(0);
		int c=(int)ceil(ppw*(double)k.re*curveLength/2.0/PI);
		c=std::max(c,10);
		if (L==0){ 
			partition.push_back(curveLength/c);
		}
		else {
			double tsigma=1;
			for (int m=0;m<L;m++) tsigma=tsigma*sigma;
			partition.push_back(tsigma*curveLength/c);
			for (int m=0;m<L;m++){
				tsigma=tsigma/sigma;
				partition.push_back(tsigma*curveLength/c);
			}
		}
		for (int j=2;j<c;j++) partition.push_back(curveLength/c*j);
		if (L==0){
			partition.push_back(curveLength);
		}
		else {
			double td=curveLength/c*(c-1);
			double tsigma=1;
			for (int m=0;m<L;m++){
				tsigma=tsigma*sigma;
				partition.push_back(td+curveLength/c*(1-tsigma));
			}
			partition.push_back(curveLength);
		}
			
        ParameterizeArc(partition);
		int n=partition.size()-1;
		elements_.resize(n);
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
        elements_[n-1]->set_prev(n-2);
        if (closed) {
            elements_[0]->set_prev(n-1);
            elements_[n-1]->set_next(0);
        }
	}
	

pGeometry AnalyticCurve::GetGeometry()
{
        return pGeometry(new bem2d::Geometry(elements_));
}



void AnalyticCurve::ParameterizeArc(dvector& partition){

  int n=partition.size()-1;
  arclengthparam.resize(n+1);
  arclengthparam[0]=0;
  double L=curve_->Length();
  const gsl_odeiv_step_type* T = gsl_odeiv_step_rk2;
  gsl_odeiv_step* s = gsl_odeiv_step_alloc(T,1);
  gsl_odeiv_control* c=gsl_odeiv_control_y_new(1E-10,0);
  gsl_odeiv_evolve* e=gsl_odeiv_evolve_alloc(1);
  gsl_odeiv_system sys={&InvAbsDerivative,NULL,1,&(*curve_)};
  double t0=0.0, t1=L;

  double y=0.0;
  double t=0;

  for (int i=1; i<=n;i++)
    {
      double ti=partition[i];
	  double h=partition[i]-partition[i-1];
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
