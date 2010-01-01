#include<cmath>
#include<iostream>
#include "bem2d_quadpoints.h"
#include "bem2d_quadrature.h"

namespace bem2d
{

Gauss1D::Gauss1D(int N) throw (ParameterException): N_(N)
{
        if (N<1) throw ParameterException();
        Gauss(x_,w_,N);
}

QuadPoints1D Gauss1D::Rescale(double a, double b) const throw (IntervalException)
{
        if (a>=b) throw IntervalException();
        QuadPoints1D p;
        p.a=a;
        p.b=b;
        p.x=x_;
        p.w=w_;
        MapPoints(p.x,p.w,a,b);
        return p;
}

Gauss1D::~Gauss1D()
{
}

Gauss2D::Gauss2D(int N) throw (ParameterException): N_(N*N)
{
        if (N<1) throw ParameterException();

        dvector xt; // 1d points
        dvector wt; // 1d weights
        Gauss(xt, wt, N);

        x_.resize(N_);
        y_.resize(N_);
        w_.resize(N_);

        for (int i=0; i<N; i++)
                for (int j=0; j<N; j++) {
                        x_[i*N+j]=xt[i];
                        y_[i*N+j]=xt[j];
                        w_[i*N+j]=wt[i]*wt[j];
                }

}

QuadPoints2D Gauss2D::Rescale(double a, double b, double c, double d) const throw (IntervalException)
{
        if (a>=b) throw IntervalException();
        if (c>=d) throw IntervalException();


        QuadPoints2D p;
        p.a=a;
        p.b=b;
        p.c=c;
        p.d=d;
        p.x=x_;
        p.y=y_;
        p.w=w_;

        MapPoints2d(p.x,p.y,p.w,a,b,c,d);
        return p;
}

Gauss2D::~Gauss2D()
{
}

AdaptedGauss1::AdaptedGauss1(int N, int L, double sigma): N_(N*N*L*L), L_(L), sigma_(sigma)
{


        dvector xl;
        for (int i=0; i<L; i++) xl.push_back(pow(sigma,L-i));
        xl.push_back(1.0);

        x_.reserve(N_);
        y_.reserve(N_);
        w_.reserve(N_);

        Gauss2D g2d(N);

        // Now iterate to create the right boxes and insert them


        for (int i=0; i<L; i++)
                for (int j=0; j<L; j++) {
                        QuadPoints2D q=g2d.Rescale(xl[i],xl[i+1],xl[j],xl[j+1]);
                        x_.insert(x_.end(),q.x.begin(),q.x.end());
                        y_.insert(y_.end(),q.y.begin(),q.y.end());
                        w_.insert(w_.end(),q.w.begin(),q.w.end());
                }

}

AdaptedGauss1::~AdaptedGauss1() {}


AdaptedGauss2::AdaptedGauss2(int N, int L, double sigma): N_(3*N*N*L), L_(L), sigma_(sigma)
{


        x_.reserve(N_);
        y_.reserve(N_);
        w_.reserve(N_);

        Gauss2D g2d(N);

        // Now iterate to create the right boxes and insert them

        double a=sigma;
        double b=1;

        while (L>0) {
                QuadPoints2D q1=g2d.Rescale(0,a,1-b,1-a);
                x_.insert(x_.end(),q1.x.begin(),q1.x.end());
                y_.insert(y_.end(),q1.y.begin(),q1.y.end());
                w_.insert(w_.end(),q1.w.begin(),q1.w.end());

                QuadPoints2D q2=g2d.Rescale(a,b,1-b,1-a);
                x_.insert(x_.end(),q2.x.begin(),q2.x.end());
                y_.insert(y_.end(),q2.y.begin(),q2.y.end());
                w_.insert(w_.end(),q2.w.begin(),q2.w.end());

                QuadPoints2D q3=g2d.Rescale(a,b,1-a,1);
                x_.insert(x_.end(),q3.x.begin(),q3.x.end());
                y_.insert(y_.end(),q3.y.begin(),q3.y.end());
                w_.insert(w_.end(),q3.w.begin(),q3.w.end());

                a*=sigma;
                b*=sigma;
                L=L-1;
        }


}

AdaptedGauss2::~AdaptedGauss2() {}


AdaptedGauss3::AdaptedGauss3(int N, int L, double sigma): N_(3*N*N*L), L_(L), sigma_(sigma)
{


        x_.reserve(N_);
        y_.reserve(N_);
        w_.reserve(N_);

        Gauss2D g2d(N);

        // Now iterate to create the right boxes and insert them

        double a=1;
        double b=sigma;

        while (L>0) {
                QuadPoints2D q1=g2d.Rescale(1-a,1-b,0,b);
                x_.insert(x_.end(),q1.x.begin(),q1.x.end());
                y_.insert(y_.end(),q1.y.begin(),q1.y.end());
                w_.insert(w_.end(),q1.w.begin(),q1.w.end());

                QuadPoints2D q2=g2d.Rescale(1-a,1-b,b,a);
                x_.insert(x_.end(),q2.x.begin(),q2.x.end());
                y_.insert(y_.end(),q2.y.begin(),q2.y.end());
                w_.insert(w_.end(),q2.w.begin(),q2.w.end());

                QuadPoints2D q3=g2d.Rescale(1-b,1,b,a);
                x_.insert(x_.end(),q3.x.begin(),q3.x.end());
                y_.insert(y_.end(),q3.y.begin(),q3.y.end());
                w_.insert(w_.end(),q3.w.begin(),q3.w.end());

                a*=sigma;
                b*=sigma;
                L=L-1;
        }

}

AdaptedGauss3::~AdaptedGauss3() {}


}

