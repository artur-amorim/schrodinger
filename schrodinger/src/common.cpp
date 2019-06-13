#ifndef _COMMON_CPP_
#define _COMMON_CPP_

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

#include "../include/common.h"

using namespace std;

template <typename T> 
int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}

// Van Wijngaarden–Dekker–Brent Method for finding root, from NR
double zbrent(function<double(double)>& func, double x1, double x2, double tol, bool silent)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;

	if (!silent && ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)))
		cout << "?";
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c  = a;			//Rename a, b, c and adjust bounding interval d.
			fc = fa;
			e  = d = b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a  = b;
			b  = c;
			c  = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol1 = 2.0*EPS*fabs(b)+0.5*tol; //Convergence check.
		xm = 0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;		//Attempt inverse quadratic interpolation.
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0-s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;		//Check whether in bounds.
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				//Accept interpolation.
				d = p/q;
			} else {
				d = xm;
				//Interpolation failed, use bisection.
				e = d;
			}
		} else {		//Bounds decreasing too slowly, use bisection.
			d = xm;
			e = d;
		}
		a = b;			//Move last best guess to a.
		fa = fb;
		if (fabs(d) > tol1)	//Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);
		fb = func(b);
	}
	cout << "[WARN] Maximum number of iterations exceeded in zbrent, returning biggest value"<< endl;
	return x2;
}

double bisection(function<double(double)> diffFunc, double min, double max){
	double E;
	int i=0;			     /* counter for iterations */
	do{
		i++;
		E=(max+min)/2.0;		     /* divide energy range */
		if (diffFunc(max)*diffFunc(E)>0)
			max=E;   /* the bisection algorithm */
		else
			min=E;
	}while(fabs(diffFunc(E))>EPS);
	return E;
}

#endif
