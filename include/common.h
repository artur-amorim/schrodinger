#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <cmath>
#include <functional>

using namespace std;

template <typename T> 
int sgn(T val);

struct Point {
    double x = 0;
    double y = 0;
    Point() {
        x = 0;
        y = 0;
    }
    Point(double xx, double yy) {
        x = xx;
        y = yy;
    }
};

struct Range {
	double eMin, eMax;
	Range(double m0, double m1) {
	    eMin = m0;
	    eMax = m1;
	}
};

struct Mode {
	double energy;
  int index;
	vector<Point> wavefunction;
	Mode() {   // default constructor for non good modes
    index = -1;
	}
	Mode(double e, vector<Point> f){
		energy = e;
		wavefunction = f;
	}
	Mode(double e, vector<Point> f, int n){
		energy = e;
		wavefunction = f;
		index = n;
	}
};

struct Spectrum {
	vector<Mode>  modes;
	vector<vector<double> > potential;

	void addMode(Mode m) {
   	modes.push_back(m);
	}

	void clear() 
	{
   		modes.clear();
	  	potential.clear();
 	}

  vector<double> getEnergies()
  {
    vector<double> energies;
    for(int i = 0; i < modes.size(); i++)
	{
      energies.push_back(modes[i].energy);
	}
    return energies;
  }

  vector<vector<Point> > getWavefunctions()
  {
    vector<vector<Point> > wfs;
    for(int i = 0; i < modes.size(); i++)
	{
      wfs.push_back(modes[i].wavefunction);
	}
    return wfs;
  }
};

// Van Wijngaarden–Dekker–Brent Method for finding root, from NR
#define ITMAX 600
#define EPS 1e-9
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double zbrent(function<double(double)>& func, double x1, double x2, double tol, bool silent = false);

// bisection method to find root
double bisection(function<double(double)> diffFunc, double min, double max);

;
#endif