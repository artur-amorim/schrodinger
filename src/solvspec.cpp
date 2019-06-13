#ifndef _SOLVSPEC_CPP_
#define _SOLVSPEC_CPP_

#include <iostream>
#include <fstream>
#include "../include/solvspec.h"

using namespace std;

SolvSpec::SolvSpec(){}

void SolvSpec::setPotential(const vector<double> &XX, const vector<double> &VV){
	if (XX.size() != VV.size())
	{
		cout << "X and V must have the same size" << endl;
		throw "error";
	}
	if(XX.size() > 1) 
	{
    	X = XX;
		PotVals = VV;
  	} 
	else 
	{
		cout << "Give more points to define the potential." << endl;
		throw "error"; 
  }

  nPoints = XX.size() ;
  xMin    = XX.front() ;
  xMax    = XX.back() ;
  h = (xMax - xMin) / nPoints;
  // Now we define the object that will be used
  // to compute the potential at any x
  double yder1 = (VV[0] - VV[1]) / (XX[0] - XX[1]) ;
  double yder2 = (VV[nPoints-1] - VV[nPoints-2]) / (XX[nPoints-1] - XX[nPoints-2]) ;
  potFunc = Spline_Interp<double>(XX, VV, yder1, yder2);

}

vector<Point> SolvSpec::getPotential()
{
  // Returns the potential attribute
	vector<Point> v;
	for(int i = 0; i < nPoints; i++)
	{
		v.push_back(Point(X[i], PotVals[i]));
	}
	return v;
}

void SolvSpec::findSpectrum(int nEigen){}

Spectrum SolvSpec::getSpectrum()
{
	// Returns the spectrum attribute
	return spectrum;
}

void SolvSpec::savePotential() 
{
  ofstream f("potential.dat");
	for (double x = xMin; x < xMax; x += 0.1) 
	{
		f << x << " " << potFunc.interp(x) << endl;
	}
	f.close();
}

SolvSpec::~SolvSpec(){}

;
#endif
