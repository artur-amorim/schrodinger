#ifndef _SOLVSPEC_H
#define _SOLVSPEC_H

#include "../src/common.hpp"
#include "../src/interp_1d.hpp"

class SolvSpec {
  protected:
  	vector<double> X, PotVals;
  	Spectrum spectrum;
  public:
  	double xMin, xMax, tol;
    // Attributes relevant for the Numerov Method
  	double h, dEmin;
    int nPoints;
    // Object with interp method that computes the potential at arbitrary x
    Spline_Interp<double> potFunc;
    // Class constructor
  	SolvSpec();
    // Setter of potential
  	virtual void setPotential(const vector<double> &XX, const vector<double> &VV);
    // Getter of potential
  	vector<Point> getPotential();
    // Computes the spectrum using a given method
  	virtual void findSpectrum(int nEigen);
    // Get the spectrum
  	Spectrum getSpectrum();
    // Save the potential in a file
  	void savePotential();
    // Class destructor
	virtual ~SolvSpec();
};

;
#endif