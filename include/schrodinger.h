#ifndef _SCHRODINGER_H_
#define _SCHRODINGER_H_

#include "chebspec.h"
#include "numerov.h"

#define NUMEROV 1
#define CHEBYSHEV 2

struct List
{
  vector<double> Es;
  vector< vector< vector<double> > > wfs;
  List(const vector<double> &es, const vector< vector< vector<double> > > &Wfs)
  {
    Es = es;
    wfs = Wfs;
  }
};

SolvSpec* setSchroMethod(std::string method);
vector< vector<double> > getPotential(SolvSpec* n);
vector<double> getEnergies(SolvSpec* n);
vector< vector< vector<double> > > getWavefunctions(SolvSpec* n);
List computeSpectrum(const vector<double> &px , const vector<double> &py,
                     int nEigen = 3, std::string method = "cheb",
                     double dE = 0.1, double tol = 1e-9);
;
#endif