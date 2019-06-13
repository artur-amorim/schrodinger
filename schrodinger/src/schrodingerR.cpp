#include <string>
#include "../include/schrodinger.h"
#include "chebspec.cpp"
#include "numerov.cpp"
#include "common.cpp"


SolvSpec* setSchroMethod(std::string method) 
{
  if(method == "numerov") {
    return new Numerov();
  }
  // use chebyshev as default method
  return new ChebSpec();
}

vector< vector<double> > getPotential(SolvSpec* n) 
{
  vector<Point> p = n->getPotential();
  int length = p.size();

  vector<double> x, y;
  for(int i = 0; i < length; i++) 
  {
    x.push_back(p[i].x);
    y.push_back(p[i].y);
  }

  vector< vector<double> > pot{x, y};
  return pot;
}

vector<double> getEnergies(SolvSpec* n) 
{
  vector<double> energies = n->getSpectrum().getEnergies();
  return energies;
}

vector< vector< vector<double> > > getWavefunctions(SolvSpec* n) 
{
  vector< vector< vector<double> > > wfs;
  vector<vector<Point> > WFs = n->getSpectrum().getWavefunctions();

  for(int i = 0; i < WFs.size(); i++) 
  {
    int length = WFs[i].size();
    vector<double> x, y;
    for(int j = 0; j < length; j++) 
    {
      x.push_back(WFs[i][j].x);
      y.push_back(WFs[i][j].y);
    }
    vector< vector<double> > wf{x, y};
    wfs.push_back(wf);
  }
  return wfs;
}

List computeSpectrum(const vector<double> &px , const vector<double> &py,
                     int nEigen, std::string method,
                     double dE, double tol) 
{
  
  SolvSpec* n = setSchroMethod(method);
  
  if(px.size() != py.size()) 
  {
    cout << "Please pass two columns with the same size for the potential" << endl;
    throw "error";
  }

  n->setPotential(px, py);
  n->dEmin = dE;
  n->tol = tol;
  n->findSpectrum(nEigen);
  vector<double> energies = getEnergies(n) ;
  vector<vector<vector<double> > > wavefuncs = getWavefunctions(n) ;
  // Free memory
  delete n ;
  return List(energies, wavefuncs) ;
}