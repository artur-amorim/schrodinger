#ifndef _CHEBSPEC_H
#define _CHEBSPEC_H

#include "solvspec.h"

class ChebSpec : public SolvSpec{
    private:
        double a, b, scal;
        vector<double> V;
        vector< vector<double> > Ehat, UE, US;
        void invMat(double*, int);
        void showMatrix(double* A, int Nr);
    public:
        static int N, L;
        static vector<double> x;
        static vector< vector<double> > Tm, TmInt, D, D2;
        ChebSpec();
        void findSpectrum(int nEigen);
        void setPotential(const vector<double> &XX, const vector<double> &VV);
        ~ChebSpec();
};

;
#endif