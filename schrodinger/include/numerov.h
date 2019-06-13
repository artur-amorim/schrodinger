#ifndef _NUMEROV_H
#define _NUMEROV_H

#include "solvspec.h"
#include "common.h"

class Numerov : public SolvSpec{
    private:
        vector<Range> zeros;
        vector<Point> solLR, solRL, sol;
        Point minPot;
        void findMinPot();
        void buildSol();
        double diff(double E);
        void scanForZeroRanges(int nZeros);
        int getSolIndex();
        function<double(double)> diffFunc = [this](double E){ return diff(E);};
    public:
        Numerov();
        void findSpectrum(int nEigen);
        ~Numerov();
};
;
#endif