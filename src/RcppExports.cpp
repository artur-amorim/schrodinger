// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// setMethod
void setMethod(std::string method);
RcppExport SEXP schrodinger_setMethod(SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    setMethod(method);
    return R_NilValue;
END_RCPP
}
// getEnergiesAndIndices
List getEnergiesAndIndices();
RcppExport SEXP schrodinger_getEnergiesAndIndices() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(getEnergiesAndIndices());
    return __result;
END_RCPP
}
// setPotential
void setPotential(NumericVector px, NumericVector py);
RcppExport SEXP schrodinger_setPotential(SEXP pxSEXP, SEXP pySEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type px(pxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type py(pySEXP);
    setPotential(px, py);
    return R_NilValue;
END_RCPP
}
// getPotential
List getPotential();
RcppExport SEXP schrodinger_getPotential() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(getPotential());
    return __result;
END_RCPP
}
// computeSpectrum
void computeSpectrum(int nEigen, double dE, double tol, int N);
RcppExport SEXP schrodinger_computeSpectrum(SEXP nEigenSEXP, SEXP dESEXP, SEXP tolSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type nEigen(nEigenSEXP);
    Rcpp::traits::input_parameter< double >::type dE(dESEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    computeSpectrum(nEigen, dE, tol, N);
    return R_NilValue;
END_RCPP
}
// getEnergies
NumericVector getEnergies();
RcppExport SEXP schrodinger_getEnergies() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(getEnergies());
    return __result;
END_RCPP
}
// getWavefunctions
List getWavefunctions();
RcppExport SEXP schrodinger_getWavefunctions() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(getWavefunctions());
    return __result;
END_RCPP
}