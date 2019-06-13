// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// chebSetN
void chebSetN(int n);
RcppExport SEXP _schrodinger_chebSetN(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    chebSetN(n);
    return R_NilValue;
END_RCPP
}
// computeSpectrum
List computeSpectrum(NumericVector px, NumericVector py, int nEigen, std::string method, double dE, double tol);
RcppExport SEXP _schrodinger_computeSpectrum(SEXP pxSEXP, SEXP pySEXP, SEXP nEigenSEXP, SEXP methodSEXP, SEXP dESEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type px(pxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type py(pySEXP);
    Rcpp::traits::input_parameter< int >::type nEigen(nEigenSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type dE(dESEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(computeSpectrum(px, py, nEigen, method, dE, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_schrodinger_chebSetN", (DL_FUNC) &_schrodinger_chebSetN, 1},
    {"_schrodinger_computeSpectrum", (DL_FUNC) &_schrodinger_computeSpectrum, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_schrodinger(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
