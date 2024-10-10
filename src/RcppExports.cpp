// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mvndrawC
arma::colvec mvndrawC(arma::colvec mu, arma::mat sig);
RcppExport SEXP _MFVART_mvndrawC(SEXP muSEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(mvndrawC(mu, sig));
    return rcpp_result_gen;
END_RCPP
}
// carterkohn
List carterkohn(arma::mat y, arma::mat Z, arma::mat Ht, arma::mat Qt, double m, double p, double t, arma::colvec B0, arma::mat V0);
RcppExport SEXP _MFVART_carterkohn(SEXP ySEXP, SEXP ZSEXP, SEXP HtSEXP, SEXP QtSEXP, SEXP mSEXP, SEXP pSEXP, SEXP tSEXP, SEXP B0SEXP, SEXP V0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ht(HtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Qt(QtSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V0(V0SEXP);
    rcpp_result_gen = Rcpp::wrap(carterkohn(y, Z, Ht, Qt, m, p, t, B0, V0));
    return rcpp_result_gen;
END_RCPP
}
// alphahelper
arma::mat alphahelper(arma::mat y, arma::mat Z, arma::mat Btdraw);
RcppExport SEXP _MFVART_alphahelper(SEXP ySEXP, SEXP ZSEXP, SEXP BtdrawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Btdraw(BtdrawSEXP);
    rcpp_result_gen = Rcpp::wrap(alphahelper(y, Z, Btdraw));
    return rcpp_result_gen;
END_RCPP
}
// sigmahelper1
arma::mat sigmahelper1(arma::mat Atdraw, double M);
RcppExport SEXP _MFVART_sigmahelper1(SEXP AtdrawSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Atdraw(AtdrawSEXP);
    Rcpp::traits::input_parameter< double >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmahelper1(Atdraw, M));
    return rcpp_result_gen;
END_RCPP
}
// sigmahelper2
List sigmahelper2(arma::mat capAt, arma::mat yhat, arma::colvec qs, arma::colvec ms, arma::colvec u2s, arma::mat Sigtdraw, arma::mat Zs, arma::mat Wdraw, arma::colvec sigma_prmean, arma::mat sigma_prvar);
RcppExport SEXP _MFVART_sigmahelper2(SEXP capAtSEXP, SEXP yhatSEXP, SEXP qsSEXP, SEXP msSEXP, SEXP u2sSEXP, SEXP SigtdrawSEXP, SEXP ZsSEXP, SEXP WdrawSEXP, SEXP sigma_prmeanSEXP, SEXP sigma_prvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type capAt(capAtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type yhat(yhatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type qs(qsSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type ms(msSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type u2s(u2sSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigtdraw(SigtdrawSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zs(ZsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Wdraw(WdrawSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type sigma_prmean(sigma_prmeanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prvar(sigma_prvarSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmahelper2(capAt, yhat, qs, ms, u2s, Sigtdraw, Zs, Wdraw, sigma_prmean, sigma_prvar));
    return rcpp_result_gen;
END_RCPP
}
// sigmahelper4
List sigmahelper4(arma::mat y2, arma::colvec qs, arma::colvec ms, arma::colvec u2s, arma::mat Sigtdraw, arma::mat Zs, arma::mat Wdraw, arma::colvec sigma_prmean, arma::mat sigma_prvar);
RcppExport SEXP _MFVART_sigmahelper4(SEXP y2SEXP, SEXP qsSEXP, SEXP msSEXP, SEXP u2sSEXP, SEXP SigtdrawSEXP, SEXP ZsSEXP, SEXP WdrawSEXP, SEXP sigma_prmeanSEXP, SEXP sigma_prvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type qs(qsSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type ms(msSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type u2s(u2sSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigtdraw(SigtdrawSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zs(ZsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Wdraw(WdrawSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type sigma_prmean(sigma_prmeanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prvar(sigma_prvarSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmahelper4(y2, qs, ms, u2s, Sigtdraw, Zs, Wdraw, sigma_prmean, sigma_prvar));
    return rcpp_result_gen;
END_RCPP
}
// sigmahelper3
List sigmahelper3(arma::mat capAt, arma::mat sigt);
RcppExport SEXP _MFVART_sigmahelper3(SEXP capAtSEXP, SEXP sigtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type capAt(capAtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigt(sigtSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmahelper3(capAt, sigt));
    return rcpp_result_gen;
END_RCPP
}
// getvc
List getvc(arma::mat Ht);
RcppExport SEXP _MFVART_getvc(SEXP HtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Ht(HtSEXP);
    rcpp_result_gen = Rcpp::wrap(getvc(Ht));
    return rcpp_result_gen;
END_RCPP
}
// makeregs_fcC
arma::mat makeregs_fcC(arma::mat ydat, double p);
RcppExport SEXP _MFVART_makeregs_fcC(SEXP ydatSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ydat(ydatSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(makeregs_fcC(ydat, p));
    return rcpp_result_gen;
END_RCPP
}
// mz
arma::colvec mz(double n);
RcppExport SEXP _MFVART_mz(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mz(n));
    return rcpp_result_gen;
END_RCPP
}
// meye
arma::mat meye(double n);
RcppExport SEXP _MFVART_meye(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(meye(n));
    return rcpp_result_gen;
END_RCPP
}
// vechC
arma::colvec vechC(arma::mat x);
RcppExport SEXP _MFVART_vechC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vechC(x));
    return rcpp_result_gen;
END_RCPP
}
// getfcsts
List getfcsts(arma::colvec Bt0, arma::colvec At0, arma::colvec Sigt0, arma::mat Qdraw, arma::mat Sdraw, arma::mat Wdraw, arma::mat ydat, double nf, double p);
RcppExport SEXP _MFVART_getfcsts(SEXP Bt0SEXP, SEXP At0SEXP, SEXP Sigt0SEXP, SEXP QdrawSEXP, SEXP SdrawSEXP, SEXP WdrawSEXP, SEXP ydatSEXP, SEXP nfSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type Bt0(Bt0SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type At0(At0SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Sigt0(Sigt0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Qdraw(QdrawSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sdraw(SdrawSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Wdraw(WdrawSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ydat(ydatSEXP);
    Rcpp::traits::input_parameter< double >::type nf(nfSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getfcsts(Bt0, At0, Sigt0, Qdraw, Sdraw, Wdraw, ydat, nf, p));
    return rcpp_result_gen;
END_RCPP
}
// wishdrawC
arma::mat wishdrawC(arma::mat h, double n);
RcppExport SEXP _MFVART_wishdrawC(SEXP hSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(wishdrawC(h, n));
    return rcpp_result_gen;
END_RCPP
}
// makeregs2_fcC
arma::rowvec makeregs2_fcC(arma::mat dat, double p);
RcppExport SEXP _MFVART_makeregs2_fcC(SEXP datSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dat(datSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(makeregs2_fcC(dat, p));
    return rcpp_result_gen;
END_RCPP
}
// matmult
arma::mat matmult(arma::mat x, double nt);
RcppExport SEXP _MFVART_matmult(SEXP xSEXP, SEXP ntSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type nt(ntSEXP);
    rcpp_result_gen = Rcpp::wrap(matmult(x, nt));
    return rcpp_result_gen;
END_RCPP
}
// varfcst
List varfcst(arma::mat b, arma::mat sig, arma::mat y, double nf);
RcppExport SEXP _MFVART_varfcst(SEXP bSEXP, SEXP sigSEXP, SEXP ySEXP, SEXP nfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type nf(nfSEXP);
    rcpp_result_gen = Rcpp::wrap(varfcst(b, sig, y, nf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MFVART_mvndrawC", (DL_FUNC) &_MFVART_mvndrawC, 2},
    {"_MFVART_carterkohn", (DL_FUNC) &_MFVART_carterkohn, 9},
    {"_MFVART_alphahelper", (DL_FUNC) &_MFVART_alphahelper, 3},
    {"_MFVART_sigmahelper1", (DL_FUNC) &_MFVART_sigmahelper1, 2},
    {"_MFVART_sigmahelper2", (DL_FUNC) &_MFVART_sigmahelper2, 10},
    {"_MFVART_sigmahelper4", (DL_FUNC) &_MFVART_sigmahelper4, 9},
    {"_MFVART_sigmahelper3", (DL_FUNC) &_MFVART_sigmahelper3, 2},
    {"_MFVART_getvc", (DL_FUNC) &_MFVART_getvc, 1},
    {"_MFVART_makeregs_fcC", (DL_FUNC) &_MFVART_makeregs_fcC, 2},
    {"_MFVART_mz", (DL_FUNC) &_MFVART_mz, 1},
    {"_MFVART_meye", (DL_FUNC) &_MFVART_meye, 1},
    {"_MFVART_vechC", (DL_FUNC) &_MFVART_vechC, 1},
    {"_MFVART_getfcsts", (DL_FUNC) &_MFVART_getfcsts, 9},
    {"_MFVART_wishdrawC", (DL_FUNC) &_MFVART_wishdrawC, 2},
    {"_MFVART_makeregs2_fcC", (DL_FUNC) &_MFVART_makeregs2_fcC, 2},
    {"_MFVART_matmult", (DL_FUNC) &_MFVART_matmult, 2},
    {"_MFVART_varfcst", (DL_FUNC) &_MFVART_varfcst, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_MFVART(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
