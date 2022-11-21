// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// DPWalkerUpdate.cpp
double min_value(doubles vec);
extern "C" SEXP _carbondate_min_value(SEXP vec) {
  BEGIN_CPP11
    return cpp11::as_sexp(min_value(cpp11::as_cpp<cpp11::decay_t<doubles>>(vec)));
  END_CPP11
}
// DPWalkerUpdate.cpp
list DPWalkerUpdate_cpp(doubles calendar_ages, doubles weightprev, doubles vprev, integers cluster_identifiers, doubles phi, doubles tau, int n_clust, double alpha, double mu_phi, double lambda, double nu1, double nu2);
extern "C" SEXP _carbondate_DPWalkerUpdate_cpp(SEXP calendar_ages, SEXP weightprev, SEXP vprev, SEXP cluster_identifiers, SEXP phi, SEXP tau, SEXP n_clust, SEXP alpha, SEXP mu_phi, SEXP lambda, SEXP nu1, SEXP nu2) {
  BEGIN_CPP11
    return cpp11::as_sexp(DPWalkerUpdate_cpp(cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<doubles>>(weightprev), cpp11::as_cpp<cpp11::decay_t<doubles>>(vprev), cpp11::as_cpp<cpp11::decay_t<integers>>(cluster_identifiers), cpp11::as_cpp<cpp11::decay_t<doubles>>(phi), cpp11::as_cpp<cpp11::decay_t<doubles>>(tau), cpp11::as_cpp<cpp11::decay_t<int>>(n_clust), cpp11::as_cpp<cpp11::decay_t<double>>(alpha), cpp11::as_cpp<cpp11::decay_t<double>>(mu_phi), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2)));
  END_CPP11
}
// update_helpers.cpp
integers which_equal(integers vec, int i);
extern "C" SEXP _carbondate_which_equal(SEXP vec, SEXP i) {
  BEGIN_CPP11
    return cpp11::as_sexp(which_equal(cpp11::as_cpp<cpp11::decay_t<integers>>(vec), cpp11::as_cpp<cpp11::decay_t<int>>(i)));
  END_CPP11
}
// update_helpers.cpp
integers which_mt(doubles vec, double i);
extern "C" SEXP _carbondate_which_mt(SEXP vec, SEXP i) {
  BEGIN_CPP11
    return cpp11::as_sexp(which_mt(cpp11::as_cpp<cpp11::decay_t<doubles>>(vec), cpp11::as_cpp<cpp11::decay_t<double>>(i)));
  END_CPP11
}
// update_helpers.cpp
integers which_mt_int(integers vec, int i);
extern "C" SEXP _carbondate_which_mt_int(SEXP vec, SEXP i) {
  BEGIN_CPP11
    return cpp11::as_sexp(which_mt_int(cpp11::as_cpp<cpp11::decay_t<integers>>(vec), cpp11::as_cpp<cpp11::decay_t<int>>(i)));
  END_CPP11
}
// UpdateCalendarAges.cpp
doubles UpdateCalendarAges_cpp(int n, doubles calendar_ages, double w, double m, integers cluster_identifiers, doubles phi, doubles tau, doubles c14_determinations, doubles c14_sigmas, doubles mucalallyr, doubles sigcalallyr);
extern "C" SEXP _carbondate_UpdateCalendarAges_cpp(SEXP n, SEXP calendar_ages, SEXP w, SEXP m, SEXP cluster_identifiers, SEXP phi, SEXP tau, SEXP c14_determinations, SEXP c14_sigmas, SEXP mucalallyr, SEXP sigcalallyr) {
  BEGIN_CPP11
    return cpp11::as_sexp(UpdateCalendarAges_cpp(cpp11::as_cpp<cpp11::decay_t<int>>(n), cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<double>>(w), cpp11::as_cpp<cpp11::decay_t<double>>(m), cpp11::as_cpp<cpp11::decay_t<integers>>(cluster_identifiers), cpp11::as_cpp<cpp11::decay_t<doubles>>(phi), cpp11::as_cpp<cpp11::decay_t<doubles>>(tau), cpp11::as_cpp<cpp11::decay_t<doubles>>(c14_determinations), cpp11::as_cpp<cpp11::decay_t<doubles>>(c14_sigmas), cpp11::as_cpp<cpp11::decay_t<doubles>>(mucalallyr), cpp11::as_cpp<cpp11::decay_t<doubles>>(sigcalallyr)));
  END_CPP11
}
// UpdatePhiTau.cpp
doubles UpdatePhiTau_cpp(doubles calendar_ages, double mu_phi, double lambda, double nu1, double nu2);
extern "C" SEXP _carbondate_UpdatePhiTau_cpp(SEXP calendar_ages, SEXP mu_phi, SEXP lambda, SEXP nu1, SEXP nu2) {
  BEGIN_CPP11
    return cpp11::as_sexp(UpdatePhiTau_cpp(cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<double>>(mu_phi), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_carbondate_DPWalkerUpdate_cpp",     (DL_FUNC) &_carbondate_DPWalkerUpdate_cpp,     12},
    {"_carbondate_UpdateCalendarAges_cpp", (DL_FUNC) &_carbondate_UpdateCalendarAges_cpp, 11},
    {"_carbondate_UpdatePhiTau_cpp",       (DL_FUNC) &_carbondate_UpdatePhiTau_cpp,        5},
    {"_carbondate_min_value",              (DL_FUNC) &_carbondate_min_value,               1},
    {"_carbondate_which_equal",            (DL_FUNC) &_carbondate_which_equal,             2},
    {"_carbondate_which_mt",               (DL_FUNC) &_carbondate_which_mt,                2},
    {"_carbondate_which_mt_int",           (DL_FUNC) &_carbondate_which_mt_int,            2},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_carbondate(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
