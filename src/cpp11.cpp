// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// dp_walker_update.cpp
list DPWalkerUpdate_cpp(doubles calendar_ages, doubles current_weight, doubles current_v, integers current_cluster_ids, int current_n_clust, double alpha, double mu_phi, double lambda, double nu1, double nu2);
extern "C" SEXP _carbondate_DPWalkerUpdate_cpp(SEXP calendar_ages, SEXP current_weight, SEXP current_v, SEXP current_cluster_ids, SEXP current_n_clust, SEXP alpha, SEXP mu_phi, SEXP lambda, SEXP nu1, SEXP nu2) {
  BEGIN_CPP11
    return cpp11::as_sexp(DPWalkerUpdate_cpp(cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<doubles>>(current_weight), cpp11::as_cpp<cpp11::decay_t<doubles>>(current_v), cpp11::as_cpp<cpp11::decay_t<integers>>(current_cluster_ids), cpp11::as_cpp<cpp11::decay_t<int>>(current_n_clust), cpp11::as_cpp<cpp11::decay_t<double>>(alpha), cpp11::as_cpp<cpp11::decay_t<double>>(mu_phi), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2)));
  END_CPP11
}
// polya_urn_update.cpp
double LogMarginalNormalGamma(double calendar_age, double lambda, double nu1, double nu2, double mu_phi);
extern "C" SEXP _carbondate_LogMarginalNormalGamma(SEXP calendar_age, SEXP lambda, SEXP nu1, SEXP nu2, SEXP mu_phi) {
  BEGIN_CPP11
    return cpp11::as_sexp(LogMarginalNormalGamma(cpp11::as_cpp<cpp11::decay_t<double>>(calendar_age), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2), cpp11::as_cpp<cpp11::decay_t<double>>(mu_phi)));
  END_CPP11
}
// polya_urn_update.cpp
double log_cpp(double val);
extern "C" SEXP _carbondate_log_cpp(SEXP val) {
  BEGIN_CPP11
    return cpp11::as_sexp(log_cpp(cpp11::as_cpp<cpp11::decay_t<double>>(val)));
  END_CPP11
}
// polya_urn_update.cpp
double lgamma_cpp(double val);
extern "C" SEXP _carbondate_lgamma_cpp(SEXP val) {
  BEGIN_CPP11
    return cpp11::as_sexp(lgamma_cpp(cpp11::as_cpp<cpp11::decay_t<double>>(val)));
  END_CPP11
}
// polya_urn_update.cpp
list PolyaUrnUpdateClusterIdentifier(doubles calendar_ages, integers current_cluster_ids, doubles current_phi, doubles current_tau, double alpha, double mu_phi, double lambda, double nu1, double nu2);
extern "C" SEXP _carbondate_PolyaUrnUpdateClusterIdentifier(SEXP calendar_ages, SEXP current_cluster_ids, SEXP current_phi, SEXP current_tau, SEXP alpha, SEXP mu_phi, SEXP lambda, SEXP nu1, SEXP nu2) {
  BEGIN_CPP11
    return cpp11::as_sexp(PolyaUrnUpdateClusterIdentifier(cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<integers>>(current_cluster_ids), cpp11::as_cpp<cpp11::decay_t<doubles>>(current_phi), cpp11::as_cpp<cpp11::decay_t<doubles>>(current_tau), cpp11::as_cpp<cpp11::decay_t<double>>(alpha), cpp11::as_cpp<cpp11::decay_t<double>>(mu_phi), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2)));
  END_CPP11
}
// update_calendar_ages.cpp
doubles UpdateCalendarAges_cpp(int n, doubles calendar_ages, double w, double m, integers cluster_identifiers, doubles phi, doubles tau, doubles c14_determinations, doubles c14_sigmas, doubles mucalallyr, doubles sigcalallyr);
extern "C" SEXP _carbondate_UpdateCalendarAges_cpp(SEXP n, SEXP calendar_ages, SEXP w, SEXP m, SEXP cluster_identifiers, SEXP phi, SEXP tau, SEXP c14_determinations, SEXP c14_sigmas, SEXP mucalallyr, SEXP sigcalallyr) {
  BEGIN_CPP11
    return cpp11::as_sexp(UpdateCalendarAges_cpp(cpp11::as_cpp<cpp11::decay_t<int>>(n), cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<double>>(w), cpp11::as_cpp<cpp11::decay_t<double>>(m), cpp11::as_cpp<cpp11::decay_t<integers>>(cluster_identifiers), cpp11::as_cpp<cpp11::decay_t<doubles>>(phi), cpp11::as_cpp<cpp11::decay_t<doubles>>(tau), cpp11::as_cpp<cpp11::decay_t<doubles>>(c14_determinations), cpp11::as_cpp<cpp11::decay_t<doubles>>(c14_sigmas), cpp11::as_cpp<cpp11::decay_t<doubles>>(mucalallyr), cpp11::as_cpp<cpp11::decay_t<doubles>>(sigcalallyr)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_carbondate_DPWalkerUpdate_cpp",              (DL_FUNC) &_carbondate_DPWalkerUpdate_cpp,              10},
    {"_carbondate_LogMarginalNormalGamma",          (DL_FUNC) &_carbondate_LogMarginalNormalGamma,           5},
    {"_carbondate_PolyaUrnUpdateClusterIdentifier", (DL_FUNC) &_carbondate_PolyaUrnUpdateClusterIdentifier,  9},
    {"_carbondate_UpdateCalendarAges_cpp",          (DL_FUNC) &_carbondate_UpdateCalendarAges_cpp,          11},
    {"_carbondate_lgamma_cpp",                      (DL_FUNC) &_carbondate_lgamma_cpp,                       1},
    {"_carbondate_log_cpp",                         (DL_FUNC) &_carbondate_log_cpp,                          1},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_carbondate(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
