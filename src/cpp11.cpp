// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// find_predictive_density.cpp
data_frame FindPredictiveDensityAndCIWalker(doubles calendar_ages, list weights, list phis, list taus, doubles mu_phis, double lambda, double nu1, double nu2, int n_posterior_samples, double quantile_edge_width);
extern "C" SEXP _carbondate_FindPredictiveDensityAndCIWalker(SEXP calendar_ages, SEXP weights, SEXP phis, SEXP taus, SEXP mu_phis, SEXP lambda, SEXP nu1, SEXP nu2, SEXP n_posterior_samples, SEXP quantile_edge_width) {
  BEGIN_CPP11
    return cpp11::as_sexp(FindPredictiveDensityAndCIWalker(cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<list>>(weights), cpp11::as_cpp<cpp11::decay_t<list>>(phis), cpp11::as_cpp<cpp11::decay_t<list>>(taus), cpp11::as_cpp<cpp11::decay_t<doubles>>(mu_phis), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2), cpp11::as_cpp<cpp11::decay_t<int>>(n_posterior_samples), cpp11::as_cpp<cpp11::decay_t<double>>(quantile_edge_width)));
  END_CPP11
}
// find_predictive_density.cpp
data_frame FindPredictiveDensityandCIPolyaUrn(doubles calendar_ages, list observations_per_clusters, list phis, list taus, doubles alphas, doubles mu_phis, double n_obs, double lambda, double nu1, double nu2, int n_posterior_samples, double quantile_edge_width);
extern "C" SEXP _carbondate_FindPredictiveDensityandCIPolyaUrn(SEXP calendar_ages, SEXP observations_per_clusters, SEXP phis, SEXP taus, SEXP alphas, SEXP mu_phis, SEXP n_obs, SEXP lambda, SEXP nu1, SEXP nu2, SEXP n_posterior_samples, SEXP quantile_edge_width) {
  BEGIN_CPP11
    return cpp11::as_sexp(FindPredictiveDensityandCIPolyaUrn(cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<list>>(observations_per_clusters), cpp11::as_cpp<cpp11::decay_t<list>>(phis), cpp11::as_cpp<cpp11::decay_t<list>>(taus), cpp11::as_cpp<cpp11::decay_t<doubles>>(alphas), cpp11::as_cpp<cpp11::decay_t<doubles>>(mu_phis), cpp11::as_cpp<cpp11::decay_t<double>>(n_obs), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2), cpp11::as_cpp<cpp11::decay_t<int>>(n_posterior_samples), cpp11::as_cpp<cpp11::decay_t<double>>(quantile_edge_width)));
  END_CPP11
}
// find_predictive_density.cpp
doubles FindInstantPredictiveDensityWalker(doubles calendar_ages, doubles weight, doubles phi, doubles tau, double mu_phi, double lambda, double nu1, double nu2);
extern "C" SEXP _carbondate_FindInstantPredictiveDensityWalker(SEXP calendar_ages, SEXP weight, SEXP phi, SEXP tau, SEXP mu_phi, SEXP lambda, SEXP nu1, SEXP nu2) {
  BEGIN_CPP11
    return cpp11::as_sexp(FindInstantPredictiveDensityWalker(cpp11::as_cpp<cpp11::decay_t<doubles>>(calendar_ages), cpp11::as_cpp<cpp11::decay_t<doubles>>(weight), cpp11::as_cpp<cpp11::decay_t<doubles>>(phi), cpp11::as_cpp<cpp11::decay_t<doubles>>(tau), cpp11::as_cpp<cpp11::decay_t<double>>(mu_phi), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2)));
  END_CPP11
}
// polya_urn_update_helpers.cpp
double PolyaUrnUpdateAlpha_test(int n, integers nci, double current_alpha, double alpha_shape, double alpha_rate);
extern "C" SEXP _carbondate_PolyaUrnUpdateAlpha_test(SEXP n, SEXP nci, SEXP current_alpha, SEXP alpha_shape, SEXP alpha_rate) {
  BEGIN_CPP11
    return cpp11::as_sexp(PolyaUrnUpdateAlpha_test(cpp11::as_cpp<cpp11::decay_t<int>>(n), cpp11::as_cpp<cpp11::decay_t<integers>>(nci), cpp11::as_cpp<cpp11::decay_t<double>>(current_alpha), cpp11::as_cpp<cpp11::decay_t<double>>(alpha_shape), cpp11::as_cpp<cpp11::decay_t<double>>(alpha_rate)));
  END_CPP11
}
// polya_urn_update_step.cpp
list PolyaUrnUpdateStep(doubles current_calendar_ages, integers current_cluster_ids, doubles current_phi, doubles current_tau, double current_alpha, double current_mu_phi, double alpha_shape, double alpha_rate, double lambda, double nu1, double nu2, double A, double B, double w, double m, doubles c14_determinations, doubles c14_sigmas, int calcurve_yr_index_offset, doubles mucalallyr, doubles sigcalallyr);
extern "C" SEXP _carbondate_PolyaUrnUpdateStep(SEXP current_calendar_ages, SEXP current_cluster_ids, SEXP current_phi, SEXP current_tau, SEXP current_alpha, SEXP current_mu_phi, SEXP alpha_shape, SEXP alpha_rate, SEXP lambda, SEXP nu1, SEXP nu2, SEXP A, SEXP B, SEXP w, SEXP m, SEXP c14_determinations, SEXP c14_sigmas, SEXP calcurve_yr_index_offset, SEXP mucalallyr, SEXP sigcalallyr) {
  BEGIN_CPP11
    return cpp11::as_sexp(PolyaUrnUpdateStep(cpp11::as_cpp<cpp11::decay_t<doubles>>(current_calendar_ages), cpp11::as_cpp<cpp11::decay_t<integers>>(current_cluster_ids), cpp11::as_cpp<cpp11::decay_t<doubles>>(current_phi), cpp11::as_cpp<cpp11::decay_t<doubles>>(current_tau), cpp11::as_cpp<cpp11::decay_t<double>>(current_alpha), cpp11::as_cpp<cpp11::decay_t<double>>(current_mu_phi), cpp11::as_cpp<cpp11::decay_t<double>>(alpha_shape), cpp11::as_cpp<cpp11::decay_t<double>>(alpha_rate), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2), cpp11::as_cpp<cpp11::decay_t<double>>(A), cpp11::as_cpp<cpp11::decay_t<double>>(B), cpp11::as_cpp<cpp11::decay_t<double>>(w), cpp11::as_cpp<cpp11::decay_t<double>>(m), cpp11::as_cpp<cpp11::decay_t<doubles>>(c14_determinations), cpp11::as_cpp<cpp11::decay_t<doubles>>(c14_sigmas), cpp11::as_cpp<cpp11::decay_t<int>>(calcurve_yr_index_offset), cpp11::as_cpp<cpp11::decay_t<doubles>>(mucalallyr), cpp11::as_cpp<cpp11::decay_t<doubles>>(sigcalallyr)));
  END_CPP11
}
// walker_update_step.cpp
list WalkerUpdateStep(doubles current_calendar_ages, doubles current_weight, doubles current_v, integers current_cluster_ids, double current_alpha, double current_mu_phi, double alpha_shape, double alpha_rate, double lambda, double nu1, double nu2, double A, double B, double w, double m, doubles c14_determinations, doubles c14_sigmas, int calcurve_yr_index_offset, doubles mucalallyr, doubles sigcalallyr);
extern "C" SEXP _carbondate_WalkerUpdateStep(SEXP current_calendar_ages, SEXP current_weight, SEXP current_v, SEXP current_cluster_ids, SEXP current_alpha, SEXP current_mu_phi, SEXP alpha_shape, SEXP alpha_rate, SEXP lambda, SEXP nu1, SEXP nu2, SEXP A, SEXP B, SEXP w, SEXP m, SEXP c14_determinations, SEXP c14_sigmas, SEXP calcurve_yr_index_offset, SEXP mucalallyr, SEXP sigcalallyr) {
  BEGIN_CPP11
    return cpp11::as_sexp(WalkerUpdateStep(cpp11::as_cpp<cpp11::decay_t<doubles>>(current_calendar_ages), cpp11::as_cpp<cpp11::decay_t<doubles>>(current_weight), cpp11::as_cpp<cpp11::decay_t<doubles>>(current_v), cpp11::as_cpp<cpp11::decay_t<integers>>(current_cluster_ids), cpp11::as_cpp<cpp11::decay_t<double>>(current_alpha), cpp11::as_cpp<cpp11::decay_t<double>>(current_mu_phi), cpp11::as_cpp<cpp11::decay_t<double>>(alpha_shape), cpp11::as_cpp<cpp11::decay_t<double>>(alpha_rate), cpp11::as_cpp<cpp11::decay_t<double>>(lambda), cpp11::as_cpp<cpp11::decay_t<double>>(nu1), cpp11::as_cpp<cpp11::decay_t<double>>(nu2), cpp11::as_cpp<cpp11::decay_t<double>>(A), cpp11::as_cpp<cpp11::decay_t<double>>(B), cpp11::as_cpp<cpp11::decay_t<double>>(w), cpp11::as_cpp<cpp11::decay_t<double>>(m), cpp11::as_cpp<cpp11::decay_t<doubles>>(c14_determinations), cpp11::as_cpp<cpp11::decay_t<doubles>>(c14_sigmas), cpp11::as_cpp<cpp11::decay_t<int>>(calcurve_yr_index_offset), cpp11::as_cpp<cpp11::decay_t<doubles>>(mucalallyr), cpp11::as_cpp<cpp11::decay_t<doubles>>(sigcalallyr)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_carbondate_FindInstantPredictiveDensityWalker", (DL_FUNC) &_carbondate_FindInstantPredictiveDensityWalker,  8},
    {"_carbondate_FindPredictiveDensityAndCIWalker",   (DL_FUNC) &_carbondate_FindPredictiveDensityAndCIWalker,   10},
    {"_carbondate_FindPredictiveDensityandCIPolyaUrn", (DL_FUNC) &_carbondate_FindPredictiveDensityandCIPolyaUrn, 12},
    {"_carbondate_PolyaUrnUpdateAlpha_test",           (DL_FUNC) &_carbondate_PolyaUrnUpdateAlpha_test,            5},
    {"_carbondate_PolyaUrnUpdateStep",                 (DL_FUNC) &_carbondate_PolyaUrnUpdateStep,                 20},
    {"_carbondate_WalkerUpdateStep",                   (DL_FUNC) &_carbondate_WalkerUpdateStep,                   20},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_carbondate(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
