#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "R.h"
#include "local_rng.h"
using namespace cpp11;

// Adapted from `ProbSampleNoReplace` in R/src/main/sample.c and from SampleNoReplace
// in Rcpp package, but substantially simplified since we know we're only ever
// going to be calling this with sz = 1. Tested that it gives the same result as
// calling sample.int from R.
int SampleInt(std::vector<double> &prob) {

  double rT, mass, sum_p = 0.;
  int i, j, n = prob.size();
  std::vector<double> p(n);
  std::vector<int> perm(n);

  for (i = 0; i < n; i++) {
    perm[i] = i;
    if (R_FINITE(prob[i]) && prob[i] > 0.0) {
      p[i] = prob[i];
      sum_p += p[i];
    } else {
      p[i] = 0.0;
    }
  }
  Rf_revsort(&p[0], &perm[0], n);

  rT = unif_rand() * sum_p;
  mass = 0.0;
  for (j = 0; j < n-1; j++) {
    mass += p[j];
    if (rT <= mass) {
      break;
    }
  }

  return perm[j];
}


[[cpp11::register]] doubles UpdateCalendarAgesGibbsCpp(
    doubles prior_calendar_ages,
    doubles calendar_age_grid,
    list likelihood_calendar_ages_from_calibration_curve,
    integers likelihood_offset) {

  local_rng rng_state;              // Ensures RNG follows R and R follows after
  int n_obs = likelihood_calendar_ages_from_calibration_curve.size();
  cpp11::writable::doubles updated_calendar_age(n_obs);

  for (int i = 0; i < n_obs; i++) {
    doubles likelihood = likelihood_calendar_ages_from_calibration_curve[i];
    std::vector<double> posterior_cal_age(likelihood.size());

    for (int j = 0; j < likelihood.size(); j++) {
      posterior_cal_age[j] = likelihood[j] * prior_calendar_ages[j + likelihood_offset[i]];
    }

    updated_calendar_age[i] = calendar_age_grid[SampleInt(posterior_cal_age) + likelihood_offset[i]];
  }

  return updated_calendar_age;
}
