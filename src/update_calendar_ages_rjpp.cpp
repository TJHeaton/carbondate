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
int SampleInt(int n, doubles prob) {

  double rT, mass, sum_p = 0.;
  int i, j;
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


[[cpp11::register]] doubles UpdateCalendarAgesGibbs(
    doubles prior_calendar_ages,
    doubles calendar_age_grid,
    list likelihood_calendar_ages_from_calibration_curve,
    integers likelihood_offsets) {

  local_rng rng_state;              // Ensures RNG follows R and R follows after
  int n_obs = likelihood_calendar_ages_from_calibration_curve.size();
  cpp11::writable::doubles updated_calendar_ages(n_obs);

  for (int i = 0; i < n_obs; i++) {
    cpp11::writable::doubles likelihood = likelihood_calendar_ages_from_calibration_curve[i];

    for (int j = 0; j < likelihood.size(); j++) likelihood[j] *= prior_calendar_ages[j + likelihood_offsets[i]];

    int new_index = SampleInt(likelihood.size(), likelihood) + likelihood_offsets[i];
    updated_calendar_ages[i] = calendar_age_grid[new_index];
  }

  return updated_calendar_ages;
}
