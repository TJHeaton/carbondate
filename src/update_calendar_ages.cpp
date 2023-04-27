#include <limits>
#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;


double ThetaLogLikelihood_cpp(
    double theta,
    double prmean,
    double prsig,
    double c14obs,
    double c14sig,
    int year_index_offset,
    const doubles& mucalallyr,
    const doubles& sigcalallyr) {

  double loglik;
  double mucal, sigcal;
  int yr_index = (int) theta - year_index_offset;

  if ((yr_index < 0) | (yr_index >= mucalallyr.size())) {
    return -std::numeric_limits<double>::infinity();
  }

  mucal = mucalallyr[yr_index];
  sigcal= sigcalallyr[yr_index];

  loglik = Rf_dnorm4(theta, prmean, prsig, 1);
  loglik += Rf_dnorm4(c14obs, mucal, sqrt(sigcal*sigcal + c14sig*c14sig), 1);

  return loglik;
}


double SliceSample_cpp(
    double x0,
    double w,
    double m,
    double prmean,
    double prsig,
    double c14obs,
    double c14sig,
    int year_index_offset,
    const doubles& mucalallyr,
    const doubles& sigcalallyr) {

  double y;     // Slice height
  double L, R;  // Each side of calendar age slice interval
  double J, K;  // Max steps on left and right hand sides
  double x1;    // Current sampled value

  //////////////////////////////////////////////
  // Slice height
  y = ThetaLogLikelihood_cpp(
    x0, prmean, prsig, c14obs, c14sig, year_index_offset, mucalallyr, sigcalallyr);
  y -= Rf_rexp(1);

  //////////////////////////////////////////////
  // Find the slice interval

  // Initialise slice
  L = x0 - w * Rf_runif(0, 1);
  R = L + w;

  // Decide how far we can extend on both sides
  J = floor(m * Rf_runif(0, 1));
  K = m - 1 - J;

  // LHS stepping out
  while ((J > 0) && (y < ThetaLogLikelihood_cpp(
      L, prmean, prsig, c14obs, c14sig, year_index_offset, mucalallyr, sigcalallyr))) {
    L -= w;
    J -= 1;
  }

  // RHS stepping out
  while ((K > 0) && (y < ThetaLogLikelihood_cpp(
      R, prmean, prsig, c14obs, c14sig, year_index_offset, mucalallyr, sigcalallyr))) {
    R += w;
    K -= 1;
  }

  //////////////////////////////////////////////
  // Get sampled value from the slice interval
  while (true) {
    x1 = L + Rf_runif(0, 1) * (R - L);

    // Break loop if we have sampled satisfactorily
    if (y < ThetaLogLikelihood_cpp(
        x1, prmean, prsig, c14obs, c14sig, year_index_offset, mucalallyr, sigcalallyr)) {
      return x1;
    }

    // Else shrink the interval
    if (x1 < x0) {
      L = x1;
    } else {
      R = x1;
    }
  }
}


std::vector<double> UpdateCalendarAges(
    int n,
    const doubles& calendar_ages,
    double w,
    double m,
    const std::vector<int>& cluster_identifiers,
    const std::vector<double>& phi,
    const std::vector<double>& tau,
    const doubles& c14_determinations,
    const doubles& c14_sigmas,
    int year_index_offset,
    const doubles& mucalallyr,
    const doubles& sigcalallyr) {

  std::vector<double> calendar_ages_new(n);
  double prmean;
  double prsig;
  int ci;

  for (int k = 0; k < n; ++k) {
    ci = cluster_identifiers[k];
    prmean = phi[ci-1];
    prsig = 1.0 / sqrt(tau[ci-1]);

    calendar_ages_new[k] = SliceSample_cpp(
      calendar_ages[k],
      w,
      m,
      prmean,
      prsig,
      c14_determinations[k],
      c14_sigmas[k],
      year_index_offset,
      mucalallyr,
      sigcalallyr);
  }
  return calendar_ages_new;
}
