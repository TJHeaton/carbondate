#include "cpp11/doubles.hpp"
#include "Rmath.h"
#include <limits>
using namespace cpp11;


[[cpp11::register]] double random_test_cpp(void) {
  return Rf_runif(0, 1);
}


[[cpp11::register]] double ThetaLogLikelihood_cpp(
    double theta,
    double prmean,
    double prsig,
    double c14obs,
    double c14sig,
    doubles mucalallyr,
    doubles sigcalallyr) {

  double loglik;
  double mucal, sigcal;
  int yr;

  yr = floor(theta) - 1;

  if (yr < 0 | yr >= mucalallyr.size()) {
    return -std::numeric_limits<double>::infinity();
  }

  mucal = mucalallyr[yr];
  sigcal= sigcalallyr[yr];

  loglik = Rf_dnorm4(theta, prmean, prsig, 1);
  loglik += Rf_dnorm4(c14obs, mucal, sqrt(sigcal*sigcal + c14sig*c14sig), 1);

  return loglik;
}


[[cpp11::register]] double SliceSample_cpp(
    double x0,
    double w,
    double m,
    double prmean,
    double prsig,
    double c14obs,
    double c14sig,
    doubles mucalallyr,
    doubles sigcalallyr) {

    double y;     // Slice height
    double L, R;  // Each side of calender age slice interval
    double J, K;  // Max steps on left and right hand sides
    double x1;    // Current sampled value

  //////////////////////////////////////////////
  // Slice height
  y = ThetaLogLikelihood_cpp(
    x0, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr);
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
  while (J > 0 & y < ThetaLogLikelihood_cpp(
      L, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr)) {
    L -= w;
    J -= 1;
  }

  // RHS stepping out
  while (K > 0 & y < ThetaLogLikelihood_cpp(
      R, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr)) {
    R += w;
    K -= 1;
  }

  //////////////////////////////////////////////
  // Get sampled value from the slice interval
  while (1) {
    x1 = L + Rf_runif(0, 1) * (R - L);

    // Break loop if we have sampled satisfactorily
    if (y < ThetaLogLikelihood_cpp(
        x1, prmean, prsig, c14obs, c14sig, mucalallyr, sigcalallyr)) {
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



