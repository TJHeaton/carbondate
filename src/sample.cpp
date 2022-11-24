#include <vector>
#include "Rmath.h"
#include "R.h"

// Adapted from `ProbSampleNoReplace` in R/src/main/sample.c and from SampleNoReplace
// in Rcpp package, but substantially simplified since we know we're only ever
// going to be calling this with sz = 1. Tested that it gives the same result as
// calling sample.int from R.
int SampleInt(int n, std::vector<double> prob, bool one_based) {

  int adj = one_based ? 0 : 1;
  double rT, mass, sum_p = 0.;
  int i, j, ans;
  std::vector<double> p(n);
  std::vector<int> perm(n);

  for (i = 0; i < n; i++) {
    perm[i] = i + 1;
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

  ans = perm[j] - adj;
  return ans;
}
