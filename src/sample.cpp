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
  int i, j;
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

  return perm[j] - adj;
}

// Adapted from do_sample in in R/src/main/sample.c and from EmpiricalSample in Rcpp package
std::vector<int> GetSampleIds(int start_index, int finish_index, int size) {

  int n = finish_index - start_index + 1;
  bool replace = size >= n;
  std::vector<int> ans(size);

  if (replace || size < 2) {
    for (int i = 0 ; i < size; i++) {
      ans[i] = static_cast<int>(R_unif_index(n)) + start_index;
    }
    return ans;
  }

  std::vector<int> x(n);
  for (int i = 0; i < n; i++) {
    x[i] = i;
  }

  for (int i = 0 ; i < size; i++) {
    int j = static_cast<int>(R_unif_index(n));
    ans[i] = x[j] + start_index;
    x[j] = x[--n];
  }

  return ans;
}
