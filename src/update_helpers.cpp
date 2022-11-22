#include <limits>
#include <algorithm>
#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;



[[cpp11::register]] integers which_equal(integers vec, int i) {
  std::vector<int> indices;
  indices.reserve(20);
  for (int k = 0; k < vec.size(); ++k) {
    if (vec[k] == i) {
      indices.push_back(k + 1);
    }
  }
  return integers(indices);
}


[[cpp11::register]] integers which_mt(doubles vec, double i) {
  std::vector<int> indices;
  indices.reserve(20);
  for (int k = 0; k < vec.size(); ++k) {
    if (vec[k] > i) {
      indices.push_back(k + 1);
    }
  }
  return integers(indices);
}

[[cpp11::register]] integers which_mt_int(integers vec, int i) {
  std::vector<int> indices;
  indices.reserve(20);
  for (int k = 0; k < vec.size(); ++k) {
    if (vec[k] > i) {
      indices.push_back(k + 1);
    }
  }
  return integers(indices);
}


double min_value(std::vector<double> vec) {
  auto pos = std::min_element(vec.begin(), vec.end());
  return *pos;
}


bool any_more_than(integers vec, int j) {
  return std::any_of(vec.cbegin(), vec.cend(), [=](int i){ return i > j; });
}


// Adapted from `ProbSampleNoReplace` in R/src/main/sample.c and from SampleNoReplace
// in Rcpp package, but substantially simplified since we know we're only ever
// going to be calling this with sz = 1. Tested that it gives the same result as
// calling sample.int from R.
[[cpp11::register]] int sample_int(int n, doubles prob, bool one_based) {

  local_rng rng_state;
  int adj = one_based ? 0 : 1;
  double rT, mass, sum_p = 0.;
  int i, j, ans;
  std::vector<double> p(n);   // Possible not not initialise?
  std::vector<int> perm(n);   // Possible not not initialise?

  for (i = 0; i < n; i++) {
    perm[i] = i + 1;
    p[i] = prob[i];
    sum_p += p[i];
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


double find_a(integers cluster_identifiers, int clust_num, double brprod, std::vector<double> u) {
  double a = 0;

  for (int k = 0; k < cluster_identifiers.size(); ++k) {
    if ((cluster_identifiers[k] == clust_num) && (u[k] > a)) {
      a = u[k];
    }
  }
  return a/brprod;
}


double find_b(
    integers cluster_identifiers,
    int clust_num,
    std::vector<double> u,
    std::vector<double> v) {

  int n = cluster_identifiers.size();
  int m = v.size();
  writable::doubles prodv(m + 1);
  double b = 1.0;
  double b_sub = 0.0;
  double b_sub_comp;
  int index;

  if (any_more_than(cluster_identifiers, clust_num)) {
    prodv[0] = 1.0 - v[0];
    for (int k = 1; k < m; ++k) {
      prodv[k] = prodv[k-1] * (1.0 - v[k]);
      prodv[k-1] /= (1.0 - v[clust_num - 1]);
    }
    prodv[m-1] /= (1.0 - v[clust_num - 1]);

    for (int k = 0; k < n; ++k) {
      if (cluster_identifiers[k] > clust_num) {
        index = cluster_identifiers[k] - 1;
        b_sub_comp = u[k] / (v[index] * prodv[index - 1]);
        if (b_sub_comp > b_sub) {
          b_sub = b_sub_comp;
        }
      }
    }
    b -= b_sub;
  }

  return b;
}

