#include <vector>
#include "Rmath.h"

double Mean(const std::vector<double>& vec) {
  double mean = 0.0;
  for (double elem : vec) {
    mean += elem;
  }
  mean /= (double) vec.size();
  return mean;
}


// Function which works out the marginal of theta when
// theta ~ N(phi, sd = sqrt(1/tau)) and (phi,tau) are NormalGamma
double LogMarginalNormalGamma(
    double calendar_age,
    double lambda,
    double nu1,
    double nu2,
    double mu_phi) {

  double logden, margprec, margdf;

  margprec = (nu1 * lambda) / (nu2 * (lambda + 1.));
  margdf = 2. * nu1;

  logden = lgamma((margdf + 1.) / 2.) - lgamma(margdf / 2.);
  logden += 0.5 * (log(margprec) - log(margdf) - log(M_PI));
  logden -= ((margdf + 1) / 2) * log(1 + margprec * pow(calendar_age - mu_phi, 2) / margdf);

  return logden;
}


// Samples from the posterior of (phi, tau) given a set of calendar ages from that cluster
// For each cluster in turn we have calendar_age ~ N(phi, sd = 1/sqrt(tau))
void UpdatePhiTau(
    const std::vector<double>& calendar_ages,  // calendar ages all belonging to the same cluster
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    double& phi,
    double& tau) {

  int nclus = calendar_ages.size();
  std::vector<double> thetadiff(nclus);
  double thetabar;
  double s;
  double mu_phi_new;
  double lambda_new, nu1_new, nu2_new;

  thetabar = Mean(calendar_ages);
  for (int i = 0; i < nclus; ++i) {
    thetadiff[i] = pow(calendar_ages[i] - thetabar, 2);
  }
  s = Mean(thetadiff);

  // Update parameters according to conjugate prior
  nu1_new = nu1 + nclus / 2.;
  nu2_new = nu2;
  nu2_new += 0.5 * nclus * (s + lambda * pow(thetabar - mu_phi, 2)/(lambda + nclus));

  lambda_new = lambda + nclus;
  mu_phi_new = (lambda * mu_phi + nclus * thetabar) / (lambda + nclus);

  tau = Rf_rgamma(nu1_new, 1./nu2_new);
  phi = Rf_rnorm(mu_phi_new, 1./sqrt(lambda_new * tau));
}


// Updates mu_phi via Gibbs sampling based upon current (phi, tau) values
double UpdateMuPhi(
    const std::vector<double>& phi,   // means of clusters
    const std::vector<double>& tau,   // precisions of clusters
    double lambda,
    double A,      // prior mean of mu_phi
    double B) {    // prior precision of mu_phi

  double mu_phi;
  double posterior_mean, posterior_precision;
  int n = phi.size();
  double sum_tau = 0.;
  double sum_tau_mult_phi = 0.;

  for (int j = 0; j < n; j++) {
    sum_tau += tau[j];
    sum_tau_mult_phi += tau[j] * phi[j];
  }

  posterior_mean = (A * B + lambda * sum_tau_mult_phi) / (B + lambda * sum_tau);
  posterior_precision = B + lambda * sum_tau;
  mu_phi = Rf_rnorm(posterior_mean, 1. / sqrt(posterior_precision));
  return mu_phi;
}
