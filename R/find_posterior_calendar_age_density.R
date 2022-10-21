# Returns a dataframe containing the `calendar_age`, the `density_mean`
# and the 95% confidence intervals for the density `density_ci_lower` and
# `density_ci_upper`.
.FindPosteriorDensityMeanAndCI <- function(
    output_data, calendar_age_sequence, n_posterior_samples) {

  posterier_density_matrix <- .FindDensityPerSampleID(
    output_data, calendar_age_sequence, n_posterior_samples)

  posterier_density_confidence_intervals <- apply(
    posterier_density_matrix, 1, stats::quantile, probs = c(0.025, 0.975))
  posterier_density_mean <- apply(posterier_density_matrix, 1, mean)

  return(data.frame(
    calendar_age=calendar_age_sequence,
    density_mean=posterier_density_mean,
    density_ci_lower=posterier_density_confidence_intervals[1, ],
    density_ci_upper=posterier_density_confidence_intervals[2, ]))
}


# Creates a matrix where each column is the density for a particular sample id
.FindDensityPerSampleID <- function(
    output_data, calendar_age_sequence, n_posterior_samples) {
  n_out <- length(output_data$alpha)
  n_burn <- floor(n_out / 2)

  posterior_sample_ids <- sample(
    x = n_burn:n_out,
    size = n_posterior_samples,
    replace = n_posterior_samples > (n_out - n_burn))

  posterior_density_matrix <- apply(
    matrix(posterior_sample_ids, 1, n_posterior_samples),
    2,
    function(i, output_data, x) {
      if (output_data$update_type == "walker") {
        .FindPredictedDensityWalker(
          x,
          weight = output_data$weight[[i]],
          phi = output_data$phi[[i]],
          tau = output_data$tau[[i]],
          mu_phi = output_data$mu_phi[i],
          lambda = output_data$lambda,
          nu1 = output_data$nu1,
          nu2 = output_data$nu2)
      } else {
        .FindPredictedDensityNeal(
          x,
          cluster_identifiers = output_data$cluster_identifiers[i, ],
          phi = output_data$phi[[i]],
          tau = output_data$tau[[i]],
          alpha = output_data$alpha[i],
          mu_phi = output_data$mu_phi[i],
          lambda = output_data$lambda,
          nu1 = output_data$nu1,
          nu2 = output_data$nu2)
      }
    },
    output_data = output_data, x = calendar_age_sequence)
  return(posterior_density_matrix)
}


.FindPredictedDensityNeal <- function(
    x, cluster_identifiers, phi, tau, alpha, mu_phi, lambda, nu1, nu2) {
  n_clust <- length(phi)
  nci <- .NumberOfObservationsInEachCluster(cluster_identifiers)

  # Find probability that new theta[i+1] is in a particular cluster or a new one
  pci <- c(nci, alpha) # Could form new cluster
  pci <- pci / sum(pci)

  dens <- .MixtureDens(x, w = pci[1:n_clust], mu = phi, sd = 1 / sqrt(tau)) +
    .PredNewDens(x, w = pci[n_clust + 1], mu_phi, lambda, nu1, nu2)
}


.FindPredictedDensityWalker <- function(
    x, weight, phi, tau, mu_phi, lambda, nu1, nu2) {
  prob_new_clust <- 1 - sum(weight)

  dens <- .MixtureDens(x, w = weight, mu = phi, sd = 1 / sqrt(tau)) +
    .PredNewDens(x, w = prob_new_clust, mu_phi, lambda, nu1, nu2)
  return(dens)
}


# Find mixture density
# Function which where you pass it a set of means, sds and weights and it
# returns the density of the corresponding mixture of normals
# Arguments:
# x - vector of values at which to evaluate mixture density
# w - vector of weights
# mu - the means
# sd - the sds
.MixtureDens <- function(x, w, mu, sd) {
  DTemp <- mapply(
    function(mu, sig, w, x) w * stats::dnorm(x, mean = mu, sd = sig),
    mu,
    sd,
    w,
    MoreArgs = list(x = x))
  apply(DTemp, 1, sum) # Sum up the various mixtures
}


# The predictive for a new observation is a scaled t-distribution
.PredNewDens <- function(x, w, mu_phi, lambda, nu1, nu2) {
  w * exp(.LogMarginalNormalGamma(x, mu_phi, lambda, nu1, nu2))
}

