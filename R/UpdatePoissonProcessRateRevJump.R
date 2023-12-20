#' Title
#'
#' @param theta Vector of the observed times of the events
#' @param rate_s The current changepoints in the rate
#' @param rate_h The current heights in the rate (corresponding to above)
#' @param integrated_rate The integral of piecewise constant rate i.e.,
#' integral_0^L nu(t) dt
#' @param prior_h_shape,prior_h_rate prior parameters on heights
#' We assume each height in piecewise rate has h ~ Gamma(shape, rate)

#' @param prior_n_internal_changepoints_lambda Prior on number of internal changepoints n ~ Po(lambda)
#' @param prob_move Dataframe with probability of each type of RJ move
#' (change_pos, change_height, birth, death)
#'
#' @return TODO
#' @export
#'
#'
#' @examples # TO DO
UpdatePoissonProcessRateRevJump <- function(
    theta,
    rate_s,
    rate_h,
    integrated_rate,
    prior_h_shape,
    prior_h_rate,
    prior_n_internal_changepoints_lambda,
    prob_move) {

  n_changepoints <- length(rate_s)
  n_heights <- length(rate_h)
  n_internal_changepoints <- n_changepoints - 2

  # TODO - CHECK ARGUMENTS OF prob_move

  if(n_heights != n_changepoints - 1) {
    stop("Error in matching dimension of rate_s and rate_h")
  }

  u <- stats::runif(1)

  if(u < prob_move$pos[n_heights]) # Propose moving position of a changepoint
  {
    update <- .ChangePos(
      theta = theta,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate)
    rate_s <- update$rate_s
    integrated_rate <- update$integrated_rate
  }
  else if(u < (
    prob_move$pos[n_heights]
    + prob_move$height[n_heights]
    )) # Propose moving height of a step
  {
    update <- .ChangeHeight(
      theta = theta,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate)
    rate_h <- update$rate_h
    integrated_rate <- update$integrated_rate
  }
  else if(u < (
    prob_move$pos[n_heights]
    + prob_move$height[n_heights]
    + prob_move$birth[n_heights]
    )) # Propose birth step
  {
    proposal_ratio <- (
      (prob_move$death[n_heights + 1] * (rate_s[n_changepoints] - rate_s[1]))
      / (prob_move$birth[n_heights] * (n_internal_changepoints + 1))
    )
    update <- .Birth(
      theta = theta,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape =  prior_h_shape,
      prior_h_rate = prior_h_rate,
      prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
      proposal_ratio = proposal_ratio)
    rate_s <- update$rate_s
    rate_h <- update$rate_h
    integrated_rate <- update$integrated_rate
  }
  else # Propose a death step
  {
    proposal_ratio <- (
      (prob_move$birth[n_heights - 1 ] * n_internal_changepoints)
      / (prob_move$death[n_heights] * (rate_s[n_changepoints] - rate_s[1]))
    )
    update <- .Death(
      theta = theta,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate,
      prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
      proposal_ratio = proposal_ratio)

    rate_s <- update$rate_s
    rate_h <- update$rate_h
    integrated_rate <- update$integrated_rate
  }

  list(rate_s = rate_s, rate_h = rate_h, integrated_rate = integrated_rate)
}



## The 4 proposal reversible jump steps
## 1: Change the position h of a changepoint
## 2: Change the height s of a section
## 3: Remove a changepoint
## 4: Add a changepoint


## Proposal 1: Alter the position of a randomly chosen internal changepoint
# Arguments:
# theta - the observed times (calendar ages) of the events
# rate_s - the current changepoints in the rate
# rate_h - the heights
# integrated_rate - integral_0^L nu(t) dt
.ChangePos <- function(
    theta,
    rate_s,
    rate_h,
    integrated_rate)
{
  n_changepoints <- length(rate_s)

  if(n_changepoints < 3) {
    stop("Internal Error (ChangePos): Proposed to move an internal changepoint when there are none")
  }

  # Select internal changepoint to move at random
  j <- .resample(2:(n_changepoints-1), 1) # Need careful version

  # Propose new changepoint position
  rate_s_new <- rate_s
  rate_s_new[j] <- stats::runif(1, min = rate_s[j - 1], max = rate_s[j + 1])

  # Find prior ratio for rate s
  prior_rate_s_new   <- (
    (rate_s_new[j + 1] - rate_s_new[j]) * (rate_s_new[j] - rate_s_new[j - 1]))
  prior_rate_s_old <- (
    (rate_s[j + 1] - rate_s[j]) * (rate_s[j] - rate_s[j - 1]))

  prior_rate_s_ratio <- prior_rate_s_new / prior_rate_s_old

  # Find the likelihood of the theta data given both sets of changepoints (only changes in period)

  # Find new integrated_rate by adjusting previous integrated_rate
  # Calculates as (Diff in changepoint location) * (Diff in height)
  adjust_integral <- (rate_h[j] - rate_h[j - 1]) * (rate_s[j] - rate_s_new[j])
  integrated_rate_new <- integrated_rate + adjust_integral

  # Find which thetas will contribute to the likelihood ratio and find their log-likelihood
  if(rate_s_new[j] < rate_s[j]) { # Have shifted new changepoint towards t = 0
    # Find the range of cal ages in which the rate has changed
    min_sj <- rate_s_new[j]
    max_sj <- rate_s[j]
    # Find the value of the rate in this interval (both existing and proposed)
    old_h_interval <- rate_h[j-1]
    new_h_interval <- rate_h[j]
  } else {  # Have shifted new changepoint away from t = 0
    # Find the range of cal ages in which the rate has changed
    min_sj <- rate_s[j]
    max_sj <- rate_s_new[j]
    # Find the value of the rate in this interval (note swap from above)
    old_h_interval <- rate_h[j]
    new_h_interval <- rate_h[j-1]
  }

  # Number of thetas that have been affected by changepoint shift
  n_theta_affected <- sum(theta > min_sj & theta < max_sj)

  log_lik_old <- (n_theta_affected * log(old_h_interval)) - integrated_rate
  log_lik_new <- (n_theta_affected * log(new_h_interval)) - integrated_rate_new
  theta_lik_ratio <- exp(log_lik_new - log_lik_old)

  # Find acceptance probability
  hastings_ratio <- theta_lik_ratio * prior_rate_s_ratio

  # Determine acceptance and return result
  if(stats::runif(1) < hastings_ratio) {
    # Accept
    retlist <- list(
      rate_s = rate_s_new,
      integrated_rate = integrated_rate_new,
      hastings_ratio = hastings_ratio)
  } else {
    # Reject
    retlist <- list(
      rate_s = rate_s,
      integrated_rate = integrated_rate,
      hastings_ratio = hastings_ratio)
  }
  return(retlist)
}


## Proposal 2: Alter the height of a randomly chosen step
# Arguments:
# theta - the observed times (calendar ages) of the events
# rate_s - the current changepoints in the rate
# rate_h - the heights
# integrated_rate - integral_0^L nu(t) dt
# prior_h_shape, prior_h_rate - prior parameters on heights h ~ Gamma(shape, rate)
# Heights h are drawn from Gamma(alpha, beta) distribution
.ChangeHeight <- function(
    theta,
    rate_s,
    rate_h,
    integrated_rate,
    prior_h_shape,
    prior_h_rate)
{
  n_heights <- length(rate_h)

  # Select step height to alter - will change height between s[j] and s[j+1]
  j <- sample(n_heights, 1)
  # Can use sample as just pass a integer so will pick from 1:n_heights

  # Propose new height of step so that log(h_j_new/h_j_old) ~ Unif[-0.5, 0.5]
  h_j_old <- rate_h[j]
  u <- stats::runif(1, min = -0.5, max = 0.5)
  h_j_new <- h_j_old * exp(u)

  # Store new set of heights
  rate_h_new <- rate_h
  rate_h_new[j] <- h_j_new

  # Find prior h ratio
  log_prior_h_ratio <- (prior_h_shape * u) - (prior_h_rate * (h_j_new - h_j_old))

  # Adjust the integrated rate to account for new height between s[j] and s[j+1]
  integrated_rate_new <- integrated_rate + (h_j_new - h_j_old)*(rate_s[j+1] - rate_s[j])

  # Number of thetas that have been affected by changepoint shift
  n_theta_affected <- sum(theta > rate_s[j] & theta < rate_s[j+1])

  log_lik_old <- (n_theta_affected * log(h_j_old)) - integrated_rate # In old these have rate h_j_old
  log_lik_new <- (n_theta_affected * log(h_j_new)) - integrated_rate_new # In new they have rate h_j_new
  log_theta_lik_ratio <- log_lik_new - log_lik_old

  # Find acceptance probability (use log rather than multiply as logphratio is simple format)
  hastings_ratio <- exp(log_prior_h_ratio + log_theta_lik_ratio)

  # Determine acceptance and return result
  if(stats::runif(1) < hastings_ratio)	{ 	# Accept
    retlist <- list(
      rate_h = rate_h_new,
      integrated_rate = integrated_rate_new,
      hastings_ratio = hastings_ratio)
  }	else { # Reject
    retlist <- list(
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      hastings_ratio = hastings_ratio)
  }
  return(retlist)
}


## Proposal 3: Giving birth to a new changepoint
# Arguments:
# theta - the observed times (calendar ages) of the events
# rate_s - the current changepoints in the rate
# rate_h - the heights
# integrated_rate - integral_0^L nu(t) dt
# prior_h_shape, prior_h_rate - prior parameters on heights h ~ Gamma(shape, rate)
# prior_n_internal_changepoints_lambda - prior on number of changepoints n ~ Po(lambda)
# proposal_ratio - proposal ratio for an additional changepoint
.Birth <- function(
    theta,
    rate_s,
    rate_h,
    integrated_rate,
    prior_h_shape,
    prior_h_rate,
    prior_n_internal_changepoints_lambda,
    proposal_ratio)
{
  n_changepoints <- length(rate_s)
  n_internal_changepoints <- n_changepoints - 2

  # Propose new changepoint
  s_star <- stats::runif(1, min = rate_s[1], max = rate_s[n_changepoints])

  # Find which interval it falls into
  j <- max(which(rate_s < s_star)) # Don't need to worry about boundary

  # Current rate height between s_j and s_{j+1}
  h_j_old <- rate_h[j]

  # Sample u = U[0,1] and adjust heights
  u <- stats::runif(1)
  h_A_new <- h_j_old * (u/(1-u))^((rate_s[j+1]- s_star)/(rate_s[j+1] - rate_s[j]))
  h_B_new <- h_A_new * (1-u) / u

  # Create new changepoint and height vectors in sorted order
  rate_s_new <- c(rate_s[1:j], s_star, rate_s[(j+1):n_changepoints])
  # No care as no change to first/last element i.e. j = 1, ns -1
  rate_h_new <- append(rate_h[-j], c(h_A_new, h_B_new), after = j-1)
  # Care as could change first or last elements

  # Find the prior ratio for dimension
  log_prior_num_change_ratio <- (
    stats::dpois(n_internal_changepoints + 1, prior_n_internal_changepoints_lambda, log = TRUE)
    - stats::dpois(n_internal_changepoints, prior_n_internal_changepoints_lambda, log = TRUE)
  )

  prior_spacing_ratio <- (
    (
      2 * (n_internal_changepoints + 1)
      * (2 * n_internal_changepoints + 3) / (rate_s[n_changepoints] - rate_s[1])^2
    )
    * (s_star - rate_s[j]) * (rate_s[j+1] - s_star) / (rate_s[j+1] - rate_s[j])
  )

  # Find the prior ratio for the heights NEED CARE WITH ROUNDING
  prior_h_ratio <- (
    (
      (prior_h_rate ^ prior_h_shape) / gamma(prior_h_shape)
    )
    * exp(
      (prior_h_shape-1)*(log(h_A_new)+log(h_B_new)-log(h_j_old))
      - prior_h_rate*(h_A_new + h_B_new - h_j_old)
    )
  )

  jacobian <- (h_A_new + h_B_new)^2 / h_j_old

  # Find the likelihood of thetas
  integrated_rate_adjustment <- (
    ((h_A_new - h_j_old) * (s_star - rate_s[j]))
    + ((h_B_new - h_j_old) * (rate_s[j+1] - s_star))
  )
  integrated_rate_new <- integrated_rate + integrated_rate_adjustment

  n_theta_affected_A <- sum(theta <= s_star & theta > rate_s[j])
  n_theta_affected_B <- sum(theta < rate_s[j+1] & theta > s_star)

  log_lik_old <- ((n_theta_affected_A + n_theta_affected_B) * log(h_j_old)) - integrated_rate # All have rate h_j_old
  log_lik_new <- (n_theta_affected_A * log(h_A_new)) + (n_theta_affected_B * log(h_B_new)) - integrated_rate_new # In new have rate h_A_new or h_B_new dependent upon if after sstar

  log_theta_lik_ratio <- log_lik_new - log_lik_old

  # Find acceptance probability
  hastings_ratio <- (
    exp(log_prior_num_change_ratio + log_theta_lik_ratio)
    * (prior_spacing_ratio * prior_h_ratio * jacobian * proposal_ratio)
  )

  # Determine acceptance and return result
  if(stats::runif(1) < hastings_ratio)	{ # Accept
    retlist <- list(
      rate_s = rate_s_new,
      rate_h = rate_h_new,
      integrated_rate = integrated_rate_new,
      hastings_ratio = hastings_ratio)
  } else { # Reject
    retlist <- list(
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      hastings_ratio = hastings_ratio)
  }
  return(retlist)
}


## Proposal 4: Killing a current changepoint
# Arguments:
# theta - the observed times (calendar ages) of the events
# rate_s - the current changepoints in the rate
# rate_h - the heights
# integrated_rate - integral_0^L nu(t) dt
# prior_h_shape, prior_h_rate - prior parameters on heights h ~ Gamma(shape, rate)
# prior_n_internal_changepoints_lambda - prior on number of internal changepoints n ~ Po(lambda)
# proposal_ratio - proposal ratio for an additional changepoint
.Death <- function(
    theta,
    rate_s,
    rate_h,
    integrated_rate,
    prior_h_shape,
    prior_h_rate,
    prior_n_internal_changepoints_lambda,
    proposal_ratio)
{
  n_changepoints <- length(rate_s)
  n_internal_changepoints <- n_changepoints - 2

  if(n_internal_changepoints <= 0) {
    stop("Error: You have called death when there are no internal changepoints")
  }

  # Select changepoint to remove
  j <- .resample(2:(n_changepoints-1), 1) # Care as 2:(ns-1) can be a single integer

  # Create new height for interval s[j-1], s[j+1] via inverse of birth step
  h_j_new <- exp(
    1/(rate_s[j+1] - rate_s[j-1]) * (
      (rate_s[j] - rate_s[j-1]) * log(rate_h[j-1]) + (rate_s[j+1]- rate_s[j]) * log(rate_h[j])
    )
  )

  rate_s_new <- rate_s[-j]
  rate_h_new <- append(rate_h[-c(j-1, j)], h_j_new, after = j-2) # Care as could change first/last element

  # Find the prior ratio for dimension
  log_prior_num_change_ratio <- (
    stats::dpois(n_internal_changepoints - 1 , prior_n_internal_changepoints_lambda, log = TRUE)
    - stats::dpois(n_internal_changepoints, prior_n_internal_changepoints_lambda, log = TRUE)
  )

  prior_spacing_ratio <- (
    (rate_s[n_changepoints] - rate_s[1])^2 / (2 * n_internal_changepoints
                                              * (2*n_internal_changepoints+1) )
    * (rate_s[j+1] - rate_s[j-1]) / ((rate_s[j+1] - rate_s[j]) * (rate_s[j] - rate_s[j-1]))
  )

  # Find the prior ratio for the heights NEED CARE WITH ROUNDING
  prior_h_ratio <- (
    (gamma(prior_h_shape) / (prior_h_rate^prior_h_shape))
    / exp(
      (prior_h_shape - 1) * (log(rate_h[j]) + log(rate_h[j-1]) - log(h_j_new))
      - prior_h_rate * (rate_h[j] + rate_h[j-1] - h_j_new)
    )
  )

  jacobian <- h_j_new / ((rate_h[j-1] + rate_h[j])^2)

  # Find likelihood of thetas
  integrated_rate_adjustment <- (
    ((h_j_new - rate_h[j-1]) * (rate_s[j] - rate_s[j-1]))
    + ((h_j_new - rate_h[j]) * (rate_s[j+1] - rate_s[j]))
  )
  integrated_rate_new <- integrated_rate + integrated_rate_adjustment

  n_theta_affected_A <- sum(theta <= rate_s[j] & theta > rate_s[j-1]) # in old will have rate h[j-1]
  n_theta_affected_B <- sum(theta < rate_s[j+1] & theta > rate_s[j])  # in old will have rate h[j]

  log_lik_old <- (n_theta_affected_A * log(rate_h[j-1])) + (n_theta_affected_B * log(rate_h[j]))  - integrated_rate # In old they have rate h[j-1] or h[j] dependent upon if after s[j]
  log_lik_new <- ((n_theta_affected_A + n_theta_affected_B) * log(h_j_new)) - integrated_rate_new # In new all have rate h_j_new

  log_theta_lik_ratio <- log_lik_new - log_lik_old


  # Find acceptance probability
  hastings_ratio <- (
    exp(log_theta_lik_ratio + log_prior_num_change_ratio)
    * (prior_spacing_ratio * prior_h_ratio * jacobian * proposal_ratio)
  )

  # Determine acceptance and return result
  if(stats::runif(1) < hastings_ratio)	{ # Accept
    retlist <- list(
      rate_s = rate_s_new,
      rate_h = rate_h_new,
      integrated_rate = integrated_rate_new,
      hastings_ratio = hastings_ratio)

  } else { # Reject
    retlist <- list(
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      hastings_ratio = hastings_ratio)
  }
  return(retlist)
}

### .FindMoveProbability
# Work out the probabilities for each move in the RJ MCMC sampler
# Return 4 vectors in list with each move probability for current n_heights = 1,..., kmax+1
# Note: Uses number of heights to avoid confusion as sometimes nk = 0
# and R does not create vectors with index p[0]
# Note: nk = nh - 1 (i.e. n_internal_changepoints = n_heights - 1)
# Arguments:
# prior_n_internal_changepoints_lambda - prior mean on number of changepoints
# k_max_internal_changepoints - maximum number of changepoints permitted
# rescale_factor - comparison of birth/death vs changing height/position
.FindMoveProbability <- function(
    prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints,
    rescale_factor = 0.9)
{
  prob_move_pos <- rep(NA, length = k_max_internal_changepoints + 1)
  prob_move_height <- rep(NA, length = k_max_internal_changepoints + 1)
  prob_move_birth <- rep(NA, length = k_max_internal_changepoints + 1)
  prob_move_death <- rep(NA, length = k_max_internal_changepoints + 1)

  # Fixed constraints
  prob_move_birth[k_max_internal_changepoints + 1] <- 0
  prob_move_death[1] <- 0

  # Now find other probabilities of birth and death ignoring constant c
  prob_move_birth[1:k_max_internal_changepoints] <- pmin(
    1,
    (stats::dpois(1:k_max_internal_changepoints, lambda = prior_n_internal_changepoints_lambda)
     / stats::dpois(0:(k_max_internal_changepoints - 1), lambda = prior_n_internal_changepoints_lambda))
  )
  prob_move_death[2:(k_max_internal_changepoints + 1)] <- pmin(
    1,
    (stats::dpois(0:(k_max_internal_changepoints - 1), lambda = prior_n_internal_changepoints_lambda)
     / stats::dpois(1:k_max_internal_changepoints, lambda = prior_n_internal_changepoints_lambda))
  )

  # Rescale to allow other moves a reasonable probability
  rescale_constant <- rescale_factor / max(prob_move_birth + prob_move_death)
  prob_move_birth <- rescale_constant * prob_move_birth
  prob_move_death <- rescale_constant * prob_move_death

  prob_move_pos <- prob_move_height <- (1 - prob_move_birth - prob_move_death)/2

  # Deal with case that if n_heights = 1 (i.e. n_internal_changes = 0)
  # then cannot move position of a changepoint
  prob_move_pos[1] <- 0
  prob_move_height[1] <- 2 * prob_move_height[1]

  list(
    pos = prob_move_pos,
    height = prob_move_height,
    birth = prob_move_birth,
    death = prob_move_death)
}




