CombineRuns <- function(output_data_list) {

  n_iter <- output_data_list[[1]]$input_parameters$n_iter
  n_thin <- output_data_list[[1]]$input_parameters$n_thin
  n_obs <- length(output_data_list[[1]]$input_data$rc_determinations)

  pre_burn_rate_s <- list()
  pre_burn_rate_h <- list()
  post_burn_rate_s <- list()
  post_burn_rate_h <- list()
  pre_burn_n_internal_changes <- c()
  post_burn_n_internal_changes <- c()
  combined_n_iter <- 0

  for(i in 1:length(output_data_list)) {

    individual_rate_s <- output_data_list[[i]]$rate_s
    individual_rate_h <- output_data_list[[i]]$rate_h
    n_internal_changes <- output_data_list[[i]]$n_internal_changes
    n_samples <- length(individual_rate_s)
    n_burn <- ceiling(n_samples/2)


    pre_burn_id <- 1:n_burn
    post_burn_id <- (n_burn + 1):n_samples

    pre_burn_rate_s <- c(pre_burn_rate_s,
                         individual_rate_s[pre_burn_id])
    post_burn_rate_s <- c(post_burn_rate_s,
                          individual_rate_s[post_burn_id])

    pre_burn_rate_h <- c(pre_burn_rate_h,
                         individual_rate_h[pre_burn_id])
    post_burn_rate_h <- c(post_burn_rate_h,
                          individual_rate_h[post_burn_id])

    pre_burn_n_internal_changes <- c(pre_burn_n_internal_changes,
                                     n_internal_changes[pre_burn_id])
    post_burn_n_internal_changes <- c(post_burn_n_internal_changes,
                                     n_internal_changes[post_burn_id])

    combined_n_iter <- combined_n_iter + output_data_list[[i]]$input_parameters$n_iter

  }

  combined_output <- combined_output <- output_data_list[[1]]
  combined_output$rate_s <- c(pre_burn_rate_s, post_burn_rate_s)
  combined_output$rate_h <- c(pre_burn_rate_h, post_burn_rate_h)
  combined_output$n_internal_changes <- c(pre_burn_n_internal_changes,
                                          post_burn_n_internal_changes)
  combined_output$input_parameters$n_iter <- combined_n_iter

  return(combined_output)
}

