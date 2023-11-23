library(devtools)
library(proftools)
load_all()

set.seed(14)

n_iter <- 1e5
n_thin <- 10

profile_walker <- function(n_iter, n_thin, calculate_density_for_convergence) {
  file_pref <- paste("Walker", n_iter, n_thin, "conv", calculate_density_for_convergence, sep="_")

  Rprof(paste("profiling/", file_pref, ".out", sep = ""))
  walker_output <- WalkerBivarDirichlet(
    kerr$c14_age,
    kerr$c14_sig,
    FALSE,
    intcal20,
    n_iter,
    n_thin,
    calculate_density_for_convergence = calculate_density_for_convergence)
  Rprof(NULL)

  pd <- readProfileData(paste("profiling/", file_pref, ".out", sep = ""))
  flat_profile <- flatProfile(pd)
  print(flat_profile)
  save(flat_profile, file=paste("profiling/", file_pref, ".rda", sep = ""))

  profileCallGraph2Dot(pd, score = "total", filename = paste("profiling/", file_pref, ".dot", sep = ""))
}

profile_walker(n_iter, n_thin, TRUE)
profile_walker(n_iter, n_thin, FALSE)
