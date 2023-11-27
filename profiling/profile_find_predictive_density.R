library(proftools)

set.seed(14)

n_iter <- 1e4
n_thin <- 10

profile_walker <- function(n_iter, n_thin, use_cpp) {
  file_pref <- paste("Density", n_iter, n_thin, "cpp", use_cpp, sep="_")

  walker_output <- WalkerBivarDirichlet(
    kerr$c14_age,
    kerr$c14_sig,
    intcal20,
    n_iter,
    n_thin,
    use_cpp = use_cpp)

  polya_urn_output <- PolyaUrnBivarDirichlet(
    kerr$c14_age,
    kerr$c14_sig,
    intcal20,
    n_iter,
    n_thin,
    use_cpp = use_cpp)

  Rprof(paste0("profiling/", file_pref, ".out"))
  PlotPredictiveCalendarAgeDensity(list(walker_output, polya_urn_output), 5000)
  Rprof(NULL)

  pd <- readProfileData(paste0("profiling/", file_pref, ".out"))
  flat_profile <- flatProfile(pd)
  print(flat_profile)
  save(flat_profile, file= paste0("profiling/", file_pref, ".rda"))

  profileCallGraph2Dot(pd, score = "total", filename = paste0("profiling/", file_pref, ".dot"))
}

profile_walker(n_iter, n_thin, TRUE)

