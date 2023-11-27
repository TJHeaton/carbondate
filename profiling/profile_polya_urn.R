library(proftools)

set.seed(14)

n_iter <- 100000
n_thin <- 10

profile_polya_urn <- function(n_iter, n_thin, use_cpp) {
  file_pref <- paste("PU", n_iter, n_thin, "cpp", use_cpp, sep="_")

  Rprof(paste0("profiling/", file_pref, ".out"))
  PolyaUrnBivarDirichlet(
    kerr$c14_age,
    kerr$c14_sig,
    intcal20,
    n_iter,
    n_thin,
    use_cpp = use_cpp)
  Rprof(NULL)

  pd <- readProfileData(paste0("profiling/", file_pref, ".out"))
  flat_profile <- flatProfile(pd)
  print(flat_profile)
  save(flat_profile, file= paste0("profiling/", file_pref, ".rda"))

  profileCallGraph2Dot(pd, score = "total", filename = paste0("profiling/", file_pref, ".dot"))
}

profile_polya_urn(n_iter, n_thin, TRUE)
profile_polya_urn(n_iter, n_thin, FALSE)
