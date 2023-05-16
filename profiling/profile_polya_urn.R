library(proftools)

set.seed(14)

n_iter = 100000
n_thin = 10

profile_polya_urn <- function(n_iter, n_thin, use_cpp) {
  file_pref = paste("PU", n_iter, n_thin, "cpp", use_cpp, sep="_")

  Rprof(paste("profiling/", file_pref, ".out", sep = ""))
  walker_output = PolyaUrnBivarDirichlet(
    kerr$c14_age,
    kerr$c14_sig,
    intcal20,
    n_iter,
    n_thin,
    use_cpp = use_cpp)
  Rprof(NULL)

  pd = readProfileData(paste("profiling/", file_pref, ".out", sep = ""))
  flat_profile = flatProfile(pd)
  print(flat_profile)
  save(flat_profile, file=paste("profiling/", file_pref, ".rda", sep = ""))

  profileCallGraph2Dot(pd, score = "total", filename = paste("profiling/", file_pref, ".dot", sep = ""))
}

profile_polya_urn(n_iter, n_thin, TRUE)
profile_polya_urn(n_iter, n_thin, FALSE)
