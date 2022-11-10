file_pref = "Walker_before_changes"

set.seed(14)
n_iter = 1000
n_thin = 10
Rprof(paste("profiling/", file_pref, ".out", sep = ""))
walker_output = WalkerBivarDirichlet(kerr$c14_ages, kerr$c14_sig, intcal20, n_iter, n_thin)
Rprof(NULL)
pd = readProfileData(paste("profiling/", file_pref, ".out", sep = ""))
flat_profile = flatProfile(pd)
print(flat_profile)
save(flat_profile, file=paste("profiling/", file_pref, ".rda", sep = ""))

profileCallGraph2Dot(pd, score = "total", filename = paste("profiling/", file_pref, ".dot", sep = ""))

