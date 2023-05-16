walker_outputs = list()

for (i in 1:10) {
  set.seed(i)
  walker_outputs[[i]] = WalkerBivarDirichlet(kerr$c14_age, kerr$c14_sig, FALSE, intcal20, n_iter=1e5, n_thin=10)
  walker_outputs[[i]]$label = paste("Seed", i)
}


for (i in 2:10) {
  for (j in 1:(i-1)) {
    .KullbackLeiblerDivergenceFinal(walker_outputs[[i]], walker_outputs[[j]])
  }
}
