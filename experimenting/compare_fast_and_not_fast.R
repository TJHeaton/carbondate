all_data <- list(
  alces = alces, bison = bison, cervus = cervus, equus = equus, human = human, mammuthus = mammuthus)

ppo_fast <- list()
ppo_normal <- list()
ppo_cpp <- list()

for (species in names(all_data)) {
  print (species)

  set.seed(14)
  start = Sys.time()
  ppo_normal[[species]] <- PPcalibrate(
    all_data[[species]]$c14_age, all_data[[species]]$c14_sig, intcal20, n_iter = 1e4, calendar_grid_resolution = 10, use_fast = FALSE)
  end = Sys.time()
  print(round(end - start), 2)

  set.seed(14)
  start = Sys.time()
  ppo_cpp[[species]] <- PPcalibrate(
    all_data[[species]]$c14_age, all_data[[species]]$c14_sig, intcal20, n_iter = 1e4, calendar_grid_resolution = 10, use_fast = TRUE)
  end = Sys.time()
  print(round(end - start), 2)

  layout.matrix <- matrix(1:8, nrow = 4, ncol = 2)
  layout(mat = layout.matrix, heights = c(1, 1), widths = c(1, 1))

  set.seed(10)
  PlotPosteriorMeanRate(ppo_normal[[species]])
  title(sub = paste("Original version", species))
  PlotPosteriorChangePoints(ppo_normal[[species]], n_changes = c(3, 4, 5))
  PlotPosteriorHeights(ppo_normal[[species]], n_changes = c(3, 4, 5))
  PlotNumberOfInternalChanges(ppo_normal[[species]])

  set.seed(10)
  PlotPosteriorMeanRate(ppo_cpp[[species]])
  title(sub = paste("CPP version", species))
  PlotPosteriorChangePoints(ppo_cpp[[species]], n_changes = c(3, 4, 5))
  PlotPosteriorHeights(ppo_cpp[[species]], n_changes = c(3, 4, 5))
  PlotNumberOfInternalChanges(ppo_cpp[[species]])

}
