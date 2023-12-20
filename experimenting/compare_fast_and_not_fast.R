load_all()

all_data <- list(
  alces = alces, bison = bison, cervus = cervus, equus = equus, human = human, mammuthus = mammuthus)

ppo_fast <- list()
ppo_normal <- list()

for (species in c("alces", "bison")) {
  print (species)

  set.seed(14)
  ppo_fast[[species]] <- PPcalibrate(
    all_data[[species]]$c14_age, all_data[[species]]$c14_sig, intcal20, n_iter = 5e4, calendar_grid_resolution = 10, use_cpp = TRUE)

  set.seed(14)
  ppo_normal[[species]] <- PPcalibrate(
    all_data[[species]]$c14_age, all_data[[species]]$c14_sig, intcal20, n_iter = 5e4, calendar_grid_resolution = 10, use_cpp = FALSE)

  layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, ncol = 2)
  layout(mat = layout.matrix, heights = c(1, 1), widths = c(1, 1))

  PlotPosteriorMeanRate(ppo_fast[[species]])
  title(sub = paste("Fast version", species))
  PlotPosteriorChangePoints(ppo_fast[[species]], n_changes = c(1, 2, 3, 4))
  PlotPosteriorHeights(ppo_fast[[species]], n_changes = c(1, 2, 3, 4))
  PlotNumberOfInternalChanges(ppo_fast[[species]])

  PlotPosteriorMeanRate(ppo_normal[[species]])
  title(sub = paste("Original version", species))
  PlotPosteriorChangePoints(ppo_normal[[species]], n_changes = c(1, 2, 3, 4))
  PlotPosteriorHeights(ppo_normal[[species]], n_changes = c(1, 2, 3, 4))
  PlotNumberOfInternalChanges(ppo_normal[[species]])
}
