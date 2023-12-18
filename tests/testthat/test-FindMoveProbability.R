test_that("FindMoveProbability gives expected outcome", {

  prior_n_internal_changepoints_lambda <- 60
  k_max_internal_changepoints <- 300

  prob_move <- .FindMoveProbability(
    prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints = k_max_internal_changepoints,
    rescale_factor = 0.9)

  expect_true(all(prob_move$pos >= 0))
  expect_true(all(prob_move$height >= 0))
  expect_true(all(prob_move$birth >= 0))
  expect_true(all(prob_move$death >= 0))

  expect_length(prob_move$pos, k_max_internal_changepoints + 1)

  expect_true(all(prob_move$pos + prob_move$height + prob_move$birth + prob_move$death == 1))

})


test_that("FindMoveProbability gives same as legacy code", {
  prior_n_internal_changepoints_lambda <- 60
  k_max_internal_changepoints <- 300

  revised_prob_move <- .FindMoveProbability(
    prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints = k_max_internal_changepoints,
    rescale_factor = 0.9)

  source(test_path("fixtures", "LegacyFindMoveProb.R"))
  legacy_prob_move <- LegacyFindMoveProb(
    lambda = prior_n_internal_changepoints_lambda,
    kmax = k_max_internal_changepoints)

  expect_identical(revised_prob_move$pos,
                   legacy_prob_move$Pos)
  expect_identical(revised_prob_move$height,
                   legacy_prob_move$He)
  expect_identical(revised_prob_move$birth,
                   legacy_prob_move$Birth)
  expect_identical(revised_prob_move$death,
                   legacy_prob_move$Death)

})
