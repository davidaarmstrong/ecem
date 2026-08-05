test_that("enumerate_gaps produces k+1 gaps for k distinct interior observed values", {
  x <- c(1, 5, 10, 15, 20)  # 5 distinct values, none at the range endpoints
  g <- enumerate_gaps(x, c(0, 25))
  expect_equal(nrow(g), 6)              # k + 1 = 6
  expect_equal(sum(g$width), 25)        # widths sum to the range's total width
})

test_that("enumerate_gaps collapses to one gap when no observed values fall inside", {
  x <- c(1, 2, 3, 4, 11, 12, 13, 14)
  g <- enumerate_gaps(x, c(5, 6))
  expect_equal(nrow(g), 1)
  expect_equal(g$rep, 5.5)
  expect_equal(g$width, 1)
})

test_that("count_achievable_configs returns K = 1 when nothing varies", {
  data <- make_toy_data()
  info <- count_achievable_configs(data, list(x = c(7)))  # fixed only
  expect_equal(info$K, 1L)
  expect_length(info$per_variable, 0)
})

test_that("count_achievable_configs multiplies across covariates", {
  data <- data.frame(x = c(1, 5, 10), y = c(1, 5, 10))
  specs <- list(x = list(c(0, 12)), y = list(c(0, 12)))
  info <- count_achievable_configs(data, specs)
  kx <- count_achievable_configs(data, list(x = specs$x))$K
  ky <- count_achievable_configs(data, list(y = specs$y))$K
  expect_equal(info$K, kx * ky)
})

test_that("enumerate_configs produces K specs with weights summing to 1", {
  data <- data.frame(x = c(1, 5, 10))
  specs <- list(x = list(c(0, 12)))
  enum <- enumerate_configs(data, specs)
  expect_equal(length(enum$specs), count_achievable_configs(data, specs)$K)
  expect_equal(sum(enum$weights), 1)
})

test_that("run_enumerated_draws / pool_rubins_rules_exact recovers the known ATT", {
  data <- make_toy_data()
  ## No observed x values fall in (5, 10), so there is exactly one
  ## achievable cutpoint (the midpoint, 7.5) -- the same clean two-strata
  ## split as the fixed-cutpoint = 7 case tested elsewhere.
  specs <- list(x = list(c(5, 10)))
  expect_equal(count_achievable_configs(data, specs)$K, 1L)

  ## estimator = "mean_diff" here since this test is about the enumeration/
  ## pooling mechanics (one achievable config -> B = 0), not about which
  ## estimator fit_effect() uses -- see test-match-and-fit.R for that.
  enum_draws <- run_enumerated_draws(data, "D", "Y", specs, estimator = "mean_diff")
  pooled_exact <- pool_rubins_rules_exact(enum_draws$results, enum_draws$weights)
  expect_equal(pooled_exact$tau_bar, 5.5)
  expect_equal(pooled_exact$B, 0)   # only one achievable config
})
