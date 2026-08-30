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

test_that("count_achievable_configs(count_only = TRUE) matches the materialized count", {
  data <- make_toy_data()
  info_full  <- count_achievable_configs(data, list(x = c(7)))                # fixed only
  info_count <- count_achievable_configs(data, list(x = c(7)), count_only = TRUE)
  expect_equal(info_count$K, info_full$K)
  expect_null(info_count$per_variable)   # only $K/$counts_by_variable, no materialized configs
})

test_that("count_achievable_configs(count_only = TRUE) matches across covariates and multiple ranges", {
  data <- data.frame(x = c(1, 5, 10), y = c(1, 5, 10))
  specs <- list(x = list(c(0, 12), c(0, 3)), y = list(c(0, 12)))
  info_full  <- count_achievable_configs(data, specs)
  info_count <- count_achievable_configs(data, specs, count_only = TRUE)
  expect_equal(info_count$K, info_full$K)
  expect_equal(unname(info_count$counts_by_variable), unname(info_full$counts_by_variable))
})

test_that("count_achievable_configs(count_only = TRUE) sums (not multiplies) across a regime()'s own arms", {
  data <- make_toy_data()
  specs <- list(
    x = regime(
      low  = list(c(5, 6)),   # 1 achievable gap
      high = list(c(8, 9))    # 1 achievable gap
    )
  )
  info_full  <- count_achievable_configs(data, specs)
  info_count <- count_achievable_configs(data, specs, count_only = TRUE)
  expect_equal(info_full$K, 2L)      # 1 + 1, not 1 * 1 -- see test-regime.R
  expect_equal(info_count$K, 2)
})

test_that("run_M_draws(exact_if_K_leq =) falls back to Monte Carlo without materializing configs when K is too large", {
  ## A stand-in for the real-world case that motivated count_only: K here
  ## is small enough to actually check quickly in a test, but the point is
  ## the *code path* -- run_M_draws() must not materialize the full
  ## achievable-config list just to compare K against exact_if_K_leq when
  ## the comparison is going to fail anyway.
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))   # K = 9 (see the recovers-attributes test elsewhere)
  draws <- run_M_draws(data, "D", "Y", specs, M = 5, exact_if_K_leq = 1)
  expect_false(attr(draws, "exact"))
  expect_length(draws, 5)
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
