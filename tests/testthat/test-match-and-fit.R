test_that("elicit_and_match retains all rows when every stratum has both arms", {
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  expect_length(m$retained_idx, 8)
  expect_equal(nlevels(m$stratum), 2)
})

test_that("elicit_and_match prunes strata lacking common support", {
  data <- make_toy_data_with_unmatched()
  m <- elicit_and_match(data, "D", list(x = c(7, 20)))
  expect_setequal(m$retained_idx, 1:8)
  expect_false(9 %in% m$retained_idx)
  expect_false(10 %in% m$retained_idx)
})

test_that("elicit_and_match reports kinds and NA regimes for non-regime covariates", {
  data <- make_toy_data()
  data$g <- rep(c("a", "b"), 4)
  m <- elicit_and_match(data, "D", list(x = c(7), g = "exact"))
  expect_equal(unname(m$kinds[c("x", "g")]), c("fixed", "exact"))
  expect_true(all(is.na(m$regimes)))
})

test_that("elicit_and_match falls back to one stratum when nothing coarsens", {
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = NULL))
  expect_equal(nlevels(m$stratum), 1)
  expect_length(m$retained_idx, 8)
})

test_that("fit_effect's mean_diff estimator recovers the known pooled ATT and variance", {
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  fit <- fit_effect(data, "D", "Y", m, estimator = "mean_diff")
  expect_equal(fit$tau_hat, 5.5)
  expect_equal(fit$var_hat, 0.625)
  expect_equal(fit$n_used, 8)
})

test_that("fit_effect's regression estimator is the default and differs from mean_diff here", {
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  fit_default <- fit_effect(data, "D", "Y", m)
  fit_reg     <- fit_effect(data, "D", "Y", m, estimator = "regression")
  expect_equal(fit_default$tau_hat, fit_reg$tau_hat)
  ## Regression additionally adjusts for x's raw (uncoarsened) value, which
  ## has real explanatory power here, so it shouldn't reproduce mean_diff's
  ## number on this toy data.
  expect_false(isTRUE(all.equal(fit_reg$tau_hat, 5.5)))
})

test_that("fit_effect's regression estimator matches a direct weighted lm() on the matched sample", {
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  fit <- fit_effect(data, "D", "Y", m, estimator = "regression")

  idx <- m$retained_idx
  w <- cem_weights(data$D[idx], m$stratum)
  ref <- stats::lm(Y ~ D + x, data = data[idx, ], weights = w)
  co <- summary(ref)$coefficients

  expect_equal(fit$tau_hat, unname(co["D", "Estimate"]))
  expect_equal(fit$var_hat, unname(co["D", "Std. Error"])^2)
  expect_equal(fit$n_used, 8)
})

test_that("fit_effect's regression estimator matches an independently hand-computed value", {
  ## The two strata here are perfectly balanced (2 treated, 2 control
  ## each), so cem_weights() is 1 for every unit -- this reduces to a
  ## plain (unweighted) lm(Y ~ D + x), independently verified via numpy.
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  fit <- fit_effect(data, "D", "Y", m, estimator = "regression")
  expect_equal(fit$tau_hat, 7.50990099009901, tolerance = 1e-6)
  expect_equal(fit$var_hat, 0.10344574061366504, tolerance = 1e-6)
})

test_that("fit_effect's regression estimator defaults covariates to the non-excluded spec entries", {
  data <- make_toy_data()
  data$z <- c(1, 1, 1, 1, 2, 2, 2, 2)
  m <- elicit_and_match(data, "D", list(x = c(7), z = NULL))  # z excluded
  fit_default <- fit_effect(data, "D", "Y", m, estimator = "regression")

  idx <- m$retained_idx
  w <- cem_weights(data$D[idx], m$stratum)
  ref <- stats::lm(Y ~ D + x, data = data[idx, ], weights = w)  # no z
  expect_equal(fit_default$tau_hat, unname(coef(ref)["D"]))
})

test_that("fit_effect's regression estimator accepts an explicit covariates override", {
  data <- make_toy_data()
  data$z <- c(1, 1, 1, 1, 2, 2, 2, 2)
  m <- elicit_and_match(data, "D", list(x = c(7), z = NULL))
  fit_z <- fit_effect(data, "D", "Y", m, estimator = "regression", covariates = c("x", "z"))

  idx <- m$retained_idx
  w <- cem_weights(data$D[idx], m$stratum)
  ref <- stats::lm(Y ~ D + x + z, data = data[idx, ], weights = w)
  expect_equal(fit_z$tau_hat, unname(coef(ref)["D"]))
})

test_that("fit_effect's regression estimator still computes unit_tau_hat via stratum mean-differences", {
  ## The covariance diagnostics need a covariate-indexed tau(x) proxy,
  ## which a single regression coefficient per draw can't give -- so this
  ## should be unaffected by estimator choice.
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  fit_md  <- fit_effect(data, "D", "Y", m, estimator = "mean_diff")
  fit_reg <- fit_effect(data, "D", "Y", m, estimator = "regression")
  expect_equal(fit_reg$unit_tau_hat, fit_md$unit_tau_hat)
  expect_equal(fit_reg$strat_diff, fit_md$strat_diff)
})

test_that("fit_effect's regression estimator returns NA gracefully when treat_var is aliased", {
  data <- make_toy_data()
  data$D <- 1  # no within-sample variation in treatment
  m <- list(retained_idx = 1:8, stratum = factor(rep("s1", 8)), kinds = c(x = "fixed"))
  fit <- fit_effect(data, "D", "Y", m, estimator = "regression", covariates = "x")
  expect_true(is.na(fit$tau_hat))
  expect_true(is.na(fit$var_hat))
})

test_that("fit_effect returns NA/empty results when nothing is retained", {
  data <- make_toy_data()
  m <- list(retained_idx = integer(0), stratum = factor(character(0)))
  fit <- fit_effect(data, "D", "Y", m)
  expect_true(is.na(fit$tau_hat))
  expect_equal(fit$n_used, 0)
})

test_that("run_draw composes elicit_and_match and fit_effect", {
  data <- make_toy_data()
  d <- run_draw(data, "D", "Y", list(x = c(7)), estimator = "mean_diff")
  expect_equal(d$tau_hat, 5.5)
  expect_equal(d$matched$kinds[["x"]], "fixed")
})

test_that("run_draw defaults to the regression estimator", {
  data <- make_toy_data()
  d <- run_draw(data, "D", "Y", list(x = c(7)))
  expect_equal(d$tau_hat, 7.50990099009901, tolerance = 1e-6)
})

test_that("run_M_draws returns M equal-weight Monte Carlo draws by default", {
  data <- make_toy_data()
  draws <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 5)
  expect_length(draws, 5)
  expect_false(attr(draws, "exact"))
  expect_equal(attr(draws, "weights"), rep(0.2, 5))
  ## x = c(7) is a fixed spec, so every draw is identical -- true
  ## regardless of which estimator produced the (unspecified) number.
  expect_length(unique(vapply(draws, function(d) d$tau_hat, numeric(1))), 1)
})

test_that("run_M_draws's estimator/covariates arguments reach fit_effect", {
  data <- make_toy_data()
  draws_md  <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 2, estimator = "mean_diff")
  draws_reg <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 2, estimator = "regression")
  expect_true(all(vapply(draws_md,  function(d) d$tau_hat, numeric(1)) == 5.5))
  expect_true(all(vapply(draws_reg, function(d) d$tau_hat, numeric(1)) != 5.5))
})

test_that("run_M_draws switches to exact enumeration when K is small enough", {
  data <- make_toy_data()
  ## No observed x values fall in (5, 6), so there is exactly one achievable
  ## cutpoint here regardless of M.
  draws <- run_M_draws(data, "D", "Y", list(x = list(c(5, 6))), M = 1,
                        exact_if_K_leq = 1000, estimator = "mean_diff")
  expect_true(attr(draws, "exact"))
  expect_length(draws, 1)
  expect_equal(draws[[1]]$tau_hat, 5.5)
})

test_that("run_M_draws falls back to Monte Carlo when K exceeds the threshold", {
  data <- make_toy_data()
  draws <- run_M_draws(data, "D", "Y", list(x = list(c(5, 6))), M = 3,
                        exact_if_K_leq = 0)
  expect_false(attr(draws, "exact"))
  expect_length(draws, 3)
})

test_that("run_M_draws carries treat_var/outcome_var/cutpoint_specs/estimator/covariates as attributes", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws_mc <- run_M_draws(data, "D", "Y", specs, M = 3)
  expect_equal(attr(draws_mc, "treat_var"), "D")
  expect_equal(attr(draws_mc, "outcome_var"), "Y")
  expect_equal(attr(draws_mc, "cutpoint_specs"), specs)
  expect_equal(attr(draws_mc, "estimator"), "regression")  # default
  expect_null(attr(draws_mc, "covariates"))                # default: auto-derive per draw

  draws_md <- run_M_draws(data, "D", "Y", specs, M = 3, estimator = "mean_diff", covariates = "x")
  expect_equal(attr(draws_md, "estimator"), "mean_diff")
  expect_equal(attr(draws_md, "covariates"), "x")

  draws_exact <- run_M_draws(data, "D", "Y", list(x = list(c(5, 6))), M = 1,
                              exact_if_K_leq = 1000)
  expect_equal(attr(draws_exact, "treat_var"), "D")
  expect_equal(attr(draws_exact, "outcome_var"), "Y")
  expect_equal(attr(draws_exact, "cutpoint_specs"), list(x = list(c(5, 6))))
  expect_equal(attr(draws_exact, "estimator"), "regression")
})
