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

test_that("elicit_and_match's min_n_per_arm defaults to CEM's usual any-unit-per-arm rule", {
  ## min_n_per_arm = 1 (the default) must reproduce exactly what
  ## elicit_and_match() did before this parameter existed: retain any
  ## stratum with at least one unit in each arm.
  data <- make_toy_data()
  m1      <- elicit_and_match(data, "D", list(x = c(7)))
  m1_default <- elicit_and_match(data, "D", list(x = c(7)), min_n_per_arm = 1)
  expect_identical(m1$retained_idx, m1_default$retained_idx)
  expect_identical(m1$stratum, m1_default$stratum)
})

test_that("elicit_and_match's min_n_per_arm prunes strata with too few units in either arm", {
  ## Hand-built so the two strata have different per-arm sample sizes:
  ## x <= 7 has 3 treated/3 control, x > 7 has only 1 treated/1 control --
  ## raising min_n_per_arm should prune the thin stratum first, then both.
  data <- data.frame(
    x = c(1, 2, 3, 4, 5, 6, 11, 12),
    D = c(1, 1, 1, 0, 0, 0,  1,  0),
    Y = c(1, 2, 3, 4, 5, 6,  7,  8)
  )
  specs <- list(x = c(7))

  m1 <- elicit_and_match(data, "D", specs, min_n_per_arm = 1)
  expect_length(m1$retained_idx, 8)
  expect_equal(nlevels(m1$stratum), 2)

  m2 <- elicit_and_match(data, "D", specs, min_n_per_arm = 2)
  ## Only the x <= 7 stratum (3 treated/3 control) clears the bar; the
  ## x > 7 stratum (1 treated/1 control) is pruned.
  expect_setequal(m2$retained_idx, 1:6)
  expect_equal(nlevels(m2$stratum), 1)

  m4 <- elicit_and_match(data, "D", specs, min_n_per_arm = 4)
  ## Neither stratum has 4 per arm -- everything is pruned.
  expect_length(m4$retained_idx, 0)
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

test_that("fit_effect's vcov_type defaults to \"classical\" and reproduces the original (pre-vcov_type) behavior", {
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  fit_default   <- fit_effect(data, "D", "Y", m, estimator = "regression")
  fit_classical <- fit_effect(data, "D", "Y", m, estimator = "regression", vcov_type = "classical")
  expect_identical(fit_default$tau_hat, fit_classical$tau_hat)
  expect_identical(fit_default$var_hat, fit_classical$var_hat)
})

test_that("fit_effect's HC vcov_type changes var_hat but never tau_hat, and matches a direct sandwich/lmtest reference", {
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  fit_classical <- fit_effect(data, "D", "Y", m, estimator = "regression", vcov_type = "classical")
  fit_hc4       <- fit_effect(data, "D", "Y", m, estimator = "regression", vcov_type = "HC4")

  ## tau_hat is the same weighted-least-squares coefficient regardless of
  ## which variance estimator is applied to it afterward.
  expect_equal(fit_hc4$tau_hat, fit_classical$tau_hat)
  ## var_hat, on the other hand, should actually differ between the
  ## classical and HC4 variance estimators on this data (not a tautology --
  ## it's a real check that vcov_type is actually doing something).
  expect_false(isTRUE(all.equal(fit_hc4$var_hat, fit_classical$var_hat)))

  idx <- m$retained_idx
  w <- cem_weights(data$D[idx], m$stratum)
  ref <- stats::lm(Y ~ D + x, data = data[idx, ], weights = w)
  co  <- lmtest::coeftest(ref, vcov = sandwich::vcovHC(ref, type = "HC4"))
  expect_equal(fit_hc4$var_hat, unname(co["D", "Std. Error"])^2)
})

test_that("fit_effect's vcov_type is ignored for mean_diff", {
  data <- make_toy_data()
  m <- elicit_and_match(data, "D", list(x = c(7)))
  fit_classical <- fit_effect(data, "D", "Y", m, estimator = "mean_diff", vcov_type = "classical")
  fit_hc4       <- fit_effect(data, "D", "Y", m, estimator = "mean_diff", vcov_type = "HC4")
  expect_identical(fit_classical$tau_hat, fit_hc4$tau_hat)
  expect_identical(fit_classical$var_hat, fit_hc4$var_hat)
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

test_that("run_draw's vcov_type argument reaches fit_effect", {
  data <- make_toy_data()
  d_classical <- run_draw(data, "D", "Y", list(x = c(7)), vcov_type = "classical")
  d_hc4       <- run_draw(data, "D", "Y", list(x = c(7)), vcov_type = "HC4")
  expect_equal(d_classical$tau_hat, d_hc4$tau_hat)
  expect_false(isTRUE(all.equal(d_classical$var_hat, d_hc4$var_hat)))
})

test_that("run_M_draws's vcov_type argument reaches every draw and defaults to classical", {
  data <- make_toy_data()
  draws_default <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 2)
  draws_hc4     <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 2, vcov_type = "HC4")

  ## x = c(7) is a fixed spec, so every draw within a run is identical.
  expect_equal(attr(draws_default, "vcov_type"), "classical")
  expect_equal(attr(draws_hc4, "vcov_type"), "HC4")
  expect_equal(draws_default[[1]]$tau_hat, draws_hc4[[1]]$tau_hat)
  expect_false(isTRUE(all.equal(draws_default[[1]]$var_hat, draws_hc4[[1]]$var_hat)))
})

test_that("run_M_draws carries min_n_per_arm as an attribute, defaulting to 1", {
  data <- make_toy_data()
  draws_default <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 2)
  draws_strict  <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 2, min_n_per_arm = 2)
  expect_equal(attr(draws_default, "min_n_per_arm"), 1)
  expect_equal(attr(draws_strict, "min_n_per_arm"), 2)
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
