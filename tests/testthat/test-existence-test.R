test_that("existence_test(data, draws) recovers treat_var/outcome_var/cutpoint_specs from draws", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  ex <- existence_test(data, draws, n_boot = 5)

  expect_s3_class(ex, "ecem_existence_test")
  expect_true(is.numeric(ex$observed_stat))
  expect_true(is.numeric(ex$p_value))
  expect_length(ex$null_stats, 5)
  expect_equal(ex$M, 3)
  expect_equal(ex$n_boot, 5)
})

test_that("existence_test still accepts explicit treat_var/outcome_var/cutpoint_specs", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  ex <- existence_test(data, draws, treat_var = "D", outcome_var = "Y",
                        cutpoint_specs = specs, n_boot = 5)
  expect_s3_class(ex, "ecem_existence_test")
})

test_that("existence_test errors informatively if draws has no recoverable attributes", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  attr(draws, "treat_var") <- NULL
  attr(draws, "outcome_var") <- NULL
  attr(draws, "cutpoint_specs") <- NULL
  expect_error(existence_test(data, draws, n_boot = 5), "could not all be recovered")
})

test_that("existence_test labels a plain stat_fun and falls back for anonymous ones", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)

  ex1 <- existence_test(data, draws, n_boot = 5)
  expect_equal(ex1$stat_label, "median")

  ex2 <- existence_test(data, draws, n_boot = 5, stat_fun = mean)
  expect_equal(ex2$stat_label, "mean")

  ex3 <- existence_test(data, draws, n_boot = 5, stat_fun = function(v, na.rm = TRUE) stats::median(v, na.rm = na.rm))
  expect_equal(ex3$stat_label, "custom")
})

test_that("existence_test(diagnostics =) reuses pooling_diagnostics()'s cached bootstrap", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  diag <- pooling_diagnostics(data, draws, n_boot = 5)

  expect_type(diag$existence_boot, "list")
  expect_equal(dim(diag$existence_boot$tau_own), c(5, 3))
  expect_equal(diag$existence_boot$tau_hat_m, vapply(draws, function(d) d$tau_hat, numeric(1)))

  ex <- existence_test(data, draws, diagnostics = diag)
  expect_true(ex$reused_bootstrap)
  expect_equal(ex$n_boot, 5)
  expect_length(ex$null_stats, 5)
  ## No bootstrapping was rerun, so this should match a manual summary of
  ## the exact same cached matrix.
  expect_equal(ex$null_stats, apply(diag$existence_boot$tau_own, 1, stats::median, na.rm = TRUE))
})

test_that("existence_test(diagnostics =) lets stat_fun vary without rerunning the bootstrap", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  diag <- pooling_diagnostics(data, draws, n_boot = 5)

  ex_mean <- existence_test(data, draws, diagnostics = diag, stat_fun = mean)
  expect_equal(ex_mean$null_stats, apply(diag$existence_boot$tau_own, 1, mean, na.rm = TRUE))
})

test_that("existence_test errors if diagnostics doesn't match draws", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  diag <- pooling_diagnostics(data, draws, n_boot = 5)

  ## Force a deterministic mismatch rather than relying on two independent
  ## draws happening to differ (which a fixed-only spec like this one, with
  ## no randomness in the matching, would not guarantee).
  draws_other <- draws
  draws_other[[1]]$tau_hat <- draws_other[[1]]$tau_hat + 1
  expect_error(existence_test(data, draws_other, diagnostics = diag), "different set")
})

test_that("existence_test errors if diagnostics isn't an ecem_pooling_diagnostics object", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  expect_error(existence_test(data, draws, diagnostics = list()), "must be the result of pooling_diagnostics")
})

## print.ecem_existence_test's branching logic, tested against hand-built
## objects with known p-values -- deterministic and independent of n_boot.
make_fake_existence <- function(observed_stat = 2, p_value = 0.5, M = 10, n_boot = 50,
                                 stat_label = "median") {
  structure(
    list(
      observed_stat = observed_stat,
      p_value       = p_value,
      null_stats    = numeric(0),
      M             = M,
      n_boot        = n_boot,
      stat_label    = stat_label
    ),
    class = "ecem_existence_test"
  )
}

test_that("print.ecem_existence_test: rejects -> existence evidence, no pooled tau_bar", {
  x <- make_fake_existence(p_value = 0.01)
  expect_output(print(x), "REJECT null of no effect")
  expect_output(print(x), "some nonzero effect exists")
  expect_output(print(x), "pooled tau_bar")
})

test_that("print.ecem_existence_test: fails to reject -> no evidence", {
  x <- make_fake_existence(p_value = 0.5)
  expect_output(print(x), "fail to reject")
  expect_output(print(x), "no evidence here")
})

test_that("print.ecem_existence_test respects the alpha argument", {
  x <- make_fake_existence(p_value = 0.03)
  expect_output(print(x, alpha = 0.01), "fail to reject")
  expect_output(print(x, alpha = 0.05), "REJECT null of no effect")
})

test_that("print.ecem_existence_test returns x invisibly", {
  x <- make_fake_existence()
  expect_invisible(print(x))
})
