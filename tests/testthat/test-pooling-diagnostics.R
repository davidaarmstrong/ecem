test_that("elicited_union_range covers all cutpoint kinds correctly", {
  expect_null(elicited_union_range(NULL))
  expect_null(elicited_union_range("exact"))
  expect_null(elicited_union_range(c(12, 16)))
  expect_equal(elicited_union_range(list(c(25, 33), c(60, 68))), c(25, 68))

  r <- regime(low = list(c(5, 6)), high = list(c(8, 9)))
  expect_equal(elicited_union_range(r), c(5, 9))
})

test_that("pooling_diagnostics(data, draws) recovers treat_var/outcome_var/cutpoint_specs from draws", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, n_boot = 10)

  expect_s3_class(diag, "ecem_pooling_diagnostics")
  expect_named(diag$flatness, "x")
  expect_true(is.numeric(diag$excess_variance$p_value))
  expect_false(is.null(diag$retention))
  expect_length(diag$retention$gap, 5)
  ## 8 distinct x values fall in [0, 20], neither endpoint coinciding with
  ## one of them, so K = 8 + 1 = 9 (see count_achievable_configs()).
  expect_equal(diag$K, 9L)
})

test_that("pooling_diagnostics still accepts explicit treat_var/outcome_var/cutpoint_specs", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, treat_var = "D", outcome_var = "Y",
                               cutpoint_specs = specs, n_boot = 10)
  expect_named(diag$flatness, "x")
})

test_that("pooling_diagnostics errors informatively if draws has no recoverable attributes", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  attr(draws, "treat_var") <- NULL
  attr(draws, "outcome_var") <- NULL
  attr(draws, "cutpoint_specs") <- NULL
  expect_error(pooling_diagnostics(data, draws, n_boot = 10), "could not all be recovered")
})

test_that("pooling_diagnostics reports no elicited covariates when the spec is entirely fixed", {
  data <- make_toy_data()
  specs <- list(x = c(7))  # fixed only, nothing elicited
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, n_boot = 10)
  expect_length(diag$flatness, 0)
})

test_that("pooling_diagnostics skips the retention diagnostic when asked to", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, n_boot = 10, run_retention = FALSE)
  expect_null(diag$retention)
})

test_that("pooling_diagnostics computes pooled itself when not supplied", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, n_boot = 10)
  expect_equal(diag$pooled, pool_draws(draws))
})

## print.ecem_pooling_diagnostics's branching logic, tested against
## hand-built objects with known p-values rather than relying on the
## stochastic tests above -- deterministic and independent of n_boot.
make_fake_diag <- function(flat_p, ev_p, ev_ratio = 1,
                            retention = list(gap = c(0.01, -0.02), mean = -0.005, sd = 0.02),
                            exact = FALSE, M = 10, K = NULL, alpha = 0.05) {
  structure(
    list(
      flatness        = list(age = list(p_value = flat_p, model = NULL, n = 100)),
      excess_variance = list(p_value = ev_p, ratio = ev_ratio, B_null_dist = numeric(0)),
      retention       = retention,
      pooled          = NULL,
      M               = M,
      K               = K,
      exact           = exact,
      alpha           = alpha
    ),
    class = "ecem_pooling_diagnostics"
  )
}

test_that("print.ecem_pooling_diagnostics: neither test rejects -> pool cleanly", {
  d <- make_fake_diag(flat_p = 0.5, ev_p = 0.5)
  expect_output(print(d), "pool by Rubin's rules with no caveat")
})

test_that("print.ecem_pooling_diagnostics: excess-variance rejects -> existence_test guidance", {
  d <- make_fake_diag(flat_p = 0.5, ev_p = 0.01)
  expect_output(print(d), "EXCESS VARIANCE DETECTED")
  expect_output(print(d), "existence_test\\(\\)")
  expect_output(print(d), "mechanism-locator")
})

test_that("print.ecem_pooling_diagnostics: flatness rejects, excess-variance doesn't -> flag discrepancy", {
  d <- make_fake_diag(flat_p = 0.01, ev_p = 0.5)
  expect_output(print(d), "HETEROGENEOUS")
  expect_output(print(d), "flag the discrepancy")
  expect_output(print(d), "complete explanation")
})

test_that("print.ecem_pooling_diagnostics adapts guidance when retention wasn't computed", {
  d <- make_fake_diag(flat_p = 0.01, ev_p = 0.5, retention = NULL)
  expect_output(print(d), "Rerun with run_retention = TRUE")
})

## The "retention not small -> what next" sub-branch: exact enumeration if
## K is within reach of M, a larger M if not, and neither if draws were
## already exact -- see power_next_step().
test_that("print.ecem_pooling_diagnostics recommends exact enumeration when K is within reach of M", {
  d <- make_fake_diag(flat_p = 0.01, ev_p = 0.5, M = 10, K = 50, exact = FALSE)
  expect_output(print(d), "exact_if_K_leq = 50")
  expect_output(print(d), "rather than a larger M")
})

test_that("print.ecem_pooling_diagnostics recommends a larger M when K is far beyond reach", {
  d <- make_fake_diag(flat_p = 0.01, ev_p = 0.5, M = 10, K = 1e7, exact = FALSE)
  expect_output(print(d), "far more than")
  expect_output(print(d), "more realistic")
})

test_that("print.ecem_pooling_diagnostics says it isn't a power problem when already exact", {
  d <- make_fake_diag(flat_p = 0.01, ev_p = 0.5, M = 50, K = 50, exact = TRUE)
  expect_output(print(d), "isn't an M/K power problem")
  expect_output(print(d), "real scope")
})

test_that("print.ecem_pooling_diagnostics falls back gracefully when K wasn't recorded", {
  d <- make_fake_diag(flat_p = 0.01, ev_p = 0.5, M = 10, K = NULL, exact = FALSE)
  expect_output(print(d), "Consider exact enumeration")
})

test_that("print.ecem_pooling_diagnostics returns x invisibly", {
  d <- make_fake_diag(flat_p = 0.5, ev_p = 0.5)
  expect_invisible(print(d))
})
