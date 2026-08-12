test_that("label_diagnostics(data, draws) recovers treat_var/cutpoint_specs from draws", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  lab <- label_diagnostics(data, draws)

  expect_s3_class(lab, "ecem_label_diagnostics")
  expect_true(is.numeric(lab$cov_tau_p))
  expect_equal(lab$att_ate_gap, lab$cov_tau_p / mean(data$D == 1))
  expect_length(lab$gap_m, 5)
  expect_equal(lab$se, sqrt(pool_draws(draws)$T))
  expect_equal(lab$covariates, "x")
})

test_that("label_diagnostics still accepts explicit treat_var/cutpoint_specs", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  lab <- label_diagnostics(data, draws, treat_var = "D", cutpoint_specs = specs)
  expect_s3_class(lab, "ecem_label_diagnostics")
})

test_that("label_diagnostics errors informatively if draws has no recoverable attributes", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  attr(draws, "treat_var") <- NULL
  attr(draws, "cutpoint_specs") <- NULL
  expect_error(label_diagnostics(data, draws), "could not both be recovered")
})

test_that("label_diagnostics defaults covariates to the non-excluded spec entries", {
  data <- make_toy_data()
  data$z <- c(1, 1, 1, 1, 2, 2, 2, 2)
  specs <- list(x = c(7), z = NULL)  # z excluded from matching
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  lab <- label_diagnostics(data, draws)
  expect_equal(lab$covariates, "x")
})

test_that("label_diagnostics accepts an explicit, larger covariates set", {
  data <- make_toy_data()
  data$z <- c(1, 1, 1, 1, 2, 2, 2, 2)
  specs <- list(x = c(7), z = NULL)
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  ## x and this hand-constructed z jointly separate treatment status
  ## perfectly on this tiny, deterministic toy sample -- fit_propensity()'s
  ## glm() correctly warns about it (fitted probabilities numerically 0 or
  ## 1), and that warning is left untouched at the source deliberately: on
  ## real data, it's a legitimate signal that the propensity model has
  ## degenerated and p_hat(X) near 0/1 could destabilize the ATT-ATE
  ## covariance term downstream, so fit_propensity() shouldn't be modified
  ## to hide it just because this particular test triggers it too. Assert
  ## on the warning explicitly instead of suppressing it, so if a future
  ## change to make_toy_data()/z stops triggering separation (or starts
  ## triggering a different warning), that shows up as a test failure
  ## rather than silently passing either way.
  expect_warning(
    lab <- label_diagnostics(data, draws, covariates = c("x", "z")),
    "fitted probabilities numerically 0 or 1 occurred"
  )
  expect_equal(lab$covariates, c("x", "z"))
})

## print.ecem_label_diagnostics's branching logic, tested against hand-built
## objects -- deterministic and independent of any actual bootstrap or fit.
make_fake_label <- function(att_ate_gap = 0.1, gap_r_bar = 0.1, se = 1,
                             cov_tau_p = att_ate_gap / 2, gap_m = c(0.1, -0.1),
                             gap_r_sd = 0.05, appreciable_frac_se = 0.5) {
  structure(
    list(
      cov_tau_p = cov_tau_p, att_ate_gap = att_ate_gap, gap_m = gap_m,
      gap_r_bar = gap_r_bar, gap_r_sd = gap_r_sd, se = se,
      pooled = NULL, covariates = "x", appreciable_frac_se = appreciable_frac_se
    ),
    class = "ecem_label_diagnostics"
  )
}

test_that("print.ecem_label_diagnostics: both terms appreciable -> FSATT", {
  x <- make_fake_label(att_ate_gap = 1, gap_r_bar = 1, se = 1)
  expect_output(print(x), "estimate of: FSATT")
})

test_that("print.ecem_label_diagnostics: retention appreciable, ATT-ATE negligible -> still FSATT", {
  ## The one combination the paper doesn't spell out -- see
  ## label_diagnostics()'s Details for why the chain stops at FSATT here.
  x <- make_fake_label(att_ate_gap = 0.01, gap_r_bar = 1, se = 1)
  expect_output(print(x), "estimate of: FSATT")
  expect_output(print(x), "promoted past FSATT")
})

test_that("print.ecem_label_diagnostics: retention negligible, ATT-ATE appreciable -> ATT", {
  x <- make_fake_label(att_ate_gap = 1, gap_r_bar = 0.01, se = 1)
  expect_output(print(x), "estimate of: ATT")
  expect_output(print(x), "most specific label")
})

test_that("print.ecem_label_diagnostics: both terms negligible -> ATE", {
  x <- make_fake_label(att_ate_gap = 0.01, gap_r_bar = 0.01, se = 1)
  expect_output(print(x), "estimate of: ATE")
  expect_output(print(x), "FSATT ~ ATT ~ ATE")
})

test_that("print.ecem_label_diagnostics treats an undefined SE conservatively as appreciable", {
  x <- make_fake_label(att_ate_gap = 0.01, gap_r_bar = 0.01, se = NA_real_)
  expect_output(print(x), "= undefined")
  expect_output(print(x), "estimate of: FSATT")
})

test_that("print.ecem_label_diagnostics returns x invisibly", {
  x <- make_fake_label()
  expect_invisible(print(x))
})
