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
  ## Nothing in this test reads $existence_boot -- skip that bootstrap
  ## (run_existence_cache = TRUE is otherwise unconditional regardless of
  ## whether anything downstream needs it, and on this tiny N = 8 toy
  ## dataset a resample can easily land on a near-zero-residual fit and
  ## throw summary.lm()'s harmless-but-noisy "essentially perfect fit"
  ## warning).
  diag <- pooling_diagnostics(data, draws, n_boot = 10, run_existence_cache = FALSE)

  expect_s3_class(diag, "ecem_pooling_diagnostics")
  expect_named(diag$flatness, "x")
  ## congeniality_test() is run by default now; on this tiny (N = 8) toy
  ## dataset the matched interaction model may or may not have enough
  ## residual df to identify -- either a numeric p-value or a graceful NA
  ## (see congeniality_test()'s degeneracy guard) is an acceptable result
  ## here, the point of this test is that pooling_diagnostics() itself
  ## doesn't error.
  expect_false(is.null(diag$congeniality))
  expect_true(is.numeric(diag$congeniality$p_value))
  ## excess_variance is no longer computed by default (superseded design).
  expect_null(diag$excess_variance)
  expect_false(is.null(diag$retention))
  expect_length(diag$retention$gap, 5)
  ## 8 distinct x values fall in [0, 20], neither endpoint coinciding with
  ## one of them, so K = 8 + 1 = 9 (see count_achievable_configs()).
  expect_equal(diag$K, 9L)
})

test_that("pooling_diagnostics(run_excess_variance = TRUE) still reproduces the superseded design", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  ## This test's whole point is the bootstrap-based excess-variance design,
  ## so its own "essentially perfect fit" warnings (see the note in the
  ## first test above) are an inherent, expected cost of exercising it --
  ## run_existence_cache = FALSE at least skips computing the *other*,
  ## unneeded half of that same bootstrap.
  diag <- pooling_diagnostics(data, draws, n_boot = 10, run_excess_variance = TRUE,
                               run_existence_cache = FALSE)
  expect_true(is.numeric(diag$excess_variance$p_value))
})

test_that("pooling_diagnostics skips congeniality gracefully when nothing is elicited", {
  data <- make_toy_data()
  specs <- list(x = c(7))  # fixed only, nothing elicited
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, n_boot = 10, run_existence_cache = FALSE)
  expect_null(diag$congeniality)
})

test_that("pooling_diagnostics(run_congeniality = FALSE) skips the congeniality test", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, n_boot = 10, run_congeniality = FALSE,
                               run_existence_cache = FALSE)
  expect_null(diag$congeniality)
})

test_that("pooling_diagnostics(run_existence_cache = FALSE, run_excess_variance = FALSE) skips the bootstrap entirely", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, run_existence_cache = FALSE, run_excess_variance = FALSE)
  expect_null(diag$existence_boot)
  expect_null(diag$excess_variance)
})

test_that("pooling_diagnostics still accepts explicit treat_var/outcome_var/cutpoint_specs", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, treat_var = "D", outcome_var = "Y",
                               cutpoint_specs = specs, n_boot = 10,
                               run_existence_cache = FALSE)
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
  diag <- pooling_diagnostics(data, draws, n_boot = 10, run_existence_cache = FALSE)
  expect_length(diag$flatness, 0)
})

test_that("pooling_diagnostics skips the retention diagnostic when asked to", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, n_boot = 10, run_retention = FALSE,
                               run_existence_cache = FALSE)
  expect_null(diag$retention)
})

test_that("pooling_diagnostics computes pooled itself when not supplied", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, n_boot = 10, run_existence_cache = FALSE)
  expect_equal(diag$pooled, pool_draws(draws))
})

## print.ecem_pooling_diagnostics's branching logic, tested against
## hand-built objects with known p-values rather than relying on the
## stochastic tests above -- deterministic and independent of n_boot.
## cg_p = NULL reproduces run_congeniality = FALSE (or nothing elicited);
## ev_p = NULL reproduces the (now off-by-default) excess-variance test
## not having been requested.
make_fake_diag <- function(flat_p, cg_p, cg_n = 100, cg_n_bins = 3, cg_position = "mid",
                            ev_p = NULL, ev_ratio = 1,
                            retention = list(gap = c(0.01, -0.02), mean = -0.005, sd = 0.02),
                            exact = FALSE, M = 10, K = NULL, alpha = 0.05) {
  congeniality <- if (is.null(cg_p)) {
    NULL
  } else {
    list(p_value = cg_p, model = NULL, n = cg_n, n_bins = cg_n_bins,
         elicited_vars = "age", fixed_cutpoints = list(age = c(29, 64)), position = cg_position)
  }
  excess_variance <- if (is.null(ev_p)) {
    NULL
  } else {
    list(p_value = ev_p, ratio = ev_ratio, B_null_dist = numeric(0))
  }
  structure(
    list(
      flatness        = list(age = list(p_value = flat_p, model = NULL, n = 100)),
      congeniality    = congeniality,
      excess_variance = excess_variance,
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
  d <- make_fake_diag(flat_p = 0.5, cg_p = 0.5)
  expect_output(print(d), "pool by Rubin's rules with no caveat")
})

test_that("print.ecem_pooling_diagnostics: congeniality rejects -> existence_test guidance", {
  d <- make_fake_diag(flat_p = 0.5, cg_p = 0.01)
  expect_output(print(d), "FSATT DRIFTS")
  expect_output(print(d), "existence_test\\(\\)")
  expect_output(print(d), "mechanism-locator")
})

test_that("print.ecem_pooling_diagnostics: flatness rejects, congeniality doesn't -> flag discrepancy", {
  d <- make_fake_diag(flat_p = 0.01, cg_p = 0.5)
  expect_output(print(d), "HETEROGENEOUS")
  expect_output(print(d), "flag the discrepancy")
  expect_output(print(d), "complete explanation")
})

test_that("print.ecem_pooling_diagnostics adapts guidance when retention wasn't computed", {
  d <- make_fake_diag(flat_p = 0.01, cg_p = 0.5, retention = NULL)
  expect_output(print(d), "Rerun with run_retention = TRUE")
})

test_that("print.ecem_pooling_diagnostics recommends trying another position when retention doesn't explain a discrepancy", {
  d <- make_fake_diag(flat_p = 0.01, cg_p = 0.5,
                       retention = list(gap = c(5, -5), mean = 4, sd = 6))
  expect_output(print(d), "may simply sit where drift")
  expect_output(print(d), "\"low\" or \"high\"")
})

test_that("print.ecem_pooling_diagnostics reports NA congeniality p-values as undecided, not a verdict", {
  d <- make_fake_diag(flat_p = 0.01, cg_p = NA_real_)
  expect_output(print(d), "p = NA")
  expect_output(print(d), "returned no verdict")
})

test_that("print.ecem_pooling_diagnostics notes when congeniality wasn't run at all", {
  d <- make_fake_diag(flat_p = 0.01, cg_p = NULL)
  expect_output(print(d), "not run")
  expect_output(print(d), "wasn't run -- rerun with run_congeniality = TRUE")
})

test_that("print.ecem_pooling_diagnostics reports the superseded excess-variance test only when present, and doesn't use it to decide anything", {
  d_without <- make_fake_diag(flat_p = 0.5, cg_p = 0.5, ev_p = NULL)
  expect_false(grepl("Excess-variance test", paste(capture.output(print(d_without)), collapse = "\n")))

  d_with <- make_fake_diag(flat_p = 0.5, cg_p = 0.5, ev_p = 0.01)
  out <- paste(capture.output(print(d_with)), collapse = "\n")
  expect_match(out, "superseded design")
  expect_match(out, "Not used below")
  ## Still pools cleanly per the (non-rejecting) congeniality test, even
  ## though the superseded excess-variance p-value would have rejected.
  expect_match(out, "pool by Rubin's rules with no caveat")
})

test_that("print.ecem_pooling_diagnostics returns x invisibly", {
  d <- make_fake_diag(flat_p = 0.5, cg_p = 0.5)
  expect_invisible(print(d))
})
