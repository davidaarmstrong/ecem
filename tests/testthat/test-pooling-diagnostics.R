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
  diag <- pooling_diagnostics(data, draws)

  expect_s3_class(diag, "ecem_pooling_diagnostics")
  expect_named(diag$flatness, "x")
  ## congeniality_test() is run at all three positions by default now; on
  ## this tiny (N = 8) toy dataset the matched interaction model may or
  ## may not have enough residual df to identify at each one (see
  ## congeniality_test()'s degeneracy guard) -- the point of this test is
  ## that pooling_diagnostics() itself doesn't error, at any position.
  expect_false(is.null(diag$congeniality))
  expect_named(diag$congeniality, c("low", "mid", "high"))
  for (pos in c("low", "mid", "high")) {
    expect_true(is.numeric(diag$congeniality[[pos]]$p_value))
  }
  expect_false(is.null(diag$retention))
  expect_length(diag$retention$gap, 5)
  ## 8 distinct x values fall in [0, 20], neither endpoint coinciding with
  ## one of them, so K = 8 + 1 = 9 (see count_achievable_configs()).
  expect_equal(diag$K, 9L)
})

test_that("pooling_diagnostics tests congeniality at low/mid/high by default", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 8), c(9, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws)
  expect_named(diag$congeniality, c("low", "mid", "high"))
  for (pos in c("low", "mid", "high")) {
    expect_equal(diag$congeniality[[pos]]$position, pos)
  }
})

test_that("pooling_diagnostics's congeniality_position restricts to just the requested position(s)", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 8), c(9, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)

  diag_one <- pooling_diagnostics(data, draws, congeniality_position = "high")
  expect_named(diag_one$congeniality, "high")
  expect_equal(diag_one$congeniality$high$position, "high")

  diag_two <- pooling_diagnostics(data, draws, congeniality_position = c("low", "high"))
  expect_named(diag_two$congeniality, c("low", "high"))
})

test_that("pooling_diagnostics's congeniality_min_n_per_arm defaults to 1 and is independent of the bootstrap's own min_n_per_arm", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 8), c(9, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)

  diag <- pooling_diagnostics(data, draws, congeniality_position = "mid")
  expect_equal(diag$congeniality$mid$min_n_per_arm, 1)

  ## Raising congeniality_min_n_per_arm past what either bin in this toy
  ## dataset can support (each has only 2 treated/2 control -- see the
  ## elicit_and_match()/congeniality_test() min_n_per_arm tests) prunes the
  ## congeniality test's matched sample to nothing (NA), mirroring
  ## congeniality_test()'s own graceful behavior, without affecting
  ## anything else pooling_diagnostics() computed.
  diag_strict <- pooling_diagnostics(data, draws, congeniality_position = "mid",
                                      congeniality_min_n_per_arm = 3)
  expect_equal(diag_strict$congeniality$mid$min_n_per_arm, 3)
  expect_true(is.na(diag_strict$congeniality$mid$p_value))
  ## Nothing else should have changed.
  expect_equal(diag_strict$K, diag$K)
})

test_that("pooling_diagnostics defaults congeniality_correction to bonferroni and returns adjusted p-values/verdict", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 8), c(9, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws)

  expect_equal(diag$congeniality_correction, "bonferroni")
  expect_named(diag$congeniality_p_adjusted, c("low", "mid", "high"))

  raw_p  <- vapply(diag$congeniality, function(r) r$p_value, numeric(1))
  non_na <- !is.na(raw_p)
  expect_equal(unname(diag$congeniality_p_adjusted[non_na]),
               unname(stats::p.adjust(raw_p[non_na], method = "bonferroni")))
  expect_equal(diag$congeniality_reject_any,
               if (!any(non_na)) NA else isTRUE(any(diag$congeniality_p_adjusted[non_na] < diag$alpha)))
})

test_that("pooling_diagnostics's congeniality_correction = 'none' reproduces the raw, unadjusted p-values", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 8), c(9, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, congeniality_correction = "none")

  raw_p <- vapply(diag$congeniality, function(r) r$p_value, numeric(1))
  expect_equal(unname(diag$congeniality_p_adjusted), unname(raw_p))
})

test_that("pooling_diagnostics's congeniality_correction is a no-op when only one position is tested", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 8), c(9, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)

  diag_b <- pooling_diagnostics(data, draws, congeniality_position = "high", congeniality_correction = "bonferroni")
  diag_h <- pooling_diagnostics(data, draws, congeniality_position = "high", congeniality_correction = "holm")
  diag_n <- pooling_diagnostics(data, draws, congeniality_position = "high", congeniality_correction = "none")

  expect_equal(diag_b$congeniality_p_adjusted[["high"]], diag_b$congeniality$high$p_value)
  expect_equal(diag_h$congeniality_p_adjusted, diag_b$congeniality_p_adjusted)
  expect_equal(diag_n$congeniality_p_adjusted, diag_b$congeniality_p_adjusted)
})

test_that("pooling_diagnostics errors on an invalid congeniality_correction", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  expect_error(pooling_diagnostics(data, draws, congeniality_correction = "sidak"), "should be one of")
})

test_that("pooling_diagnostics skips congeniality gracefully when nothing is elicited", {
  data <- make_toy_data()
  specs <- list(x = c(7))  # fixed only, nothing elicited
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws)
  expect_null(diag$congeniality)
})

test_that("pooling_diagnostics(run_congeniality = FALSE) skips the congeniality test", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, run_congeniality = FALSE)
  expect_null(diag$congeniality)
})

test_that("pooling_diagnostics still accepts explicit treat_var/outcome_var/cutpoint_specs", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, treat_var = "D", outcome_var = "Y",
                               cutpoint_specs = specs)
  expect_named(diag$flatness, "x")
})

test_that("pooling_diagnostics errors informatively if draws has no recoverable attributes", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 3)
  attr(draws, "treat_var") <- NULL
  attr(draws, "outcome_var") <- NULL
  attr(draws, "cutpoint_specs") <- NULL
  expect_error(pooling_diagnostics(data, draws), "could not all be recovered")
})

test_that("pooling_diagnostics reports no elicited covariates when the spec is entirely fixed", {
  data <- make_toy_data()
  specs <- list(x = c(7))  # fixed only, nothing elicited
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws)
  expect_length(diag$flatness, 0)
})

test_that("pooling_diagnostics skips the retention diagnostic when asked to", {
  data <- make_toy_data()
  specs <- list(x = list(c(0, 20)))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws, run_retention = FALSE)
  expect_null(diag$retention)
})

test_that("pooling_diagnostics computes pooled itself when not supplied", {
  data <- make_toy_data()
  specs <- list(x = c(7))
  draws <- run_M_draws(data, "D", "Y", specs, M = 5)
  diag <- pooling_diagnostics(data, draws)
  expect_equal(diag$pooled, pool_draws(draws))
})

## print.ecem_pooling_diagnostics's branching logic, tested against
## hand-built objects with known p-values rather than relying on the
## stochastic tests above -- deterministic.
## cg_p = NULL reproduces run_congeniality = FALSE (or nothing elicited).
## cg_p is normally a single (unnamed) p-value, tested at "mid" only; pass
## a named vector/list keyed by position (e.g. c(low = 0.3, mid = 0.1,
## high = 0.005)) to build a multi-position fixture instead.
## cg_correction mirrors pooling_diagnostics()'s congeniality_correction
## ("bonferroni" by default, matching the package default) -- this fixture
## builder computes congeniality_p_adjusted/congeniality_reject_any the
## same way pooling_diagnostics() itself does (via stats::p.adjust()), so
## fixtures stay consistent with real output without duplicating that
## logic differently here.
make_fake_diag <- function(flat_p, cg_p, cg_n = 100, cg_n_bins = 3,
                            cg_vcov_type = NULL, cg_min_n_per_arm = NULL,
                            cg_correction = "bonferroni",
                            retention = list(gap = c(0.01, -0.02), mean = -0.005, sd = 0.02),
                            exact = FALSE, M = 10, K = NULL, alpha = 0.05) {
  congeniality <- NULL
  congeniality_p_adjusted <- NULL
  congeniality_reject_any <- NULL
  if (!is.null(cg_p)) {
    if (is.null(names(cg_p))) {
      cg_p <- stats::setNames(list(cg_p), "mid")
    }
    res <- lapply(names(cg_p), function(pos) {
      list(p_value = cg_p[[pos]], model = NULL, n = cg_n, n_bins = cg_n_bins,
           elicited_vars = "age", fixed_cutpoints = list(age = c(29, 64)), position = pos,
           vcov_type = cg_vcov_type, min_n_per_arm = cg_min_n_per_arm)
    })
    names(res) <- names(cg_p)
    congeniality <- res

    raw_p <- vapply(congeniality, function(r) r$p_value, numeric(1))
    congeniality_p_adjusted <- stats::setNames(rep(NA_real_, length(raw_p)), names(raw_p))
    non_na <- !is.na(raw_p)
    if (any(non_na)) {
      congeniality_p_adjusted[non_na] <- stats::p.adjust(raw_p[non_na], method = cg_correction)
    }
    congeniality_reject_any <- if (!any(non_na)) {
      NA
    } else {
      isTRUE(any(congeniality_p_adjusted[non_na] < alpha))
    }
  }
  structure(
    list(
      flatness                = list(age = list(p_value = flat_p, model = NULL, n = 100)),
      congeniality             = congeniality,
      congeniality_correction  = cg_correction,
      congeniality_p_adjusted  = congeniality_p_adjusted,
      congeniality_reject_any  = congeniality_reject_any,
      retention    = retention,
      pooled       = NULL,
      M            = M,
      K            = K,
      exact        = exact,
      alpha        = alpha
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

test_that("print.ecem_pooling_diagnostics recommends testing the untested positions when retention doesn't explain a discrepancy", {
  ## Only "mid" was tested in this fixture -- "low"/"high" remain untested,
  ## so the guidance should point at covering them rather than asserting
  ## the discrepancy is already fully explained.
  d <- make_fake_diag(flat_p = 0.01, cg_p = 0.5,
                       retention = list(gap = c(5, -5), mean = 4, sd = 6))
  expect_output(print(d), "untested position\\(s\\) may simply be where drift")
  expect_output(print(d), "congeniality_position =")
})

test_that("print.ecem_pooling_diagnostics doesn't suggest more positions once low/mid/high were all already tested", {
  ## All three positions tested, none reject -- the discrepancy isn't
  ## explained by only having checked a weak-drift position, so the
  ## guidance shouldn't suggest testing "more" positions that don't exist.
  d <- make_fake_diag(flat_p = 0.01, cg_p = c(low = 0.4, mid = 0.5, high = 0.3),
                       retention = list(gap = c(5, -5), mean = 4, sd = 6))
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(out, "already tested at \"low\", \"mid\", and\n\"high\"")
  expect_false(grepl("congeniality_position =", out))
})

test_that("print.ecem_pooling_diagnostics: congeniality rejects at only one of several tested positions", {
  ## "low" and "mid" pass, "high" rejects -- the overall verdict should
  ## still be REJECT (per Section 6.3's "any position rejects" rule), and
  ## the summary line should name which position(s) rejected.
  d <- make_fake_diag(flat_p = 0.5, cg_p = c(low = 0.4, mid = 0.3, high = 0.006))
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(out, "no drift detected")
  expect_match(out, "FSATT DRIFTS across the elicited range")
  expect_match(out, "Overall: REJECTS -- drift detected at position\\(s\\) \"high\"")
  expect_match(out, "Congeniality test rejects: congeniality has failed")
})

test_that("print.ecem_pooling_diagnostics shows the adjusted p-value alongside the raw one by default", {
  d <- make_fake_diag(flat_p = 0.5, cg_p = c(low = 0.4, mid = 0.3, high = 0.006))
  out <- paste(capture.output(print(d)), collapse = "\n")
  ## high: raw p = 0.006, bonferroni-adjusted (x3) = 0.018.
  expect_match(out, "p = 0.006.*bonferroni-adjusted: p = 0.018", perl = TRUE)
})

test_that("print.ecem_pooling_diagnostics with congeniality_correction = 'none' reproduces the old uncorrected behavior", {
  d <- make_fake_diag(flat_p = 0.5, cg_p = c(low = 0.4, mid = 0.3, high = 0.006), cg_correction = "none")
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_false(grepl("adjusted", out))
  expect_match(out, "Overall: REJECTS -- drift detected at position\\(s\\) \"high\"")
})

test_that("print.ecem_pooling_diagnostics: holm rejects more positions than bonferroni at the same p-values", {
  ## low = .001, mid = .02, high = .03 -- bonferroni (x3) only clears low
  ## (.003 < .05; mid -> .06, high -> .09, neither < .05), while holm's
  ## step-down rejects all three (adjusted: .003, .04, .04, all < .05).
  ## This is exactly the power difference documented in
  ## pooling_diagnostics()'s congeniality_correction argument.
  cg_p <- c(low = 0.001, mid = 0.02, high = 0.03)

  d_bonf <- make_fake_diag(flat_p = 0.5, cg_p = cg_p, cg_correction = "bonferroni")
  out_bonf <- paste(capture.output(print(d_bonf)), collapse = "\n")
  expect_match(out_bonf, "Overall: REJECTS -- drift detected at position\\(s\\) \"low\"")

  d_holm <- make_fake_diag(flat_p = 0.5, cg_p = cg_p, cg_correction = "holm")
  out_holm <- paste(capture.output(print(d_holm)), collapse = "\n")
  expect_match(out_holm, "\"low\", \"mid\", \"high\"")
})

test_that("print.ecem_pooling_diagnostics: all tested positions degenerate to NA -> no verdict, not either outcome", {
  d <- make_fake_diag(flat_p = 0.01, cg_p = c(low = NA_real_, mid = NA_real_, high = NA_real_))
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(out, "returned no verdict at any tested position")
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

test_that("print.ecem_pooling_diagnostics falls back to HC4 and hides min_n_per_arm when the fake diag predates both fields", {
  ## Mirrors an older cached diagnostics object / hand-built fixture from
  ## before vcov_type and min_n_per_arm were tracked on congeniality_test()'s
  ## result -- both should degrade gracefully rather than printing "NULL".
  d <- make_fake_diag(flat_p = 0.5, cg_p = 0.5, cg_vcov_type = NULL, cg_min_n_per_arm = NULL)
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(out, "HC4-robust")
  expect_false(grepl("min_n_per_arm", out))
})

test_that("print.ecem_pooling_diagnostics reports vcov_type and min_n_per_arm in the congeniality header when non-default", {
  d <- make_fake_diag(flat_p = 0.5, cg_p = 0.5, cg_vcov_type = "HC4m", cg_min_n_per_arm = 5)
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(out, "HC4m-robust")
  expect_match(out, "min_n_per_arm = 5")
})

test_that("print.ecem_pooling_diagnostics hides min_n_per_arm from the header when it's just the default of 1", {
  d <- make_fake_diag(flat_p = 0.5, cg_p = 0.5, cg_vcov_type = "HC4", cg_min_n_per_arm = 1)
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_false(grepl("min_n_per_arm", out))
})

test_that("print.ecem_pooling_diagnostics returns x invisibly", {
  d <- make_fake_diag(flat_p = 0.5, cg_p = 0.5)
  expect_invisible(print(d))
})
