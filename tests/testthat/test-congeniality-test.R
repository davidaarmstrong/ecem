## Direct sanity tests for congeniality_test(), separate from the
## pooling_diagnostics()-level tests in test-pooling-diagnostics.R (which
## exercise it indirectly and hand-build fake results to test
## print.ecem_pooling_diagnostics()'s branching). These use
## simulate_population() rather than the tiny make_toy_data() toy dataset,
## since a real (co)variance structure and adequate N are needed for the
## interaction model to actually identify -- make_toy_data() is reserved
## for the degeneracy-path test below, where under-identification is the
## point.

test_that("congeniality_test runs on a single elicited covariate and returns a well-formed result", {
  set.seed(1)
  pop <- simulate_population(N = 1000, heterogeneous = FALSE)
  specs <- list(age = list(c(25, 33), c(60, 68)))
  cg <- congeniality_test(pop, "D", "Y", specs)

  expect_true(is.numeric(cg$p_value))
  expect_false(is.na(cg$p_value))
  expect_true(cg$p_value >= 0 && cg$p_value <= 1)
  expect_s3_class(cg$model, "lm")
  expect_equal(cg$elicited_vars, "age")
  expect_equal(cg$position, "mid")
  expect_true(cg$n > 0)
  expect_true(cg$n_bins >= 2)
  ## mid position: fixed cutpoint per range is that range's mean.
  expect_equal(cg$fixed_cutpoints$age, c(mean(c(25, 33)), mean(c(60, 68))))
})

test_that("congeniality_test has power against real heterogeneity and stays quiet without it", {
  ## Not a formal calibration/power check (that's what the simulations/
  ## interflex-hc4-check.R and simulations/sim-congeniality.R scripts are
  ## for, at N/reps large enough to say something about size and power
  ## properly) -- just confirms the wiring behaves in the expected
  ## direction on one dataset large enough that a real age-varying effect
  ## should be easy to detect.
  set.seed(2)
  specs <- list(age = list(c(25, 33), c(60, 68)))

  pop_flat <- simulate_population(N = 2000, tau0 = 3, heterogeneous = FALSE)
  cg_flat <- congeniality_test(pop_flat, "D", "Y", specs)

  pop_het <- simulate_population(N = 2000, tau0 = 3, heterogeneous = TRUE)
  cg_het <- congeniality_test(pop_het, "D", "Y", specs)

  expect_true(cg_flat$p_value > 0.05)
  expect_true(cg_het$p_value < 0.05)
})

test_that("congeniality_test tests multiple elicited covariates jointly via a crossed bin factor", {
  set.seed(3)
  pop <- simulate_population(N = 1500, heterogeneous = FALSE)
  specs <- list(age = list(c(25, 33), c(60, 68)), educ = list(c(10, 12), c(14, 16)))
  cg <- congeniality_test(pop, "D", "Y", specs)

  expect_equal(sort(cg$elicited_vars), c("age", "educ"))
  ## age alone has 3 bins (2 cutpoints), educ alone has 3 bins -- crossed
  ## (before dropping empty combinations) is at most 3 * 3 = 9.
  expect_true(cg$n_bins >= 3)
  expect_true(cg$n_bins <= 9)
  expect_true(is.numeric(cg$p_value))
})

test_that("congeniality_test's position argument moves the fixed cutpoint to the range's bound", {
  set.seed(4)
  pop <- simulate_population(N = 1000, heterogeneous = FALSE)
  specs <- list(age = list(c(25, 33), c(60, 68)))

  cg_low  <- congeniality_test(pop, "D", "Y", specs, position = "low")
  cg_high <- congeniality_test(pop, "D", "Y", specs, position = "high")

  expect_equal(cg_low$fixed_cutpoints$age,  c(25, 60))
  expect_equal(cg_high$fixed_cutpoints$age, c(33, 68))
})

test_that("congeniality_test's fixed_cutpoints override takes precedence over position", {
  set.seed(5)
  pop <- simulate_population(N = 1000, heterogeneous = FALSE)
  specs <- list(age = list(c(25, 33), c(60, 68)))

  cg <- congeniality_test(pop, "D", "Y", specs, position = "low",
                           fixed_cutpoints = list(age = c(29, 64)))
  expect_equal(cg$fixed_cutpoints$age, c(29, 64))
})

test_that("congeniality_test passes fixed/exact/excluded covariates through as additive controls", {
  set.seed(6)
  pop <- simulate_population(N = 1000, heterogeneous = FALSE)
  specs <- list(age = list(c(25, 33), c(60, 68)), educ = c(12, 16), region = "exact", noise = NULL)
  cg <- congeniality_test(pop, "D", "Y", specs)

  expect_equal(cg$elicited_vars, "age")
  expect_true(is.numeric(cg$p_value))
  ## educ and region entered the matched model as linear controls; noise
  ## was excluded entirely -- both should show up (or not) in the fitted
  ## interaction model's terms accordingly.
  term_labels <- attr(stats::terms(cg$model), "term.labels")
  expect_true(any(grepl("^educ$", term_labels)))
  expect_true(any(grepl("^region$", term_labels)))
  expect_false(any(grepl("noise", term_labels)))
})

test_that("congeniality_test errors on a regime() covariate without an explicit fixed_cutpoints override", {
  set.seed(7)
  pop <- simulate_population(N = 500, heterogeneous = FALSE)
  specs <- list(age = regime(low = list(c(25, 33)), high = list(c(60, 68))))
  expect_error(congeniality_test(pop, "D", "Y", specs), "regime\\(\\) covariates")
})

test_that("congeniality_test accepts a regime() covariate once fixed_cutpoints supplies its cutpoint", {
  set.seed(8)
  pop <- simulate_population(N = 500, heterogeneous = FALSE)
  specs <- list(age = regime(low = list(c(25, 33)), high = list(c(60, 68))))
  ## regime() covariates aren't elicited (list-of-ranges) in the sense
  ## congeniality_test() tests -- with only a regime() entry and no other
  ## elicited covariate, there is still nothing to test congeniality of,
  ## so this should hit the "needs at least one elicited covariate" error
  ## rather than silently testing nothing.
  expect_error(
    congeniality_test(pop, "D", "Y", specs, fixed_cutpoints = list(age = c(29))),
    "needs at least one elicited"
  )
})

test_that("congeniality_test errors when nothing in cutpoint_specs is elicited", {
  pop <- simulate_population(N = 200, heterogeneous = FALSE)
  specs <- list(age = c(45))
  expect_error(congeniality_test(pop, "D", "Y", specs), "needs at least one elicited")
})

test_that("congeniality_test identifies fine on a small but adequately-covered matched sample", {
  ## make_toy_data() is deliberately tiny (N = 8), but a single interior
  ## cutpoint per range keeps both bins populated with both treatment arms
  ## after matching and pruning -- the interaction model still identifies
  ## (this is a sanity check that a small, non-degenerate case is NOT
  ## flagged as degenerate; see the tests below for cases that actually
  ## are).
  data <- make_toy_data()
  specs <- list(x = list(c(0, 8), c(9, 20)))
  cg <- congeniality_test(data, "D", "Y", specs)

  expect_true(is.numeric(cg$p_value))
  expect_false(is.na(cg$p_value))
  expect_true(cg$p_value >= 0 && cg$p_value <= 1)
  expect_equal(cg$elicited_vars, "x")
})

test_that("congeniality_test fails gracefully (NA, not an error) when pruning collapses the matched sample to a single bin", {
  ## Each of x = 1, 2, 3 is its own singleton bin here, and each is a
  ## single treatment arm on its own (x = 1, 2 are both D = 1 with no
  ## control; x = 3 is D = 0 with no treated) -- elicit_and_match() prunes
  ## all three for lacking a counterfactual, leaving only the x = 4..14
  ## bin retained. That collapses the matched sample onto a single
  ## populated .bin level, which R's own contrasts machinery can't build a
  ## treatment-by-bin interaction from -- congeniality_test() should catch
  ## this and return NA rather than letting that error propagate.
  data <- make_toy_data()
  specs <- list(x = list(c(1, 2), c(2, 3), c(3, 4)))
  cg <- congeniality_test(data, "D", "Y", specs)

  expect_true(is.na(cg$p_value))
  expect_true(cg$n_bins < 2)
})

test_that("congeniality_test returns an empty-match result gracefully when nothing is retained", {
  ## A range with no data on one side (no controls above x = 4) leaves
  ## elicit_and_match() nothing to retain in that stratum; force complete
  ## non-retention by asking for a cutpoint entirely outside the treated
  ## group's support so every stratum lacks one arm.
  data <- make_toy_data_with_unmatched()
  specs <- list(x = list(c(24, 27)))
  cg <- congeniality_test(data, "D", "Y", specs)

  ## An empty match (n = 0), a single-populated-bin match (n_bins < 2), or
  ## a degenerate one (too little residual df) are all acceptable outcomes
  ## here -- the point of this test is that congeniality_test() doesn't
  ## error regardless of which of those three it hits.
  expect_true(is.na(cg$p_value))
})
