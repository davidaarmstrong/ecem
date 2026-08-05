## These tests target the bug found and fixed during development: exact
## enumeration bakes each config's realized cutpoints into a plain numeric
## vector before handing it to run_draw() (see enumerate_var_configs()),
## which strips the cem_regime_spec class. Without the "regime" attribute
## fix in enumerate_var_configs() / draw_cutpoints_for_var(), an exactly-
## enumerated draw would silently report kind = "fixed" with no regime
## recorded, even though it really was one regime or the other -- and
## $matched$regimes would come back NA for every such draw. Monte Carlo
## draws never hit this, since they resolve the actual regime() object
## each time, which is why it wasn't caught until exact enumeration was
## exercised directly.

test_that("draw_cutpoints_for_var recovers a regime tag from a realized cutpoint vector", {
  cuts <- c(5.5)
  attr(cuts, "regime") <- "low"
  r <- draw_cutpoints_for_var(cuts)
  expect_equal(r$kind, "regime")
  expect_equal(r$regime, "low")
  expect_equal(as.numeric(r$cutpoints), 5.5)
})

test_that("draw_cutpoints_for_var treats an untagged numeric vector as plain fixed", {
  r <- draw_cutpoints_for_var(c(12, 16))
  expect_equal(r$kind, "fixed")
  expect_null(r$regime)
})

test_that("Monte Carlo regime draws never mix cutpoints from different regimes", {
  set.seed(42)
  data <- make_toy_data()
  specs <- list(
    x = regime(
      low  = list(c(0, 2)),
      high = list(c(12, 14))
    )
  )
  draws <- run_M_draws(data, "D", "Y", specs, M = 30)
  regimes <- vapply(draws, function(d) d$matched$regimes[["x"]], character(1))
  expect_false(any(is.na(regimes)))
  expect_true(all(regimes %in% c("low", "high")))
  ## sanity check that sampling actually varies -- astronomically unlikely
  ## to fail by chance with M = 30 and default (uniform) regime weights
  expect_equal(length(unique(regimes)), 2)
})

test_that("regime information survives exact enumeration", {
  data <- make_toy_data()
  ## No observed x values fall in (5, 6) or (8, 9), so each regime has
  ## exactly one achievable cutpoint, giving K = 1 + 1 = 2.
  specs <- list(
    x = regime(
      low  = list(c(5, 6)),
      high = list(c(8, 9))
    )
  )
  info <- count_achievable_configs(data, specs)
  expect_equal(info$K, 2L)

  draws <- run_M_draws(data, "D", "Y", specs, M = 1, exact_if_K_leq = 1000)
  expect_true(attr(draws, "exact"))
  expect_length(draws, 2)

  regimes <- vapply(draws, function(d) d$matched$regimes[["x"]], character(1))
  expect_false(any(is.na(regimes)))
  expect_setequal(regimes, c("low", "high"))

  kinds <- vapply(draws, function(d) d$matched$kinds[["x"]], character(1))
  expect_true(all(kinds == "regime"))
})

test_that("K sums within a regime rather than multiplying across regimes", {
  data <- make_toy_data()
  specs <- list(
    x = regime(
      low  = list(c(5, 6)),   # 1 achievable gap
      high = list(c(8, 9))    # 1 achievable gap
    )
  )
  info <- count_achievable_configs(data, specs)
  expect_equal(info$K, 2L)   # 1 + 1, not 1 * 1
})
