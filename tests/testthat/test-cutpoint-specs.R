test_that("regime() requires at least two named regimes", {
  expect_error(regime(a = list(c(1, 2))), "at least two")
})

test_that("regime() requires every regime to be named", {
  expect_error(regime(a = list(c(1, 2)), list(c(3, 4))), "named")
})

test_that("regime() defaults to uniform weights that sum to 1", {
  r <- regime(a = list(c(1, 2)), b = list(c(3, 4)))
  expect_s3_class(r, "cem_regime_spec")
  expect_equal(unname(r$weights), c(0.5, 0.5))
  expect_equal(names(r$weights), c("a", "b"))
})

test_that("regime() renormalizes supplied weights", {
  r <- regime(a = list(c(1, 2)), b = list(c(3, 4)), weights = c(3, 1))
  expect_equal(unname(r$weights), c(0.75, 0.25))
})

test_that("regime() errors if weights length doesn't match the number of regimes", {
  expect_error(regime(a = list(c(1, 2)), b = list(c(3, 4)), weights = c(1, 2, 3)))
})

test_that("draw_cutpoints_for_var treats NULL as excluded", {
  r <- draw_cutpoints_for_var(NULL)
  expect_equal(r$kind, "excluded")
  expect_null(r$cutpoints)
})

test_that("draw_cutpoints_for_var treats \"exact\" as exact", {
  r <- draw_cutpoints_for_var("exact")
  expect_equal(r$kind, "exact")
  expect_null(r$cutpoints)
})

test_that("draw_cutpoints_for_var treats a numeric vector as fixed, and sorts it", {
  r <- draw_cutpoints_for_var(c(16, 12))
  expect_equal(r$kind, "fixed")
  expect_equal(r$cutpoints, c(12, 16))
  expect_null(r$regime)
})

test_that("draw_cutpoints_for_var draws within range for elicited specs, sorted", {
  set.seed(1)
  r <- draw_cutpoints_for_var(list(c(25, 33), c(60, 68)))
  expect_equal(r$kind, "elicited")
  expect_length(r$cutpoints, 2)
  expect_true(r$cutpoints[1] >= 25 && r$cutpoints[1] <= 33)
  expect_true(r$cutpoints[2] >= 60 && r$cutpoints[2] <= 68)
  expect_true(r$cutpoints[1] <= r$cutpoints[2])
})

test_that("draw_cutpoints_for_var errors on an unrecognized spec", {
  expect_error(draw_cutpoints_for_var(TRUE), "Unrecognized")
})

test_that("coarsen_one returns NULL for an excluded covariate", {
  expect_null(coarsen_one(c(1, 2, 3), list(kind = "excluded", cutpoints = NULL)))
})

test_that("coarsen_one returns a factor of raw values for an exact-match covariate", {
  x <- c("a", "b", "a")
  out <- coarsen_one(x, list(kind = "exact", cutpoints = NULL))
  expect_s3_class(out, "factor")
  expect_equal(as.character(out), x)
})

test_that("coarsen_one bins correctly for fixed/elicited/regime cutpoints alike", {
  x <- c(1, 5, 10, 15)
  expected <- c(1L, 1L, 2L, 2L)
  expect_equal(coarsen_one(x, list(kind = "fixed",    cutpoints = 7)), expected)
  expect_equal(coarsen_one(x, list(kind = "elicited", cutpoints = 7)), expected)
  expect_equal(coarsen_one(x, list(kind = "regime",   cutpoints = 7)), expected)
})
