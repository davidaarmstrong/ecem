test_that("pool_rubins_rules matches a hand-computed example", {
  draws <- list(
    list(tau_hat = 1, var_hat = 0.1, n_used = 20),
    list(tau_hat = 2, var_hat = 0.2, n_used = 20),
    list(tau_hat = 3, var_hat = 0.3, n_used = 20)
  )
  out <- pool_rubins_rules(draws)
  expect_equal(out$tau_bar, 2)
  expect_equal(out$Wbar, 0.2)
  expect_equal(out$B, 1)              # var(c(1, 2, 3)) == 1
  expect_equal(out$T, 0.2 + (1 + 1 / 3) * 1)
  expect_equal(out$M, 3)
})

test_that("pool_rubins_rules drops draws with NA tau_hat or var_hat", {
  draws <- list(
    list(tau_hat = 1,          var_hat = 0.1, n_used = 20),
    list(tau_hat = NA_real_,   var_hat = 0.2, n_used = 20),
    list(tau_hat = 3,          var_hat = 0.3, n_used = 20)
  )
  out <- pool_rubins_rules(draws)
  expect_equal(out$M, 2)
  expect_equal(out$tau_bar, 2)
})

test_that("pool_draws dispatches to Monte Carlo pooling when not flagged exact", {
  draws <- list(
    list(tau_hat = 1, var_hat = 0.1, n_used = 20),
    list(tau_hat = 2, var_hat = 0.2, n_used = 20)
  )
  attr(draws, "exact") <- FALSE
  out <- pool_draws(draws)
  expect_false(out$exact)
  expect_equal(out$tau_bar, 1.5)
})

test_that("pool_draws dispatches to exact, weighted pooling when flagged", {
  draws <- list(
    list(tau_hat = 0,  var_hat = 0.1),
    list(tau_hat = 10, var_hat = 0.1)
  )
  attr(draws, "exact") <- TRUE
  attr(draws, "weights") <- c(0.9, 0.1)
  out <- pool_draws(draws)
  expect_true(out$exact)
  expect_equal(out$tau_bar, 1)   # 0.9 * 0 + 0.1 * 10
  expect_true(is.na(out$df))
  expect_equal(out$lambda_hat, out$B / out$T)
})

test_that("pool_rubins_rules_exact reproduces an unweighted mean under equal weights", {
  results <- list(
    list(tau_hat = 4, var_hat = 0.5),
    list(tau_hat = 6, var_hat = 0.5)
  )
  out <- pool_rubins_rules_exact(results, weights = c(1, 1))
  expect_equal(out$tau_bar, 5)
  expect_equal(out$Wbar, 0.5)
  expect_equal(out$B, 1)   # 0.5*(4-5)^2 + 0.5*(6-5)^2
  expect_equal(out$K, 2)
})
