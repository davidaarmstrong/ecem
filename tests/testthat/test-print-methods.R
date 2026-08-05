test_that("as.data.frame.ecem_draws returns one row per draw with known values", {
  data <- make_toy_data()
  ## estimator = "mean_diff" so the hand-verified numbers in
  ## helper-data.R's make_toy_data() apply directly -- this test is about
  ## as.data.frame()'s structure, not which estimator produced tau_hat.
  draws <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 3, estimator = "mean_diff")
  tab <- as.data.frame(draws)
  expect_equal(nrow(tab), 3)
  expect_equal(tab$tau_hat, rep(5.5, 3))
  expect_equal(tab$var_hat, rep(0.625, 3))
  expect_equal(tab$n_used, rep(8, 3))
  expect_equal(tab$weight, rep(1 / 3, 3))
  expect_false(any(grepl("^regime_", names(tab))))
})

test_that("as.data.frame.ecem_draws adds a regime column only for regime covariates", {
  data <- make_toy_data()
  specs <- list(x = regime(low = list(c(5, 6)), high = list(c(8, 9))))
  draws <- run_M_draws(data, "D", "Y", specs, M = 1, exact_if_K_leq = 1000)
  tab <- as.data.frame(draws)
  expect_true("regime_x" %in% names(tab))
  expect_setequal(tab$regime_x, c("low", "high"))
  expect_equal(tab$weight, attr(draws, "weights"))
})

test_that("print.ecem_draws prints a header and table, and returns the table invisibly", {
  data <- make_toy_data()
  draws <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 2)
  expect_output(out <- print(draws), "ecem_draws")
  expect_output(print(draws), "Monte Carlo")
  expect_equal(out, as.data.frame(draws))
  expect_invisible(print(draws))
})

test_that("print.ecem_draws reports exact enumeration correctly", {
  data <- make_toy_data()
  draws <- run_M_draws(data, "D", "Y", list(x = list(c(5, 6))), M = 1,
                        exact_if_K_leq = 1000)
  expect_output(print(draws), "exactly-enumerated")
})

test_that("summary.ecem_draws matches pool_draws() directly", {
  data <- make_toy_data()
  draws <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 4)
  expect_equal(summary(draws), pool_draws(draws))
})

test_that("pool_rubins_rules, pool_rubins_rules_exact, and pool_draws all return class ecem_pooled", {
  draws <- list(list(tau_hat = 1, var_hat = 0.1, n_used = 20),
                list(tau_hat = 2, var_hat = 0.2, n_used = 20))
  expect_s3_class(pool_rubins_rules(draws), "ecem_pooled")

  results <- list(list(tau_hat = 1, var_hat = 0.1), list(tau_hat = 2, var_hat = 0.2))
  expect_s3_class(pool_rubins_rules_exact(results, c(0.5, 0.5)), "ecem_pooled")

  data <- make_toy_data()
  mc_draws <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 3)
  expect_s3_class(pool_draws(mc_draws), "ecem_pooled")
  exact_draws <- run_M_draws(data, "D", "Y", list(x = list(c(5, 6))), M = 1,
                              exact_if_K_leq = 1000)
  expect_s3_class(pool_draws(exact_draws), "ecem_pooled")
})

test_that("print.ecem_pooled prints Monte Carlo pooling with a t-based CI, and returns x invisibly", {
  data <- make_toy_data()
  draws <- run_M_draws(data, "D", "Y", list(x = c(7)), M = 4)
  pooled <- pool_draws(draws)
  expect_output(out <- print(pooled), "Monte Carlo")
  expect_output(print(pooled), "Barnard-Rubin")
  expect_equal(out, pooled)
  expect_invisible(print(pooled))
})

test_that("print.ecem_pooled prints exact pooling with a normal-based CI", {
  data <- make_toy_data()
  draws <- run_M_draws(data, "D", "Y", list(x = list(c(5, 6))), M = 1,
                        exact_if_K_leq = 1000)
  pooled <- pool_draws(draws)
  expect_output(print(pooled), "exact pooling")
  expect_output(print(pooled), "normal reference distribution")
  ## df = NA under exact pooling is by design, not a bug -- the printout
  ## should say so rather than leaving a bare "NA" that looks like one.
  expect_output(print(pooled), "df = NA \\(exact enumeration\\)")
})

test_that("print.ecem_pooled works directly on pool_rubins_rules_exact() output", {
  ## Regardless of entry point -- not just via pool_draws() -- print should
  ## correctly infer exact pooling from the presence of $K and fill in
  ## lambda_hat/df at print time.
  results <- list(list(tau_hat = 4, var_hat = 0.5), list(tau_hat = 6, var_hat = 0.5))
  out <- pool_rubins_rules_exact(results, weights = c(1, 1))
  expect_output(print(out), "exact pooling")
  expect_output(print(out), "normal reference distribution")
})
