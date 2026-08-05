## ecem demo: regime-competition
##
## age is elicited under two competing theories rather than one: a
## lifecycle account (three breakpoints around early adulthood, midlife,
## retirement) and a generational account (four breakpoints at cohort
## boundaries). A draw is always entirely one or the other -- never, e.g.,
## the first lifecycle breakpoint combined with the last generational one.

set.seed(20260803)
pop <- simulate_population(N = 3000, tau0 = 3, heterogeneous = TRUE)

regime_specs <- list(
  age    = regime(
    lifecycle    = list(c(22, 28), c(40, 48), c(62, 68)),
    generational = list(c(24, 28), c(42, 46), c(58, 62), c(76, 80))
  ),
  educ   = c(12, 16),
  region = "exact",
  income = NULL,
  noise  = NULL
)

## Age rounded to years (as in the ceiling demo) so K is small enough for
## exact_if_K_leq below to actually trigger exact enumeration rather than
## falling back to Monte Carlo draws.
pop_int <- pop
pop_int$age <- round(pop_int$age)

set.seed(20260804)
M_reg <- 20
draws_reg <- run_M_draws(pop_int, "D", "Y", regime_specs, M_reg, exact_if_K_leq = 1000)
cat("\nexact =", attr(draws_reg, "exact"), "   K/M =", length(draws_reg), "\n")

regime_hits <- vapply(draws_reg, function(d) d$matched$regimes[["age"]], character(1))
cat("Regime drawn, by frequency across", length(draws_reg), "draws:\n")
print(table(regime_hits))

## pool_draws(), not pool_rubins_rules() directly -- draws_reg may be the
## exact-enumeration result above, whose draws carry unequal weights that
## pool_rubins_rules() would silently ignore.
pooled_reg <- pool_draws(draws_reg)
cat(sprintf("\nPooled across both regimes: tau_bar = %.3f   B = %.4f   T = %.3f\n",
            pooled_reg$tau_bar, pooled_reg$B, pooled_reg$T))

## Does it matter which theory you use? Compare tau_hat within each regime.
tau_by_regime <- tapply(
  vapply(draws_reg, function(d) d$tau_hat, numeric(1)),
  regime_hits, mean, na.rm = TRUE
)
cat("Mean tau_hat within each regime (informal check, not a formal test):\n")
print(round(tau_by_regime, 3))

## K sums within a regime and multiplies across covariates: K = K_lifecycle
## + K_generational here, since age is the only varying covariate.
K_regime <- count_achievable_configs(pop_int, regime_specs)
cat("\nAchievable configurations (age rounded), by regime:\n")
per_regime_K <- vapply(
  split(K_regime$per_variable$age, vapply(K_regime$per_variable$age, `[[`, character(1), "regime")),
  length, integer(1)
)
print(per_regime_K)
cat("Total K =", K_regime$K, "(= sum of the above, not their product)\n")
