## ecem demo: achievable-configuration-ceiling
##
## Continuous covariates admit essentially unlimited distinct coarsenings,
## so Monte Carlo draws are fine. Real survey age is reported in whole
## years; rounding it shows the finite ceiling K on achievable
## coarsenings, and what exact enumeration looks like once K is small
## enough to replace Monte Carlo draws outright.

cutpoint_specs <- list(
  age    = list(c(25, 33), c(60, 68)),
  educ   = c(12, 16),
  region = "exact",
  income = NULL,
  noise  = NULL
)

set.seed(20260803)
pop <- simulate_population(N = 3000, tau0 = 3, heterogeneous = TRUE)

K_continuous <- count_achievable_configs(pop, cutpoint_specs)
cat("\nAchievable configurations with continuous age: K =", K_continuous$K, "\n")

pop_int <- pop
pop_int$age <- round(pop_int$age)

K_int <- count_achievable_configs(pop_int, cutpoint_specs)
cat("Achievable configurations with age rounded to years: K =", K_int$K, "\n")
cat("Achievable positions per covariate:\n")
print(K_int$counts_by_variable)

cat("\nEnumerating all", K_int$K, "configurations exactly (age rounded)...\n")
enum_draws <- run_enumerated_draws(pop_int, "D", "Y", cutpoint_specs)
pooled_exact <- pool_rubins_rules_exact(enum_draws$results, enum_draws$weights)
cat(sprintf("Exact pooling: tau_bar = %.3f   Wbar = %.3f   B = %.4f   T = %.3f   (K = %d, no MC error)\n",
            pooled_exact$tau_bar, pooled_exact$Wbar, pooled_exact$B, pooled_exact$T, pooled_exact$K))

M_over <- 20 * K_int$K
cat(sprintf("\nFor comparison, M = %d Monte Carlo draws on the same (rounded) data\n", M_over))
draws_over <- run_M_draws(pop_int, "D", "Y", cutpoint_specs, M_over)
pooled_over <- pool_rubins_rules(draws_over)
cat(sprintf("would give: tau_bar = %.3f   B = %.4f -- should match the exact values above\n",
            pooled_over$tau_bar, pooled_over$B))
cat("up to Monte Carlo noise, at far higher computational cost for no extra power.\n")

## The easy way to get the same exact result as above: pass
## exact_if_K_leq to run_M_draws() and call pool_draws() instead of
## pool_rubins_rules(). M is irrelevant here since K = 81 <= 200, so it's
## enumerated exactly rather than Monte Carlo sampled -- no need to call
## count_achievable_configs()/run_enumerated_draws() by hand first.
cat("\nSame thing the easy way (exact_if_K_leq = 200):\n")
draws_easy  <- run_M_draws(pop_int, "D", "Y", cutpoint_specs, M = 10, exact_if_K_leq = 200)
pooled_easy <- pool_draws(draws_easy)
cat(sprintf("exact = %s   K = %d   tau_bar = %.3f   B = %.4f   T = %.3f\n",
            attr(draws_easy, "exact"), length(draws_easy),
            pooled_easy$tau_bar, pooled_easy$B, pooled_easy$T))
cat("The M argument above is ignored whenever it enumerates exactly -- pass\n")
cat("whatever M you'd use as a Monte Carlo fallback for larger K.\n")
