## ecem demo: basic-workflow
##
## Elicit, draw, match, estimate, pool, and diagnose -- the full workflow
## on simulated data with one elicited covariate (age), one fixed (educ),
## one exact-match (region), and two excluded (income, noise).

## One spec exercising four of the five cutpoint kinds (see ?regime for
## the fifth, competing regimes -- demo("regime-competition", "ecem")).
cutpoint_specs <- list(
  age    = list(c(25, 33), c(60, 68)),  # elicited (X_E): two internal cutpoints
  educ   = c(12, 16),                   # fixed (X_F)
  region = "exact",                     # exact match on raw category
  income = NULL,                        # excluded
  noise  = NULL                         # excluded
)

set.seed(20260803)

pop <- simulate_population(N = 500, tau0 = 3, heterogeneous = TRUE)
M <- 100
draws <- run_M_draws(pop, "D", "Y", cutpoint_specs, M)
pooled <- pool_draws(draws)

diag <- pooling_diagnostics(pop, draws, pooled = pooled)
diag

if (diag$excess_variance$p_value < diag$alpha) {
  existence_test(pop, draws, diagnostics = diag)
}

label_diagnostics(pop, draws)



cat("\n--- Rubin's rules pooling across", M, "draws ---\n")
cat(sprintf("tau_bar = %.3f   Wbar = %.3f   B = %.4f   T = %.3f   lambda_hat = %.3f   df = %.1f\n",
            pooled$tau_bar, pooled$Wbar, pooled$B, pooled$T, pooled$lambda_hat, pooled$df))

treat_idx <- which(pop$D == 1)
R <- retention_matrix(draws, treat_idx)
cat("\nRetention rate per draw (share of treated units retained):\n")
print(round(colMeans(R), 3))

## label_diagnostics() bundles the ATT-ATE (Cov(tau_hat(X), p_hat(X))) and
## FSATT-ATT (Gap_r, per draw) covariance diagnostics from Section 4 and
## suggests which of FSATT/ATT/ATE the pooled tau_bar should be reported
## as. Its default propensity model uses only the covariates entering the
## matching (age, educ) -- pass `covariates` explicitly to include income,
## which was excluded from matching but is still worth conditioning on
## here since it plausibly drives treatment assignment.
lab <- label_diagnostics(pop, draws, covariates = c("age", "income", "educ"))
print(lab)

flat <- flatness_test_XE(pop, "D", "Y", "age", c(20, 75))
cat("\nPre-matching flatness test on age [illustrative lm-interaction version]: p =",
    round(flat$p_value, 4), "\n")

## n_boot is kept small below purely so the demo runs quickly; use
## something like 500-1000 for actual inference. pooling_diagnostics()
## bundles the flatness/excess-variance/retention checks above into one
## call and, as a side effect, caches its bootstrap matches for
## existence_test() to reuse below -- see its `diagnostics` argument.
cat("\nRunning pooling_diagnostics() (n_boot = 50)...\n")
diag <- pooling_diagnostics(pop, draws, n_boot = 50)
print(diag)

## Reuses diag's cached bootstrap instead of rerunning n_boot x M matches
## from scratch -- essentially instant.
cat("\nRunning Simonsohn-style existence test, reusing diag's bootstrap...\n")
ex <- existence_test(pop, draws, diagnostics = diag)
print(ex)
