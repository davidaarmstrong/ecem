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

## diag$congeniality is a named list, one result per position tested
## (low/mid/high by default) -- congeniality is treated as having failed
## if *any* tested position rejects, after a multiple-testing correction
## across however many positions were tested (congeniality_correction =
## "bonferroni" by default; the raw, uncorrected rule runs at roughly
## twice its nominal size -- see pooling_diagnostics()'s Details).
## diag$congeniality_reject_any is that already-corrected verdict, so
## there's no need to recompute it from the individual p-values here.
if (isTRUE(diag$congeniality_reject_any)) {
  existence_test(pop, draws)
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

## pooling_diagnostics() bundles the flatness/congeniality/retention
## checks above into one call -- it never bootstraps, so it's cheap
## regardless of how it's called.
cat("\nRunning pooling_diagnostics()...\n")
diag <- pooling_diagnostics(pop, draws)
print(diag)

## existence_test() always bootstraps fresh. n_boot is kept small here
## purely so the demo runs quickly; use something like 500-1000 for actual
## inference.
cat("\nRunning Simonsohn-style existence test (n_boot = 50)...\n")
ex <- existence_test(pop, draws, n_boot = 50)
print(ex)
