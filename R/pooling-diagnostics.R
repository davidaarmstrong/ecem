#' Run the congeniality diagnostics on a set of draws, together
#'
#' Runs the pre-matching flatness test ([flatness_test_XE()], once per
#' elicited or [regime()] covariate, over that covariate's full union
#' range), the excess-variance test on the realized draws
#' ([excess_variance_test()]), and, by default, the retention-interaction
#' diagnostic ([cov_retention_per_draw()]) -- the three checks Section 6.3
#' of the paper uses to decide whether pooling via Rubin's rules is
#' trustworthy, and what to report if it is not. [print.ecem_pooling_diagnostics()]
#' reports the results together with that decision logic spelled out.
#'
#' `treat_var`, `outcome_var`, and `cutpoint_specs` are recovered from
#' `draws` itself by default -- [run_M_draws()] stores them as attributes
#' for exactly this purpose -- so the common case is just
#' `pooling_diagnostics(data, draws)`. Pass them explicitly only if `draws`
#' didn't come from [run_M_draws()] (e.g. it was assembled some other way)
#' or you deliberately want to diagnose against a different spec than the
#' one that produced `draws`.
#'
#' This function does not itself run [existence_test()], even when the
#' excess-variance test rejects and the printed guidance recommends it --
#' that decision is left as a deliberate next step rather than something to
#' run automatically inside a diagnostics call. It does, however, do the
#' expensive part of existence_test()'s work as a side effect: the same `M`
#' matches per bootstrap replicate used here for the excess-variance test
#' are also refit under each draw's own null and cached (as
#' `existence_boot`), since matching never depends on the outcome and so
#' those matches are reusable for free. Pass this function's result as
#' `existence_test()`'s `diagnostics` argument to use that cache instead of
#' bootstrapping a second time.
#'
#' @param data A data frame, the same one used to produce `draws`. Needed
#'   again here because both the flatness and excess-variance tests rerun
#'   parts of the pipeline from scratch -- the flatness test on the full,
#'   pre-matching sample, and the excess-variance test on resampled data
#'   redrawn from `cutpoint_specs`. (Unlike `treat_var`/`outcome_var`/
#'   `cutpoint_specs`, `data` is not stored on `draws` -- keeping a full
#'   copy of the data as an attribute would be wasteful and could go
#'   stale, so it's always passed fresh.)
#' @param draws An object of class `"ecem_draws"`, as returned by
#'   [run_M_draws()].
#' @param treat_var,outcome_var,cutpoint_specs `NULL` (the default) to
#'   recover these from `draws`'s attributes; supply them explicitly to
#'   override, or if `draws` doesn't carry them (see Details).
#' @param estimator,covariates `NULL` (the default) to recover these from
#'   `draws`'s attributes too (see [run_M_draws()]), so the excess-variance
#'   bootstrap refits draws the same way `draws` was actually computed.
#' @param pooled Optional; the result of [pool_draws(draws)][pool_draws()].
#'   Computed automatically if not supplied.
#' @param n_boot Integer; number of bootstrap replicates for the
#'   excess-variance test. Kept modest by default for interactive use; use
#'   at least a few hundred for actual inference.
#' @param run_retention Logical; whether to compute the retention-
#'   interaction (`FSATT_m - ATT`) diagnostic. Defaults to `TRUE`; it is
#'   cheap (no bootstrapping) and is exactly what the printed guidance
#'   points to when flatness and excess-variance disagree.
#' @param alpha Significance level used to label each test's verdict in
#'   the printed output (and to choose which guidance
#'   [print.ecem_pooling_diagnostics()] prints). Defaults to `0.05`.
#' @param progress Logical; show a text progress bar over the excess-
#'   variance test's `n_boot` bootstrap replicates. Passed straight through
#'   to [excess_variance_test()]; see its documentation. Defaults to
#'   `interactive()`.
#'
#' @return An object of class `"ecem_pooling_diagnostics"`, a list with
#'   elements `flatness` (named list of [flatness_test_XE()] results, one
#'   per elicited/regime covariate), `excess_variance` (same shape as
#'   [excess_variance_test()]'s result), `retention` (`NULL` if
#'   `run_retention = FALSE`, else a list with `gap`, `mean`, `sd`),
#'   `pooled`, `M`, `K` (from [count_achievable_configs()], used by
#'   [print.ecem_pooling_diagnostics()] to decide whether recommending
#'   exact enumeration is actually practical here), `exact`, `alpha`, and
#'   `existence_boot` (a cache of this call's bootstrap matches, refit
#'   under each draw's own null, for [existence_test()]'s `diagnostics`
#'   argument to reuse -- see Details).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- simulate_population(N = 1000, heterogeneous = TRUE)
#' specs <- list(age = list(c(25, 33), c(60, 68)), educ = c(12, 16))
#' draws <- run_M_draws(pop, "D", "Y", specs, M = 15)
#' diag <- pooling_diagnostics(pop, draws, n_boot = 30)
#' diag
#' }
#'
#' @export
pooling_diagnostics <- function(data, draws, treat_var = NULL, outcome_var = NULL,
                                 cutpoint_specs = NULL, estimator = NULL, covariates = NULL,
                                 pooled = NULL, n_boot = 200,
                                 run_retention = TRUE, alpha = 0.05, progress = interactive()) {
  if (is.null(treat_var))      treat_var      <- attr(draws, "treat_var")
  if (is.null(outcome_var))    outcome_var    <- attr(draws, "outcome_var")
  if (is.null(cutpoint_specs)) cutpoint_specs <- attr(draws, "cutpoint_specs")

  if (is.null(treat_var) || is.null(outcome_var) || is.null(cutpoint_specs)) {
    stop(
      "treat_var, outcome_var, and cutpoint_specs could not all be recovered ",
      "from `draws`. This happens if `draws` didn't come from run_M_draws() ",
      "(e.g. it was assembled by hand or predates this package's tracking of ",
      "these attributes) -- pass whichever of treat_var/outcome_var/",
      "cutpoint_specs is missing explicitly."
    )
  }

  ## estimator has no NULL-means-"missing" ambiguity to guard against the
  ## way treat_var/outcome_var/cutpoint_specs do -- draws predating this
  ## attribute (there shouldn't be any outside development) fall back to
  ## the current default rather than erroring.
  if (is.null(estimator))  estimator  <- attr(draws, "estimator")
  if (is.null(estimator))  estimator  <- "regression"
  if (is.null(covariates)) covariates <- attr(draws, "covariates")

  if (is.null(pooled)) {
    pooled <- pool_draws(draws)
  }

  xe_ranges <- Filter(Negate(is.null), lapply(cutpoint_specs, elicited_union_range))
  flatness <- lapply(names(xe_ranges), function(v) {
    flatness_test_XE(data, treat_var, outcome_var, v, xe_ranges[[v]])
  })
  names(flatness) <- names(xe_ranges)

  M <- length(draws)
  tau_hat_m <- vapply(draws, function(d) d$tau_hat, numeric(1))

  ## One combined bootstrap serves both the excess-variance test here and
  ## existence_test()'s own-null refit, since elicit_and_match() never
  ## looks at the outcome -- see bootstrap_congeniality()'s comment.
  boot <- bootstrap_congeniality(
    data, treat_var, outcome_var, cutpoint_specs, M,
    tau_bar = pooled$tau_bar, tau_hat_m = tau_hat_m, n_boot = n_boot,
    estimator = estimator, covariates = covariates,
    progress = progress
  )
  excess_variance <- excess_variance_from_B_null(boot$B_null, B_obs = pooled$B)
  existence_boot <- list(tau_own = boot$tau_own, tau_hat_m = tau_hat_m, M = M, n_boot = n_boot)

  retention <- NULL
  if (run_retention) {
    treat_idx <- which(data[[treat_var]] == 1)
    tau_i <- pooled_unit_tau(draws, nrow(data))
    gap <- cov_retention_per_draw(draws, tau_i, treat_idx)
    retention <- list(gap = gap, mean = mean(gap, na.rm = TRUE), sd = stats::sd(gap, na.rm = TRUE))
  }

  ## Cheap (no bootstrapping, no matching) and lets print.ecem_pooling_diagnostics()
  ## give concrete, K-aware advice -- try exact enumeration if K turns out to
  ## be within reach, increase M otherwise -- rather than mentioning both
  ## options without saying which one actually applies here.
  K <- count_achievable_configs(data, cutpoint_specs)$K

  out <- list(
    flatness        = flatness,
    excess_variance = excess_variance,
    retention       = retention,
    pooled          = pooled,
    M               = M,
    K               = K,
    exact           = isTRUE(attr(draws, "exact")),
    alpha           = alpha,
    existence_boot  = existence_boot
  )
  class(out) <- "ecem_pooling_diagnostics"
  out
}
