#' Run the congeniality diagnostics on a set of draws, together
#'
#' Runs the pre-matching flatness test ([flatness_test_XE()], once per
#' elicited or [regime()] covariate, over that covariate's full union
#' range), the congeniality test ([congeniality_test()], at a single fixed
#' representative specification -- see its documentation for why this
#' replaced the earlier resampling-based excess-variance design), and, by
#' default, the retention-interaction diagnostic ([cov_retention_per_draw()])
#' -- the three checks Section 6.3 of the paper uses to decide whether
#' pooling via Rubin's rules is trustworthy, and what to report if it is
#' not. [print.ecem_pooling_diagnostics()] reports the results together
#' with that decision logic spelled out.
#'
#' `treat_var`, `outcome_var`, and `cutpoint_specs` are recovered from
#' `draws` itself by default -- [run_M_draws()] stores them as attributes
#' for exactly this purpose -- so the common case is just
#' `pooling_diagnostics(data, draws)`. Pass them explicitly only if `draws`
#' didn't come from [run_M_draws()] (e.g. it was assembled some other way)
#' or you deliberately want to diagnose against a different spec than the
#' one that produced `draws`.
#'
#' Neither the flatness test nor the congeniality test needs any
#' covariate in `cutpoint_specs` to actually be elicited -- if none is
#' (an all-fixed/exact spec), both are skipped gracefully (`flatness`
#' comes back an empty list and `congeniality` comes back `NULL`) rather
#' than erroring, since there is nothing to test congeniality of.
#'
#' This function does not itself run [existence_test()], even when the
#' congeniality test rejects and the printed guidance recommends it --
#' that decision is left as a deliberate next step. It can, however,
#' precompute the expensive part of `existence_test()`'s own work as a
#' side effect (`run_existence_cache = TRUE`, the default): a bootstrap
#' that resamples the retained-and-treated units and refits each draw's
#' own null, cached as `existence_boot` for `existence_test()`'s
#' `diagnostics` argument to reuse instead of bootstrapping a second time.
#' Set `run_existence_cache = FALSE` if you only want the cheap (no
#' resampling) congeniality/retention diagnostics and don't plan to call
#' `existence_test()` afterward -- that skips the bootstrap entirely,
#' which otherwise runs unconditionally regardless of whether anything
#' downstream needs it.
#'
#' `run_excess_variance = FALSE` by default: [excess_variance_test()] is
#' the earlier, resampling-based congeniality design the paper's Appendix
#' documents as invalid (near-zero power regardless of true heterogeneity,
#' for reasons rooted in the bootstrap's inconsistency for discrete
#' matching estimators). It's kept computable here, off by default, only
#' for reproducing that earlier design or comparing against it directly --
#' it is no longer part of the recommended workflow, and
#' [print.ecem_pooling_diagnostics()] does not use it to make any
#' decision.
#'
#' @param data A data frame, the same one used to produce `draws`. Needed
#'   again here because the flatness and congeniality tests both rerun
#'   parts of the pipeline from scratch -- the flatness test on the full,
#'   pre-matching sample, and the congeniality test on a freshly matched
#'   fixed specification. (Unlike `treat_var`/`outcome_var`/
#'   `cutpoint_specs`, `data` is not stored on `draws` -- keeping a full
#'   copy of the data as an attribute would be wasteful and could go
#'   stale, so it's always passed fresh.)
#' @param draws An object of class `"ecem_draws"`, as returned by
#'   [run_M_draws()].
#' @param treat_var,outcome_var,cutpoint_specs `NULL` (the default) to
#'   recover these from `draws`'s attributes; supply them explicitly to
#'   override, or if `draws` doesn't carry them (see Details).
#' @param estimator,covariates `NULL` (the default) to recover these from
#'   `draws`'s attributes too (see [run_M_draws()]), so any bootstrap
#'   refits draws the same way `draws` was actually computed.
#' @param pooled Optional; the result of [pool_draws(draws)][pool_draws()].
#'   Computed automatically if not supplied.
#' @param run_congeniality Logical; whether to run [congeniality_test()].
#'   Defaults to `TRUE`.
#' @param congeniality_position `"mid"` (default), `"low"`, or `"high"` --
#'   passed straight through to [congeniality_test()]'s `position`.
#' @param run_retention Logical; whether to compute the retention-
#'   interaction (`FSATT_m - ATT`) diagnostic. Defaults to `TRUE`; it is
#'   cheap (no bootstrapping, no rematching) and is exactly what the
#'   printed guidance points to when flatness and the congeniality test
#'   disagree.
#' @param run_existence_cache Logical; whether to precompute and cache the
#'   bootstrap [existence_test()] would otherwise redo itself (see
#'   Details). Defaults to `TRUE`.
#' @param run_excess_variance Logical; whether to also compute the
#'   superseded resampling-based excess-variance test (see Details).
#'   Defaults to `FALSE`.
#' @param n_boot Integer; number of bootstrap replicates, used if either
#'   `run_existence_cache` or `run_excess_variance` is `TRUE` (ignored
#'   otherwise). Kept modest by default for interactive use; use at least
#'   a few hundred for actual inference.
#' @param alpha Significance level used to label each test's verdict in
#'   the printed output (and to choose which guidance
#'   [print.ecem_pooling_diagnostics()] prints). Defaults to `0.05`.
#' @param progress Logical; show a text progress bar over the bootstrap's
#'   `n_boot` replicates, if it runs at all (see `run_existence_cache`/
#'   `run_excess_variance`). Defaults to `interactive()`.
#'
#' @return An object of class `"ecem_pooling_diagnostics"`, a list with
#'   elements `flatness` (named list of [flatness_test_XE()] results, one
#'   per elicited/regime covariate), `congeniality` (`NULL` if
#'   `run_congeniality = FALSE` or nothing in `cutpoint_specs` is
#'   elicited, else [congeniality_test()]'s result), `excess_variance`
#'   (`NULL` unless `run_excess_variance = TRUE`), `retention` (`NULL` if
#'   `run_retention = FALSE`, else a list with `gap`, `mean`, `sd`),
#'   `pooled`, `M`, `K` (from [count_achievable_configs()], used by
#'   [print.ecem_pooling_diagnostics()] to decide whether recommending
#'   exact enumeration is actually practical for pooling precision),
#'   `exact`, `alpha`, and `existence_boot` (`NULL` unless
#'   `run_existence_cache = TRUE`; a cache of this call's bootstrap
#'   matches, refit under each draw's own null, for `existence_test()`'s
#'   `diagnostics` argument to reuse -- see Details).
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
                                 pooled = NULL,
                                 run_congeniality = TRUE, congeniality_position = "mid",
                                 run_retention = TRUE,
                                 run_existence_cache = TRUE, run_excess_variance = FALSE,
                                 n_boot = 200, alpha = 0.05, progress = interactive()) {
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

  congeniality <- NULL
  if (run_congeniality && length(xe_ranges) > 0) {
    congeniality <- congeniality_test(data, treat_var, outcome_var, cutpoint_specs,
                                       position = congeniality_position)
  }

  retention <- NULL
  if (run_retention) {
    treat_idx <- which(data[[treat_var]] == 1)
    tau_i <- pooled_unit_tau(draws, nrow(data))
    gap <- cov_retention_per_draw(draws, tau_i, treat_idx)
    retention <- list(gap = gap, mean = mean(gap, na.rm = TRUE), sd = stats::sd(gap, na.rm = TRUE))
  }

  ## The bootstrap below is the expensive part of this function (it reruns
  ## CEM matching n_boot x M times) -- only pay for it if something
  ## actually needs it. run_existence_cache is what most callers want
  ## (existence_test() reuse, for free, if the congeniality test rejects);
  ## run_excess_variance reproduces the superseded design (see Details)
  ## and is off by default. One combined bootstrap serves both when both
  ## are requested, since elicit_and_match() never looks at the outcome --
  ## see bootstrap_congeniality()'s comment.
  excess_variance <- NULL
  existence_boot  <- NULL
  if (run_existence_cache || run_excess_variance) {
    tau_hat_m <- vapply(draws, function(d) d$tau_hat, numeric(1))
    boot <- bootstrap_congeniality(
      data, treat_var, outcome_var, cutpoint_specs, M,
      tau_bar   = if (run_excess_variance)  pooled$tau_bar else NULL,
      tau_hat_m = if (run_existence_cache) tau_hat_m       else NULL,
      n_boot = n_boot, estimator = estimator, covariates = covariates,
      progress = progress
    )
    if (run_excess_variance) {
      excess_variance <- excess_variance_from_B_null(boot$B_null, B_obs = pooled$B)
    }
    if (run_existence_cache) {
      existence_boot <- list(tau_own = boot$tau_own, tau_hat_m = tau_hat_m, M = M, n_boot = n_boot)
    }
  }

  ## Cheap (no bootstrapping, no matching) and lets print.ecem_pooling_diagnostics()
  ## give concrete, K-aware advice about POOLING precision (exact
  ## enumeration vs. a larger M) -- the congeniality test itself no longer
  ## depends on K or M at all, so this is reported for that separate
  ## question, not as a lever on the congeniality verdict.
  K <- count_achievable_configs(data, cutpoint_specs)$K

  out <- list(
    flatness        = flatness,
    congeniality    = congeniality,
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
