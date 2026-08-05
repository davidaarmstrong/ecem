#' Pre-matching flatness test on one elicited covariate (illustrative)
#'
#' A simple linear-model interaction F-test for whether the treatment
#' effect is flat over a covariate's elicited range, restricted to the
#' union of that range. **This is a placeholder** for the Crump et al.
#' (2008) / generic machine-learning heterogeneity tests the paper actually
#' recommends -- it is fine for exercising the rest of the pipeline, but
#' should be swapped out before using this function for real inference on
#' whether `tau(x)` is flat.
#'
#' @param data A data frame containing `treat_var`, `outcome_var`, and
#'   `xe_var`.
#' @param treat_var Character; name of the 0/1 treatment indicator column.
#' @param outcome_var Character; name of the outcome column.
#' @param xe_var Character; name of the elicited covariate to test.
#' @param xe_range Numeric length-2 vector, the union of the elicited range
#'   for `xe_var` (e.g. spanning all of that covariate's elicited
#'   cutpoint ranges).
#'
#' @return A list with elements `p_value` (from the interaction F-test),
#'   `model` (the fitted interaction model), and `n` (the restricted
#'   sample size used).
#'
#' @export
flatness_test_XE <- function(data, treat_var, outcome_var, xe_var, xe_range) {
  in_range <- data[[xe_var]] >= xe_range[1] & data[[xe_var]] <= xe_range[2]
  sub <- data[in_range, , drop = FALSE]

  f0 <- stats::as.formula(paste(outcome_var, "~", treat_var, "+", xe_var))
  f1 <- stats::as.formula(paste(outcome_var, "~", treat_var, "*", xe_var))
  fit0 <- stats::lm(f0, data = sub)
  fit1 <- stats::lm(f1, data = sub)
  cmp  <- stats::anova(fit0, fit1)

  list(p_value = cmp[["Pr(>F)"]][2], model = fit1, n = nrow(sub))
}

## Shared bootstrap engine for excess_variance_test() and existence_test().
##
## Both tests resample rows and rerun M matches per replicate, differing
## only in which null they substitute into the outcome before refitting:
## excess_variance_test() forces one shared tau_bar onto every draw, while
## existence_test() forces each draw's own tau_hat_m. elicit_and_match()
## never looks at the outcome column at all -- only fit_effect() does -- so
## the M matches formed within a replicate are identical either way, and
## only the (cheap) effect-fitting step needs rerunning for the second
## null. This computes both nulls' tau's off the *same* M matches per
## replicate, instead of matching 2 x n_boot x M times across two separate
## bootstraps.
##
## Set compute_shared/compute_own = FALSE (the default when tau_bar/
## tau_hat_m aren't supplied) to skip the half you don't need -- e.g.
## excess_variance_test() run standalone has no tau_hat_m to give an own
## null, and existence_test() run standalone has no tau_bar.
##
## estimator/covariates are passed straight through to fit_effect() and
## should match whatever `draws` was actually computed with (callers with a
## `draws` object recover these from its attributes; see run_M_draws()) --
## refitting bootstrap replicates with a different estimator than the one
## that produced the real draws would make B_null/tau_own answer the wrong
## question.
##
## Returns a list with `B_null` (`NULL` if !compute_shared) and `tau_own`
## (an n_boot x M matrix, `NULL` if !compute_own) -- deliberately left
## unsummarized, so any stat_fun can be applied to `tau_own` later without
## rerunning this bootstrap (see existence_test()'s `diagnostics` argument).
bootstrap_congeniality <- function(data, treat_var, outcome_var, cutpoint_specs, M,
                                    tau_bar = NULL, tau_hat_m = NULL, n_boot = 200,
                                    compute_shared = !is.null(tau_bar),
                                    compute_own = !is.null(tau_hat_m),
                                    estimator = "regression", covariates = NULL,
                                    progress = interactive()) {
  n <- nrow(data)
  D <- data[[treat_var]]
  Y <- data[[outcome_var]]

  Y_null_shared <- if (compute_shared) Y - tau_bar * D else NULL
  Y_null_own    <- if (compute_own) matrix(Y, nrow = n, ncol = M) - outer(D, tau_hat_m) else NULL

  pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0, max = n_boot, style = 3) else NULL
  if (!is.null(pb)) on.exit(close(pb), add = TRUE)

  B_null  <- if (compute_shared) numeric(n_boot) else NULL
  tau_own <- if (compute_own) matrix(NA_real_, nrow = n_boot, ncol = M) else NULL

  for (b in seq_len(n_boot)) {
    rows      <- sample.int(n, n, replace = TRUE)
    boot_base <- data[rows, , drop = FALSE]
    Yb_shared <- if (compute_shared) Y_null_shared[rows] else NULL
    Yb_own    <- if (compute_own) Y_null_own[rows, , drop = FALSE] else NULL

    tau_shared_b <- if (compute_shared) numeric(M) else NULL
    for (m in seq_len(M)) {
      matched_m <- elicit_and_match(boot_base, treat_var, cutpoint_specs)

      if (compute_shared) {
        boot_data <- boot_base
        boot_data[[outcome_var]] <- Yb_shared
        tau_shared_b[m] <- fit_effect(boot_data, treat_var, outcome_var, matched_m,
                                       estimator = estimator, covariates = covariates)$tau_hat
      }
      if (compute_own) {
        boot_data <- boot_base
        boot_data[[outcome_var]] <- Yb_own[, m]
        tau_own[b, m] <- fit_effect(boot_data, treat_var, outcome_var, matched_m,
                                     estimator = estimator, covariates = covariates)$tau_hat
      }
    }
    if (compute_shared) B_null[b] <- stats::var(tau_shared_b, na.rm = TRUE)

    if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
  }

  list(B_null = B_null, tau_own = tau_own)
}

## Excess-variance test's p_value/ratio from a B_null distribution --
## shared by excess_variance_test() (standalone) and pooling_diagnostics()
## (which computes B_null itself via the combined bootstrap above, so it
## can also cache tau_own for existence_test() to reuse).
excess_variance_from_B_null <- function(B_null, B_obs) {
  list(
    p_value     = mean(B_null >= B_obs, na.rm = TRUE),
    ratio       = B_obs / stats::median(B_null, na.rm = TRUE),
    B_null_dist = B_null
  )
}

#' Excess-variance bootstrap test on the realized draws
#'
#' Tests whether the observed between-draw variance `B` exceeds what
#' sampling noise alone predicts, computed directly on the realized draws
#' (as opposed to [flatness_test_XE()], which is a proxy computed on the
#' full retained sample before matching). Forces the null
#' \eqn{FSATT_1 = \dots = FSATT_M} onto the data by substituting the pooled
#' `tau_bar` for every unit's draw-specific contribution
#' (\eqn{Y^* = Y - \bar\tau D}), resamples rows with replacement, reruns
#' all `M` draws on that same resample, and records the resulting
#' between-draw variance `B*`. Repeating this traces out `B`'s null
#' distribution under exact congeniality -- correlation across draws is
#' preserved, since the resample is shared across all draws within a given
#' replicate -- without ever estimating or inverting a covariance matrix.
#'
#' @inheritParams run_draw
#' @param M Integer; number of draws to rerun within each bootstrap
#'   replicate (should match the `M` used to compute `B_obs`).
#' @param tau_bar The pooled point estimate from [pool_rubins_rules()] or
#'   [pool_draws()].
#' @param B_obs The observed between-draw variance from the same pooling
#'   call.
#' @param n_boot Integer; number of bootstrap replicates. Kept small in
#'   examples for speed; use at least a few hundred for actual inference.
#' @param estimator,covariates Passed straight through to [fit_effect()]
#'   for every refit inside the bootstrap. Should match whatever `draws`
#'   (or whatever produced `tau_bar`/`B_obs`) was actually computed with --
#'   refitting under a different estimator than the real draws used would
#'   make this test answer a different question than the one its `B_obs`
#'   was computed for. If you have a `draws` object, prefer calling this
#'   indirectly through [pooling_diagnostics()], which recovers these from
#'   `draws`'s attributes automatically rather than requiring you to repeat
#'   them here.
#' @param progress Logical; show a text progress bar over the `n_boot`
#'   bootstrap replicates (each of which reruns all `M` draws). Defaults to
#'   `interactive()`, so it stays quiet in scripts, `R CMD check`, and
#'   `testthat` runs. Written to `stderr()` via [utils::txtProgressBar()],
#'   so it never interferes with captured `stdout` output.
#'
#' @return A list with elements `p_value`, `ratio` (`B_obs` divided by the
#'   median of the null distribution), and `B_null_dist` (the full vector
#'   of bootstrapped `B*` values).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- simulate_population(N = 500)
#' specs <- list(age = list(c(25, 33), c(60, 68)), educ = c(12, 16))
#' draws <- run_M_draws(pop, "D", "Y", specs, M = 10)
#' pooled <- pool_rubins_rules(draws)
#' excess_variance_test(pop, "D", "Y", specs, M = 10,
#'                       tau_bar = pooled$tau_bar, B_obs = pooled$B,
#'                       n_boot = 20)
#' }
#'
#' @export
excess_variance_test <- function(data, treat_var, outcome_var, cutpoint_specs,
                                  M, tau_bar, B_obs, n_boot = 200,
                                  estimator = c("regression", "mean_diff"), covariates = NULL,
                                  progress = interactive()) {
  estimator <- match.arg(estimator)
  boot <- bootstrap_congeniality(
    data, treat_var, outcome_var, cutpoint_specs, M,
    tau_bar = tau_bar, n_boot = n_boot,
    estimator = estimator, covariates = covariates,
    progress = progress
  )
  excess_variance_from_B_null(boot$B_null, B_obs)
}

#' Simonsohn-style existence test
#'
#' A fallback for when [excess_variance_test()] rejects and pooling via
#' Rubin's rules is no longer trusted. For each draw `m` separately, forces
#' that draw's *own* null (zero effect: \eqn{Y^*_m = Y - \hat\tau_m D}) --
#' deliberately different from [excess_variance_test()], which forces one
#' shared null across all draws. Resamples rows once per replicate (the
#' same resampled rows used for every draw within that replicate) and
#' reruns the full pipeline for each draw on its own null-consistent,
#' resampled outcome.
#'
#' `treat_var`, `outcome_var`, and `cutpoint_specs` are recovered from
#' `draws`'s attributes by default (see [run_M_draws()]), and each draw's
#' `tau_hat` is read directly from `draws` -- so the common case is
#' `existence_test(data, draws)`. Pass any of the three explicitly to
#' override, or if `draws` doesn't carry them.
#'
#' If you already ran [pooling_diagnostics()] on this `draws` (e.g. because
#' its excess-variance test rejected and this is the recommended fallback),
#' pass that result as `diagnostics` to reuse its bootstrap instead of
#' rerunning one from scratch -- [pooling_diagnostics()] and
#' [excess_variance_test()] resample rows and rerun `M` matches per
#' replicate too, and matching doesn't depend on the outcome at all, so the
#' matches it already did can be refit under this test's own-draw null for
#' free instead of matching a second time. This turns `existence_test()`
#' from another `n_boot`-replicate bootstrap into an essentially instant
#' lookup.
#'
#' @param data A data frame, the same one used to produce `draws`. Needed
#'   again here because the test reruns the full pipeline on resampled
#'   data redrawn from `cutpoint_specs` -- unless `diagnostics` is
#'   supplied, in which case its cached bootstrap is used instead and
#'   `data` is not touched.
#' @param draws An object of class `"ecem_draws"`, as returned by
#'   [run_M_draws()].
#' @param treat_var,outcome_var,cutpoint_specs `NULL` (the default) to
#'   recover these from `draws`'s attributes; supply them explicitly to
#'   override, or if `draws` doesn't carry them. Ignored if `diagnostics`
#'   is supplied.
#' @param estimator,covariates `NULL` (the default) to recover these from
#'   `draws`'s attributes too (see [run_M_draws()]), so the bootstrap
#'   refits draws the same way `draws` was actually computed. Ignored if
#'   `diagnostics` is supplied (the cached bootstrap was already computed
#'   consistently with `draws` when [pooling_diagnostics()] made it).
#' @param diagnostics `NULL` (the default) to bootstrap fresh, or the
#'   result of [pooling_diagnostics(data, draws)][pooling_diagnostics()] to
#'   reuse its cached bootstrap (see Details). Must have been run on this
#'   same `draws` -- checked by comparing cached and current `tau_hat`
#'   values.
#' @param n_boot Integer; number of bootstrap replicates. Kept small in
#'   examples for speed; use at least a few hundred for actual inference.
#'   Ignored if `diagnostics` is supplied (its `n_boot` is used instead).
#' @param stat_fun Summary statistic applied to the draws' `tau_hat`
#'   values, both observed and under each bootstrap replicate's null.
#'   Defaults to [stats::median()]. Can be freely changed even when reusing
#'   a cached bootstrap via `diagnostics`, since the cache stores each
#'   replicate's raw per-draw `tau_hat`s, unsummarized.
#' @param progress Logical; show a text progress bar over the `n_boot`
#'   bootstrap replicates (each of which reruns all `M` draws, one per
#'   draw's own null). Defaults to `interactive()`, so it stays quiet in
#'   scripts, `R CMD check`, and `testthat` runs. Written to `stderr()` via
#'   [utils::txtProgressBar()], so it never interferes with captured
#'   `stdout` output. Ignored if `diagnostics` is supplied (nothing to
#'   show progress on).
#'
#' @return An object of class `"ecem_existence_test"` (see
#'   [print.ecem_existence_test()]): a list with elements `observed_stat`,
#'   `p_value`, `null_stats` (the full vector of the statistic under the
#'   null), `M` (number of draws), `n_boot`, `stat_label` (a short
#'   description of `stat_fun`, for printing), and `reused_bootstrap`
#'   (logical, whether `diagnostics`'s cached bootstrap was used).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- simulate_population(N = 500)
#' specs <- list(age = list(c(25, 33), c(60, 68)), educ = c(12, 16))
#' draws <- run_M_draws(pop, "D", "Y", specs, M = 10)
#' existence_test(pop, draws, n_boot = 20)
#'
#' ## Reusing an already-run pooling_diagnostics() bootstrap instead:
#' diag <- pooling_diagnostics(pop, draws, n_boot = 20)
#' existence_test(pop, draws, diagnostics = diag)
#' }
#'
#' @export
existence_test <- function(data, draws, treat_var = NULL, outcome_var = NULL,
                            cutpoint_specs = NULL, estimator = NULL, covariates = NULL,
                            diagnostics = NULL,
                            n_boot = 200, stat_fun = stats::median,
                            progress = interactive()) {
  stat_label <- deparse(substitute(stat_fun))
  if (nchar(stat_label) > 30 || grepl("^function", stat_label)) {
    stat_label <- "custom"
  }

  tau_hat_m <- vapply(draws, function(d) d$tau_hat, numeric(1))
  M <- length(tau_hat_m)
  observed_stat <- stat_fun(tau_hat_m, na.rm = TRUE)

  if (!is.null(diagnostics)) {
    if (!inherits(diagnostics, "ecem_pooling_diagnostics") || is.null(diagnostics$existence_boot)) {
      stop(
        "`diagnostics` must be the result of pooling_diagnostics() run on ",
        "this same `draws` (so it has a cached $existence_boot to reuse)."
      )
    }
    cached <- diagnostics$existence_boot
    if (!identical(cached$tau_hat_m, tau_hat_m)) {
      stop(
        "`diagnostics`'s cached bootstrap was computed from a different set ",
        "of draws (its tau_hat values don't match `draws`'s) -- rerun ",
        "pooling_diagnostics() on this `draws`, or omit `diagnostics` to let ",
        "existence_test() bootstrap fresh."
      )
    }

    null_stats <- apply(cached$tau_own, 1, stat_fun, na.rm = TRUE)
    n_boot <- cached$n_boot
    reused_bootstrap <- TRUE

  } else {
    if (is.null(treat_var))      treat_var      <- attr(draws, "treat_var")
    if (is.null(outcome_var))    outcome_var    <- attr(draws, "outcome_var")
    if (is.null(cutpoint_specs)) cutpoint_specs <- attr(draws, "cutpoint_specs")

    if (is.null(treat_var) || is.null(outcome_var) || is.null(cutpoint_specs)) {
      stop(
        "treat_var, outcome_var, and cutpoint_specs could not all be recovered ",
        "from `draws`. This happens if `draws` didn't come from run_M_draws() ",
        "-- pass whichever of treat_var/outcome_var/cutpoint_specs is missing ",
        "explicitly."
      )
    }

    ## estimator has no NULL-means-"missing" ambiguity to guard against the
    ## way treat_var/outcome_var/cutpoint_specs do -- draws predating this
    ## attribute (there shouldn't be any outside development) fall back to
    ## the current default rather than erroring.
    if (is.null(estimator))  estimator  <- attr(draws, "estimator")
    if (is.null(estimator))  estimator  <- "regression"
    if (is.null(covariates)) covariates <- attr(draws, "covariates")

    boot <- bootstrap_congeniality(
      data, treat_var, outcome_var, cutpoint_specs, M,
      tau_hat_m = tau_hat_m, n_boot = n_boot,
      estimator = estimator, covariates = covariates,
      progress = progress
    )
    null_stats <- apply(boot$tau_own, 1, stat_fun, na.rm = TRUE)
    reused_bootstrap <- FALSE
  }

  out <- list(
    observed_stat    = observed_stat,
    p_value          = mean(abs(null_stats) >= abs(observed_stat), na.rm = TRUE),
    null_stats       = null_stats,
    M                = M,
    n_boot           = n_boot,
    stat_label       = stat_label,
    reused_bootstrap = reused_bootstrap
  )
  class(out) <- "ecem_existence_test"
  out
}
