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

## Bootstrap engine behind existence_test(). Resamples rows with
## replacement, substitutes each draw's own null (tau_hat_m) into the
## outcome, and reruns all M draws on that resample, recording the
## resulting per-draw tau under the null. Repeating this traces out each
## draw's own null distribution of tau.
##
## Previously also drove excess_variance_test() (a second, shared-null
## bootstrap computing a between-draw null variance B_null) -- removed
## along with that test, which the paper's Appendix documents as invalid
## (CEM's discrete matching solution is destabilized by resampling with
## replacement; see congeniality_test()'s documentation for the citation).
## Also previously cached by pooling_diagnostics() (run_existence_cache)
## for existence_test() to reuse -- removed too, since with only one
## resampling test left there was nothing left to amortize; see
## pooling_diagnostics()'s documentation for why.
##
## Set compute_own = FALSE (the default when tau_hat_m isn't supplied) to
## skip the bootstrap entirely and just return an empty result.
##
## estimator/covariates/min_n_per_arm are passed straight through to
## elicit_and_match()/fit_effect() and should match whatever `draws` was
## actually computed with (callers with a `draws` object recover these
## from its attributes; see run_M_draws()) -- refitting bootstrap
## replicates with a different estimator, covariate set, or retention
## threshold than the one that produced the real draws would make tau_own
## answer the wrong question.
##
## Returns a list with `tau_own` (an n_boot x M matrix, `NULL` if
## !compute_own) -- deliberately left unsummarized, so any stat_fun can be
## applied to it later.
bootstrap_congeniality <- function(data, treat_var, outcome_var, cutpoint_specs, M,
                                    tau_hat_m = NULL, n_boot = 200,
                                    compute_own = !is.null(tau_hat_m),
                                    estimator = "regression", covariates = NULL,
                                    min_n_per_arm = 1,
                                    progress = interactive()) {
  n <- nrow(data)
  D <- data[[treat_var]]
  Y <- data[[outcome_var]]

  Y_null_own <- if (compute_own) matrix(Y, nrow = n, ncol = M) - outer(D, tau_hat_m) else NULL

  pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0, max = n_boot, style = 3) else NULL
  if (!is.null(pb)) on.exit(close(pb), add = TRUE)

  tau_own <- if (compute_own) matrix(NA_real_, nrow = n_boot, ncol = M) else NULL

  for (b in seq_len(n_boot)) {
    rows      <- sample.int(n, n, replace = TRUE)
    boot_base <- data[rows, , drop = FALSE]
    Yb_own    <- if (compute_own) Y_null_own[rows, , drop = FALSE] else NULL

    for (m in seq_len(M)) {
      matched_m <- elicit_and_match(boot_base, treat_var, cutpoint_specs, min_n_per_arm = min_n_per_arm)

      if (compute_own) {
        boot_data <- boot_base
        boot_data[[outcome_var]] <- Yb_own[, m]
        tau_own[b, m] <- fit_effect(boot_data, treat_var, outcome_var, matched_m,
                                     estimator = estimator, covariates = covariates)$tau_hat
      }
    }

    if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
  }

  list(tau_own = tau_own)
}

#' Congeniality test: does FSATT drift within the elicited range?
#'
#' Tests the condition Rubin's-rules pooling actually needs (the paper's
#' Section on congeniality): whether the FSATT varies across the elicited
#' range of the coarsened covariate(s), rather than
#' [flatness_test_XE()]'s weaker pre-matching proxy for the same question
#' (run before any coarsening is applied, on the full retained-eligible
#' sample).
#'
#' This supersedes an earlier, resampling-based excess-variance-bootstrap
#' design the paper's Appendix documents as invalid -- CEM's discrete
#' matching solution is destabilized by resampling with replacement
#' (Abadie & Imbens 2008's bootstrap-inconsistency result for matching
#' estimators), which produced near-zero power regardless of true
#' heterogeneity. This test is not run on the `M` realized coarsening
#' draws at all, and needs no resampling and no null-imposed outcome
#' construction: it fixes ONE representative specification within each
#' elicited covariate's range, matches at that specification with CEM's
#' own pruning rule, and tests a treatment-by-bin interaction on the
#' resulting matched, common-support-restricted, CEM-weighted sample -- a
#' design adapted from the binning estimator of Hainmueller, Mummolo & Xu
#' (2019) for diagnosing linearity assumptions in interaction models.
#' Because CEM's weights are design weights (built to recover the FSATT),
#' not inverse-variance weights, the interaction is tested by a Wald test
#' using a heteroskedasticity-consistent sandwich covariance matrix (the
#' leverage-adaptive HC4 correction of Cribari-Neto 2004) rather than the
#' classical weighted-least-squares F-test, which assumes the latter.
#'
#' Every **elicited** (list-of-ranges) entry in `cutpoint_specs` gets a
#' single fixed cutpoint per range, taken from `position` (or from
#' `fixed_cutpoints`, if supplied for that covariate) rather than drawn.
#' **Fixed**, **exact**, and excluded (`NULL`) entries pass through to
#' [elicit_and_match()] unchanged and enter the test as additive controls.
#' When more than one covariate is elicited, they are tested *jointly*:
#' bins from every elicited covariate are crossed into one combined
#' factor, and the interaction is tested against that joint factor,
#' matching the joint-flatness requirement in the paper (marginal
#' flatness in each elicited covariate is necessary but not sufficient
#' for congeniality when more than one is elicited) rather than testing
#' each covariate's flatness separately. [regime()] entries are not
#' currently supported here -- a "fixed representative specification" has
#' no single natural meaning across competing regimes -- supply an
#' explicit entry in `fixed_cutpoints` for any regime covariate you want
#' included.
#'
#' Requires the \pkg{sandwich} and \pkg{lmtest} packages (both Imports of
#' this package, so ordinarily already installed).
#'
#' @param data A data frame containing `treat_var`, `outcome_var`, and
#'   every covariate in `cutpoint_specs`.
#' @param treat_var Character; name of the 0/1 treatment indicator column.
#' @param outcome_var Character; name of the outcome column.
#' @param cutpoint_specs A named list, the same format [elicit_and_match()]
#'   takes (see [draw_cutpoints_for_var()] for the five kinds of entry).
#' @param position `"mid"` (default), `"low"`, or `"high"`: which point
#'   within each elicited range becomes that range's fixed cutpoint -- the
#'   range's mean, lower bound, or upper bound respectively. Ignored for
#'   any covariate with an explicit override in `fixed_cutpoints`. Running
#'   the test at more than one `position` and comparing is a reasonable
#'   sensitivity check; this function itself only ever tests one.
#' @param fixed_cutpoints Optional named list, one entry per elicited
#'   covariate, of explicit numeric cutpoint vectors to use instead of
#'   deriving them from `position`.
#' @param min_n_per_arm Integer, default `1`; passed straight through to
#'   [elicit_and_match()] -- see its documentation for the full rationale.
#'   Raising this is the more targeted fix for the same numerically
#'   unstable-leverage problem `vcov_type` addresses when a treatment-by-
#'   bin cell is a near-singleton: it removes the offending cell from the
#'   matched sample entirely (a common-support restriction) rather than
#'   asking a different sandwich correction to cope with it. If you also
#'   want the headline point estimate (from [run_M_draws()]) to describe
#'   the *same* sample this test does, pass the same `min_n_per_arm` to
#'   both rather than only here.
#' @param vcov_type `"HC4"` (default), `"HC4m"`, `"HC3"`, `"HC5"`,
#'   `"HC1"`, or `"HC0"` -- which heteroskedasticity-consistent sandwich
#'   correction [sandwich::vcovHC()] applies. HC4 is the validated default
#'   (see the paper's Appendix for the comparison against HC1/HC3/HC4m/HC5
#'   that led to it), clean across every simulated condition checked. On
#'   real, unrounded data, though, a treatment-by-bin cell can end up with
#'   only one or two units -- rare, but common enough in application-sized
#'   samples -- pushing that observation's leverage (its hat value) close
#'   to, or exactly at, 1. HC4's correction divides by `(1 - h_i)^delta_i`
#'   for each observation, which is numerically unstable (and undefined at
#'   `h_i = 1`) regardless of the `delta_i` leverage-adaptivity cap --
#'   `sandwich::vcovHC()` warns exactly when this happens
#'   ("HC4 covariances are numerically unstable for hat values close to
#'   1..."). HC4m (Cribari-Neto & Souza-Vasconcellos, 2007) was
#'   developed specifically to stay bounded in this regime and is the
#'   principled thing to switch to if you see that warning on your own
#'   data, rather than a workaround -- it's already part of the same
#'   HC1/HC3/HC4/HC4m/HC5 family comparison the paper's appendix reports.
#'   Worth also just looking at which bin those high-leverage observations
#'   fall into (e.g. `hatvalues(cg$model)`) -- a near-singleton
#'   treatment-by-bin cell is itself worth knowing about independent of
#'   which correction handles it numerically.
#'
#' @return A list with elements `p_value` (from the `vcov_type`-robust
#'   Wald F-test; `NA` if the matched sample is empty or too small to
#'   identify the interaction model -- the same graceful-failure
#'   convention [fit_effect()] uses for an aliased coefficient, rather
#'   than erroring), `model` (the fitted interaction model), `n` (matched
#'   sample size), `n_bins` (levels of the joint bin factor tested),
#'   `elicited_vars` (which covariates were tested), `fixed_cutpoints`
#'   (the cutpoints actually used, named by covariate), `position`,
#'   `vcov_type` (which correction was actually used), and
#'   `min_n_per_arm` (the retention threshold actually used).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- simulate_population(N = 1000, heterogeneous = TRUE)
#' specs <- list(age = list(c(25, 33), c(60, 68)), educ = c(12, 16))
#' congeniality_test(pop, "D", "Y", specs)
#' congeniality_test(pop, "D", "Y", specs, position = "low")
#' }
#'
#' @export
congeniality_test <- function(data, treat_var, outcome_var, cutpoint_specs,
                               position = c("mid", "low", "high"),
                               fixed_cutpoints = NULL,
                               min_n_per_arm = 1,
                               vcov_type = c("HC4", "HC4m", "HC3", "HC5", "HC1", "HC0")) {
  if (!requireNamespace("sandwich", quietly = TRUE) || !requireNamespace("lmtest", quietly = TRUE)) {
    stop("congeniality_test() requires the 'sandwich' and 'lmtest' packages ",
         "(both Imports of ecem; install via install.packages(c('sandwich','lmtest')) ",
         "if somehow missing).")
  }
  position  <- match.arg(position)
  vcov_type <- match.arg(vcov_type)

  var_names   <- names(cutpoint_specs)
  is_regime   <- vapply(cutpoint_specs, inherits, logical(1), what = "cem_regime_spec")
  is_elicited <- vapply(cutpoint_specs, function(spec) is.list(spec) && !inherits(spec, "cem_regime_spec"),
                         logical(1))

  missing_regime_fix <- var_names[is_regime & !(var_names %in% names(fixed_cutpoints))]
  if (length(missing_regime_fix) > 0) {
    stop(
      "congeniality_test() does not derive a fixed representative specification ",
      "for regime() covariates automatically (there's no single natural choice ",
      "across competing regimes). Supply an explicit entry in fixed_cutpoints for: ",
      paste(missing_regime_fix, collapse = ", ")
    )
  }

  elicited_vars <- var_names[is_elicited]
  if (length(elicited_vars) == 0) {
    stop("congeniality_test() needs at least one elicited (list-of-ranges) ",
         "covariate in cutpoint_specs -- nothing to test congeniality of otherwise.")
  }

  fixed_specs   <- cutpoint_specs
  resolved_cuts <- list()
  for (v in elicited_vars) {
    if (!is.null(fixed_cutpoints[[v]])) {
      cuts <- sort(fixed_cutpoints[[v]])
    } else {
      ranges <- cutpoint_specs[[v]]
      cuts <- sort(vapply(ranges, function(r) {
        switch(position, mid = mean(r), low = r[1], high = r[2])
      }, numeric(1)))
    }
    resolved_cuts[[v]] <- cuts
    fixed_specs[[v]]   <- cuts
  }

  matched <- elicit_and_match(data, treat_var, fixed_specs, min_n_per_arm = min_n_per_arm)

  if (length(matched$retained_idx) == 0) {
    return(list(p_value = NA_real_, model = NULL, n = 0, n_bins = NA_integer_,
                elicited_vars = elicited_vars, fixed_cutpoints = resolved_cuts, position = position,
                vcov_type = vcov_type, min_n_per_arm = min_n_per_arm))
  }

  sub <- data[matched$retained_idx, , drop = FALSE]

  bin_list <- lapply(elicited_vars, function(v) {
    cut(sub[[v]], breaks = c(-Inf, resolved_cuts[[v]], Inf), include.lowest = TRUE)
  })
  names(bin_list) <- elicited_vars
  sub$.bin <- if (length(bin_list) == 1) {
    bin_list[[1]]
  } else {
    do.call(interaction, c(bin_list, list(drop = TRUE)))
  }

  ## Drop any elicited bin (or, for >1 elicited covariate, joint bin) the
  ## matched-and-retained sample doesn't actually populate. CEM's own
  ## pruning (a bin with only one treatment arm has no counterfactual and
  ## is dropped entirely) or a position/fixed_cutpoints choice that puts a
  ## cutpoint outside the retained sample's range can both leave a
  ## declared factor level with zero rows behind it -- contributing an
  ## identically-zero column to the design matrix at best, and, if it
  ## collapses .bin to a single remaining level, making R's own contrasts
  ## machinery error outright (a factor needs at least 2 levels to get a
  ## contrast matrix) rather than fitting at all. Fail gracefully instead,
  ## the same convention used below for a rank-deficient/zero-residual-df
  ## fit, rather than letting that error propagate.
  sub$.bin <- droplevels(sub$.bin)
  if (nlevels(sub$.bin) < 2) {
    return(list(p_value = NA_real_, model = NULL, n = nrow(sub), n_bins = nlevels(sub$.bin),
                elicited_vars = elicited_vars, fixed_cutpoints = resolved_cuts, position = position,
                vcov_type = vcov_type, min_n_per_arm = min_n_per_arm))
  }

  w <- cem_weights(sub[[treat_var]], matched$stratum)

  ## Every matched covariate (elicited ones included, at their raw
  ## uncoarsened values) enters as an additive linear control, alongside
  ## the joint bin factor -- generalizes the "+ X_E + X_F" structure the
  ## design was validated with to an arbitrary number of elicited/fixed
  ## covariates.
  linear_controls <- names(matched$kinds)[matched$kinds != "excluded"]
  rhs_extra <- if (length(linear_controls) > 0) {
    paste("+", paste(linear_controls, collapse = " + "))
  } else {
    ""
  }
  f0 <- stats::as.formula(paste(outcome_var, "~", treat_var, "+ .bin", rhs_extra))
  f1 <- stats::as.formula(paste(outcome_var, "~", treat_var, "* .bin", rhs_extra))
  fit0 <- stats::lm(f0, data = sub, weights = w)
  fit1 <- stats::lm(f1, data = sub, weights = w)

  ## Small/degenerate matched samples (too few units per treatment x bin
  ## cell) can leave the interaction model rank-deficient or without
  ## residual degrees of freedom to test against -- fail gracefully with
  ## an NA p-value, the same convention fit_effect() uses for an aliased
  ## D coefficient, rather than letting lmtest::waldtest() error out.
  rank_ok     <- fit1$rank == length(stats::coef(fit1)) && !anyNA(stats::coef(fit1))
  df_resid_ok <- stats::df.residual(fit1) > 0
  if (!rank_ok || !df_resid_ok) {
    return(list(p_value = NA_real_, model = fit1, n = nrow(sub), n_bins = nlevels(sub$.bin),
                elicited_vars = elicited_vars, fixed_cutpoints = resolved_cuts, position = position,
                vcov_type = vcov_type, min_n_per_arm = min_n_per_arm))
  }

  wt <- tryCatch(
    lmtest::waldtest(fit0, fit1, vcov = function(x) sandwich::vcovHC(x, type = vcov_type), test = "F"),
    error = function(e) NULL
  )
  p_value <- if (is.null(wt)) NA_real_ else wt[["Pr(>F)"]][2]

  list(
    p_value         = p_value,
    model           = fit1,
    n               = nrow(sub),
    n_bins          = nlevels(sub$.bin),
    elicited_vars   = elicited_vars,
    fixed_cutpoints = resolved_cuts,
    position        = position,
    vcov_type       = vcov_type,
    min_n_per_arm   = min_n_per_arm
  )
}

#' Simonsohn-style existence test
#'
#' A fallback for when [congeniality_test()] rejects and pooling via
#' Rubin's rules is no longer trusted. For each draw `m` separately, forces
#' that draw's *own* null (zero effect: \eqn{Y^*_m = Y - \hat\tau_m D}).
#' Resamples rows once per replicate (the same resampled rows used for
#' every draw within that replicate) and reruns the full pipeline for each
#' draw on its own null-consistent, resampled outcome.
#'
#' `treat_var`, `outcome_var`, and `cutpoint_specs` are recovered from
#' `draws`'s attributes by default (see [run_M_draws()]), and each draw's
#' `tau_hat` is read directly from `draws` -- so the common case is
#' `existence_test(data, draws)`. Pass any of the three explicitly to
#' override, or if `draws` doesn't carry them.
#'
#' Each call bootstraps fresh (`n_boot` replicates, each rerunning all `M`
#' draws) -- there is nothing to precompute or reuse from
#' [pooling_diagnostics()], which does not itself resample (see its
#' documentation for why).
#'
#' @param data A data frame, the same one used to produce `draws`. Needed
#'   again here because the test reruns the full pipeline on resampled
#'   data redrawn from `cutpoint_specs`.
#' @param draws An object of class `"ecem_draws"`, as returned by
#'   [run_M_draws()].
#' @param treat_var,outcome_var,cutpoint_specs `NULL` (the default) to
#'   recover these from `draws`'s attributes; supply them explicitly to
#'   override, or if `draws` doesn't carry them.
#' @param estimator,covariates `NULL` (the default) to recover these from
#'   `draws`'s attributes too (see [run_M_draws()]), so the bootstrap
#'   refits draws the same way `draws` was actually computed.
#' @param min_n_per_arm `NULL` (the default) to recover this from `draws`'s
#'   `"min_n_per_arm"` attribute too (falling back to `1`, i.e. CEM's usual
#'   any-unit-per-arm retention rule, if `draws` doesn't carry one -- e.g.
#'   older cached `draws` objects predating this parameter), so the
#'   bootstrap prunes matched strata the same way `draws` was actually
#'   computed.
#' @param n_boot Integer; number of bootstrap replicates. Kept small in
#'   examples for speed; use at least a few hundred for actual inference.
#' @param stat_fun Summary statistic applied to the draws' `tau_hat`
#'   values, both observed and under each bootstrap replicate's null.
#'   Defaults to [stats::median()].
#' @param progress Logical; show a text progress bar over the `n_boot`
#'   bootstrap replicates (each of which reruns all `M` draws, one per
#'   draw's own null). Defaults to `interactive()`, so it stays quiet in
#'   scripts, `R CMD check`, and `testthat` runs. Written to `stderr()` via
#'   [utils::txtProgressBar()], so it never interferes with captured
#'   `stdout` output.
#'
#' @return An object of class `"ecem_existence_test"` (see
#'   [print.ecem_existence_test()]): a list with elements `observed_stat`,
#'   `p_value`, `null_stats` (the full vector of the statistic under the
#'   null), `M` (number of draws), `n_boot`, and `stat_label` (a short
#'   description of `stat_fun`, for printing).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- simulate_population(N = 500)
#' specs <- list(age = list(c(25, 33), c(60, 68)), educ = c(12, 16))
#' draws <- run_M_draws(pop, "D", "Y", specs, M = 10)
#' existence_test(pop, draws, n_boot = 20)
#' }
#'
#' @export
existence_test <- function(data, draws, treat_var = NULL, outcome_var = NULL,
                            cutpoint_specs = NULL, estimator = NULL, covariates = NULL,
                            min_n_per_arm = NULL,
                            n_boot = 200, stat_fun = stats::median,
                            progress = interactive()) {
  stat_label <- deparse(substitute(stat_fun))
  ## An unmodified default (stat_fun = stats::median) substitutes the
  ## literal, namespace-qualified expression from the function signature
  ## rather than the bare name -- strip any "pkg::" qualifier so the
  ## printed label reads "median" whether the caller relied on the
  ## default or wrote stats::median explicitly themselves.
  stat_label <- sub("^[A-Za-z0-9_.]+::", "", stat_label)
  if (nchar(stat_label) > 30 || grepl("^function", stat_label)) {
    stat_label <- "custom"
  }

  tau_hat_m <- vapply(draws, function(d) d$tau_hat, numeric(1))
  M <- length(tau_hat_m)
  observed_stat <- stat_fun(tau_hat_m, na.rm = TRUE)

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
  if (is.null(min_n_per_arm)) min_n_per_arm <- attr(draws, "min_n_per_arm")
  if (is.null(min_n_per_arm)) min_n_per_arm <- 1

  boot <- bootstrap_congeniality(
    data, treat_var, outcome_var, cutpoint_specs, M,
    tau_hat_m = tau_hat_m, n_boot = n_boot,
    estimator = estimator, covariates = covariates,
    min_n_per_arm = min_n_per_arm,
    progress = progress
  )
  null_stats <- apply(boot$tau_own, 1, stat_fun, na.rm = TRUE)

  out <- list(
    observed_stat = observed_stat,
    p_value       = mean(abs(null_stats) >= abs(observed_stat), na.rm = TRUE),
    null_stats    = null_stats,
    M             = M,
    n_boot        = n_boot,
    stat_label    = stat_label
  )
  class(out) <- "ecem_existence_test"
  out
}
