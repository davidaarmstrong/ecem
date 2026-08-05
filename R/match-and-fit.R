#' Resolve cutpoints, coarsen, match, and prune for one draw
#'
#' Resolves cutpoints for every covariate in `cutpoint_specs` (drawing a
#' value for elicited/regime covariates, passing fixed covariates through
#' unchanged), coarsens each one via [coarsen_one()], forms CEM strata from
#' their interaction, and prunes strata that lack common support (no
#' treated or no control units). Does not estimate any effect -- see
#' [fit_effect()] for that.
#'
#' @param data A data frame containing `treat_var` and every covariate
#'   named in `cutpoint_specs`.
#' @param treat_var Character; name of the 0/1 treatment indicator column
#'   in `data`.
#' @param cutpoint_specs A named list, one entry per covariate to include
#'   in the matching, in the format described in
#'   [draw_cutpoints_for_var()].
#'
#' @return A list with elements: `retained_idx` (row indices of `data`
#'   retained after pruning), `stratum` (a factor giving each retained
#'   unit's matched stratum), `resolved` (the per-covariate resolution from
#'   [draw_cutpoints_for_var()]), `kinds` (named character vector of each
#'   covariate's `kind`), and `regimes` (named character vector, `NA` for
#'   non-regime covariates and the drawn regime name for any covariate
#'   specified via [regime()] -- lets a run of draws be tabulated by which
#'   theory each one used).
#'
#' @export
elicit_and_match <- function(data, treat_var, cutpoint_specs) {
  var_names <- names(cutpoint_specs)
  stopifnot(all(var_names %in% names(data)), treat_var %in% names(data))

  resolved <- lapply(cutpoint_specs, draw_cutpoints_for_var)
  names(resolved) <- var_names

  coarsened <- vector("list", length(var_names))
  names(coarsened) <- var_names
  for (v in var_names) {
    coarsened[[v]] <- coarsen_one(data[[v]], resolved[[v]])
  }
  coarsened <- coarsened[!vapply(coarsened, is.null, logical(1))]

  n <- nrow(data)
  if (length(coarsened) == 0) {
    strata <- factor(rep("all", n))
  } else {
    strata <- do.call(interaction, c(coarsened, list(drop = TRUE)))
  }

  D <- data[[treat_var]]
  has_treat   <- tapply(D, strata, function(d) any(d == 1))
  has_control <- tapply(D, strata, function(d) any(d == 0))
  keep_levels <- levels(strata)[has_treat & has_control]

  retained <- which(strata %in% keep_levels)

  list(
    retained_idx = retained,
    stratum      = droplevels(strata[retained]),
    resolved     = resolved,
    kinds        = vapply(resolved, function(r) r$kind, character(1)),
    regimes      = vapply(resolved, function(r) {
      if (is.null(r$regime)) NA_character_ else r$regime
    }, character(1))
  )
}

## Iacus, King & Porro's CEM control-weighting scheme (2011/2012): treated
## units get weight 1; each control's weight is scaled so that, within its
## own stratum, the reweighted control mass matches the treated units'
## relative stratum sizes. A weighted difference of means using these
## weights recovers exactly the treated-count-weighted stratum
## mean-difference (see fit_effect()'s "mean_diff" branch), which is why
## these weights are the right ones to hand to a regression too -- doing
## so just lets the regression additionally adjust for covariates instead
## of relying on coarsening/stratification alone to do so.
##
## @param D 0/1 treatment vector, restricted to the retained/matched units.
## @param strat A factor giving each of those units' matched stratum.
## @return A numeric vector of weights, same length and order as `D`.
cem_weights <- function(D, strat) {
  T_s <- tapply(D, strat, sum)
  C_s <- tapply(1 - D, strat, sum)
  T_tot <- sum(D)
  C_tot <- sum(1 - D)

  w <- rep(1, length(D))
  is_control <- D == 0
  s_ctrl <- as.character(strat[is_control])
  w[is_control] <- (T_s[s_ctrl] / C_s[s_ctrl]) * (C_tot / T_tot)
  unname(w)
}

#' Estimate the treatment effect on a matched sample
#'
#' Two estimators are available, both operating on the matched object
#' returned by [elicit_and_match()]:
#'
#' - `"regression"` (the default): a weighted regression of `outcome_var`
#'   on `treat_var` plus `covariates`, using CEM's own control-weighting
#'   scheme ([cem_weights()]) on the retained sample -- the "run whatever
#'   model you would have run anyway, on the matched data" approach the
#'   CEM literature itself recommends (Iacus, King & Porro 2011/2012,
#'   building on Ho, Imai, King & Stuart's (2007) matching-as-
#'   nonparametric-preprocessing rationale). `tau_hat`/`var_hat` are `treat_var`'s
#'   coefficient and squared standard error. Unlike `"mean_diff"`, this can
#'   adjust for covariates' *uncoarsened* values, mopping up residual
#'   imbalance that coarsening-into-bins leaves behind -- at the cost of
#'   the usual regression caveats (functional-form assumptions on
#'   `covariates`, classical rather than robust standard errors here).
#' - `"mean_diff"`: a treated-count-weighted stratum mean-difference (no
#'   regression, no covariate adjustment beyond the matching strata
#'   themselves) -- simpler and always well-defined, but forgoes the
#'   efficiency and residual-imbalance adjustment a regression gives.
#'
#' Either way, the per-unit `tau(x)` surface used by the covariance
#' diagnostics ([cov_att_ate()], [cov_retention_per_draw()]) -- returned
#' here as `unit_tau_hat` -- is always the stratum mean-difference: those
#' diagnostics need a covariate-indexed CATE proxy, which a single
#' regression coefficient per draw can't provide, so the estimator choice
#' only changes what's reported as the headline `tau_hat`/`var_hat`, not
#' what feeds the diagnostics.
#'
#' @param data The same data frame passed to [elicit_and_match()].
#' @param treat_var Character; name of the 0/1 treatment indicator column.
#' @param outcome_var Character; name of the outcome column.
#' @param matched The list returned by [elicit_and_match()].
#' @param estimator `"regression"` (the default) or `"mean_diff"`; see
#'   Details.
#' @param covariates Character vector of covariate names for the
#'   `"regression"` estimator, using their raw (uncoarsened) values.
#'   `NULL` (the default) uses every non-excluded covariate in `matched`
#'   (i.e. `names(matched$kinds)` minus any `"excluded"` entries) -- the
#'   same covariates that entered the matching. Ignored for
#'   `"mean_diff"`.
#'
#' @return A list with elements: `tau_hat` (the within-draw effect
#'   estimate), `var_hat` (its estimated sampling variance), `n_used`
#'   (retained sample size), `strat_diff` (per-stratum mean differences),
#'   `unit_tau_hat` (each retained unit's own stratum diff -- the
#'   operationalized `tau(x)` surface used by the covariance diagnostics,
#'   computed the same way regardless of `estimator`), and `unit_idx` (the
#'   row indices these correspond to).
#'
#' @export
fit_effect <- function(data, treat_var, outcome_var, matched,
                        estimator = c("regression", "mean_diff"), covariates = NULL) {
  estimator <- match.arg(estimator)

  idx   <- matched$retained_idx
  D     <- data[[treat_var]][idx]
  Y     <- data[[outcome_var]][idx]
  strat <- matched$stratum

  if (length(idx) == 0) {
    return(list(tau_hat = NA_real_, var_hat = NA_real_, n_used = 0,
                strat_diff = numeric(0), unit_tau_hat = numeric(0),
                unit_idx = integer(0)))
  }

  strat_levels <- levels(strat)
  strat_stats <- lapply(strat_levels, function(lv) {
    rows <- which(strat == lv)
    d <- D[rows]; y <- Y[rows]
    n1 <- sum(d == 1); n0 <- sum(d == 0)
    y1 <- mean(y[d == 1]); y0 <- mean(y[d == 0])
    v1 <- if (n1 > 1) stats::var(y[d == 1]) / n1 else 0
    v0 <- if (n0 > 1) stats::var(y[d == 0]) / n0 else 0
    c(n1 = n1, n0 = n0, diff = y1 - y0, var = v1 + v0)
  })
  names(strat_stats) <- strat_levels
  agg <- do.call(rbind, strat_stats)
  unit_tau_hat <- unname(agg[as.character(strat), "diff"])

  if (estimator == "mean_diff") {
    w <- agg[, "n1"] / sum(agg[, "n1"])
    tau_hat <- unname(sum(w * agg[, "diff"]))
    var_hat <- unname(sum(w^2 * agg[, "var"]))

  } else {
    if (is.null(covariates)) {
      covariates <- names(matched$kinds)[matched$kinds != "excluded"]
    }

    cw <- cem_weights(D, strat)
    reg_data <- data[idx, covariates, drop = FALSE]
    reg_data[[treat_var]]   <- D
    reg_data[[outcome_var]] <- Y
    form <- stats::as.formula(paste(
      outcome_var, "~", treat_var,
      if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else ""
    ))
    reg_fit <- stats::lm(form, data = reg_data, weights = cw)
    co <- summary(reg_fit)$coefficients

    if (treat_var %in% rownames(co) && is.finite(co[treat_var, "Estimate"])) {
      tau_hat <- unname(co[treat_var, "Estimate"])
      var_hat <- unname(co[treat_var, "Std. Error"])^2
    } else {
      ## Rank-deficient/aliased D coefficient (e.g. a covariate perfectly
      ## collinear with D on this particular retained sample) -- fail
      ## gracefully rather than silently reporting a bogus number.
      tau_hat <- NA_real_
      var_hat <- NA_real_
    }
  }

  list(
    tau_hat      = tau_hat,
    var_hat      = var_hat,
    n_used       = length(idx),
    strat_diff   = agg[, "diff"],
    unit_tau_hat = unit_tau_hat,
    unit_idx     = idx
  )
}

#' Run one draw: elicit, match, and estimate
#'
#' Thin composition of [elicit_and_match()] and [fit_effect()] for a single
#' draw. This is the atomic unit reused, unmodified, by [run_M_draws()] and
#' by the resampling-based diagnostics ([excess_variance_test()],
#' [existence_test()]).
#'
#' @inheritParams elicit_and_match
#' @param outcome_var Character; name of the outcome column in `data`.
#' @inheritParams fit_effect
#'
#' @return A list combining the elements of [elicit_and_match()]'s output
#'   (under `matched`) and [fit_effect()]'s output (`tau_hat`, `var_hat`,
#'   `n_used`, `strat_diff`, `unit_tau_hat`, `unit_idx`).
#'
#' @export
run_draw <- function(data, treat_var, outcome_var, cutpoint_specs,
                      estimator = c("regression", "mean_diff"), covariates = NULL) {
  estimator <- match.arg(estimator)
  matched <- elicit_and_match(data, treat_var, cutpoint_specs)
  fit <- fit_effect(data, treat_var, outcome_var, matched, estimator = estimator, covariates = covariates)
  c(list(matched = matched), fit)
}

#' Run M draws, or every achievable configuration if that's small enough
#'
#' Runs `M` Monte Carlo draws via repeated calls to [run_draw()], unless
#' the achievable-configuration ceiling `K` (see [count_achievable_configs()];
#' this already sums within a [regime()] and multiplies across covariates,
#' so it accounts for every regime) is small enough to enumerate exactly
#' (`K <= exact_if_K_leq`), in which case it enumerates all `K` of them
#' instead and ignores `M` entirely.
#'
#' Exactness is signaled via attributes on the returned list, not a
#' different return shape: either way this returns a plain list of
#' per-draw records with the usual `$tau_hat`/`$var_hat`/`$matched`/...
#' fields, so [retention_matrix()], [cov_att_ate()],
#' [cov_retention_per_draw()], [pooled_unit_tau()], etc. all keep working
#' unmodified on the result. **Use [pool_draws()], not [pool_rubins_rules()]
#' directly, to pool it** -- `pool_rubins_rules()` assumes equal-weight
#' Monte Carlo draws and will silently compute the wrong between-draw
#' variance if handed weighted exact draws.
#'
#' `exact_if_K_leq` defaults to `NULL`, i.e. always Monte Carlo -- the
#' original behavior -- so existing callers are unaffected unless they opt
#' in. [excess_variance_test()] and [existence_test()] call `run_M_draws()`
#' internally and summarize it with an unweighted variance; turning
#' exactness on for them would also require switching that summary to a
#' weighted variance, which this package does not currently do, so those
#' two functions do not set `exact_if_K_leq`.
#'
#' @inheritParams run_draw
#' @param M Integer; number of Monte Carlo draws to run (ignored if exact
#'   enumeration is used instead -- see `exact_if_K_leq`).
#' @param exact_if_K_leq `NULL` (the default) to always run `M` Monte Carlo
#'   draws, or an integer threshold: if the achievable-configuration count
#'   `K` is at most this value, every achievable configuration is
#'   enumerated exactly instead of drawing `M` times.
#'
#' @return An object of class `"ecem_draws"`: a list of `M` (or `K`)
#'   per-draw records, as returned by [run_draw()], with attributes
#'   `exact` (logical), `weights` (numeric, summing to 1), and, when
#'   exact, `regimes` and `K`. Also carries `treat_var`, `outcome_var`,
#'   `cutpoint_specs`, `estimator`, and `covariates` as attributes, so that
#'   [pooling_diagnostics()], [existence_test()] (and, in principle, your
#'   own code) can recover them from `draws` alone rather than having to
#'   pass them again -- which matters more than it might seem for
#'   `estimator`/`covariates`, since the bootstrap-based diagnostics need
#'   to refit draws the *same* way `draws` was originally computed to give
#'   a meaningful null distribution. Has [print.ecem_draws()],
#'   [summary.ecem_draws()], and [as.data.frame.ecem_draws()] methods;
#'   still usable as a plain list everywhere else in this package (e.g.
#'   [retention_matrix()], [pool_draws()]).
#'
#' @export
run_M_draws <- function(data, treat_var, outcome_var, cutpoint_specs, M,
                         exact_if_K_leq = NULL,
                         estimator = c("regression", "mean_diff"), covariates = NULL) {
  estimator <- match.arg(estimator)
  exact <- FALSE
  if (!is.null(exact_if_K_leq)) {
    K <- count_achievable_configs(data, cutpoint_specs)$K
    if (K <= exact_if_K_leq) {
      enum  <- enumerate_configs(data, cutpoint_specs)
      draws <- lapply(enum$specs, function(spec_i) {
        run_draw(data, treat_var, outcome_var, spec_i, estimator = estimator, covariates = covariates)
      })
      attr(draws, "weights") <- enum$weights
      attr(draws, "regimes") <- enum$regimes
      attr(draws, "K")       <- K
      exact <- TRUE
    }
  }
  if (!exact) {
    draws <- lapply(seq_len(M), function(m) {
      run_draw(data, treat_var, outcome_var, cutpoint_specs, estimator = estimator, covariates = covariates)
    })
    attr(draws, "weights") <- rep(1 / M, M)
  }

  attr(draws, "exact")          <- exact
  attr(draws, "treat_var")      <- treat_var
  attr(draws, "outcome_var")    <- outcome_var
  attr(draws, "cutpoint_specs") <- cutpoint_specs
  attr(draws, "estimator")      <- estimator
  attr(draws, "covariates")     <- covariates
  class(draws) <- c("ecem_draws", "list")
  draws
}
