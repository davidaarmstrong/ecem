#' Tidy an ecem_draws object into a data frame
#'
#' Extracts one row per draw from an [ecem_draws][run_M_draws] object: the
#' individual `tau_hat` and `var_hat` estimates that [pool_draws()] would
#' otherwise combine into a single pooled number. Useful when you want the
#' draw-level detail directly -- e.g. to plot the distribution of
#' `tau_hat` across draws, or to hand to your own pooling/diagnostic code.
#'
#' @param x An object of class `"ecem_draws"`, as returned by
#'   [run_M_draws()].
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return A data frame with one row per draw and columns `draw` (index),
#'   `tau_hat`, `var_hat`, `n_used`, `weight` (uniform `1/M` for Monte
#'   Carlo draws, generally unequal for exact enumeration -- see
#'   [run_M_draws()]'s `exact_if_K_leq`), and, for any covariate specified
#'   via [regime()], a `regime_<covariate>` column giving which regime
#'   that draw used.
#'
#' @examples
#' pop <- simulate_population(N = 500)
#' draws <- run_M_draws(pop, "D", "Y", list(age = c(45)), M = 5)
#' as.data.frame(draws)
#'
#' @export
as.data.frame.ecem_draws <- function(x, ...) {
  tau_hat <- vapply(x, function(d) d$tau_hat, numeric(1))
  var_hat <- vapply(x, function(d) d$var_hat, numeric(1))
  n_used  <- vapply(x, function(d) d$n_used,  numeric(1))

  weight <- attr(x, "weights")
  if (is.null(weight)) weight <- rep(1 / length(x), length(x))

  out <- data.frame(
    draw    = seq_along(x),
    tau_hat = tau_hat,
    var_hat = var_hat,
    n_used  = n_used,
    weight  = weight,
    row.names = NULL
  )

  ## One regime_<covariate> column for every covariate that used regime()
  ## -- everything else in $matched$regimes is NA for every draw and is
  ## left out rather than cluttering the table with all-NA columns.
  regime_mat <- do.call(rbind, lapply(x, function(d) d$matched$regimes))
  has_regime <- apply(regime_mat, 2, function(col) any(!is.na(col)))
  for (v in colnames(regime_mat)[has_regime]) {
    out[[paste0("regime_", v)]] <- unname(regime_mat[, v])
  }

  out
}

#' Print an ecem_draws object
#'
#' Prints a short header (number of draws, Monte Carlo vs. exact, and each
#' covariate's kind) followed by the per-draw table from
#' [as.data.frame.ecem_draws()]. Returns that table invisibly, so
#' `tab <- print(draws)` (or, more directly, `as.data.frame(draws)`)
#' captures it.
#'
#' There is deliberately no separate "summary" of these numbers beyond
#' [summary.ecem_draws()], which just pools them -- see [pool_draws()].
#' The individual draws are the thing worth seeing raw; anything you'd
#' want to summarize about them (a pooled point estimate and variance) is
#' exactly what pooling already computes.
#'
#' @param x An object of class `"ecem_draws"`, as returned by
#'   [run_M_draws()].
#' @param digits Integer; number of significant digits to print for
#'   `tau_hat`, `var_hat`, and `weight`.
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return The per-draw data frame from [as.data.frame.ecem_draws()],
#'   invisibly.
#'
#' @examples
#' pop <- simulate_population(N = 500)
#' draws <- run_M_draws(pop, "D", "Y", list(age = c(45)), M = 5)
#' draws
#'
#' @export
print.ecem_draws <- function(x, digits = 4, ...) {
  exact <- isTRUE(attr(x, "exact"))
  kinds <- x[[1]]$matched$kinds

  cat(sprintf("<ecem_draws> %d %s\n", length(x),
              if (exact) "exactly-enumerated configuration(s)" else "Monte Carlo draw(s)"))
  cat("covariates: ", paste(sprintf("%s [%s]", names(kinds), kinds), collapse = ", "), "\n\n", sep = "")

  tab <- as.data.frame(x)
  print(tab, row.names = FALSE, digits = digits)

  invisible(tab)
}

#' Summarize an ecem_draws object by pooling it
#'
#' A thin wrapper around [pool_draws()]: the individual draws in an
#' [ecem_draws][run_M_draws] object don't have a separate "summary" beyond
#' their pooled point estimate and variance, so `summary()` just returns
#' that. See [print.ecem_draws()] / [as.data.frame.ecem_draws()] for the
#' individual `tau_hat`/`var_hat` values this pools together.
#'
#' @param object An object of class `"ecem_draws"`, as returned by
#'   [run_M_draws()].
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return The list returned by [pool_draws()]: `tau_bar`, `Wbar`, `B`,
#'   `T`, `lambda_hat`, `df`, `exact`.
#'
#' @examples
#' pop <- simulate_population(N = 500)
#' draws <- run_M_draws(pop, "D", "Y", list(age = c(45)), M = 5)
#' summary(draws)
#'
#' @export
summary.ecem_draws <- function(object, ...) {
  pool_draws(object)
}

#' Print a pooled result
#'
#' Prints the pooled point estimate, its standard error, and a confidence
#' interval, followed by the components (`Wbar`, `B`, `T`, `lambda_hat`,
#' `df`) that produced them. Works on the output of [pool_draws()],
#' [pool_rubins_rules()], or [pool_rubins_rules_exact()] -- whichever
#' produced `x`, "exact" pooling is detected from the presence of `K`
#' (`pool_rubins_rules_exact()`'s draw count) rather than requiring the
#' `exact` field that only [pool_draws()] sets, and `lambda_hat`/`df` are
#' filled in at print time if the object being printed didn't already have
#' them.
#'
#' The confidence interval uses a `t` reference distribution with the
#' Barnard-Rubin `df` under Monte Carlo pooling, as classical Rubin's
#' rules would. Under exact pooling `df` is `NA` -- there is no finite-`M`
#' uncertainty left to adjust for, since every achievable configuration
#' was actually run -- so the interval instead uses a normal reference
#' distribution, reflecting only the ordinary sampling uncertainty in
#' `Wbar` and `B` (see the paper's appendix remark on the achievable-
#' configuration ceiling for why a confidence interval remains meaningful
#' here even though there is no draw-selection uncertainty left).
#'
#' @param x An object of class `"ecem_pooled"`.
#' @param level Confidence level for the interval. Defaults to `0.95`.
#' @param digits Integer; number of significant digits to print.
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' pop <- simulate_population(N = 500)
#' draws <- run_M_draws(pop, "D", "Y", list(age = c(45)), M = 10)
#' pool_draws(draws)
#'
#' @export
print.ecem_pooled <- function(x, level = 0.95, digits = 4, ...) {
  exact <- if (!is.null(x$exact)) isTRUE(x$exact) else !is.null(x$K)
  n     <- if (!is.null(x$K)) x$K else x$M
  n_lab <- if (exact) "exactly-enumerated configuration(s)" else "Monte Carlo draw(s)"

  lambda_hat <- if (!is.null(x$lambda_hat)) x$lambda_hat else x$B / x$T
  df         <- if (!is.null(x$df)) x$df else NA_real_

  se    <- sqrt(x$T)
  alpha <- 1 - level
  if (exact || is.na(df)) {
    crit    <- stats::qnorm(1 - alpha / 2)
    ci_note <- "normal reference distribution (exact enumeration: no draw-selection uncertainty left to adjust for)"
  } else {
    crit    <- stats::qt(1 - alpha / 2, df = df)
    ci_note <- sprintf("t distribution, df = %.1f (Barnard-Rubin)", df)
  }
  ci <- x$tau_bar + c(-1, 1) * crit * se

  cat(sprintf("<ecem_pooled> %s pooling across %s %s\n",
              if (exact) "exact" else "Monte Carlo", n, n_lab))
  cat(sprintf("  tau_bar = %s   SE = %s   %d%% CI = [%s, %s]\n",
              formatC(x$tau_bar, digits = digits, format = "g"),
              formatC(se,        digits = digits, format = "g"),
              round(level * 100),
              formatC(ci[1], digits = digits, format = "g"),
              formatC(ci[2], digits = digits, format = "g")))
  cat(sprintf("  (%s)\n\n", ci_note))
  cat(sprintf("  Wbar = %s   B = %s   T = %s\n",
              formatC(x$Wbar, digits = digits, format = "g"),
              formatC(x$B,    digits = digits, format = "g"),
              formatC(x$T,    digits = digits, format = "g")))
  ## NA df has two different causes worth distinguishing in the printout:
  ## under exact enumeration it's NA by design (no finite-M uncertainty
  ## left to adjust for -- not a bug), whereas an NA under Monte Carlo
  ## pooling would mean the Barnard-Rubin formula actually degenerated
  ## (e.g. M too small), which is worth flagging differently.
  df_label <- if (!is.na(df)) {
    formatC(df, digits = 1, format = "f")
  } else if (exact) {
    "NA (exact enumeration)"
  } else {
    "NA (undefined)"
  }
  cat(sprintf("  lambda_hat = %s (share of T from between-draw variation)   df = %s\n",
              formatC(lambda_hat, digits = digits, format = "g"), df_label))

  invisible(x)
}

## No longer called by print.ecem_pooling_diagnostics() (the congeniality
## test it used to accompany doesn't depend on M/K at all, so this M/K
## guidance is no longer causally tied to its verdict -- see the comment
## on K in pooling_diagnostics()). Kept here, unused, as a documented
## building block for the separate exact-enumeration-vs-Monte-Carlo
## pooling-precision question, which K/M still answer.
##
## Originally: when flatness rejects, excess-variance doesn't, and the
## retention diagnostic doesn't explain it away, the next step depended
## on whether exact enumeration was actually within reach:
##   - already exact (every achievable config was run): not a power
##     problem at all -- there's no simulation error left to add power
##     against, so point at a real scope difference instead.
##   - K within a generous multiple of the M already run: exact
##     enumeration costs about the same as (or less than) what was already
##     paid for and removes simulation error entirely, so recommend it by
##     name with the actual call to make.
##   - K far beyond that: exact enumeration isn't practical, so a larger M
##     is the more realistic (if still imperfect) fix.
## The "20x" multiplier is a judgment call, not a formal rule -- it mirrors
## the oversampling factor used in demo("achievable-configuration-ceiling",
## "ecem") to show Monte Carlo converging to the exact answer, i.e. a
## budget most users would already consider generous.
power_next_step <- function(x) {
  if (isTRUE(x$exact)) {
    paste0(
      "  This was already exact enumeration (every achievable configuration was\n",
      "  run), so it isn't an M/K power problem -- there's no simulation error left\n",
      "  to add power against. The discrepancy more likely reflects a real scope\n",
      "  difference: the flatness test uses the full retained sample, before\n",
      "  matching prunes it, while B is computed only on the matched, retained draws.\n"
    )
  } else if (is.null(x$K)) {
    paste0(
      "  Consider exact enumeration (see count_achievable_configs()) if the\n",
      "  achievable-configuration count K is small, or a larger M otherwise.\n"
    )
  } else if (x$K <= 20 * x$M) {
    sprintf(paste0(
      "  K = %d achievable configurations is within reach of the %d draws you\n",
      "  already ran -- rerun with exact_if_K_leq = %d (or higher) in run_M_draws()\n",
      "  for an exact answer with no simulation error, rather than a larger M.\n"
    ), x$K, x$M, x$K)
  } else {
    sprintf(paste0(
      "  K = %d achievable configurations is far more than the %d draws you ran,\n",
      "  so full enumeration isn't practical here. A larger M is the more realistic\n",
      "  fix, though power is still capped by how close M gets to K.\n"
    ), x$K, x$M)
  }
}

#' Print a pooling-diagnostics result
#'
#' Reports the flatness and congeniality diagnostics from
#' [pooling_diagnostics()] together (and the retention-interaction and
#' excess-variance diagnostics, if computed), and spells out which of the
#' paper's Section 6.3 workflow outcomes they imply: pool cleanly, flag a
#' discrepancy (and check whether the retention diagnostic explains it --
#' and, if not, whether testing at a different position in the elicited
#' range clarifies things), or treat congeniality as having failed and
#' fall back on [existence_test()].
#'
#' @param x An object of class `"ecem_pooling_diagnostics"`.
#' @param digits Integer; number of significant digits to print.
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return `x`, invisibly.
#'
#' @export
print.ecem_pooling_diagnostics <- function(x, digits = 4, ...) {
  g <- function(v) formatC(v, digits = digits, format = "g")

  cat("<ecem_pooling_diagnostics>\n\n")

  flatness_reject <- FALSE
  if (length(x$flatness) == 0) {
    cat("Flatness pre-check: no elicited or regime covariates in this spec -- skipped.\n\n")
  } else {
    cat("Pre-matching flatness test (per elicited covariate, full retained sample):\n")
    flat_p <- vapply(x$flatness, function(f) f$p_value, numeric(1))
    for (v in names(x$flatness)) {
      verdict <- if (flat_p[[v]] < x$alpha) "HETEROGENEOUS" else "homogeneous"
      cat(sprintf("  %-12s p = %s   [%s]\n", v, g(flat_p[[v]]), verdict))
    }
    cat("\n")
    flatness_reject <- any(flat_p < x$alpha, na.rm = TRUE)
  }

  cg <- x$congeniality
  cg_reject <- FALSE
  if (is.null(cg)) {
    cat("Congeniality test: not run (run_congeniality = FALSE, or nothing in this\n")
    cat("spec is elicited -- there's nothing to test congeniality of).\n\n")
  } else {
    cg_reject <- isTRUE(cg$p_value < x$alpha)
    cat(sprintf("Congeniality test (fixed specification, position = \"%s\", matched sample):\n",
                cg$position))
    if (is.na(cg$p_value)) {
      cat(sprintf("  p = NA   n = %d   [matched sample too small/degenerate to test -- see congeniality_test()]\n",
                  cg$n))
    } else {
      cat(sprintf("  p = %s   n = %d   n_bins = %d   [%s]\n",
                  g(cg$p_value), cg$n, cg$n_bins,
                  if (cg_reject) "FSATT DRIFTS across the elicited range" else "no drift detected"))
    }
    cat("\n")
  }

  if (!is.null(x$excess_variance)) {
    ev <- x$excess_variance
    ev_reject <- ev$p_value < x$alpha
    cat(sprintf("Excess-variance test (superseded design, reproduced on request -- on the %d realized %s):\n",
                x$M, if (x$exact) "exactly-enumerated configuration(s)" else "Monte Carlo draw(s)"))
    cat(sprintf("  p = %s   B_obs / median(B*) = %s   [%s]\n",
                g(ev$p_value), g(ev$ratio),
                if (ev_reject) "excess variance detected" else "no excess variance detected"))
    cat("  (Not used below -- see the paper's Appendix for why this design is invalid,\n")
    cat("  and congeniality_test() for what replaced it.)\n\n")
  }

  if (!is.null(x$retention)) {
    cat("Retention-interaction diagnostic (Gap_r = FSATT_m - ATT, per draw):\n")
    cat(sprintf("  mean = %s   sd = %s\n", g(x$retention$mean), g(x$retention$sd)))
    cat("\n")
  }

  cat("--- What this implies (Section 6.3 workflow) ---\n")
  if (is.null(cg)) {
    cat("Congeniality test wasn't run -- rerun with run_congeniality = TRUE (and at\n")
    cat("least one elicited covariate) for a verdict; flatness alone is a screen, not\n")
    cat("a certificate.\n")
  } else if (is.na(cg$p_value)) {
    cat("Congeniality test returned no verdict (matched sample too small/degenerate at\n")
    cat("this specification) -- try a different position, or treat this as missing\n")
    cat("rather than as either outcome below.\n")
  } else if (cg_reject) {
    cat("Congeniality test rejects: congeniality has failed regardless of the flatness\n")
    cat("result. Report the draws descriptively and fall back on existence_test()\n")
    cat("rather than trusting the pooled T -- pass this object as diagnostics = to\n")
    cat("reuse its cached bootstrap instead of rerunning one from scratch (requires\n")
    cat("run_existence_cache = TRUE, the default).\n")
    if (!is.null(x$retention)) {
      cat("The retention diagnostic above is worth reading now as a mechanism-locator:\n")
      cat("it isolates how much of the drift routes through the retention channel.\n")
    } else {
      cat("Rerun with run_retention = TRUE to help locate the mechanism.\n")
    }
  } else if (!flatness_reject) {
    cat("Neither test rejects: pool by Rubin's rules with no caveat.\n")
  } else {
    cat("Flatness rejects but the congeniality test does not: pool by Rubin's rules,\n")
    cat("but flag the discrepancy.\n")
    if (!is.null(x$retention)) {
      cat("If the retention diagnostic above is small and stable across draws, that's a\n")
      cat("complete explanation for why the congeniality test stayed quiet despite real\n")
      cat("heterogeneity in tau(x).\n")
      cat("If it is not small, the fixed specification tested may simply sit where drift\n")
      cat("is weak -- rerun congeniality_test() (or pooling_diagnostics() with\n")
      cat("congeniality_position =) at \"low\" or \"high\" rather than \"mid\" before\n")
      cat("concluding the discrepancy is benign.\n")
    } else {
      cat("Rerun with run_retention = TRUE to check whether the retention channel\n")
      cat("explains the discrepancy, and consider testing at another position in the\n")
      cat("elicited range before concluding it's benign.\n")
    }
  }

  invisible(x)
}

#' Print an existence-test result
#'
#' Reports the observed statistic, its bootstrap p-value under each draw's
#' own null, and a verdict at the given `alpha`. See [existence_test()] for
#' when this test is the right fallback (excess-variance has rejected and
#' pooling via Rubin's rules is no longer trusted).
#'
#' Rejecting the null here is evidence that *some* nonzero effect exists
#' somewhere among the draws -- it is deliberately not evidence for any
#' single pooled point estimate, since congeniality has already failed and
#' the draws are no longer treated as interchangeable estimates of one
#' quantity. Report the individual draws (e.g. [as.data.frame.ecem_draws()])
#' descriptively rather than a pooled `tau_bar` from here on.
#'
#' @param x An object of class `"ecem_existence_test"`.
#' @param alpha Significance level for the printed verdict. Defaults to
#'   `0.05`.
#' @param digits Integer; number of significant digits to print.
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return `x`, invisibly.
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
print.ecem_existence_test <- function(x, alpha = 0.05, digits = 4, ...) {
  g <- function(v) formatC(v, digits = digits, format = "g")
  reject <- x$p_value < alpha

  cat(sprintf("<ecem_existence_test> Simonsohn-style existence test (%s, %d bootstrap replicate(s)%s)\n",
              x$stat_label, x$n_boot,
              if (isTRUE(x$reused_bootstrap)) ", reused from pooling_diagnostics()" else ""))
  cat(sprintf("  observed %s(tau_hat) = %s   p = %s   [%s]\n",
              x$stat_label, g(x$observed_stat), g(x$p_value),
              if (reject) "REJECT null of no effect" else "fail to reject"))
  cat(sprintf("  across %d draw(s)\n\n", x$M))

  if (reject) {
    cat("Rejecting is evidence that some nonzero effect exists among the draws --\n")
    cat("not validation of any single pooled point estimate. Congeniality has already\n")
    cat("failed here, so report the individual draws descriptively rather than a\n")
    cat("pooled tau_bar.\n")
  } else {
    cat("Fails to reject even the draw-specific nulls -- no evidence here that any\n")
    cat("draw's effect differs from zero.\n")
  }

  invisible(x)
}

## Ratio of a diagnostic term's magnitude to SE(tau_bar), and whether it
## counts as "appreciable" under the appreciable_frac_se heuristic. NA/0/
## non-finite SE (or an NA term) is treated conservatively as appreciable
## -- undecidable is not the same as negligible, and defaulting to
## negligible here would understate how far the estimate might be from
## the label being suggested.
label_term_status <- function(term, se, frac) {
  ratio <- if (is.na(term) || is.na(se) || !is.finite(se) || se <= 0) {
    NA_real_
  } else {
    abs(term) / se
  }
  appreciable <- if (is.na(ratio)) TRUE else ratio >= frac
  list(ratio = ratio, appreciable = appreciable)
}

#' Print a label-diagnostics result
#'
#' Reports the ATT-ATE term (`Cov(tau_hat(X), p_hat(X))`, and the
#' equivalent ATT - ATE gap) and the FSATT-ATT term
#' (\eqn{\overline{\mathrm{Gap}}_r}, mean and spread across draws) from
#' [label_diagnostics()], followed by a *suggested* FSATT/ATT/ATE label --
#' explicitly flagged here as a heuristic this package proposes (see
#' [label_diagnostics()]'s Details), not a threshold the paper itself
#' specifies; the paper deliberately leaves "appreciable" to the analyst.
#'
#' @param x An object of class `"ecem_label_diagnostics"`.
#' @param digits Integer; number of significant digits to print.
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- simulate_population(N = 1000, heterogeneous = TRUE)
#' specs <- list(age = list(c(25, 33), c(60, 68)), educ = c(12, 16))
#' draws <- run_M_draws(pop, "D", "Y", specs, M = 15)
#' label_diagnostics(pop, draws)
#' }
#'
#' @export
print.ecem_label_diagnostics <- function(x, digits = 4, ...) {
  g <- function(v) formatC(v, digits = digits, format = "g")
  fmt_ratio <- function(ratio) if (is.na(ratio)) "undefined" else g(ratio)

  cat("<ecem_label_diagnostics>\n\n")

  att_ate <- label_term_status(x$att_ate_gap, x$se, x$appreciable_frac_se)
  cat("ATT-ATE term (Proposition 1):\n")
  cat(sprintf("  Cov(tau_hat(X), p_hat(X)) = %s\n", g(x$cov_tau_p)))
  cat(sprintf("  implied ATT - ATE = %s   |gap| / SE(tau_bar) = %s   [%s]\n",
              g(x$att_ate_gap), fmt_ratio(att_ate$ratio),
              if (att_ate$appreciable) "appreciable" else "negligible"))
  cat("\n")

  gap_r <- label_term_status(x$gap_r_bar, x$se, x$appreciable_frac_se)
  cat("FSATT-ATT term (Proposition 2, retention channel):\n")
  cat(sprintf("  Gap_r_bar (mean across %d draw(s)) = %s   sd = %s\n",
              length(x$gap_m), g(x$gap_r_bar), g(x$gap_r_sd)))
  cat(sprintf("  |Gap_r_bar| / SE(tau_bar) = %s   [%s]\n",
              fmt_ratio(gap_r$ratio), if (gap_r$appreciable) "appreciable" else "negligible"))
  cat("\n")

  label <- if (gap_r$appreciable) {
    "FSATT"
  } else if (att_ate$appreciable) {
    "ATT"
  } else {
    "ATE"
  }

  cat(sprintf("--- Suggested label (heuristic, threshold = %s x SE, not from the paper) ---\n",
              g(x$appreciable_frac_se)))
  cat(sprintf("Report tau_bar as an estimate of: %s\n", label))
  if (gap_r$appreciable) {
    cat("The retention term is appreciable (or undecidable), so tau_bar can't be\n")
    cat("promoted past FSATT -- this holds regardless of the ATT-ATE term (see\n")
    cat("label_diagnostics()'s Details for why the chain stops here).\n")
  } else if (att_ate$appreciable) {
    cat("The retention term is negligible (FSATT ~ ATT), but the ATT-ATE term is not,\n")
    cat("so ATT is the most specific label the diagnostics support.\n")
  } else {
    cat("Both terms are negligible relative to SE(tau_bar): FSATT ~ ATT ~ ATE.\n")
  }
  cat("\n")
  cat("Per the paper (Sec 5.2), report these numbers as interpretive aids and label\n")
  cat("the estimate accordingly -- never subtract them off as a correction.\n")

  invisible(x)
}
