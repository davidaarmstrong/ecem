#' Run the congeniality diagnostics on a set of draws, together
#'
#' Runs the pre-matching flatness test ([flatness_test_XE()], once per
#' elicited or [regime()] covariate, over that covariate's full union
#' range), the congeniality test ([congeniality_test()] -- see its
#' documentation for why this replaced the earlier resampling-based
#' excess-variance design), and, by default, the retention-interaction
#' diagnostic ([cov_retention_per_draw()]) -- the three checks Section 6.3
#' of the paper uses to decide whether pooling via Rubin's rules is
#' trustworthy, and what to report if it is not.
#' [print.ecem_pooling_diagnostics()] reports the results together with
#' that decision logic spelled out.
#'
#' The congeniality test is run at `"low"`, `"mid"`, and `"high"` by
#' default (see `congeniality_position`), not just one fixed
#' specification: a test that stays quiet at a single position can simply
#' be sitting where drift is weak. Congeniality is treated as having
#' failed if *any* tested position rejects -- catching drift anywhere in
#' the elicited range is the whole point of testing more than one
#' position -- but testing three positions and rejecting if any one clears
#' `alpha` inflates the false-rejection rate under the true null relative
#' to `alpha`, since it's an uncorrected union of (correlated, but not
#' redundant) tests. Simulation checked this directly: at the null, the
#' raw three-position "any rejects" rule rejected roughly twice as often
#' as its nominal rate; a Bonferroni correction across however many
#' positions were actually tested brought the empirical rate back within
#' Monte Carlo noise of nominal. `congeniality_correction` (default
#' `"bonferroni"`) applies that correction to the combination rule; see
#' its own documentation for the alternatives.
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
#' that decision is left as a deliberate next step, and [existence_test()]
#' bootstraps fresh when you take it. There is nothing to precompute or
#' cache here: an earlier design let this function's bootstrap be shared
#' with a second, now-removed test, but [existence_test()] is the only
#' test left that resamples, so running its bootstrap inside this function
#' first would cost exactly the same `n_boot x M` rematches as running it
#' later inside [existence_test()] itself -- with the downside that it
#' would pay that cost on every call, including the (usually common)
#' specifications where congeniality doesn't reject and [existence_test()]
#' is never needed. This function is therefore cheap unconditionally: no
#' resampling, no rematching.
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
#' @param pooled Optional; the result of [pool_draws(draws)][pool_draws()].
#'   Computed automatically if not supplied.
#' @param run_congeniality Logical; whether to run [congeniality_test()].
#'   Defaults to `TRUE`.
#' @param congeniality_position `c("low", "mid", "high")` (the default) to
#'   run [congeniality_test()] at all three positions, each passed straight
#'   through to its `position` argument; supply a subset (even a single
#'   value, e.g. `"mid"`) to restrict which positions are tested. Testing
#'   all three is the point of the default -- drift that's weak at one
#'   fixed specification can be strong at another, so a single position
#'   silently passing is not by itself evidence of congeniality.
#' @param congeniality_vcov_type `"HC4"` (default) -- passed straight
#'   through to [congeniality_test()]'s `vcov_type`. Switch to `"HC4m"` if
#'   you see a `sandwich::vcovHC()` warning about numerically unstable
#'   covariances at high-leverage observations (common on real,
#'   unrounded data with a sparse treatment-by-bin cell) -- see
#'   [congeniality_test()]'s documentation for why.
#' @param congeniality_min_n_per_arm `1` (default, CEM's usual any-unit-
#'   per-arm retention rule) -- passed straight through to
#'   [congeniality_test()]'s `min_n_per_arm`. Raise this to also require a
#'   minimum sample size per treatment arm within each matched stratum
#'   before it's used for the congeniality test, e.g. to avoid
#'   `sandwich::vcovHC()` instability at singleton/near-singleton cells --
#'   see [congeniality_test()]'s documentation for why.
#' @param congeniality_correction `"bonferroni"` (default), `"holm"`, or
#'   `"none"` -- the multiple-testing correction applied, via
#'   [stats::p.adjust()], to each tested position's p-value before the
#'   "reject if any position rejects" combination rule is evaluated (see
#'   Details for why this matters: the uncorrected rule runs at roughly
#'   twice its nominal size). `"bonferroni"` divides `alpha` by the number
#'   of positions actually tested (so a single-position call is
#'   unaffected -- there's nothing to correct for). `"holm"` is uniformly
#'   at least as powerful (step-down rather than a flat threshold) while
#'   still controlling the family-wise rate, at the cost of a slightly
#'   less transparent per-position threshold. `"none"` restores the raw,
#'   uncorrected rule (not recommended when testing more than one
#'   position, but available for comparison). Ignored when only one
#'   position is tested, or when congeniality isn't run at all.
#' @param run_retention Logical; whether to compute the retention-
#'   interaction (`FSATT_m - ATT`) diagnostic. Defaults to `TRUE`; it is
#'   cheap (no bootstrapping, no rematching) and is exactly what the
#'   printed guidance points to when flatness and the congeniality test
#'   disagree.
#' @param alpha Significance level used to label each test's verdict in
#'   the printed output (and to choose which guidance
#'   [print.ecem_pooling_diagnostics()] prints). Defaults to `0.05`.
#'
#' @return An object of class `"ecem_pooling_diagnostics"`, a list with
#'   elements `flatness` (named list of [flatness_test_XE()] results, one
#'   per elicited/regime covariate), `congeniality` (`NULL` if
#'   `run_congeniality = FALSE` or nothing in `cutpoint_specs` is
#'   elicited, else a named list of [congeniality_test()] results, one per
#'   position actually tested, named by that position -- e.g.
#'   `congeniality$high` -- even when only one position was requested),
#'   `congeniality_correction` (the method actually used, echoing the
#'   argument), `congeniality_p_adjusted` (`NULL` alongside `congeniality`,
#'   else a named numeric vector, one adjusted p-value per tested
#'   position, `NA` wherever the raw p-value was `NA`), `congeniality_reject_any`
#'   (`NULL` alongside `congeniality`; else `NA` if every tested position's
#'   p-value was `NA`, else `TRUE`/`FALSE` -- whether any position's
#'   *adjusted* p-value cleared `alpha`, i.e. the actual combination-rule
#'   verdict), `retention` (`NULL` if
#'   `run_retention = FALSE`, else a list with `gap`, `mean`, `sd`),
#'   `pooled`, `M`, `K` (from [count_achievable_configs()], used by
#'   [print.ecem_pooling_diagnostics()] to decide whether recommending
#'   exact enumeration is actually practical for pooling precision),
#'   `exact`, and `alpha`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- simulate_population(N = 1000, heterogeneous = TRUE)
#' specs <- list(age = list(c(25, 33), c(60, 68)), educ = c(12, 16))
#' draws <- run_M_draws(pop, "D", "Y", specs, M = 15)
#' diag <- pooling_diagnostics(pop, draws)
#' diag
#' }
#'
#' @export
pooling_diagnostics <- function(data, draws, treat_var = NULL, outcome_var = NULL,
                                 cutpoint_specs = NULL,
                                 pooled = NULL,
                                 run_congeniality = TRUE,
                                 congeniality_position = c("low", "mid", "high"),
                                 congeniality_vcov_type = "HC4",
                                 congeniality_min_n_per_arm = 1,
                                 congeniality_correction = c("bonferroni", "holm", "none"),
                                 run_retention = TRUE,
                                 alpha = 0.05) {
  congeniality_position <- match.arg(congeniality_position, choices = c("low", "mid", "high"),
                                      several.ok = TRUE)
  congeniality_correction <- match.arg(congeniality_correction, choices = c("bonferroni", "holm", "none"))
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
  congeniality_p_adjusted <- NULL
  congeniality_reject_any <- NULL
  if (run_congeniality && length(xe_ranges) > 0) {
    congeniality <- lapply(congeniality_position, function(pos) {
      congeniality_test(data, treat_var, outcome_var, cutpoint_specs,
                         position = pos,
                         vcov_type = congeniality_vcov_type,
                         min_n_per_arm = congeniality_min_n_per_arm)
    })
    names(congeniality) <- congeniality_position

    ## The "reject if any tested position rejects" combination rule is an
    ## uncorrected union of (correlated, but not redundant) tests unless
    ## adjusted here -- see this function's Details for the simulation
    ## evidence. p.adjust(method = "none") is a documented no-op, so this
    ## same code path handles all three congeniality_correction choices
    ## without a special case. NA raw p-values (degenerate matched sample
    ## at that position) are left NA rather than passed to p.adjust().
    cg_p <- vapply(congeniality, function(r) r$p_value, numeric(1))
    congeniality_p_adjusted <- stats::setNames(rep(NA_real_, length(cg_p)), names(cg_p))
    non_na <- !is.na(cg_p)
    if (any(non_na)) {
      congeniality_p_adjusted[non_na] <- stats::p.adjust(cg_p[non_na], method = congeniality_correction)
    }
    congeniality_reject_any <- if (!any(non_na)) {
      NA
    } else {
      isTRUE(any(congeniality_p_adjusted[non_na] < alpha))
    }
  }

  retention <- NULL
  if (run_retention) {
    treat_idx <- which(data[[treat_var]] == 1)
    tau_i <- pooled_unit_tau(draws, nrow(data))
    gap <- cov_retention_per_draw(draws, tau_i, treat_idx)
    retention <- list(gap = gap, mean = mean(gap, na.rm = TRUE), sd = stats::sd(gap, na.rm = TRUE))
  }

  ## Cheap (no bootstrapping, no matching) and lets print.ecem_pooling_diagnostics()
  ## give concrete, K-aware advice about POOLING precision (exact
  ## enumeration vs. a larger M) -- the congeniality test itself no longer
  ## depends on K or M at all, so this is reported for that separate
  ## question, not as a lever on the congeniality verdict. count_only =
  ## TRUE is essential here, not just an optimization: this function only
  ## ever reads $K, but count_achievable_configs()'s default path
  ## materializes every achievable config's cutpoints/weight to get it --
  ## on real, unrounded covariates with many distinct values, that list
  ## can run into the millions or billions of entries and make an
  ## otherwise-cheap diagnostic call hang indefinitely. count_only = TRUE
  ## computes K directly from gap counts instead.
  K <- count_achievable_configs(data, cutpoint_specs, count_only = TRUE)$K

  out <- list(
    flatness                = flatness,
    congeniality            = congeniality,
    congeniality_correction = congeniality_correction,
    congeniality_p_adjusted = congeniality_p_adjusted,
    congeniality_reject_any = congeniality_reject_any,
    retention               = retention,
    pooled                  = pooled,
    M                       = M,
    K                       = K,
    exact                   = isTRUE(attr(draws, "exact")),
    alpha                   = alpha
  )
  class(out) <- "ecem_pooling_diagnostics"
  out
}
