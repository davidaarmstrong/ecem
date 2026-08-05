#' Label diagnostics: FSATT, ATT, or ATE?
#'
#' Computes the two covariance diagnostics from Section 4 that separate the
#' three estimands CEM's estimator can be read as targeting -- FSATT
#' (feasible sample ATT, defined on the retained sample), ATT, and ATE --
#' and reports them together. Mirrors [pooling_diagnostics()] in shape, but
#' for the FSATT -> ATT -> ATE labeling question rather than pooling
#' validity; Section 5.3 treats the two as parallel, always-run checks
#' rather than one gating the other, so this function does not depend on
#' [pooling_diagnostics()] or vice versa, and either can be run first.
#'
#' The ATT-ATE term is `Cov(tau_hat(X), p_hat(X))` ([cov_att_ate()]), and
#' the equivalent ATT - ATE gap (Proposition 1) is that covariance divided
#' by `Pr(D = 1)`. The FSATT-ATT term is `Gap_m`
#' ([cov_retention_per_draw()]), the full ratio from Proposition 2 (not the
#' covariance alone -- the covariance understates the gap whenever
#' retention is appreciably below complete), computed once per draw and
#' reported as its mean and spread across the `M` (or `K`) draws: the paper
#' treats stability across draws as itself informative -- a stable, small
#' \eqn{\overline{\mathrm{Gap}}_r} means the retention channel isn't
#' driving any unexplained between-draw variance, while a volatile one is
#' the retention-interaction channel's own signature.
#'
#' **On the suggested label.** The paper is explicit that it gives no
#' numeric threshold for "appreciable" versus "negligible" (Section 5.2)
#' and recommends against ever using these diagnostics to numerically
#' correct the point estimate -- they are for labeling and disclosure only,
#' never subtraction. [print.ecem_label_diagnostics()] nonetheless prints a
#' *suggested* label, clearly flagged there as a heuristic this package
#' proposes rather than something the paper specifies: a term is called
#' "appreciable" if its magnitude is at least `appreciable_frac_se` times
#' the pooled estimate's own standard error (`sqrt(T)`, from
#' [pool_draws()]) -- i.e., large enough that it could plausibly be
#' distinguished from the sampling noise already in `tau_bar`, not large in
#' any absolute sense.
#'
#' The paper's own rule covers three of the four possible combinations: an
#' FSATT label if both terms are appreciable, an ATT label if only the
#' retention term is negligible, an ATE label if both are negligible. It
#' does not say what to do if the retention term is appreciable but the
#' ATT-ATE term is negligible. This function extends the rule to that case
#' by treating FSATT -> ATT -> ATE as a chain rather than two independent
#' conditions: you cannot claim ATT -- and so there is nothing to gain by
#' also asking whether ATT and ATE agree -- unless the FSATT-ATT link is
#' itself negligible first. So an appreciable retention term labels the
#' estimate FSATT regardless of the ATT-ATE term.
#'
#' @param data A data frame, the same one used to produce `draws`.
#' @param draws An object of class `"ecem_draws"`, as returned by
#'   [run_M_draws()].
#' @param treat_var,cutpoint_specs `NULL` (the default) to recover these
#'   from `draws`'s attributes; supply them explicitly to override, or if
#'   `draws` doesn't carry them. Unlike [pooling_diagnostics()], this
#'   function has no use for `outcome_var` -- everything here is computed
#'   from `draws`'s already-fitted `tau_hat`/retention information, not
#'   from the raw outcome column.
#' @param covariates Character vector of covariate names for the
#'   propensity model ([fit_propensity()]). `NULL` (the default) uses every
#'   non-excluded covariate in `cutpoint_specs` -- i.e., exactly the
#'   covariates entering the matching. Pass a different (larger) set
#'   explicitly to condition the propensity model on covariates that were
#'   excluded from matching (e.g. `X`s elicited as `NULL` in the spec).
#' @param pooled Optional; the result of [pool_draws(draws)][pool_draws()].
#'   Computed automatically if not supplied. Used only for its `T`
#'   (`sqrt(T)` is the yardstick the suggested label compares each term
#'   to).
#' @param appreciable_frac_se Numeric; a term is flagged "appreciable" in
#'   the suggested label if its magnitude is at least this fraction of
#'   `sqrt(pooled$T)`. Defaults to `0.5`. Purely a printing heuristic (see
#'   Details) -- change it freely, or ignore the suggested label entirely
#'   and read the reported numbers yourself.
#'
#' @return An object of class `"ecem_label_diagnostics"` (see
#'   [print.ecem_label_diagnostics()]): a list with elements `cov_tau_p`
#'   (`Cov(tau_hat(X), p_hat(X))`), `att_ate_gap` (`cov_tau_p / Pr(D=1)`,
#'   Proposition 1's ATT - ATE estimate), `gap_m` (per-draw
#'   \eqn{\widehat{\mathrm{Gap}}_m}, Proposition 2), `gap_r_bar` and
#'   `gap_r_sd` (its mean and spread across draws), `se` (`sqrt(pooled$T)`),
#'   `pooled`, `covariates` (the set actually used), and
#'   `appreciable_frac_se`.
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
label_diagnostics <- function(data, draws, treat_var = NULL, cutpoint_specs = NULL,
                               covariates = NULL, pooled = NULL,
                               appreciable_frac_se = 0.5) {
  if (is.null(treat_var))      treat_var      <- attr(draws, "treat_var")
  if (is.null(cutpoint_specs)) cutpoint_specs <- attr(draws, "cutpoint_specs")

  if (is.null(treat_var) || is.null(cutpoint_specs)) {
    stop(
      "treat_var and cutpoint_specs could not both be recovered from ",
      "`draws`. This happens if `draws` didn't come from run_M_draws() ",
      "(e.g. it was assembled by hand or predates this package's tracking ",
      "of these attributes) -- pass whichever of treat_var/cutpoint_specs ",
      "is missing explicitly."
    )
  }

  if (is.null(covariates)) {
    covariates <- names(Filter(Negate(is.null), cutpoint_specs))
  }

  if (is.null(pooled)) {
    pooled <- pool_draws(draws)
  }

  treat_idx <- which(data[[treat_var]] == 1)
  tau_i <- pooled_unit_tau(draws, nrow(data))
  p_hat <- fit_propensity(data, treat_var, covariates)

  cov_tau_p <- cov_att_ate(tau_i, p_hat, treat_idx)
  pr_treat  <- mean(data[[treat_var]] == 1)
  att_ate_gap <- cov_tau_p / pr_treat

  gap_m     <- cov_retention_per_draw(draws, tau_i, treat_idx)
  gap_r_bar <- mean(gap_m, na.rm = TRUE)
  gap_r_sd  <- stats::sd(gap_m, na.rm = TRUE)

  out <- list(
    cov_tau_p           = cov_tau_p,
    att_ate_gap         = att_ate_gap,
    gap_m               = gap_m,
    gap_r_bar           = gap_r_bar,
    gap_r_sd            = gap_r_sd,
    se                  = sqrt(pooled$T),
    pooled              = pooled,
    covariates          = covariates,
    appreciable_frac_se = appreciable_frac_se
  )
  class(out) <- "ecem_label_diagnostics"
  out
}
