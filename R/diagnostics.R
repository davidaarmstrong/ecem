#' Build the retention matrix R_m(X_i)
#'
#' For every treated unit and every draw, records whether that unit was
#' retained (matched to at least one control in its stratum) in that draw.
#'
#' @param draws A list of per-draw records, as returned by [run_M_draws()].
#' @param treat_idx Integer vector of row indices for the treated units
#'   (e.g. `which(data$treat_var == 1)`).
#'
#' @return An integer matrix with one row per treated unit and one column
#'   per draw, entries 0/1.
#'
#' @export
retention_matrix <- function(draws, treat_idx) {
  M <- length(draws)
  R <- matrix(0L, nrow = length(treat_idx), ncol = M,
              dimnames = list(as.character(treat_idx), paste0("m", seq_len(M))))
  for (m in seq_len(M)) {
    retained_treated <- intersect(draws[[m]]$matched$retained_idx, treat_idx)
    if (length(retained_treated) > 0) {
      R[as.character(retained_treated), m] <- 1L
    }
  }
  R
}

#' Fit a propensity score on the full (uncoarsened) data
#'
#' Estimates `p(x)` once on the full data via logistic regression. Does not
#' depend on the coarsening draw, unlike `r_m(x)`.
#'
#' @param data A data frame containing `treat_var` and every covariate in
#'   `covariates`.
#' @param treat_var Character; name of the 0/1 treatment indicator column.
#' @param covariates Character vector of covariate names to include in the
#'   propensity model.
#'
#' @return A numeric vector of fitted propensity scores, one per row of
#'   `data`.
#'
#' @export
fit_propensity <- function(data, treat_var, covariates) {
  form <- stats::as.formula(paste(treat_var, "~", paste(covariates, collapse = " + ")))
  fit  <- stats::glm(form, data = data, family = stats::binomial())
  as.numeric(stats::fitted(fit))
}

#' Pool each treated unit's tau_hat across the draws that retained it
#'
#' Each treated unit's `tau_hat`, averaged across the draws in which it was
#' retained. Used as the `tau(x)` surface for both covariance diagnostics
#' ([cov_att_ate()] and [cov_retention_per_draw()]) -- deliberately the
#' *pooled* surface, not a per-draw one, since `tau(x)` is a fixed
#' population quantity while `r_m(x)` is what actually varies by draw.
#'
#' @param draws A list of per-draw records, as returned by [run_M_draws()].
#' @param n Total number of rows in the original data (i.e. `nrow(data)`).
#'
#' @return A numeric vector of length `n`, `NA` for units never retained in
#'   any draw.
#'
#' @export
pooled_unit_tau <- function(draws, n) {
  sums   <- numeric(n)
  counts <- integer(n)
  for (d in draws) {
    if (length(d$unit_idx) > 0) {
      sums[d$unit_idx]   <- sums[d$unit_idx] + d$unit_tau_hat
      counts[d$unit_idx] <- counts[d$unit_idx] + 1L
    }
  }
  ifelse(counts > 0, sums / counts, NA_real_)
}

#' ATT-ATE labeling diagnostic: Cov(tau_hat(X), p_hat(X))
#'
#' A single number (unlike [cov_retention_per_draw()]), since `p(x)` does
#' not depend on the coarsening draw. A nonzero estimate does not
#' invalidate the pooled estimate; it indicates the pooled FSATT should be
#' labeled ATT rather than ATE (see the paper's Proposition on the ATT-ATE
#' gap).
#'
#' @param tau_i Numeric vector, the pooled per-unit `tau(x)` surface from
#'   [pooled_unit_tau()].
#' @param p_hat Numeric vector of propensity scores from
#'   [fit_propensity()], same length and indexing as `tau_i`.
#' @param treat_idx Integer vector of row indices for the treated units.
#'
#' @return A single numeric covariance estimate.
#'
#' @export
cov_att_ate <- function(tau_i, p_hat, treat_idx) {
  use <- treat_idx[!is.na(tau_i[treat_idx])]
  stats::cov(tau_i[use], p_hat[use])
}

#' FSATT-ATT labeling diagnostic: the per-draw retention gap
#'
#' Computes \eqn{\widehat{Gap}_m = Cov(\hat\tau(X), \hat r_m(X) \mid D=1) /
#' \widehat{Pr}(R_m = 1 \mid D = 1)} for each draw -- the full identity for
#' \eqn{FSATT_m - ATT} (see the paper's Proposition on this gap), not just
#' its numerator. Dividing by the draw's own treated-side retention rate
#' matters: a small, stable covariance can still hide a large, volatile gap
#' if retention rates themselves swing across draws, which is exactly what
#' the retention-interaction channel predicts.
#'
#' Report the mean across draws alongside its spread: under stability, the
#' `M` (or `K`) per-draw estimates should cluster tightly, and substantial
#' spread among them is itself informative -- it is the observable
#' signature of the retention-interaction channel.
#'
#' @inheritParams cov_att_ate
#' @param draws A list of per-draw records, as returned by [run_M_draws()].
#'
#' @return A numeric vector, one estimate per draw, `NA` for any draw whose
#'   treated-side retention rate among used units is zero or undefined.
#'
#' @export
cov_retention_per_draw <- function(draws, tau_i, treat_idx) {
  vapply(draws, function(d) {
    r_m <- as.integer(treat_idx %in% d$matched$retained_idx)
    use <- !is.na(tau_i[treat_idx])
    cov_r <- stats::cov(tau_i[treat_idx][use], r_m[use])
    pr_retain <- mean(r_m[use])
    if (!is.finite(pr_retain) || pr_retain <= 0) {
      return(NA_real_)
    }
    cov_r / pr_retain
  }, numeric(1))
}
