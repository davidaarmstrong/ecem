#' Pool Monte Carlo draws via Rubin's rules
#'
#' Flat, single-level Rubin's rules pooling across a set of Monte Carlo
#' draws produced by [run_M_draws()]. Assumes equal weight across draws --
#' if `draws` came from exact enumeration (unequal weights), use
#' [pool_draws()] instead, which dispatches correctly either way.
#'
#' @param draws A list of per-draw records, as returned by [run_M_draws()]
#'   (with `exact_if_K_leq` left at its default `NULL`) or by repeated
#'   calls to [run_draw()].
#'
#' @return An object of class `"ecem_pooled"` (see [print.ecem_pooled()]):
#'   a list with elements `tau_bar`, `Wbar`, `B`, `T`, `lambda_hat`, `df`,
#'   and `M`.
#'
#' @section Caveat on `df`:
#' The Barnard-Rubin small-sample degrees-of-freedom adjustment needs a
#' "complete-data" degrees of freedom (`df_com`) from the estimator used
#' within each draw. This function uses `mean(n_used) - 2` as a
#' placeholder; replace it with something specific to your estimator before
#' relying on `df` for actual inference.
#'
#' @export
pool_rubins_rules <- function(draws) {
  tau    <- vapply(draws, function(d) d$tau_hat, numeric(1))
  v      <- vapply(draws, function(d) d$var_hat, numeric(1))
  n_used <- vapply(draws, function(d) d$n_used,  numeric(1))

  ok <- !is.na(tau) & !is.na(v)
  tau <- tau[ok]; v <- v[ok]; n_used <- n_used[ok]
  M <- length(tau)

  tau_bar    <- mean(tau)
  Wbar       <- mean(v)
  B          <- stats::var(tau)
  Tvar       <- Wbar + (1 + 1 / M) * B
  lambda_hat <- (1 + 1 / M) * B / Tvar

  df_com <- mean(n_used) - 2
  gamma  <- (1 + 1 / M) * B / Tvar
  df_old <- (M - 1) / gamma^2
  df_obs <- (df_com + 1) / (df_com + 3) * df_com * (1 - gamma)
  df_adj <- 1 / (1 / df_old + 1 / df_obs)

  out <- list(tau_bar = tau_bar, Wbar = Wbar, B = B, T = Tvar,
              lambda_hat = lambda_hat, df = df_adj, M = M)
  class(out) <- c("ecem_pooled", "list")
  out
}

#' Pool a set of draws, dispatching to the correct formula
#'
#' Pools a `draws` object from [run_M_draws()], dispatching to the exact
#' ([pool_rubins_rules_exact()]) or Monte Carlo ([pool_rubins_rules()])
#' formula depending on how the draws were produced. This is the function
#' to call by default -- it returns the same set of fields either way, so
#' nothing downstream needs to know or care which path `run_M_draws()`
#' took.
#'
#' `df` is a Barnard-Rubin quantity motivated by finite-`M`
#' multiple-imputation uncertainty: how much is left to learn from drawing
#' more configurations. Once every achievable configuration has actually
#' been enumerated there is none left, so `df` comes back `NA` under exact
#' pooling -- `lambda_hat` (still computed as `B / T`, the share of total
#' variance coming from between-draw variation) and `B`/`T` themselves
#' (exact rather than Monte Carlo estimates in that case) remain
#' meaningful and are what [print.ecem_pooled()] falls back to a normal,
#' rather than a `t`, reference distribution for.
#'
#' @param draws A list of per-draw records, as returned by [run_M_draws()].
#'
#' @return An object of class `"ecem_pooled"` (see [print.ecem_pooled()]):
#'   a list with elements `tau_bar`, `Wbar`, `B`, `T`, `lambda_hat`, `df`,
#'   and `exact` (logical, whether exact pooling was used).
#'
#' @export
pool_draws <- function(draws) {
  if (isTRUE(attr(draws, "exact"))) {
    out <- pool_rubins_rules_exact(draws, attr(draws, "weights"))
    out$lambda_hat <- out$B / out$T
    out$df         <- NA_real_
    out$exact      <- TRUE
    class(out)     <- c("ecem_pooled", "list")
    out
  } else {
    out <- pool_rubins_rules(draws)
    out$exact  <- FALSE
    out
  }
}
