#' Exact enumeration of achievable coarsenings
#'
#' @description
#' `cut()`-based coarsening only changes at a covariate's own observed
#' values, so the number of practically distinct outcomes achievable by
#' varying a cutpoint within an elicited range is finite: if k distinct
#' observed values fall in `[lo, hi]`, there are exactly k+1 intervals that
#' produce a distinct binning (before the first value, between each pair,
#' after the last). Two draws landing in the same interval produce an
#' *identical* matched sample, not merely a similar one, since the data
#' doesn't change between draws. Increasing `M` past this ceiling doesn't
#' buy more power -- it just re-estimates a fixed quantity more precisely.
#' [count_achievable_configs()] reports the ceiling before you commit to an
#' `M`; [enumerate_configs()] / [run_enumerated_draws()] replace Monte
#' Carlo drawing with exact computation once that ceiling is small enough
#' to be worth enumerating directly. [run_M_draws()]'s `exact_if_K_leq`
#' argument does this automatically.
#'
#' @name exact-enumeration
NULL

## Achievable gap midpoints and widths for one range of one covariate.
## Internal to enumerate_var_configs().
enumerate_gaps <- function(x, range_) {
  vals   <- sort(unique(x[x >= range_[1] & x <= range_[2]]))
  bounds <- unique(c(range_[1], vals, range_[2]))
  n      <- length(bounds) - 1
  if (n < 1) {
    return(data.frame(rep = mean(range_), width = diff(range_)))
  }
  reps   <- (bounds[-length(bounds)] + bounds[-1]) / 2
  widths <- diff(bounds)
  data.frame(rep = reps, width = widths)
}

## Enumerates every achievable coarsening for ONE covariate, given its
## spec. Returns NULL if the covariate does not vary across draws
## (excluded, exact, or fixed) -- such covariates contribute no
## combinatorial factor to K and are left untouched in the spec passed
## through to run_draw().
##
## For a plain elicited covariate (list of ranges), this is the full
## Cartesian product across its own cutpoint slots: one achievable config
## per combination of gaps, one per slot.
##
## For a regime() covariate, this is the UNION of each regime's own
## Cartesian product -- regimes are never crossed with each other, only
## weighted (by the regime's own prior weight times the interval widths
## within it) and concatenated. A single achievable config is therefore
## always entirely one regime, matching how draw_cutpoints_for_var() draws
## it: pick the regime first, then draw cutpoints only within it. The
## covariate's own achievable-config count is thus a SUM across regimes,
## not a product -- crossing happens within a regime, never across
## regimes. Unrelated covariates still cross against this one as a normal
## factor in the product computed by count_achievable_configs().
##
## Internal to count_achievable_configs() / enumerate_configs().
enumerate_var_configs <- function(data, var, spec) {
  x <- data[[var]]

  if (is.null(spec) || identical(spec, "exact")) {
    return(NULL)                                  # excluded / exact: fixed across draws
  }

  cartesian_product <- function(ranges) {
    gaps <- lapply(ranges, function(r) enumerate_gaps(x, r))
    grid <- expand.grid(lapply(gaps, function(g) seq_len(nrow(g))), KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(grid)), function(i) {
      cuts <- vapply(seq_along(gaps), function(j) gaps[[j]]$rep[grid[i, j]], numeric(1))
      w    <- prod(vapply(seq_along(gaps), function(j) gaps[[j]]$width[grid[i, j]], numeric(1)))
      list(cutpoints = sort(cuts), weight = w)
    })
  }

  if (inherits(spec, "cem_regime_spec")) {
    out <- list()
    for (rname in names(spec$regimes)) {
      within <- cartesian_product(spec$regimes[[rname]])
      for (cfg in within) {
        cfg$weight <- cfg$weight * spec$weights[[rname]]
        cfg$regime <- rname
        ## Tag the realized cutpoints themselves with which regime
        ## produced them. enumerate_configs() assigns cfg$cutpoints
        ## straight into spec_i[[var]] as a plain numeric vector, which
        ## strips the cem_regime_spec class -- without this attribute,
        ## draw_cutpoints_for_var() would re-resolve that vector as an
        ## ordinary "fixed" spec and $matched$regimes would come back NA
        ## for every exactly-enumerated draw, even though it really is one
        ## regime or the other.
        attr(cfg$cutpoints, "regime") <- rname
        out[[length(out) + 1]] <- cfg
      }
    }
    return(out)                                    # sum across regimes, not a product
  }

  if (is.list(spec)) {                              # ordinary elicited covariate
    return(lapply(cartesian_product(spec), function(cfg) {
      cfg$regime <- NA_character_
      cfg
    }))
  }

  NULL                                              # fixed (X_F): doesn't vary
}

## Achievable-config COUNT for one covariate, without materializing the
## per-config cutpoints/weight list enumerate_var_configs() builds --
## mirrors its logic exactly (product of per-range gap counts for a plain
## elicited covariate; sum across regimes of each regime's own product for
## a regime() covariate) but only ever computes integers. enumerate_gaps()
## itself is already cheap (O(N log N) per range, just sorting the
## observed values in range); what's expensive in enumerate_var_configs()
## is the per-row lapply() that follows it to build one list entry per
## achievable combination, which this skips entirely. Real, unrounded
## covariates can put the achievable count in the millions or billions per
## variable -- fine to know as a number, infeasible to ever materialize as
## a list -- so count_achievable_configs()'s count_only path (the one
## pooling_diagnostics() uses, since it only ever reads $K) calls this
## instead of enumerate_var_configs() for exactly that reason.
##
## Internal to count_achievable_configs().
count_var_configs <- function(data, var, spec) {
  x <- data[[var]]

  if (is.null(spec) || identical(spec, "exact")) {
    return(NULL)                                  # excluded / exact: fixed across draws
  }

  range_gap_count <- function(r) nrow(enumerate_gaps(x, r))

  if (inherits(spec, "cem_regime_spec")) {
    total <- 0
    for (rname in names(spec$regimes)) {
      total <- total + prod(vapply(spec$regimes[[rname]], range_gap_count, integer(1)))
    }
    return(total)                                  # sum across regimes, not a product
  }

  if (is.list(spec)) {                              # ordinary elicited covariate
    return(prod(vapply(spec, range_gap_count, integer(1))))
  }

  NULL                                              # fixed (X_F): doesn't vary
}

#' Count the achievable-configuration ceiling K
#'
#' See `?"exact-enumeration"` for the underlying argument. Reports how many
#' practically distinct coarsenings are achievable given `cutpoint_specs`
#' and the observed values in `data` -- the ceiling past which a larger `M`
#' stops buying additional power and starts just re-estimating the same
#' fixed quantity more precisely.
#'
#' `count_only = FALSE` (the default) additionally materializes and
#' returns every achievable config's realized cutpoints/weight in
#' `per_variable` -- this is what [enumerate_configs()] needs to actually
#' build them, but it means the *count* isn't free: real, unrounded
#' covariates can have thousands of distinct values per elicited range, and
#' the achievable count for even one covariate with more than one elicited
#' range is the *product* of those per-range counts, easily reaching
#' millions or billions. Materializing a list that large is what's slow,
#' not computing its length. Set `count_only = TRUE` to get just `K` (and
#' `counts_by_variable`) computed directly from each range's gap count,
#' without ever building the per-config list -- this stays fast regardless
#' of how large `K` turns out to be, which is why
#' [pooling_diagnostics()] uses it: it only ever reads `$K`, to report
#' pooling precision, never the enumerated configs themselves.
#'
#' @inheritParams elicit_and_match
#' @param count_only Logical; if `TRUE`, compute `K` (and
#'   `counts_by_variable`) directly from gap counts without materializing
#'   `per_variable`. Defaults to `FALSE` for backward compatibility (this
#'   function's output has always included `per_variable`), but callers
#'   who only need `K` should pass `TRUE` -- see Details.
#'
#' @return A list with elements `K` (the total achievable-configuration
#'   count, the product across covariates of each covariate's own count),
#'   `per_variable` (named list, one entry per varying covariate, each a
#'   list of achievable configs with their cutpoints/weight/regime --
#'   omitted when `count_only = TRUE`), and, if any covariate varies,
#'   `counts_by_variable` (named numeric vector).
#'
#' @export
count_achievable_configs <- function(data, cutpoint_specs, count_only = FALSE) {
  if (count_only) {
    counts <- list()
    for (v in names(cutpoint_specs)) {
      n <- count_var_configs(data, v, cutpoint_specs[[v]])
      if (!is.null(n)) counts[[v]] <- n
    }
    if (length(counts) == 0) {
      return(list(K = 1L))
    }
    counts <- unlist(counts)
    return(list(K = prod(counts), counts_by_variable = counts))
  }

  per_var <- list()
  for (v in names(cutpoint_specs)) {
    cfgs <- enumerate_var_configs(data, v, cutpoint_specs[[v]])
    if (!is.null(cfgs)) per_var[[v]] <- cfgs
  }
  if (length(per_var) == 0) {
    return(list(K = 1L, per_variable = per_var))
  }
  counts <- vapply(per_var, length, integer(1))
  list(K = prod(counts), per_variable = per_var, counts_by_variable = counts)
}

#' Enumerate every achievable configuration exactly
#'
#' See `?"exact-enumeration"` for the underlying argument. Returns `K`
#' exact `cutpoint_specs` (each with elicited ranges replaced by a single
#' realized, fixed cutpoint vector -- [draw_cutpoints_for_var()] already
#' treats a plain numeric vector as fixed, so [run_draw()] needs no changes
#' to consume these) plus their exact probability weights, proportional to
#' interval width under the continuous elicited prior (and, for a
#' [regime()] covariate, to the regime's own prior weight as well). Also
#' returns, per config, which regime (if any) each varying regime
#' covariate landed in.
#'
#' @inheritParams elicit_and_match
#'
#' @return A list with elements `specs` (list of `K` realized
#'   `cutpoint_specs`), `weights` (numeric, summing to 1), and `regimes`
#'   (list of named lists giving the drawn regime, if any, per config).
#'
#' @export
enumerate_configs <- function(data, cutpoint_specs) {
  info    <- count_achievable_configs(data, cutpoint_specs)
  per_var <- info$per_variable
  if (length(per_var) == 0) {
    return(list(specs = list(cutpoint_specs), weights = 1, regimes = list(list())))
  }

  grids  <- lapply(per_var, function(cfgs) seq_along(cfgs))
  combos <- expand.grid(grids, KEEP.OUT.ATTRS = FALSE)
  names(combos) <- names(per_var)

  n_configs <- nrow(combos)
  specs   <- vector("list", n_configs)
  weights <- numeric(n_configs)
  regimes <- vector("list", n_configs)

  for (i in seq_len(n_configs)) {
    spec_i <- cutpoint_specs
    w <- 1
    reg_i <- list()
    for (v in names(per_var)) {
      cfg <- per_var[[v]][[combos[[v]][i]]]
      spec_i[[v]] <- cfg$cutpoints
      w <- w * cfg$weight
      if (!is.na(cfg$regime)) reg_i[[v]] <- cfg$regime
    }
    specs[[i]]   <- spec_i
    weights[i]   <- w
    regimes[[i]] <- reg_i
  }

  weights <- weights / sum(weights)
  list(specs = specs, weights = weights, regimes = regimes)
}

#' Run every achievable configuration exactly
#'
#' See `?"exact-enumeration"`. Runs [run_draw()] once per config returned
#' by [enumerate_configs()]. Prefer [run_M_draws()] with `exact_if_K_leq`
#' set for the common case of "enumerate exactly if small enough, else fall
#' back to Monte Carlo" -- call this directly only if you specifically want
#' exact enumeration regardless of how large `K` is.
#'
#' @inheritParams elicit_and_match
#' @param outcome_var Character; name of the outcome column in `data`.
#' @inheritParams fit_effect
#'
#' @return A list with elements `results` (list of `K` per-draw records, as
#'   from [run_draw()]), `weights`, `regimes`, and `K`.
#'
#' @export
run_enumerated_draws <- function(data, treat_var, outcome_var, cutpoint_specs,
                                  estimator = c("regression", "mean_diff"), covariates = NULL) {
  estimator <- match.arg(estimator)
  enum <- enumerate_configs(data, cutpoint_specs)
  results <- lapply(enum$specs, function(spec_i) {
    run_draw(data, treat_var, outcome_var, spec_i, estimator = estimator, covariates = covariates)
  })
  list(results = results, weights = enum$weights, regimes = enum$regimes, K = length(enum$specs))
}

#' Pool exactly-enumerated configurations
#'
#' Exact-weight pooling: no `(1 + 1/K)` inflation, because that correction
#' exists specifically to account for not having sampled the full
#' population of draws -- which, under exact enumeration, has now been
#' sampled in full. Both `tau_bar` and `B` are exact, not Monte Carlo
#' estimates, so `T = Wbar + B` with no further adjustment. Prefer
#' [pool_draws()], which dispatches to this automatically when appropriate.
#'
#' @param results A list of per-draw records, as returned by
#'   [run_enumerated_draws()]'s `results` element (or [run_M_draws()]'s
#'   output when it enumerated exactly).
#' @param weights Numeric vector of weights, one per element of `results`,
#'   as returned by [enumerate_configs()]'s `weights` element.
#'
#' @return An object of class `"ecem_pooled"` (see [print.ecem_pooled()]):
#'   a list with elements `tau_bar`, `Wbar`, `B`, `T`, and `K`. Unlike
#'   [pool_draws()]'s output, `lambda_hat`, `df`, and `exact` are not set
#'   here; [print.ecem_pooled()] fills them in at print time (`lambda_hat
#'   = B / T`, `df = NA`, exact inferred from the presence of `K`).
#'
#' @export
pool_rubins_rules_exact <- function(results, weights) {
  tau <- vapply(results, function(d) d$tau_hat, numeric(1))
  v   <- vapply(results, function(d) d$var_hat, numeric(1))
  ok  <- !is.na(tau) & !is.na(v)
  tau <- tau[ok]; v <- v[ok]; w <- weights[ok]
  w <- w / sum(w)

  tau_bar <- sum(w * tau)
  Wbar    <- sum(w * v)
  B       <- sum(w * (tau - tau_bar)^2)
  Tvar    <- Wbar + B

  out <- list(tau_bar = tau_bar, Wbar = Wbar, B = B, T = Tvar, K = length(tau))
  class(out) <- c("ecem_pooled", "list")
  out
}
