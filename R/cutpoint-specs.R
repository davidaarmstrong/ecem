#' Specify competing coarsening regimes for one covariate
#'
#' Bundles two or more named, mutually exclusive elicited cutpoint sets for
#' a single covariate -- e.g. a lifecycle-stage theory and a
#' generational-cohort theory of where age should be coarsened. Pass the
#' result as that covariate's entry in `cutpoint_specs`.
#'
#' Each draw first picks one whole regime (by `weights`, uniform by
#' default), then draws cutpoints only within that regime, exactly as it
#' would for a plain elicited covariate. Cutpoints from different regimes
#' never appear together in the same draw -- a draw is entirely one regime
#' or entirely another, never a mix -- but other, unrelated covariates
#' still vary independently against whichever regime was drawn.
#'
#' @param ... Two or more named arguments, each a list of numeric ranges
#'   (the same format used for a plain elicited covariate), e.g.
#'   `lifecycle = list(c(18, 25), c(35, 45), c(60, 70))`.
#' @param weights Optional numeric vector of prior weights, one per regime,
#'   in the same order as `...`. Defaults to uniform. Renormalized to sum
#'   to 1 regardless of what is supplied.
#'
#' @return An object of class `"cem_regime_spec"`, for use as one entry in
#'   a `cutpoint_specs` list.
#'
#' @examples
#' age_regimes <- regime(
#'   lifecycle    = list(c(22, 28), c(40, 48), c(62, 68)),
#'   generational = list(c(24, 28), c(42, 46), c(58, 62), c(76, 80))
#' )
#'
#' @export
regime <- function(..., weights = NULL) {
  regimes <- list(...)
  if (length(regimes) < 2) {
    stop("regime() needs at least two named regimes to choose between.")
  }
  if (is.null(names(regimes)) || any(names(regimes) == "")) {
    stop("Every regime passed to regime() must be named.")
  }
  if (is.null(weights)) {
    weights <- rep(1 / length(regimes), length(regimes))
  }
  stopifnot(length(weights) == length(regimes))
  names(weights) <- names(regimes)
  structure(list(regimes = regimes, weights = weights / sum(weights)),
            class = "cem_regime_spec")
}

#' Resolve one covariate's cutpoint spec for a single draw
#'
#' Internal building block of [elicit_and_match()]. A cutpoint spec for one
#' covariate is one of five things:
#'
#' - a list of numeric ranges, e.g. `list(c(25, 33), c(60, 68))`:
#'   **elicited** (`X_E`). Draws one value uniformly, continuously, within
#'   each range (a list with k ranges produces k cutpoints, i.e. k+1 bins).
#'   This is a genuine continuous draw, not snapped to any observed data
#'   value -- snapping to observed values would distort the prior toward
#'   wherever the data happens to be dense.
#' - a plain numeric vector, e.g. `c(12, 16)`: **fixed** (`X_F`). Used as
#'   is, identical across every draw. Also used internally by exact
#'   enumeration to pass a single realized cutpoint vector through
#'   unchanged (see [enumerate_configs()]).
#' - the string `"exact"`: matches on the raw, uncoarsened value of this
#'   covariate (no binning at all).
#' - `NULL`: **excluded**. This covariate does not enter the matching
#'   specification at all.
#' - a [regime()] object: two or more competing elicited cutpoint sets;
#'   see [regime()] for the semantics.
#'
#' @param spec One covariate's cutpoint spec, as described above.
#'
#' @return A list with element `kind` (one of `"excluded"`, `"exact"`,
#'   `"elicited"`, `"fixed"`, `"regime"`), `cutpoints` (the realized numeric
#'   cutpoints, or `NULL` for `"excluded"`/`"exact"`), and, for `kind ==
#'   "regime"`, `regime` giving the name of the regime that was drawn.
#'
#' @keywords internal
draw_cutpoints_for_var <- function(spec) {
  if (is.null(spec)) {
    return(list(kind = "excluded", cutpoints = NULL))
  }
  if (identical(spec, "exact")) {
    return(list(kind = "exact", cutpoints = NULL))
  }
  if (inherits(spec, "cem_regime_spec")) {
    regime_name <- sample(names(spec$regimes), 1, prob = spec$weights)
    ranges <- spec$regimes[[regime_name]]
    cuts <- vapply(ranges, function(r) stats::runif(1, r[1], r[2]), numeric(1))
    return(list(kind = "regime", regime = regime_name, cutpoints = sort(cuts)))
  }
  if (is.list(spec)) {
    cuts <- vapply(spec, function(r) stats::runif(1, r[1], r[2]), numeric(1))
    return(list(kind = "elicited", cutpoints = sort(cuts)))
  }
  if (is.numeric(spec)) {
    ## enumerate_configs() realizes a regime() covariate's drawn cutpoints
    ## into a plain numeric vector for run_draw() to consume (see the
    ## comment in enumerate_var_configs()), tagged with which regime they
    ## came from via a "regime" attribute since the class itself doesn't
    ## survive that flattening. Recover it here so a config produced by
    ## exact enumeration reports kind = "regime" exactly like a Monte Carlo
    ## draw of the same covariate would, rather than silently downgrading
    ## to a plain "fixed" spec with no regime recorded.
    from_regime <- attr(spec, "regime")
    if (is.null(from_regime)) {
      return(list(kind = "fixed", cutpoints = sort(spec)))
    }
    return(list(kind = "regime", regime = from_regime, cutpoints = sort(spec)))
  }
  stop("Unrecognized cutpoint spec of class: ", paste(class(spec), collapse = "/"))
}

#' Coarsen one covariate given its resolved cutpoints
#'
#' Internal building block of [elicit_and_match()]. Applies the cutpoints
#' resolved by [draw_cutpoints_for_var()] to a single covariate vector.
#'
#' @param x A numeric (or, for `kind == "exact"`, any) vector: one
#'   covariate's column from the data.
#' @param resolved The list returned by [draw_cutpoints_for_var()] for this
#'   covariate.
#'
#' @return `NULL` for an excluded covariate; a factor of the raw values for
#'   an exact-match covariate; otherwise an integer vector of bin
#'   memberships from [cut()].
#'
#' @keywords internal
coarsen_one <- function(x, resolved) {
  switch(
    resolved$kind,
    excluded = NULL,
    exact    = factor(x),
    elicited = ,
    fixed    = ,
    regime   = cut(x, breaks = c(-Inf, resolved$cutpoints, Inf), labels = FALSE)
  )
}

## The union range spanned by an elicited or regime() covariate's spec --
## every range across every cutpoint slot (and, for a regime, across every
## regime), collapsed to a single c(lo, hi). Used by pooling_diagnostics()
## to auto-detect which covariates the pre-matching flatness test
## (flatness_test_XE()) should run on, and over what range: the flatness
## question is "is tau(x) flat anywhere this covariate could plausibly be
## coarsened," which spans every regime, not just whichever one a given
## draw happened to use. Returns NULL for excluded/exact/fixed covariates,
## which have no uncertainty to test flatness over.
elicited_union_range <- function(spec) {
  if (inherits(spec, "cem_regime_spec")) {
    ranges <- unlist(spec$regimes, recursive = FALSE)
  } else if (is.list(spec)) {
    ranges <- spec
  } else {
    return(NULL)
  }
  bounds <- unlist(ranges)
  c(min(bounds), max(bounds))
}
