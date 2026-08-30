#' ecem: Elicited Priors for Coarsening Regimes in Coarsened Exact Matching
#'
#' @description
#' Coarsened exact matching (CEM) requires the analyst to fix a coarsening
#' (a set of cutpoints) for every covariate before matching. **ecem**
#' replaces that single fixed choice with an elicited *prior* over
#' cutpoints, draws repeatedly from it, matches and estimates on each draw,
#' and pools the results -- along with diagnostics for when that pooling is
#' (and is not) trustworthy.
#'
#' @section Specifying a covariate's coarsening:
#' Every covariate in a `cutpoint_specs` list is one of five things -- see
#' [regime()] and [draw_cutpoints_for_var()] for the full description:
#' a list of ranges (elicited, `X_E`), a plain numeric vector (fixed,
#' `X_F`), the string `"exact"`, `NULL` (excluded), or a [regime()] object
#' (competing elicited cutpoint sets for one covariate).
#'
#' @section The core workflow:
#' [run_M_draws()] draws (or, when the achievable-configuration count is
#' small enough, exactly enumerates -- see `exact_if_K_leq`) repeated
#' coarsenings, matches, and estimates on each; [pool_draws()] combines
#' them via Rubin's rules. [pooling_diagnostics()] bundles the checks for
#' whether that pooling is valid ([flatness_test_XE()],
#' [congeniality_test()], and the retention-interaction diagnostic) into
#' one call; if pooling is not valid, [existence_test()] provides a
#' fallback (and can reuse [pooling_diagnostics()]'s bootstrap rather than
#' rerunning one). [label_diagnostics()] bundles [cov_att_ate()] and
#' [cov_retention_per_draw()], the covariance-based diagnostics for
#' labeling the pooled estimate as an FSATT, ATT, or ATE.
#'
#' @section The within-draw estimator:
#' [fit_effect()] estimates each draw's effect on its own matched sample
#' via a CEM-weighted regression by default (`estimator = "regression"`,
#' [run_M_draws()]'s and [run_draw()]'s default too) -- the "run whatever
#' model you would have run anyway, on the matched data" approach the CEM
#' literature itself recommends -- or, via `estimator = "mean_diff"`, a
#' simpler treated-count-weighted stratum mean-difference with no
#' regression adjustment. Either way, the covariate-indexed `tau(x)`
#' surface the covariance diagnostics need is always the stratum
#' mean-difference version, regardless of which estimator produces the
#' headline `tau_hat`.
#'
#' @section Exact enumeration:
#' [count_achievable_configs()] reports how many practically distinct
#' coarsenings are actually achievable given a covariate's elicited range
#' and the data's own observed values; when that number is small,
#' [enumerate_configs()] / [run_enumerated_draws()] / [pool_rubins_rules_exact()]
#' replace Monte Carlo draws with exact computation, eliminating simulation
#' error entirely. `run_M_draws(..., exact_if_K_leq = )` does this
#' automatically.
#'
#' @keywords internal
"_PACKAGE"
