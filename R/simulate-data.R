#' Simulate a population for examples and demos
#'
#' Generates a simple population with a treatment assignment mechanism
#' depending on `age`, `income`, `educ`, and `region`, and an outcome with
#' either a constant or age-varying treatment effect. Used throughout this
#' package's examples, demos, and tests; not intended as a general-purpose
#' data-generating function for real analyses.
#'
#' @param N Integer; population size.
#' @param tau0 Numeric; baseline treatment effect (the effect at `age =
#'   45` when `heterogeneous = TRUE`, or the constant effect otherwise).
#' @param heterogeneous Logical; if `TRUE`, the treatment effect varies
#'   linearly in `age` (`tau0 + 0.15 * (age - 45)`); if `FALSE`, the effect
#'   is constant at `tau0`.
#'
#' @return A data frame with columns `age`, `income`, `educ`, `region`
#'   (categorical, 3 levels, for exact-match examples), `noise` (unrelated
#'   to treatment or outcome, for excluded-covariate examples), `D` (0/1
#'   treatment), and `Y` (outcome).
#'
#' @examples
#' pop <- simulate_population(N = 500, tau0 = 3, heterogeneous = TRUE)
#' head(pop)
#'
#' @export
simulate_population <- function(N = 3000, tau0 = 3, heterogeneous = FALSE) {
  age    <- stats::runif(N, 20, 75)
  income <- pmax(0, 10000 + 1200 * age + stats::rnorm(N, 0, 10000))
  educ   <- pmin(20, pmax(6, 9 + income / 40000 * 4 + stats::rnorm(N, 0, 2)))
  region <- sample(1:3, N, replace = TRUE)   # categorical, exact-match demo
  noise  <- stats::rnorm(N)                   # excluded-covariate demo

  lin_ps <- -3 + 0.03 * age + 0.00002 * income + 0.05 * educ + 0.10 * (region == 2)
  p <- stats::plogis(lin_ps)
  D <- stats::rbinom(N, 1, p)

  tau_i <- if (heterogeneous) tau0 + 0.15 * (age - 45) else rep(tau0, N)
  Y <- 5 + 0.10 * age + 0.00005 * income + 0.30 * educ + tau_i * D + stats::rnorm(N, 0, 5)

  data.frame(age, income, educ, region, noise, D, Y)
}
