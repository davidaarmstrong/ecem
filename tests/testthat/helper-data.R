## Deterministic toy dataset for exact, hand-verified expectations. A
## cutpoint at x = 7 splits it into two perfectly balanced strata, each
## with 2 treated and 2 control units and an identical within-stratum mean
## difference of 5.5, so the pooled ATT is exactly 5.5 and its variance is
## exactly 0.625 regardless of how the two strata are weighted (they have
## equal n1, so the weights are equal too). Worked by hand:
##   stratum 1 (x in {1,2,3,4}):   D=1 mean(10,12)=11, D=0 mean(5,6)=5.5,  diff=5.5
##   stratum 2 (x in {11..14}):    D=1 mean(20,22)=21, D=0 mean(15,16)=15.5, diff=5.5
##   var: v1 = var(c(10,12))/2 = 1,    v0 = var(c(5,6))/2  = 0.25 -> stratum var 1.25
##        v1 = var(c(20,22))/2 = 1,    v0 = var(c(15,16))/2 = 0.25 -> stratum var 1.25
##   var_hat = 0.5^2 * 1.25 + 0.5^2 * 1.25 = 0.625
make_toy_data <- function() {
  data.frame(
    x = c(1, 2, 3, 4, 11, 12, 13, 14),
    D = c(1, 1, 0, 0, 1, 1, 0, 0),
    Y = c(10, 12, 5, 6, 20, 22, 15, 16)
  )
}

## Adds a third stratum (x > 20) containing only treated units, which
## elicit_and_match() should prune for lacking control support.
make_toy_data_with_unmatched <- function() {
  base <- make_toy_data()
  extra <- data.frame(x = c(25, 26), D = c(1, 1), Y = c(30, 31))
  rbind(base, extra)
}
