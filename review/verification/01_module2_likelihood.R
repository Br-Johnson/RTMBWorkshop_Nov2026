## Module 2 (modules/2_likelihoods.qmd): two checks.
##  A. The function called `negloglik` returns the log-likelihood (sign flipped twice).
##  B. The "first taste of RTMB" benchmark times a tape that contains one multiplication.

source(here::here("review", "verification", "00_setup.R"))

section("A. `negloglik` on the 'Writing a Likelihood in RTMB' slide")
set.seed(1)
x <- rnorm(50, 2, 5)

## --- verbatim from the slide --------------------------------------------------
negloglik <- \(theta){
  getAll(theta)
  sigma <- exp(logsigma)
  ll <- -sum(dnorm(x, mu, sigma, log = TRUE))
  return(-ll)
}
pars <- list(mu = 0, logsigma = 0)
obj <- MakeADFun(negloglik, pars, silent = TRUE)
## ------------------------------------------------------------------------------

cat("value returned at (mu = 0, logsigma = 0):   ", obj$fn(c(0, 0)), "\n")
cat("true NEGATIVE log-likelihood at that point:", -sum(dnorm(x, 0, 1, log = TRUE)), "\n")
cat("=> the function returns the log-likelihood itself; minimising it maximises the NLL.\n\n")

opt <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr))
cat("nlminb on the slide's function:\n")
print(c(mu = unname(opt$par[1]), sigma = unname(exp(opt$par[2])),
        objective = opt$objective, convergence = opt$convergence))
cat("message:", opt$message, "\n")

nll_ok <- \(theta){
  getAll(theta)
  -sum(dnorm(x, mu, exp(logsigma), log = TRUE))
}
obj2 <- MakeADFun(nll_ok, pars, silent = TRUE)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
cat("\nnlminb on the corrected function (single minus sign):\n")
print(c(mu = unname(opt2$par[1]), sigma = unname(exp(opt2$par[2])),
        objective = opt2$objective, convergence = opt2$convergence))
cat("closed form: mean(x) =", mean(x), "  ML sd =", sqrt(mean((x - mean(x))^2)), "\n")

section("B. What does the taped 'looped' benchmark function contain?")
x <- rnorm(1000, 0, 5)
## --- verbatim from the 'Fitting a Likelihood in R' slide ----------------------
ll_loop <- \(a){
  ll <- 0
  for( i in 1:1000){
    ll <- ll + dnorm(x[i], 0, 5, log = TRUE)
  }
  ll*a
}
ll_loopc <- RTMB::MakeTape(ll_loop, 1)
## ------------------------------------------------------------------------------
cat("Tape of the 1000-iteration loop:\n")
ll_loopc$print()
cat("\nThe loop never touches the tape input `a`, so all 1000 density evaluations are\n",
    "done once in R while taping and stored as a single constant. The benchmark\n",
    "therefore times `constant * a`, not a compiled loop.\n", sep = "")

## A version in which the loop really depends on the input:
ll_loop_dep <- \(mu){
  ll <- 0
  for (i in 1:1000) ll <- ll + dnorm(x[i], mu, 5, log = TRUE)
  ll
}
tp <- RTMB::MakeTape(ll_loop_dep, 0)
n_ops <- nrow(tp$data.frame())
cat("\nOperations on the tape when the loop depends on the parameter:", n_ops, "\n")
t_r    <- system.time(for (k in 1:200) ll_loop_dep(rnorm(1)))["elapsed"] / 200
t_tape <- system.time(for (k in 1:2000) tp(rnorm(1)))["elapsed"] / 2000
cat(sprintf("mean time per call: plain R loop %.1f microsec | taped loop %.1f microsec\n",
            1e6 * t_r, 1e6 * t_tape))
