## Module 7 (modules/7_misc_coding.qmd) and 7a (modules/7a_custom_distribution.qmd,
## duplicated as examples/custom_distribution.qmd).
##  A. `make_adfun(pars, data)` ignores its `pars` argument.
##  B. The ad_context() example cannot run as written.
##  C. rdirmulti() does not draw from a Dirichlet-multinomial unless all alpha = 1.
##  D. The Dirichlet-multinomial "negloglik" sums densities, not log-densities.
##  E. ddirmulti_osa(): what the unfinished `for` statement actually does.
##  F. %~% : what the function returns when a statement follows the last `%~%`.
##  G. A uniform prior needs no density at all.

source(here::here("review", "verification", "00_setup.R"))
sox <- load_sox()

section("A. 'Passing Data to RTMB (Option 2)'")
param <- list(logsigma = 0, alpha = 0, logbeta = 0)        # the global the slide relies on
## --- verbatim from the slide ---------------------------------------------------
objective_function <- \(pars){
  getAll(pars, data)
  mu <- alpha - exp(logbeta)*spawners
  sigma <- exp(logsigma)

  nll <- -sum(dnorm(logRS, mu, sigma, log = TRUE))
  return(nll)
}
make_adfun <- function(pars, data){
  data_env <- new.env()

  ## The entire function's environment
  # environment(objective_function) <- environment()

  ## Just data:
  data_env <- local({data <- data; environment()})
  environment(objective_function) <- data_env

  obj <- MakeADFun(objective_function, param, silent = TRUE)
  obj
}
## ------------------------------------------------------------------------------
nadina <- sox %>% filter(PopName == "Nadina")
my_start <- list(logsigma = log(2), alpha = 1.5, logbeta = -1)
obj_nadina <- make_adfun(my_start, nadina)
cat("start values passed to make_adfun():", paste(unlist(my_start), collapse = ", "), "\n")
cat("start values the object actually uses:", paste(obj_nadina$par, collapse = ", "),
    " <- the global `param`, not the `pars` argument\n")
rm(param)
res <- try(make_adfun(my_start, nadina), silent = TRUE)
cat("same call once the global `param` is gone:",
    if (inherits(res, "try-error")) conditionMessage(attr(res, "condition")) else "no error", "\n")

section("B. 'ADcontext' slide")
data(cars)
pred_func <- function(pars) {
  getAll(pars, cars, warn = FALSE)
  mu <- beta[1] + beta[2] * speed
  sigma <- exp(logsigma)
  dist %~% dnorm(mu, sd = sigma)
  if (!RTMB:::ad_context()) {
    sim_matrix <- replicate(10000, rnorm(length(mu), mean = mu, sd = sigma))
    pred_intervals <- apply(sim_matrix, 1, quantile, probs = c(0.025, 0.975))
    REPORT(mu)
    REPORT(pred_intervals)
  }
}
res <- try(MakeADFun(pred_func, list(beta = c(0,0), logsigma = log(sigma)), silent = TRUE), silent = TRUE)
cat("MakeADFun(pred_func, list(beta = c(0,0), logsigma = log(sigma))) ->\n  ",
    if (inherits(res, "try-error")) conditionMessage(attr(res, "condition")) else "no error",
    "\n  (`sigma` is not defined at that point, so R finds the function stats::sigma)\n")
cat("is ad_context exported by RTMB?", "ad_context" %in% getNamespaceExports("RTMB"), "\n")
## With a valid start value the example runs; compare its plug-in interval with lm():
obj_p <- MakeADFun(pred_func, list(beta = c(0,0), logsigma = log(15)), silent = TRUE)
fit_p <- nlminb(obj_p$par, obj_p$fn, obj_p$gr)
set.seed(1)
pi_sim <- t(obj_p$report(fit_p$par)$pred_intervals)
pi_lm  <- predict(lm(dist ~ speed, cars), interval = "prediction")
cat("mean width of the 95% prediction interval: simulated at the MLE (slide) =",
    round(mean(pi_sim[, 2] - pi_sim[, 1]), 1), "| predict.lm =", round(mean(pi_lm[, "upr"] - pi_lm[, "lwr"]), 1),
    "\n  (the plug-in version ignores parameter uncertainty and uses the ML sigma, so it is narrower)\n")

section("C. rdirmulti(): Dirichlet draws")
## --- verbatim from 7a ----------------------------------------------------------
rdirmulti <- function(n, size, alpha) {
  nalpha <- length(alpha)
  counts <- matrix(0, nrow = n, ncol = nalpha)
  for( i in 1:n ){
    theta <- rgamma(length(alpha), 1, alpha)
    theta <- theta/sum(theta)
    counts[i,] <- rmultinom(1, size = size, prob = theta)
  }
  if(n==1) return(drop(counts)) # return as vector
  counts
}
## ------------------------------------------------------------------------------
rdirmulti_ok <- function(n, size, alpha) {
  t(replicate(n, { theta <- rgamma(length(alpha), shape = alpha, rate = 1)
                   drop(rmultinom(1, size = size, prob = theta / sum(theta))) }))
}
alpha <- c(1, 5, 20)
set.seed(1)
m_slide <- colMeans(rdirmulti(20000, 100, alpha)) / 100
m_ok    <- colMeans(rdirmulti_ok(20000, 100, alpha)) / 100
print(round(rbind("theory: alpha / sum(alpha)"           = alpha / sum(alpha),
                  "slide: rgamma(K, 1, alpha)"           = m_slide,
                  "correct: rgamma(K, shape = alpha, 1)" = m_ok), 3))
cat("With alpha = c(1,1,1) the two coincide, which is why the slide's own check passes.\n")

section("D. Dirichlet-multinomial likelihood")
## --- verbatim from 7a (density function and objective) --------------------------
ddirmulti <- function(x, size, alpha, log = FALSE) {
  if(inherits(x, "simref")) {
    x[] <- rdirmulti(1, size=size, alpha=alpha)
    return(0)
  }
  if(inherits(x, "osa")) {
    return(ddirmulti_osa(x, alpha, log))
  }
  alpha0 <- sum(alpha)
  logdens <- lgamma(alpha0) + lgamma(size + 1) - lgamma(size + alpha0)
  logdens <- logdens + sum(lgamma(x + alpha) - lgamma(alpha) - lgamma(x + 1))
  if(log) return(logdens)
  return(exp(logdens))
}
set.seed(2)
x <- rdirmulti_ok(200, 5, c(2, 4, 8))
negloglik <- function(pars){
  getAll(pars)
  alpha <- exp(logalpha)
  x <- OBS(x)

  ll <- 0
  for( i in 1:nrow(x) ){
    ll <- ll + ddirmulti(x[i,], size = 5, alpha = alpha)
  }
  -ll
}
## ------------------------------------------------------------------------------
negloglik_log <- function(pars){
  getAll(pars)
  alpha <- exp(logalpha)
  ll <- 0
  for (i in 1:nrow(x)) ll <- ll + ddirmulti(x[i,], size = 5, alpha = alpha, log = TRUE)
  -ll
}
o1 <- MakeADFun(negloglik,     list(logalpha = c(0, 0, 0)), silent = TRUE)
o2 <- MakeADFun(negloglik_log, list(logalpha = c(0, 0, 0)), silent = TRUE)
f1 <- suppressWarnings(nlminb(o1$par, o1$fn, o1$gr, control = list(eval.max = 1e4, iter.max = 1e4)))
f2 <- nlminb(o2$par, o2$fn, o2$gr)
cat("objective as written at logalpha = 0:", o1$fn(), " = -(sum of densities);  with log = TRUE:", o2$fn(), "\n")
print(round(rbind("true alpha"                     = c(2, 4, 8),
                  "as written (log = FALSE)"       = exp(f1$par),
                  "with log = TRUE"                = exp(f2$par)), 2))
cat("as-written fit: convergence =", f1$convergence, "|", f1$message, "\n")

section("E. ddirmulti_osa() as written")
f_txt <- 'for( a in 1:length(alpha-1) )

  # Calculate conditional density metrics
  alpha0_sub <- sum(alpha_sub)
  size_sub <- 1'
e <- parse(text = f_txt)
cat("The `for (...)` line has no body, so R takes the NEXT statement as the loop body. Parsed as:\n")
print(e[[1]])
alpha <- c(2, 4, 8)
cat("and 1:length(alpha-1) is", paste(1:length(alpha - 1), collapse = " "),
    "(length of alpha, not 1:(length(alpha)-1) =", paste(1:(length(alpha) - 1), collapse = " "), ")\n")

section("F. What a %~% model returns")
set.seed(3); y <- rnorm(30, 5, 2)
f_tilde <- function(pars) { getAll(pars); y %~% dnorm(mu, exp(logsigma)); ADREPORT(exp(logsigma)) }
o <- MakeADFun(f_tilde, list(mu = 0, logsigma = 0), silent = TRUE)
cat("objective:", o$fn(), "| direct:", -sum(dnorm(y, 0, 1, log = TRUE)),
    " -> the hidden .nll is returned even though ADREPORT() is the last statement\n")

## Mixing the two styles: a hand-accumulated term plus %~%, returned explicitly.
f_mixed <- function(pars) {
  getAll(pars)
  nll <- -dnorm(mu, 0, 0.1, log = TRUE)          # e.g. a prior/penalty added by hand
  y %~% dnorm(mu, exp(logsigma))
  return(nll)                                    # looks as if the hand-built total is returned
}
o_mixed <- MakeADFun(f_mixed, list(mu = 0.5, logsigma = 0), silent = TRUE)
by_hand <- -dnorm(0.5, 0, 0.1, log = TRUE)
by_tilde <- -sum(dnorm(y, 0.5, 1, log = TRUE))
cat("mixed function returns:", o_mixed$fn(), "| hand term only:", by_hand, "| %~% term only:", by_tilde,
    "| both:", by_hand + by_tilde, "\n")
cat("=> ", if (isTRUE(all.equal(o_mixed$fn(), by_tilde))) "the explicit return() is overridden: the hand-accumulated term is silently dropped."
           else if (isTRUE(all.equal(o_mixed$fn(), by_hand))) "only the hand-accumulated term is used: the %~% term is silently dropped."
           else if (isTRUE(all.equal(o_mixed$fn(), by_hand + by_tilde))) "both terms are included."
           else "unexpected value.", "\n", sep = "")

section("G. Uniform prior = bounds, not a density")
## A Uniform(0, 1) prior on a probability p contributes a constant, so nothing is added to the NLL.
## Enforce the support by transformation (and, for MCMC, add the log-Jacobian of that transformation).
k <- 7; n <- 20
f_unif <- function(pars) { p <- plogis(pars$logitp); -dbinom(k, n, p, log = TRUE) }
o <- MakeADFun(f_unif, list(logitp = 0), silent = TRUE)
cat("MLE of p with implicit Uniform(0,1) prior via logit transform:",
    round(plogis(nlminb(o$par, o$fn, o$gr)$par), 4), "| k/n =", k / n, "\n")
