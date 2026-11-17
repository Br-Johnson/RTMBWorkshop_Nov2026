## Module 5 (modules/5_random_effects.qmd).
##  A. First code slide: MakeADFun() without `random=` maximises the JOINT likelihood.
##  B. The GMRF slide: `ll` is used before it is defined.
##  C. Do the three AR(1) codings agree once B is fixed?
##  D. obj$report() uses the LAST evaluated parameters, not the best ones.
##  E. ML vs REML for the between-stock standard deviation.
##  F. checkConsistency() on the hierarchical model.
##  G. Sparsity: state ("centered") vs innovation ("non-centered") coding of a time series.

source(here::here("review", "verification", "00_setup.R"))
sox <- load_sox()
nstocks <- max(sox$stock)
nyears  <- max(sox$time)
big <- list(eval.max = 1e4, iter.max = 1e4)

## Data-informed starting values, used wherever we are not testing the slides' own inits.
init_hier <- function() list(
  alpha0 = mean(sox$logRS), re = rep(0, nstocks),
  logbeta = rep(log(0.5), nstocks), logsigma = rep(log(sd(sox$logRS)), nstocks),
  logsigma_re = log(0.5))

## --- verbatim from the 'Coding in RTMB' slide ---------------------------------
negll <- function(pars){

  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta)

  sigma <- exp(logsigma)
  sigma_re <- exp(logsigma_re)
  alpha <- re + alpha0

  ADREPORT(alpha)

  mu <- alpha[stock] - beta[stock]*spawners_scaled
  logRS <- OBS(logRS)

  ll <- sum(dnorm(re, 0, sigma_re, log = TRUE))
  ll <- ll + sum(dnorm(logRS, mu, sigma[stock], log = TRUE))

  return(-ll)
}
## ------------------------------------------------------------------------------

section("A. Joint likelihood (no `random=`), as on the first code slide")
for (seed in 1:3) {
  set.seed(seed)
  parameters <- list(alpha0 = runif(1, 2, 5), re = rnorm(nstocks), logbeta = runif(nstocks, 0, 1),
                     logsigma = runif(nstocks, 3, 5), logsigma_re = runif(1, 0.5, 2))
  obj <- MakeADFun(func = negll, parameters = parameters, silent = TRUE)
  fit <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr, control = big))
  cat(sprintf("seed %d: objective = %.4g | sigma_re = %.3g | max|re| = %.3g | conv = %d (%s)\n",
              seed, fit$objective, exp(fit$par["logsigma_re"]),
              max(abs(fit$par[names(fit$par) == "re"])), fit$convergence, fit$message))
}
cat("The joint density is unbounded as sigma_re -> 0 with re -> 0, so this 'fit' is degenerate.\n")

obj_ml <- MakeADFun(func = negll, parameters = init_hier(), random = "re", silent = TRUE)
fit_ml <- nlminb(obj_ml$par, obj_ml$fn, obj_ml$gr, control = big)
cat(sprintf("With random = \"re\" (Laplace): objective = %.4f | sigma_re = %.4f | conv = %d (%s)\n",
            fit_ml$objective, exp(fit_ml$par["logsigma_re"]), fit_ml$convergence, fit_ml$message))

section("B. GMRF slide as written")
init_ar <- function() c(init_hier(), list(re_time = rep(0, nyears), logitphi = 0, logsigma_time = log(0.5)))
## --- verbatim from the 'As a Gauss Markov Random Field Distribution' slide ----
negll_gmrf <- function(pars){

  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta)

  sigma <- exp(logsigma)
  sigma_re <- exp(logsigma_re)
  sigma_time <- exp(logsigma_time)
  phi <- 2*plogis(logitphi)-1

  alpha_stock <- re + alpha0
  ADREPORT(alpha_stock)

  mu <- alpha_stock[stock] + re_time[time] - beta[stock]*spawners_scaled
  logRS <- OBS(logRS)

  Q<-Matrix::spMatrix(nyears,nyears)
  diag(Q)<-c(1,rep(1+phi^2, nyears-2),1)
  Q[cbind(2:nyears,1:(nyears-1))] <- -phi
  Q[cbind(1:(nyears-1), 2:nyears)] <- -phi
  Q <- Q/(sigma_time^2)
  ll <- ll + dgmrf(re_time, Q=Q, log=TRUE)

  ll <- ll + sum(dnorm(re, 0, sigma_re, log = TRUE))
  ll <- ll + sum(dnorm(logRS, mu, sigma[stock], log = TRUE))

  REPORT(mu)
  REPORT(re_time)

  return(-ll)
}
## ------------------------------------------------------------------------------
if (exists("ll")) rm(ll)
res <- try(MakeADFun(func = negll_gmrf, parameters = init_ar(), random = c("re", "re_time"), silent = TRUE),
           silent = TRUE)
cat("MakeADFun(negll_gmrf, ...) in a clean session ->\n  ",
    if (inherits(res, "try-error")) conditionMessage(attr(res, "condition")) else "no error", "\n")

## Same function with the missing initialisation added (one line: `ll <- 0`).
negll_gmrf_fixed <- negll_gmrf
body(negll_gmrf_fixed) <- as.call(append(as.list(body(negll_gmrf)), quote(ll <- 0), after = 1))
res2 <- try(MakeADFun(func = negll_gmrf_fixed, parameters = init_ar(), random = c("re", "re_time"), silent = TRUE),
            silent = TRUE)
cat("After adding `ll <- 0` ->\n  ",
    if (inherits(res2, "try-error")) conditionMessage(attr(res2, "condition")) else "tapes without error", "\n")

section("C. Three codings of the AR(1) year effect")
ar_template <- function(kind) function(pars){
  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta)
  sigma <- exp(logsigma)
  sigma_re <- exp(logsigma_re)
  sigma_time <- exp(logsigma_time)
  phi <- 2*plogis(logitphi)-1
  alpha_stock <- re + alpha0
  mu <- alpha_stock[stock] + re_time[time] - beta[stock]*spawners_scaled
  logRS <- OBS(logRS)
  if (kind == "loop") {
    ll <- dnorm(re_time[1], 0, sqrt(sigma_time^2/(1-phi^2)), log = TRUE)
    for( t in 2:nyears )
      ll <- ll + dnorm(re_time[t], phi*re_time[t-1], sigma_time, log = TRUE)
  } else if (kind == "dautoreg") {
    ll <- dautoreg(re_time, phi = phi, scale = sqrt(sigma_time^2/(1-phi^2)), log = TRUE)
  } else {
    ## Sparse precision built from constant sparse matrices times AD scalars.
    I0 <- Matrix::Diagonal(nyears)
    I1 <- Matrix::Diagonal(nyears, x = c(0, rep(1, nyears - 2), 0))
    B  <- Matrix::bandSparse(nyears, k = c(-1, 1), diagonals = list(rep(1, nyears - 1), rep(1, nyears - 1)))
    Q <- (I0 + phi^2 * I1 - phi * B) / sigma_time^2
    ll <- dgmrf(re_time, Q = Q, log = TRUE)
  }
  ll <- ll + sum(dnorm(re, 0, sigma_re, log = TRUE))
  ll <- ll + sum(dnorm(logRS, mu, sigma[stock], log = TRUE))
  REPORT(mu)
  REPORT(re_time)
  return(-ll)
}
kinds <- c("loop", "dautoreg", "gmrf")
objs <- lapply(kinds, function(k) MakeADFun(ar_template(k), init_ar(), random = c("re", "re_time"), silent = TRUE))
names(objs) <- kinds
set.seed(11)
p_test <- objs$loop$par + rnorm(length(objs$loop$par), 0, 0.1)
vals <- sapply(objs, function(o) o$fn(p_test))
cat("Marginal NLL at one random parameter vector:\n"); print(vals, digits = 12)
tim <- sapply(objs, function(o) system.time(for (i in 1:20) o$fn(p_test + rnorm(length(p_test), 0, 1e-3)))["elapsed"] / 20)
cat("seconds per marginal-likelihood evaluation:\n"); print(round(tim, 4))
fit_ar <- nlminb(objs$dautoreg$par, objs$dautoreg$fn, objs$dautoreg$gr, control = big)
sdr_ar <- sdreport(objs$dautoreg)
cat(sprintf("AR(1) fit: NLL = %.3f | phi = %.3f | sigma_time = %.3f | sigma_re = %.3f | conv = %d | pdHess = %s | max|grad| = %.2g\n",
            fit_ar$objective, 2*plogis(fit_ar$par["logitphi"]) - 1, exp(fit_ar$par["logsigma_time"]),
            exp(fit_ar$par["logsigma_re"]), fit_ar$convergence, sdr_ar$pdHess, max(abs(sdr_ar$gradient.fixed))))
cat(sprintf("Hierarchical model without year effect: NLL = %.3f  => year effect improves the marginal log-likelihood by %.1f for 3 parameters\n",
            fit_ml$objective, fit_ml$objective - fit_ar$objective))

section("D. obj$report() after nlminb(): last.par vs last.par.best")
obj_d <- MakeADFun(ar_template("dautoreg"), init_ar(), random = c("re", "re_time"), silent = TRUE)
fit_d <- nlminb(obj_d$par, obj_d$fn, obj_d$gr, control = big)
lp <- obj_d$env$last.par; lpb <- obj_d$env$last.par.best
cat("max |last.par - last.par.best| immediately after nlminb():", signif(max(abs(lp - lpb)), 3), "\n")
mu_last <- obj_d$report()$mu
mu_best <- obj_d$report(lpb)$mu
cat("max |mu from obj$report() - mu from obj$report(last.par.best)|:", signif(max(abs(mu_last - mu_best)), 3), "\n")
## Any later evaluation moves last.par - e.g. a likelihood profile, a simulation check, or simply:
invisible(obj_d$fn(fit_d$par + 0.5))
cat("after one further obj$fn() call elsewhere: max |mu from obj$report() - mu at the optimum| =",
    signif(max(abs(obj_d$report()$mu - mu_best)), 3), "\n")
cat("=> always pass the parameters explicitly: obj$report(obj$env$last.par.best).\n")

section("E. ML vs REML")
## E1. The textbook case from the pre-course material: sd of a normal sample.
set.seed(1); x <- rnorm(50)
f <- function(p) -sum(dnorm(x, p$mu, exp(p$logsd), log = TRUE))
o_ml   <- MakeADFun(f, list(mu = 0, logsd = 0), silent = TRUE)
o_reml <- MakeADFun(f, list(mu = 0, logsd = 0), random = "mu", silent = TRUE)
s_ml   <- exp(nlminb(o_ml$par,   o_ml$fn,   o_ml$gr)$par["logsd"])
s_reml <- exp(nlminb(o_reml$par, o_reml$fn, o_reml$gr)$par["logsd"])
print(c(ML = unname(s_ml), "ML * sqrt(n/(n-1))" = unname(s_ml) * sqrt(50 / 49),
        "REML (random = 'mu')" = unname(s_reml), "sd(x)" = sd(x)))

## E2. The hierarchical Ricker: integrate the mean parameter alpha0 out as well.
obj_reml <- MakeADFun(func = negll, parameters = init_hier(), random = c("re", "alpha0"), silent = TRUE)
fit_reml <- nlminb(obj_reml$par, obj_reml$fn, obj_reml$gr, control = big)
print(c(sigma_re_ML = unname(exp(fit_ml$par["logsigma_re"])),
        sigma_re_REML = unname(exp(fit_reml$par["logsigma_re"])),
        ratio = unname(exp(fit_reml$par["logsigma_re"] - fit_ml$par["logsigma_re"])),
        "sqrt(19/18)" = sqrt(19 / 18)))

section("F. checkConsistency() on the hierarchical model (Laplace is EXACT for this normal-normal model)")
## F1. Is obj$simulate() drawing from the right model?  (re is used to build mu BEFORE dnorm(re, ...))
set.seed(42)
p_hat <- obj_ml$env$parList(fit_ml$par, obj_ml$env$last.par.best)
sims <- replicate(400, {
  s <- obj_ml$simulate(obj_ml$env$last.par.best)
  mu_s <- p_hat$alpha0 + s$re[sox$stock] - exp(p_hat$logbeta)[sox$stock] * sox$spawners_scaled
  z <- (s$logRS - mu_s) / exp(p_hat$logsigma)[sox$stock]
  c(sd_re = sd(s$re), mean_z = mean(z), sd_z = sd(z))
})
cat("simulated random effects and observations (mean over 400 simulations):\n")
print(round(c(rowMeans(sims), "fitted sigma_re" = exp(p_hat$logsigma_re)), 4))
cat("=> sd(re) matches sigma_re and the standardised observation errors are N(0,1): the simulation is right.\n\n")

## F2. The p-value reported by summary(checkConsistency()) refers n * m' S^-1 m to a chi-square.
##     That statistic is Hotelling's T^2, whose exact null law is ((n-1)p/(n-p)) * F(p, n-p).
hotelling_p <- function(grad) {
  n <- ncol(grad); p <- nrow(grad); m <- rowMeans(grad)
  q <- n * drop(t(m) %*% solve(var(t(grad)), m))
  c(q = q, p_chisq = 1 - pchisq(q, p), p_hotelling = 1 - pf(q * (n - p) / ((n - 1) * p), p, n - p))
}
cat("number of fixed effects tested:", length(fit_ml$par), "\n")
cat("Default call on the slide, checkConsistency(obj) with n = 100 simulations, five seeds:\n")
resF <- t(sapply(1:5, function(seed) {
  set.seed(seed)
  s <- summary(checkConsistency(obj_ml, par = fit_ml$par))
  h <- hotelling_p(s$gradient)
  c(seed = seed, "s$marginal$p.value" = s$marginal$p.value, "Hotelling-corrected p" = unname(h["p_hotelling"]))
}))
print(signif(resF, 3))
set.seed(99)
s_big <- summary(checkConsistency(obj_ml, par = fit_ml$par, n = 4000))
cat("With n = 4000 simulations: s$marginal$p.value =", signif(s_big$marginal$p.value, 3),
    "| Hotelling-corrected p =", signif(hotelling_p(s_big$gradient)["p_hotelling"], 3), "\n")
cat("largest |estimated bias| / standard error over the fixed effects (n = 4000):",
    signif(max(abs(s_big$marginal$bias) / sqrt(diag(sdreport(obj_ml)$cov.fixed))), 3), "\n")

section("G. Hessian sparsity: states vs innovations for a random-walk year effect")
n <- nyears
set.seed(5)
y <- cumsum(rnorm(n, 0, 0.3)) + rnorm(n, 0, 0.5)
f_state <- function(p) {
  u <- p$u
  nll <- -sum(dnorm(diff(u), 0, exp(p$ls), log = TRUE)) - dnorm(u[1], 0, 10, log = TRUE)
  nll - sum(dnorm(y, u, exp(p$lo), log = TRUE))
}
f_innov <- function(p) {
  u <- cumsum(p$z) * exp(p$ls)             # states rebuilt from standard-normal innovations
  nll <- -sum(dnorm(p$z[-1], 0, 1, log = TRUE)) - dnorm(u[1], 0, 10, log = TRUE)
  nll - sum(dnorm(y, u, exp(p$lo), log = TRUE))
}
o_s <- MakeADFun(f_state, list(u = rep(0, n), ls = log(0.3), lo = log(0.5)), random = "u", silent = TRUE)
o_i <- MakeADFun(f_innov, list(z = rep(0, n), ls = log(0.3), lo = log(0.5)), random = "z", silent = TRUE)
invisible(o_s$fn()); invisible(o_i$fn())
nnz <- function(o) { h <- o$env$spHess(random = TRUE); length(h@x) }
cat("non-zeros in the random-effect Hessian (lower triangle), n =", n, "years:\n")
print(c(states = nnz(o_s), innovations = nnz(o_i), dense_would_be = n * (n + 1) / 2))
