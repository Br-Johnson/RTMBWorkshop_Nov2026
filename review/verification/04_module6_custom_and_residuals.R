## Module 6 (modules/6_custom_functions.qmd, modules/6_residuals.qmd).
##  A. Global on/off flags: obj$simulate() follows the flag's CURRENT value, not the taped one.
##  B. The slide's `dlgamma` masks RTMB's exported `dlgamma` (different parameterisation).
##  C. RTMB already ships dlnorm() with simulation and OSA support.
##  D. Additive-log-ratio Jacobian: closed form vs tape + log-determinant.
##  E. Poisson-gamma: what leaving out the Jacobian does to the estimates.
##  F. OSA residual methods for the Gaussian sockeye model: agreement and cost.

source(here::here("review", "verification", "00_setup.R"))
sox <- load_sox()
nstocks <- max(sox$stock); nyears <- max(sox$time)
big <- list(eval.max = 1e4, iter.max = 1e4)

section("A. 'Log Normal Distribution Cont...' + 'Checking Simulation' slides")
soxch <- sox %>% filter(PopName == "Chilko")
calls <- c(dlognorm_branch = 0, dnorm_branch = 0)
## --- verbatim from the slides (only the two counters are added) ----------------
dlognorm <- function(x, mean, sd, log = FALSE){
  y <- log(x)
  logpdens <- dnorm(y, mean, sd, log = TRUE) - y # -y == log(1/x)

  if (inherits(x, "simref")) {
    return(RTMB:::dGenericSim(
      "dlognorm",
      x = x,
      mean = mean,
      sd = sd,
      log = log
    ))
  }

  if(log) return(logpdens)
  return(exp(logpdens))
}
rlognorm <- function(n, mean, sd){
  exp(rnorm(n, mean, sd))
}
param <- list(logsigma = 0, alpha = 0, logbeta = 0)
soxch$RS <- soxch$recruits/soxch$spawners
negll <- function(pars){
  getAll(pars, soxch)
  mu <- alpha - exp(logbeta)*spawners
  sigma <- exp(logsigma)
  RS <- OBS(RS)

  if(lognormal == TRUE){
    calls["dlognorm_branch"] <<- calls["dlognorm_branch"] + 1     # counter (added)
    nll <- -sum(dlognorm(RS, mu, sigma, log = TRUE))
  }else{
    calls["dnorm_branch"] <<- calls["dnorm_branch"] + 1           # counter (added)
    nll <- -sum(dnorm(log(RS), mu, sigma, log = TRUE))
  }
  return(nll)
}
lognormal <- TRUE
obj <- MakeADFun(negll, param, silent = TRUE)
lognormal <- FALSE
obj2 <- MakeADFun(negll, param, silent = TRUE)
## ------------------------------------------------------------------------------
cat("objective at the start values:  obj (taped with lognormal = TRUE) =", obj$fn(),
    "| obj2 (taped with lognormal = FALSE) =", obj2$fn(), "\n")
cat("difference = sum(log(RS)) =", sum(log(soxch$RS)), " (the Jacobian term, constant in the parameters)\n")
calls[] <- 0
invisible(obj$simulate())          # as on the 'Checking Simulation' slide: sim1 uses `obj`
cat("Branch executed by obj$simulate() while the global flag is FALSE:\n"); print(calls)
cat("=> `sim1` and `sim2` on the slide are both drawn through the dnorm(log(RS)) branch;\n",
    "   the custom dlognorm()/rlognorm() simulation code is never exercised by that check.\n", sep = "")
## Now make the simulation actually go through dlognorm(), as the slide intended:
lognormal <- TRUE; calls[] <- 0
err <- try(obj$simulate(), silent = TRUE)
cat("\nWith the flag set back to TRUE, obj$simulate() gives:\n  ",
    if (inherits(err, "try-error")) paste("ERROR:", conditionMessage(attr(err, "condition"))) else "no error", "\n")
cat("Cause: dlognorm() calls dnorm(log(x), ...) BEFORE testing inherits(x, 'simref'). RTMB's dnorm()\n",
    "already simulates x through the log(); the explicit dGenericSim() call then tries to assign\n",
    "the same observation a second time.\n", sep = "")

## Corrected custom density: branch on 'simref' first, do nothing else in simulation mode.
dlognorm2 <- function(x, mean, sd, log = FALSE){
  if (inherits(x, "simref"))
    return(RTMB:::dGenericSim("dlognorm2", x = x, mean = mean, sd = sd, log = log))
  y <- log(x)
  logpdens <- dnorm(y, mean, sd, log = TRUE) - y
  if (log) logpdens else exp(logpdens)
}
rlognorm2 <- function(n, mean, sd) exp(rnorm(n, mean, sd))
negll2 <- function(pars){
  getAll(pars, soxch)
  mu <- alpha - exp(logbeta)*spawners
  RS <- OBS(RS)
  -sum(dlognorm2(RS, mu, exp(logsigma), log = TRUE))
}
obj_fix <- MakeADFun(negll2, param, silent = TRUE)
lognormal <- FALSE
set.seed(1); s_fix <- obj_fix$simulate()$RS
set.seed(1); s_ref <- obj2$simulate()$RS
cat("Corrected dlognorm2(): simulates without error; same draws as the dnorm(log(RS)) model with the\n",
    "same seed: max |difference| =", signif(max(abs(s_fix - s_ref)), 3), "\n")

## The same thing without any unexported RTMB function: in simulation mode, assign the draw to
## the observation and return 0 (the pattern RTMB's maintainer suggests for custom densities).
dlognorm3 <- function(x, mean, sd, log = FALSE){
  if (inherits(x, "simref")) {
    x[] <- exp(rnorm(length(x), mean, sd))
    return(0)
  }
  y <- log(x)
  logpdens <- dnorm(y, mean, sd, log = TRUE) - y
  if (log) logpdens else exp(logpdens)
}
negll3 <- function(pars){
  getAll(pars, soxch)
  mu <- alpha - exp(logbeta)*spawners
  RS <- OBS(RS)
  -sum(dlognorm3(RS, mu, exp(logsigma), log = TRUE))
}
obj_fix3 <- MakeADFun(negll3, param, silent = TRUE)
set.seed(1); s_fix3 <- obj_fix3$simulate()$RS
cat("dlognorm3() (no RTMB:::): objective", obj_fix3$fn(), "| same draws as above: max |difference| =",
    signif(max(abs(s_fix3 - s_ref)), 3), "\n")

section("B. Name clash: dlgamma")
cat("RTMB exports dlgamma with arguments:", paste(names(formals(RTMB::dlgamma)), collapse = ", "), "\n")
dlgamma_slide <- function(x, shape, rate, log = FALSE){     # as defined on the 'When it Matters' slide
  y <- exp(x)
  logpdens <- dgamma(y, shape = shape, rate = rate, log = TRUE) + x
  if(log) return(logpdens)
  return(exp(logpdens))
}
x0 <- log(3.7)
cat("slide's dlgamma(x, shape = 5, rate = 2)        =", dlgamma_slide(x0, 5, 2, log = TRUE), "\n")
cat("RTMB::dlgamma(x, shape = 5, scale = 1/2)        =", RTMB::dlgamma(x0, 5, 1 / 2, log = TRUE), "\n")
cat("RTMB::dlgamma(x, 5, 2) (positional, = scale 2)  =", RTMB::dlgamma(x0, 5, 2, log = TRUE),
    " <- what a reader gets if the slide's definition is not in scope\n")

section("C. Built-in dlnorm() with OBS(): likelihood, simulation and OSA")
negll_ln <- function(pars){
  getAll(pars, soxch)
  mu <- alpha - exp(logbeta)*spawners
  RS <- OBS(RS)
  -sum(dlnorm(RS, mu, exp(logsigma), log = TRUE))
}
obj_ln <- MakeADFun(negll_ln, param, silent = TRUE)
cat("NLL with built-in dlnorm at the start values:", obj_ln$fn(), "(custom dlognorm:", obj$fn(), ")\n")
set.seed(2); sim_ln <- obj_ln$simulate()$RS
cat("obj$simulate() works with dlnorm: all simulated RS positive =", all(sim_ln > 0),
    "| n =", length(sim_ln), "\n")
fit_ln <- nlminb(obj_ln$par, obj_ln$fn, obj_ln$gr)
try_osa <- function(o, ...) {
  r <- try(suppressWarnings(oneStepPredict(o, trace = FALSE, ...)), silent = TRUE)
  if (inherits(r, "try-error")) paste("ERROR:", conditionMessage(attr(r, "condition")))
  else sprintf("ok (residual mean %.3f, sd %.3f)", mean(r$residual), sd(r$residual))
}
cat("OSA with dlnorm, method = 'oneStepGeneric', range = c(0, Inf):",
    try_osa(obj_ln, method = "oneStepGeneric", range = c(0, Inf)), "\n")
cat("OSA with dlnorm, method = 'cdf':",
    try_osa(obj_ln, method = "cdf"), "\n")
cat("  (RTMB 2.0 has no AD version of plnorm(), which the 'cdf' method needs.)\n")
## The same model with the observation defined on the log scale supports every OSA method:
soxch$logRS_ch <- log(soxch$RS)
negll_lg <- function(pars){
  getAll(pars, soxch)
  mu <- alpha - exp(logbeta)*spawners
  logRS_ch <- OBS(logRS_ch)
  -sum(dnorm(logRS_ch, mu, exp(logsigma), log = TRUE))
}
obj_lg <- MakeADFun(negll_lg, param, silent = TRUE)
fit_lg <- nlminb(obj_lg$par, obj_lg$fn, obj_lg$gr)
cat("OSA with dnorm on OBS(log RS), method = 'cdf':         ", try_osa(obj_lg, method = "cdf"), "\n")
cat("OSA with dnorm on OBS(log RS), method = 'fullGaussian':", try_osa(obj_lg, method = "fullGaussian"), "\n")
cat("same parameter estimates either way: max |diff| =", signif(max(abs(fit_ln$par - fit_lg$par)), 3),
    "| NLL(dlnorm) - NLL(dnorm on logs) =", signif(fit_ln$objective - fit_lg$objective, 6),
    "= sum(log RS) =", signif(sum(log(soxch$RS)), 6), "\n")

section("D. Jacobian of the additive log-ratio (multivariate logit) transform")
mlogit <- \(x){ n <- length(x); log(x[-n]/x[n]) }
p <- c(0.1, 0.5, 0.2, 0.05, 0.15)
## alr acts on the first n-1 free proportions, with p_n = 1 - sum(p_1..p_{n-1})
alr_free <- function(pf) { pn <- 1 - sum(pf); log(pf / pn) }
J <- MakeTape(alr_free, p[-5])$jacobian(p[-5])
cat("log|det J| from a tape + determinant():", as.numeric(determinant(J, logarithm = TRUE)$modulus), "\n")
cat("closed form  -sum(log(p)):              ", -sum(log(p)), "\n")
## determinant() also works directly on AD matrices - no ADjoint needed:
Fd <- MakeTape(function(x) determinant(matrix(x, 2, 2), logarithm = TRUE)$modulus, c(2, 0.3, 0.3, 1))
cat("determinant() taped natively on an AD matrix: value =", Fd(c(2, 0.3, 0.3, 1)),
    "| base R =", log(det(matrix(c(2, 0.3, 0.3, 1), 2))), "\n")

section("E. Poisson counts with a gamma-distributed rate (true shape = 5, rate = 1)")
set.seed(10)
lambda <- rgamma(100, 5, 1)
count <- rpois(100, lambda)
fit_pg <- function(kind) {
  f <- function(pars){
    getAll(pars)
    a <- exp(logshape); b <- exp(lograte)
    if (kind == "nb") return(-sum(dnbinom(count, a, b/(1+b), log = TRUE)))
    mu <- exp(logre)
    nll <- if (kind == "with_jacobian") -sum(dgamma(mu, a, b, log = TRUE) + logre)   # density of log(lambda)
           else                         -sum(dgamma(mu, a, b, log = TRUE))           # Jacobian omitted
    nll - sum(dpois(count, mu, log = TRUE))
  }
  start <- list(logshape = log(3), lograte = 0)
  if (kind == "nb") {
    o <- MakeADFun(f, start, silent = TRUE)
  } else {
    o <- MakeADFun(f, c(start, list(logre = log(count + 0.5))), random = "logre", silent = TRUE)
  }
  r <- try(suppressWarnings(nlminb(o$par, o$fn, o$gr, control = big)), silent = TRUE)
  if (inherits(r, "try-error")) return(c(shape = NA, rate = NA, mean = NA, nll = NA, convergence = NA))
  c(shape = unname(exp(r$par["logshape"])), rate = unname(exp(r$par["lograte"])),
    mean = unname(exp(r$par["logshape"] - r$par["lograte"])),
    nll = r$objective, convergence = r$convergence)
}
print(round(rbind("exact negative binomial"            = fit_pg("nb"),
                  "Laplace, Jacobian included"         = fit_pg("with_jacobian"),
                  "Laplace, Jacobian omitted (wrong)"  = fit_pg("without_jacobian")), 3))
cat("(NA = the optimiser stopped with a NaN gradient. Without the Jacobian the implied density of\n",
    " log(lambda) is not a gamma on lambda; it is improper as log(lambda) -> -Inf when shape <= 1.)\n", sep = "")

section("F. OSA residuals for the Gaussian AR(1) sockeye model")
init <- list(alpha0 = mean(sox$logRS), re = rep(0, nstocks), logbeta = rep(log(0.5), nstocks),
             logsigma = rep(log(sd(sox$logRS)), nstocks), logsigma_re = log(0.5),
             re_time = rep(0, nyears), logitphi = 0, logsigma_time = log(0.5))
make_obj <- function(dat) {
  f <- function(pars){
    getAll(pars, dat, warn = FALSE)
    beta <- exp(logbeta); sigma <- exp(logsigma)
    sigma_re <- exp(logsigma_re); sigma_time <- exp(logsigma_time)
    phi <- 2*plogis(logitphi)-1
    logRS <- OBS(logRS)
    alpha_stock <- re + alpha0
    mu <- alpha_stock[stock] + re_time[time] - beta[stock]*spawners_scaled
    ll <- dautoreg(re_time, phi = phi, scale = sqrt(sigma_time^2/(1-phi^2)), log = TRUE)
    ll <- ll + sum(dnorm(re, 0, sigma_re, log = TRUE))
    ll <- ll + sum(dnorm(logRS, mu, sigma[stock], log = TRUE))
    REPORT(mu)
    -ll
  }
  o <- MakeADFun(f, init, random = c("re", "re_time"), silent = TRUE)
  fit <- nlminb(o$par, o$fn, o$gr, control = big)
  list(obj = o, fit = fit)
}
m <- make_obj(sox)
cat("data order in the slides: sorted by stock, then year (first rows):\n")
print(head(sox[, c("PopName", "yr")], 3))
tm <- function(expr) { t0 <- proc.time()["elapsed"]; v <- suppressWarnings(expr); list(v = v, t = proc.time()["elapsed"] - t0) }
r_full <- tm(oneStepPredict(m$obj, method = "fullGaussian", trace = FALSE))
r_osg  <- tm(oneStepPredict(m$obj, method = "oneStepGaussian", trace = FALSE))
r_cdf  <- tm(oneStepPredict(m$obj, method = "cdf", trace = FALSE))
r_gen  <- tm(oneStepPredict(m$obj, method = "oneStepGeneric", trace = FALSE))   # the slide's choice
cmp <- data.frame(
  method  = c("fullGaussian", "oneStepGaussian", "cdf", "oneStepGeneric (as on the slide)"),
  seconds = round(c(r_full$t, r_osg$t, r_cdf$t, r_gen$t), 1),
  n       = c(nrow(r_full$v), nrow(r_osg$v), nrow(r_cdf$v), nrow(r_gen$v)),
  max_abs_diff_vs_fullGaussian = signif(c(0,
        max(abs(r_osg$v$residual - r_full$v$residual)),
        max(abs(r_cdf$v$residual - r_full$v$residual)),
        max(abs(r_gen$v$residual - r_full$v$residual))), 3))
print(cmp, row.names = FALSE)
dg <- abs(r_gen$v$residual - r_full$v$residual)
cat("oneStepGeneric vs the exact residuals: |difference| > 0.01 for", sum(dg > 0.01), "of", length(dg),
    "observations, > 0.1 for", sum(dg > 0.1), "; non-finite:", sum(!is.finite(r_gen$v$residual)), "\n")
## ?oneStepPredict: "Pay special attention to the integration domain which must be set ... via `range`"
r_gen2 <- tm(oneStepPredict(m$obj, method = "oneStepGeneric", trace = FALSE,
                            range = range(sox$logRS) + c(-10, 10)))
dg2 <- abs(r_gen2$v$residual - r_full$v$residual)
cat(sprintf("oneStepGeneric with an explicit finite `range`: %.1f s | max |difference| %.3g | > 0.01 for %d observations\n",
            r_gen2$t, max(dg2), sum(dg2 > 0.01)))
res_full <- r_full$v$residual
cat(sprintf("fullGaussian residuals: mean %.3f, sd %.3f, Shapiro-Wilk p = %.3g\n",
            mean(res_full), sd(res_full), shapiro.test(res_full)$p.value))
cat("excess kurtosis of the OSA residuals:", round(mean((res_full - mean(res_full))^4) / var(res_full)^2 - 3, 2),
    "| share with |z| > 3:", round(mean(abs(res_full) > 3), 4), "(normal: 0.0027)\n")

## Ordering: OSA residuals condition on the observations that come EARLIER in the data vector.
sox_t <- sox[order(sox$yr, sox$stock), ]
m_t <- make_obj(sox_t)
res_t <- suppressWarnings(oneStepPredict(m_t$obj, method = "fullGaussian", trace = FALSE))$residual
cat("same model, data sorted by year then stock: NLL difference =", signif(m_t$fit$objective - m$fit$objective, 3), "\n")
key <- paste(sox$PopName, sox$yr); key_t <- paste(sox_t$PopName, sox_t$yr)
cat("correlation between the two sets of OSA residuals (matched by stock-year):",
    round(cor(res_full, res_t[match(key, key_t)]), 3), "\n")
cat("=> both orderings give valid N(0,1) residuals, but they are different residuals; only the\n",
    "   time-ordered version has the 'predict this year from earlier years' interpretation.\n", sep = "")
## Serial dependence left in the time-ordered OSA residuals, within stock:
rt <- data.frame(stock = sox_t$stock, yr = sox_t$yr, r = res_t)
lag1 <- merge(rt, transform(rt, yr = yr + 1), by = c("stock", "yr"), suffixes = c("", "_prev"))
ct <- cor.test(lag1$r, lag1$r_prev)
cat(sprintf("within-stock lag-1 correlation of the time-ordered OSA residuals: %.3f (95%% CI %.3f to %.3f, n = %d pairs)\n",
            ct$estimate, ct$conf.int[1], ct$conf.int[2], nrow(lag1)))
by_stock <- sapply(split(lag1, lag1$stock), function(d) cor(d$r, d$r_prev))
names(by_stock) <- levels(factor(sox$PopName))
cat("per-stock lag-1 correlations (largest five):\n"); print(round(sort(by_stock, decreasing = TRUE)[1:5], 2))
cat("=> the shared year effect leaves stock-specific serial correlation that a QQ-plot cannot show.\n")
cat("conditional (plug-in) standardized residuals, as on the 'Classic Residual Analysis' slide: sd =",
    round(sd((sox$logRS - m$obj$report(m$obj$env$last.par.best)$mu) /
               exp(m$obj$env$parList(m$fit$par, m$obj$env$last.par.best)$logsigma)[sox$stock]), 3),
    " (below 1 because the random effects are fitted to the same data)\n")
