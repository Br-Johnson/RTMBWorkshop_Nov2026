## examples/GAM.qmd: penalised thin-plate spline for the year effect.
##  A. As written: lambda is estimated as an ordinary parameter next to the spline
##     coefficients, so the "penalised" fit collapses to the unpenalised one.
##  B. The mixed-model form (mgcv::smooth2random) reproduces mgcv's ML fit.
##  C. The same recipe inside the workshop's hierarchical Ricker model.

source(here::here("review", "verification", "00_setup.R"))
suppressPackageStartupMessages(library(mgcv))
sox <- load_sox()
nstocks <- max(sox$stock); nyears <- max(sox$time); nobs <- nrow(sox)
big <- list(eval.max = 1e4, iter.max = 1e4)
cat("years in the smooth:", nyears, "\n")

sm <- smoothCon(s(time, bs = "tp", k = 10), data = data.frame(time = 1:nyears), absorb.cons = TRUE)[[1]]
B <- sm$X
S <- sm$S[[1]]
K <- ncol(B)
cat("basis columns:", K, "| rank of penalty matrix S:", qr(S)$rank, "(the unpenalised null space has dimension", K - qr(S)$rank, ")\n")

init <- function(extra = list()) c(list(
  alpha0 = mean(sox$logRS), re = rep(0, nstocks), logbeta = rep(log(0.5), nstocks),
  logsigma = rep(log(sd(sox$logRS)), nstocks), logsigma_re = log(0.5)), extra)

section("A. The penalised model exactly as written in examples/GAM.qmd")
## --- verbatim objective (second model in the document) --------------------------
negll <- function(pars){

  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta)
  lambda <- exp(loglambda)

  logRS <- OBS(logRS)

  sigma <- exp(logsigma)
  sigma_re <- exp(logsigma_re)
  sigma_obs <- sigma[stock]

  alpha_stock <- re + alpha0
  ADREPORT(alpha_stock)

  time_effect <- (B %*% b)

  mu <- numeric(nobs)
  for( i in 1:nobs )  {
    mu[i] <- alpha_stock[stock[i]] - beta[stock[i]]*spawners_scaled[i] + time_effect[time[i],1]
  }
  ll <- sum(dnorm(re, 0, sigma_re, log = TRUE))
  ll <- ll + sum(dnorm(logRS, mu, sigma_obs, log = TRUE))

  REPORT(mu)
  REPORT(sigma_obs)
  REPORT(time_effect)

  pen <- lambda/2*t(b)%*%S%*%b

  return(-ll + pen)
}
## ------------------------------------------------------------------------------
obj_w <- MakeADFun(func = negll, parameters = init(list(b = rep(0, K), loglambda = 0)), random = "re", silent = TRUE)
fit_w <- suppressWarnings(nlminb(obj_w$par, obj_w$fn, obj_w$gr, control = big))
cat("loglambda at the reported optimum:", round(fit_w$par["loglambda"], 2),
    "-> lambda =", signif(exp(fit_w$par["loglambda"]), 3), "\n")
cat("convergence:", fit_w$convergence, "|", fit_w$message, "| objective", round(fit_w$objective, 4), "\n")
b_w <- fit_w$par[names(fit_w$par) == "b"]
cat("penalty term lambda/2 * b'Sb at that point:", signif(exp(fit_w$par["loglambda"]) / 2 * drop(t(b_w) %*% S %*% b_w), 3), "\n")

## The same model with no penalty at all:
obj_u <- MakeADFun(func = negll, parameters = init(list(b = rep(0, K), loglambda = -30)), random = "re",
                   map = list(loglambda = factor(NA)), silent = TRUE)
fit_u <- suppressWarnings(nlminb(obj_u$par, obj_u$fn, obj_u$gr, control = big))
te_w <- obj_w$report(obj_w$env$last.par.best)$time_effect
te_u <- obj_u$report(obj_u$env$last.par.best)$time_effect
cat("objective with lambda fixed at exp(-30) (no smoothing):", round(fit_u$objective, 4), "\n")
cat("max |fitted year effect: 'penalised' - unpenalised| =", signif(max(abs(te_w - te_u)), 3), "\n")
cat("=> the penalty is >= 0 and vanishes as lambda -> 0, so minimising over lambda removes it.\n",
    "   lambda is only estimable after the coefficients are integrated out (and the\n",
    "   normalising term of their Gaussian 'prior' is kept).\n", sep = "")

section("B. Mixed-model form checked against mgcv (common residual SD so the models are identical)")
rs <- smooth2random(sm, "", type = 2)
Xf <- rs$Xf                 # unpenalised part (here: the linear trend)
Zs <- rs$rand[[1]]          # penalised part: coefficients ~ iid N(0, sigma_s^2)
cat("smooth2random: fixed columns =", ncol(Xf), "| random columns =", ncol(Zs), "\n")
sox$PopNamef <- factor(sox$PopName)
gam_ml <- gam(logRS ~ 0 + PopNamef + PopNamef:spawners_scaled + s(time, bs = "tp", k = 10),
              data = sox, method = "ML")
Xlin <- model.matrix(~ 0 + PopNamef + PopNamef:spawners_scaled, data = sox)
f_b <- function(pars){
  getAll(pars)
  smooth_t <- drop(Xf %*% bf + Zs %*% bs)
  mu <- drop(Xlin %*% coefs) + smooth_t[sox$time]
  nll <- -sum(dnorm(bs, 0, exp(logsd_s), log = TRUE))
  nll <- nll - sum(dnorm(sox$logRS, mu, exp(logsd), log = TRUE))
  REPORT(smooth_t)
  nll
}
obj_b <- MakeADFun(f_b, list(coefs = rep(0, ncol(Xlin)), bf = rep(0, ncol(Xf)), bs = rep(0, ncol(Zs)),
                             logsd_s = 0, logsd = 0), random = "bs", silent = TRUE)
fit_b <- nlminb(obj_b$par, obj_b$fn, obj_b$gr, control = big)
sm_rtmb <- obj_b$report(obj_b$env$last.par.best)$smooth_t
sm_gam  <- drop(predict(gam_ml, newdata = data.frame(time = 1:nyears, PopNamef = sox$PopNamef[1], spawners_scaled = 0),
                        type = "terms", terms = "s(time)"))
## The two smooths are centred differently (over the 75 years here, over the 1245 observations in
## gam), which is absorbed by the stock intercepts; compare them after removing that constant.
ctr <- function(v) v - mean(v)
cat("max |smooth from RTMB - smooth from mgcv::gam(method = 'ML')| after centring both =",
    signif(max(abs(ctr(sm_rtmb) - ctr(sm_gam))), 3),
    "(range of the smooth:", paste(round(range(ctr(sm_gam)), 2), collapse = " to "), ")\n")
cat("residual SD: RTMB", round(exp(fit_b$par["logsd"]), 4), "| mgcv", round(sqrt(gam_ml$sig2 * df.residual(gam_ml) / nobs), 4),
    "(ML scale)   | effective df of s(time) in mgcv:", round(sum(gam_ml$edf[grep("s\\(time\\)", names(gam_ml$edf))]), 2), "\n")
cat("log-likelihood: RTMB (Laplace, ML)", round(-fit_b$objective, 3), "| mgcv ML score", round(-gam_ml$gcv.ubre, 3), "\n")

section("C. The workshop's model with the smooth done properly")
f_c <- function(pars){
  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta); sigma <- exp(logsigma)
  alpha_stock <- re + alpha0
  smooth_t <- drop(Xf %*% bf + Zs %*% bs)
  mu <- alpha_stock[stock] - beta[stock]*spawners_scaled + smooth_t[time]   # vectorised: no loop
  nll <- -sum(dnorm(re, 0, exp(logsigma_re), log = TRUE))
  nll <- nll - sum(dnorm(bs, 0, exp(logsd_s), log = TRUE))
  nll <- nll - sum(dnorm(logRS, mu, sigma[stock], log = TRUE))
  REPORT(smooth_t)
  nll
}
t_loop <- system.time(MakeADFun(func = negll, parameters = init(list(b = rep(0, K), loglambda = 0)), random = "re", silent = TRUE))["elapsed"]
t_vec  <- system.time(obj_c <- MakeADFun(f_c, init(list(bf = rep(0, ncol(Xf)), bs = rep(0, ncol(Zs)), logsd_s = 0)),
                                         random = c("re", "bs"), silent = TRUE))["elapsed"]
cat("time to build the model object: element-by-element loop (as written)", round(t_loop, 2), "s | vectorised", round(t_vec, 2), "s\n")
fit_c <- nlminb(obj_c$par, obj_c$fn, obj_c$gr, control = big)
sdr_c <- sdreport(obj_c)
cat("convergence", fit_c$convergence, "| pdHess", sdr_c$pdHess, "| marginal NLL", round(fit_c$objective, 3),
    "| SD of penalised spline coefficients", round(exp(fit_c$par["logsd_s"]), 3), "\n")
sm_c <- obj_c$report(obj_c$env$last.par.best)$smooth_t
cat("roughness (sum of squared second differences) of the fitted year effect:\n")
print(round(c(as_written_unpenalised = sum(diff(drop(te_w), differences = 2)^2),
              smooth2random_ML = sum(diff(sm_c, differences = 2)^2)), 4))
cat("for reference, marginal NLL of the AR(1) year-effect model from Module 5 on the same data: 1650.68\n",
    "(see 03_module5_random_effects.R) - a smooth trend cannot capture year-to-year shocks shared by all stocks.\n", sep = "")
