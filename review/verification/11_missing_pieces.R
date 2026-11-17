## Working examples of things the workshop does not cover, run on the workshop's own data.
##  A. Biological reference points with standard errors (Lambert W via ADjoint; S_gen via AD uniroot()).
##  B. Epsilon bias correction for a non-linear function of random effects.
##  C. Ricker vs Larkin (delayed density dependence) by stock.
##  D. Leave-future-out comparison of stationary vs random-walk productivity for all 19 stocks.

source(here::here("review", "verification", "00_setup.R"))
sox <- load_sox()
nstocks <- max(sox$stock); nyears <- max(sox$time)
stock_names <- levels(factor(sox$PopName))
big <- list(eval.max = 1e4, iter.max = 1e4)

section("A. Reference points for one stock (Chilko), with delta-method standard errors")
## Lambert W exactly as in ?ADjoint (and on the workshop's 'Lambert W Example' slide)
W <- function(x) {
  logx <- log(x)
  y <- pmax(logx, 0)
  while (any(abs(logx - log(y) - y) > 1e-9, na.rm = TRUE)) {
      y <- y - (y - exp(logx - y)) / (1 + y)
  }
  y
}
dW <- function(x, y, dy) dy / (exp(y) * (1 + y))
LamW <- ADjoint(W, dW)

ch <- sox %>% filter(PopName == "Chilko")
f_rp <- function(pars) {
  getAll(pars, ch, warn = FALSE)
  beta <- exp(logbeta); sigma <- exp(logsigma)
  nll <- -sum(dnorm(logRS, alpha - beta * spawners, sigma, log = TRUE))
  Smax <- 1 / beta                                   # spawners giving maximum recruits
  Seq  <- alpha / beta                               # unfished equilibrium ("carrying capacity")
  Umsy <- 1 - LamW(exp(1 - alpha))                   # Scheuerell (2016)
  Smsy <- Umsy / beta
  ## S_gen: spawners that return S_msy in one generation. Solve it numerically; uniroot() is
  ## differentiable in RTMB >= 1.9. NOTE: in AD mode `lower`/`upper` must be plain numbers, so the
  ## equation is written for the ratio x = S_gen / S_msy, which always lies in (0, 1).
  x_gen <- uniroot(function(x) log(x) + alpha - beta * Smsy * x,
                   lower = 1e-10, upper = 1, tol = 1e-12)$root
  Sgen <- x_gen * Smsy
  ADREPORT(Smax); ADREPORT(Seq); ADREPORT(Smsy); ADREPORT(Umsy); ADREPORT(Sgen)
  nll
}
o_rp <- MakeADFun(f_rp, list(alpha = 1.5, logbeta = 0, logsigma = 0), silent = TRUE)
fit_rp <- nlminb(o_rp$par, o_rp$fn, o_rp$gr, control = big)
sdr_rp <- sdreport(o_rp)
rp <- summary(sdr_rp, "report")
print(round(cbind(rp, CV = rp[, 2] / rp[, 1]), 4))
a <- fit_rp$par["alpha"]; b <- exp(fit_rp$par["logbeta"])
cat("check: S_gen * exp(alpha - beta * S_gen) =", round(rp["Sgen", 1] * exp(a - b * rp["Sgen", 1]), 5),
    "| S_msy =", round(rp["Smsy", 1], 5), "\n")
cat("Hilborn (1985) approximation S_msy ~ (alpha/beta)(0.5 - 0.07 alpha) =", round(a / b * (0.5 - 0.07 * a), 4),
    "| beta itself =", round(b, 3), "per million spawners (a rate, not a capacity)\n")

section("B. Bias correction of exp(random effect): yearly productivity index in the AR(1) model")
init <- list(alpha0 = mean(sox$logRS), re = rep(0, nstocks), logbeta = rep(log(0.5), nstocks),
             logsigma = rep(log(sd(sox$logRS)), nstocks), logsigma_re = log(0.5),
             re_time = rep(0, nyears), logitphi = 0, logsigma_time = log(0.5))
f_ar <- function(pars){
  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta); sigma <- exp(logsigma)
  sigma_time <- exp(logsigma_time); phi <- 2*plogis(logitphi)-1
  mu <- (re + alpha0)[stock] + re_time[time] - beta[stock]*spawners_scaled
  nll <- -dautoreg(re_time, phi = phi, scale = sqrt(sigma_time^2/(1-phi^2)), log = TRUE)
  nll <- nll - sum(dnorm(re, 0, exp(logsigma_re), log = TRUE))
  nll <- nll - sum(dnorm(logRS, mu, sigma[stock], log = TRUE))
  prod_index <- exp(alpha0 + re_time)        # recruits per spawner at low density, average stock
  ADREPORT(prod_index)
  nll
}
o_ar <- MakeADFun(f_ar, init, random = c("re", "re_time"), silent = TRUE)
fit_ar <- nlminb(o_ar$par, o_ar$fn, o_ar$gr, control = big)
sdr_bc <- sdreport(o_ar, bias.correct = TRUE, bias.correct.control = list(sd = FALSE))
plug <- sdr_bc$value; corr <- sdr_bc$unbiased$value
cat("ratio (bias-corrected / plug-in) of the yearly productivity index: mean", round(mean(corr / plug), 3),
    "| range", paste(round(range(corr / plug), 3), collapse = " to "), "\n")
yrs <- min(sox$yr) + seq_len(nyears) - 1
show <- c(1, 20, 40, 60, nyears)
print(data.frame(brood_year = yrs[show], plug_in = round(plug[show], 2), bias_corrected = round(corr[show], 2)),
      row.names = FALSE)

section("C. Ricker vs Larkin by stock (ordinary least squares on the same rows)")
raw <- read.csv(here::here("data", "sockeye", "SRDATA2026 1.csv")) %>%
  left_join(read.csv(here::here("data", "sockeye", "FC_StockInfo 2.csv")), by = "PopID")
lagS <- function(d, k) d$spawners[match(d$yr - k, d$yr)]
lk <- do.call(rbind, lapply(split(raw, raw$PopName), function(d) {
  d$S1 <- lagS(d, 1); d$S2 <- lagS(d, 2); d$S3 <- lagS(d, 3)
  d$logRS <- log(d$recruits / d$spawners)
  d <- d[complete.cases(d[, c("logRS", "spawners", "S1", "S2", "S3")]), ]
  m0 <- lm(logRS ~ spawners, data = d); m1 <- lm(logRS ~ spawners + S1 + S2 + S3, data = d)
  cyc <- tapply(raw$spawners[raw$PopName == d$PopName[1]], raw$cycle[raw$PopName == d$PopName[1]], mean, na.rm = TRUE)
  data.frame(stock = d$PopName[1], n = nrow(d), dominance = round(max(cyc) / median(cyc), 1),
             dAIC_Larkin_minus_Ricker = round(AIC(m1) - AIC(m0), 1),
             lag_coefs_negative = sum(coef(m1)[c("S1", "S2", "S3")] < 0))
}))
lk <- lk[order(lk$dAIC_Larkin_minus_Ricker), ]
cat("dominance = mean spawners on the largest cycle line / median over the four lines\n")
print(lk, row.names = FALSE)
cat("stocks where the Larkin form lowers AIC by more than 2:", sum(lk$dAIC_Larkin_minus_Ricker < -2), "of", nrow(lk), "\n")

## Lagged spawners and a drifting productivity can explain the same pattern. Fit all four
## combinations by (marginal) maximum likelihood on identical rows and compare AIC.
nll_rw <- function(d, X) {
  f <- function(pars) {
    getAll(pars)
    a <- c(a_init, a_rest)                       # level in the first year is a fixed effect, so the
    nll <- -sum(dnorm(diff(a), 0, exp(logtau), log = TRUE))   # model nests the stationary one (tau -> 0)
    nll - sum(dnorm(d$logRS, a[d$idx] - drop(X %*% b), exp(logsigma), log = TRUE))
  }
  m <- lm(d$logRS ~ X)
  n_state <- max(d$idx)
  o <- MakeADFun(f, list(a_init = unname(coef(m)[1]), a_rest = rep(unname(coef(m)[1]), n_state - 1),
                         b = -unname(coef(m)[-1]), logsigma = log(sd(resid(m))), logtau = log(0.1)),
                 random = "a_rest", silent = TRUE)
  fit <- suppressWarnings(nlminb(o$par, o$fn, o$gr, control = big, lower = c(rep(-Inf, 1 + ncol(X)), -6, -9),
                                 upper = c(rep(Inf, 1 + ncol(X)), 3, 1)))
  c(nll = fit$objective, tau = unname(exp(fit$par["logtau"])))
}
four <- do.call(rbind, lapply(split(raw, raw$PopName), function(d) {
  d$S1 <- lagS(d, 1); d$S2 <- lagS(d, 2); d$S3 <- lagS(d, 3)
  d$logRS <- log(d$recruits / d$spawners)
  d <- d[complete.cases(d[, c("logRS", "spawners", "S1", "S2", "S3")]), ]
  sc <- sd(d$spawners)                            # scale spawners for conditioning
  X0 <- as.matrix(d[, "spawners", drop = FALSE]) / sc
  X1 <- as.matrix(d[, c("spawners", "S1", "S2", "S3")]) / sc
  d$idx <- d$yr - min(d$yr) + 1
  r0 <- nll_rw(d, X0); r1 <- nll_rw(d, X1)
  aic <- c(Ricker = 2 * -as.numeric(logLik(lm(d$logRS ~ X0))) + 2 * 3,
           Larkin = 2 * -as.numeric(logLik(lm(d$logRS ~ X1))) + 2 * 6,
           Ricker_RW = 2 * r0["nll"] + 2 * 4, Larkin_RW = 2 * r1["nll"] + 2 * 7)
  data.frame(stock = d$PopName[1], t(round(aic - min(aic), 1)), best = names(aic)[which.min(aic)],
             tau_Ricker_RW = round(r0["tau"], 2), check.names = FALSE)
}))
names(four)[2:5] <- c("Ricker", "Larkin", "Ricker+RW", "Larkin+RW")
four$best <- sub("\\.nll$", "", four$best)
cat("\nAIC difference from the best of four models, by stock (0 = best):\n")
print(four[order(four$best, four$stock), ], row.names = FALSE)
cat("best model counts:\n"); print(table(sub("_", "+", four$best)))
cat("stocks where the plain Ricker is within 2 AIC units of the best model:", sum(four$Ricker <= 2), "of", nrow(four), "\n")

section("D. Leave-future-out check: stationary vs random-walk productivity, all stocks")
## Forecast target: log(R/S) of brood year t, using only brood years <= t - 2 (when year t is
## forecast pre-season, the age-5 returns of brood t - 1 have not happened yet).
f_rw <- function(d, S_target) function(pars) {
  getAll(pars)
  beta <- exp(logbeta); sigma <- exp(logsigma); tau <- exp(logtau)
  nll <- -dnorm(a[1], a0, 5, log = TRUE) - sum(dnorm(diff(a), 0, tau, log = TRUE))
  nll <- nll - sum(dnorm(d$logRS, a[d$idx] - beta * d$S, sigma, log = TRUE))
  pred <- a[length(a)] - beta * S_target
  ADREPORT(pred)
  REPORT(sigma)
  nll
}
lfo_stock <- function(nm, n_test = 20) {
  d_all <- sox[sox$PopName == nm, c("yr", "logRS", "spawners_scaled")]
  names(d_all)[3] <- "S"
  d_all <- d_all[order(d_all$yr), ]
  targets <- tail(d_all$yr, n_test)
  do.call(rbind, lapply(targets, function(t) {
    d <- d_all[d_all$yr <= t - 2, ]
    if (nrow(d) < 15) return(NULL)
    tgt <- d_all[d_all$yr == t, ]
    ## stationary Ricker (closed form)
    m <- lm(logRS ~ S, data = d)
    ps <- predict(m, newdata = tgt, se.fit = TRUE)
    sd_s <- sqrt(ps$se.fit^2 + sum(resid(m)^2) / nrow(d))
    ## random-walk productivity: states for every year from the first up to t
    d$idx <- d$yr - min(d$yr) + 1
    n_state <- t - min(d$yr) + 1
    o <- MakeADFun(f_rw(d, tgt$S),
                   list(a = rep(coef(m)[1], n_state), a0 = unname(coef(m)[1]),
                        logbeta = log(max(-coef(m)[2], 0.05)), logsigma = log(sd(resid(m))), logtau = log(0.1)),
                   random = "a", map = list(a0 = factor(NA)), silent = TRUE)
    f <- suppressWarnings(nlminb(o$par, o$fn, o$gr, control = big,
                                 lower = c(-Inf, -6, -7), upper = c(Inf, 3, 1)))
    sdr <- suppressWarnings(try(sdreport(o), silent = TRUE))
    okv <- !inherits(sdr, "try-error") && is.finite(sdr$sd[1])
    sig <- exp(f$par["logsigma"])
    pm <- if (okv) unname(sdr$value[1]) else NA
    psd <- if (okv) sqrt(sdr$sd[1]^2 + sig^2) else NA
    data.frame(stock = nm, yr = t, obs = tgt$logRS, stat = unname(ps$fit), stat_sd = unname(sd_s),
               rw = pm, rw_sd = unname(psd), tau = unname(exp(f$par["logtau"])))
  }))
}
t0 <- proc.time()["elapsed"]
lfo <- do.call(rbind, lapply(stock_names, lfo_stock))
cat("forecasts made:", nrow(lfo), "| random-walk fits without a usable sdreport:", sum(is.na(lfo$rw)),
    "| elapsed", round(proc.time()["elapsed"] - t0), "s\n")
lfo <- lfo[!is.na(lfo$rw), ]
sc <- function(pred, sd, obs) c(bias = mean(pred - obs), RMSE = sqrt(mean((pred - obs)^2)),
                                mean_log_score = mean(dnorm(obs, pred, sd, log = TRUE)),
                                cover80 = mean(abs(obs - pred) <= qnorm(0.9) * sd))
cat("\nPooled over stocks (", nrow(lfo), " forecasts of log(R/S), two brood years ahead):\n", sep = "")
print(round(rbind("stationary Ricker"        = sc(lfo$stat, lfo$stat_sd, lfo$obs),
                  "random-walk productivity" = sc(lfo$rw, lfo$rw_sd, lfo$obs)), 3))
cat("(bias > 0 means over-forecasting; mean_log_score: higher is better; cover80: share of observations\n",
    " inside the 80% forecast interval)\n", sep = "")
nh <- lfo$stock != "Harrison"
cat("\nHarrison (ocean-type, ages 3-4, escapements up to 400x apart) dominates the pooled RMSE. Without it:\n")
print(round(rbind("stationary Ricker"        = sc(lfo$stat[nh], lfo$stat_sd[nh], lfo$obs[nh]),
                  "random-walk productivity" = sc(lfo$rw[nh], lfo$rw_sd[nh], lfo$obs[nh])), 3))
cat("median absolute error, all stocks: stationary", round(median(abs(lfo$stat - lfo$obs)), 3),
    "| random walk", round(median(abs(lfo$rw - lfo$obs)), 3), "\n")
cat("share of the", nrow(lfo), "forecasts that were too high: stationary", round(mean(lfo$stat > lfo$obs), 3),
    "| random walk", round(mean(lfo$rw > lfo$obs), 3), "(0.5 if unbiased in the median)\n")
by <- do.call(rbind, lapply(split(lfo, lfo$stock), function(d) data.frame(
  stock = d$stock[1], n = nrow(d),
  bias_stat = mean(d$stat - d$obs), bias_rw = mean(d$rw - d$obs),
  rmse_stat = sqrt(mean((d$stat - d$obs)^2)), rmse_rw = sqrt(mean((d$rw - d$obs)^2)),
  median_tau = median(d$tau))))
by[, -1] <- round(by[, -1], 2)
print(by[order(by$rmse_rw - by$rmse_stat), ], row.names = FALSE)
cat("stocks where the random-walk model has the lower RMSE:", sum(by$rmse_rw < by$rmse_stat), "of", nrow(by),
    "| stocks where the stationary model over-forecasts on average (bias > 0):", sum(by$bias_stat > 0), "of", nrow(by), "\n")
