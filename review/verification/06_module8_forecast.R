## Module 8 (modules/8_practical.qmd): the Stellako age-structured forecast.
##  A. The covariate CSV starts with a byte-order mark: the join on "yr" is locale dependent.
##  B. Data check: `recruits` vs the sum of the age columns.
##  C. Covariate model as written: intercept counted twice -> singular Hessian.
##  D. Sensitivity to the "+ 1e-4" constant.
##  E. Median vs mean forecast.
##  F. Retrospective skill: the slides' models vs naive benchmarks, plus interval coverage.

source(here::here("review", "verification", "00_setup.R"))
suppressPackageStartupMessages(library(MASS))
select <- dplyr::select
big <- list(eval.max = 1e4, iter.max = 1e4)
cov_file <- here::here("data", "sockeye", "FC_Environmental_Data_2026_Std_Offset_Real.csv")

section("A. Byte-order mark in the covariate file")
cat("first three bytes of the file:", paste(readBin(cov_file, "raw", 3), collapse = " "), "(ef bb bf = UTF-8 BOM)\n")
old <- Sys.getlocale("LC_CTYPE")
for (loc in c("C", "C.UTF-8")) {
  ok <- suppressWarnings(Sys.setlocale("LC_CTYPE", loc))
  if (identical(ok, "")) next
  nm <- names(read.csv(cov_file))[1]
  cat(sprintf("  LC_CTYPE = %-8s first column is named '%s' -> left_join(covariates, by = \"yr\") %s\n",
              loc, nm, if (nm == "yr") "works" else "FAILS"))
}
invisible(Sys.setlocale("LC_CTYPE", old))
cat("  read.csv(..., fileEncoding = \"UTF-8-BOM\") gives '",
    names(read.csv(cov_file, fileEncoding = "UTF-8-BOM"))[1], "' in every locale.\n", sep = "")

## Data exactly as prepared on the 'Load in Data' slide (BOM-safe read)
sox_raw <- read.csv(here::here("data", "sockeye", "SRDATA2026 1.csv"))
stock   <- read.csv(here::here("data", "sockeye", "FC_StockInfo 2.csv"))
covariates <- load_covariates()
sox <- sox_raw %>% left_join(stock, by = "PopID") %>% left_join(covariates, by = "yr") %>%
  filter(!is.na(spawners), !is.na(recruits)) %>% as.data.frame()

section("B. Is `recruits` consistent with the age columns?")
## FC_StockInfo gives each stock's main age classes: c(4,5) for 18 stocks, c(3,4) for Harrison.
chk <- sox %>% filter(!is.na(age2 + age3 + age4 + age5)) %>%
  mutate(expected = ifelse(PopAgeClasses == "c(3,4)", age3 + age4, age4 + age5),
         rel = recruits / expected - 1)
cat("main age classes by stock:\n"); print(table(stock$PopAgeClasses))
cat("rows compared:", nrow(chk), "| rows where `recruits` is within 1% of the sum of the main age classes:",
    sum(abs(chk$rel) <= 0.01), "| rows that differ by more than 1%:", sum(abs(chk$rel) > 0.01), "\n")
odd <- chk %>% filter(abs(rel) > 0.01)
cat("brood years of the rows that differ:\n"); print(table(odd$yr))
cat("number of stocks affected in brood year 2018:", length(unique(odd$PopName[odd$yr == 2018])), "of", length(unique(chk$PopName)), "\n")
cat("largest discrepancies (rel = recruits / sum of main age classes - 1):\n")
print(odd %>% arrange(rel) %>% select(PopName, yr, age3, age4, age5, recruits, rel) %>% head(8) %>%
        mutate(across(c(age3, age4, age5, recruits), ~ round(.x, 4)), rel = round(rel, 2)), row.names = FALSE)
cat("Stellako, brood year 2018 (used in Module 8): age4 =", sox$age4[sox$PopName == "Stellako" & sox$yr == 2018],
    "| age5 =", sox$age5[sox$PopName == "Stellako" & sox$yr == 2018],
    "| recruits =", sox$recruits[sox$PopName == "Stellako" & sox$yr == 2018], "\n")
cat("=> Modules 3-7 model log(recruits/spawners); Module 8 models age4 and age5. For brood year 2018\n",
    "   the two are not the same data.\n", sep = "")
cat("covariate columns are in raw units (not z-scores):\n")
print(round(sapply(covariates[c("peak", "aflow", "mapsst", "jnesst", "pdo", "stn2js")],
                   function(v) c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE))), 2))

## ---- model functions ---------------------------------------------------------
## verbatim from 'Fitting it in RTMB'
negll <- function(pars, data){
  getAll(pars, data, warn = FALSE)

  phi <- -1 + 2*plogis(logitphi)
  sigma <- exp(logsigma)
  beta <- exp(logbeta)
  mu4 <- alpha[1] - beta*spawners
  mu5 <- alpha[2] - beta*spawners + phi*sigma[2]/sigma[1]*(logRS4 - mu4)
  mis4 <- is.na(logRS4)
  mis5 <- is.na(logRS5)
  sigma5 <- sqrt(sigma[2]^2*(1-phi^2))

  ll <- sum(dnorm(logRS4[!mis4], mu4[!mis4], sigma[1], log = TRUE))
  ll <- ll + sum(dnorm(logRS5[!mis5], mu5[!mis5], sigma5, log = TRUE))

  ## Prediction is last year of 4 plus previous year of 5.
  n <- length(mu4)
  Eage4 <- exp(mu4)*spawners
  Eage5 <- exp(mu5)*spawners
  pred <- Eage4[n] + Eage5[n-1]
  REPORT(pred)
  REPORT(Eage4)
  REPORT(Eage5)
  return(-ll)
}
cmb <- function(f, d) function(p) f(p, d)
## verbatim from 'Adding Covariates'
negll_cov <- function(pars, data, X){
  getAll(pars, data, warn = FALSE)

  phi <- -1 + 2*plogis(logitphi)
  sigma <- exp(logsigma)
  beta <- exp(logbeta)
  delta_alpha <- drop(X %*% a)

  mu4 <- delta_alpha + alpha[1] - beta*spawners
  mu5 <- delta_alpha + alpha[2] - beta*spawners + phi*sigma[2]/sigma[1]*(logRS4 - mu4)
  mis4 <- is.na(logRS4)
  mis5 <- is.na(logRS5)
  sigma5 <- sqrt(sigma[2]^2*(1-phi^2))

  ll <- sum(dnorm(logRS4[!mis4], mu4[!mis4], sigma[1], log = TRUE))
  ll <- ll + sum(dnorm(logRS5[!mis5], mu5[!mis5], sigma5, log = TRUE))

  n <- length(mu4)
  Eage4 <- exp(mu4)*spawners
  Eage5 <- exp(mu5)*spawners
  pred <- Eage4[n] + Eage5[n-1]
  REPORT(pred)
  REPORT(Eage4)
  REPORT(Eage5)
  return(-ll)
}
cmb3 <- function(f, d, X) function(p) f(p, d, X)
## General version used for the alternatives below: optional covariates Z (no intercept
## column), optional ridge shrinkage of their coefficients, and the pieces of the
## forecast distribution reported separately.
negll_gen <- function(pars, data, Z = NULL, ridge = FALSE, tv = FALSE){
  getAll(pars, data, warn = FALSE)
  phi <- -1 + 2*plogis(logitphi)
  sigma <- exp(logsigma)
  beta <- exp(logbeta)
  delta <- if (is.null(Z)) 0 else drop(Z %*% a)
  if (tv) delta <- delta + a_t               # time-varying productivity: random walk in log(alpha)
  mu4 <- delta + alpha[1] - beta*spawners
  mu5 <- delta + alpha[2] - beta*spawners + phi*sigma[2]/sigma[1]*(logRS4 - mu4)
  mis4 <- is.na(logRS4); mis5 <- is.na(logRS5)
  sigma5 <- sqrt(sigma[2]^2*(1-phi^2))
  nll <- 0
  if (tv) {                                   # states, not innovations: keeps the Hessian sparse
    tau_a <- exp(logtau_a)
    nll <- nll - dnorm(a_t[1], 0, tau_a, log = TRUE) - sum(dnorm(diff(a_t), 0, tau_a, log = TRUE))
  }
  if (ridge) nll <- nll - sum(dnorm(a, 0, exp(logtau), log = TRUE))
  nll <- nll - sum(dnorm(logRS4[!mis4], mu4[!mis4], sigma[1], log = TRUE))
  nll <- nll - sum(dnorm(logRS5[!mis5], mu5[!mis5], sigma5, log = TRUE))
  n <- length(mu4)
  m4 <- mu4[n] + log(spawners[n]); m5 <- mu5[n-1] + log(spawners[n-1])   # log-scale means of the two pieces
  s4 <- sigma[1]; s5 <- sigma5
  pred <- exp(m4) + exp(m5)
  fc <- c(m4, m5)
  ADREPORT(fc)                                # so sdreport() returns their joint covariance
  REPORT(pred); REPORT(m4); REPORT(m5); REPORT(s4); REPORT(s5)
  nll
}
gen <- function(d, Z = NULL, ridge = FALSE, tv = FALSE) function(p) negll_gen(p, d, Z, ridge, tv)

prep_stock <- function(name, const = 1e-4) {
  sox %>% filter(PopName == name) %>%
    mutate(logRS4 = log((age4 + const)/spawners), logRS5 = log((age5 + const)/spawners))
}
sox_ch <- prep_stock("Stellako")
params <- list(logsigma = c(0,0), alpha = c(0,0), logbeta = 0, logitphi = 0)
covs <- c("peak", "aflow", "jnesst", "pdo", "stn2js")

section("C. 'Adding Covariates' as written")
X <- model.matrix(~peak + aflow + mapsst + pdo + stn2js, data = sox_ch)
cat("columns of X:", paste(colnames(X), collapse = ", "), "\n")
set.seed(1)
params_cov <- list(logsigma = c(0,0), alpha = c(0,0), a = rnorm(ncol(X)), logbeta = 0, logitphi = 0)
obj_c <- MakeADFun(cmb3(negll_cov, sox_ch, X), params_cov, silent = TRUE)
fit_c <- suppressWarnings(nlminb(obj_c$par, obj_c$fn, obj_c$gr))
cat("nlminb as on the slide (default control): convergence =", fit_c$convergence, "|", fit_c$message,
    "| iterations", fit_c$iterations, "\n")
fit_c2 <- suppressWarnings(nlminb(obj_c$par, obj_c$fn, obj_c$gr, control = big))
cat("with eval.max = iter.max = 1e4:           convergence =", fit_c2$convergence, "|", fit_c2$message,
    "| iterations", fit_c2$iterations, "| NLL", round(fit_c2$objective, 4), "\n")
H <- obj_c$he(fit_c2$par); ev <- eigen(H, symmetric = TRUE)
cat("eigenvalues of the Hessian (smallest three):", paste(signif(tail(ev$values, 3), 3), collapse = ", "), "\n")
v <- ev$vectors[, ncol(H)]; names(v) <- names(fit_c2$par)
cat("direction of the (near-)zero eigenvalue:\n"); print(round(v[abs(v) > 0.05], 3))
sdr_c <- suppressWarnings(try(sdreport(obj_c), silent = TRUE))
cat("sdreport(obj)$pdHess:", if (inherits(sdr_c, "try-error")) "sdreport failed" else sdr_c$pdHess, "\n")
se <- if (inherits(sdr_c, "try-error")) NA else suppressWarnings(sqrt(diag(sdr_c$cov.fixed)))
cat("standard errors reported for alpha[1], alpha[2], a[1] (intercept):",
    paste(signif(se[c(which(names(fit_c2$par) == "alpha"), which(names(fit_c2$par) == "a")[1])], 3), collapse = ", "), "\n")
cat("=> a[1] (the '(Intercept)' column of X) and alpha[1], alpha[2] enter only through a[1]+alpha[1] and\n",
    "   a[1]+alpha[2]: three parameters for two identifiable quantities.\n", sep = "")

## Same model, identifiable: no intercept column, covariates centred and scaled.
Zs <- scale(as.matrix(sox_ch[, c("peak", "aflow", "mapsst", "pdo", "stn2js")]))
p_gen <- list(logsigma = c(0,0), alpha = c(0,0), a = rep(0, ncol(Zs)), logbeta = 0, logitphi = 0)
obj_g <- MakeADFun(gen(sox_ch, Zs), p_gen, silent = TRUE)
fit_g <- nlminb(obj_g$par, obj_g$fn, obj_g$gr, control = big)
sdr_g <- sdreport(obj_g)
cat("\nIdentifiable version (no intercept in the covariate matrix, covariates standardised):\n")
cat("  NLL", round(fit_g$objective, 4), "(as written:", round(fit_c2$objective, 4), ") | pdHess", sdr_g$pdHess,
    "| iterations", fit_g$iterations, "| condition number of Hessian",
    signif(kappa(obj_g$he(fit_g$par), exact = TRUE), 3), "\n")
est <- summary(sdr_g)[rownames(summary(sdr_g)) == "a", ]
rownames(est) <- colnames(Zs)
cat("  covariate effects per 1 SD (estimate, SE, z):\n")
print(round(cbind(est, z = est[, 1] / est[, 2]), 3))
obj_0 <- MakeADFun(gen(sox_ch), params, silent = TRUE)
fit_0 <- nlminb(obj_0$par, obj_0$fn, obj_0$gr, control = big)
cat(sprintf("  no-covariate model NLL %.3f; 5 covariates lower it by %.3f -> AIC changes by %+.2f\n",
            fit_0$objective, fit_0$objective - fit_g$objective, 2 * 5 - 2 * (fit_0$objective - fit_g$objective)))

section("D. Sensitivity to the constant added before taking logs")
sens <- function(name, consts = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) {
  out <- t(sapply(consts, function(cc) {
    d <- prep_stock(name, cc)
    o <- MakeADFun(gen(d), params, silent = TRUE)
    f <- suppressWarnings(nlminb(o$par, o$fn, o$gr, control = big))
    r <- o$report(f$par)
    c(const = cc, alpha5 = unname(f$par[names(f$par) == "alpha"][2]), sigma5_marginal = unname(exp(f$par[names(f$par) == "logsigma"][2])),
      phi = unname(2 * plogis(f$par["logitphi"]) - 1), age5_forecast = exp(r$m5), total_forecast = r$pred)
  }))
  d0 <- sox %>% filter(PopName == name)
  cat(sprintf("%s: %d zero age-5 values of %d; smallest positive age-5 = %.6f (millions)\n",
              name, sum(d0$age5 == 0, na.rm = TRUE), sum(!is.na(d0$age5)), min(d0$age5[d0$age5 > 0], na.rm = TRUE)))
  print(signif(out, 3))
}
sens("Stellako")
sens("Quesnel")
sens("Late Shuswap")

section("E. Median vs mean forecast (Stellako, no covariates, as on the slide)")
r0 <- obj_0$report(fit_0$par)
cat(sprintf("age-4 piece: median %.4f, mean %.4f (x %.3f) | age-5 piece: median %.4f, mean %.4f (x %.3f)\n",
            exp(r0$m4), exp(r0$m4 + r0$s4^2 / 2), exp(r0$s4^2 / 2),
            exp(r0$m5), exp(r0$m5 + r0$s5^2 / 2), exp(r0$s5^2 / 2)))
cat(sprintf("slide's `pred` (sum of two medians) = %.4f | expected return = %.4f | ratio %.2f\n",
            r0$pred, exp(r0$m4 + r0$s4^2 / 2) + exp(r0$m5 + r0$s5^2 / 2),
            (exp(r0$m4 + r0$s4^2 / 2) + exp(r0$m5 + r0$s5^2 / 2)) / r0$pred))

section("F. Retrospective evaluation, return years 2004-2025 (22 forecasts)")
## Observed returns by return year = age-4 fish from brood Y-4 + age-5 fish from brood Y-5
ret <- merge(data.frame(year = sox_ch$yr + 4, a4 = sox_ch$age4), data.frame(year = sox_ch$yr + 5, a5 = sox_ch$age5))
ret$obs <- ret$a4 + ret$a5
ret <- ret[!is.na(ret$obs), c("year", "obs")]
years <- 2000:2021

hindcast <- function(kind) {
  out <- lapply(years, function(by) {
    dat <- sox_ch %>% filter(yr <= by)
    dat$logRS4[dat$yr == by] <- NA
    dat$logRS5[dat$yr == by] <- NA
    dat$logRS5[dat$yr == by - 1] <- NA
    if (kind == "slide_base") {
      o <- MakeADFun(cmb(negll, dat), params, silent = TRUE)
      f <- suppressWarnings(nlminb(o$par, o$fn, o$gr))
      return(data.frame(year = by + 4, pred = o$report()$pred, conv = f$convergence))
    }
    if (kind == "slide_cov") {
      Xd <- model.matrix(~ peak + aflow + jnesst + pdo + stn2js, data = dat)
      o <- MakeADFun(cmb3(negll_cov, dat, Xd), params_cov, silent = TRUE)
      f <- suppressWarnings(nlminb(o$par, o$fn, o$gr))
      return(data.frame(year = by + 4, pred = o$report()$pred, conv = f$convergence))
    }
    ## alternatives: covariates centred/scaled with TRAINING-window statistics only
    Zd <- scale(as.matrix(dat[, covs]))
    o <- switch(kind,
      cov_fixed = MakeADFun(gen(dat, Zd), p_gen, silent = TRUE),
      cov_ridge = MakeADFun(gen(dat, Zd, ridge = TRUE), c(p_gen, list(logtau = log(0.2))),
                            random = "a", silent = TRUE),
      base      = MakeADFun(gen(dat), params, silent = TRUE),
      rw        = MakeADFun(gen(dat, tv = TRUE), c(params, list(a_t = rep(0, nrow(dat)), logtau_a = log(0.2))),
                            random = "a_t", silent = TRUE))
    f <- suppressWarnings(nlminb(o$par, o$fn, o$gr, control = big))
    r <- o$report(o$env$last.par.best)
    q <- rep(NA_real_, 5)
    if (kind %in% c("base", "rw")) {
      ## Forecast distribution = estimation uncertainty in the two log-scale means (delta method,
      ## including the random effects for the random-walk model) + lognormal process error.
      sdr <- sdreport(o)
      mm  <- MASS::mvrnorm(4000, sdr$value, sdr$cov)
      tot <- exp(rnorm(4000, mm[, 1], r$s4)) + exp(rnorm(4000, mm[, 2], r$s5))
      q <- quantile(tot, c(0.1, 0.25, 0.5, 0.75, 0.9))
    }
    data.frame(year = by + 4, pred = r$pred, conv = f$convergence,
               mean_pred = exp(r$m4 + r$s4^2 / 2) + exp(r$m5 + r$s5^2 / 2),
               q10 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q90 = q[5],
               tau_a = if (kind == "rw") unname(exp(f$par["logtau_a"])) else NA_real_)
  })
  do.call(rbind, out)
}
set.seed(2026)
h <- list(slide_base = hindcast("slide_base"), slide_cov = hindcast("slide_cov"),
          cov_fixed = hindcast("cov_fixed"), cov_ridge = hindcast("cov_ridge"),
          base_dist = hindcast("base"), rw = hindcast("rw"))
cat("non-converged refits (nlminb code != 0) out of 22:  ",
    paste(names(h), sapply(h, function(d) sum(d$conv != 0)), sep = " = ", collapse = " | "), "\n")
cat("largest |difference| between the slide's covariate forecasts and the identifiable version of the same model:",
    signif(max(abs(h$slide_cov$pred - h$cov_fixed$pred)), 3), "(millions)\n")

## Naive benchmarks built only from returns observed before the forecast year
obs_of <- function(y) ret$obs[match(y, ret$year)]
naive <- data.frame(year = years + 4)
naive$R1C <- obs_of(naive$year - 4)
naive$R2C <- (obs_of(naive$year - 4) + obs_of(naive$year - 8)) / 2
naive$RAC <- sapply(naive$year, function(y) mean(ret$obs[ret$year < y & (y - ret$year) %% 4 == 0]))
naive$TSA <- sapply(naive$year, function(y) mean(ret$obs[ret$year < y]))
naive$obs <- obs_of(naive$year)

score <- function(pred, obs) c(MRE = mean(pred - obs), MAE = mean(abs(pred - obs)), RMSE = sqrt(mean((pred - obs)^2)),
                               MPE = mean((pred - obs) / obs), log_RMSE = sqrt(mean(log(pred / obs)^2)))
tab <- rbind(
  "slide: age-structured Ricker"                   = score(h$slide_base$pred, naive$obs),
  "slide: + 5 covariates"                          = score(h$slide_cov$pred, naive$obs),
  "5 covariates, ridge (random-effect) shrinkage"  = score(h$cov_ridge$pred, naive$obs),
  "age-structured Ricker, mean instead of median"  = score(h$base_dist$mean_pred, naive$obs),
  "random-walk productivity (state-space)"         = score(h$rw$pred, naive$obs),
  "naive R1C (return 4 years earlier)"             = score(naive$R1C, naive$obs),
  "naive R2C (mean of 4 and 8 years earlier)"      = score(naive$R2C, naive$obs),
  "naive RAC (mean of earlier same-cycle returns)" = score(naive$RAC, naive$obs),
  "naive TSA (mean of all earlier returns)"        = score(naive$TSA, naive$obs))
tab <- cbind(tab, MASE_vs_R1C = tab[, "MAE"] / tab["naive R1C (return 4 years earlier)", "MAE"])
cat("\nForecast skill, Stellako, 22 return years (units: millions of fish; MRE = mean raw error = bias):\n")
print(round(tab, 3))
cat("mean observed return over the 22 years:", round(mean(naive$obs), 3), "\n")

cat("\nCalibration of the forecast distributions (estimation + process uncertainty), 22 years:\n")
calib <- function(d) c("inside 50% interval (expect 11)"   = sum(naive$obs >= d$q25 & naive$obs <= d$q75),
                       "inside 80% interval (expect 17.6)" = sum(naive$obs >= d$q10 & naive$obs <= d$q90),
                       "below 10% point (expect 2.2)"      = sum(naive$obs < d$q10),
                       "below median (expect 11)"          = sum(naive$obs < d$q50),
                       "above 90% point (expect 2.2)"      = sum(naive$obs > d$q90))
print(rbind("stationary age-structured Ricker" = calib(h$base_dist),
            "random-walk productivity"         = calib(h$rw)))
## Interval score for the central 80% interval (Gneiting & Raftery 2007): lower is better
int_score <- function(d, a = 0.2) mean((d$q90 - d$q10) + (2 / a) * pmax(d$q10 - naive$obs, 0) +
                                         (2 / a) * pmax(naive$obs - d$q90, 0))
cat("mean 80% interval score (lower is better): stationary =", round(int_score(h$base_dist), 3),
    "| random walk =", round(int_score(h$rw), 3), "\n")
cat("random-walk SD of log-productivity (tau_a) across the 22 refits: range",
    paste(round(range(h$rw$tau_a), 3), collapse = " to "), "\n")
cat("\nForecast vs observed by year (millions):\n")
print(data.frame(year = naive$year, observed = round(naive$obs, 3), slide_ricker = round(h$slide_base$pred, 3),
                 random_walk = round(h$rw$pred, 3), naive_R1C = round(naive$R1C, 3)), row.names = FALSE)
