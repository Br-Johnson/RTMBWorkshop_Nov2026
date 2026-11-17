## Tested versions of the workflow recommended in the review document.
##  A. A fitting helper: generous limits, optional Newton steps, and the checks that matter.
##  B. Process residuals from ONE joint draw of the random effects (checks the random-effect model).
##  C. Simulate-and-refit self-test: is the between-stock SD estimated without bias? (ML vs REML)

source(here::here("review", "verification", "00_setup.R"))
suppressPackageStartupMessages(library(Matrix))
sox <- load_sox()
nstocks <- max(sox$stock); nyears <- max(sox$time)

## ---- A. fitting helper --------------------------------------------------------
fit_rtmb <- function(obj, newton_steps = 2, control = list(eval.max = 1e4, iter.max = 1e4), quiet = FALSE) {
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
  for (i in seq_len(newton_steps)) {                 # polish: Newton steps on the (marginal) likelihood
    g <- as.numeric(obj$gr(opt$par))
    h <- optimHess(opt$par, obj$fn, obj$gr)
    step <- try(solve(h, g), silent = TRUE)
    if (inherits(step, "try-error")) break
    new_par <- opt$par - step
    new_obj <- obj$fn(new_par)
    if (is.finite(new_obj) && new_obj <= opt$objective) {   # keep the step only if it does not make things worse
      opt$par <- new_par; opt$objective <- new_obj
    } else break
  }
  invisible(obj$fn(opt$par))                         # make sure last.par.best holds the final point
  sdr <- sdreport(obj, par.fixed = opt$par)
  ev  <- eigen(optimHess(opt$par, obj$fn, obj$gr), symmetric = TRUE, only.values = TRUE)$values
  se  <- suppressWarnings(sqrt(diag(sdr$cov.fixed)))
  checks <- c(nlminb_code_0     = opt$convergence == 0,
              max_gradient_ok   = max(abs(sdr$gradient.fixed)) < 1e-3,
              hessian_pos_def   = all(ev > 0),
              all_se_finite     = all(is.finite(se)),
              no_huge_se        = all(se[is.finite(se)] < 100))
  if (!quiet) {
    cat(sprintf("objective %.4f | nlminb: %s | max |gradient| %.2e | smallest Hessian eigenvalue %.3g\n",
                opt$objective, opt$message, max(abs(sdr$gradient.fixed)), min(ev)))
    cat("checks passed: ", paste(names(checks)[checks], collapse = ", "), "\n")
    if (any(!checks)) cat("CHECKS FAILED: ", paste(names(checks)[!checks], collapse = ", "), "\n")
  }
  list(opt = opt, sdr = sdr, checks = checks, par_full = obj$env$last.par.best)
}

init <- list(alpha0 = mean(sox$logRS), re = rep(0, nstocks), logbeta = rep(log(0.5), nstocks),
             logsigma = rep(log(sd(sox$logRS)), nstocks), logsigma_re = log(0.5),
             re_time = rep(0, nyears), logitphi = 0, logsigma_time = log(0.5))
f_ar <- function(pars){
  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta); sigma <- exp(logsigma)
  sigma_time <- exp(logsigma_time); phi <- 2*plogis(logitphi)-1
  logRS <- OBS(logRS)
  mu <- (re + alpha0)[stock] + re_time[time] - beta[stock]*spawners_scaled
  nll <- -dautoreg(re_time, phi = phi, scale = sqrt(sigma_time^2/(1-phi^2)), log = TRUE)
  nll <- nll - sum(dnorm(re, 0, exp(logsigma_re), log = TRUE))
  nll - sum(dnorm(logRS, mu, sigma[stock], log = TRUE))
}
section("A. Fitting helper on the AR(1) hierarchical Ricker")
obj <- MakeADFun(f_ar, init, random = c("re", "re_time"), silent = TRUE)
plain <- nlminb(obj$par, obj$fn, obj$gr)
cat("single default nlminb() call as on the slides:", plain$message, "| iterations", plain$iterations,
    "| max |gradient|", signif(max(abs(obj$gr(plain$par))), 2), "\n")
obj <- MakeADFun(f_ar, init, random = c("re", "re_time"), silent = TRUE)
fit <- fit_rtmb(obj)

section("B. Process residuals from one joint sample of the random effects")
## Joint precision of (fixed, random) at the optimum; draw the random effects given the fixed effects.
sdr_j <- sdreport(obj, par.fixed = fit$opt$par, getJointPrecision = TRUE)
r_idx <- obj$env$random
Quu <- sdr_j$jointPrecision[r_idx, r_idx]
L <- Matrix::Cholesky(Quu, LDL = FALSE)
pl <- obj$env$parList(fit$opt$par, fit$par_full)
nm <- names(fit$par_full)[r_idx]
phi <- 2 * plogis(pl$logitphi) - 1; s_t <- exp(pl$logsigma_time)
std_innov <- function(a) c(a[1] * sqrt(1 - phi^2), a[-1] - phi * a[-length(a)]) / s_t   # iid N(0,1) if the AR(1) is right
draw_once <- function() {
  z <- rnorm(length(r_idx))
  u <- fit$par_full[r_idx] + as.numeric(Matrix::solve(L, Matrix::solve(L, z, system = "Lt"), system = "Pt"))
  innov <- std_innov(u[nm == "re_time"]); stock_z <- u[nm == "re"] / exp(pl$logsigma_re)
  half <- seq_along(innov) <= length(innov) / 2
  c(sd_innov = sd(innov), shapiro_p_innov = shapiro.test(innov)$p.value,
    lag1_innov = acf(innov, plot = FALSE)$acf[2],
    first_half_minus_second = mean(innov[half]) - mean(innov[!half]),
    min_innov = min(innov), sd_stock = sd(stock_z), shapiro_p_stock = shapiro.test(stock_z)$p.value)
}
set.seed(11)
one <- draw_once()
cat("one draw (this is what you would plot):\n"); print(round(one, 3))
many <- t(replicate(500, draw_once()))
cat("\nthe same summaries over 500 independent draws (to show the single-draw conclusions are not luck):\n")
print(round(rbind(mean = colMeans(many), sd = apply(many, 2, sd)), 3))
cat("share of draws in which Shapiro-Wilk rejects normality at 5%: year-effect innovations",
    round(mean(many[, "shapiro_p_innov"] < 0.05), 2), "| stock effects", round(mean(many[, "shapiro_p_stock"] < 0.05), 2), "\n")
cat("expected SD of 'first half minus second half' if the innovations were iid N(0,1):",
    round(sqrt(2 / (nyears / 2)), 2), "\n")
## The same check with the conditional modes instead of a draw is too optimistic:
a_hat <- fit$par_full[r_idx][nm == "re_time"]
cat("SD of the innovations computed from the conditional modes (shrunk towards zero):", round(sd(std_innov(a_hat)), 2), "\n")
yrs <- min(sox$yr) + seq_len(nyears) - 1
cat("brood year with the most negative year effect (conditional mode):", yrs[which.min(a_hat)],
    "| value", round(min(a_hat), 2), "\n")

section("C. How well is the between-stock SD identified? Profile + simulate-and-refit (ML vs REML)")
## Hierarchical Ricker from Module 5 (random stock intercepts, stock-specific slopes and residual SDs).
## `linear = TRUE` writes the slopes without exp(), so that they can be integrated out for full REML.
f_h <- function(dat, linear = FALSE) function(pars){
  getAll(pars, dat, warn = FALSE)
  slope <- if (linear) b else exp(logbeta)
  sigma <- exp(logsigma)
  mu <- (re + alpha0)[stock] - slope[stock]*spawners_scaled
  nll <- -sum(dnorm(re, 0, exp(logsigma_re), log = TRUE))
  nll - sum(dnorm(logRS, mu, sigma[stock], log = TRUE))
}
init_h <- init[c("alpha0", "re", "logbeta", "logsigma", "logsigma_re")]
big <- list(eval.max = 1e4, iter.max = 1e4)
o_h <- MakeADFun(f_h(sox), init_h, random = "re", silent = TRUE)
f0 <- nlminb(o_h$par, o_h$fn, o_h$gr, control = big)
truth <- o_h$env$parList(f0$par)
sd_hat <- exp(truth$logsigma_re)
se_log <- sqrt(diag(sdreport(o_h)$cov.fixed))[names(f0$par) == "logsigma_re"]
cat(sprintf("ML estimate on the real data: sigma_re = %.4f | Wald 95%% interval on the log scale: %.3f to %.3f\n",
            sd_hat, exp(log(sd_hat) - 1.96 * se_log), exp(log(sd_hat) + 1.96 * se_log)))
grid <- c(1e-4, 0.02, 0.05, sd_hat, 0.15, 0.2, 0.3)
prof <- sapply(grid, function(sv) {
  st <- truth; st$logsigma_re <- log(sv)
  o <- MakeADFun(f_h(sox), st, random = "re", map = list(logsigma_re = factor(NA)), silent = TRUE)
  nlminb(o$par, o$fn, o$gr, control = big)$objective })
cat("profile of the marginal NLL over sigma_re (difference from its minimum; 1.92 = edge of a 95% interval):\n")
print(rbind(sigma_re = signif(grid, 3), delta_NLL = round(prof - min(prof), 2)))
cat("=> sigma_re = 0 is only", round(prof[1] - min(prof), 2), "log-likelihood units worse than the estimate:",
    "the data barely distinguish 'no between-stock variation' from the fitted value.\n\n")

one_sim <- function(seed, sd_gen) {
  set.seed(seed)
  d <- sox
  re_sim <- rnorm(nstocks, 0, sd_gen)
  mu <- (truth$alpha0 + re_sim)[d$stock] - exp(truth$logbeta)[d$stock] * d$spawners_scaled
  d$logRS <- rnorm(nrow(d), mu, exp(truth$logsigma)[d$stock])
  st <- truth; st$re <- rep(0, nstocks); st$logsigma_re <- log(max(sd_gen, 0.05))
  st_lin <- st; st_lin$b <- exp(st$logbeta); st_lin$logbeta <- NULL
  est <- function(o) unname(exp(suppressWarnings(nlminb(o$par, o$fn, o$gr, control = big))$par["logsigma_re"]))
  c(ML          = est(MakeADFun(f_h(d), st, random = "re", silent = TRUE)),
    REML_alpha0 = est(MakeADFun(f_h(d), st, random = c("re", "alpha0"), silent = TRUE)),
    REML_full   = est(MakeADFun(f_h(d, linear = TRUE), st_lin, random = c("re", "alpha0", "b"), silent = TRUE)))
}
t0 <- proc.time()["elapsed"]
for (sd_gen in c(sd_hat, 0.3)) {
  sims <- t(sapply(1:150, one_sim, sd_gen = sd_gen))
  cat(sprintf("data generated with sigma_re = %.3f, 150 simulated data sets:\n", sd_gen))
  print(round(rbind(mean_estimate = colMeans(sims), median_estimate = apply(sims, 2, median),
                    share_estimated_at_zero = colMeans(sims < 0.01)), 3))
}
cat("simulation time:", round(proc.time()["elapsed"] - t0), "s\n")
cat("REML_alpha0 integrates the common intercept only (random = c('re', 'alpha0'));\n",
    "REML_full also integrates the 19 slopes, written without exp() so the model is linear in them.\n", sep = "")
