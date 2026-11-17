## Module 5 says "With more data, the likelihood is smooth and has a sharp peak" when discussing
## when the Laplace approximation works. What matters is information PER RANDOM EFFECT.
## Binary responses, random intercept per group, true SD = 1.5, 300 groups:
##   Laplace (RTMB, random = "u")  vs  the exact marginal likelihood (adaptive Gauss-Hermite quadrature).
## The total number of observations is varied only through the number of observations per group.

source(here::here("review", "verification", "00_setup.R"))

## Gauss-Hermite nodes/weights for integrals against the standard normal (Golub-Welsch)
gauss_hermite <- function(n) {
  i <- seq_len(n - 1)
  e <- eigen(diag(0, n) + rbind(cbind(0, diag(sqrt(i / 2))), 0) + t(rbind(cbind(0, diag(sqrt(i / 2))), 0)),
             symmetric = TRUE)
  list(z = sqrt(2) * e$values, w = e$vectors[1, ]^2)      # nodes for N(0,1), weights sum to 1
}
gh <- gauss_hermite(41)
stopifnot(abs(sum(gh$w) - 1) < 1e-10, abs(sum(gh$w * gh$z^2) - 1) < 1e-8)

fit_both <- function(seed, n_groups, n_per, sd_true = 1.5, b0_true = -1) {
  set.seed(seed)
  u <- rnorm(n_groups, 0, sd_true)
  s <- rbinom(n_groups, n_per, plogis(b0_true + u))       # successes per group
  ## Laplace
  f_lap <- function(p) {
    nll <- -sum(dnorm(p$u, 0, 1, log = TRUE))
    nll - sum(dbinom(s, n_per, plogis(p$b0 + exp(p$logsd) * p$u), log = TRUE))
  }
  o1 <- MakeADFun(f_lap, list(b0 = 0, logsd = 0, u = rep(0, n_groups)), random = "u", silent = TRUE)
  r1 <- suppressWarnings(nlminb(o1$par, o1$fn, o1$gr, control = list(eval.max = 1e4, iter.max = 1e4)))
  ## Exact marginal likelihood by ADAPTIVE Gauss-Hermite quadrature (nodes centred and scaled at
  ## each group's own mode), in plain R. Groups with the same count share one integral.
  tab <- table(factor(s, levels = 0:n_per)); cnt <- as.numeric(tab); val <- 0:n_per
  log_int <- function(j, b0, sdv) {
    lf <- function(z) dbinom(j, n_per, plogis(b0 + sdv * z), log = TRUE) + dnorm(z, log = TRUE)
    zh <- optimize(lf, c(-15, 15), maximum = TRUE, tol = 1e-12)$maximum
    ph <- plogis(b0 + sdv * zh)
    sh <- 1 / sqrt(n_per * sdv^2 * ph * (1 - ph) + 1)          # exact curvature at the mode
    zk <- zh + sh * gh$z                                        # shifted and scaled nodes
    lf(zh) + log(sh) + log(sum(gh$w * exp(lf(zk) - lf(zh) - dnorm(gh$z, log = TRUE))))   # derivation below
  }
  ## Derivation: int f(z) dz = sh * int f(zh + sh x) dx = sh * E_{x ~ N(0,1)}[ f(zh + sh x) / dnorm(x) ],
  ## and the expectation is the Gauss-Hermite sum with weights gh$w at nodes gh$z.
  f_ex <- function(par) {
    -sum(sapply(which(cnt > 0), function(i) cnt[i] * log_int(val[i], par[1], exp(par[2]))))
  }
  r2 <- suppressWarnings(nlminb(c(0, 0), f_ex, control = list(eval.max = 1e4, iter.max = 1e4)))
  names(r2$par) <- c("b0", "logsd")
  c(sd_laplace = unname(exp(r1$par["logsd"])), sd_exact = unname(exp(r2$par["logsd"])),
    b0_laplace = unname(r1$par["b0"]), b0_exact = unname(r2$par["b0"]))
}

## Self-test of the adaptive quadrature against stats::integrate() for one sharply peaked integrand
local({
  n_per <- 100; j <- 30; b0 <- -1; sdv <- 1.5
  lf <- function(z) dbinom(j, n_per, plogis(b0 + sdv * z), log = TRUE) + dnorm(z, log = TRUE)
  zh <- optimize(lf, c(-15, 15), maximum = TRUE, tol = 1e-12)$maximum
  ph <- plogis(b0 + sdv * zh); sh <- 1 / sqrt(n_per * sdv^2 * ph * (1 - ph) + 1)
  aghq <- lf(zh) + log(sh) + log(sum(gh$w * exp(lf(zh + sh * gh$z) - lf(zh) - dnorm(gh$z, log = TRUE))))
  ref  <- log(integrate(function(z) exp(lf(z)), zh - 12 * sh, zh + 12 * sh, rel.tol = 1e-12)$value)
  fixed <- log(sum(gh$w * dbinom(j, n_per, plogis(b0 + sdv * gh$z))))          # non-adaptive nodes
  cat(sprintf("quadrature self-test (100 obs in the group): adaptive %.8f | integrate() %.8f | fixed nodes %.8f\n",
              aghq, ref, fixed))
})

section("Laplace vs exact ML for a binary random-intercept model (true SD 1.5, intercept -1, 300 groups)")
res <- lapply(c(2, 5, 20, 100), function(n_per) {
  r <- t(sapply(1:30, fit_both, n_groups = 300, n_per = n_per))
  data.frame(obs_per_group = n_per, total_obs = 300 * n_per,
             sd_exact = mean(r[, "sd_exact"]), sd_laplace = mean(r[, "sd_laplace"]),
             ratio_sd = mean(r[, "sd_laplace"] / r[, "sd_exact"]),
             b0_exact = mean(r[, "b0_exact"]), b0_laplace = mean(r[, "b0_laplace"]))
})
res <- do.call(rbind, res)
print(round(res, 3), row.names = FALSE)
cat("(means over 30 simulated data sets; ratio_sd = Laplace estimate / exact ML estimate of the SD)\n")

section("Same total number of observations (6000), arranged differently")
a <- t(sapply(1:30, fit_both, n_groups = 3000, n_per = 2))
b <- t(sapply(1:30, fit_both, n_groups = 60, n_per = 100))
print(round(rbind("3000 groups x 2 obs"  = c(ratio_sd = mean(a[, "sd_laplace"] / a[, "sd_exact"]),
                                              sd_exact = mean(a[, "sd_exact"]), sd_laplace = mean(a[, "sd_laplace"])),
                  "60 groups x 100 obs"   = c(ratio_sd = mean(b[, "sd_laplace"] / b[, "sd_exact"]),
                                              sd_exact = mean(b[, "sd_exact"]), sd_laplace = mean(b[, "sd_laplace"]))), 3))
cat("=> adding groups does not help the approximation; adding observations per group does.\n")
