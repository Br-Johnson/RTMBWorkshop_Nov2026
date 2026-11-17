## Modules 3-4 (modules/3_gradients.qmd, modules/4_inference.qmd).
##  A. `beta_real <- beta * transform_beta` is not a back-transformation.
##  B. How do the random starting values (logsigma ~ U(3, 5)) behave?
##  C. What does confint(tmbprofile) actually compute?
##  D. Does the DataEval bootstrap refit correctly, and when does it go stale?
##  E. solve(H) succeeding does not show H is positive definite.

source(here::here("review", "verification", "00_setup.R"))
sox <- load_sox()
nstocks <- max(sox$stock)
sd_by_stock <- tapply(sox$sd_spawners, sox$PopName, unique)
cat("per-stock SD of spawners (millions): smallest", signif(min(sd_by_stock), 3), "(", names(which.min(sd_by_stock)),
    ") largest", signif(max(sd_by_stock), 3), "(", names(which.max(sd_by_stock)), ") ratio", round(max(sd_by_stock) / min(sd_by_stock)), "\n")

## Objective exactly as on the Module 3/4 slides (data read from the workspace).
negll <- function(pars){
  getAll(pars)
  beta <- exp(logbeta)
  sigma <- exp(logsigma)
  mu <- alpha[sox$stock] - beta[sox$stock]*sox$spawners_scaled
  -sum(dnorm(sox$logRS, mu, sigma[sox$stock], log = TRUE))
}
workshop_inits <- function() list(
  alpha = runif(nstocks, 2, 5),
  logbeta = runif(nstocks, 0, 1),
  logsigma = runif(nstocks, 3, 5)
)

## Reference solution: each stock is a separate linear regression.
ref <- lapply(split(sox, sox$stock), function(d) {
  m <- lm(logRS ~ spawners_scaled, data = d)
  c(alpha = unname(coef(m)[1]), beta = -unname(coef(m)[2]),
    sigma_ml = sqrt(mean(resid(m)^2)), nll = -as.numeric(logLik(m)))
})
ref <- do.call(rbind, ref)
nll_ref <- sum(ref[, "nll"])
cat("Reference NLL from 19 separate lm() fits:", nll_ref, "\n")
cat("Stocks whose unconstrained OLS slope implies beta <= 0 (no compensation):",
    paste(levels(factor(sox$PopName))[ref[, "beta"] <= 0], collapse = ", "), "\n")

section("A. 'Transform Parameters' slide: beta_real <- beta * transform_beta")
set.seed(101)
parameters <- workshop_inits()
transform_beta <- sox$spawners/sox$sd_spawners          # verbatim from the slide
negll_A <- function(pars){
  getAll(pars)
  beta <- exp(logbeta)
  sigma <- exp(logsigma)
  beta_real <- beta*transform_beta                      # verbatim from the slide
  REPORT(beta_real)                                     # added so we can look at it
  ADREPORT(beta)
  mu <- alpha[sox$stock] - beta[sox$stock]*sox$spawners_scaled
  -sum(dnorm(sox$logRS, mu, sigma[sox$stock], log = TRUE))
}
objA <- MakeADFun(negll_A, parameters, silent = TRUE)
fitA <- nlminb(objA$par, objA$fn, objA$gr, control = list(eval.max = 1e4, iter.max = 1e4))
cat("length(beta) =", nstocks, "| length(transform_beta) =", length(transform_beta),
    "| length(beta_real) =", length(objA$report(fitA$par)$beta_real), "\n")
cat("transform_beta is identical to spawners_scaled:",
    isTRUE(all.equal(transform_beta, sox$spawners_scaled)), "\n")
cat("=> a length-19 vector is silently recycled against one value per observation;\n",
    "   `beta_real` is never reported, and ADREPORT(beta) is still on the scaled-spawner scale.\n", sep = "")

## Correct back-transformation: beta_real[s] = beta_scaled[s] / sd_spawners[s]
sd_s <- as.numeric(tapply(sox$sd_spawners, sox$stock, unique))
negll_A2 <- function(pars){
  getAll(pars)
  beta <- exp(logbeta)
  sigma <- exp(logsigma)
  beta_real <- beta / sd_s
  ADREPORT(beta_real)
  mu <- alpha[sox$stock] - beta[sox$stock]*sox$spawners_scaled
  -sum(dnorm(sox$logRS, mu, sigma[sox$stock], log = TRUE))
}
objA2 <- MakeADFun(negll_A2, parameters, silent = TRUE)
fitA2 <- nlminb(objA2$par, objA2$fn, objA2$gr, control = list(eval.max = 1e4, iter.max = 1e4))
sdrA2 <- sdreport(objA2)
beta_real_hat <- summary(sdrA2, "report")[, "Estimate"]
lm_raw <- sapply(split(sox, sox$stock), function(d) -coef(lm(logRS ~ spawners, data = d))[2])
ok <- ref[, "beta"] > 1e-6
cat("max relative difference between corrected beta_real and -slope of lm(logRS ~ spawners)\n",
    "  (stocks with positive OLS beta only):",
    signif(max(abs(beta_real_hat[ok] / lm_raw[ok] - 1)), 3), "\n")

section("B. Starting values used on the slides")
run_start <- function(seed, control = list()) {
  set.seed(seed)
  p <- workshop_inits()
  obj <- MakeADFun(negll, p, silent = TRUE)
  fit <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr, control = control))
  data.frame(seed = seed, conv = fit$convergence, iters = fit$iterations,
             dnll = fit$objective - nll_ref,
             maxgrad = max(abs(obj$gr(fit$par))), msg = fit$message)
}
resB <- do.call(rbind, lapply(1:50, run_start))
cat("Workshop inits (alpha~U(2,5), logbeta~U(0,1), logsigma~U(3,5)), default nlminb control, 50 seeds:\n")
cat("  sigma implied by logsigma in (3,5): from", round(exp(3), 1), "to", round(exp(5), 1),
    " (observed sd of logRS =", round(sd(sox$logRS), 2), ")\n")
print(table(convergence_code = resB$conv, message = resB$msg))
cat("  NLL above the lm() reference: "); print(summary(resB$dnll))
cat("  max |gradient| at the reported optimum: "); print(summary(resB$maxgrad))

resB2 <- do.call(rbind, lapply(1:50, run_start, control = list(eval.max = 1e4, iter.max = 1e4)))
cat("\nSame starts with eval.max = iter.max = 1e4:\n")
print(table(convergence_code = resB2$conv, message = resB2$msg))
cat("  NLL above the lm() reference: "); print(summary(resB2$dnll))

## Data-informed starts
p_good <- list(alpha = rep(mean(sox$logRS), nstocks), logbeta = rep(log(0.5), nstocks),
               logsigma = rep(log(sd(sox$logRS)), nstocks))
obj_g <- MakeADFun(negll, p_good, silent = TRUE)
fit_g <- nlminb(obj_g$par, obj_g$fn, obj_g$gr)
cat("\nData-informed starts (alpha = mean(logRS), beta = 0.5, sigma = sd(logRS)), default control:\n")
cat("  convergence", fit_g$convergence, "|", fit_g$message, "| iterations", fit_g$iterations,
    "| NLL above reference", signif(fit_g$objective - nll_ref, 3),
    "| max|grad|", signif(max(abs(obj_g$gr(fit_g$par))), 3), "\n")

## The slide's multi-start loop discards non-converged fits silently.
multi_start <- function(seed, control = list()) {
  set.seed(seed)
  obj_ms <- MakeADFun(negll, workshop_inits(), silent = TRUE)
  out <- lapply(1:50, function(i) {
    fit <- suppressWarnings(nlminb(rnorm(length(obj_ms$par), 0, 3), obj_ms$fn, obj_ms$gr,
                                   control = control))
    data.frame(conv = fit$convergence, dnll = fit$objective - nll_ref, msg = fit$message)
  })
  do.call(rbind, out)
}
cat("\nMulti-start loop from the 'Convergence' slide (50 starts ~ N(0, 3^2), default control).\n",
    "The slide keeps a run only `if(mle$convergence == 0)` and reports max(nll) - min(nll)\n",
    "over the survivors. Runs that did NOT converge, for three different seeds:\n", sep = "")
for (s in c(202, 203, 204)) {
  ms <- multi_start(s)
  cat(sprintf("  seed %d: %d of 50 dropped", s, sum(ms$conv != 0)), "|",
      paste(names(table(ms$msg[ms$conv != 0])), table(ms$msg[ms$conv != 0]), sep = ": ", collapse = "; "), "\n")
}
ms2 <- multi_start(202, control = list(eval.max = 1e4, iter.max = 1e4))
cat("  seed 202 with eval.max = iter.max = 1e4:", sum(ms2$conv != 0), "of 50 dropped;",
    "range of NLL above reference among all 50:", paste(signif(range(ms2$dnll), 3), collapse = " to "), "\n")

section("C. Profile-likelihood interval for parameter 31 (= logbeta[12])")
obj_c <- MakeADFun(negll, p_good, silent = TRUE)
fit_c <- nlminb(obj_c$par, obj_c$fn, obj_c$gr, control = list(eval.max = 1e4, iter.max = 1e4))
cat("parameter 31 is", names(obj_c$par)[31], "for stock", levels(factor(sox$PopName))[31 - nstocks], "\n")
prof <- TMB::tmbprofile(obj_c, name = 31, trace = FALSE)
ci <- confint(prof)
cat("confint(prof):", round(ci, 4), "\n")
cat("How confint.tmbprofile is defined (TMB source):\n")
print(body(getS3method("confint", "tmbprofile")))
## (i) likelihood-ratio definition: NLL rises by qchisq(0.95, 1)/2 = 1.92 above its minimum
d <- prof$value - min(prof$value)
i <- which.min(d)
lr <- suppressWarnings(
  c(approx(d[1:i], prof[[1]][1:i], qchisq(0.95, 1) / 2)$y,
    approx(d[i:nrow(prof)], prof[[1]][i:nrow(prof)], qchisq(0.95, 1) / 2)$y))
## (ii) the slide's description: central 95% of the area under the profile likelihood
o <- order(prof[[1]]); xg <- prof[[1]][o]; lik <- exp(-d[o])
cdf <- cumsum(c(0, diff(xg) * (head(lik, -1) + tail(lik, -1)) / 2)); cdf <- cdf / max(cdf)
area <- approx(cdf, xg, c(0.025, 0.975))$y
sdr_c <- sdreport(obj_c)
wald <- fit_c$par[31] + c(-1, 1) * qnorm(0.975) * sqrt(diag(sdr_c$cov.fixed))[31]
print(rbind("likelihood-ratio cut (1.92 units)" = lr,
            "95% of area under profile"         = area,
            "Wald on log scale"                 = wald))
cat("confint() reproduces the likelihood-ratio cut, not the area definition.\n")

section("D. Bootstrap with DataEval (verbatim structure from the slide)")
data <- list(stock = sox$stock,
  spawners_scaled = sox$spawners_scaled,
  transform_beta = sox$spawners/sox$sd_spawners)
boot <- list(logRS = sox$logRS)
negll_boot <- function(pars){
  getAll(pars, data, warn=FALSE)
  logRS <- DataEval(\(x){boot$logRS}) ## Updates the data.
  logRS <- OBS(logRS)
  beta <- exp(logbeta)
  sigma <- exp(logsigma)
  REPORT(beta)
  mu <- alpha[stock] - beta[stock]*spawners_scaled
  logRS %~% dnorm(mu, sigma[stock])
}
obj_b <- MakeADFun(negll_boot, p_good, silent = TRUE)
fit_b <- nlminb(obj_b$par, obj_b$fn, obj_b$gr, control = list(eval.max = 1e4, iter.max = 1e4))
set.seed(303)
y_sim1 <- obj_b$simulate(fit_b$par)$logRS
y_sim2 <- obj_b$simulate(fit_b$par)$logRS
cat("simulated data differ from the observed data:", !isTRUE(all.equal(y_sim1, sox$logRS)), "\n")

## D1. The slide's loop body (no force.update): swap the data, refit from obj$par.
boot$logRS <- y_sim1
refit_slide <- nlminb(obj_b$par, obj_b$fn, obj_b$gr, control = list(eval.max = 1e4, iter.max = 1e4))
sox_sim <- sox; sox_sim$logRS <- y_sim1
negll_sim <- function(pars){
  getAll(pars)
  beta <- exp(logbeta); sigma <- exp(logsigma)
  mu <- alpha[sox_sim$stock] - beta[sox_sim$stock]*sox_sim$spawners_scaled
  -sum(dnorm(sox_sim$logRS, mu, sigma[sox_sim$stock], log = TRUE))
}
obj_r <- MakeADFun(negll_sim, p_good, silent = TRUE)
refit_retape <- nlminb(obj_r$par, obj_r$fn, obj_r$gr, control = list(eval.max = 1e4, iter.max = 1e4))
cat("D1. slide's loop body vs a freshly taped object on the same simulated data:\n",
    "    objective difference =", signif(refit_slide$objective - refit_retape$objective, 3),
    "| max |par diff| =", signif(max(abs(refit_slide$par - refit_retape$par)), 3), "\n",
    "    (correct here only because nlminb's first evaluation is at a different point\n",
    "     from the tape's previous evaluation)\n")

## D2. Stale evaluation: same parameter vector, new data, no force.update()
nll_direct <- function(p, y) {
  a <- p[names(p) == "alpha"]; b <- exp(p[names(p) == "logbeta"]); s <- exp(p[names(p) == "logsigma"])
  -sum(dnorm(y, a[sox$stock] - b[sox$stock] * sox$spawners_scaled, s[sox$stock], log = TRUE))
}
p0 <- refit_slide$par
v_old <- obj_b$fn(p0)                 # tape evaluated at p0 with y_sim1
boot$logRS <- y_sim2                  # data change
v_stale <- obj_b$fn(p0)               # same p0 again
obj_b$force.update()
v_fresh <- obj_b$fn(p0)
cat("D2. evaluate twice at the SAME parameter vector, changing the data in between:\n")
print(c(before_change = v_old, after_change_no_update = v_stale,
        after_force.update = v_fresh, direct_R_value_new_data = nll_direct(p0, y_sim2)))
cat("=> the tape returns the OLD data's value until obj$force.update() is called\n",
    "   (the remedy documented in ?DataEval since RTMB 1.8).\n", sep = "")

section("E. An invertible Hessian is not necessarily positive definite")
H <- matrix(c(-1/3, 2/3, 2/3, -1/3), 2)
cat("H =\n"); print(H)
cat("solve(H) works and its diagonal is positive:", paste(round(diag(solve(H)), 3), collapse = ", "), "\n")
cat("but the eigenvalues of H are:", paste(round(eigen(H)$values, 3), collapse = ", "),
    " -> indefinite (a saddle point).\n")
cat("Use eigen(H)$values > 0, or sdreport(obj)$pdHess, rather than 'solve() did not fail'.\n")
