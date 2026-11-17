## examples/latent_variable_model.qmd: two latent factors by year with stock-specific loadings.
##  A. Fit as written; which residual SD collapses (the "Heywood case" the text notices).
##  B. Post-processing bug: second loading column is multiplied by the first.
##  C. The shared year effect is left out of the "covariance between stocks".
##  D. The "penalised" refit: the penalty weight is estimated, so it goes to zero.
##  E. Loadings depend on which stocks happen to be listed first; rotate before interpreting.

source(here::here("review", "verification", "00_setup.R"))
sox <- load_sox()
nstocks <- max(sox$stock); nyears <- max(sox$time); nobs <- nrow(sox)
stock_names <- levels(factor(sox$PopName))
big <- list(eval.max = 1e4, iter.max = 1e4)

## --- verbatim objective from the document (first model) -------------------------
negll <- function(pars){

  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta)

  logRS <- OBS(logRS)

  sigma <- exp(logsigma)
  sigma_re <- exp(logsigma_re)
  sigma_time <- exp(logsigma_time)
  sigma_obs <- sigma[stock]
  scale_z <- exp(logscale_z)
  phi <- 2*plogis(logitphi)-1

  alpha_stock <- re + alpha0
  ADREPORT(alpha_stock)

  mu <- AD(numeric(nobs))
  for( i in 1:nobs ){
    mu[i] <- alpha_stock[stock[i]] + re_time[time[i]] - beta[stock[i]]*spawners_scaled[i] +
      z[time[i],1]*lambda1[stock[i]]*scale_z[1] + z[time[i],2]*lambda2[stock[i]]*scale_z[2]
  }

  ## Use this for built in residuals
  ll <- dautoreg(re_time,phi=phi,scale=sqrt(sigma_time^2/(1-phi^2)), log=TRUE)
  ## Add in latent variables
  for( i in 1:2 ) ll <- ll + sum(dnorm(z[,i], 0, 1, log = TRUE))


  ll <- ll + sum(dnorm(re, 0, sigma_re, log = TRUE))
  ll <- ll + sum(dnorm(logRS, mu, sigma_obs, log = TRUE))

  REPORT(mu)
  REPORT(re_time)
  REPORT(sigma_obs)

  return(-ll)
}
## ------------------------------------------------------------------------------
workshop_inits <- function() list(
  alpha0 = runif(1, 2, 5),
  re = rnorm(nstocks),
  logbeta = runif(nstocks, 0, 1),
  logsigma = runif(nstocks, 3, 5),
  logsigma_re = runif(1, 0.5, 2),
  re_time = rnorm(nyears),
  logitphi = 0,
  logsigma_time = 0,
  z = matrix(0, nrow = nyears, ncol = 2),
  logscale_z = c(0,0),
  lambda1 = c(1, rnorm(nstocks-1)),
  lambda2 = c(0,1, rnorm(nstocks-2))
)
map_lam <- list(lambda1 = factor(c(NA, 2:nstocks)), lambda2 = factor(c(NA, NA, 3:nstocks)))

section("A. Fit as written (default nlminb control, the document's random starting values)")
fits <- lapply(1:3, function(seed) {
  set.seed(seed)
  obj <- MakeADFun(func = negll, parameters = workshop_inits(), random = c("re", "re_time", "z"),
                   map = map_lam, silent = TRUE)
  fit <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr))
  list(obj = obj, fit = fit)
})
for (i in 1:3) with(fits[[i]], cat(sprintf(
  "seed %d: NLL %.3f | convergence %d (%s) | iterations %d | smallest residual SD %.2g (%s)\n",
  i, fit$objective, fit$convergence, fit$message, fit$iterations,
  min(exp(fit$par[names(fit$par) == "logsigma"])),
  stock_names[which.min(fit$par[names(fit$par) == "logsigma"])])))
## Continue the best of them with a larger iteration budget
best <- which.min(sapply(fits, function(f) f$fit$objective))
obj <- fits[[best]]$obj
fit <- suppressWarnings(nlminb(fits[[best]]$fit$par, obj$fn, obj$gr, control = big))
cat(sprintf("continued (seed %d) with eval.max = iter.max = 1e4: NLL %.3f | convergence %d (%s)\n",
            best, fit$objective, fit$convergence, fit$message))
sdr <- suppressWarnings(sdreport(obj))
cat("pdHess:", sdr$pdHess, "| max |gradient|:", signif(max(abs(sdr$gradient.fixed)), 3),
    "| number of NaN standard errors:", sum(is.nan(suppressWarnings(sqrt(diag(sdr$cov.fixed))))), "of", length(fit$par), "\n")
pList <- obj$env$parList(fit$par, obj$env$last.par.best)
sig <- exp(pList$logsigma); names(sig) <- stock_names
cat("three smallest residual SDs:\n"); print(signif(sort(sig)[1:3], 3))

section("B. Factor loadings in the post-processing chunk")
cat("`pList$logscale` (the name used in the document) works only through partial matching of `$`:",
    !is.null(pList$logscale), "| exact name is `logscale_z`\n")
Lambda_doc <- cbind(pList$lambda1*exp(pList$logscale[1]),
                    pList$lambda2*pList$lambda1*exp(pList$logscale[2]))     # verbatim
Lambda_ok  <- cbind(pList$lambda1*exp(pList$logscale_z[1]),
                    pList$lambda2*exp(pList$logscale_z[2]))                 # what the model uses
cat("second column, document vs model (first 6 stocks):\n")
print(round(rbind(document = Lambda_doc[1:6, 2], model = Lambda_ok[1:6, 2]), 3))
C_doc <- cov2cor(Lambda_doc %*% t(Lambda_doc) + diag(sig^2))
C_ok  <- cov2cor(Lambda_ok  %*% t(Lambda_ok)  + diag(sig^2))
off <- upper.tri(C_ok)
cat("pairwise stock correlations: largest |document - model| =", round(max(abs(C_doc[off] - C_ok[off])), 3),
    "| pairs whose sign differs:", sum(sign(C_doc[off]) != sign(C_ok[off])), "of", sum(off), "\n")

section("C. What the 'covariance between stocks' leaves out")
phi <- 2*plogis(pList$logitphi) - 1
var_year <- exp(pList$logsigma_time)^2 / (1 - phi^2)       # stationary variance of the shared AR(1) effect
C_full <- cov2cor(Lambda_ok %*% t(Lambda_ok) + var_year + diag(sig^2))
cat(sprintf("shared year effect: stationary SD %.3f (phi = %.2f). It adds the same covariance %.3f to every pair of stocks.\n",
            sqrt(var_year), phi, var_year))
cat("mean pairwise correlation: factors only (as in the document, corrected) =", round(mean(C_ok[off]), 3),
    "| factors + shared year effect =", round(mean(C_full[off]), 3), "\n")

section("D. The 'penalised' refit with logpen estimated")
negll_pen <- negll
body(negll_pen)[[length(body(negll_pen))]] <- quote(return(-ll + exp(logpen)/2*sum(lambda1^2 + lambda2^2)))
start_pen <- c(pList, list(logpen = 0))          # start from the fitted values, penalty weight exp(0) = 1
obj_p <- MakeADFun(func = negll_pen, parameters = start_pen, random = c("re", "re_time", "z"),
                   map = map_lam, silent = TRUE)
cat("objective at the unpenalised optimum with penalty weight 1:", round(obj_p$fn(), 3), "\n")
fit_p <- suppressWarnings(nlminb(obj_p$par, obj_p$fn, obj_p$gr, control = big))
cat(sprintf("after optimisation: logpen = %.1f (weight %.2g) | objective %.3f (unpenalised NLL %.3f) | convergence %d (%s)\n",
            fit_p$par["logpen"], exp(fit_p$par["logpen"]), fit_p$objective, fit$objective, fit_p$convergence, fit_p$message))
cat("largest change in any loading relative to the unpenalised fit:",
    signif(max(abs(fit_p$par[grepl("lambda", names(fit_p$par))] - fit$par[grepl("lambda", names(fit$par))])), 3), "\n")
cat("=> as in the GAM example, a penalty weight that is itself minimised over goes to zero.\n")

## A penalty that does something: fix the strength (here a N(0, 1) prior on the free loadings)
## and keep the residual SDs away from zero with a weak boundary-avoiding prior.
negll_fix <- negll
body(negll_fix)[[length(body(negll_fix))]] <- quote(return(
  -ll - sum(dnorm(lambda1, 0, 1, log = TRUE)) - sum(dnorm(lambda2, 0, 1, log = TRUE)) -
    sum(dgamma(sigma, shape = 2, rate = 2, log = TRUE))))
obj_f <- MakeADFun(func = negll_fix, parameters = pList, random = c("re", "re_time", "z"), map = map_lam, silent = TRUE)
fit_f <- suppressWarnings(nlminb(obj_f$par, obj_f$fn, obj_f$gr, control = big))
sdr_f <- suppressWarnings(sdreport(obj_f))
pF <- obj_f$env$parList(fit_f$par, obj_f$env$last.par.best)
cat(sprintf("fixed-strength penalties: convergence %d | pdHess %s | smallest residual SD %.3f (%s) | NaN SEs: %d\n",
            fit_f$convergence, sdr_f$pdHess, min(exp(pF$logsigma)), stock_names[which.min(pF$logsigma)],
            sum(is.nan(suppressWarnings(sqrt(diag(sdr_f$cov.fixed)))))))

section("E. Loadings are identified only up to rotation")
L <- cbind(pF$lambda1*exp(pF$logscale_z[1]), pF$lambda2*exp(pF$logscale_z[2])); rownames(L) <- stock_names
cat("The constraint lambda2[1] = 0 ties the axes to whichever stock is first alphabetically (",
    stock_names[1], ").\n", sep = "")
vm <- varimax(L)
cat("share of loading variance on factor 1 / factor 2: as estimated",
    paste(round(colSums(L^2) / sum(L^2), 2), collapse = " / "),
    "| after varimax rotation", paste(round(colSums(vm$loadings^2) / sum(L^2), 2), collapse = " / "), "\n")
cat("L %*% t(L) is unchanged by rotation: max |difference| =",
    signif(max(abs(L %*% t(L) - vm$loadings %*% t(vm$loadings))), 3), "\n")
