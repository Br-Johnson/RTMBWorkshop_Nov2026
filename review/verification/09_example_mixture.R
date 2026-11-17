## examples/mixture_model.qmd: four-component Gaussian mixture.
## The document says direct maximisation "is not possible due to the sum within a log" and that
## EM is "way better for the same starting values".
##  A. What the document's starting values are, relative to the data.
##  B. Direct maximisation vs EM from the SAME sensible starting values, over 50 simulated data sets.
##  C. The heteroscedastic Gaussian mixture likelihood is unbounded (true for EM as well).
##  D. Log-sum-exp: why `log(sum(prob * dnorm(...)))` is fragile.

source(here::here("review", "verification", "00_setup.R"))
suppressPackageStartupMessages(library(mixtools))
big <- list(eval.max = 1e4, iter.max = 1e4)

K <- 4
mu_true <- c(35,50,59,72); sigma_true <- c(1.5, 2, 2.5, 7); prop_true <- c(0.1, 0.45, 0.3, 0.15)
simulate_x <- function(nobs = 100) {
  z <- sample(x = 1:K, size = nobs, replace = TRUE, prob = prop_true)
  rnorm(nobs, mu_true[z], sigma_true[z])
}
itrans <- function(lp){ p <- c(exp(lp), 1); p/sum(p) }
## mixtools::normalmixEM prints its iteration count; run it quietly and return a try-error on failure.
quiet_em <- function(...) {
  em <- NULL
  invisible(capture.output(em <- try(suppressWarnings(normalmixEM(..., verb = FALSE)), silent = TRUE)))
  em
}

## Direct marginal likelihood, written with logspace_add() so nothing under- or overflows.
make_nll <- function(x) function(pars){
  getAll(pars)
  lprob <- log(itrans(logitp))
  sigma <- exp(logsigma)
  acc <- lprob[1] + dnorm(x, mu[1], sigma[1], log = TRUE)
  for (k in 2:K) acc <- logspace_add(acc, lprob[k] + dnorm(x, mu[k], sigma[k], log = TRUE))
  -sum(acc)
}
## The document's version of the same likelihood
make_nll_doc <- function(x) function(pars){
  getAll(pars)
  prob = itrans(logitp)
  sigma <- exp(logsigma)
  ll <- 0
  for( i in 1:length(x) ){
    pobs <- sum(prob*dnorm(x[i], mu, sigma))  ## Numerical Issues?
    ll <- ll + log(pobs)
  }
  -ll
}

section("A. The document's starting values")
set.seed(1)
x <- simulate_x()
cat("data range:", paste(round(range(x), 1), collapse = " to "), "\n")
set.seed(2)
start_doc <- replicate(5, sort(round(exp(rnorm(K, 3, 0.25)), 1)))
cat("five draws of the document's starting means, exp(rnorm(4, 3, 0.25)), sorted:\n"); print(t(start_doc))
cat("starting SD for every component: exp(2) =", round(exp(2), 2), "\n")
cat("=> all four components start below the smallest observation, 10-55 units from most of the data.\n")
res_doc <- t(sapply(1:20, function(seed) {
  set.seed(seed)
  param <- list(logitp = c(0,0,0), mu = exp(rnorm(K, 3, 0.25)), logsigma = numeric(K)+2)
  o <- MakeADFun(make_nll_doc(x), param, silent = TRUE)
  f <- suppressWarnings(nlminb(o$par, o$fn, o$gr))
  em <- quiet_em(x, lambda = itrans(param$logitp), mu = param$mu, sigma = exp(param$logsigma), maxit = 5000)
  c(direct_nll = f$objective, em_nll = if (inherits(em, "try-error")) NA else -em$loglik)
}))
cat("From the document's starting values, 20 random draws (one data set):\n")
cat("  direct NLL: min", round(min(res_doc[, 1]), 2), "| median", round(median(res_doc[, 1]), 2), "| max", round(max(res_doc[, 1]), 2), "\n")
cat("  EM NLL:     min", round(min(res_doc[, 2], na.rm = TRUE), 2), "| median", round(median(res_doc[, 2], na.rm = TRUE), 2),
    "| max", round(max(res_doc[, 2], na.rm = TRUE), 2), "| EM failures:", sum(is.na(res_doc[, 2])), "\n")
cat("  NLL of a single normal fitted to all the data (what 'one component took everything' means):",
    round(-sum(dnorm(x, mean(x), sqrt(mean((x - mean(x))^2)), log = TRUE)), 2), "\n")

section("B. Same sensible starting values for both algorithms, 50 simulated data sets")
sensible_start <- function(x) {
  km <- kmeans(x, centers = quantile(x, (1:K - 0.5) / K))          # deterministic k-means start
  o <- order(km$centers)
  list(logitp = log(km$size[o][1:3] / km$size[o][4]), mu = as.numeric(km$centers[o]),
       logsigma = log(pmax(tapply(x, match(km$cluster, o), sd), 0.5)))
}
cmp <- t(sapply(1:50, function(seed) {
  set.seed(seed); x <- simulate_x()
  st <- sensible_start(x)
  o <- MakeADFun(make_nll(x), st, silent = TRUE)
  t_d <- system.time(f <- suppressWarnings(nlminb(o$par, o$fn, o$gr, control = big)))["elapsed"]
  t_e <- system.time(em <- quiet_em(x, lambda = itrans(st$logitp), mu = st$mu, sigma = exp(st$logsigma),
                                    maxit = 10000, epsilon = 1e-10))["elapsed"]
  truth <- -sum(log(sapply(x, function(xi) sum(prop_true * dnorm(xi, mu_true, sigma_true)))))
  c(direct = f$objective, em = if (inherits(em, "try-error")) NA else -em$loglik, truth = truth,
    min_sigma_direct = min(exp(f$par[names(f$par) == "logsigma"])),
    min_sigma_em = if (inherits(em, "try-error")) NA else min(em$sigma),
    em_iter = if (inherits(em, "try-error")) NA else length(em$all.loglik) - 1,
    direct_iter = f$iterations, t_direct = t_d, t_em = t_e)
}))
ok <- !is.na(cmp[, "em"])
d <- cmp[ok, "direct"] - cmp[ok, "em"]
cat("EM failed outright (mixtools error) in", sum(!ok), "of 50 data sets\n")
cat("direct NLL minus EM NLL over the rest:  equal within 0.01:", sum(abs(d) <= 0.01),
    "| direct lower (better):", sum(d < -0.01), "| EM lower (better):", sum(d > 0.01), "\n")
cat("both end at or below the NLL of the TRUE parameters:  direct", sum(cmp[, "direct"] <= cmp[, "truth"]),
    "of 50 | EM", sum(cmp[ok, "em"] <= cmp[ok, "truth"]), "of", sum(ok), "\n")
cat("median iterations: direct (quasi-Newton with AD gradient)", median(cmp[, "direct_iter"]),
    "| EM", median(cmp[ok, "em_iter"]), "\n")
cat("fits ending with a component SD below 0.3 (true SDs are 1.5-7): direct", sum(cmp[, "min_sigma_direct"] < 0.3),
    "| EM", sum(cmp[ok, "min_sigma_em"] < 0.3), "\n")

section("C. The likelihood has no maximum: a component can collapse onto one observation")
set.seed(1); x <- simulate_x()
st <- sensible_start(x)
o <- MakeADFun(make_nll(x), st, silent = TRUE)
f <- suppressWarnings(nlminb(o$par, o$fn, o$gr, control = big))
p_spike <- o$env$parList(f$par)
p_spike$mu[1] <- x[1]
for (s in c(1, 1e-4, 1e-8, 1e-30, 1e-100)) {
  p_spike$logsigma[1] <- log(s)
  cat(sprintf("  component 1 centred on x[1] with SD %-6g  -> NLL %.2f\n", s, o$fn(unlist(p_spike))))
}
cat("  NLL at the fitted (interior) solution:", round(f$objective, 2), "\n")
cat("=> the NLL falls without limit as that SD shrinks (about 2.3 units per factor of 10), so the\n",
    "   'maximum likelihood estimate' does not exist; optimisers (EM included) find local maxima.\n",
    "   The rendered document's own EM answer has a component with SD exp(-1.52) = 0.22.\n", sep = "")

section("D. Numerical stability of log(sum(prob * dnorm()))")
x_far <- 500                                  # an observation far from every component
lp <- log(prop_true) + dnorm(x_far, mu_true, sigma_true, log = TRUE)
cat("naive  log(sum(prob * dnorm(x))) =", log(sum(prop_true * dnorm(x_far, mu_true, sigma_true))), "\n")
cat("stable log-sum-exp of the log terms =", max(lp) + log(sum(exp(lp - max(lp)))), "\n")
