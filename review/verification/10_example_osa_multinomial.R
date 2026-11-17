## examples/osa_multinomial.qmd: OSA residuals for multinomial observations.
##  A. `as.matrix(r$residual, nrow = np, ncol = 100)` does not reshape; `osa[-np, ]` then drops
##     one value instead of the 100 uninformative last-category residuals.
##  B. What those last-category residuals are, and what the QQ plot should be based on.
##  C. Would the diagnostic detect a wrong model? (over-dispersed data, as in real age compositions)

source(here::here("review", "verification", "00_setup.R"))

mlogit  <- function(p){ np <- length(p); log(p[1:(np-1)]/p[np]) }
imlogit <- function(alpha){ expalpha <- exp(c(alpha, 0)); expalpha/sum(expalpha) }
np <- 5; size <- 15; p_true <- c(0.1,0.5,0.2,0.05,0.15)

fit_and_osa <- function(X) {
  dat <- list(X = X)
  nll <- function (par){                    # as in the document (Pearson residual code omitted)
    getAll(par, dat)
    X <- OBS(X)
    p <- imlogit(alpha)
    ll <- 0
    for(i in 1:100) ll <- ll + dmultinom(X[,i], prob=p, log=TRUE)
    -ll
  }
  obj <- MakeADFun(nll, list(alpha = numeric(np-1)), silent = TRUE)
  fit <- nlminb(obj$par, obj$fn, obj$gr)
  r <- oneStepPredict(obj, method = "oneStepGeneric", observation.name = "X",
                      discrete = TRUE, discreteSupport = 0:size, trace = FALSE)
  list(fit = fit, r = r)
}

section("A. Shape of the residual object")
set.seed(1)
X <- rmultinom(100, size = size, prob = p_true)
out <- fit_and_osa(X)
r <- out$r
osa_doc <- as.matrix(r$residual, nrow = np, ncol = 100)       # verbatim from the document
cat("dim(as.matrix(r$residual, nrow = np, ncol = 100)):", paste(dim(osa_doc), collapse = " x "),
    " (as.matrix() has no nrow/ncol arguments; they are silently ignored)\n")
cat("number of residuals left by osa[-np, ]:", length(osa_doc[-np, ]), " (only element", np, "is dropped)\n")
osa_ok <- matrix(r$residual, nrow = np, ncol = 100)           # categories in rows, samples in columns
cat("dim(matrix(r$residual, nrow = np, ncol = 100)):", paste(dim(osa_ok), collapse = " x "),
    "| residuals left after dropping the last category:", length(osa_ok[-np, ]), "\n")

section("B. The last-category residuals")
last <- osa_ok[np, ]
cat("Given the first four counts and the total, the fifth count is known exactly, so its\n",
    "'residual' carries no information about fit (WHAM, for example, removes these).\n",
    "Values returned for the 100 last-category cells:\n", sep = "")
print(summary(last))
cat("number of non-finite values among them:", sum(!is.finite(last)), "\n")
doc_set <- osa_doc[-np, ]; ok_set <- as.numeric(osa_ok[-np, ])
cat(sprintf("document's set: n = %d, non-finite = %d | corrected set: n = %d, non-finite = %d\n",
            length(doc_set), sum(!is.finite(doc_set)), length(ok_set), sum(!is.finite(ok_set))))
cat(sprintf("corrected set: mean %.3f, sd %.3f, Shapiro-Wilk p = %.3f (the model is true here)\n",
            mean(ok_set), sd(ok_set), shapiro.test(ok_set)$p.value))

section("C. Power to detect over-dispersion (Dirichlet-multinomial data fitted with a multinomial)")
rdirmult <- function(n, size, alpha) sapply(1:n, function(i) {
  th <- rgamma(length(alpha), shape = alpha, rate = 1); rmultinom(1, size, th / sum(th)) })
sim_once <- function(seed) {
  set.seed(seed)
  Xo <- rdirmult(100, size, 10 * p_true)                 # effective sample size well below 15
  ro <- fit_and_osa(Xo)$r$residual
  m_ok  <- matrix(ro, nrow = np)[-np, ]
  m_doc <- as.matrix(ro)[-np, ]
  fin <- function(v) v[is.finite(v)]
  c(sd_corrected = sd(fin(m_ok)), sd_document = sd(fin(m_doc)),
    p_corrected = shapiro.test(fin(m_ok))$p.value, p_document = shapiro.test(fin(m_doc))$p.value)
}
pw <- t(sapply(1:20, sim_once))
cat("20 over-dispersed data sets; residual SD should be clearly above 1 if the check has power:\n")
print(round(colMeans(pw[, 1:2]), 3))
cat("share of the 20 data sets in which Shapiro-Wilk rejects at 5%: corrected",
    mean(pw[, "p_corrected"] < 0.05), "| document", mean(pw[, "p_document"] < 0.05), "\n")
cat("=> here the reshaping bug changes the diagnostic only slightly, because most last-category\n",
    "   residuals are returned as NaN and qqnorm() drops NaN silently; a handful of spurious\n",
    "   finite values remain in the document's plot. The code is still wrong and would matter\n",
    "   for any summary that does not drop NaN (mean(), sd(), acf()).\n", sep = "")
