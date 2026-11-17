## The README promises "Compare output with glmmTMB" for the hierarchical model; the modules do not
## contain that comparison. Here it is, for the Module 5 model:
##   log(R/S) = alpha0 + re[stock] - beta[stock] * S_scaled + eps,   eps ~ N(0, sigma[stock]^2)
## This script needs only glmmTMB + dplyr (it does not load RTMB), so that it can be run with
## whatever TMB version your glmmTMB was built against.

suppressPackageStartupMessages({library(glmmTMB); library(dplyr)})
cat(sprintf("[setup] R %s | glmmTMB %s | TMB %s\n", getRversion(), packageVersion("glmmTMB"), packageVersion("TMB")))

sox <- read.csv(here::here("data", "sockeye", "SRDATA2026 1.csv")) %>%
  left_join(read.csv(here::here("data", "sockeye", "FC_StockInfo 2.csv")), by = "PopID") %>%
  filter(!is.na(spawners), !is.na(recruits)) %>%
  group_by(PopName) %>% mutate(spawners_scaled = spawners / sd(spawners)) %>% ungroup() %>%
  mutate(logRS = log(recruits / spawners), PopName = factor(PopName)) %>% as.data.frame()

fit <- function(REML) glmmTMB(logRS ~ 1 + PopName:spawners_scaled + (1 | PopName),
                              dispformula = ~ 0 + PopName, data = sox, REML = REML)
m_ml   <- fit(FALSE)
m_reml <- fit(TRUE)
sd_of <- function(m) attr(VarCorr(m)$cond$PopName, "stddev")[[1]]

cat("\nglmmTMB, maximum likelihood:\n")
cat("  -logLik =", round(-as.numeric(logLik(m_ml)), 4), "| between-stock SD =", round(sd_of(m_ml), 4),
    "| intercept (alpha0) =", round(fixef(m_ml)$cond[["(Intercept)"]], 4), "\n")
cat("  RTMB values from 03_module5_random_effects.R / 13_recommended_workflow.R: -logLik = 1883.8099, SD = 0.1056\n")
cat("  all 19 fitted slopes negative (so the exp() constraint in the slides is not binding):",
    all(fixef(m_ml)$cond[-1] < 0), "\n")
cat("\nglmmTMB, REML (all fixed effects integrated out):\n")
cat("  between-stock SD =", round(sd_of(m_reml), 4), " (ML:", round(sd_of(m_ml), 4), "-> ratio",
    round(sd_of(m_reml) / sd_of(m_ml), 2), ")\n")
ci <- suppressWarnings(try(confint(m_ml, parm = "theta_", method = "profile"), silent = TRUE))
if (!inherits(ci, "try-error")) {
  ## glmmTMB profiles the log of the SD ("theta"); NA for the lower limit means the profile never
  ## rises by 1.92 units on that side, i.e. the interval extends to SD = 0.
  lo <- ci[1, 1]; hi <- ci[1, 2]
  cat("\nprofile-likelihood 95% interval for the between-stock SD (ML fit):",
      if (is.na(lo)) "0 (no lower limit found)" else signif(exp(lo), 3), "to", signif(exp(hi), 3), "\n")
} else cat("\n(profile interval not available in this glmmTMB version)\n")
