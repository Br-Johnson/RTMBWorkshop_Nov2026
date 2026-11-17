## Shared setup for the verification scripts in this folder.
##
## Run every script from the project root, e.g.
##   Rscript review/verification/01_module2_likelihood.R
##
## The data preparation below mirrors the hidden setup chunk that modules 3-7
## all share, so the objects have the same names and meaning as in the slides.

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
})

cat(sprintf("[setup] R %s | RTMB %s | TMB %s | locale %s\n",
            getRversion(), packageVersion("RTMB"), packageVersion("TMB"),
            Sys.getlocale("LC_CTYPE")))

load_sox <- function() {
  sox   <- read.csv(here::here("data", "sockeye", "SRDATA2026 1.csv"))
  stock <- read.csv(here::here("data", "sockeye", "FC_StockInfo 2.csv"))
  sox <- sox %>%
    left_join(stock, by = "PopID") %>%
    filter(!is.na(spawners), !is.na(recruits)) %>%
    group_by(PopName) %>%
    mutate(sd_spawners = sd(spawners), spawners_scaled = spawners / sd_spawners) %>%
    ungroup() %>%
    as.data.frame()
  sox$stock <- as.numeric(factor(sox$PopName))
  sox$logRS <- log(sox$recruits / sox$spawners)
  sox$time  <- sox$yr - min(sox$yr) + 1
  sox
}

## The environmental covariate file starts with a UTF-8 byte-order mark; reading
## it this way works in every locale (see 06_module8_forecast.R).
load_covariates <- function() {
  read.csv(here::here("data", "sockeye", "FC_Environmental_Data_2026_Std_Offset_Real.csv"),
           fileEncoding = "UTF-8-BOM")
}

section <- function(...) cat("\n", strrep("=", 78), "\n", ..., "\n", strrep("=", 78), "\n", sep = "")
