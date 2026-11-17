## Run every verification script and save its console output next to it.
##   Rscript review/verification/run_all.R
## Needs: RTMB (>= 1.9), TMB, dplyr, here, MASS, Matrix, mgcv, mixtools, and glmmTMB (script 14 only).
## Total run time is roughly 20 minutes; scripts 03 and 04 account for most of it.

dir_v   <- here::here("review", "verification")
dir_out <- file.path(dir_v, "output")
dir.create(dir_out, showWarnings = FALSE)

scripts <- sort(list.files(dir_v, pattern = "^(0[1-9]|1[0-9])_.*\\.R$", full.names = TRUE))
rscript <- file.path(R.home("bin"), "Rscript")

for (s in scripts) {
  out <- file.path(dir_out, sub("\\.R$", ".txt", basename(s)))
  t0 <- proc.time()["elapsed"]
  status <- system2(rscript, shQuote(s), stdout = out, stderr = out)
  cat(sprintf("%-45s exit %d  %5.0f s\n", basename(s), status, proc.time()["elapsed"] - t0))
}
