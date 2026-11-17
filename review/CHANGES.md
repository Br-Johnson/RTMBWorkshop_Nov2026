# Changes made to the workshop files

[`RTMB_workshop_review.md`](RTMB_workshop_review.md) says what is wrong or missing in the workshop material. This file lists what was then changed in the workshop sources. IDs such as M4-1 or G-1 refer to the tables in section 3 of the review, and "§4.8" refers to a section of it. Line numbers in the review refer to the files before these changes (commit `2321223`).

## Summary

- **Errors.** Every finding the review rates High is fixed, and nearly all the Medium and Low ones. The exceptions are under [Not changed](#not-changed).
- **Additions.** New slides cover the gaps the review found most costly: a fitting routine with convergence checks (Module 4); a check against `glmmTMB`, REML with a profile interval for the between-stock variance, and model comparison (Module 5); salmon reference points (Module 6); residual checks for independence and for the random effects (Module 6); and forecast evaluation with a random-walk productivity model (Module 8).
- **Builds.** All 18 documents build with Quarto 1.12.2. Before the changes three did not: `3_gradients.qmd` and `6_residuals.qmd` stopped on malformed chunk options (`#| echo:true`, `#| echo: ture`), and `8_practical.qmd` stopped in a C locale on the byte-order mark in the covariate file.
- **HTML.** The 15 HTML files in the repository were rebuilt from the new sources. They are in a separate commit so that you can drop it and rebuild on your own machine.

## How the changes were checked

- Each document was built with Quarto (R 4.3.3, RTMB 2.0, TMB 1.9.25) and its printed output read. Where new slide text states a result, the result is in the output of that slide.
- Every slide deck was opened in a browser and screenshotted to check layout, and every formula was typeset with MathJax: no errors in any document.
- Hand-coded models were checked against packages where one exists. The Module 5 model matches `glmmTMB` (negative log-likelihood 1883.81 in both). The smoother recipe in `GAM.qmd` reproduces `mgcv::gam()` to 1.5 × 10⁻⁵ when the two models are identical (verification script 07); in the document itself the two curves are close but not identical, because `gam()` has one residual SD where the workshop model has one per stock. Direct maximisation of the mixture likelihood matches `mixtools` and the document's own EM (346.47 in all three).
- The verification scripts in `verification/` still describe the original code. They quote the original slides and were not changed.

## Corrections, by file

| File | What changed | Findings |
|---|---|---|
| `modules/0_precourse_material.qmd` | "Carrying capacity" replaced by density dependence, with `S = 1/γ` and the equilibrium `α/γ` stated. `$\beta_2$` and `nlminb` typos. One sentence linking the `n/(n − a)` correction to REML. | P-1, P-4 |
| `modules/1_introduction.qmd` | RTMB described as recording a tape and replaying it, with exact derivatives. Missing bracket in the outline. | P-2 |
| `modules/2_likelihoods.qmd` | `negloglik` now returns the negative log-likelihood (three places). The benchmark loop depends on the parameter, so it times the tape and not a stored constant. A check of the RTMB value against plain R. Seed for the simulated data. Typos. | M2-1, M2-2, M2-3, P-1 |
| `modules/3_gradients.qmd` | Seeded starting values in a plausible range. `par_init` is used. `nlminb` limits raised and its message printed. Wording on `optim` and on classifying a stationary point. `∂` in the denominators. `#\| echo: true`. | M3-1 to M3-4 |
| `modules/4_inference.qmd` | `beta_real <- beta/sd_by_stock`, reported with `ADREPORT`. The multi-start loop counts failed starts and no longer drops them. Hessian check by eigenvalues. Profile interval described as a likelihood-ratio cut, parameter chosen by name. Bootstrap: `obj$force.update()` after the data change, refits start from the estimate, flagged refits are counted and kept, the real data are restored. Typos. | M4-1 to M4-6 |
| `modules/5_random_effects.qmd` | The joint-likelihood fit is labelled and its collapse shown. `ll <- dgmrf(...)`. `checkConsistency(obj, n = 2000)` with the bias in standard-error units. `obj$report(obj$env$last.par.best)`. Formula fixes (sign of β, zero-inflation indicator, Laplace σ² and limits). Wording on when Laplace is accurate. The six time-varying chunks are now evaluated. | M5-1 to M5-7 |
| `modules/6_custom_functions.qmd` | `dlognorm` tests for simulation first and no longer uses an internal RTMB function. Global on/off flags replaced by a function factory, so the simulation check compares two different models. The Jacobian constant is printed. `dlgamma` renamed `dloggamma`. The Poisson–gamma comparison is fitted and printed. The multivariate logit-normal is finished with its closed-form Jacobian and a numerical check. Built-in distributions and the status of `RTMB:::` functions mentioned. | M6-1 to M6-6 |
| `modules/6_residuals.qmd` | Pearson formula. Sign of β. `logitphi`. Duplicate `OBS()` removed. Data sorted by year before the fit. Simulated quantiles kept away from 0 and 1. `method = "fullGaussian"`. `#\| echo: true`. | M6-7 to M6-10, M5-5 |
| `modules/7_misc_coding.qmd` | `make_adfun()` uses its `pars` argument. The `ad_context()` example runs, and its interval is named a plug-in interval and compared with `lm`. Wording on uniform priors, `ADoverload("if")`, `%~%` and prediction intervals. "Option 2" label on the second tmbstan model. | M7-1 to M7-6 |
| `modules/7a_custom_distribution.qmd`, `examples/custom_distribution.qmd` | Dirichlet draws use `shape = alpha`. `log = TRUE` in the likelihood. Checks with unequal `alpha` that need no other package. The fit recovers the simulated values. The unfinished `ddirmulti_osa()` is removed and the density stops with a message if residuals are requested. `M` in the formula. The two files are still identical. | M7a-1 to M7a-4 |
| `modules/7b_tmbstan.qmd` | Only `mu[i] <- ...`. | part of M7b-1 |
| `modules/8_practical.qmd` | See [Module 8](#module-8) below. | M8-1, M8-3 to M8-8, P-1 |
| `examples/GAM.qmd` | The estimated penalty weight is kept as a first attempt that fails (λ goes to zero). The working version uses `mgcv::smooth2random()` with the spline coefficients as random effects, and is compared with `gam()`. Vectorised mean. Seeded starts. | G-1, G-2 |
| `examples/latent_variable_model.qmd` | Loading matrix corrected. Four seeded starts, summarised in a table. Penalties of fixed strength on the loadings and the residual SDs in place of the estimated one. The shared year effect is named in the model and in what the correlation plot shows. Loadings rotated before plotting. | L-1 to L-5 |
| `examples/mixture_model.qmd` | The statement that direct maximisation is not possible is corrected, and the cause of the poor first fit (starting values) is shown. Log-sum-exp with `logspace_add()`. A k-means start, from which direct maximisation, `mixtools` and the document's EM agree. A paragraph on the unbounded likelihood. Two formula typos. | X-1 to X-3 |
| `examples/osa_multinomial.qmd` | `matrix(r$residual, nrow = np)`, so the last category is really dropped. `log(p_i/p_5)`. `ncol`. Seed. | O-1, O-2 |
| `examples/ricker_example.qmd` | The restart passes log σ. Wording. | R-1, P-1 |

## New material

| Where | Slide or section | What it adds | Review |
|---|---|---|---|
| Module 4 | Putting the Checks Together | `fit_rtmb()`: raised limits, Newton steps, and four checks (optimizer code, gradient, Hessian eigenvalues, finite standard errors) | §4.1 |
| Module 5 | Check Against a Package | The same model in `glmmTMB`. Runs only if `glmmTMB` is installed | P-3, §5.7 |
| Module 5 | How Well Do We Know the Variance? | REML by adding the intercept to `random`. Profile interval for the between-stock SD, which has no lower limit | §4.3 |
| Module 5 | Is the Year Effect Worth It? | Marginal AIC, with a caution about variances at zero | §4.5 |
| Module 6 | Using It: Reference Points | S_max, S_eq, U_MSY, S_MSY and S_gen for Chilko with standard errors, using the Lambert W function the module already defines and AD `uniroot()` | §4.8(a) |
| Module 6 residuals | Which Method? | How to choose the `oneStepPredict()` method | §5.4 |
| Module 6 residuals | Are the Residuals Independent in Time? | Lag-1 correlation of the residuals within stock (0.16 overall, 0.49 for Cultus) | M6-9 |
| Module 6 residuals | Checking the Random Effects, One Draw of the Random Effects | One draw of the year effects from their distribution given the data, checked against the AR(1) model | §4.6 |
| Module 8 | The deck grew from 9 slides to 23 | See below | §4.4, §4.8 |
| `GAM.qmd` | The smoothing parameter as a variance, Compare with mgcv | Why the penalty needs its normalising term, and the check | §5.2 |

### Module 8

Module 8 changed the most, so it gets its own list. The model on the slides is the same age-structured Ricker. The reference solution around it is new.

- **Reading the data.** `fileEncoding = "UTF-8-BOM"` (M8-1).
- **Zeros.** `log((age + 1e-4)/spawners)` is replaced by treating a zero return as missing. For Stellako this changes the age-5 SD from 1.18 to 1.02 and the age-4/age-5 correlation from 0.33 to 0.49, because the single zero no longer dominates. The likelihood now handles a missing age-4 value. This is a stop-gap for stocks where zeros are rare. The slides say so and point to modelling total recruits and age composition separately, which is not implemented (M8-5).
- **Median and mean.** `Eage4` is renamed `med4`, and the forecast is reported both ways: 0.359 against 0.558 million (M8-6).
- **Forecast interval.** `forecast_dist()` simulates estimation uncertainty and process variation (§4.4).
- **Hindcast.** One `hindcast()` function for all models. It reports at the best parameters, raises the optimizer limits and records convergence (M5-5, M8-7).
- **Evaluation.** Bias, mean absolute error, and error relative to the naive forecast "the return four years earlier", plus coverage of the 80% interval (M8-7, §4.8(e)).
- **Covariates.** No intercept column, covariates standardised with the data available at each forecast, one covariate set throughout, z-values and AIC (M8-3, M8-4).
- **Random-walk productivity.** One more random effect and two lines of likelihood (§4.8(d)).

What the new slides report for Stellako, 22 forecasts (millions of fish, mean return 0.27):

| Forecast | Bias | Mean absolute error | Too high in | Error relative to naive | Inside 80% interval |
|---|---|---|---|---|---|
| Stationary Ricker | +0.142 | 0.238 | 21 of 22 | 1.13 | 14 of 22 |
| Ricker with 5 covariates | +0.106 | 0.216 | 17 of 22 | 1.02 | 14 of 22 |
| Random-walk productivity | −0.012 | 0.174 | 12 of 22 | 0.82 | 16 of 22 |
| Naive | +0.051 | 0.211 | 11 of 22 | 1.00 | |

These differ a little from the table in §4.8(e) of the review, which kept the original `+1e-4`.

## Not changed

- **`modules/7b_tmbstan.qmd`** is still a copy of the pre-course notes under the title "Custom Distribution" (M7b-1). A Bayesian module has to be written (§4.9), and `tmbstan` could not be installed where this work was done, so nothing here could have been run. Module 7's tmbstan slide has two new bullets and a corrected label; its code is unchanged and still `eval: false`.
- **The data files.** The disagreement between `recruits` and the age columns for brood year 2018 needs whoever maintains the file (M8-2). File names with spaces and the byte-order mark are as they were; the reading code now copes with the mark.
- **`examples/custom_distribution.qmd`** duplicates Module 7a. It was kept, and kept identical, in case something links to it.
- **Modelling zeros properly** (M8-5), the **Larkin model** (§4.8(c)), **bias correction** with `sdreport(bias.correct = TRUE)` (§4.4) and **conditional AIC or leave-future-out scoring** beyond the hindcast (§4.5) are named on the slides as things to try but are not worked examples.
- **`README.md`, `homework.qmd`, `examples/ricker_lm.R`** are untouched.

## The rebuilt HTML

- Built on Linux with Quarto 1.12.2, R 4.3.3, RTMB 2.0, TMB 1.9.25, ggplot2 3.4.4 and glmmTMB 1.1.8. The previous files were built with Quarto 1.3.353, so the HTML diff is large for reasons that have nothing to do with content.
- `RTMBdist` was not installed on the build machine (it needs CRAN, which that machine could not reach). In Module 7a the comparison with `RTMBdist::ddirmult()` sits behind `requireNamespace()`, so the committed HTML shows the call without output. It runs where `RTMBdist` is installed.
- The six pages with `format: html` (pre-course notes, Module 7a and the four examples) link to one script on `cdnjs.cloudflare.com` and do not embed it, because the build machine could not reach that host. It is a compatibility script for old browsers that Pandoc adds next to MathJax. The pages display normally. MathJax itself was loaded from the web in the old files and still is, so math needs a connection, as before.
- Timings on the benchmark slides are from the build machine.
- The build machine's `glmmTMB` (1.1.8, from Ubuntu) was compiled against TMB 1.9.10 and warns about that when loaded next to TMB 1.9.25. The warning is in the build log only, and the fit still matches RTMB.
- The same set of 15 HTML files is committed as before. `7b_tmbstan`, `ricker_example` and `examples/custom_distribution` have no HTML in the repository and were not added.

To rebuild everything on your own machine, from the project root:

```sh
quarto render modules/8_practical.qmd      # one document
for f in modules/*.qmd examples/*.qmd; do quarto render "$f"; done
```

Packages the documents load: RTMB, TMB, dplyr, tidyr, ggplot2, gganimate, here, MASS, Matrix, mgcv, microbenchmark, pracma, numDeriv, igraph, splines2, mixtools, corrplot, gclus, optimx, and optionally glmmTMB, RTMBdist and tmbstan. The changes add no required package: `MASS` and `mgcv` come with R, and `glmmTMB` is used only if it is installed.
