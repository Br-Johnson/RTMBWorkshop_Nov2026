# Review of the PSC RTMB workshop materials (November 2026)

Prepared 17 November 2026 with Claude Code (an AI coding assistant) at the repository owner's request, for the fork `Br-Johnson/RTMBWorkshop_Nov2026`. It covers every `.qmd` and `.R` file in `modules/` and `examples/` at commit `2321223`. `homework.qmd` is your own work and is not reviewed here. "I" below is the assistant.

The review asks three questions of the material:

1. Is anything missing?
2. Is anything in it that should not be, because it is wrong or poor practice?
3. Are more modern or more sophisticated approaches available that it does not mention?

**How to read this.** Section 1 answers the three questions in about two pages. Sections 3 to 5 give the detail. Section 6 says what I would do with it. Every claim carries a tag saying how it was checked.

| Tag | Meaning |
|---|---|
| **[run]** | I executed the code. The script and its saved output are in `review/verification/`. "V06-F" means script `06_…R`, section F. |
| **[doc]** | Stated in current package documentation or source that I read directly (RTMB 2.0, TMB 1.9.25, and others named in the text). |
| **[read]** | Visible from reading the workshop source. No execution needed. |
| **[lit]** | Published literature. Network limits meant most papers were confirmed from abstracts or search snippets, not full text. Section 7 says what could not be confirmed. |

Contents: [1 Short answers](#1-short-answers) · [2 What the workshop does well](#2-what-the-workshop-does-well) · [3 Errors and poor practice](#3-errors-and-poor-practice-by-module) · [4 What is missing](#4-what-is-missing) · [5 More modern approaches](#5-more-modern-or-sophisticated-approaches) · [6 What to do with this](#6-what-to-do-with-this) · [7 Limits of this review](#7-limits-of-this-review) · [8 References](#8-references) · [Appendix: verification scripts](#appendix-verification-scripts)

---

## 1. Short answers

### Overall

The workshop is a sound introduction to the mechanics of RTMB: write a likelihood in R, get exact gradients, declare random effects, report quantities with standard errors. Its central hierarchical model is correct. An independent fit in glmmTMB reproduces its likelihood (−logLik 1883.81) and its between-stock standard deviation (0.106) **[run V14]**.

The problems fall into three groups.

- **Code on the slides that is wrong, crashes, or works only by accident.** I confirmed more than twenty cases by running them. Many sit in chunks marked `eval: false`, so they were never executed when the slides were built. One is hidden by a second bug, and another by a test case that happens not to exercise it.
- **Three worked examples reach conclusions that do not hold.** The GAM example's smoothing penalty is estimated in a way that removes it. The latent-variable example did not converge in any of my replications, and its correlation plot is computed from the wrong loadings. The mixture example says direct maximisation "is not possible" when it works as well as EM.
- **The applied end of the course stops short.** Module 8 produces a point forecast with no uncertainty, no skill score and no benchmark. Run with those added, the workshop's own hindcast shows the forecast is biased high and, by mean absolute error, no better than "same as four years ago". The two things these data most clearly call for, time-varying productivity and delayed density dependence, are not used.

Use the workshop to learn how RTMB works. Do not carry Module 8 or the three examples into real analysis without the fixes in section 3.

### Is anything missing?

Yes. In rough order of importance for stock–recruit and forecasting work:

- **Forecast evaluation.** No skill metric, naive benchmark, forecast interval or calibration check (section 4.8).
- **Time-varying productivity in the forecast.** Module 5 teaches an AR(1) year effect, but Module 8 forecasts from a stationary Ricker. A random-walk version removes the bias on these data **[run V06-F, V11-D]**.
- **Delayed density dependence (the Larkin model) and cyclic dominance.** The data file carries a `cycle` column that no module uses. By AIC, the plain Ricker is within 2 units of the best of four simple alternatives for only 4 of 19 stocks **[run V11-C]**.
- **Biological reference points** (S_MSY, U_MSY, S_gen, S_max) with uncertainty. The Lambert W function is introduced but never applied to anything.
- **A fitting and convergence routine.** The slides make one default `nlminb()` call. No raised iteration limits, Newton steps, eigenvalue check or restart (section 4.1).
- **REML.** The pre-course notes discuss the bias of ML variance estimates but never give the general remedy. It matters here: full REML raises the between-stock SD from 0.106 to 0.18 **[run V14]**.
- **Checks on the random-effects part of a model.** Every residual shown tests the observation model only.
- **Model comparison.** AIC, likelihood-ratio tests and cross-validation do not appear anywhere in the material.
- **Bias correction of derived quantities, and prediction intervals** that include process error.
- **Known biases of spawner–recruit regression** (time-series bias, spawner measurement error) and the state-space models built to deal with them.
- **Simulation self-tests** (simulate from the fitted model, refit, compare).
- **A real Bayesian section.** `modules/7b_tmbstan.qmd` is a mis-titled copy of the pre-course notes.

### Is anything in it that should not be?

Yes. The items I would remove or rewrite:

- `log((age + 1e-4)/spawners)` to handle zero recruits (Module 8). The constant drives the answer. For Quesnel the forecast changes 2.7-fold as the constant goes from 10⁻⁶ to 10⁻², and the age-4/age-5 correlation changes sign **[run V06-D]**.
- Estimating a penalty weight (`loglambda`, `logpen`) as an ordinary parameter (GAM and latent-variable examples). The optimiser sends it to zero, so there is no penalty **[run V07-A, V08-D]**.
- A covariate design matrix with an intercept added to a model that already has intercepts (Module 8). The Hessian is singular and the standard errors are `NaN` **[run V06-C]**.
- The custom `dlognorm()`, `rdirmulti()`, `ddirmulti()` likelihood call and `ddirmulti_osa()` (Modules 6 and 7a). One crashes in simulation mode, one samples the wrong distribution, one sums densities where it should sum log-densities, and one is unfinished **[run V04-A, V05-C/D/E]**.
- Global on/off flags that switch branches inside the objective function (Module 6). `obj$simulate()` and `obj$report()` follow the flag's current value, not the value the model was built with **[run V04-A]**.
- The multi-start "convergence check" that silently discards runs that did not converge. In my replication it discards 46 to 49 of 50 **[run V02-B]**.
- "The Hessian is positive definite if the variance can be computed." It is not a reliable test **[run V02-E, V06-C]**.
- `checkConsistency(obj)` with default settings read as a test of the Laplace approximation. With 40 parameters and 100 simulations it rejects in 4 of 5 runs for a model where the approximation is exact **[run V03-F]**.
- `oneStepGeneric` residuals for a Gaussian model. It is more than 50 times slower than `fullGaussian`, which is exact here **[run V04-F]**.
- The statement that direct maximisation of a mixture likelihood "is not possible" **[run V09-B]**.
- "Carrying capacity" as the name for the Ricker β, and `Eage4` ("expected") as the name for a median.
- Calls that rely on default arguments: `obj$report()`, `obj$simulate()`. TMB's own help page says not to **[doc]**.

### Are more modern or sophisticated approaches available?

Yes. Section 5 has a table with the maturity of each. The ones most relevant to you:

- **Recent RTMB features** (versions 1.8 to 2.0, October 2025 to September 2026): differentiable `uniroot()` and `integrate()`, `obj$force.update()`, `dflat()`, a "tips" vignette, and a published paper to cite **[doc]**.
- **`mgcv::smooth2random()`** for smoothers inside an RTMB likelihood, as sdmTMB, glmmTMB and brms do. In my test it reproduces `mgcv::gam(method = "ML")` to 10⁻⁵ **[run V07-B]**.
- **State-space spawner–recruit models** with random-walk productivity. The `samEst` package from the same organisation that hosts this workshop already implements them, with leave-future-out cross-validation **[doc]**.
- **Leave-future-out validation and proper scoring** for forecasts.
- **SparseNUTS and StanEstimators** as newer routes to Bayesian sampling from RTMB objects, alongside tmbstan **[doc]**.
- **RTMBdist** (on CRAN) for extra distributions, an AD Lambert W, and an importance-sampling check of the Laplace approximation **[doc]**.
- **One-sample "process residuals"** to check random effects, and the `compResidual` approach for composition data.

### The findings that matter most

| # | Finding | Where | Evidence |
|---|---|---|---|
| 1 | The Module 8 Ricker forecast over-forecasts Stellako in 21 of 22 hindcast years (mean error +0.14 million on a mean return of 0.27) and does not beat "same as 4 years ago" on mean absolute error. A random-walk productivity version is unbiased and has the lowest mean absolute error and RMSE of everything tried. | `8_practical.qmd` | V06-F |
| 2 | Across all 19 stocks, the stationary Ricker over-forecasts log(R/S) in 69% of 379 leave-future-out forecasts, by a factor of about 1.8 on average. The random-walk version: 49%. | — | V11-D |
| 3 | The `+1e-4` zero fix: changing the constant between 10⁻⁶ and 10⁻² changes a forecast by up to 2.7-fold. | `8_practical.qmd:106` | V06-D |
| 4 | The covariate model is not identifiable. Only 1 of 5 covariates has \|z\| > 2, and AIC is worse than with no covariates. | `8_practical.qmd:189` | V06-C |
| 5 | GAM example: λ goes to 2.5 × 10⁻¹⁰. The "penalised" fit equals the unpenalised one to 3 × 10⁻⁶. | `GAM.qmd:193–231` | V07-A |
| 6 | Latent-variable example: no convergence from any of 3 seeds, the penalty weight goes to 10⁻⁸, and a post-processing bug changes the sign of 29 of 171 stock correlations. | `latent_variable_model.qmd` | V08 |
| 7 | Mixture example: with sensible starting values, direct optimisation and EM reach the same likelihood in 50 of 50 simulated data sets, and direct needs half the iterations. | `mixture_model.qmd:27` | V09-B |
| 8 | The first three likelihood slides return the log-likelihood from a function named `negloglik`. `nlminb()` then "converges" (code 0) to σ ≈ 10⁻¹⁴⁸. | `2_likelihoods.qmd:305–348` | V01-A |
| 9 | ML drives the between-stock SD to zero in 79% of data sets simulated from the fitted model. Full REML: 28%. On the real data the 95% profile interval for that SD runs from 0 to 0.25. | Module 5 model | V13-C, V14 |
| 10 | `checkConsistency()` used as on the slide gives false alarms. | `5_random_effects.qmd:311` | V03-F |
| 11 | The fitted AR(1) model's own residuals show stock-specific autocorrelation (0.16) and heavy tails, and its year-effect innovations fail a normality check in 99% of joint draws. None of this is visible in the slides' QQ-plot. | `6_residuals.qmd` | V04-F, V13-B |
| 12 | Data: the covariate CSV starts with a byte-order mark, so the join on `yr` fails outside UTF-8 locales. `recruits` disagrees with the age columns in 22 rows, 18 of them brood year 2018. | `data/sockeye/` | V06-A, V06-B |

---

## 2. What the workshop does well

These are worth keeping, and they are why the criticisms below are fixable.

- **The progression.** `lm()` → a hand-written likelihood → the same likelihood in RTMB, on one data set throughout.
- **It shows a real numerical failure and its cure.** `examples/ricker_example.qmd` shows `nlminb()` missing the `lm()` answer and traces it to the scale of the spawner numbers. Per-stock spawner SDs span a factor of 170 **[run V02]**, so scaling is the right fix.
- **AD is compared with finite differences by counting function calls** (115 for a gradient, 6,556 for a Hessian in the rendered slides).
- **The Laplace approximation is checked against closed forms** (Poisson–gamma against the negative binomial, normal–normal against the multivariate normal).
- **Three equivalent codings of an AR(1) effect** (a loop, `dautoreg()`, a sparse precision matrix). Once one typo is fixed they agree to 12 significant digits **[run V03-C]**.
- **Four kinds of confidence interval on one parameter** (Wald on the log scale, Wald on the natural scale, profile, bootstrap).
- **Jacobians are explained and used correctly** in the prior for the tmbstan example.
- **The right toolbox is introduced**: `map`, `OBS()`, `%~%`, `ADREPORT()`, prediction by augmenting the data with `NA`, `checkConsistency()`, one-step-ahead residuals, and the closure pattern for passing data that RTMB's documentation recommends.
- **The examples are candid.** "This is certainly not working that great." "Very suspicious… Likely overparameterized." Those instincts were right.

---

## 3. Errors and poor practice, by module

Severity: **High** means wrong numbers, a crash, or a wrong conclusion. **Medium** means it works here but misleads, or breaks under small changes. **Low** means cosmetic. Line numbers refer to the `.qmd` source.

A caveat that applies throughout: I reviewed the files, not the live sessions. Some of this may have been corrected or explained aloud.

### 3.1 Pre-course notes, README and Module 1

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| P-1 | Medium | `0_precourse_material.qmd:39,43`; `2_likelihoods.qmd:365`; `8_practical.qmd:52` | The Ricker γ (later β) is called "carrying capacity". In `R = S·exp(α − βS)`, β is the strength of density dependence, in units of 1/spawners. Recruitment peaks at `S_max = 1/β`, and the unfished equilibrium, which is what "carrying capacity" normally means, is `α/β`. For Chilko: β = 2.58 per million spawners, S_max = 0.39 million, equilibrium = 0.94 million **[run V11-A]**. `samEst` reports `Smax <- 1/beta` **[doc]**. | Call β "density dependence (1/S_max)". |
| P-2 | Medium | `1_introduction.qmd:32` | "Run R code in C++ for speed" gives the wrong mental model. RTMB runs your R function once, records the arithmetic on a "tape", and replays that tape. It does not translate R to C++. This is why an `if`, a loop length or an array size cannot depend on a parameter. The speed gain is mostly in gradients: the rendered Module 2 benchmark shows no gain for a vectorised likelihood (12.7 vs 12.4 µs) **[read]**. | Say "records and replays the calculation, with exact derivatives". |
| P-3 | Low | `README.md:30` | The schedule promises "Compare output with `glmmTMB`". No module does. | Script V14 does it. The two agree. |
| P-4 | Low | `0_precourse_material.qmd:67,147` | `$beta_2$` (missing backslash), `nlmbin`. | — |

### 3.2 Module 2: writing a likelihood

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| M2-1 | High | `2_likelihoods.qmd:305–310, 328–333, 343–348` | `ll <- -sum(dnorm(..., log = TRUE)); return(-ll)` flips the sign twice, so the function called `negloglik` returns the log-likelihood. The rendered slide prints −704.37 for it. Minimising it does the opposite of fitting: `nlminb()` reports convergence code 0 at μ = −28, σ ≈ 10⁻¹⁴⁸ **[run V01-A]**. The slide asks "Problems with this code?", but the answer given on the next slide is the missing `exp()` on σ, and the sign error survives into both later versions. | One minus sign (below). |
| M2-2 | Medium | `2_likelihoods.qmd:268–294` | The "first taste of RTMB" benchmark does not time a compiled loop. The loop never touches the tape's input `a`, so all 1,000 density evaluations happen once, while taping, and are stored as one constant. The tape has a single multiplication **[run V01-B]**. The conclusion (loops are cheap on the tape) is true, but this benchmark does not show it. | Make the loop depend on the parameter. Then the tape has 10,001 operations and is about 175 times faster than the R loop **[run V01-B]**. |
| M2-3 | Low | `2_likelihoods.qmd:90, 319` | `f(x_1, …, f_n)`; "log maps a positive number to {∞, ∞}". | `x_n`; (−∞, ∞). |

```r
negloglik <- \(theta){
  getAll(theta)
  -sum(dnorm(x, mu, exp(logsigma), log = TRUE))
}
```

### 3.3 Module 3: gradients and optimisation

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| M3-1 | Medium | `3_gradients.qmd:290–294` (repeated in Modules 4–6 and the examples) | Starting values are drawn at random with no seed, so the slides are not reproducible. `logsigma = runif(nstocks, 3, 5)` starts σ between 20 and 148 when the SD of the data is 1.2. From these starts 1 run in 50 hits the iteration limit **[run V02-B]**. | Set a seed. Start from the data: `alpha = mean(logRS)`, `logsigma = log(sd(logRS))`. |
| M3-2 | Low | `3_gradients.qmd:440–441` | `par_init` is computed and then `nlminb(obj$par, …)` ignores it. | Pass `par_init`. |
| M3-3 | Low | `3_gradients.qmd:370–371` | "Newton based functions to use `optim` or `nlminb`… If gradient not provided, then finite differences are used." `optim()`'s default method is Nelder–Mead, which uses no derivatives at all. | Name `nlminb`, or `optim(method = "BFGS")`. |
| M3-4 | Low | `3_gradients.qmd:124, 130, 144, 171` | "Identify inflection point" is used for classifying a maximum or minimum. An inflection point is where curvature changes sign. The partial-derivative denominators lack ∂. | "Classify the stationary point". |

### 3.4 Module 4: inference

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| M4-1 | High | `4_inference.qmd:142, 147`; again at `:317, 325` | `beta_real <- beta*transform_beta` is meant to put β back on the scale of unscaled spawners. But `transform_beta <- sox$spawners/sox$sd_spawners` is the scaled spawner vector itself, one value per observation. A length-19 vector is recycled against 1,245 values, the result is never reported, and `ADREPORT(beta)` is still on the scaled scale. So the "real scale" standard errors and intervals on the next slides are not on the real scale **[run V02-A]**. | Divide by the per-stock SD (below). That matches `lm()` on unscaled spawners to 4 × 10⁻⁴ **[run V02-A]**. |
| M4-2 | High | `4_inference.qmd:114–124` | The multi-start check keeps a run only `if(mle$convergence == 0)` and reports the spread among the survivors. With starts drawn from N(0, 3²) and default `nlminb()` limits, 48, 49 and 46 of 50 runs failed for three seeds, almost all on the 150-iteration limit **[run V02-B]**. The slide's "max − min = 2.6 × 10⁻⁸" therefore describes a handful of runs. With limits raised all 50 return code 0, yet at least one sits 0.43 log-likelihood units above the optimum. | Raise the limits, report how many runs failed, and compare all of them. |
| M4-3 | Medium | `4_inference.qmd:98–103` | "Hessian should be positive definite. This one will be obvious if the variance can be computed." Not so. A matrix can be inverted, and give positive "variances", while having a negative eigenvalue **[run V02-E]**. In Module 8's covariate model `solve()` works, `sdreport()` reports `pdHess = TRUE`, and the standard errors are `NaN` **[run V06-C]**. | Look at the eigenvalues, and at whether every standard error is finite (section 4.1). |
| M4-4 | Medium | `4_inference.qmd:279` | "Find 95% area under profile likelihood." A profile-likelihood interval is the set of values whose negative log-likelihood is within `qchisq(0.95, 1)/2 = 1.92` of the minimum. That is how `confint.tmbprofile` is written **[doc]**. For the parameter on the slide the two definitions give (−1.92, −0.25) and (−1.82, −0.35) **[run V02-C]**. | Describe it as a likelihood-ratio cut. |
| M4-5 | Medium | `4_inference.qmd:311–339` | The bootstrap swaps the data through `DataEval()`. A tape evaluated twice at the same parameter values does not notice that the data changed: it returned 1829.64 where the right value was 1854.20 **[run V02-D]**. The loop as written gives correct refits only because each `nlminb()` starts somewhere new. Refits also start from the poor `obj$par`, and non-converged refits are kept. | Call `obj$force.update()` after changing the data (documented in `?DataEval` since RTMB 1.8 **[doc]**). Start refits at `fit$par`. Drop or count failures. |
| M4-6 | Low | `4_inference.qmd:132, 219, 223, 285` | `sdeport`. The Wald formula mixes β̂ and sd(α̂). "Confidence intervals are transformation invariant" needs care: you may transform the end points of an interval, but the Wald recipe gives different intervals on different scales, which is the slide's own next point. `name = 31` is hard-coded. | Select the parameter by name. |

```r
sd_s <- as.numeric(tapply(sox$sd_spawners, sox$stock, unique))   # one SD per stock
# inside the objective function:
beta_real <- beta / sd_s          # because spawners_scaled = spawners / sd_s
ADREPORT(beta_real)
```

### 3.5 Module 5: random effects

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| M5-1 | Medium | `5_random_effects.qmd:103–104` | The first model is built without `random =`, so `nlminb()` maximises the joint density of data and random effects. That density has no maximum: it grows without limit as `sigma_re → 0` with every `re → 0`. The run ends with `sigma_re ≈ 10⁻⁸` and "false convergence" **[run V03-A]**. The next slide explains why this is wrong, but the code slide is not labelled as a counter-example. | Label it, or add `random = "re"`. |
| M5-2 | High | `5_random_effects.qmd:521` | The GMRF version uses `ll <- ll + dgmrf(...)` before `ll` exists: `object 'll' not found` **[run V03-B]**. If a variable called `ll` happens to exist in the workspace, the code runs and silently adds it. | Start with `ll <- dgmrf(...)`. Building `Q` from constant sparse matrices times AD scalars is the documented idiom **[doc]** (V03-C has a version). |
| M5-3 | High | `5_random_effects.qmd:311–315` | `checkConsistency(obj)` followed by `s$marginal$p.value` is presented as the check on the Laplace approximation. The p-value refers `n·m′S⁻¹m` to a χ² distribution, which is only valid when the number of simulations is much larger than the number of parameters **[doc: TMB source]**. This model has 40 fixed effects and the default is 100 simulations. For this normal–normal model the Laplace approximation is exact, and yet five default runs gave p = 0.004, 0.24, 0.002, 0.027 and 0.000005. With 4,000 simulations p = 0.25 **[run V03-F]**. | Use many more simulations than parameters, and read the estimated bias in units of standard error (below). |
| M5-4 | Medium | `5_random_effects.qmd:209` | "With more data, the likelihood is smooth and has a sharp peak." What the Laplace approximation needs is information about each random effect, not more data in total. For binary data with a random intercept, the approximation underestimates the SD by 22% with 2 observations per group, 5% with 5, 1% with 20. With 6,000 observations arranged as 3,000 groups of 2 the error is 23%; arranged as 60 groups of 100 it is 0.2% **[run V12]**. This matches Joe (2008) **[lit]**. | "With more data per random effect." |
| M5-5 | Medium | throughout Modules 5, 6 and 8 | `obj$report()` and `obj$simulate()` are called with no argument. They use the last parameters evaluated, not the best ones. Straight after `nlminb()` these coincided in my runs, but after one more evaluation `obj$report()$mu` was off by up to 1.36 **[run V03-D]**. Profiling, `checkConsistency()` and `oneStepPredict()` all evaluate the function elsewhere. `?TMB::MakeADFun` says: "Do not rely upon the default arguments of any of the functions in the model object" **[doc]**. | `obj$report(obj$env$last.par.best)`. |
| M5-6 | Low | `5_random_effects.qmd:26, 61, 130, 163–172, 240` | The first equation has `+β` where the code has `−β`. "All models are numerically the same" should be "mathematically equivalent". The zero-inflation formula lacks the indicator: `p(x) = p(z=0)·1[x=0] + p(z=1)·p(x\|z=1)`. The Laplace slide defines σ as `1/(−(log f)″)`, which is the variance, and writes ∫ from ∞ to ∞. `optimize(…, lower = 0)` restricts log λ to be positive. | — |
| M5-7 | Low | `5_random_effects.qmd:365–590` | Every time-varying chunk is `eval: false`, so none was run when the slides were built. That is how M5-2 got through. | Evaluate them, with `cache: true` if speed is the worry. |

```r
chk <- checkConsistency(obj, n = 2000)       # far more simulations than fixed effects
s   <- summary(chk)
s$marginal$p.value
max(abs(s$marginal$bias) / sqrt(diag(sdreport(obj)$cov.fixed)))   # 0.04 here: negligible
```

One more point on the "ways to code random effects" slide (`:46–63`). For the Laplace approximation the centred and non-centred forms give the same likelihood, because the change of variable is linear. The choice still matters in TMB, for a reason the slide does not give: speed. Writing a random walk or AR(1) through its innovations makes the Hessian dense. For 75 years the state form has 149 non-zero entries and the innovation form 2,850 **[run V03-G]**. NOAA's FIMS documentation reports the same contrast **[doc]**. The loop and `dautoreg()` codings on the later slides are the state form, which is the right choice.

### 3.6 Module 6: custom functions

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| M6-1 | High | `6_custom_functions.qmd:181–197` | The custom `dlognorm()` calls `dnorm(log(x), …)` before it tests for simulation mode. RTMB's `dnorm()` already simulates through the `log()`, so the explicit simulation call that follows tries to assign the same observation twice: `A simulation can only be assigned to once!` **[run V04-A]**. | Test for simulation mode first (below), or use the built-in `dlnorm()`. |
| M6-2 | High | `6_custom_functions.qmd:228–233, 248–255, 307–317` | Branches are chosen by global flags (`lognormal`, `nb`, `lgamma`). The tape is fixed when the object is built, but `obj$simulate()` and `obj$report()` re-run the R function and follow the flag's current value. On the "Checking Simulation" slide the flag is `FALSE` by then, so both `sim1` and `sim2` come from the `dnorm(log(RS))` branch **[run V04-A]**. The plot compares a thing with itself, and that is why M6-1 never surfaced. | Pass switches through the data closure, like the data. |
| M6-3 | Medium | `6_custom_functions.qmd:274–280` | The slide defines `dlgamma(x, shape, rate)`. RTMB already exports `dlgamma(x, shape, scale)` **[doc]**. The same positional call means two different things depending on which definition is in scope **[run V04-B]**. | Rename it, e.g. `dloggamma_rate()`. |
| M6-4 | Medium | `6_custom_functions.qmd:398–432` | The multivariate logit-normal section stops before it defines a density, and `mlogit_logdet()` is never used. The tape-inside-a-tape plus `ADjoint` log-determinant is unnecessary. `determinant()` works on AD matrices directly, and the Jacobian has a closed form, `log\|J\| = −Σ log p_i` **[run V04-D]**. | Finish it with the closed form, or cut it. |
| M6-5 | Medium | `6_custom_functions.qmd:186`; `7_misc_coding.qmd:123`; `7a_custom_distribution.qmd:27`; `mixture_model.qmd:133` | `RTMB:::dGenericSim`, `RTMB:::ad_context`, `RTMB:::getValues` are not exported. Their export lines are commented out in RTMB's NAMESPACE **[doc]**. They can change without notice, and CRAN's checks flag `:::` calls into another package. RTMBdist does use the same helpers, through `asNamespace()` **[doc]**, so this is a convention with a known risk. | Say so on the slide. For simulation, the pattern below needs no internals. |
| M6-6 | Low | `6_custom_functions.qmd:153–202` | The slides never mention that `dlnorm()` is built in. It supports simulation and `oneStepGeneric` residuals, though not the `cdf` method, because RTMB 2.0 has no AD `plnorm()` **[run V04-C]**. Putting the observation on the log scale with `dnorm()` supports every method. | Mention it. |

```r
dlognorm <- function(x, mean, sd, log = FALSE){
  if (inherits(x, "simref")) {                 # simulation mode: assign the draw and stop
    x[] <- exp(rnorm(length(x), mean, sd))
    return(0)
  }
  y <- log(x)
  logpdens <- dnorm(y, mean, sd, log = TRUE) - y
  if (log) logpdens else exp(logpdens)
}
```

This version gives the same draws as the `dnorm(log(RS))` model with the same seed **[run V04-A]**.

The slide's question, "Does it matter if I include the Jacobian for my observed distribution?", has a two-part answer that is worth stating. For the estimates, no: the two objectives differ by the constant Σ log(RS) = 135.81 and give the same parameters **[run V04-C]**. For AIC, yes: models fitted to `log(R/S)` and to `R/S` can only be compared if that constant is included. For a random effect it always matters. Leaving the Jacobian out of the Poisson–gamma example moves the estimated shape from 4.1 to 1.8, where the exact negative-binomial likelihood gives 4.0 **[run V04-E]**. The slide sets this up but comments out the lines that would show it.

### 3.7 Module 6: residuals

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| M6-7 | Medium | `6_residuals.qmd:250` | `method = "oneStepGeneric"` for a linear Gaussian model. It took about 150 s against under 3 s for `fullGaussian`, which is exact for this model. Three of the 1,245 residuals differ from the exact ones by more than 0.01, the worst by 0.73, and setting `range` did not change that **[run V04-F]**. TMB's help does call `oneStepGeneric` "the most accurate" in general **[doc]**. Here it is neither needed nor the most accurate. | `oneStepPredict(obj, method = "fullGaussian")`. |
| M6-8 | Medium | `6_residuals.qmd:232–259` | One-step-ahead residuals condition each observation on those before it in the data vector. The data are sorted by stock, then year, so "the past" of a 1950 Late Stuart observation includes Early Stuart in 2021. Any order gives valid N(0,1) residuals, but they are different residuals: sorting by year first changes them (correlation 0.94) **[run V04-F]**. | Sort by year before building the object if "predict this year from earlier years" is what you want to check. |
| M6-9 | Medium | `6_residuals.qmd:140–260` | The only diagnostic is a QQ-plot, read as "looks okay". The exact residuals of this model have within-stock lag-1 autocorrelation of 0.16 (95% CI 0.11 to 0.22; Cultus 0.49), a Shapiro–Wilk p of 0.0015 and twice the normal share beyond ±3 **[run V04-F]**. The module also states, correctly, that these residuals say nothing about the random-effect structure, then offers no tool that does (section 4.6). | Add an autocorrelation check by stock and a check on the random effects. |
| M6-10 | Low | `6_residuals.qmd:30–32, 95, 104, 116, 202, 246` | The "Pearson residual" formula divides by √expected, which is right only for Poisson data. `logtphi`. `OBS(logRS)` is called twice. Simulated CDF values of exactly 0 or 1 give ±∞. `#\| echo: ture`. | Divide by √Var. |

### 3.8 Modules 7, 7a and 7b: miscellaneous coding, custom distribution, tmbstan

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| M7-1 | Medium | `7_misc_coding.qmd:83–95` | `make_adfun(pars, data)` calls `MakeADFun(objective_function, param, …)`. It ignores its own `pars` argument and reads a global `param`. Start values passed in are discarded, and the function fails when `param` does not exist **[run V05-A]**. | Use `pars`. |
| M7-2 | Medium | `7_misc_coding.qmd:115–133` | The `ad_context()` example cannot run: `logsigma = log(sigma)` with no `sigma` defined finds the function `stats::sigma` **[run V05-B]**. With that fixed, its "prediction interval" is simulated at the point estimate, so it ignores parameter uncertainty and uses the ML σ. Mean width 59.1 against 63.1 from `predict.lm()` **[run V05-B]**. | Give a numeric start value. Name the interval "plug-in". |
| M7-3 | Low | `7_misc_coding.qmd:154` | "You can't use uniform distributions in RTMB." A uniform prior is a constant. It needs no density at all, only a transformation or bounds to keep the parameter in range **[run V05-G]**. What you cannot do is estimate the end points of a uniform observation model. | Reword. |
| M7-4 | Low | `7_misc_coding.qmd:178–201` | `ADoverload("if")` is shown as routine. RTMB's NEWS for 1.9 lists it as "Experimental atomic branching", and it takes a scalar condition only **[doc]**. | Say "experimental". |
| M7-5 | Medium | `7_misc_coding.qmd:204–227` | The `%~%` slide says you "don't need to declare we are returning anything". It should also say that you cannot. `%~%` returns its own hidden total on exit and overrides an explicit `return()`. A prior added by hand next to a `%~%` statement is silently dropped: the function returned 309.71 where the full objective was 320.83 **[run V05-F]**. | Use one style per function. |
| M7-6 | Low | `7_misc_coding.qmd:320–350` | The tmbstan slide is correct, including the Jacobian. It is also the whole of the Bayesian content (section 4.9). | — |
| M7a-1 | High | `7a_custom_distribution.qmd:75` | `rgamma(length(alpha), 1, alpha)` draws gamma variables with shape 1 and rate α. A Dirichlet(α) needs shape α and a common rate. For α = (1, 5, 20) the simulated mean proportions are (0.70, 0.23, 0.07) where theory gives (0.04, 0.19, 0.77) **[run V05-C]**. The slide's own check uses α = (1, 1, 1), the one case where the two agree. | `rgamma(length(alpha), shape = alpha, rate = 1)`. |
| M7a-2 | High | `7a_custom_distribution.qmd:132` | `ll <- ll + ddirmulti(x[i,], size = 5, alpha = alpha)` leaves out `log = TRUE`, so the objective is minus the sum of the densities. Fitted to data simulated with α = (2, 4, 8) it returns (0.01, 2.5 × 10⁷, 7.4 × 10⁷) with "false convergence". With `log = TRUE`: (2.6, 4.8, 9.1) **[run V05-D]**. | Add `log = TRUE`. |
| M7a-3 | Medium | `7a_custom_distribution.qmd:153–172` | `ddirmulti_osa()` is unfinished. `for( a in 1:length(alpha-1) )` has no body, so R takes the next statement as the loop body, and `1:length(alpha-1)` is `1:length(alpha)` **[run V05-E]**. One-step-ahead residuals for composition data are hard to get right (Trijoulet et al. 2023 **[lit]**). RTMBdist's own `ddirmult` refuses them **[doc]**. | Remove it until it is finished and tested. |
| M7a-4 | Low | `7a_custom_distribution.qmd:53`; `examples/custom_distribution.qmd` | The density formula uses `n` where it means `M`. The file is duplicated byte-for-byte in `examples/`. | Keep one copy. |
| M7b-1 | Medium | `7b_tmbstan.qmd` | Titled "Custom Distribution", it is a copy of the pre-course notes with nothing about tmbstan. It still has a bug that Module 0 has since fixed (`:93`, `mu <- sum(X[i,] * beta)` overwrites `mu` on every pass). | Delete it or write it. |

### 3.9 Module 8: the practical forecast

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| M8-1 | Medium | `8_practical.qmd:43–45`; `data/sockeye/FC_Environmental_Data_2026_Std_Offset_Real.csv` | The covariate file starts with a UTF-8 byte-order mark. In a UTF-8 locale R drops it. In a C locale (many containers and CI runners) the first column is read as `X...yr` and `left_join(covariates, by = "yr")` stops with an error **[run V06-A]**. I expect the same on Windows with R before 4.2, but did not test that. | `read.csv(..., fileEncoding = "UTF-8-BOM")`, or re-save the file without the mark. |
| M8-2 | Medium | `data/sockeye/SRDATA2026 1.csv` | `recruits` equals the sum of each stock's main age classes (ages 4+5, or 3+4 for Harrison) in 1,194 of 1,216 rows. It does not in 22 rows, 18 of which are brood year 2018, one for each of 18 stocks. Stellako 2018 has age4 = 0.200 and age5 = 0.042 but recruits = 0.129 **[run V06-B]**. Modules 3 to 7 model `recruits`, and Module 8 models the age columns, so for 2018 they use different numbers. The README does warn that these are "not necessarily official numbers". | Check brood year 2018 with whoever maintains the file. |
| M8-3 | High | `8_practical.qmd:189, 217` | `model.matrix(~ peak + aflow + …)` includes an intercept column, and the model already has `alpha[1]` and `alpha[2]`. Three parameters carry two identifiable quantities. The Hessian has an eigenvalue of −3 × 10⁻¹², along (α₁, α₂, −intercept). `nlminb()` returns code 0 and `sdreport()` reports `pdHess = TRUE`, but the standard errors of all three are `NaN` **[run V06-C]**. Point forecasts are unaffected (they match an identifiable version to 6 × 10⁻⁶), so the plot gives no warning. | Drop the intercept and standardise the covariates (below). |
| M8-4 | Medium | `8_practical.qmd:155–201` | Five covariates are added without any check that they help. In the identifiable version only one (`mapsst`, −0.32 per SD, z = −2.4) is distinguishable from zero. Together they lower the negative log-likelihood by 3.8 for 5 parameters, which makes AIC worse by 2.5 **[run V06-C]**. The first fit uses `mapsst` and the loop uses `jnesst`. Covariates enter in raw units (salinity is 30.9 ± 0.4), which is the conditioning problem the course teaches for spawners. | Section 4.8. |
| M8-5 | High | `8_practical.qmd:106–107, 207–208` | `log((age5 + 1e-4)/spawners)`. Values are in millions, so the constant is 100 fish, chosen arbitrarily. Stellako has one zero in 73 years, and moving the constant from 10⁻² to 10⁻⁶ doubles the age-5 SD (0.77 to 1.54). Quesnel has 24 zeros of 72: the forecast runs from 0.42 to 1.14 million and the age-4/age-5 correlation from −0.15 to +0.80. Late Shuswap, 19 zeros: 0.10 to 0.23 million **[run V06-D]**. The practical invites students to try other stocks. | Model total recruits and the age split separately (section 4.8). |
| M8-6 | Medium | `8_practical.qmd:56, 62, 87–89` | `Eage4 <- exp(mu4)*spawners` is a median, not an expectation. With σ₄ = 0.94 the mean is 1.56 times larger, and the mean total return is 1.57 times the slide's `pred` **[run V06-E]**. The median is a legitimate thing to forecast, and it is the reference point DFO's Fraser forecasts use **[lit]**, but the slide should say which it is. The conditional variance is written `Σ22 − Σ12²Σ11` and should be `Σ22 − Σ12²/Σ11`; `m_4t` should be `μ_4t`; the first line needs `~` in place of `=`. The code is right. | Rename to `med4`; fix the formula. |
| M8-7 | High | `8_practical.qmd:123–153, 203–234` | The hindcast ends with a line plot. No error measure, no benchmark, no forecast interval. Section 4.8 reports what those show. | Section 4.8. |
| M8-8 | Low | `8_practical.qmd:26, 128, 194` | "Run a retrospective analysis over the last 12 years"; the code does 22. `a = rnorm(ncol(X))` is an unseeded random start. | — |

```r
covs <- c("peak", "aflow", "jnesst", "pdo", "stn2js")
Z <- scale(as.matrix(dat[, covs]))        # inside the hindcast loop: training-window mean and SD only
# in the objective:  delta_alpha <- drop(Z %*% a)      (no intercept column)
```

One thing Module 8 gets right deserves a mention. The hindcast loop blanks exactly the returns that would not yet be known when the forecast is made (age 4 from the forecast brood year, age 5 from that year and the one before). I checked the alignment and found no look-ahead in the responses **[read]**.

### 3.10 The example documents

**`examples/GAM.qmd`**

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| G-1 | High | `GAM.qmd:186–231` | The penalty `lambda/2 * b′Sb` is added to the objective and `loglambda` is estimated as an ordinary parameter, with the spline coefficients `b` also ordinary parameters (`random = "re"` only). The penalty is never negative and vanishes as λ → 0, so the optimiser takes λ there. I got λ = 2.5 × 10⁻¹⁰ with convergence code 0, an objective identical to the unpenalised fit (1745.687), and the same curve to 3 × 10⁻⁶ **[run V07-A]**. A smoothing parameter can be estimated only after the coefficients are integrated out, with their Gaussian normalising term kept. That is what REML or marginal-likelihood smoothing means (Wood 2011 **[lit]**). | Section 5.2. |
| G-2 | Low | `GAM.qmd:170, 213–216` | `fit_gam` is computed and never compared with anything. The mean is filled element by element in a loop; a vectorised line builds the object several times faster **[run V07-C]**. | Compare with it. |

**`examples/latent_variable_model.qmd`**

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| L-1 | High | `latent_variable_model.qmd:121–124` | With the document's random starts and default `nlminb()`, none of three seeds converged (two hit the iteration limit, one "false convergence"). They ended at different likelihoods (1598.8, 1599.1, 1610.1) and with a different stock's residual SD at zero: Scotch Creek twice (1.8 × 10⁻⁵ after continuing), Seymour once **[run V08-A]**. A residual variance collapsing to zero is a "Heywood case", a well-known failure of factor models. R's `factanal` help calls them "much more common than most texts… would lead one to believe" and bounds the variances below **[doc]**. | A lower bound or weak prior on the residual SDs (V08-D), fewer factors, and several starts. |
| L-2 | High | `latent_variable_model.qmd:152–153` | `Lambda <- cbind(lambda1*…, lambda2*lambda1*…)`: the second column of loadings is multiplied by the first. The correlations in the plot differ from the model's by up to 0.43, and 29 of 171 have the wrong sign **[run V08-B]**. `pList$logscale` also works only because `$` partially matches `logscale_z`. | `cbind(lambda1*exp(logscale_z[1]), lambda2*exp(logscale_z[2]))`. |
| L-3 | Medium | `latent_variable_model.qmd:158` | "Covariance between stocks" is computed from the two factors and the residual variances. It leaves out the shared AR(1) year effect, which is in the model and adds the same covariance to every pair. The mean pairwise correlation is 0.05 from the factors alone and 0.39 with the year effect **[run V08-C]**. What the plot shows is covariation left after the common year effect is removed. | Add the year-effect variance, or relabel the plot. |
| L-4 | High | `latent_variable_model.qmd:203, 242` | The "penalised" refit estimates the penalty weight `logpen`. As in G-1 it goes to zero (10⁻⁸); the objective equals the unpenalised one and the loadings do not move **[run V08-D]**. With the strength fixed (a N(0,1) penalty on loadings and a weak boundary-avoiding penalty on the residual SDs) the fit converges with a positive-definite Hessian and the smallest residual SD is 0.19 **[run V08-D]**. | Fix the strength, or treat the loadings as random with their normalising term. |
| L-5 | Medium | `latent_variable_model.qmd:174–176` | Loadings are identified only up to rotation. Fixing `lambda2[1] = 0` ties the axes to whichever stock comes first alphabetically (Birkenhead), so "Latent Variable 1 vs 2" has no meaning of its own. A varimax rotation moves the split of loading variance from 62/38 to 48/52 and leaves the implied covariance unchanged **[run V08-E]**. bayesdfa and gllvm rotate before plotting **[doc]**. | Rotate first, or interpret only `ΛΛ′`. |

**`examples/mixture_model.qmd`**

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| X-1 | High | `mixture_model.qmd:27, 80` | "Direct maximization of the log likelihood… is not possible due to the sum within a log." It is possible and it is standard: the sum inside the log rules out a closed-form solution, not numerical optimisation. Stan's guide and the RTMB-based LaMa package both fit mixtures this way **[doc]**. From the same k-means starting values, direct optimisation with AD gradients and `mixtools::normalmixEM` reached the same likelihood in 50 of 50 simulated data sets, in a median of 27 iterations against 52 **[run V09-B]**. | Correct the statement. |
| X-2 | Medium | `mixture_model.qmd:62` | What failed was the start. `mu = exp(rnorm(K, 3, 0.25))` puts all four means between about 11 and 33 when the data run from 33 to 87. From there the direct fit collapses to a single normal (NLL 381.9, identical to fitting one normal) while EM finds something better (346 to 356) **[run V09-A]**. | Start from quantiles or k-means. |
| X-3 | Medium | `mixture_model.qmd:56`, rendered output | Neither method has a true maximum to find. With a separate variance per component the likelihood grows without limit as a component narrows onto one observation **[run V09-C]**; Stan's guide makes the same point **[doc]**. The rendered document's own EM answer has a component with SD `exp(−1.52) = 0.22` when every true SD is at least 1.5. The line `log(sum(prob*dnorm(…)))`, flagged "Numerical Issues?", underflows to `−Inf` for a distant observation where log-sum-exp gives −1874 **[run V09-D]**. | `logspace_add()` (below); bound the SDs below or share one variance; order the means. |

```r
lprob <- log(itrans(logitp)); sigma <- exp(logsigma)
acc <- lprob[1] + dnorm(x, mu[1], sigma[1], log = TRUE)
for (k in 2:K) acc <- logspace_add(acc, lprob[k] + dnorm(x, mu[k], sigma[k], log = TRUE))
nll <- -sum(acc)
```

**`examples/osa_multinomial.qmd` and `examples/ricker_example.qmd`**

| ID | Severity | Where | Problem | Fix |
|---|---|---|---|---|
| O-1 | Medium | `osa_multinomial.qmd:93, 141` | `as.matrix(r$residual, nrow = np, ncol = 100)` does not reshape: `as.matrix()` has no such arguments and returns 500 × 1. `osa[-np,]` then drops one value where it should drop the 100 last-category residuals **[run V10-A]**. The effect on the QQ-plot is small here, because 92 of those 100 come back as `NaN` and `qqnorm()` drops `NaN` silently **[run V10-B]**. It would matter for `mean()`, `sd()` or `acf()`. | `matrix(r$residual, nrow = np)`. |
| O-2 | Low | `osa_multinomial.qmd:25, 117` | The text writes `log(p_i/(1 − p_5))`; it should be `log(p_i/p_5)`, as the code has it. `numeric(nrow(dat$X))` should be `ncol`; it works because R recycles zeros. | — |
| R-1 | Low | `ricker_example.qmd:74` | The restart passes `sigma(fit_lm)/n*(n-npars)` where the function expects log σ. | `log(sigma(fit_lm) * sqrt((n - npars)/n))`. |

---

## 4. What is missing

Each gap below says why it matters for these data, with evidence where I could produce it.

### 4.1 A fitting and convergence routine

Every fit on the slides is one call, `nlminb(obj$par, obj$fn, obj$gr)`, and the only checks offered are the convergence code, the largest gradient and "can the variance be computed".

That is thin for three reasons.

- **The default limits are low.** `nlminb()` stops after 150 iterations or 200 function evaluations **[doc]**. From good starting values the 57-parameter Ricker of Module 3 needed 125 iterations **[run V02-B]**. The latent-variable example hit the limit **[run V08-A]**.
- **Code 0 does not mean a usable fit.** It was returned with a singular Hessian (M8-3), with a penalty weight at 10⁻¹⁰ (G-1), with a wrong-signed objective (M2-1), and 0.43 log-likelihood units short of the optimum (M4-2).
- **Positive-definiteness is not checked properly** (M4-3).

Every established TMB package wraps the optimiser in more than this **[doc: package sources]**:

| Package | What it does around `nlminb()` |
|---|---|
| glmmTMB | `iter.max = 300, eval.max = 400`; warns when the Hessian is not positive definite; a troubleshooting vignette covers "false convergence" |
| sdmTMB | `eval.max = 2000, iter.max = 1000`; one Newton step by default; `sanity()` flags gradients above 0.001, non-finite or large standard errors |
| WHAM | `iter.max = eval.max = 1000`; three Newton steps |
| SAM | three Newton steps; `rel.tol = 1e-10`; convergence redefined as largest gradient below 10⁻⁴ |
| TMBhelper `fit_tmb` | `eval.max = iter.max = 1e4`; restarts; optional Newton steps. Last commit 2022, so treat it as unmaintained |

A version for RTMB that I tested on the AR(1) model **[run V13-A]**:

```r
fit_rtmb <- function(obj, newton_steps = 2, control = list(eval.max = 1e4, iter.max = 1e4)) {
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
  for (i in seq_len(newton_steps)) {                    # polish with Newton steps
    g <- as.numeric(obj$gr(opt$par))
    h <- optimHess(opt$par, obj$fn, obj$gr)
    step <- try(solve(h, g), silent = TRUE)
    if (inherits(step, "try-error")) break
    new_par <- opt$par - step
    new_obj <- obj$fn(new_par)
    if (is.finite(new_obj) && new_obj <= opt$objective) {   # keep a step only if it does not make things worse
      opt$par <- new_par; opt$objective <- new_obj
    } else break
  }
  invisible(obj$fn(opt$par))                            # so that last.par.best is the final point
  sdr <- sdreport(obj, par.fixed = opt$par)
  ev  <- eigen(optimHess(opt$par, obj$fn, obj$gr), symmetric = TRUE, only.values = TRUE)$values
  se  <- suppressWarnings(sqrt(diag(sdr$cov.fixed)))
  checks <- c(nlminb_code_0   = opt$convergence == 0,
              max_gradient_ok = max(abs(sdr$gradient.fixed)) < 1e-3,
              hessian_pos_def = all(ev > 0),
              all_se_finite   = all(is.finite(se)),
              no_huge_se      = all(se[is.finite(se)] < 100))
  list(opt = opt, sdr = sdr, checks = checks, par_full = obj$env$last.par.best)
}
```

On the AR(1) model a single default call stops with a largest gradient of 10⁻³. The helper takes it to 5 × 10⁻⁸ and all five checks pass **[run V13-A]**. `optimHess()` is used because `obj$he()` is not available for models with random effects **[doc]**. Then use `obj$report(fit$par_full)`, never `obj$report()`.

Two more habits: set a seed and start from data-based values (M3-1), and refit from two or three perturbed starts and report all of them (M4-2).

### 4.2 How to tell whether the Laplace approximation is good enough

The slide's statement is corrected in M5-4. The useful rule is: the approximation is exact for Gaussian random effects with Gaussian data (every model in this workshop), very good when each random effect is informed by many observations, and poor for binary or small-count data with few observations per random effect (Joe 2008; Ogden 2021 **[lit]**).

When you move beyond Gaussian data you need a check. In rising order of cost:

- **`checkConsistency()` with enough simulations** (M5-3). Read the estimated bias in units of standard error, not only the p-value. TMB's help notes that a formal test "would thus fail eventually for large n" even when the bias is negligible **[doc]**.
- **Exact integration where it is cheap.** RTMB 1.9 made `integrate()` differentiable **[doc]**, so a model with one random effect per group can have each group integrated exactly inside the likelihood. Compare the estimates with the Laplace fit.
- **`RTMBdist::laplace_check()`**, an importance-sampling check **[doc]**. I did not run it.
- **MCMC without the approximation.** tmbstan samples the random effects by default (`laplace = FALSE`). Comparing fixed-effect posteriors with and without it is the test proposed by Monnahan and Kristensen (2018) **[lit]**.

One correction to a common belief: `aghq::marginal_laplace_tmb()` improves the integration over the *fixed* effects. It keeps TMB's Laplace approximation for the random effects **[doc]**.

### 4.3 REML, and how far to trust a variance component

The pre-course notes show that the ML estimate of a standard deviation is too small and patch it with `sqrt(n/(n − a))`. The general remedy is REML: integrate the mean parameters out, the same way as random effects. In RTMB that is one argument, and RTMB's own help page shows it as the "REML trick" **[doc]**.

```r
obj_ml   <- MakeADFun(f, list(mu = 0, logsd = 0))
obj_reml <- MakeADFun(f, list(mu = 0, logsd = 0), random = "mu")   # REML
```

For a normal sample this returns `sd(x)` exactly (0.8314 against 0.8230 from ML) **[run V03-E]**.

It matters for the workshop's own hierarchical Ricker, where 19 slopes and an intercept are estimated alongside the between-stock SD.

| Estimator of the between-stock SD of log-productivity | Real data | Simulated with true SD = 0.106: mean estimate | Share of simulations estimated at zero |
|---|---|---|---|
| ML (as on the slides) | 0.106 | 0.021 | 79% |
| REML, common intercept integrated out | 0.127 | 0.028 | 74% |
| REML, intercept and all slopes integrated out | 0.18 (glmmTMB) | 0.080 | 28% |

Sources: **[run V13-C, V14]**, 150 simulated data sets per row.

On the real data the profile likelihood for this SD is almost flat down to zero. A value of 0 is only 0.28 log-likelihood units worse than the estimate, and the 95% profile interval is 0 to 0.25 **[run V13-C, V14]**. So the hierarchical variance, the centrepiece of Module 5, is barely identified in this model. That is not a reason to drop the model. It is a reason to report the interval, use REML, and not read much into the point estimate.

To get full REML in RTMB, write the slopes without `exp()` so the model is linear in them, and add them to `random` (V13-C shows it). Two cautions from the package documentation **[doc]**: REML likelihoods cannot be compared between models with different fixed effects (sdmTMB), and for non-Gaussian responses "simulations should be done first to verify" (glmmTMB). Where a variance estimate collapses to zero, a weak boundary-avoiding penalty is the standard remedy (Chung et al. 2013 **[lit]**).

### 4.4 Uncertainty for derived quantities and predictions

`ADREPORT()` gives delta-method standard errors, and Module 4 explains that well. Three things are missing.

- **Bias correction.** A non-linear function of random effects, such as `exp()` of a year effect, is biased when you plug in the estimated effects. `sdreport(obj, bias.correct = TRUE)` corrects it (Thorson and Kristensen 2016 **[lit]**; **[doc]**). For a yearly productivity index `exp(alpha0 + re_time)` in the AR(1) model the correction is +2% in years with data and +14% in a year without **[run V11-B]**. Only the estimate is corrected by default, not its standard error.
- **Prediction intervals.** A forecast interval needs estimation uncertainty plus process error. Module 7's simulated interval leaves out the first (M7-2). Module 8 has no interval. V06-F builds one from `sdreport()` covariances plus lognormal process error in about ten lines.
- **Median or mean** (M8-6).

A widely used alternative for complicated derived quantities is to sample fixed and random effects from their joint normal approximation (`sdreport(getJointPrecision = TRUE)`), which is what sdmTMB's `predict(nsim = )` does **[doc]**. It assumes normality on the internal scale and does not capture estimator bias.

### 4.5 Model comparison

The material never compares two models formally. "AIC" does not occur in it **[read]**. What a short section would need:

- **Marginal AIC** from any fit: `2 * opt$objective + 2 * length(opt$par)`. Here the AR(1) year effect improves the marginal log-likelihood by 233 units for 3 parameters **[run V03-C]**.
- **Its limits.** Testing a variance against zero is a boundary problem, so the usual likelihood-ratio p-value is about twice too large (Self and Liang 1987; Stram and Lee 1994 **[lit]**). Marginal AIC targets prediction for new groups. Conditional AIC targets prediction for the same groups and is implemented as `cAIC()` in sdmTMB and tinyVAST, following Zheng, Cadigan and Thorson (2024) **[doc]**.
- **Comparability.** Likelihoods compare only on one response scale (the Jacobian point in section 3.6), and REML likelihoods only between models with the same fixed effects.
- **For forecasting, out-of-sample tests beat information criteria.** Leave-future-out cross-validation refits on data up to each origin and scores the next value (Bürkner et al. 2020 **[lit]**). `samEst::tmb_mod_lfo_cv()` does this for Ricker models **[doc]**. Holt and Michielsens (2020) report that time-varying models fitted poorly by AICc and BIC "despite being less biased" **[lit]**, which is a warning against relying on AIC alone here.

### 4.6 Checking the random-effects part of a model

Module 6 says, correctly, that ordinary residuals are "not about random effect structure". It then offers nothing that is.

The standard tool is to take **one draw** of the random effects from their joint distribution given the data, and check that draw against what the model assumes: innovations should be independent standard normals. Using the estimated random effects themselves is too optimistic, because they are shrunk (Waagepetersen 2006; Thygesen et al. 2017 **[lit]**; TMB's documentation describes the same single-sample method **[doc]**). sdmTMB's residual vignette makes the same point about observation residuals and recommends a draw **[doc]**.

Applied to the AR(1) model **[run V13-B]**:

| Check on the year-effect innovations | One draw | Mean over 500 draws | If the model were right |
|---|---|---|---|
| Standard deviation | 1.00 | 1.01 | 1 |
| Shapiro–Wilk p-value | 0.000 | rejects in 99% of draws | rejects in 5% |
| Lag-1 autocorrelation | −0.21 | −0.18 | 0 ± 0.12 |
| Mean of first half of years minus second half | 0.55 | 0.58 | 0 ± 0.23 |
| Most extreme innovation | −3.8 | −3.5 | about ±2.5 |

The stock effects pass (normality rejected in 2% of draws). The year effect does not: there is at least one extreme year (the largest negative year effect is brood year 2015), negative serial correlation, and a drift from positive innovations early to negative ones late. A stationary Gaussian AR(1) is the wrong shape for these year effects. The drift is what a random walk, or a trend, would absorb.

### 4.7 Simulation self-tests

RTMB makes it nearly free to simulate from a fitted model and refit. The material uses simulation for a bootstrap and for `checkConsistency()`, but never as a test of the estimator. TMB's documentation says "it is generally recommended to test the implementation on simulated data" **[doc]**.

The table in section 4.3 is what such a test finds: an estimator that returns zero four times in five when the truth is 0.106. I would not have guessed that from the fit, and neither the convergence code nor the Hessian flagged it. `checkConsistency(obj, estimate = TRUE)` runs the same kind of study in one call **[doc]**.

### 4.8 Salmon stock–recruit and forecasting content

**(a) Reference points.** Module 6 introduces the Lambert W function as an `ADjoint` example and then never uses it. Its use in salmon work is the exact S_MSY and U_MSY of the Ricker model (Scheuerell 2016 **[lit]**). With delta-method standard errors, for Chilko **[run V11-A]**:

```r
Smax <- 1 / beta                                    # spawners giving maximum recruits
Seq  <- alpha / beta                                # unfished equilibrium
Umsy <- 1 - LamW(exp(1 - alpha))                    # LamW as defined on the Module 6 slide
Smsy <- Umsy / beta
x_gen <- uniroot(function(x) log(x) + alpha - beta * Smsy * x,    # S_gen / S_msy, always in (0, 1)
                 lower = 1e-10, upper = 1, tol = 1e-12)$root
Sgen <- x_gen * Smsy
ADREPORT(Smax); ADREPORT(Seq); ADREPORT(Smsy); ADREPORT(Umsy); ADREPORT(Sgen)
```

| | Estimate (millions) | CV |
|---|---|---|
| S_max | 0.39 | 19% |
| S_eq | 0.94 | 15% |
| S_MSY | 0.31 | 17% |
| U_MSY | 0.80 (a rate) | 3% |
| S_gen | 0.030 | 29% |

A trap I fell into and that is worth knowing: in AD mode `uniroot()` needs plain numbers for `lower` and `upper`. Bounds that depend on parameters are silently coerced and give a wrong root. Rescale the unknown so that the bounds are constants, as above **[run V11-A; doc: RTMB source]**. RTMBdist also exports a differentiable `lambertW()` **[doc]**. Whether to apply the lognormal correction `α + σ²/2` before computing reference points is a choice the course should state; `samSim` applies it when `biasCor` is on **[doc]**.

**(b) Known biases of regressing log(R/S) on S.** Two are long documented **[lit]**. Time-series bias arises because spawners in one generation are the survivors of earlier recruitment (Walters 1985; Korman et al. 1995). Measurement error in spawners appears on both sides of the regression and distorts the apparent density dependence (Walters and Ludwig 1981; Kehler et al. 2002). The modern response is a state-space model in which recruitment, age at maturity and observation error are modelled together (Fleischman et al. 2013). Staton et al. (2020) report that conclusions "were largely consistent among state-space models but differed strongly from regression-based approaches". The per-stock regression is fine as a teaching likelihood. It should carry this caveat, and a state-space Ricker is a natural capstone for an RTMB course, because it is exactly what the Laplace approximation is for.

**(c) Delayed density dependence and cyclic dominance.** Several of these stocks have one dominant cycle line. Mean spawners on the largest line are 15 times the median line for Late Shuswap, 14 for Scotch Creek, 8 for Late Stuart **[run V11-C]**. The data file has a `cycle` column that no module uses. The standard alternative is the Larkin model, which adds spawners from the three previous brood years (Larkin 1971 **[lit]**). DFO's forecasts rank Ricker, power and Larkin forms stock by stock **[lit]**.

A direct test on the workshop data, by AIC on identical rows **[run V11-C]**:

| Best of four models | Stocks |
|---|---|
| Larkin | Chilko, Gates, Late Shuswap, Portage, Scotch Creek, Seymour |
| Larkin + random-walk productivity | Quesnel, Raft, Stellako, Upper Pitt River |
| Ricker + random-walk productivity | Birkenhead, Bowron, Early Stuart, Fennel Creek, Late Stuart |
| Plain Ricker | Cultus, Harrison, Nadina, Weaver Creek |

The plain Ricker is within 2 AIC units of the best model for 4 of 19 stocks. Read this as "the stationary Ricker is mis-specified for most of these stocks", not as a ranking of mechanisms. Lagged spawners and drifting productivity can explain the same pattern, these are in-sample comparisons with unconstrained lag coefficients, and the evidence in the literature for the Larkin mechanism is mixed (Peterman and Dorner 2012 found it sufficient only for Quesnel **[lit]**).

**(d) Time-varying productivity in the forecast.** Module 5 fits an AR(1) year effect, and Module 8 then forecasts from a model with constant productivity. A random walk in log-productivity is the conventional tool for persistent change (Peterman et al. 2003; Peterman and Dorner 2012 **[lit]**). Wor et al. (2025), from the group behind `samEst`, conclude that random-walk models "tended to perform equally well or outperform regime shift and stationary models" **[lit]**. Adding it to the Module 8 model takes two lines of likelihood and one more random effect (`a_t` in V06-F):

```r
mu4 <- a_t + alpha[1] - beta*spawners               # a_t: one value per brood year
nll <- nll - dnorm(a_t[1], 0, tau_a, log = TRUE) - sum(dnorm(diff(a_t), 0, tau_a, log = TRUE))
# MakeADFun(..., random = "a_t")
```

Two cautions from the literature **[lit]**. Benchmarks based on S_MSY drift downward when productivity declines, so DFO's 2025 guidance is to use time-varying models only when several lines of evidence support them (Holt et al. 2025). And for cyclic stocks a random walk can absorb cycle-line effects.

**(e) Forecast evaluation.** This is the gap with the largest practical consequence. I reran the Module 8 hindcast for Stellako (22 return years, 2004 to 2025, mean observed return 0.266 million) and scored it against the naive models DFO's forecast documents use **[run V06-F]**.

| Forecast | Mean error | Mean absolute error | RMSE | MAE relative to R1C |
|---|---|---|---|---|
| Slide: age-structured Ricker | +0.144 | 0.241 | 0.329 | 1.14 |
| Slide: Ricker + 5 covariates | +0.120 | 0.223 | 0.326 | 1.06 |
| Same, covariate effects shrunk as random effects | +0.145 | 0.241 | 0.329 | 1.14 |
| Ricker with random-walk productivity | −0.009 | 0.177 | 0.318 | 0.84 |
| Naive R1C (return 4 years earlier) | +0.051 | 0.211 | 0.379 | 1.00 |
| Naive R2C (mean of 4 and 8 years earlier) | +0.082 | 0.213 | 0.361 | 1.01 |
| Naive RAC (mean of earlier same-cycle returns) | +0.185 | 0.296 | 0.362 | 1.40 |
| Naive TSA (mean of all earlier returns) | +0.190 | 0.308 | 0.368 | 1.46 |

Units are millions of fish. The last column is the mean absolute scaled error (Hyndman and Koehler 2006 **[lit]**): below 1 beats the naive forecast.

Four things follow.

- The slide's forecast is biased high by about half the mean return and was too high in 21 of 22 years.
- It does not beat "the same as four years ago" on mean absolute error. It does on RMSE, because the naive forecast repeats the exceptional 2014 return four years later.
- The five covariates buy a small improvement here (MAE 0.223 against 0.241). Estimated with shrinkage, their effects shrink away and the forecasts match the no-covariate model.
- Letting productivity drift removes the bias and gives the lowest error.

The same contrast shows in the forecast distributions (estimation plus process uncertainty) **[run V06-F]**:

| | Stationary Ricker | Random-walk productivity | Expected |
|---|---|---|---|
| Observed inside the 50% interval | 6 of 22 | 12 of 22 | 11 |
| Observed inside the 80% interval | 14 of 22 | 15 of 22 | 17.6 |
| Observed below the 10% point | 8 | 6 | 2.2 |
| Observed below the median | 21 | 14 | 11 |

Stellako is one stock and 22 years is a small sample, so I repeated the comparison for all 19 stocks: forecast log(R/S) for each of the last 20 brood years from data two or more brood years older **[run V11-D]**.

| 379 forecasts | Share too high | Mean error (log scale) | RMSE | RMSE without Harrison | Mean log score (higher is better) |
|---|---|---|---|---|---|
| Stationary Ricker | 69% | +0.58 | 1.99 | 1.57 | −2.32 |
| Random-walk productivity | 49% | −0.03 | 2.02 | 1.48 | −2.16 |

The stationary model over-forecasts on average in 18 of 19 stocks, by a factor of about 1.8. The random walk removes that and improves the probabilistic score. Its gain in point accuracy is modest: lower RMSE in 10 of 19 stocks. Both models' 80% intervals cover only about 63% of outcomes, which fits the heavy tails found in the residuals (M6-9) and argues for a heavier-tailed error distribution.

A forecast section should therefore teach: a rolling-origin hindcast that uses only data available at the time (Module 8's loop does this), a naive benchmark, an error measure relative to it, forecast intervals with their empirical coverage, and a proper score for the distribution (Gneiting and Raftery 2007 **[lit]**).

**(f) Zeros and age structure.** The alternative to `+1e-4` (M8-5) is to stop taking logs of age-specific recruits. Model total brood-year recruits with the Ricker, and the split among ages with a composition likelihood (multinomial, Dirichlet-multinomial or logistic-normal), which takes zeros in its stride. That is the design of the standard state-space model (Fleischman et al. 2013) and, per its documents, of DFO's forecast, which models total recruits and applies age proportions or a sibling regression **[lit]**. Adding a constant before a log is a known source of arbitrary results (O'Hara and Kotze 2010 **[lit]**). Module 8's conditional age-5-given-age-4 term is a sibling relationship in disguise, and it would be worth saying so.

**(g) Environmental covariates.** The literature is cautious **[lit]**. Published environment–recruitment correlations rarely survive a retest (Myers 1998), and climate–salmon relationships in the North Pacific have changed over time (Litzow et al. 2018, 2019, 2020). DFO's own documents report that covariate models "do not rank high in terms of model performance" (Science Response 2018/034, from a search snippet I could not read in context). The Stellako result above agrees. If covariates are used: choose few on prior grounds, standardise with training-window statistics, test out of sample against a naive model, and shrink the coefficients, which in RTMB means `random = "a"` with a normal distribution on them (V06-F).

### 4.9 A Bayesian section

One slide passes an object to `tmbstan()`. A section would need: priors on every parameter (a flat prior on an unbounded log-SD is improper); when a Jacobian is needed (for sampling, when the prior is written for a transformed parameter; not for a penalised-likelihood mode, as sdmTMB's documentation notes **[doc]**); the difference between a posterior mode and a posterior mean; and convergence diagnostics. `7b_tmbstan.qmd` was evidently meant to be this (M7b-1). Section 5.6 lists the current samplers.

### 4.10 Reproducibility

- **Seeds and versions.** Random starts are unseeded and no package versions are recorded. RTMB changed materially between 1.7 and 2.0 (section 5.1).
- **Unevaluated code.** `eval: false` chunks are where most of the broken code lives (M5-2, M5-3, M7-2). The rendered HTML is committed, and for Module 5 it cannot show the time-varying results because they were never run.
- **Files.** Names with spaces and download suffixes (`SRDATA2026 1.csv`, `FC_StockInfo 2.csv`), a byte-order mark (M8-1), a duplicated example and a placeholder module (M7a-4, M7b-1).

---

## 5. More modern or sophisticated approaches

"Status" below is as of November 2026, from package sources and documentation **[doc]** unless marked otherwise.

| Area | What the workshop uses | A more current option | Status | Worth it when |
|---|---|---|---|---|
| RTMB itself | The core interface, plus `ADoverload("if")` | AD `uniroot()`, `integrate()`, `Vectorize()`; `obj$force.update()`; `dflat()`; tips vignette | In RTMB 1.8 to 2.0, on CRAN | Now. Section 5.1 |
| Fitting | One `nlminb()` call | Raised limits, Newton steps, eigenvalue and gradient checks | Standard in glmmTMB, sdmTMB, WHAM, SAM | Always. Section 4.1 |
| Variance components | ML | REML by integrating mean parameters; profile intervals | One argument in RTMB; `REML = TRUE` in glmmTMB | Few groups, or many mean parameters. Section 4.3 |
| Smoothers | Penalty with estimated weight | `mgcv::smooth2random()`, coefficients as random effects | Used by sdmTMB, glmmTMB, brms | Any smooth term. Section 5.2 |
| Time-varying productivity | AR(1) year effect, not used for forecasting | Random-walk or regime-switching Ricker; state-space run reconstruction | `samEst` (TMB and Stan; RTMB versions in progress) | Forecasting; reference points under change. Section 5.3 |
| Residuals | `oneStepGeneric`; QQ-plot | `fullGaussian` or `cdf` where valid; one-sample process residuals; composition-specific residuals | In TMB and RTMB; `compResidual` on GitHub | Always check both levels. Section 5.4 |
| Forecast evaluation | A plot | Leave-future-out hindcast, naive benchmark, scaled error, interval coverage, proper scores | Standard practice; `samEst::tmb_mod_lfo_cv` | Any forecast. Section 4.8(e) |
| Model comparison | None | Marginal and conditional AIC; leave-future-out | `cAIC()` in sdmTMB and tinyVAST | Choosing structure. Section 4.5 |
| Bayesian sampling | `tmbstan` | SparseNUTS; StanEstimators | SparseNUTS 1.0.2 and StanEstimators 0.4.0, neither on CRAN; tmbstan 1.2.1 on CRAN | Correlated hierarchical posteriors. Section 5.6 |
| Extra distributions | Hand-written, using RTMB internals | RTMBdist | On CRAN (1.2.1) | Before writing your own. Section 5.5 |
| Covariance among stocks | Hand-built two-factor model | `glmmTMB::rr()`, gllvm, dynamic factor analysis (MARSS, bayesdfa); `unstructured()` + `dmvnorm()` in RTMB for a few stocks | All on CRAN | Section 5.7 |
| Age composition | `log(x + 1e-4)` per age | Dirichlet-multinomial, logistic-normal, multivariate Tweedie | In WHAM and SAM; Dirichlet-multinomial also in RTMBdist | Any age-structured model. Section 4.8(f) |

### 5.1 RTMB has moved

The material was written against a recent RTMB (it uses `ADoverload("if")`, which arrived in 1.9), but it leaves out several additions from the last year that suit it. Current CRAN versions are RTMB 2.0 (17 September 2026) and TMB 1.9.25. From the RTMB NEWS file **[doc]**:

- **1.8 (October 2025).** `force.update()` for tapes and model objects, for use with `DataEval()` (M4-5). `dgamma(rate = )` fixed. AD versions of `sort()`, `findInterval()` and `%%`. More sparse-matrix methods.
- **1.9 (March 2026).** "AD adaptive integration via `integrate()`", "AD root finding via `uniroot()`", a fast AD `Vectorize()`, and "Experimental atomic branching via `ADoverload('if')`". A new "tips" vignette by Ben Bolker, which is the nearest thing to an official do-and-don't list. Special functions now keep the attributes of their arguments.
- **2.0 (September 2026).** `dflat()`, an experimental `obj$vectorize()`, and a fix to `dautoreg()` simulation for orders above 1.
- **A paper to cite**: Kristensen, Nielsen, Albertsen, Berg and Trijoulet (2026), *Ecology and Evolution* 16(11): e74572.

Things in the tips vignette that the course does not say **[doc]**: no `abs()`, `round()`, `min()` or `max()` of parameters; `[[` indexing can be much faster than `[`; vectorised code builds the model object much faster; `cat()` and `print()` of AD values do not show what you expect.

One open defect to know about: RTMBdist's NEWS reports that RTMB 2.0's `pbeta()` and `qbeta()` have wrong derivatives in the shape parameters when a shape is an integer or half-integer, fixed in development TMB **[doc]**.

RTMB also has tools the course does not reach that suit this material: `unstructured()` with `dmvnorm()` for a full correlation matrix, `dseparable()` for stock-by-year effects, `F$laplace(indices, SPA = TRUE)` for a saddle-point variant, and `MakeADFun(..., integrate = )` for alternatives to Laplace **[doc]**.

### 5.2 Smoothers through the mixed-model form

This replaces the penalised model in `examples/GAM.qmd` (G-1). `mgcv` builds the basis. `smooth2random()` splits it into an unpenalised part and a part whose coefficients are independent normals with one SD. That SD is the smoothing parameter, and declaring the coefficients `random` lets the Laplace approximation estimate it.

```r
sm <- mgcv::smoothCon(mgcv::s(time, bs = "tp", k = 10),
                      data = data.frame(time = 1:nyears), absorb.cons = TRUE)[[1]]
rs <- mgcv::smooth2random(sm, "", type = 2)
Xf <- rs$Xf                 # unpenalised part (here the linear trend)
Zs <- rs$rand[[1]]          # penalised part

# inside the objective function:
smooth_t <- drop(Xf %*% bf + Zs %*% bs)
mu  <- alpha_stock[stock] - beta[stock]*spawners_scaled + smooth_t[time]
nll <- nll - sum(dnorm(bs, 0, exp(logsd_s), log = TRUE))

# and when building the object:
obj <- MakeADFun(f, parameters, random = c("re", "bs"))
```

With a common residual SD this gives the same curve as `mgcv::gam(method = "ML")` to 1.5 × 10⁻⁵ and the same log-likelihood to two decimals **[run V07-B]**. sdmTMB, glmmTMB (`s()`, labelled experimental) and brms all build smooths this way **[doc]**. tinyVAST does not: it uses the penalty matrix directly with a normalising term, so check that pattern before copying it **[doc, as read by a research agent; I did not check it]**.

A note on what the smooth shows here. mgcv gives this smooth 8.6 effective degrees of freedom out of 9, so it is barely penalised. The year-to-year signal is rougher than a 10-knot spline can follow: the AR(1) year effect fits far better (marginal NLL 1650.7 against 1771.1) **[run V07-C, V03-C]**.

### 5.3 State-space and time-varying spawner–recruit models

- **`samEst`** (Pacific-salmon-assess, the organisation hosting this workshop) has stationary, autocorrelated, random-walk (`ricker_rw_TMB`), Kalman-filter and hidden-Markov regime (`ricker_hmm_TMB`) Ricker models in TMB and Stan, reference-point functions, and leave-future-out cross-validation (`tmb_mod_lfo_cv`). It exports RTMB versions of the stationary and autocorrelated models **[doc]**. The accompanying paper is Wor et al. (2025) **[lit]**. This is the natural next step from Module 5, and the workshop does not point to it.
- **Age-structured state-space models** (Fleischman et al. 2013) are the standard for escapement-goal analysis in Alaska and are used in Canada **[lit]**. Published implementations are in JAGS or Stan. The literature search found no published RTMB port, though it reports that the `salmonMSE` package conditions its operating model with RTMB **[lit, not verified by me]**.
- **`LaMa`** (CRAN) fits hidden-Markov and related models by direct likelihood maximisation on top of RTMB **[doc]**, which is the route to regime-shift Ricker models in RTMB.

### 5.4 Model validation

- **Pick the one-step-ahead method deliberately.** RTMB's help says "please carefully read the documentation, especially the `method` argument" **[doc]**. For linear Gaussian models `fullGaussian` is exact and fast. RTMB supports `method = "cdf"` automatically for built-in densities that have an AD distribution function **[doc]**; here it matched the exact residuals to 0.006 **[run V04-F]**.
- **Check the random effects with one joint draw** (section 4.6).
- **DHARMa** changed its default to conditional simulations in version 0.5.0 **[doc]**. If you use it with a random-effects model, draw the random effects once from their joint distribution, as sdmTMB's `simulate(type = "mle-mvn")` does **[doc]**.
- **Composition data** need their own residuals (Trijoulet et al. 2023 **[lit]**). The last category carries no information and should be dropped; WHAM removes it **[doc]**.

### 5.5 RTMBdist before writing a distribution

RTMBdist (CRAN 1.2.1, November 2026) adds many distributions with simulation and one-step-ahead support, a differentiable `lambertW()`, and `laplace_check()` **[doc]**. Module 7a already uses its `ddirmult()` to check the hand-written density, and its `lambertW()` would replace the `ADjoint` version in Module 6. Two cautions: its `ddirmult` does not support one-step-ahead residuals, and it defines `dnbinom2(x, mu, size)` where RTMB has `dnbinom2(x, mu, var)` **[doc]**.

### 5.6 Bayesian sampling from an RTMB object

| Package | Status | Notes |
|---|---|---|
| tmbstan | CRAN 1.2.1 (September 2026). Needs R ≥ 4.5 | The route shown in the workshop. Its README says of `laplace = TRUE`: "This is generally not recommended" **[doc]** |
| SparseNUTS | 1.0.2, GitHub only | From the adnuts author. Decorrelates and rescales the posterior before sampling (the "sparse NUTS" algorithm), and runs through StanEstimators. Preprint: Monnahan et al. (2026) **[doc]** |
| StanEstimators | 0.4.0, r-universe | Runs Stan's NUTS, Pathfinder, ADVI and Laplace on any R function. Supports `grad_fun = "RTMB"` **[doc]** |
| adnuts | CRAN 1.1.2 (2021) | Frozen. Its experimental TMB features moved to SparseNUTS **[doc]** |

### 5.7 Packages that already do what the examples hand-code

Hand-coding is the point of the course. It is still worth knowing what you would check against.

- **glmmTMB** fits the Module 5 model in one line (V14), does REML, has `ar1()` and `rr()` (reduced-rank, which is the latent-factor covariance), and `s()` smooths. It has an experimental RTMB back-end (`glmmTMBControl(use_rtmb = TRUE)`) **[doc]**.
- **sdmTMB** made RTMB its default back-end in 1.2.0 (October 2026), and its development version has removed the C++ one **[doc]**. Its source is now a large, maintained example of production RTMB code: smoothers, priors, residuals, cross-validation.
- **gllvm**, **MARSS** and **bayesdfa** fit factor and dynamic-factor models, with rotation handled for you.

### 5.8 Forecasting methods beyond the Ricker

For completeness, from the literature **[lit]**. Empirical dynamic modelling was demonstrated on nine Fraser sockeye stocks by Ye et al. (2015), who reported better forecasts than Ricker models. Machine-learning and ensemble forecasts for Bristol Bay sockeye gave gains of about 15% in some systems and losses in others (Ovando et al. 2022). Dynamic linear models with time-varying covariate effects go back to Scheuerell and Williams (2005). The search also found a DFO repository exploring EDM for pre-season forecasting (`Pacific-salmon-assess/salmon_edm`). None of these changes the first priority, which is to evaluate whatever is used against a naive benchmark.

---

## 6. What to do with this

**If you use this material for your own analyses**

1. Keep the mechanics: the closure pattern, `getAll()`, `OBS()`, `ADREPORT()`, `map`, `random`.
2. Wrap every fit in something like `fit_rtmb()` (section 4.1) and never call `obj$report()` without parameters.
3. Report variance components with profile intervals, and use REML when the number of groups is small relative to the mean parameters (section 4.3).
4. Check residuals at both levels: `fullGaussian` or `cdf` one-step-ahead residuals with an autocorrelation check, and one joint draw of the random effects (sections 4.6, 5.4).
5. Simulate from the fitted model and refit before trusting an estimator (section 4.7).
6. For forecasts: allow productivity to change, drop `+1e-4`, and score every candidate against R1C with intervals (section 4.8).
7. For smooths and penalties: never estimate a penalty weight as a free parameter (section 5.2).
8. Look at `samEst` before writing a time-varying Ricker from scratch (section 5.3).

**If you want to send fixes upstream**, these are small, self-contained and uncontroversial, in order of value:

1. `2_likelihoods.qmd`: the double sign flip (M2-1).
2. `7a_custom_distribution.qmd`: `rgamma(shape = alpha)` and `log = TRUE` (M7a-1, M7a-2).
3. `5_random_effects.qmd:521`: `ll <- dgmrf(...)` (M5-2).
4. `4_inference.qmd`: `beta / sd_s`, and the multi-start loop (M4-1, M4-2).
5. `8_practical.qmd`: `fileEncoding = "UTF-8-BOM"`, no intercept in the covariate matrix (M8-1, M8-3).
6. `latent_variable_model.qmd:152–153`: the loadings (L-2).
7. `osa_multinomial.qmd:93`: `matrix()` for `as.matrix()` (O-1).
8. `6_custom_functions.qmd`: test for `simref` first (M6-1).
9. `GAM.qmd`: the `smooth2random` version (G-1).
10. `mixture_model.qmd:27`: the claim about direct maximisation (X-1).

I have not changed any workshop file. Everything I added is under `review/`.

---

## 7. Limits of this review

- **What was run, and with what.** R 4.3.3, RTMB 2.0 and TMB 1.9.25. This environment could not reach CRAN (`cloud.r-project.org` is blocked by its network policy), so I built RTMB and TMB from their public GitHub sources (commits `3d264fb` and `0008857`, both 12 November 2026) and took other packages from Ubuntu (dplyr, mgcv, mixtools, glmmTMB 1.1.8 with TMB 1.9.10). Results could differ slightly on your versions.
- **What was not run.** `tmbstan` (not installable here), the first, B-spline half of `GAM.qmd` (`splines2`), `optimx`, `DHARMa`, `RTMBdist`, and all plotting code. I read the rendered HTML only for its printed output.
- **Literature.** I delegated the literature and documentation search to three parallel research agents and then spot-checked what they returned. Journal and agency sites were blocked, so papers were confirmed from abstracts, index records or search snippets, and quotations from papers and DFO documents come from snippets that nobody could read in context. The web-search budget ran out before DFO's 2021 to 2023 forecast documents could be checked. Package claims are firmer: they come from source code and help files, and I re-checked the ones this review leans on against the cloned repositories.
- **Details that could not be confirmed, and that nothing above relies on.** How DFO's ranking of forecast models is cross-validated (leave-one-out or leave-future-out); whether DFO's forecast method changed after 2019; whether Zheng et al. (2024) has appeared in a journal; volume and page details of some references marked †.
- **Statistical limits of my own checks.** The Stellako hindcast is 22 forecasts of one stock. The four-model comparison is in-sample AIC with unconstrained lag coefficients. Differences between methods with similar scores are within noise. The clear results are the direction and size of the bias, which the 19-stock check confirms.
- **The data.** I can show that `recruits` and the age columns disagree for brood year 2018. I cannot say which is right.
- **The live course.** I saw files, not the sessions.
- **My own code.** The verification scripts are mine and may contain mistakes. They are in the repository so you can check them. Two of my first attempts were wrong, and the scripts keep the evidence: a fixed-node quadrature reference in V12 (its self-test line shows the error) and parameter-dependent `uniroot()` bounds in V11 (section 4.8a).

---

## 8. References

Entries marked † were confirmed at abstract or snippet level only, or have a bibliographic detail that could not be confirmed.

**Software and documentation (read directly)**

- Kristensen K, Nielsen A, Albertsen CM, Berg CW, Trijoulet V (2026). RTMB: An R package for fitting latent variable models with Template Model Builder. *Ecology and Evolution* 16(11): e74572. https://doi.org/10.1002/ece3.74572 (citation as given in the RTMB package)
- Kristensen K, Nielsen A, Berg CW, Skaug H, Bell BM (2016). TMB: automatic differentiation and Laplace approximation. *Journal of Statistical Software* 70(5). https://doi.org/10.18637/jss.v070.i05
- RTMB 2.0: `NEWS`, `NAMESPACE`, help pages `ADoverload`, `Tape`, `Simulation`, `OSA-residuals`, `TMB-interface`, `Distributions`, `MVgauss`, and the "RTMB tips" vignette. https://github.com/kaskr/RTMB
- TMB 1.9.25: help pages `MakeADFun`, `oneStepPredict`, `checkConsistency`, `sdreport`; `R/checker.R`, `R/sdreport.R`; the TMB book chapters on sparsity, simulation and validation. https://github.com/kaskr/adcomp
- Package sources: glmmTMB, sdmTMB, tinyVAST, WHAM, SAM (`stockassessment`), TMBhelper, RTMBdist, LaMa, SparseNUTS, StanEstimators, tmbstan, DHARMa, mgcv, FIMS, samEst and samSim (https://github.com/Pacific-salmon-assess), the Stan User's Guide, and R's `factanal` and `nlminb` help.

**Statistics**

- Bürkner P-C, Gabry J, Vehtari A (2020). Approximate leave-future-out cross-validation for Bayesian time series models. *J Stat Comput Simul* 90:2499–2523. https://doi.org/10.1080/00949655.2020.1783262 †
- Chung Y, Rabe-Hesketh S, Dorie V, Gelman A, Liu J (2013). A nondegenerate penalized likelihood estimator for variance parameters in multilevel models. *Psychometrika* 78:685–709. https://doi.org/10.1007/s11336-013-9328-2 †
- Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction, and estimation. *JASA* 102:359–378. https://doi.org/10.1198/016214506000001437 †
- Hyndman RJ, Koehler AB (2006). Another look at measures of forecast accuracy. *Int J Forecast* 22:679–688. https://doi.org/10.1016/j.ijforecast.2006.03.001 †
- Joe H (2008). Accuracy of Laplace approximation for discrete response mixed models. *Comput Stat Data Anal* 52:5066–5074. https://doi.org/10.1016/j.csda.2008.05.002 †
- Monnahan CC, Kristensen K (2018). No-U-turn sampling for fast Bayesian inference in ADMB and TMB: introducing the adnuts and tmbstan R packages. *PLoS ONE* 13:e0197954. https://doi.org/10.1371/journal.pone.0197954 †
- Monnahan CC, Thorson JT, Kristensen K, Carpenter B (2026). Leveraging sparsity to improve no-U-turn sampling efficiency for hierarchical Bayesian models. arXiv:2603.02437 †
- Ogden HE (2021). On the error in Laplace approximations of high-dimensional integrals. *Stat* 10:e380. https://doi.org/10.1002/sta4.380 †
- O'Hara RB, Kotze DJ (2010). Do not log-transform count data. *Methods Ecol Evol* 1:118–122. https://doi.org/10.1111/j.2041-210X.2010.00021.x †
- Self SG, Liang K-Y (1987). Asymptotic properties of maximum likelihood estimators and likelihood ratio tests under nonstandard conditions. *JASA* 82:605–610. †
- Stram DO, Lee JW (1994). Variance components testing in the longitudinal mixed effects model. *Biometrics* 50:1171–1177. †
- Thorson JT, Kristensen K (2016). Implementing a generic method for bias correction in statistical models using random effects. *Fish Res* 175:66–74. https://doi.org/10.1016/j.fishres.2015.11.016 †
- Thygesen UH, Albertsen CM, Berg CW, Kristensen K, Nielsen A (2017). Validation of ecological state space models using the Laplace approximation. *Environ Ecol Stat* 24:317–339. https://doi.org/10.1007/s10651-017-0372-4 †
- Trijoulet V, Albertsen CM, Kristensen K, Legault CM, Miller TJ, Nielsen A (2023). Model validation for compositional data in stock assessment models: calculating residuals with correct properties. *Fish Res* 257:106487. https://doi.org/10.1016/j.fishres.2022.106487 †
- Waagepetersen R (2006). A simulation-based goodness-of-fit test for random effects in generalized linear mixed models. *Scand J Stat* 33:721–731. https://doi.org/10.1111/j.1467-9469.2006.00504.x †
- Wood SN (2011). Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. *JRSS B* 73:3–36. https://doi.org/10.1111/j.1467-9868.2010.00749.x †
- Zheng N, Cadigan N, Thorson JT (2024). A note on numerical evaluation of conditional Akaike information for nonlinear mixed-effects models. arXiv:2411.14185 †

**Salmon stock–recruit analysis and forecasting**

- Fleischman SJ, Catalano MJ, Clark RA, Bernard DR (2013). An age-structured state-space stock–recruit model for Pacific salmon. *CJFAS* 70:401–414. https://doi.org/10.1139/cjfas-2012-0112 †
- Holt CA, Connors B, Greenberg D, Wor C (2025). *Can. Tech. Rep. Fish. Aquat. Sci.* 3653. https://doi.org/10.60825/002c-as40 †
- Holt CA, Michielsens CGJ (2020). Impact of time-varying productivity on estimated stock–recruitment parameters and biological reference points. *CJFAS* 77:836–847. https://doi.org/10.1139/cjfas-2019-0104 †
- Kehler DG, Myers RA, Field CA (2002). Measurement error and bias in the maximum reproductive rate for the Ricker model. *CJFAS* 59:854–864. †
- Korman J, Peterman RM, Walters CJ (1995). Empirical and theoretical analyses of correction of time-series bias in stock–recruitment relationships of sockeye salmon. *CJFAS* 52:2174–2189. https://doi.org/10.1139/f95-810 †
- Larkin PA (1971). Simulation studies of the Adams River sockeye salmon. *J Fish Res Board Can* 28:1493–1502. https://doi.org/10.1139/f71-230 †
- Litzow MA et al. (2018) *Proc R Soc B* 285:20181855; (2019) *Ecology* 100:e02760; (2020) *PNAS* 117:7665–7671. †
- Myers RA (1998). When do environment–recruitment correlations work? *Rev Fish Biol Fish* 8:285–305. https://doi.org/10.1023/A:1008828730759 †
- Ovando D, Cunningham C, Kuriyama P, Boatright C, Hilborn R (2022). Improving forecasts of sockeye salmon with parametric and nonparametric models. *CJFAS* 79:1198–1210. https://doi.org/10.1139/cjfas-2021-0287 †
- Peterman RM, Dorner B (2012). A widespread decrease in productivity of sockeye salmon populations in western North America. *CJFAS* 69:1255–1260. https://doi.org/10.1139/f2012-063 †
- Peterman RM, Pyper BJ, MacGregor BW (2003). Use of the Kalman filter to reconstruct historical trends in productivity of Bristol Bay sockeye salmon. *CJFAS* 60:809–824. https://doi.org/10.1139/f03-069 †
- Scheuerell MD (2016). An explicit solution for calculating optimum spawning stock size from Ricker's stock recruitment model. *PeerJ* 4:e1623. https://doi.org/10.7717/peerj.1623 †
- Scheuerell MD, Williams JG (2005). Forecasting climate-induced changes in the survival of Snake River spring/summer Chinook salmon. *Fish Oceanogr* 14:448–457. †
- Staton BA, Catalano MJ, Connors BM, et al. (2020). Evaluation of methods for spawner–recruit analysis in mixed-stock Pacific salmon fisheries. *CJFAS* 77:1149–1162. https://doi.org/10.1139/cjfas-2019-0281 †
- Walters CJ (1985). Bias in the estimation of functional relationships from time series data. *CJFAS* 42:147–149. †
- Walters CJ, Ludwig D (1981). Effects of measurement errors on the assessment of stock–recruitment relationships. *CJFAS* 38:704–710. https://doi.org/10.1139/f81-093 †
- Wor C, Greenberg D, Holt C, Connors B, Feddern M, Freshwater C, Britten G, Mazur M (2025). Recommendations for estimating and detecting time-varying spawner–recruit dynamics in fish populations. *Ecological Modelling* 507:111159. †
- Ye H, Beamish RJ, Glaser SM, Grant SCH, Hsieh C-H, Richards LJ, Schnute JT, Sugihara G (2015). Equation-free mechanistic ecosystem forecasting using empirical dynamic modeling. *PNAS* 112:E1569–E1576. https://doi.org/10.1073/pnas.1417063112 †
- DFO pre-season Fraser sockeye forecast documents, including Grant et al. (2010) CSAS Res. Doc. 2010/042, MacDonald and Grant (2012) CSAS Res. Doc. 2012/011, and Science Response 2018/034. †

---

## Appendix: verification scripts

All in `review/verification/`. Run from the project root, for example `Rscript review/verification/06_module8_forecast.R`, or all of them with `Rscript review/verification/run_all.R` (about 20 minutes; scripts 03 and 04 take most of it). Saved console output from my runs is in `review/verification/output/`.

They need RTMB 1.9 or later, TMB, dplyr, here, MASS, Matrix, mgcv and mixtools. Script 14 needs glmmTMB instead of RTMB.

| Script | Covers | Findings it supports |
|---|---|---|
| `00_setup.R` | Shared data preparation, mirroring the slides' hidden setup chunk | — |
| `01_module2_likelihood.R` | Sign of `negloglik`; what the "looped" tape contains | M2-1, M2-2 |
| `02_module4_inference.R` | `beta_real`; starting values; multi-start; profile interval; `DataEval` bootstrap; Hessian check | M3-1, M4-1 to M4-5 |
| `03_module5_random_effects.R` | Joint likelihood; GMRF slide; three AR(1) codings; `last.par`; ML and REML; `checkConsistency`; sparsity | M5-1 to M5-5, section 4.3 |
| `04_module6_custom_and_residuals.R` | `dlognorm` and global flags; `dlgamma`; `dlnorm`; Jacobians; one-step-ahead methods, ordering and what the residuals show | M6-1 to M6-9 |
| `05_module7_misc_and_custom_dist.R` | `make_adfun`; `ad_context` example; Dirichlet sampler; missing `log = TRUE`; `%~%` | M7-1 to M7-5, M7a-1 to M7a-3 |
| `06_module8_forecast.R` | Byte-order mark; `recruits` consistency; covariate model; the `+1e-4` constant; median and mean; hindcast skill and calibration | M8-1 to M8-7, section 4.8(e) |
| `07_example_gam.R` | Penalty weight; `smooth2random` against mgcv | G-1, G-2, section 5.2 |
| `08_example_latent_variable.R` | Convergence; loadings bug; shared year effect; penalty; rotation | L-1 to L-5 |
| `09_example_mixture.R` | Starting values; direct optimisation against EM; unbounded likelihood; log-sum-exp | X-1 to X-3 |
| `10_example_osa_multinomial.R` | Residual reshaping; last-category residuals | O-1 |
| `11_missing_pieces.R` | Reference points; bias correction; Ricker and Larkin with and without random-walk productivity; leave-future-out for 19 stocks | Sections 4.4, 4.8 |
| `12_laplace_accuracy.R` | Laplace against exact likelihood for binary data | M5-4, section 4.2 |
| `13_recommended_workflow.R` | `fit_rtmb()`; one-sample process residuals; profile and simulate-and-refit for the between-stock SD | Sections 4.1, 4.3, 4.6, 4.7 |
| `14_crosscheck_glmmTMB.R` | The Module 5 model in glmmTMB, ML and REML | P-3, section 4.3 |

