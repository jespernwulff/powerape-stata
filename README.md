# powerape for Stata

**Power analysis for average partial effects (APEs) and average interaction
effects (AIEs) of probit and logit models.**

**Documentation, a getting-started guide and workflows, in Stata and R:
<https://jespernwulff.github.io/powerape/>**

`powerape` is the Stata command of the R package of the same name
(<https://github.com/jespernwulff/powerape>). In a probit or logit model the
coefficient is not the quantity anyone interprets: conclusions are stated in
probability points, as the APE, or, for a moderation question, as the AIE.
`powerape` takes the effect in those units, simulates the whole study (draw
the data, fit the model, compute the effect with a delta-method confidence
interval) and reports power for three claims on that interval (Riesthuis,
2024): that the effect exceeds the smallest effect size of interest (SESOI),
that it differs from zero, or that it is too small to matter.

## Installation

```stata
net install powerape, from("https://raw.githubusercontent.com/jespernwulff/powerape-stata/main/")
mata: mata mlib index
powerape version
```

`mata: mata mlib index` makes the newly installed Mata library known to the
running session (needed once, right after installing). `powerape version`
checks that the command and its library belong together. To update, add the
`replace` option to `net install`.

`net get powerape, from(...)` with the same address copies two files into the
current folder: `showcase.do`, a commented tour of every subcommand, and
`build.do`, which recompiles the Mata library.

## Quick start

```stata
* 1. Describe the world: a baseline outcome rate and a covariate signal
powerape dgp d, focal(treat, binary p(0.5))                            ///
    covariates((age, normal mean(45) sd(12)) (female, binary p(0.55))) ///
    correlation(0.2) baseline(0.30) signal(0.15)

* 2. Pin the true APE at 10 percentage points
powerape setape d, target(0.10)

* 3. Power for showing the effect exceeds a SESOI of 5 points, at n = 1,500
powerape power d, n(1500) claim(minimum) sesoi(0.05) seed(1)

* 4. Required n, the minimum detectable effect, robustness, the statement
powerape n d, power(0.80) claim(minimum) sesoi(0.05) seed(2) engine(mata)
powerape mde d, n(1500) claim(minimum) sesoi(0.05) seed(3) engine(mata)
powerape robust d, n(2200) claim(minimum) sesoi(0.05) vary(baseline(0.20 0.40)) engine(mata)
powerape statement
```

The site walks through these steps with their output:
[Get started in Stata](https://jespernwulff.github.io/powerape/articles/stata.html).

## What powerape does

* **Worlds.** Parametric covariates with a Gaussian copula; covariates
  resampled from pilot data in memory (`empirical`); the world of a fitted
  `probit` or `logit` pilot model (`fromfit`); correlated random effects
  panel probit with Mundlak means, unit-clustered standard errors and
  unbalanced panels (`panel`); an endogenous treatment estimated by control
  functions, the estimator of `cfprobit` (`iv`).
* **Effects.** APEs of binary and continuous variables and AIEs (the
  estimand of `ginteff`), pinned in probability points by `setape` and
  `setaie`; targets that cannot exist are refused with the attainable range.
* **Power and searches.** Power with the distribution of study outcomes,
  power curves, required sample sizes and minimum detectable effects with a
  high-precision confirmation stage, and robustness sweeps with n_max
  (Hancock & Feng, 2026).
* **Error rates and reporting.** Alpha as a sample-size rule or as the
  error-cost optimum (Maier & Lakens, 2022), from one stored simulation;
  a citable methods paragraph for any result.

## Engines

Every simulated study is fitted, by default, by Stata's own `probit`,
`logit` or `cfprobit` (`engine(stata)`), with the effects computed from
`e(b)` and `e(V)` as `margins` computes them; `effects(margins)` computes
them with `margins` or `ginteff` themselves. `engine(mata)` fits each study
with a Mata port of the R package's routine and is 10 to 30 times faster:
use it for searches, curves and sweeps, and confirm the final design with
`engine(stata)`. Every simulation command prints a `fit:` line that names
the engine.

## Requirements

Stata 16 or later. The control-function route under the default engine uses
`cfprobit` (Stata 18 or later). `effects(margins)` for interaction effects
uses `ginteff` (Radean, 2023):
`net install st0711, from(http://www.stata-journal.com/software/sj23-2)`.
Tested on StataNow 19.5. The Mata library `lpowerape.mlib` was compiled by
StataNow 19.5; on Stata 16 to 18, run `net get powerape` and `do build.do`
to recompile it, then `mata: mata mlib index`.

## Validation

The command is checked against the R package in four tiers (exact power by
enumeration, the estimators on fixed data, the calibration constants, and
simulated power, sample sizes and detectable effects in matched designs) and
against Stata's `probit`, `logit`, `margins`, `ginteff` and `cfprobit`. See
[validation/RESULTS.md](validation/RESULTS.md) and
[How powerape is validated](https://jespernwulff.github.io/powerape/articles/validation.html).

Everything needed to rerun it is in [validation/](validation/): one script
per tier (`v_tierA.do` to `v_tierD.do`) with its results
(`tier?-results.txt`) and log (`v_tier?.log`), the fixed data sets and R
reference values in `validation/data/`, and the R scripts that produced them
(`export_fixtures.R` and `tierD_r.R`, with the R package powerape 1.11.3).
Start Stata in the root of a clone of this repository and run a tier:

```stata
cd "path/to/powerape-stata"
do validation/v_tierA.do
```

The scripts put the clone's `ado/` folder first on the adopath, so they check
the command in the clone rather than an installed copy, and they stop if
started from another folder. Tier B takes a minute or two, Tier C about
three, Tier A about half an hour and Tier D 20 to 40 minutes. Tiers B and D
also use `cfprobit` and `ginteff` (see Requirements); the scripts declare
`version 19` and were run on StataNow 19.5. The R scripts and the regression
tests run from the root too: `Rscript validation/tierD_r.R`,
`do tests/run-tests.do`.

## Version history

* 0.1.3 (2026-10-17): panel designs calibrate about as fast as in the R
  package. A panel `powerape dgp`, `setape`, `setaie`, `trueape` or
  `trueaie` at the default `nint(1e5)` took about a minute and now takes one
  to three seconds, so robustness sweeps over panel scenarios and `mde`
  searches, which recalibrate at every step, take seconds instead of
  minutes. Every result is bit-identical to 0.1.2.
* 0.1.2 (2026-10-16): every simulation command prints a `fit:` line naming
  who fitted the simulated studies and stores it in `r(fit)`, `r(engine)`
  and `r(effects)`; `engine()` and `effects()` given on one command apply to
  that command only.
* 0.1.1 (2026-10-16): `net install`, `powerape version`, Stata 16 or later,
  `separation(keep)` under `engine(stata)`, `effects(margins)` for the AIE
  on the control-function route.
* 0.1.0 (2026-10-15): first release, a port of the R package powerape
  1.11.3.

## References

Hancock, G. R., & Feng, Y. (2026). nmax and the quest to restore caution,
integrity, and practicality to the sample size planning process.
*Psychological Methods*, *31*(1), 174--199. <https://doi.org/10.1037/met0000776>

Maier, M., & Lakens, D. (2022). Justify your alpha: A primer on two
practical approaches. *Advances in Methods and Practices in Psychological
Science*, 5(2). <https://doi.org/10.1177/25152459221080396>

Radean, M. (2023). ginteff: A generalized command for computing interaction
effects. *The Stata Journal*, 23(2), 301--335.
<https://doi.org/10.1177/1536867X231175253>

Riesthuis, P. (2024). Simulation-based power analyses for the smallest
effect size of interest: A confidence-interval approach for minimum-effect
and equivalence testing. *Advances in Methods and Practices in Psychological
Science*, 7(2). <https://doi.org/10.1177/25152459241240722>

## Author

Jesper N. Wulff, Aarhus University. Report bugs at
<https://github.com/jespernwulff/powerape-stata/issues>.

## License

MIT; see [LICENSE](LICENSE).
