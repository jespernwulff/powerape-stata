# powerape for Stata: validation against the R package

Generated 2026-10-15 from `v_tierA.do`, `v_tierB.do`, `v_tierC.do`, `v_tierD.do`
(StataNow/SE 19.5, update 14 Oct 2026) against powerape 1.11.3 (R 4.5.2); all four
tiers regenerated 2026-10-16 under version 0.1.1 (sections "2026-10-16 additions" and
"2026-10-16 Tier C/D rerun"); version 0.1.2 (the evening of 2026-10-16) changed only what
the commands print and store; Tier A, the one tier with a record gap from the engine leak,
was rerun under it the same night (12/12 PASS, result file byte-identical; section
"2026-10-16 Tier A rerun"). On 2026-10-17 the scripts became runnable from a clone of the
repository (start Stata in its root, then `do validation/v_tierA.do`) and all four tiers
were rerun under version 0.1.3, every result file byte-identical (last two sections).
Reference values come from `export_fixtures.R` (fixed datasets and R internals) and
`tierD_r.R` (matched simulations), also run from the repository root. Raw result lines:
`tier?-results.txt`; logs `v_tier?.log`.

**Updated 2026-10-16 (version 0.1.1):** all four tiers were regenerated after the
engine changes of 0.1.1 (the Stata engine's `separation(keep)`, `effects(margins)`
for the IV AIE, `version 16`). Tiers A and B gained rows; their pre-existing rows are
bit-identical to the 2026-10-15 run. Tiers C and D, rerun later the same day, are
byte-identical to the 2026-10-15 result files (30/30 and 45/45 PASS). Details in the
sections "2026-10-16 additions (version 0.1.1)" and "2026-10-16 Tier C/D rerun under
0.1.1" below; Tier A was rerun once more under 0.1.2 (section "2026-10-16 Tier A
rerun under 0.1.2").

Two engines are validated throughout. `engine(stata)` (the default) fits every
simulated study with Stata's own `probit`/`logit` (`vce(robust)`, `vce(cluster)`,
`asis`) or `cfprobit`, and computes the APE/AIE by closed-form delta-method
Jacobians from e(b)/e(V) (or by `margins`/`ginteff` with `effects(margins)`).
`engine(mata)` is the self-contained Mata port of the R engine (Fisher scoring,
expected information, HC0 and G/(G-1) sandwiches, stacked control-function vcov).

## Tier A: exact finite-sample power by binomial enumeration (12/12 PASS)

The saturated design (binary focal, no covariates) has an exact power for the
engine's decision rule (R battery V1; `exact_power_sat` ported to Mata). Criterion:
|simulated - exact| <= 4 MCSE + 1e-4.

| design | claim | engine | nsim | exact | simulated | MCSE | diff |
|---|---|---|---|---|---|---|---|
| .30 -> .40, n = 712 | detect | mata | 40,000 | 0.8018 | 0.8059 | 0.0020 | 0.0041 |
| .30 -> .40, n = 712 | minimum (sesoi .04) | mata | 40,000 | 0.5175 | 0.5163 | 0.0025 | 0.0012 |
| truth 0, n = 712 | equivalence (sesoi .05) | mata | 40,000 | 0.0000 | 0.0000 | 0 | 0 |
| truth 0, n = 712 | equivalence (sesoi .07) | mata | 40,000 | 0.3060 | 0.3068 | 0.0023 | 0.0008 |
| .30 -> .40, n = 712 | detect | stata (probit) | 4,000 | 0.8018 | 0.8030 | 0.0063 | 0.0012 |
| .30 -> .40, n = 712 | minimum (sesoi .04) | stata (probit) | 4,000 | 0.5175 | 0.5102 | 0.0079 | 0.0073 |
| .30 -> .40, n = 712 | detect | stata + margins | 1,000 | 0.8018 | 0.8040 | 0.0126 | 0.0022 |
| .05 -> .01, n = 300, separation(fail) | detect | mata | 20,000 | 0.3503 | 0.3496 | 0.0034 | 0.0007 |
| .05 -> .01, n = 300, separation(keep) | detect | mata | 20,000 | 0.5602 | 0.5619 | 0.0035 | 0.0017 |
| .05 -> .01, n = 300, separation(fail) | detect | stata (probit) | 4,000 | 0.3503 | 0.3490 | 0.0075 | 0.0013 |
| .05 -> .01, n = 300, separation(keep) (0.1.1; rerun 0.1.2) | detect | stata (probit, asis) | 4,000 | 0.5602 | 0.5705 | 0.0078 | 0.0103 |
| probit = logit on the same seed, n = 500 | detect | mata | 2,000 | 0.6640 | 0.6640 | - | 0 |

The rare-outcome rows certify both separation conventions under both engines: `fail`
(Stata's default `probit` behaviour, the R default since 1.11.0) and `keep` (R's glm
scoring). The Mata engine reproduces R's Fisher-scoring fit; since 0.1.1 the Stata
engine keeps the `asis` fit's separated coefficient together with its zero variance
row, so the delta-method interval collapses onto the identified cells, the same limit.
On one seed the two engines' `keep` powers coincide exactly (.555 and .555 at nsim
400, no failed fits, separated share .215; `fail` gives .345 with 86 failed
replications under both engines; values from `tests/numbers_probe_2026-10-16.do` on
the run-tests seeds, which run-tests.do asserts with `nodisplay`).

Tier A was run three times on 2026-10-15/16: the first run (2026-10-15) had eleven rows,
the 0.1.1 run (2026-10-16 09:53) added the `keep` row under `engine(stata)`, and the
0.1.2 run (2026-10-16 20:36-21:07, section "2026-10-16 Tier A rerun under 0.1.2") repeated
the whole file under the leak-fixed ado-file. The 09:53 run had a labelling gap: the
`effects(margins)` row (seed 7) precedes the two rare-outcome `engine(stata)` rows (seed
10), and on the ado-file of that morning the `margins` setting leaked into them, so the
`keep` row of 0.5705 ran through `pape_effects_margins` with `asis` and not through
`pape_fit_stata`'s kept-fit logic that its label names. The 0.1.2 rerun exercises the
path as labelled (shown by `tests/tierA_path_probe_2026-10-16.do`) and returns the
same 0.5705, so the table is unchanged and now certifies the kept-fit logic.

## Tier B: fixed-data parity (171/171 PASS)

Eight fixed datasets exported from R (`data/fix?.csv`): A probit binary focal + 2
covariates; B logit continuous focal; C-F AIE for the four type pairs
(binary x binary, continuous x binary, binary x continuous under logit, continuous x
continuous); G clustered probit with a Mundlak mean (100 units x 5); H the
control-function fixture of `test-cf.R` (linear and probit first stages, cluster
SEs, the AIE with `interact(m)` and instrument-by-moderator first stage; 0.1.1 adds
the continuous-focal AIE, `linear interact(m)`, and two separated saturated designs).

Against the R internals (`fit_index_model`, `ape_*_est`, `aie_est`,
`vcov_from_glmfit`, `vcov_sandwich`, `vcov_cluster`, `cf_fit`, `cf_ape_est`,
`cf_aie_est`): all coefficients, Jacobian elements, estimates and standard errors
agree to 1e-6 or better (typical |diff| 1e-9 to 1e-8; coefficients to ML
convergence tolerance).

Against Stata's own commands on the same data:

| quantity | Stata command | agreement |
|---|---|---|
| APE, AIE point estimates | `probit`/`logit` + `margins, dydx()`, `ginteff, dydxs()` | 1e-7 to 6e-7 (ML tolerance) |
| model SEs, logit | `margins` / `ginteff` (OIM = expected information) | 1e-9 |
| model SEs, probit | `margins` / `ginteff` (OIM vs expected information) | 2e-6 to 6e-5 (documented convention difference) |
| robust SEs | `vce(robust)` (Stata's N/(N-1) factor, OIM bread) | 1e-6 to 3e-6 |
| cluster SE | `probit, vce(cluster id)` + `margins` | 2e-5 (bread convention) |
| CF coefficients and stacked SEs | `cfprobit, vce(robust)` / `vce(cluster)` | 2e-8 to 9e-7 |
| CF APE and SE | `cfprobit` + `margins, dydx()` | 9e-9 / 2e-9 (linear), 5e-7 / 7e-8 (probit first stage) |
| CF AIE via the engine's cfprobit mapping (`interact()`, `mainonly()`) | `cfprobit` | 1e-5 |
| IV AIE by `margins, dydx(d) at(m=(0 1))` after `cfprobit` (0.1.1) | R `cf_aie_est` on the same data | 4e-8 (binary focal), 1e-8 (continuous); SEs 1.3e-7, 4e-9 |
| `separation(keep)` fit under the Stata engine (0.1.1): `probit`/`logit, asis` on a separated saturated design | Mata engine (R glm) and the other arm's binomial SE | APE 1e-5 probit (the `asis` intercept stops 1e-4 short of the closed form), 1e-8 logit; SE 2e-6; with a covariate APE 5e-6, SE 1e-5 |
| closed-form helpers (req_delta, req_margin, alpha_frontier, mz_r2) | R functions | 1e-9 |

## Tier C: calibration and inversion constants (30/30 PASS)

Closed-form worlds (no covariates; binary x binary AIE) agree to 1e-7 or better.
Worlds calibrated on a Monte Carlo integration draw differ by integration noise,
since the two packages draw their own sets: at n_int = 4e5 the nuisance
coefficients and intercepts differ by 0.02-0.3% and the pinned focal coefficients by
<0.03%; at n_int = 4e6 the differences shrink to 0.004-0.03%. Panel (n_int 1e5
units) and IV (1e5) constants differ by 0.003-0.3%; the realized true effects equal
the targets to 1e-7 on each package's own draw. The empirical route (shared pilot
file, Tier D7) is exact: beta0, gamma and beta_focal agree to 1e-15.

## Tier D: stochastic agreement on matched designs (45/45 PASS)

R: 4,000 replications per power (`tierD_r.R`); Stata: 4,000 with `engine(mata)`
and 1,000 (400 for `cfprobit`, 300 for `ginteff`, 200 for `cfprobit` + `margins`)
with `engine(stata)`. Criterion: |difference| <= 3 sqrt(MCSE_R^2 + MCSE_Stata^2).

| design | R | Stata mata (4,000) | Stata engine(stata) |
|---|---|---|---|
| D1 Design A world, minimum claim, n 2000 | 0.4415 | 0.4392 | 0.4490 (probit, 1,000) |
| D1 outcome shares detect-only / inconclusive | .262 / .296 | .270 / .291 | - |
| D1 detect claim | 0.7015 | 0.7167 | - |
| D2 logit, continuous focal, robust SE, n 800 | 0.4920 | 0.4978 | 0.4880 (vce(robust)) |
| D2 same, model SE | 0.4870 | 0.4988 | - |
| D3 AIE binary x binary, n 4000 | 0.2537 | 0.2587 | 0.2450 (probit); 0.2433 (ginteff, 300) |
| D3b AIE continuous x binary, logit, n 3000 | 0.4395 | 0.4190 | - |
| D4 equivalence at truth 0, sesoi .05, n 2000 | 0.5753 | 0.5803 | - |
| D5 panel CRE probit, 150 units x 4 | 0.5430 | 0.5423 | 0.5370 (vce(cluster)) |
| D5u unbalanced panel, retention .8 | 0.4883 | 0.4828 | - |
| D6 IV, binary focal (probit first stage), n 1500 | 0.2622 | 0.2670 | 0.2750 (cfprobit, 400) |
| D6c IV, continuous focal (linear first stage), n 900 | 0.7857 | 0.7915 | 0.7850 (cfprobit); 0.7800 (+ margins, 200) |
| D7 empirical route (shared pilot), minimum claim | 0.3623 | 0.3588 | - |
| D8 required n, detect 80%, Design A world | 2,494 | 2,495 (recorded as 2,493 until 2026-10-16; see the final section) | - |
| D8 required n, minimum 80% | 5,451 | 5,461 | - |
| D8 required n, saturated .30 -> .40 (analytic 712) | 707 | 726 | 783 (nsim 600; confirmation rule holds) |
| D9 MDE at n 712, detect 80% | 0.0996 | 0.0996 | - |
| D9 equivalence margin at n 2000 | 0.0600 | 0.0599 | - |
| D10 power under alpha rule 0.05 sqrt(500/n) at n 712 | 0.7802 (alpha .0419) | 0.7825 (alpha .0419) | - |
| D10 robustness sweep, baseline .2/.3/.4 | .393/.321/.287 | .387/.329/.272 | - |
| D11 power at alpha .01/.025/.05/.10 (stored draws) | .198/.323/.443/.582 | .199/.324/.446/.584 | - |
| D11 justified alpha (4:1), flat objective | 0.0485 (wcer .1510) | 0.0490 (wcer .1509) | - |

All power comparisons pass; every required-n answer satisfies the confirmation rule
(confirmed power >= goal - max(.005, 1.5 MCSE)) and lies within 1% of R's for the
Mata engine (2.7% on the saturated design, where R's own answer sits 0.7% below the
analytic 712). The Stata-engine search at 600 replications overshoots by 11% (783),
the documented behaviour of the search at low replication counts (overshoot is never
trimmed); it passes the confirmation rule and the 15% band used for that row.

## Regression tests and smoke tests

`tests/run-tests.do` (94 checks, all PASS on 2026-10-16 under 0.1.2; the authoritative log is `tests/run-tests.log`): assert-based checks of calibration, inversion, feasibility
errors, claim coherence guards, the outcome table, both engines on the same seed,
separation conventions, RNG-state restoration, curve monotonicity, required n near
the analytic anchor, MDE direction, equivalence margins, justified alpha, the
frontier closed forms, the empirical, pilot-model, panel and IV routes, and the
robustness sweep (factorial grid, infeasible scenarios, coefficient pinning, n_max).
Version 0.1.1 added four checks: the Stata engine's `separation(keep)` reproduces the
Mata engine's keep power on the same seed with no failed fits; `effects(margins)` for
the IV AIE equals the closed form on the same seed (binary focal, .0667 = .0667 at
nsim 30; a continuous-focal pair gives .1333 = .1333); `powerape version` finds the
ado-file and the Mata library in agreement. The engine-leak fix of 2026-10-16 (ado-file
only; next-to-last section) added two: a command and a sweep with a per-command
`engine(stata)` leave a Mata session default intact (79 checks). Version 0.1.2 added 15:
`r(fit)`, `r(engine)` and `r(effects)` on `power` under both engines, with `vce(robust)`,
`asis` (keep), `vce(cluster)` (panel), `cfprobit` with `effects(closed)` and
`effects(margins)`, on `curve`, `n`, `mde` and `robust`, the per-command value against the
session default, and the printed `fit:` line read back from a named log. The checks run with `nodisplay`; the
quoted values come from `tests/numbers_probe_2026-10-16.do` on the same seeds. The
other probes behind this record (asis fits, margins forms after cfprobit, the scratch
net install and the rebuild path) are `tests/*_2026-10-16.do` with their logs.
`tests/smoke.do` and `tests/smoke2.do` exercise every subcommand on every route with
both engines (rerun green under 0.1.1 on 2026-10-16, and again under 0.1.2 the same
evening together with `examples/showcase.do`; last section).

## Conventions that differ between the engines (documented, not defects)

* Standard errors: `engine(stata)` uses Stata's observed-information (OIM) e(V)
  and its small-sample factors (`vce(robust)`: N/(N-1); `vce(cluster)`: G/(G-1));
  `engine(mata)` uses R's expected information, HC0 without a factor, and G/(G-1).
  Differences are at the 1e-4 relative level for probit SEs and vanish for logit.
* `separation(keep)`: both engines keep the degenerate Wald interval (R's glm).
  `engine(mata)` keeps the Fisher-scoring fit; `engine(stata)` (since 0.1.1) keeps
  the `asis` fit, whose separated coefficient has an exactly zero variance row in
  e(V) (its SE prints as missing), and the delta method collapses onto the
  identified cells (for a binary focal with an empty arm: minus the other arm's
  rate, that arm's binomial SE). Stata's `asis` iteration stops with the intercept
  about 1e-4 off the closed form (the likelihood is flat along the separated
  coefficient), so the kept estimate differs from R's at the 1e-5 level, far inside
  Monte Carlo error. Terms omitted for collinearity (o.) still fail under either
  engine; a kept zero-variance fit counts as separated in the reported share.
* `effects(margins)` on the IV route (since 0.1.1): the AIE is the difference of
  `margins, dydx(d) at(m=(0 1))` after `cfprobit` with the focal-by-moderator term
  written as `1.d#c.m` for a binary focal (`c.d#c.m` for a continuous one): the same
  fit as before, only the stripe name changes, and `margins` can now differentiate
  the focal. The control function is held at its estimated value (cfprobit's
  predict(pr) is the average structural function), which is what `pape_cf_aie_est`
  computes. `margins, at(d=(0 1) m=(0 1))` is NOT the AIE: setting the endogenous
  variable through `at()` lets the generalized residual move (fixture H cells
  .526/.718/.788/.903 instead of the ASF cells .328/.518/.601/.763).

## 2026-10-16 additions (version 0.1.1)

Four items closed the port (CLAUDE.md session log 2026-10-16); pre-edit sources in
`legacy/powerape-stata-ado_pre-0.1.1_2026-10-16/`.

1. **Net-install packaging.** `stata.toc` and `powerape.pkg` (v 3, Distribution-Date
   20261016) install `ado/powerape.ado`, `ado/powerape.sthlp`, `ado/lpowerape.mlib`
   and the eight `.mata` sources (Stata treats `.mata` files as installation files);
   `build.do` is the ancillary file (`net get`). Tested by `net install powerape,
   from(<this directory>)` into a scratch PLUS directory (`sysdir set PLUS`; never
   `~/ado/plus`): `which powerape` and `findfile lpowerape.mlib` resolve to the scratch
   tree, `powerape version` reports both files with r(match) = 1, both engines run from
   the installed copy, `ado uninstall powerape` removes it. A second run exercised
   build.do's installed-copy path (project tree absent): it found the sources with
   `findfile` and replaced the installed library in place. `powerape version` compares
   the ado-file's version constant with the first word of the library's
   `pape_version()`, prints both locations, warns on a mismatch and errors with a
   pointer to `mata mlib index` when the library is missing.
2. **Stata versions.** `version 19` became `version 16` in the ado-file, the eight
   `.mata` sources and build.do. Feature inventory: frames and `st_frame*()` (Stata
   16); `cfprobit` (Stata 18, IV route under `engine(stata)` only); the community
   `ginteff` (`effects(margins)` for the AIE). Only StataNow 19.5 (update 14 Oct 2026)
   is installed here, so Stata 16-18 compatibility is a documented assumption, and the
   shipped library, compiled by 19.5, may not load in an older Stata (rebuild path
   documented in the help, README and pkg, and tested). Evidence that version control
   changed nothing on 19.5: the 153 pre-existing Tier B rows and the 11 pre-existing
   Tier A rows are bit-identical to the 2026-10-15 results, and the 73 pre-existing
   tests pass unchanged.
3. **`separation(keep)` under `engine(stata)`: implemented** (conventions above).
   Tier A row `stata (probit, asis)`: 0.5705 at 4,000 replications against the exact
   0.5602 (|diff| 0.0103 = 1.3 MCSE; 2 of the 4,000 `asis` fits did not converge and
   count against power; on the pre-fix ado-file of 2026-10-16 this row ran through
   `margins` (see the Tier C/D rerun section); the rerun under the leak-fixed 0.1.2 library
   gives the same 0.5705 through `pape_fit_stata`, see the last section). Tier B rows `sep_sat`/`sep_cov` (the asis
   fit against the Mata engine and against the other arm's binomial SE). run-tests.do:
   same-seed equality with the Mata engine. Limits: the IV route and the panel route
   (`vce(cluster)`) accept zero-variance coefficients under `asis` by the same rule but
   have no exact-value check; with
   `effects(margins)` the separated share comes from the cell check only (no
   Hauck-Donner or zero-row flag on that path), the scoring is unaffected.
4. **`effects(margins)` for the IV AIE: implemented** (conventions above). Tier B rows
   `H_aie`/`H_aie_c` (margins against the R closed form on the same data); run-tests.do:
   same-seed equality with the closed form. The `margins, at(d=(0 1) m=(0 1))` form was
   tried first and rejected because it does not hold the control function fixed.

Tiers C and D were rerun later the same day; see the next section.

## 2026-10-16 Tier C/D rerun under 0.1.1

Tiers C and D were rerun on 2026-10-16 against the 0.1.1 library (`v_tierC.do` 7 min,
`v_tierD.do` 51 min against about 28 min on 2026-10-15, with six other Stata batch
jobs running; `tests/version_probe_2026-10-16.do` records `powerape version` = 0.1.1
with the ado-file and the library in agreement). Both result files are byte-identical
to the 2026-10-15 files (Tier C 30/30 PASS, Tier D 45/45 PASS), and so are the two
logs; the 2026-10-15 copies are in
`legacy/powerape-stata-validation_pre-0.1.1-rerun_2026-10-16/`. The R reference
(`data/tierD-reference.csv`, powerape 1.11.3) was not regenerated.

Nothing moved, and the sources say why. Tier C's integration draws are seeded
(`seedint()`, default 20260812), so the calibration constants do not depend on the
session's RNG state. On the simulation path the 0.1.1 changes are inert under the
default `fail` convention: the new `rank < k` separation flag in `pape_sim_ci` can only
fire on a fit that passed `ok`, which requires full rank for the Mata engine and for
the Stata engine without `asis`; `pape_fit_stata` keeps zero-variance rows only under
`asis`; the `cfprobit` changes apply only with a moderator; and `effects(margins)` on
the IV route without a moderator (row D6c) already went through `margins` in 0.1.0.
Tier D has no row with `separation(keep)`.

One correction to the record. The Tier D results line for "D8 required n, detect 80%"
prints n in `%8.4f` as 2.5e+03 with rel.diff 4.0e-04, which fits 2,493 and 2,495
alike, so byte-identity does not pin that value. `tests/tierD_D8_n_probe_2026-10-16.do`
reruns the search (same DGP, seed 221, nsim 1,500) under 0.1.1 and returns n = 2,495
(confirmed power .7997, MCSE .0051); the same search on the pre-0.1.1 sources compiled
in memory (`pape_version()` 0.1.0, no library on the adopath;
`tests/tierD_D8_n_probe_v010_2026-10-16.*`) also returns 2,495 with the same confirmed
power. The search is seed-driven throughout (pilot, steps and confirmation all derive
their seeds from `seed()`), so the value inside the tier run is the same. The "2,493"
in the Tier D table (and in the 2026-10-15 CLAUDE.md entry) was a recording error,
corrected above; the minimum-claim search gives 5,461 under both versions, as recorded.

A second finding follows from the engine-leak fix that a sibling session made the same
afternoon (`powerape.ado` 16:56, ado-file only, `lpowerape.mlib` unchanged; CLAUDE.md
2026-10-16). Before that fix a per-command `engine()`/`effects()` leaked into the
session default, and an empty `effects()` kept the leaked value. `v_tierD.do` calls
`effects(margins)` once, in D3 (nsim 300), and afterwards runs four `engine(stata)` rows
without an `effects()` option: D5 panel probit `vce(cluster)`, D6 IV binary focal
`cfprobit` (nsim 400), the D6c `cfprobit` row at nsim 400, and the D8 Stata-engine
search. On the leaky ado-file those four rows therefore went through `margins` (or
`cfprobit` + `margins`), not through the closed-form Jacobians their labels name. The
rerun described above started at 16:47 with the pre-fix ado-file loaded, so it has that
exposure. Whether the 2026-10-15 run had it the record does not settle: the
end-of-2026-10-15 ado-file (legacy copy) carries the same `pape_setengine` logic, but
that logic arrived in a review fix late that day, and the original run's 28 minutes sit
closer to the 23 minutes of the fixed-ado-file run below than to the 51 minutes of the
exposed rerun. The values stand as power estimates for the Stata engine (`margins` and the
closed form agree to 1e-7 on every fixed dataset in Tier B, and both read the same
e(b)/e(V)), but the four labels overstated which effects path the rows exercised, and
the Tier D table's distinction between the D6c `cfprobit` row and the `cfprobit` +
`margins` row was void. Tier D was therefore run once more under the fixed ado-file
(17:54-18:17). Tier C has no simulation command and is unaffected. Tier A's two
rare-outcome `engine(stata)` rows (fail and keep) follow its `effects(margins)` row in
the same way and were last run (2026-10-16 09:53) on the pre-fix ado-file, so its 0.1.1
`keep` row came from `pape_effects_margins` with `asis` rather than from
`pape_fit_stata`'s kept-fit logic; the Tier A rerun under the fixed library was done
later that night (last section) and closes this. The
kept-fit logic itself was covered meanwhile: run-tests.do's same-seed keep check runs through
`pape_fit_stata` (it precedes the suite's only `effects(margins)` call) and passed under
the fixed ado-file (79/79, 17:22).

Under the fixed ado-file (`powerape.ado` of 16:56, md5 32fb31b5...; `lpowerape.mlib`
unchanged) Tier D again gives 45/45 PASS, and `tierD-results.txt` and `v_tierD.log` are
byte-identical to both the first rerun and the 2026-10-15 files: on these designs
`margins` and the closed-form Jacobians produced the same decision in every
replication, and the D8 Stata-engine search returned the same 783. The files on disk
are from this run, so every Tier D row now exercised the path its label names; the
first rerun's files are kept in the legacy folder as
`tierD-results_rerun1-pre-leak-fix-ado_2026-10-16-1739.txt` and the matching log. The
runtime fell from 51 to 23 minutes, consistent with the four rows no longer calling
`margins`. Open: the `recr` lines print required n with `%8.4f`, so a required-n row
cannot be pinned from the results file alone; printing `r(n)` in full (a format change
that would break bit-comparison with the earlier files) is left for the next revision
of the tier do-files.

## 2026-10-16 version 0.1.2: the `fit:` line; the engine-leak fix folded in

Version 0.1.2 (ado-file, Mata library and help bumped together; `powerape version`
reports r(match) = 1) makes every simulation command say who fitted the simulated
studies. `power`, `curve`, `n`, `mde` and `robust` print one line after the design line,
of the form

    fit: probit via engine(stata), effects(closed)
    fit: probit, vce(cluster) via engine(stata), effects(closed)     (panel route)
    fit: logit, vce(robust) via engine(stata), effects(closed)       (se(robust))
    fit: probit, asis via engine(stata), effects(closed)             (separation(keep))
    fit: cfprobit via engine(stata), effects(margins)                (IV route)
    fit: Mata engine (Fisher scoring)
    fit: Mata engine (Fisher scoring; control-function probit)       (IV route)

and store the same text in `r(fit)`, with its parts in `r(engine)` and `r(effects)`
(`robust` too). The label is built when the result is created (`pape_fit_lab` and
`pape_stamp_fit` in `powerape_stata_engine.mata`, called from `pape_power_once`,
`pape_curve`, `pape_n` and `pape_finish_mde`; the sweep builds its own), inside the
Mata call of the command, while a per-command `engine()`/`effects()` is still in
force; the print routines read the stored fields, so a result printed later (or a
`powerape alpha` on it, which sets the engine itself) cannot be mislabelled by the
session default. The existing output lines are unchanged. Nothing on the simulation
or estimation path changed; the struct `pape_res` gained the three string fields.

The 0.1.2 library also carries the engine-leak fix made in the ado-file alone earlier
on 2026-10-16: before it, a per-command `engine()`/`effects()` was written into the
Mata session externals and never restored, so one `effects(margins)` call made every
later command without an `effects()` option run through `margins`. `pape_setengine`
now saves the previous externals and `pape_runmata` restores them after each
simulation command (not after the nested `powerape dgp` calls a sweep issues). The
fix added two checks to `tests/run-tests.do` (77 -> 79); its exposure of the Tier A
and Tier D records is described in the previous section.

Runs under 0.1.2 (all logs kept; the pre-0.1.2 sources, library, documentation and
the logs they overwrite are in `legacy/powerape-stata-ado_pre-0.1.2_2026-10-16/`):

* `build.do`: library rebuilt in 5 s (`build.log`, 18:36) after the Tier D rerun's
  log had ended (18:17); the six StataSE processes running at the time were JW's
  unrelated batch jobs in other folders (checked by their command lines).
* `tests/run-tests.do`: 94/94 PASS (79 + 15 new); check 7b confirms
  r(match) = 1. A first attempt stopped at the named-log check because Mata's
  `fileread()` is not found on this Stata (r(3499)); the check now uses `cat()`
  (first-attempt log kept in the legacy folder).
* `tests/smoke.do`: green (no r() codes); fit lines `probit via engine(stata),
  effects(closed)` x2, `probit via engine(stata), effects(margins)` x2 (one the
  `ginteff` AIE), `logit, vce(robust) via engine(stata), effects(closed)`,
  `Mata engine (Fisher scoring)` x12.
* `tests/smoke2.do`: green (no r() codes; 22 min alongside the other jobs); fit lines `cfprobit via engine(stata), effects(closed)` x3 and `effects(margins)` x2, `Mata engine (Fisher scoring; control-function probit)` x3 (IV), `probit, vce(cluster) via engine(stata), effects(closed)` x3 and `effects(margins)` x1 (panel), `logit via engine(stata), effects(closed)`, `probit via engine(stata), effects(closed)`, `Mata engine (Fisher scoring)` x14; the sweep header carries the line after its range line.
* `examples/showcase.do`: green in 8.5 min; fit lines `probit via engine(stata),
  effects(closed)` x2, `probit via engine(stata), effects(margins)` (section 6,
  `ginteff`), `probit, vce(cluster) via engine(stata), effects(closed)` (panel),
  `cfprobit via engine(stata), effects(closed)` and `Mata engine (Fisher scoring;
  control-function probit)` (IV), `Mata engine (Fisher scoring)` x13.
* `tests/help_pkg_probe_v012_2026-10-16.do`: `powerape version` r(match) = 1 with
  r(ado_version) = r(mlib_version) = 0.1.2; `net describe powerape` reads the edited
  `.pkg` (Distribution-Date 20261016, the version paragraph, eleven installation
  files, `build.do` ancillary). `tests/help_render_probe_v012_2026-10-16.do` renders
  the edited help file through `type, smcl` without error (`help` itself is ignored
  in batch mode).

Tiers A-D were not rerun by the 0.1.2 session: the change touches printing and r() only,
so `tier?-results.txt` would be unchanged, while `v_tier?.log` would gain one `fit:`
line per simulation command (where the call does not pass `nodisplay`) and no longer be
byte-identical to the 2026-10-16 logs. The one rerun that was needed, Tier A under the
fixed ado-file, was done later that night (next section).

## 2026-10-16 Tier A rerun under 0.1.2

Tier A was rerun under the 0.1.2 library on 2026-10-16 to close the gap described in the
Tier C/D section: its two rare-outcome `engine(stata)` rows follow the `effects(margins)`
row, and on the pre-fix ado-file the `margins` setting leaked into them, so the 0.1.1
`separation(keep)` row (0.5705) had gone through `pape_effects_margins` with `asis`
instead of `pape_fit_stata`'s kept-fit logic.

* **Run:** `stata-do powerape-stata/validation/v_tierA.do`, StataNow 19.5 (update 14 Oct
  2026), started 20:36, finished 21:07 (31 min; JW's six `ltz_*` batch jobs ran
  throughout, no other process from `powerape-stata/` was running; one `end of do-file`,
  no r() code). Library: `powerape.ado` md5 b7c164074ee286798daa54bbb3a8fbec
  (18:35, 0.1.2 with the engine-leak fix), `lpowerape.mlib` rebuilt 18:36;
  `powerape version` reports 0.1.2 for both, r(match) = 1. A second session had been
  started for the same rerun; it was told to stand down before it launched anything
  (only one `v_tierA.do` Stata process existed), so both files on disk come from this
  single run.
* **Result: 12/12 PASS. `tierA-results.txt` and `v_tierA.log` are byte-identical to the
  09:53 files** (md5 9abb7096... for the results file); the 09:53 copies are in
  `legacy/powerape-stata-validation_pre-0.1.2-tierA-rerun_2026-10-16/`. The log does not
  gain `fit:` lines because `v_tierA.do` passes `nodisplay` to every simulation call and the
  `fit:` line is part of the printed output; `r(fit)` is stored regardless.
* **The keep row:** `.05 -> .01, n = 300, separation(keep), engine(stata)` gives
  **0.5705** (exact 0.5602, MCSE 0.0078, |diff| 0.0103 = 1.3 MCSE, well inside the 4 MCSE
  criterion), with 2 of the 4,000 `asis` fits not converged, as on 2026-10-16 09:53. The
  separated share of the Mata `fail` run is .22395, also unchanged. The same power from the
  two code paths means the same number of rejections in the 4,000 replications. That the
  09:53 run went through `margins` is inferred from the leak analysis, not observed in its
  log; since `margins` and the closed-form Jacobians read the same `e(b)`/`e(V)` and agree
  to 1e-7 in Tier B, they almost surely took the same decision in every replication (only
  one within 1e-7 of a boundary could differ).
* **Which path ran** (not visible in the log): `tests/tierA_path_probe_2026-10-16.do`
  replays the Tier A order at small nsim (100/200/200) and prints `r(fit)`: after the
  `effects(margins)` row (`effects=margins`) the `fail` row reports `probit via
  engine(stata), effects(closed)` and the `keep` row `probit, asis via engine(stata),
  effects(closed)`, with `assert`s on all three (PASS; keep power at nsim 200 is .595,
  no failed fits). So in the 0.1.2 rerun the `keep` row exercised `pape_fit_stata`'s
  kept-fit path with the closed-form effects, as its label says, and the `fail` row the
  plain closed-form `probit` path; the earlier exposure is closed. Tier A needs no
  further rerun for the leak.
* **Tiers B-D** are unaffected as recorded above (B: no simulation commands; C: none; D
  rerun under the fixed ado-file, byte-identical). Their `v_tier?.log` files have
  not been regenerated under 0.1.2.

## 2026-10-17 version 0.1.3: panel calibration speed

Version 0.1.3 changes one statement of the panel route. `pape_panel_reg_draw`
(`powerape_panel.mata`) built the unit index of a panel draw by appending one unit at
a time, `rd.unit = rd.unit \ J(T_i[j], 1, j)` for every unit `j`, which copies the
growing vector at each step, so its cost grew with the square of the number of units.
At the default `nint(1e5)` (100,000 units, 400,000 rows over four periods) that loop
was essentially the whole cost of an integration draw, about a minute, and every panel
`dgp`, `setape`, `setaie`, `trueape` and `trueaie` draws the integration set once, each
scenario of a robustness sweep twice (rebuild and re-pin) and an `mde` search once per
candidate effect. The new `pape_unit_index()` puts a 1 on each unit's first row and
takes a running sum (R's `rep(seq_len(n_units), times = T_i)`). The index holds the
same integers and draws no random numbers, so every number downstream is bit-identical.
The same function builds the unit index of every simulated panel study
(`pape_draw_x_panel`), so the simulated results are unchanged too. The rest of the
library is untouched apart from the version string (`pape_version()`), and the
ado-file changed only its version.

**Where the time went** (`tests/panel_speed_2026-10-17.do`, part B: the library's
draw timed whole, and a step-timed copy of its balanced branch). On 0.1.2 the draw took
121.8 s; in the copy, the unit-index loop took 155.8 s while everything else together
took under 0.6 s (random draws 0.08 s, `panelsetup` 0.44 s, transforms and Mundlak
means 0.06 s), and the two calibration root-finders 0.7 s and 0.6 s (14 objective
evaluations for `beta0`). The integration draw is not rebuilt inside the root-finders.
On 0.1.3 the draw takes 0.36 s and the vectorized index 0.02 s.

**Timings** in seconds at `nint(1e5)` (`tests/panel_speed_before_timings.txt`;
`panel_speed_after_quiet_timings.txt` and `panel_speed_r_timings.txt`, the in-place
0.1.3 library and R run one after the other at 19:32-19:36 with no other Stata job on
the machine; `panel_speed_after_timings.txt` and `panel_speed_r_timings_loaded.txt`,
the scratch build of the same sources and R at 18:42-18:47 with five to six other
Stata jobs). Two sibling sessions ran up to nine Stata jobs on this 8-thread machine
during part of the 0.1.2 run: q = no other Stata job running, L = five to nine
(`trueaie`: the other jobs started near its end), so the L rows overstate the speed-up
against the quiet columns. The first command of a session also loads the Mata library
(the first row). The `power` row shows the cost of simulation for scale: the change
barely touches it (with the loop, the unit index of a simulated study of 300 units cost
about a millisecond). Single runs vary by tens of percent (R's `mde`: 18.3 s and
24.3 s).

| Step | 0.1.2 | 0.1.3 quiet | R quiet | 0.1.3 load | R load |
|---|---:|---:|---:|---:|---:|
| `dgp`, balanced, 4 periods | 60.4 q | 5.2 | 1.8 | 3.2 | 1.8 |
| `setape` | 62.9 q | 3.0 | 1.3 | 2.4 | 1.4 |
| `trueape` | 61.2 q | 1.1 | 0.5 | 1.2 | 0.5 |
| `dgp`, AIE world | 68.0 q | 1.8 | 1.9 | 1.6 | 1.7 |
| `setaie` | 58.6 q | 2.6 | 3.3 | 3.2 | 3.0 |
| `trueaie` | 58.4 q/L | 1.3 | 0.8 | 1.7 | 0.8 |
| `dgp`, retention .85, 6 periods | 148.1 L | 1.1 | 1.7 | 1.3 | 1.6 |
| `setape` (retention world) | 153.4 L | 1.4 | 1.4 | 2.1 | 1.4 |
| `power`, n 300, nsim 100, Mata engine | 10.8 L | 6.0 | 0.9 | 3.6 | 1.0 |
| `robust`, per scenario (2 scenarios, nsim 100) | 323.3 L | 3.2 | 3.9 | 4.0 | 3.9 |
| `mde`, n 400, nsim 100, Mata engine | 675.8 L | 14.6 | 24.3 | 18.3 | 18.3 |
| integration draw (`pape_panel_reg_draw`) | 121.8 L | 0.34 | | 0.36 | |

The report that prompted the change measured 70 s for `dgp` and 63 s for `setape` on
0.1.2. A same-load pair (`tests/panel_tierC_hex_2026-10-17.do`, 0.1.2 then 0.1.3 back
to back with eight to nine other Stata jobs running) gives `dgp` 163.7 -> 7.5 s,
`setape` 183.8 -> 7.1 s, `trueape` 175.2 -> 2.9 s, and on the Tier D retention world
`dgp` 129.5 -> 2.9 s, `setape` 144.8 -> 1.9 s. Sweeps and `mde` searches now run as fast
as in R or faster; single calibration calls take 1.1 to 3.0 s against R's 0.5 to 3.3 s
(the first command of a session longer, as it loads the library). The quiet 0.1.3 run
also reproduced the 166-line constants file below byte for byte, this time on the
in-place library.

**Results unchanged** (all files under `tests/` unless stated; the pre-0.1.3 sources,
library, documentation and every log or result file the reruns overwrite are in
`legacy/powerape-stata-ado_pre-0.1.3_2026-10-17/`):

* `panel_speed_before_constants.txt` (0.1.2), `panel_speed_after_constants.txt` (the
  0.1.3 scratch build) and `panel_speed_after_quiet_constants.txt` (the in-place 0.1.3
  library) are byte-identical (166 lines, md5 226818a5...): every stored constant
  (`beta0`, `gamma`, `xi`, `xbar_center`, `s2`, `beta_focal`, `p1`, `beta_mod`,
  `beta_int`, `var_a`, `scale_a`, `expected_T`) of six panel worlds (balanced; AIE;
  retention .85; a 2/4/6 length mixture; a continuous focal; no covariates) in `%21x`
  (exact), the true APE/AIE, seeded power, MCSE and separation counts under both
  engines (the Stata engine's `probit, vce(cluster)` included), the two-scenario sweep
  matrix and the `mde`. The probe also checks the old loop against the new index on
  the full draw and on balanced, retention, length-mixture and edge-case supports
  (exact equality), and that re-solving `beta0` and `beta_focal` on the draw returns
  the stored values.
* `panel_tierC_hex_old.txt`, `_new.txt` and `_inplace.txt` (0.1.2; the 0.1.3 scratch
  build; the in-place 0.1.3 library): byte-identical constants of the Tier C world
  `calp` (= Tier D's D5) and of Tier D's D5u in `%21x`; the seven-decimal values
  equal `tierC-results.txt`.
* Tier C (`v_tierC.do`, rerun from the repository root under 0.1.3, 18:50-18:53):
  30/30 PASS; `tierC-results.txt` byte-identical to the 2026-10-16 file (md5
  7cb1e9b5...). The log `v_tierC.log` is now the script's own log (opened after the
  set-up, since the scripts became root-relative the same day); from the first command
  on it is line-for-line the 2026-10-16 log, apart from its last line.
* Tier D (`v_tierD.do`, from the root under 0.1.3, 18:50-19:26, 36 min under load):
  45/45 PASS; `tierD-results.txt` byte-identical to the 2026-10-16 file (md5
  3cc8cd2c...), so the panel rows (D5: .5423 with the Mata engine, .5370 with
  `probit, vce(cluster)`; D5u: .4828) and every other row are unchanged; from the first
  command on, `v_tierD.log` is line-for-line the 2026-10-16 log apart from its last
  line. The 2026-10-16 C and D files compared against were produced under 0.1.1 with
  the leak-fixed ado-file (0.1.2 changed only printing; see above).
* `tests/run-tests.do` (from the repository root, 18:50-18:59): 94/94 PASS. Its log
  (`tests/run-tests.log`) equals the 2026-10-16 one apart from the version strings, the
  path separators in the `powerape version` lines (the script now builds its paths from
  `c(pwd)`) and the temporary file name and time stamps of the named-log check.
* `tests/smoke2.do` (from the root, 18:50-18:56; 5.5 min against 22 on 2026-10-16):
  green, no r() code. Its output equals the 2026-10-16 log apart from three version
  strings in printed statements and its four simulation timers; the `fit:` lines are
  those listed for 0.1.2.
* `examples/showcase.do` (from the root, 18:53-18:59; 6.0 min under load, so its header
  now says five to ten minutes instead of ten to fifteen): green; output identical to
  the 2026-10-16 log apart from the version strings and the `ginteff` installation
  comment that was corrected in the source earlier that day.
* `build.do` rebuilt `ado/lpowerape.mlib` in place at 18:48:50 (`build.log`), after
  the same sources had been compiled in a scratch copy for the timings above;
  `powerape version` reports 0.1.3 for the ado-file and the library, r(match) = 1
  (`tests/help_pkg_probe_v013_2026-10-17.do`, which also renders the edited help
  through `type, smcl` and reads the edited `.pkg` with `net describe`).

## 2026-10-17 rerun from the repository root with the published scripts

The validation now runs from a clone of this repository. The tier do-files and
`tests/run-tests.do` (and the local smoke tests) take the repository root from
`c(pwd)` and stop with r(601), "run it from the repository root", unless Stata was
started there: `do validation/v_tierA.do`. A batch wrapper that starts Stata in the
do-file's own folder is refused that way (checked with `v_tierC.do`).
`export_fixtures.R` and `tierD_r.R` write to `validation/data` relative to the working
directory, also the repository root (`Rscript validation/tierD_r.R >
validation/tierD_r.log`). Each tier first runs `powerape version` and stops unless the
ado-file and the Mata library agree. It then opens its own log, `validation/v_tier?.log`,
after the set-up and with `nomsg`, so the logs carry neither Stata's licence banner nor a
local path. Their first output line names the versions (`powerape 0.1.3 (Stata port of
powerape 1.11.3); Stata 19.5 (14 Oct 2026)`). From the first command after the set-up
on, each log equals the earlier batch log line for line; only the last line differs
(`log close tier` instead of `exit, clear` and `end of do-file`).

* **The path change alone, under 0.1.2.** The tiers ran in a scratch copy holding exactly
  the 41 files the repository's `.gitignore` then published (listed by `git status
  --porcelain --untracked-files=all` in a `git init -b main` copy; nothing committed), so
  they used nothing outside the published set. Each tier was started by a one-line driver
  that `cd`s to the copy's root and runs `do validation/v_tier?.do`, as a user of a clone
  would. Library 0.1.2: `powerape.ado` md5 b7c164074ee286798daa54bbb3a8fbec,
  `lpowerape.mlib` md5 ad9ed7d0db1412e019f56acf334a9586, both loaded from the copy's `ado/`
  (`powerape version` in the drivers' batch logs). **A 12/12, B 171/171, C 30/30 and D
  45/45 PASS, and every `tier?-results.txt` is byte-identical to the 2026-10-16 file**
  (md5 A 9abb7096..., B 7b6876b2..., C 7cb1e9b5..., D 3cc8cd2c...). The runs shared the
  machine with up to eight other Stata jobs: A 44 min, B 1.7, C 10.9, D 56.
  `tests/run-tests.do` from the copy's root gave 94/94 PASS (33 min under that load), and
  `tests/smoke.do` and `tests/smoke2.do` ran green (no r() code).
* **The published logs, under 0.1.3.** Tiers A and B were rerun under 0.1.3
  (`powerape.ado` md5 dbb27083d2c81f98d86282f0a43f2feb, `lpowerape.mlib` md5
  fb57d686a1aeb9b2c64489b49121f3ac) in a second scratch copy of the published files:
  A 12/12 PASS (28 min) and B 171/171 PASS (29 s), results byte-identical. The Tier C
  and D runs under 0.1.3 are the ones in the previous section. All four logs in this
  folder come from these 0.1.3 runs, and each equals its 0.1.2 counterpart apart from
  the version line.
* **R side.** From the root of a separate copy (R 4.5.2, powerape 1.11.3),
  `export_fixtures.R` rewrites its nine files (`fixA.csv` to `fixH.csv`,
  `reference.csv`) byte for byte (2 min), and `tierD_r.R` its two (`tierD-reference.csv`,
  `pilot600.csv`; 20 min under load). `tierD_r.log` is the output of that run.
* **What the repository holds.** The allowlist now publishes 46 files: the package, the
  examples and this record; the four tier scripts with their results and logs;
  `export_fixtures.R`, and `tierD_r.R` with its log; the 11 files of `data/`; and
  `tests/run-tests.do`. A search of all 46 for user and folder names, drive-letter paths,
  `AppData` and licence details finds nothing; the only drive-letter-shaped hit is SMCL
  markup in the help file. Kept local: `tests/run-tests.log` (it prints a temporary
  file's path and the library's location), the smoke tests, the probes
  `tests/*_2026-10-1?.do` with their logs, and the `legacy/` folders that earlier
  sections cite.
