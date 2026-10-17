* showcase.do -- a guided tour of powerape for Stata (version 0.1.3)
*
* Run it top to bottom (five to ten minutes on this machine, depending on
* what else is running) or step through it in the Do-file Editor. Every
* command is documented in -help powerape-.
* Replication counts are kept small so the tour runs quickly; for a real
* analysis use the defaults (nsim 1000-1500) or more.

version 16
clear all
set more off

* ---------------------------------------------------------------------------
* 0. Make the package visible
* ---------------------------------------------------------------------------
* Install it once (Stata 16 or later):
*     net install powerape, from("https://raw.githubusercontent.com/jespernwulff/powerape-stata/main/")
*     mata: mata mlib index
* Run from a copy of the repository instead, this file uses the copy's ado
* folder for this session only.
capture which powerape
if _rc & fileexists("ado/powerape.ado") {
    adopath ++ "`c(pwd)'/ado"
    mata: mata mlib index
}
powerape version

* ---------------------------------------------------------------------------
* 1. Describe the world (a DGP), then set the effect in probability points
* ---------------------------------------------------------------------------
* Probit, a binary treatment with P(treat = 1) = .5, one normal covariate,
* an outcome rate of .30 in the untreated and a covariate signal (latent
* pseudo-R2) of .15. DGPs are stored by name, like -estimates store-.
powerape dgp design_a, focal(treat, binary p(0.5))              ///
    covariates((age, normal mean(45) sd(12)))                  ///
    baseline(0.30) signal(0.15)

* The effect is the average partial effect: treatment raises the outcome
* probability by 5 percentage points. setape inverts this into coefficients.
powerape setape design_a, target(0.05)
powerape trueape design_a
powerape describe design_a

* ---------------------------------------------------------------------------
* 2. Power at one sample size
* ---------------------------------------------------------------------------
* The default claim is "minimum": the CI must exclude 0 AND the smallest
* effect size of interest (sesoi). Each simulated study is fit by Stata's
* own -probit- (the default engine); the APE (for a binary focal the average
* discrete change) and its delta-method SE equal what -margins- reports.
* The "fit:" line of the output names who fitted the studies (the Stata
* command and the effects path, or the Mata engine); r(fit) stores it.
powerape power design_a, n(2000) claim(minimum) sesoi(0.02) nsim(200) seed(1)
return list
matrix list r(outcomes)          // how the replications were classified
powerape statement               // a citable paragraph for the last result

* A plain detection claim (CI excludes 0), then the same seed under the
* Mata engine: a self-contained port of the R engine, 10-30 times faster.
powerape power design_a, n(2000) claim(detect) nsim(200) seed(1)
powerape power design_a, n(2000) claim(detect) nsim(200) seed(1) engine(mata)

* ---------------------------------------------------------------------------
* 3. Power curve, required n, minimum detectable effect
* ---------------------------------------------------------------------------
* Searches run thousands of replications, so make the fast engine the
* session default; confirm a final design with engine(stata) afterwards.
* (engine() or effects() given on a single command apply to that command
* only and leave the session default alone.)
powerape engine                      // shows the current defaults
powerape engine, engine(mata)        // sets them for the session

powerape curve design_a, n(1000 2000 3000 4000) claim(detect) nsim(300) seed(2)
matrix list r(results)

* Required n for 80% power: a search followed by a high-precision
* confirmation stage at the proposed n.
powerape n design_a, power(0.8) claim(detect) nsim(500) seed(3)
display "Required n: " r(n) "   power at that n (confirmation stage): " %5.3f r(power)

* The smallest APE detectable with 80% power at n = 2000.
powerape mde design_a, n(2000) claim(detect) nsim(500) seed(4)
display "MDE (probability points): " %6.4f r(mde)
powerape statement

* ---------------------------------------------------------------------------
* 4. Robustness: how much do the nuisance assumptions matter?
* ---------------------------------------------------------------------------
* Vary the baseline rate over .20-.40 (three grid points), keep the APE
* pinned at .05, and report the "insurance premium" n_max: the sample size
* that reaches 90% power in the least favourable scenario.
powerape robust design_a, n(2000) claim(detect) vary(baseline(0.20 0.40)) ///
    gridpoints(3) nsim(300) seed(5)
matrix list r(scenarios)
display "n_max = " r(nmax) " (scenario " r(nmax_scenario) ")"

* Two nuisance inputs at once, as MDEs instead of power, without n_max.
powerape robust design_a, n(2000) claim(detect)                  ///
    vary(baseline(0.20 0.40) signal(0.05 0.25)) gridpoints(2)    ///
    mode(mde) nonmax nsim(300) seed(6)
matrix list r(scenarios)

* ---------------------------------------------------------------------------
* 5. Justified alpha
* ---------------------------------------------------------------------------
* Store the replication draws of one power run, then re-evaluate them.
powerape power design_a, n(2000) claim(minimum) sesoi(0.02) nsim(400) seed(7)
powerape powerat, alpha(0.01 0.025 0.05 0.10)        // power at other alphas
powerape alpha, cost(4) prior(1)                     // Maier & Lakens optimum
powerape statement _lastalpha
powerape frontier, cost(4) power(0.8)                // closed-form frontier

* Or let alpha fall with n by a sample-size rule: .05 at n = 1000, shrinking
* as sqrt(1000/n); @n stands for the sample size.
powerape power design_a, n(2000) claim(detect) alpharule(0.05*sqrt(1000/@n)) ///
    nsim(300) seed(8)

* ---------------------------------------------------------------------------
* 6. Interactions: the average interaction effect (AIE)
* ---------------------------------------------------------------------------
* Treatment by a binary moderator. setaie fixes the AIE (.04) together with
* the two main effects (.10 for treat, .05 for female), all in probability
* points. The estimator equals -ginteff-.
powerape dgp design_b, focal(treat, binary p(0.5))                ///
    moderator(female, binary p(0.55))                             ///
    covariates((age, normal mean(45) sd(12)))                     ///
    baseline(0.30) signal(0.15)
powerape setaie design_b, target(0.04) mainfocal(0.10) mainmoderator(0.05)
powerape trueaie design_b
powerape power design_b, n(3000) claim(detect) nsim(300) seed(9)

* The same world fit by -probit- with the AIE computed by -ginteff- itself
* (about one to two seconds per replication, hence nsim(30)). Needs the
* community package ginteff (Radean, 2023), from the Stata Journal:
*     net install st0711, from(http://www.stata-journal.com/software/sj23-2)
powerape power design_b, n(3000) claim(detect) nsim(30) seed(9) ///
    engine(stata) effects(margins)
powerape engine                      // still the session default set above

* ---------------------------------------------------------------------------
* 7. Equivalence: power to show the APE is within +/- .05 of zero
* ---------------------------------------------------------------------------
powerape dgp design_null, focal(treat, binary p(0.5)) baseline(0.30)
powerape setape design_null, target(0)
powerape power design_null, n(3000) claim(equivalence) sesoi(0.05) nsim(300) seed(10)

* ---------------------------------------------------------------------------
* 8. Pilot data in memory: resampled covariates, or a fitted pilot model
* ---------------------------------------------------------------------------
* A small artificial pilot (in practice: your own pilot or a public file).
set seed 20261016
set obs 500
gen treat  = runiform() < .5
gen age    = rnormal(45, 12)
gen female = runiform() < .55
gen y      = runiform() < normal(-0.9 + 0.3*treat + 0.01*age + 0.2*female)

* empirical: draw covariate rows from the data, calibrate to baseline/signal
powerape dgp pilot_rows, empirical focal(treat) covariates(age female) ///
    baseline(0.30) signal(0.15)
powerape setape pilot_rows, target(0.05)
powerape power pilot_rows, n(1500) claim(minimum) sesoi(0.02) nsim(200) seed(11)

* fromfit: take the world from a fitted probit (its coefficients and data)
probit y treat age female
powerape dgp pilot_fit, fromfit focal(treat)
powerape setape pilot_fit, target(0.08)
powerape power pilot_fit, n(1000) claim(detect) nsim(200) seed(12)

* ---------------------------------------------------------------------------
* 9. Panel data: correlated random effects probit, unit-clustered SEs
* ---------------------------------------------------------------------------
* 300 units observed four times; treatment is time-invariant (icc(1)), the
* covariate is 60% between-unit; rho is the unit-effect share, creshare
* the share of it that is correlated with the regressors (Mundlak terms).
powerape dgp firms, panel focal(treat, binary p(0.5) icc(1))       ///
    covariates((size, normal icc(0.6)))                            ///
    nperiods(4) rho(0.3) creshare(0.5) baseline(0.30) signal(0.10)
powerape setape firms, target(0.08)
powerape trueape firms
powerape power firms, n(300) claim(detect) nsim(200) seed(13)
* with -probit ..., vce(cluster id)- doing the fitting:
powerape power firms, n(300) claim(detect) nsim(50) seed(13) engine(stata)

* ---------------------------------------------------------------------------
* 10. Endogenous treatment: control-function probit (Stata's -cfprobit-)
* ---------------------------------------------------------------------------
* A binary instrument (subsidy) shifts take-up of training; endogeneity is
* the error correlation, ivstrength the instrument's first-stage R2.
powerape dgp iv_design, iv focal(training, binary p(0.4))          ///
    covariates((size, normal)) instruments((subsidy, binary p(0.5))) ///
    endogeneity(0.4) ivstrength(0.2) baseline(0.30) signal(0.10)
powerape setape iv_design, target(0.10)
powerape trueape iv_design
powerape power iv_design, n(1500) claim(detect) nsim(200) seed(14)
* with -cfprobit- fitting each study (about a third of a second each):
powerape power iv_design, n(1500) claim(detect) nsim(30) seed(14) engine(stata)

* ---------------------------------------------------------------------------
* 11. Housekeeping
* ---------------------------------------------------------------------------
powerape engine, engine(stata) effects(closed)     // back to the default
powerape list                                      // stored DGPs and results
powerape drop design_null
powerape list

exit, clear
