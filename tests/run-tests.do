* run-tests.do -- regression tests for the powerape Stata port (assert-based).
* Fast: every simulation uses the Mata engine at small nsim; the Stata engine
* is exercised once per route. Error paths are checked with rcof.
* Run from the repository root (the folder that holds ado/ and tests/):
*     do tests/run-tests.do
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "run-tests.do: run it from the repository root:  do tests/run-tests.do"
    exit 601
}
adopath ++ "`root'/ado"
mata: mata mlib index
local npass 0
local nfail 0
program define check
    gettoken ok 0 : 0, match(par)
    gettoken label 0 : 0
    if (`ok') {
        di as txt "PASS  `label'"
        global PAPE_NPASS = $PAPE_NPASS + 1
    }
    else {
        di as err "FAIL  `label'"
        global PAPE_NFAIL = $PAPE_NFAIL + 1
    }
end
global PAPE_NPASS 0
global PAPE_NFAIL 0

* ---- 1. DGP construction, calibration, inversion ----------------------------
powerape dgp t0, focal(t, binary p(0.5)) baseline(0.30) nodisplay
check (abs(r(beta0) - invnormal(0.30)) < 1e-12) "beta0 = invnormal(baseline) without covariates"
powerape setape t0, target(0.10) nodisplay
check (abs(r(p1) - 0.40) < 1e-7) "binary focal: implied P(Y=1|focal=1) = baseline + APE"
check (abs(r(beta_focal) - (invnormal(0.40) - invnormal(0.30))) < 1e-7) "beta_focal closed form"
powerape trueape t0
check (abs(r(ape) - 0.10) < 1e-7) "true APE equals the target"
capture noisily powerape setape t0, target(0.75)
check (_rc == 498) "infeasible APE target is refused"
powerape dgp t1, focal(treat, binary p(0.5)) covariates((age, normal mean(45) sd(12)) (female, binary p(0.55))) ///
    correlation(0.2) baseline(0.30) signal(0.15) nodisplay
check (abs(r(signal) - 0.15) < 1e-6) "signal realized on the integration draw"
check (abs(r(baseline) - 0.30) < 1e-12) "baseline stored"
powerape setape t1, target(0.05) nodisplay
powerape trueape t1
check (abs(r(ape) - 0.05) < 1e-6) "true APE with covariates equals the target"
powerape dgp t1g, focal(treat, binary p(0.5)) covariates((age, normal mean(45) sd(12))) baseline(0.30) gamma(0.02) nodisplay
check (abs(el(r(gamma), 1, 1) - 0.02) < 1e-12) "explicit gamma() is used"
capture noisily powerape dgp bad, focal(t, binary p(1.5)) baseline(0.3)
check (_rc == 498) "binary p outside (0,1) is refused"
capture noisily powerape dgp bad, focal(t, binary p(0.5)) baseline(0.3) signal(0.2)
check (_rc == 498) "signal without covariates is refused"
capture noisily powerape dgp bad, focal(t, binary p(0.5)) covariates((z, normal)) baseline(0.3) correlation(1.5)
check (_rc == 498) "non-positive-definite correlation is refused"
capture noisily powerape dgp t0, focal(t, binary p(0.5)) baseline(0.3)
check (_rc == 110) "existing DGP name needs replace"
* continuous focal: APE plateaus
powerape dgp tc, focal(x, normal mean(0) sd(1)) baseline(0.30) nodisplay
powerape setape tc, target(0.05) nodisplay
powerape trueape tc
check (abs(r(ape) - 0.05) < 1e-6) "continuous focal: true APE equals the target"
capture noisily powerape setape tc, target(0.9)
check (_rc == 498) "continuous focal: target beyond the plateau is refused"
* AIE
powerape dgp ta, focal(treat, binary p(0.5)) moderator(female, binary p(0.55)) baseline(0.30) nodisplay
powerape setaie ta, target(0.04) mainfocal(0.10) mainmoderator(0.05) nodisplay
powerape trueaie ta
check (abs(r(aie) - 0.04) < 1e-6) "AIE binary x binary pinned exactly"
capture noisily powerape setape ta, target(0.1)
check (_rc == 498) "setape on a moderated DGP is refused"
capture noisily powerape setaie t0, target(0.1) mainfocal(0.1) mainmoderator(0.1)
check (_rc == 498) "setaie without a moderator is refused"
foreach pair in "x normal|m binary p(0.5)" "treat binary p(0.5)|w normal" "x normal|w normal" {
    gettoken f m : pair, parse("|")
    local m = subinstr("`m'", "|", "", 1)
    gettoken fn fs : f
    gettoken mn ms : m
    powerape dgp tp, focal(`fn', `fs') moderator(`mn', `ms') covariates((z, normal)) baseline(0.35) signal(0.1) replace nodisplay
    powerape setaie tp, target(0.03) mainfocal(0.05) mainmoderator(0.06) nodisplay
    powerape trueaie tp
    check (abs(r(aie) - 0.03) < 1e-5) "AIE pinned: `fn' x `mn'"
}

* ---- 2. power: claims, coherence, outcome table, both engines -----------------
powerape power t0, n(712) claim(detect) nsim(300) seed(1) engine(mata) nodisplay
check (r(power) > 0.70 & r(power) < 0.90) "detect power near the four-tool anchor (.80 at n = 712)"
check (abs(r(conf) - 0.95) < 1e-12 & abs(r(alpha) - 0.05) < 1e-12) "detect claim uses the 95% interval"
check ("`r(engine)'" == "mata" & "`r(effects)'" == "closed" & "`r(fit)'" == "Mata engine (Fisher scoring)") "r(fit) names the Mata engine"
local pm = r(power)
powerape power t0, n(712) claim(detect) nsim(300) seed(1) engine(stata) nodisplay
check (abs(r(power) - `pm') <= 0.03) "stata and mata engines agree on the same seed"
check ("`r(engine)'" == "stata" & "`r(effects)'" == "closed" & "`r(fit)'" == "probit via engine(stata), effects(closed)") "r(fit) names probit under engine(stata), effects(closed)"
* the printed output carries the fit line (read back through a named log)
tempfile fitlog
log using "`fitlog'.log", text replace name(pape_fitchk)
powerape power t0, n(300) claim(detect) nsim(30) seed(1) engine(stata)
log close pape_fitchk
mata: st_local("hasfit", strofreal(any(strpos(cat(st_local("fitlog") + ".log"), "  fit: probit via engine(stata), effects(closed)") :> 0)))
check (`hasfit') "the power output prints the fit line"
capture erase "`fitlog'.log"
powerape power t1, n(2000) claim(minimum) sesoi(0.02) nsim(200) seed(2) engine(mata) nodisplay
check (abs(r(conf) - 0.90) < 1e-12) "minimum claim uses the 90% interval"
mata: o = st_matrix("r(outcomes)"); st_local("sumo", strofreal(sum(o), "%12.0g"))
check (abs(`sumo' - 1) < 1e-9) "outcome distribution sums to one"
capture noisily powerape power t1, n(2000) claim(minimum) sesoi(0.05) nsim(100)
check (_rc == 498) "minimum claim with target <= sesoi is refused"
capture noisily powerape power t1, n(2000) claim(minimum) nsim(100)
check (_rc == 498) "minimum claim without sesoi is refused"
capture noisily powerape power t1, n(2000) claim(equivalence) sesoi(0.02) nsim(100)
check (_rc == 498) "equivalence with |target| >= sesoi is refused"
powerape power t1, n(2000) claim(detect) conf(0.99) nsim(100) seed(3) engine(mata) nodisplay
check (abs(r(alpha) - 0.01) < 1e-12) "conf() override sets alpha"
powerape power t1, n(2000) claim(detect) alpharule(0.05*sqrt(1000/@n)) nsim(100) seed(3) engine(mata) nodisplay
check (abs(r(alpha) - 0.05*sqrt(1000/2000)) < 1e-9) "alpha rule evaluated at n"
powerape power t1, n(2000) claim(detect) nsim(100) seed(3) se(robust) engine(mata) nodisplay
check ("`r(se)'" == "robust") "robust SEs selected"
powerape power t1, n(2000) claim(detect) nsim(40) seed(3) se(robust) engine(stata) nodisplay
check ("`r(fit)'" == "probit, vce(robust) via engine(stata), effects(closed)") "r(fit) shows vce(robust) under engine(stata)"
* separation handling: rare outcome, protective effect
powerape dgp tr, focal(t, binary p(0.5)) baseline(0.05) nodisplay
powerape setape tr, target(-0.04) nodisplay
powerape power tr, n(300) claim(detect) nsim(400) seed(4) engine(mata) nodisplay
check (r(separated) > 0.10) "separation detected in the rare-outcome design"
local pf = r(power)
powerape power tr, n(300) claim(detect) nsim(400) seed(4) engine(mata) separation(keep) nodisplay
check (r(power) > `pf' + 0.10) "keep convention raises power by the separated share"
local pk = r(power)
local sk = r(separated)
* the Stata engine under separation(keep) fits with asis and keeps the separated
* coefficient's zero variance row: the same degenerate intervals as the Mata engine
powerape power tr, n(300) claim(detect) nsim(400) seed(4) engine(stata) separation(keep) nodisplay
check (abs(r(power) - `pk') <= 0.0125 & abs(r(separated) - `sk') < 1e-12) "keep under engine(stata) (asis) reproduces the Mata engine's keep power on the same seed"
check (r(n_failed) <= 2) "keep under engine(stata): separated replications are kept, not failed"
check ("`r(fit)'" == "probit, asis via engine(stata), effects(closed)") "r(fit) shows asis under separation(keep), engine(stata)"
* logit = probit on the saturated design, same seed
powerape dgp tl, model(logit) focal(t, binary p(0.5)) baseline(0.30) nodisplay
powerape setape tl, target(0.10) nodisplay
powerape power tl, n(500) claim(detect) nsim(300) seed(5) engine(mata) nodisplay
local pl = r(power)
powerape power t0, n(500) claim(detect) nsim(300) seed(5) engine(mata) nodisplay
check (abs(r(power) - `pl') < 1e-12) "probit and logit identical in the saturated case"
* seeded calls leave the caller's RNG stream untouched
set seed 777
local u1 = runiform()
set seed 777
powerape power t0, n(300) claim(detect) nsim(50) seed(9) engine(mata) nodisplay
local u2 = runiform()
check (abs(`u1' - `u2') < 1e-15) "seeded run restores the ambient RNG state"

* ---- 3. curve, n, mde ---------------------------------------------------------
powerape curve t0, n(300 600 900) claim(detect) nsim(200) seed(6) engine(mata) nodisplay
check ("`r(fit)'" == "Mata engine (Fisher scoring)") "curve stores r(fit)"
mata: R = st_matrix("r(results)"); st_local("mono", strofreal(R[1,3] < R[2,3] & R[2,3] < R[3,3]))
check (`mono') "power curve increases in n"
powerape n t0, power(0.8) claim(detect) nsim(500) seed(7) engine(mata) nodisplay
check ("`r(fit)'" == "Mata engine (Fisher scoring)") "n stores r(fit)"
check (r(n) > 600 & r(n) < 820) "required n near the analytic 712"
check (r(power) >= 0.8 - max(0.005, 1.5*r(mcse))) "confirmation rule holds at the returned n"
powerape mde t0, n(712) claim(detect) nsim(400) seed(8) engine(mata) nodisplay
check ("`r(fit)'" == "Mata engine (Fisher scoring)") "mde stores r(fit)"
check (abs(r(mde) - 0.10) < 0.012) "MDE at n = 712 near 0.10"
powerape mde t0, n(712) claim(detect) nsim(400) seed(8) engine(mata) direction(negative) nodisplay
check (r(mde) < 0) "negative direction returns a decrease"
powerape dgp te, focal(t, binary p(0.5)) baseline(0.30) nodisplay
powerape setape te, target(0) nodisplay
powerape mde te, n(2000) claim(equivalence) nsim(300) seed(9) engine(mata) nodisplay
check (r(is_margin) == 1 & r(mde) > 0.04 & r(mde) < 0.08) "equivalence MDE returns a margin"
powerape n te, power(0.8) claim(equivalence) sesoi(0.05) nsim(400) seed(10) engine(mata) noconfirm nodisplay
check (r(confirmed) == 0 & r(n) > 1500 & r(n) < 4500) "equivalence n without confirmation"

* ---- 4. alpha / powerat / frontier / statement ---------------------------------
powerape power t1, n(2000) claim(minimum) sesoi(0.02) nsim(400) seed(11) engine(mata) nodisplay
powerape powerat, alpha(0.01 0.05 0.10)
mata: T = st_matrix("r(table)"); st_local("monoa", strofreal(T[1,3] <= T[2,3] & T[2,3] <= T[3,3]))
check (`monoa') "power is monotone in alpha on stored draws"
powerape alpha, cost(4) prior(1) nodisplay
check (r(alpha) > 0 & r(alpha) <= 0.05) "justified alpha within the cap"
powerape frontier, cost(4) power(0.8)
check (abs(r(alpha) - 0.0274) < 5e-4) "frontier alpha 0.0274 at 4:1 and 80%"
powerape frontier, cost(1) power(0.8) claim(minimum)
check (abs(r(alpha) - 0.20) < 1e-9) "one-sided frontier at equal costs is 1 - power"
powerape statement
check (strpos(`"`r(statement)'"', "Riesthuis") > 0) "statement produced"

* ---- 5. empirical and pilot-model routes ---------------------------------------
clear
set seed 21
set obs 800
gen treat = runiform() < .5
gen age = rnormal(45, 12)
gen female = runiform() < .55
gen y = runiform() < normal(-0.9 + 0.3*treat + 0.01*age + 0.2*female)
powerape dgp e1, empirical focal(treat) covariates(age female) baseline(0.30) signal(0.15) nodisplay
powerape setape e1, target(0.05) nodisplay
powerape trueape e1
check (abs(r(ape) - 0.05) < 1e-6) "empirical route pins the APE"
powerape power e1, n(1500) claim(detect) nsim(100) seed(12) engine(mata) nodisplay
check (r(power) > 0 & r(power) <= 1) "empirical route simulates"
probit y treat age female
powerape dgp f1, fromfit focal(treat) nodisplay
check (abs(r(baseline) - normal(_b[_cons])) < 0.2) "fromfit baseline near the pilot's"
powerape setape f1, target(0.10) nodisplay
powerape power f1, n(800) claim(detect) nsim(100) seed(13) engine(stata) nodisplay
check (r(power) > 0.3) "fromfit route simulates with the Stata engine"
capture noisily powerape dgp f2, fromfit focal(nothere)
check (_rc == 498) "fromfit with a non-regressor focal is refused"

* ---- 6. panel and IV routes --------------------------------------------------------
powerape dgp p1, panel focal(treat, binary p(0.5) icc(1)) covariates((size, normal icc(0.6))) ///
    nperiods(4) rho(0.3) creshare(0.5) baseline(0.30) signal(0.10) nodisplay
powerape setape p1, target(0.08) nodisplay
powerape trueape p1
check (abs(r(ape) - 0.08) < 1e-6) "panel route pins the APE on the ASF"
powerape power p1, n(150) claim(detect) nsim(150) seed(14) engine(mata) nodisplay
local pp = r(power)
check (`pp' > 0.3 & `pp' < 0.9) "panel power in range"
powerape power p1, n(150) claim(detect) nsim(60) seed(14) engine(stata) nodisplay
check (abs(r(power) - `pp') < 0.2) "panel: Stata engine (vce(cluster)) agrees roughly"
check ("`r(fit)'" == "probit, vce(cluster) via engine(stata), effects(closed)") "r(fit) shows vce(cluster) on the panel route"
capture noisily powerape dgp pl, panel model(logit) focal(t, binary p(0.5)) nperiods(4) baseline(0.3)
check (_rc == 498) "panel logit is refused with a teaching error"
powerape dgp p2, panel focal(treat, binary p(0.5) icc(1)) nperiods(6) retention(0.8) rho(0.3) baseline(0.30) nodisplay
powerape setape p2, target(0.08) nodisplay
powerape power p2, n(150) claim(detect) nsim(100) seed(15) engine(mata) nodisplay
check (r(power) > 0 & r(power) <= 1) "unbalanced panel (retention) simulates"
powerape dgp p4, panel focal(treat, binary p(0.5) icc(1)) moderator(m, binary p(0.5) icc(0.7)) ///
    nperiods(4) rho(0.3) creshare(0) baseline(0.30) nodisplay
powerape setaie p4, target(0.05) mainfocal(0.10) mainmoderator(0.05) nodisplay
powerape trueaie p4
check (abs(r(aie) - 0.05) < 1e-6) "panel AIE pinned"
powerape dgp i1, iv focal(training, binary p(0.4)) covariates((size, normal)) instruments((subsidy, binary p(0.5))) ///
    endogeneity(0.4) ivstrength(0.2) baseline(0.30) signal(0.10) nodisplay
powerape setape i1, target(0.10) nodisplay
powerape trueape i1
check (abs(r(ape) - 0.10) < 1e-6) "IV route pins the ASF APE"
powerape power i1, n(1500) claim(detect) nsim(150) seed(16) engine(mata) nodisplay
local pi = r(power)
check ("`r(fit)'" == "Mata engine (Fisher scoring; control-function probit)") "r(fit) names the Mata control-function engine on the IV route"
powerape power i1, n(1500) claim(detect) nsim(40) seed(16) engine(stata) nodisplay
check (abs(r(power) - `pi') < 0.25) "IV: cfprobit engine agrees roughly"
check ("`r(fit)'" == "cfprobit via engine(stata), effects(closed)") "r(fit) names cfprobit on the IV route"
powerape dgp i2, iv focal(rd, normal) instruments((z1, normal)) endogeneity(0.5) ivstrength(0.3) baseline(0.30) nodisplay
powerape setape i2, target(0.08) nodisplay
powerape power i2, n(900) claim(detect) nsim(100) seed(17) engine(mata) nodisplay
check (r(power) > 0.3) "IV continuous focal simulates"
powerape dgp i3, iv focal(training, binary p(0.4)) moderator(m, binary p(0.5)) instruments((z1, normal)) ///
    endogeneity(0.4) ivstrength(0.3) baseline(0.30) nodisplay
powerape setaie i3, target(0.08) mainfocal(0.10) mainmoderator(0.05) nodisplay
powerape trueaie i3
check (abs(r(aie) - 0.08) < 1e-6) "IV AIE pinned"
powerape power i3, n(2000) claim(detect) nsim(60) seed(18) engine(mata) nodisplay
local pa = r(power)
powerape power i3, n(2000) claim(detect) nsim(30) seed(18) engine(stata) nodisplay
check (abs(r(power) - `pa') < 0.3) "IV AIE: cfprobit interact()/mainonly() engine agrees roughly"
local pc = r(power)
* effects(margins) for the IV AIE: margins, dydx(d) at(m = 0 1) after cfprobit with the
* factor-form interaction; per-replication estimates equal the closed form to 1e-7
powerape power i3, n(2000) claim(detect) nsim(30) seed(18) engine(stata) effects(margins) nodisplay
check (abs(r(power) - `pc') <= 1/30 + 1e-9 & r(n_failed) == 0) "IV AIE: effects(margins) agrees with the closed form on the same seed"
check ("`r(effects)'" == "margins" & "`r(fit)'" == "cfprobit via engine(stata), effects(margins)") "r(fit) shows effects(margins) after cfprobit"

* ---- 7. robustness sweep -------------------------------------------------------------
powerape robust t1, n(1500) claim(detect) vary(baseline(0.2 0.4) signal(0 0.3)) nsim(100) seed(19) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("nsc", strofreal(rows(S))); st_local("allfeas", strofreal(all(S[., cols(S)] :== 1)))
check (`nsc' == 9 & `allfeas') "3 x 3 factorial grid, all feasible"
powerape robust t1, n(1500) claim(detect) vary(baseline(0.05 0.95)) gridpoints(3) nsim(100) seed(20) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("nfeas", strofreal(sum(S[., cols(S)])))
check (`nfeas' == 2) "infeasible scenario (baseline .95 with APE .05) reported, not fatal"
powerape robust t1, n(1500) claim(detect) vary(baseline(0.2 0.4)) nsim(100) seed(21) pin(coefficients) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("drift", strofreal(S[1, 2] != S[3, 2]))
check (`drift') "pin(coefficients) lets the implied effect drift"
powerape robust t1, n(1500) claim(detect) vary(baseline(0.2 0.4)) nsim(200) seed(22) nmaxpower(0.8) engine(mata) nodisplay
check (r(nmax) > 0 & r(nmax) < .) "n_max reported"
capture noisily powerape robust t1, n(1500) claim(detect) vary(rho(0.1 0.3)) nsim(50)
check (_rc == 498) "varying a non-varyable input is refused"

* ---- 7b. version check (ado-file vs Mata library) ------------------------------------
powerape version
check (r(match) == 1 & "`r(version)'" != "" & "`r(mlib_version)'" == "`r(ado_version)'") "powerape version: the ado-file and lpowerape.mlib agree"

* ---- 8. session engine default and e() preservation (review items) ---------------
powerape engine, engine(mata)
powerape power t0, n(300) claim(detect) nsim(50) seed(23) nodisplay
check ("`r(engine)'" == "mata" & "`r(fit)'" == "Mata engine (Fisher scoring)") "without engine(), r(fit) reports the session default"
mata: st_local("eng", pape_engine())
check ("`eng'" == "mata") "powerape engine sets a session default that empty options keep"
* a per-command engine() applies to that command only (fixed 2026-10-16): the
* session default must survive a command and a sweep (whose scenario rebuilds
* go through nested -powerape dgp- calls) that name another engine
powerape power t0, n(300) claim(detect) nsim(50) seed(23) engine(stata) nodisplay
check ("`r(engine)'" == "stata" & "`r(fit)'" == "probit via engine(stata), effects(closed)") "r(fit) records the per-command engine, not the session default"
mata: st_local("eng", pape_engine())
check ("`eng'" == "mata") "a per-command engine() does not change the session default"
powerape robust t0, n(300) claim(detect) vary(baseline(0.25 0.40)) gridpoints(2) nsim(30) seed(23) nonmax engine(stata) nodisplay
check ("`r(engine)'" == "stata" & "`r(fit)'" == "probit via engine(stata), effects(closed)") "a sweep stores r(fit) with the per-command engine"
mata: st_local("eng", pape_engine())
check ("`eng'" == "mata") "a sweep with a per-command engine() restores the session default"
powerape engine, engine(stata) effects(closed)
clear
set seed 24
set obs 600
gen treat = runiform() < .5
gen age = rnormal(45, 12)
gen female = runiform() < .55
gen y = runiform() < normal(-0.9 + 0.3*treat + 0.01*age + 0.2*female)
probit y treat age female
local cmd0 `"`e(cmdline)'"'
local ll0 = e(ll)
powerape dgp f3, fromfit focal(treat) nodisplay
powerape setape f3, target(0.10) nodisplay
powerape power f3, n(500) claim(detect) nsim(40) seed(25) engine(stata) nodisplay
check (`"`e(cmdline)'"' == `"`cmd0'"' & abs(e(ll) - `ll0') < 1e-12) "the Stata engine leaves the user's e() untouched"
powerape robust f3, n(500) claim(detect) vary(baseline(0.25 0.40)) nsim(40) seed(26) nonmax engine(stata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("allfeas", strofreal(all(S[., cols(S)] :== 1)))
check (`allfeas') "sweep on a pilot-model DGP with the Stata engine: all scenarios feasible"
powerape robust e1, n(1500) claim(detect) vary(baseline(0.2 0.4) signal(0 0.3)) nsim(40) seed(27) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("allfeas", strofreal(all(S[., cols(S)] :== 1) & rows(S) == 9))
check (`allfeas') "sweep on an empirical DGP recalibrates from the stored rows"

powerape drop _all
di as txt _n "run-tests: $PAPE_NPASS PASS, $PAPE_NFAIL FAIL"
if $PAPE_NFAIL > 0 {
    di as err "some tests failed"
    exit 9
}
exit, clear
