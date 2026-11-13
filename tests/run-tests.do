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

* ---- 9. difference-in-differences designs -------------------------------------
clear
* design, pin, readings (closed forms without covariates)
powerape dgp k1, did model(logit) baseline(0.20) gap(0.20) trend(0.10) rho(0.4) nodisplay
check ("`r(route)'" == "did" & "`r(design)'" == "panel" & "`r(trends)'" == "index" & "`r(analysis)'" == "index") "did: route, design, trends and analysis are returned"
local cf = invlogit(logit(0.40) + logit(0.30) - logit(0.20))
check (abs(el(r(rates), 1, 4) - `cf') < 1e-12) "did: index world, the untreated rate follows the control group on the logit index"
powerape setatt k1, target(0.07) nodisplay
check (abs(r(target) - 0.07) < 1e-12 & abs(el(r(rates), 1, 5) - (`cf' + 0.07)) < 1e-12 & "`r(estimand)'" == "att") "did: setatt pins the effect on the treated"
check (abs(r(reading_index) - 0.07) < 1e-10 & abs(r(reading_linear) - (`cf' + 0.07 - 0.50)) < 1e-12) "did: both readings of the four rates"
powerape trueatt k1
check (abs(r(att) - 0.07) < 1e-9) "did: true ATT equals the target"
powerape dgp k1l, did model(logit) baseline(0.20) gap(0.20) trend(0.10) trends(linear) analysis(linear) nodisplay
powerape setatt k1l, target(0.10) nodisplay
check (abs(el(r(rates), 1, 4) - 0.50) < 1e-12 & abs(r(reading_linear) - 0.10) < 1e-12) "did: linear world, the untreated rate follows the control group in probabilities"
check (abs(r(reading_index) - (0.60 - `cf')) < 1e-10 & abs(r(reading_index) - 0.0667) < 5e-4) "did: rates .20/.30/.40/.60 read on the logit index: .067"
powerape dgp k1p, did model(probit) baseline(0.20) gap(0.20) trend(0.10) trends(linear) nodisplay
powerape setatt k1p, target(0.10) nodisplay
check (abs(r(reading_index) - (0.60 - normal(invnormal(0.40) + invnormal(0.30) - invnormal(0.20)))) < 1e-10 & abs(r(reading_index) - 0.074) < 1e-3) "did: the same rates read on the probit index: .074"
* without a baseline gap, or without a change in the control group, the two readings coincide
foreach lk in probit logit {
    powerape dgp kz, did model(`lk') baseline(0.20) gap(0) trend(0.10) trends(linear) design(cross-sections) nodisplay replace
    powerape setatt kz, target(0.10) nodisplay
    check (abs(r(reading_index) - 0.10) < 1e-10 & abs(r(reading_linear) - 0.10) < 1e-10) "did: no baseline gap, the readings coincide (`lk')"
    powerape dgp kz, did model(`lk') baseline(0.20) gap(0.20) trend(0) trends(linear) design(cross-sections) nodisplay replace
    powerape setatt kz, target(0.10) nodisplay
    check (abs(r(reading_index) - 0.10) < 1e-10 & abs(r(reading_linear) - 0.10) < 1e-10) "did: no control-group change, the readings coincide (`lk')"
}
* refusals of the builder and the pins
capture noisily powerape dgp bad, did baseline(0.90) gap(0.20)
check (_rc == 498) "did: baseline + gap outside (0, 1) is refused"
capture noisily powerape dgp bad, did baseline(0.60) gap(0.30) trend(0.20) trends(linear)
check (_rc == 498) "did: linear trends that leave (0, 1) are refused"
capture noisily powerape dgp bad, did baseline(0.30 0.20)
check (_rc == 498) "did: two baselines need subgroup()"
capture noisily powerape dgp bad, did baseline(0.30) design(wide)
check (_rc == 498) "did: unknown design() is refused"
capture noisily powerape dgp bad, did baseline(0.30) focal(t, binary p(0.5))
check (_rc == 198) "did: focal() is refused"
capture noisily powerape setatt k1, target(0.60)
check (_rc == 498) "did: infeasible ATT is refused"
capture noisily powerape setatt k1, target(0.05 0.02)
check (_rc == 498) "did: two targets need a subgroup design"
capture noisily powerape setatt t0, target(0.05)
check (_rc == 498) "setatt is refused for a design that is not difference-in-differences"
capture noisily powerape setape k1, target(0.05)
check (_rc == 498) "setape is refused for a difference-in-differences design"
capture noisily powerape trueape k1
check (_rc == 498) "trueape is refused for a difference-in-differences design"
powerape dgp k0, did model(logit) baseline(0.20) gap(0.20) trend(0.10) nodisplay
capture noisily powerape power k0, n(500) claim(detect) nsim(20)
check (_rc == 498) "did: power needs a pinned effect"
* several periods and covariates: the coefficient vector and the pin
powerape dgp k3, did model(probit) baseline(0.25) gap(0.10) trend(0.08) design(cross-sections) npre(2) npost(3) ///
    covariates((z, normal)) signal(0.2) nint(50000) nodisplay
powerape setatt k3, target(0.08) nodisplay
mata: bt = pape_did_beta(pape_get_dgp("k3")); st_local("okb", strofreal(rows(bt) == 8 & bt[3] == 0 & bt[4] == bt[5] & bt[5] == bt[6]))
check (`okb') "did: 2 + 3 periods, coefficients 0 on the pre-period dummy and equal on the post-period dummies"
check (abs(r(p_post) - 0.6) < 1e-12 & r(n_pre) == 2 & r(n_post) == 3) "did: the default share after treatment is n_post / periods"
powerape trueatt k3
check (abs(r(att) - 0.08) < 1e-7) "did: true ATT with covariates equals the target"
* a binary subgroup
powerape dgp ks, did model(logit) baseline(0.30 0.20) gap(0.10 0.15) trend(0.10 0.05) rho(0.5) subgroup(female, p(0.4)) nodisplay
capture noisily powerape setatt ks, target(0.05)
check (_rc == 498) "did: a subgroup design needs the two effects"
powerape setatt ks, target(0.10 0.03) nodisplay
check (abs(r(target) + 0.07) < 1e-12 & "`r(estimand)'" == "datt" & abs(r(att_ref) - 0.10) < 1e-12 & abs(r(att_sub) - 0.03) < 1e-12) "did: the subgroup estimand is the difference of the two effects"
powerape trueatt ks
check (abs(r(att) + 0.07) < 1e-9 & abs(r(att_ref) - 0.10) < 1e-9 & abs(r(att_sub) - 0.03) < 1e-9) "did: true subgroup effects and their difference"
mata: st_local("okb", strofreal(rows(pape_did_beta(pape_get_dgp("ks"))) == 8))
check (`okb') "did: the subgroup design matrix has eight columns with two periods"

* the saturated identities on one fixed study (both links, both analyses, both engines)
clear
set seed 31
set obs 1200
generate double d = runiform() < 0.5
generate double post = runiform() < 0.5
generate double w = d * post
generate double y = runiform() < invlogit(-1 + 0.5 * d + 0.4 * post + 0.6 * w)
forvalues a = 0/1 {
    forvalues b = 0/1 {
        quietly summarize y if d == `a' & post == `b'
        local m`a'`b' = r(mean)
    }
}
mata: X = (J(st_nobs(), 1, 1), st_data(., ("d", "post", "w"))); y = st_data(., "y")
mata: est = .; se = .; b = .; jac = .; by = .
foreach lk in logit probit {
    local sat = cond("`lk'" == "logit", `m11' - invlogit(logit(`m10') + logit(`m01') - logit(`m00')), ///
        `m11' - normal(invnormal(`m10') + invnormal(`m01') - invnormal(`m00')))
    mata: ok = pape_did_analyze(X, y, J(0, 1, .), J(0, 1, .), "`lk'", 2, "index", "model", "mata", "closed", est, se, b, jac, by)
    mata: st_local("okb", strofreal(ok & abs(est - `sat') < 1e-7))
    check (`okb') "did: saturated identity of the ATT estimate, `lk', engine(mata)"
    mata: e_m = est; s_m = se
    mata: ok = pape_did_analyze(X, y, J(0, 1, .), J(0, 1, .), "`lk'", 2, "index", "model", "stata", "closed", est, se, b, jac, by)
    mata: st_local("okb", strofreal(ok & abs(est - `sat') < 1e-6 & abs(se / s_m - 1) < 1e-5))
    check (`okb') "did: saturated design, `lk': Stata's fit gives the same estimate and SE"
    mata: ok = pape_did_analyze(X, y, J(0, 1, .), J(0, 1, .), "`lk'", 2, "index", "model", "stata", "margins", est, se, b, jac, by)
    mata: st_local("okb", strofreal(ok & abs(est - `sat') < 1e-6 & abs(se / s_m - 1) < 1e-5))
    check (`okb') "did: saturated design, `lk': margins, dydx() subpop() gives the same estimate and SE"
}
local cd = `m11' - `m10' - `m01' + `m00'
mata: ok = pape_did_analyze(X, y, J(0, 1, .), J(0, 1, .), "logit", 2, "linear", "model", "mata", "closed", est, se, b, jac, by)
mata: st_local("okb", strofreal(ok & abs(est - `cd') < 1e-12))
check (`okb') "did: the linear analysis is the cross difference of the four cell means"
quietly regress y d post w, vce(robust)
scalar sb = _b[w]
scalar sse = _se[w]
mata: st_local("okb", strofreal(abs(est - st_numscalar("sb")) < 1e-12 & abs(se - st_numscalar("sse")) < 1e-12))
check (`okb') "did: OLS coefficient and HC1 standard error equal regress, vce(robust)"
clear

* both engines on one seed (the same draws, the same decisions for logit and OLS)
powerape power k1, n(1200) claim(detect) nsim(60) seed(41) engine(mata) nodisplay
scalar pm = r(power)
check ("`r(fit)'" == "Mata engine (Fisher scoring)" & r(valid) == 1 & "`r(estimand)'" == "att") "did: r(fit), r(valid) and r(estimand) under engine(mata)"
check (abs(r(est_mean) - 0.07) < 0.02 & r(est_mcse) < 0.01) "did: r(est_mean) and r(est_mcse) are returned"
powerape power k1, n(1200) claim(detect) nsim(60) seed(41) engine(stata) nodisplay
check (r(power) == scalar(pm) & "`r(fit)'" == "logit, vce(cluster) via engine(stata), effects(closed)") "did: logit, vce(cluster) in the loop gives the same power on the same seed"
powerape power k1, n(1200) claim(detect) nsim(60) seed(41) engine(stata) effects(margins) nodisplay
check (r(power) == scalar(pm) & "`r(fit)'" == "logit, vce(cluster) via engine(stata), effects(margins)") "did: margins in the loop gives the same power on the same seed"
powerape power k1l, n(1200) claim(detect) nsim(60) seed(42) engine(mata) nodisplay
scalar pm = r(power)
check ("`r(fit)'" == "Mata engine (OLS)" & "`r(analysis)'" == "linear") "did: linear analysis under engine(mata)"
powerape power k1l, n(1200) claim(detect) nsim(60) seed(42) engine(stata) nodisplay
check (r(power) == scalar(pm) & "`r(fit)'" == "regress, vce(cluster) via engine(stata)") "did: regress, vce(cluster) in the loop gives the same power on the same seed"
powerape power k3, n(1500) claim(detect) nsim(40) seed(43) se(robust) engine(stata) nodisplay
check ("`r(fit)'" == "probit, vce(robust) via engine(stata), effects(closed)" & "`r(se)'" == "robust") "did: cross-sections accept se(robust)"
powerape power k1, n(1200) claim(detect) nsim(40) seed(44) se(robust) engine(mata) nodisplay
check ("`r(se)'" == "model") "did: se() is fixed in a panel (unit-clustered)"
powerape power ks, n(2500) claim(detect) nsim(40) seed(45) engine(mata) nodisplay
scalar pm = r(power)
check ("`r(estimand)'" == "datt" & r(target) < 0) "did: subgroup design runs, estimand datt"
powerape power ks, n(2500) claim(detect) nsim(40) seed(45) engine(stata) nodisplay
scalar ps = r(power)
powerape power ks, n(2500) claim(detect) nsim(40) seed(45) engine(stata) effects(margins) nodisplay
check (scalar(ps) == scalar(pm) & r(power) == scalar(pm)) "did: subgroup design, three ways of fitting agree on one seed"
* an old route returns none of the new results
powerape power t0, n(300) claim(detect) nsim(30) seed(46) engine(mata) nodisplay
check ("`r(est_mean)'" == "" & "`r(analysis)'" == "" & "`r(valid)'" == "") "a standard design returns no difference-in-differences results"

* an analysis that assumes another scale than the world has
powerape dgp km, did model(logit) baseline(0.20) gap(0.20) trend(0.10) rho(0.4) trends(linear) analysis(index) nodisplay
powerape setatt km, target(0.10) nodisplay
powerape power km, n(1500) claim(detect) nsim(200) seed(47) engine(mata) nodisplay
check (r(valid) == 0 & abs(r(est_mean) - 0.0667) < 0.012 & r(est_mean) < 0.09) "did: logit analysis in a linear world centres on the index reading (.067), not on the ATT (.100)"
capture noisily powerape n km, power(0.8) claim(detect) nsim(50)
check (_rc == 498) "did: n refuses an analysis that is not valid in the simulated world"
capture noisily powerape mde km, n(1500) claim(detect) nsim(50)
check (_rc == 498) "did: mde refuses it too"
capture noisily powerape powerat, alpha(0.01)
check (_rc == 498) "did: powerat refuses a rejection rate"
capture noisily powerape alpha
check (_rc == 498) "did: alpha refuses a rejection rate"
capture noisily powerape statement
check (_rc == 498) "did: statement refuses a rejection rate"
capture noisily powerape robust km, n(1500) claim(detect) vary(rho(0.2 0.6))
check (_rc == 498) "did: robust refuses a base design whose analysis is not valid"
powerape power k1l, n(900) claim(detect) nsim(300) seed(48) engine(mata) nodisplay
check (abs(r(est_mean) - 0.10) < 4 * r(est_mcse) + 0.002) "did: the linear analysis centres on the ATT in a linear world"

* curve, n, mde
powerape curve k1, n(800 1500) claim(detect) nsim(60) seed(49) engine(mata) nodisplay
check (colsof(r(results)) == 8 & el(r(results), 2, 3) >= el(r(results), 1, 3) - 0.15) "did: the curve carries the mean estimate in column 8"
powerape curve t0, n(300 500) claim(detect) nsim(40) seed(49) engine(mata) nodisplay
check (colsof(r(results)) == 7) "a standard design's curve keeps seven columns"
powerape n k1l, power(0.8) claim(detect) nsim(200) seed(50) engine(mata) noconfirm nodisplay
check (r(n) > 800 & r(n) < 2200 & "`r(kind)'" == "n") "did: required number of units, linear analysis (closed form about 1,350)"
powerape mde k1, n(1500) claim(detect) nsim(200) seed(51) engine(mata) noconfirm nodisplay
check (r(mde) > 0.07 & r(mde) < 0.13 & "`r(estimand)'" == "att") "did: minimum detectable ATT"
powerape mde ks, n(2500) claim(detect) nsim(100) seed(52) engine(mata) noconfirm nodisplay
check (r(mde) < 0 & "`r(estimand)'" == "datt") "did: the MDE of a subgroup design is the difference, in the pinned direction"
powerape dgp ks0, did model(logit) baseline(0.30 0.20) gap(0.10 0.15) trend(0.10 0.05) rho(0.5) subgroup(female, p(0.4)) nodisplay
capture noisily powerape mde ks0, n(2500) claim(detect) nsim(50)
check (_rc == 498) "did: the MDE of a subgroup design needs the reference effect pinned"

* robustness sweep
powerape robust k1, n(1200) claim(detect) vary(rho(0.2 0.6)) nsim(60) seed(1) nonmax engine(mata) nodisplay
matrix S = r(scenarios)
check (rowsof(S) == 3 & colsof(S) == 6 & abs(S[1, 2] - 0.07) < 1e-12 & abs(S[3, 2] - 0.07) < 1e-12) "did: sweep over rho re-pins the effect in every scenario"
powerape power k1, n(1200) claim(detect) nsim(60) seed(3) engine(mata) nodisplay
check (r(power) == S[2, 3]) "did: the middle scenario equals a direct run at its seed"
powerape robust k1, n(1200) claim(detect) vary(baseline(0.15 0.35) gap(0 0.2) trend(0.04 0.12) ptreated(0.3 0.5)) ///
    gridpoints(2) nsim(30) seed(2) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("okb", strofreal(rows(S) == 16 & all(S[., cols(S)] :== 1)))
check (`okb') "did: baseline, gap, trend and the treated share can be varied"
powerape dgp kx, did model(logit) baseline(0.25) gap(0.10) trend(0.08) design(cross-sections) nodisplay
powerape setatt kx, target(0.07) nodisplay
capture noisily powerape robust kx, n(1200) claim(detect) vary(rho(0.2 0.6))
check (_rc == 498) "did: rho cannot be varied in cross-sections"
capture noisily powerape robust k1, n(1200) claim(detect) vary(rho(0.2 0.6)) pin(coefficients)
check (_rc == 498) "did: pin(coefficients) is refused"
capture noisily powerape robust ks, n(2500) claim(detect) vary(baseline(0.2 0.3))
check (_rc == 498) "did: baseline cannot be varied in a subgroup design"
capture noisily powerape robust k1, n(1200) claim(detect) vary(trends(index probit))
check (_rc == 498) "did: vary(trends()) takes index and linear"
powerape robust k1, n(1500) claim(detect) vary(trends(index linear)) nsim(150) seed(2) nmaxpower(0.8) engine(mata) nodisplay
matrix S = r(scenarios)
check (colsof(S) == 8 & S[1, 1] == 1 & S[2, 1] == 2 & S[1, 8] == 1 & S[2, 8] == 0) "did: a sweep over trends flags the scenario the analysis is biased in"
check (abs(S[1, 2] - 0.07) < 1e-12 & abs(S[2, 2] - 0.07) < 1e-12 & S[2, 7] < S[1, 7]) "did: the effect is held and the biased scenario centres lower"
check (r(worst) == 1 & r(nmax_scenario) == 1 & r(nmax) > 1500) "did: the biased scenario is shown but never ranked"
powerape robust k1, n(1500) claim(detect) vary(baseline(0.15 0.35)) mode(mde) gridpoints(2) nsim(80) seed(3) engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("okb", strofreal(all(S[., 2] :> 0.04) & all(S[., 2] :< 0.2)))
check (`okb') "did: MDE mode of the sweep"
powerape robust ks, n(2500) claim(detect) vary(rho(0.3 0.7)) gridpoints(2) nsim(40) seed(4) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("okb", strofreal(all(abs(S[., 2] :+ 0.07) :< 1e-12)))
check (`okb') "did: sweep on a subgroup design holds both effects"
powerape robust k1l, n(900) claim(detect) vary(gap(0 0.2)) gridpoints(2) nsim(40) seed(5) nonmax engine(mata) nodisplay
check ("`r(analysis)'" == "linear" & "`r(fit)'" == "Mata engine (OLS)") "did: the sweep carries the linear analysis into every scenario"

* justified alpha and the methods paragraph
powerape power k1, n(1500) claim(minimum) sesoi(0.02) nsim(300) seed(5) engine(mata) nodisplay
local pw = r(power)
powerape powerat, alpha(0.01 0.05)
check (abs(el(r(table), 2, 3) - `pw') < 1e-12 & el(r(table), 1, 3) <= el(r(table), 2, 3)) "did: powerat re-reads the stored draws"
powerape alpha, cost(4) size nsimsize(300) seed(6) engine(mata) nodisplay
check (abs(r(size_truth) - 0.02) < 1e-12 & abs(r(size) - r(alpha)) < 4 * r(size_mcse) + 0.01) "did: the realized size is measured in the boundary world (true ATT = sesoi)"
powerape statement _last
local st `"`r(statement)'"'
check (strpos(`"`st'"', "average treatment effect on the treated (ATT) in a difference-in-differences design") > 0) "did statement: names the estimand"
check (strpos(`"`st'"', "a treated group (50% of the units) and a control group") > 0 & strpos(`"`st'"', "observed in one period before and one period after one common treatment date") > 0) "did statement: groups and periods"
check (strpos(`"`st'"', "as a panel of the same units, with a latent within-unit correlation of 0.40") > 0) "did statement: data structure"
check (strpos(`"`st'"', "untreated trends parallel on the logit index") > 0 & strpos(`"`st'"', "The assumed true ATT was 0.070 (7.0 percentage points).") > 0) "did statement: scale of the trends and the effect"
check (strpos(`"`st'"', "n = 1500 units observed in 2 periods (3,000 unit-period observations)") > 0) "did statement: sample size"
check (strpos(`"`st'"', "a pooled logit regression of the outcome on group, period indicators and the treatment indicator") > 0 & strpos(`"`st'"', "unit-clustered standard errors") > 0) "did statement: the analysis"
check (strpos(`"`st'"', "a group-by-before/after cell with no events or no non-events") > 0 & strpos(`"`st'"', "focal") == 0) "did statement: separation named by its cells, no focal variable"
powerape statement _lastalpha
check (strpos(`"`r(statement)'"', "the planning value (ATT = 0.070)") > 0 & strpos(`"`r(statement)'"', "Each simulated study was analysed by a pooled logit regression") > 0 & strpos(`"`r(statement)'"', "At n = 1500 units the optimum") > 0) "did: the justification paragraph names the estimand and the analysis"
powerape dgp kxl, did model(probit) baseline(0.20) gap(0.20) trend(0.10) design(cross-sections) npre(2) npost(3) ///
    trends(linear) analysis(linear) nodisplay
powerape setatt kxl, target(0.08) nodisplay
powerape power kxl, n(2000) claim(detect) nsim(60) seed(7) engine(mata) nodisplay
powerape statement
local st `"`r(statement)'"'
check (strpos(`"`st'"', "observed in 2 periods before and 3 periods after") > 0 & strpos(`"`st'"', "as repeated cross-sections, with 60% of the observations after the treatment date") > 0) "did statement: several periods, cross-sections"
check (strpos(`"`st'"', "untreated trends parallel in probabilities") > 0 & strpos(`"`st'"', "n = 2000 observations") > 0) "did statement: linear world"
check (strpos(`"`st'"', "ordinary least squares of the outcome on group, period indicators and the treatment indicator (a linear probability model)") > 0 & strpos(`"`st'"', "heteroskedasticity-robust (HC1) standard errors") > 0) "did statement: the linear analysis"
check (strpos(`"`st'"', "maximum-likelihood") == 0 & strpos(`"`st'"', "separation") == 0) "did statement: no separation language for OLS"
powerape power ks, n(2500) claim(detect) nsim(40) seed(8) engine(mata) nodisplay
powerape statement
local st `"`r(statement)'"'
check (strpos(`"`st'"', "difference between two subgroups' average treatment effects on the treated (by female)") > 0 & strpos(`"`st'"', "were 0.100 (female = 0) and 0.030 (female = 1), a true difference of -0.070") > 0) "did statement: subgroup difference"
check (strpos(`"`st'"', "fully interacted with female") > 0 & strpos(`"`st'"', "within a subgroup") > 0) "did statement: the interacted analysis"
powerape power t0, n(300) claim(detect) nsim(40) seed(9) engine(mata) nodisplay
powerape statement
check (strpos(`"`r(statement)'"', "average partial effect (APE) of t") > 0 & strpos(`"`r(statement)'"', "difference-in-differences") == 0) "a standard design keeps its paragraph"

* a design from pre-treatment pilot data (20,000 units drawn from the route's copula)
clear
set seed 61
set obs 20000
generate long firm = _n
generate double a = rnormal()
generate byte treated = runiform() < 0.4
expand 2
bysort firm: generate byte year = _n
generate double p = cond(treated, 0.34, 0.22)
generate byte y = runiform() < normal((invnormal(p) - sqrt(0.35) * a) / sqrt(1 - 0.35))
powerape dgp kp, didpilot outcome(y) group(treated) id(firm) period(year) trend(0.05) model(logit) nodisplay
matrix PL = r(pilot)
check (abs(r(baseline) - 0.22) < 0.01 & abs(el(r(rates), 1, 3) - 0.34) < 0.012 & abs(r(p_treated) - 0.40) < 0.012) "didpilot: recovers the rates and the treated share"
check (abs(r(rho) - 0.35) < 0.04 & "`r(design)'" == "panel" & PL[1, 4] == 20000) "didpilot: recovers the latent within-unit correlation"
check (abs(el(r(rates), 1, 2) - (r(baseline) + 0.05)) < 1e-12) "didpilot: trend() is the control group's change"
powerape setatt kp, target(0.06) nodisplay
powerape robust kp, n(1500) claim(detect) vary(rho(0.2 0.5)) gridpoints(2) nsim(30) seed(10) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)"); st_local("okb", strofreal(all(S[., cols(S)] :== 1)))
check (`okb') "didpilot: a pilot-built design can be swept"
capture noisily powerape dgp kp2, didpilot outcome(y) group(treated) id(firm) period(year)
check (_rc == 498) "didpilot: trend() is required"
powerape dgp kp3, didpilot outcome(y) group(treated) trend(0.05) baseline(0.25) nodisplay
check ("`r(design)'" == "cross-sections" & abs(r(baseline) - 0.25) < 1e-12 & r(rho) == 0) "didpilot: no id() gives cross-sections; a named input overrides the estimate"
capture noisily powerape dgp kp4, didpilot outcome(y) group(treated) id(firm) trend(0.05) subgroup(female, p(0.4))
check (_rc == 198) "didpilot: subgroup designs are refused"
quietly keep if year == 1
capture noisily powerape dgp kp5, didpilot outcome(y) group(treated) id(firm) trend(0.05)
check (_rc == 498) "didpilot: one period per unit needs rho()"
powerape dgp kp5, didpilot outcome(y) group(treated) id(firm) trend(0.05) rho(0.3) nodisplay
check (abs(r(rho) - 0.3) < 1e-12) "didpilot: rho() supplied for a one-period pilot"
clear

powerape drop _all
di as txt _n "run-tests: $PAPE_NPASS PASS, $PAPE_NFAIL FAIL"
if $PAPE_NFAIL > 0 {
    di as err "some tests failed"
    exit 9
}
exit, clear
