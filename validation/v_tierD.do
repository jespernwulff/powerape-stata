* v_tierD.do -- stochastic parity with the R package on matched designs.
* R reference values (validation/tierD_r.R, powerape 1.11.3, 4,000
* replications per power) are compared with the Stata port at 4,000
* replications (engine(mata)) and, for a subset, 1,000 replications of the
* Stata engine (probit/logit/cfprobit in the loop). RNG streams differ, so
* the criterion is |difference| <= 3 x sqrt(MCSE_R^2 + MCSE_Stata^2).
* Run from the repository root (the folder that holds ado/ and validation/):
*     do validation/v_tierD.do
* Writes validation/tierD-results.txt and the log validation/v_tierD.log.
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "v_tierD.do: run it from the repository root:  do validation/v_tierD.do"
    exit 601
}
adopath ++ "`root'/ado"
mata: mata mlib index
powerape version
if r(match) != 1 exit 459
local pav "`r(version)'"
cd "`root'/validation"
* the log opens after the set-up, whose output names local folders
capture log close tier
log using v_tierD.log, text replace nomsg name(tier)
di as txt "powerape `pav'; Stata `c(stata_version)' (`c(born_date)')"
import delimited data/tierD-reference.csv, clear varnames(1) stringcols(1 2) asdouble
capture destring value mcse, replace force
mata:
RD = st_sdata(., "design"); RQ = st_sdata(., "quantity"); RV = st_data(., "value"); RM = st_data(., "mcse")
real rowvector refd(string scalar d, string scalar q)
{
    external string colvector RD, RQ
    external real colvector RV, RM
    real scalar j
    for (j = 1; j <= rows(RD); j++) if (RD[j] == d & RQ[j] == q) return((RV[j], RM[j]))
    return((., .))
}
RES = J(0, 1, ""); NP = 0; NF = 0
void recd(string scalar label, real scalar s, real scalar s_mcse, real rowvector r, real scalar k)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar tol, ok, rm
    rm = (r[2] == . ? 0 : r[2])
    tol = k * sqrt(rm^2 + (s_mcse == . ? 0 : s_mcse)^2)
    ok = (abs(s - r[1]) <= tol)
    RES = RES \ sprintf("%-58s Stata %8.4f (MCSE %6.4f)  R %8.4f (MCSE %6.4f)  |diff| %6.4f  tol %6.4f  %s", label, s, s_mcse, r[1], rm, abs(s - r[1]), tol, (ok ? "PASS" : "FAIL"))
    if (ok) NP++
    else NF++
}
void recr(string scalar label, real scalar s, real scalar r, real scalar reltol)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar ok
    ok = (abs(s - r) <= reltol * max((abs(r), 1e-8)))
    RES = RES \ sprintf("%-58s Stata %8.4f  R %8.4f  rel.diff %8.2e  tol %6.1e  %s", label, s, r, abs(s - r) / max((abs(r), 1e-8)), reltol, (ok ? "PASS" : "FAIL"))
    if (ok) NP++
    else NF++
}
end
local NS 4000
local NSS 1000

* D1: Design A world
powerape dgp d1, focal(treat, binary p(0.5)) covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.15) nodisplay
powerape setape d1, target(0.05) nodisplay
powerape power d1, n(2000) claim(minimum) sesoi(0.02) nsim(`NS') seed(201) engine(mata) nodisplay
mata: recd("D1 minimum power, n 2000, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D1", "power"), 3)
mata: o = st_matrix("r(outcomes)")
mata: recd("D1 detect-only share", o[2], sqrt(o[2] * (1 - o[2]) / `NS'), (refd("D1", "detect_only")[1], sqrt(refd("D1", "detect_only")[1] * (1 - refd("D1", "detect_only")[1]) / 4000)), 3)
mata: recd("D1 inconclusive share", o[3], sqrt(o[3] * (1 - o[3]) / `NS'), (refd("D1", "inconclusive")[1], sqrt(refd("D1", "inconclusive")[1] * (1 - refd("D1", "inconclusive")[1]) / 4000)), 3)
powerape power d1, n(2000) claim(minimum) sesoi(0.02) nsim(`NSS') seed(202) engine(stata) nodisplay
mata: recd("D1 minimum power, engine(stata) probit, nsim 1000", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D1", "power"), 3)
powerape power d1, n(2000) claim(detect) nsim(`NS') seed(203) engine(mata) nodisplay
mata: recd("D1 detect power, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D1", "power_detect"), 3)

* D2: logit, continuous focal, robust and model SEs
powerape dgp d3, model(logit) focal(dose, normal mean(2) sd(1)) covariates((z1, normal) (z2, binary p(0.4))) ///
    correlation(0.2) baseline(0.25) signal(0.10) nodisplay
powerape setape d3, target(0.03) nodisplay
powerape power d3, n(800) claim(detect) nsim(`NS') seed(204) se(robust) engine(mata) nodisplay
mata: recd("D2 logit continuous, robust SE, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D2", "power"), 3)
powerape power d3, n(800) claim(detect) nsim(`NS') seed(205) engine(mata) nodisplay
mata: recd("D2 logit continuous, model SE, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D2", "power_model"), 3)
powerape power d3, n(800) claim(detect) nsim(`NSS') seed(206) se(robust) engine(stata) nodisplay
mata: recd("D2 logit continuous, vce(robust), engine(stata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D2", "power"), 3)

* D3: AIE
powerape dgp d2, focal(treat, binary p(0.5)) moderator(female, binary p(0.55)) baseline(0.30) nodisplay
powerape setaie d2, target(0.04) mainfocal(0.10) mainmoderator(0.05) nodisplay
powerape power d2, n(4000) claim(detect) nsim(`NS') seed(207) engine(mata) nodisplay
mata: recd("D3 AIE binary x binary, n 4000, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D3", "power"), 3)
powerape power d2, n(4000) claim(detect) nsim(`NSS') seed(208) engine(stata) nodisplay
mata: recd("D3 AIE binary x binary, engine(stata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D3", "power"), 3)
powerape power d2, n(4000) claim(detect) nsim(300) seed(209) engine(stata) effects(margins) nodisplay
mata: recd("D3 AIE, engine(stata) effects(margins)=ginteff, nsim 300", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D3", "power"), 3)
powerape dgp d2b, model(logit) focal(x, normal) moderator(m, binary p(0.5)) covariates((z, normal)) baseline(0.35) signal(0.1) nodisplay
powerape setaie d2b, target(0.03) mainfocal(0.05) mainmoderator(0.08) nodisplay
powerape power d2b, n(3000) claim(detect) nsim(`NS') seed(210) engine(mata) nodisplay
mata: recd("D3b AIE continuous x binary (logit), engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D3b", "power"), 3)

* D4: equivalence
powerape dgp d4, focal(t, binary p(0.5)) baseline(0.30) nodisplay
powerape setape d4, target(0) nodisplay
powerape power d4, n(2000) claim(equivalence) sesoi(0.05) nsim(`NS') seed(211) engine(mata) nodisplay
mata: recd("D4 equivalence power at truth 0, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D4", "power"), 3)

* D5: panel
powerape dgp dp, panel focal(treat, binary p(0.5) icc(1)) covariates((size, normal icc(0.6))) ///
    nperiods(4) rho(0.3) creshare(0.5) baseline(0.30) signal(0.10) nodisplay
powerape setape dp, target(0.08) nodisplay
powerape power dp, n(150) claim(detect) nsim(`NS') seed(212) engine(mata) nodisplay
mata: recd("D5 panel CRE probit, 150 units x 4, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D5", "power"), 3)
powerape power dp, n(150) claim(detect) nsim(`NSS') seed(213) engine(stata) nodisplay
mata: recd("D5 panel, engine(stata) probit vce(cluster)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D5", "power"), 3)
powerape dgp dpu, panel focal(treat, binary p(0.5) icc(1)) covariates((size, normal icc(0.6))) ///
    nperiods(6) retention(0.8) rho(0.3) creshare(0.5) baseline(0.30) signal(0.10) nodisplay
powerape setape dpu, target(0.08) nodisplay
powerape power dpu, n(150) claim(detect) nsim(`NS') seed(214) engine(mata) nodisplay
mata: recd("D5u unbalanced panel (retention .8), engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D5u", "power"), 3)

* D6: IV
powerape dgp di, iv focal(training, binary p(0.4)) covariates((size, normal)) instruments((subsidy, binary p(0.5))) ///
    endogeneity(0.4) ivstrength(0.2) baseline(0.30) signal(0.10) nodisplay
powerape setape di, target(0.10) nodisplay
powerape power di, n(1500) claim(detect) nsim(`NS') seed(215) engine(mata) nodisplay
mata: recd("D6 IV binary focal (probit first stage), engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D6", "power"), 3)
powerape power di, n(1500) claim(detect) nsim(400) seed(216) engine(stata) nodisplay
mata: recd("D6 IV binary focal, engine(stata) cfprobit, nsim 400", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D6", "power"), 3)
powerape dgp di2, iv focal(rd, normal mean(0) sd(1)) instruments((z1, normal)) endogeneity(0.5) ivstrength(0.3) baseline(0.30) nodisplay
powerape setape di2, target(0.08) nodisplay
powerape power di2, n(900) claim(detect) nsim(`NS') seed(217) engine(mata) nodisplay
mata: recd("D6c IV continuous focal (linear first stage), engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D6c", "power"), 3)
powerape power di2, n(900) claim(detect) nsim(400) seed(218) engine(stata) nodisplay
mata: recd("D6c IV continuous, engine(stata) cfprobit, nsim 400", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D6c", "power"), 3)
powerape power di2, n(900) claim(detect) nsim(200) seed(219) engine(stata) effects(margins) nodisplay
mata: recd("D6c IV continuous, cfprobit + margins, nsim 200", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D6c", "power"), 3)

* D7: empirical route on the shared pilot file (calibration is deterministic here)
import delimited data/pilot600.csv, clear varnames(1) asdouble
powerape dgp de, empirical focal(treat) covariates(age female) baseline(0.30) signal(0.15) nodisplay
powerape setape de, target(0.05) nodisplay
mata: recr("D7 empirical beta0 (exact: same rows)", st_numscalar("r(beta0)"), refd("D7", "beta0")[1], 1e-7)
mata: recr("D7 empirical gamma1", st_matrix("r(gamma)")[1], refd("D7", "gamma1")[1], 1e-7)
mata: recr("D7 empirical gamma2", st_matrix("r(gamma)")[2], refd("D7", "gamma2")[1], 1e-7)
mata: recr("D7 empirical beta_focal", st_numscalar("r(beta_focal)"), refd("D7", "beta_focal")[1], 1e-7)
powerape power de, n(1500) claim(minimum) sesoi(0.02) nsim(`NS') seed(220) engine(mata) nodisplay
mata: recd("D7 empirical route minimum power, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D7", "power"), 3)
clear

* D8: required n (answers differ by search noise; each must satisfy the confirmation rule,
*     and the two answers must lie within 6% of each other)
powerape n d1, power(0.8) claim(detect) nsim(1500) seed(221) engine(mata) nodisplay
mata: recr("D8 required n, detect 80%, Design A world (R 2494)", st_numscalar("r(n)"), refd("D8", "n_detect80")[1], 0.06)
mata: st_local("okc", strofreal(st_numscalar("r(power)") >= 0.8 - max((0.005, 1.5 * st_numscalar("r(mcse)")))))
mata: recr("D8 confirmation rule satisfied (1 = yes)", `okc', 1, 0)
powerape n d1, power(0.8) claim(minimum) sesoi(0.02) nsim(1500) seed(222) engine(mata) nodisplay
mata: recr("D8 required n, minimum 80% (R 5451)", st_numscalar("r(n)"), refd("D8", "n_minimum80")[1], 0.06)
powerape dgp d0, focal(t, binary p(0.5)) baseline(0.30) nodisplay
powerape setape d0, target(0.10) nodisplay
powerape n d0, power(0.8) claim(detect) nsim(1500) seed(223) engine(mata) nodisplay
mata: recr("D8 required n, saturated .30->.40 (R 707; analytic 712)", st_numscalar("r(n)"), refd("D8", "n_sat80")[1], 0.06)
powerape n d0, power(0.8) claim(detect) nsim(600) seed(224) engine(stata) nodisplay
mata: recr("D8 required n, saturated, engine(stata) nsim 600 (search noise)", st_numscalar("r(n)"), refd("D8", "n_sat80")[1], 0.15)
mata: st_local("okc2", strofreal(st_numscalar("r(power)") >= 0.8 - max((0.005, 1.5 * st_numscalar("r(mcse)")))))
mata: recr("D8 engine(stata) confirmation rule satisfied (1 = yes)", `okc2', 1, 0)

* D9: MDE
powerape mde d0, n(712) claim(detect) nsim(1000) seed(225) engine(mata) nodisplay
mata: recr("D9 MDE at n 712 (R .0996)", st_numscalar("r(mde)"), refd("D9", "mde")[1], 0.06)
powerape mde d4, n(2000) claim(equivalence) nsim(1000) seed(226) engine(mata) nodisplay
mata: recr("D9 equivalence margin at n 2000 (R .0600)", st_numscalar("r(mde)"), refd("D9", "margin")[1], 0.06)

* D10: alpha rule and robustness sweep
powerape power d0, n(712) claim(detect) alpharule(0.05*sqrt(500/@n)) nsim(`NS') seed(227) engine(mata) nodisplay
mata: recd("D10 power under the alpha rule", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("D10", "power_rule"), 3)
mata: recr("D10 alpha from the rule at n 712", st_numscalar("r(alpha)"), refd("D10", "alpha_rule")[1], 1e-9)
powerape robust d1, n(700) claim(detect) vary(baseline(0.20 0.40)) nsim(2000) seed(228) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)")
forvalues i = 1/3 {
    mata: recd("D10 robust scenario `i' power (baseline sweep)", S[`i', 3], S[`i', 4], refd("D10", "robust_power_`i'"), 3)
}

* D11: justified alpha
powerape power d1, n(2000) claim(minimum) sesoi(0.02) nsim(`NS') seed(229) engine(mata) nodisplay
powerape powerat, alpha(0.01 0.025 0.05 0.10)
mata: T = st_matrix("r(table)")
mata: recd("D11 power at alpha .01", T[1, 3], T[1, 4], refd("D11", "power_at_0.01"), 3)
mata: recd("D11 power at alpha .025", T[2, 3], T[2, 4], refd("D11", "power_at_0.025"), 3)
mata: recd("D11 power at alpha .05", T[3, 3], T[3, 4], refd("D11", "power_at_0.05"), 3)
mata: recd("D11 power at alpha .10", T[4, 3], T[4, 4], refd("D11", "power_at_0.1"), 3)
powerape alpha, cost(4) prior(1) nodisplay
mata: recr("D11 alpha* (flat objective; R .0485)", st_numscalar("r(alpha)"), refd("D11", "alpha_star")[1], 0.5)
mata: recr("D11 weighted combined error at alpha*", st_numscalar("r(wcer)"), refd("D11", "wcer")[1], 0.05)

mata:
printf("\n")
for (j = 1; j <= rows(RES); j++) printf("{txt}%s\n", RES[j])
printf("\n{txt}Tier D: %f PASS, %f FAIL\n", NP, NF)
unlink("tierD-results.txt")
fh = fopen("tierD-results.txt", "w")
for (j = 1; j <= rows(RES); j++) fput(fh, RES[j])
fput(fh, sprintf("Tier D: %f PASS, %f FAIL", NP, NF))
fclose(fh)
end
log close tier
cd "`root'"
exit, clear
