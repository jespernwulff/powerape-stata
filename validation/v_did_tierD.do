* v_did_tierD.do -- difference-in-differences route, stochastic parity with
* the R package on matched designs. R reference values
* (validation/did_tierD_r.R, version 1.12.2 of the R package, 4,000
* replications per power) are compared with the Stata port at 4,000
* replications (engine(mata)) and, for a subset, fewer replications of the
* Stata engine (logit / probit / regress in the loop, and margins). RNG
* streams differ, so the criterion is
* |difference| <= 3 x sqrt(MCSE_R^2 + MCSE_Stata^2).
* Run from the repository root (the folder that holds ado/ and validation/):
*     do validation/v_did_tierD.do
* Writes validation/did-tierD-results.txt and the log validation/v_did_tierD.log.
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "v_did_tierD.do: run it from the repository root:  do validation/v_did_tierD.do"
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
log using v_did_tierD.log, text replace nomsg name(tier)
di as txt "powerape `pav'; Stata `c(stata_version)' (`c(born_date)')"
import delimited data/did-tierD-reference.csv, clear varnames(1) stringcols(1 2) asdouble
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
    ok = (s != . & r[1] != . & abs(s - r[1]) <= tol)
    RES = RES \ sprintf("%-62s Stata %8.4f (MCSE %6.4f)  R %8.4f (MCSE %6.4f)  |diff| %6.4f  tol %6.4f  %s", label, s, s_mcse, r[1], rm, abs(s - r[1]), tol, (ok ? "PASS" : "FAIL"))
    if (ok) NP++
    else NF++
}
void recr(string scalar label, real scalar s, real scalar r, real scalar reltol)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar ok
    ok = (s != . & r != . & abs(s - r) <= reltol * max((abs(r), 1e-8)))
    RES = RES \ sprintf("%-62s Stata %9.4f  R %9.4f  rel.diff %8.2e  tol %6.1e  %s", label, s, r, abs(s - r) / max((abs(r), 1e-8)), reltol, (ok ? "PASS" : "FAIL"))
    if (ok) NP++
    else NF++
}
// a share against R's, each with its binomial Monte Carlo error at nsim
void recs(string scalar label, real scalar s, real scalar r, real scalar nsim)
{
    recd(label, s, sqrt(s * (1 - s) / nsim), (r, sqrt(r * (1 - r) / 4000)), 3)
}
end
local NS 4000
local NSS 1000

* P1: logit panel, index world and analysis
powerape dgp p1, did model(logit) baseline(0.25) gap(0.10) trend(0.08) rho(0.4) nodisplay
powerape setatt p1, target(0.07) nodisplay
powerape power p1, n(1500) claim(detect) nsim(`NS') seed(401) engine(mata) nodisplay
mata: recd("P1 logit panel, detect, n 1500, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P1", "detect"), 3)
mata: recd("P1 mean of the estimates", st_numscalar("r(est_mean)"), st_numscalar("r(est_mcse)"), refd("P1", "detect_estmean"), 3)
powerape power p1, n(1500) claim(minimum) sesoi(0.02) nsim(`NS') seed(402) engine(mata) nodisplay
mata: recd("P1 minimum (sesoi .02), engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P1", "minimum"), 3)
mata: o = st_matrix("r(outcomes)")
mata: recs("P1 detect-only share", o[2], refd("P1", "detect_only")[1], `NS')
mata: recs("P1 inconclusive share", o[3], refd("P1", "inconclusive")[1], `NS')
powerape power p1, n(1500) claim(detect) nsim(`NSS') seed(403) engine(stata) nodisplay
mata: recd("P1 detect, engine(stata) logit vce(cluster), nsim 1000", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P1", "detect"), 3)
powerape power p1, n(1500) claim(detect) nsim(300) seed(404) engine(stata) effects(margins) nodisplay
mata: recd("P1 detect, engine(stata) effects(margins), nsim 300", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P1", "detect"), 3)

* P2: probit cross-sections with a covariate, model and robust SEs
powerape dgp p2, did model(probit) baseline(0.25) gap(0.10) trend(0.08) design(cross-sections) ///
    covariates((z, normal)) signal(0.2) nodisplay
powerape setatt p2, target(0.08) nodisplay
powerape power p2, n(2500) claim(detect) nsim(`NS') seed(405) engine(mata) nodisplay
mata: recd("P2 probit cross-sections + covariate, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P2", "detect"), 3)
mata: recd("P2 mean of the estimates", st_numscalar("r(est_mean)"), st_numscalar("r(est_mcse)"), refd("P2", "detect_estmean"), 3)
powerape power p2, n(2500) claim(detect) nsim(`NS') seed(406) se(robust) engine(mata) nodisplay
mata: recd("P2 robust SEs, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P2", "detect_robust"), 3)
powerape power p2, n(2500) claim(detect) nsim(`NSS') seed(407) engine(stata) nodisplay
mata: recd("P2 engine(stata) probit, nsim 1000", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P2", "detect"), 3)
powerape power p2, n(2500) claim(detect) nsim(`NSS') seed(408) se(robust) engine(stata) nodisplay
mata: recd("P2 engine(stata) probit vce(robust), nsim 1000", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P2", "detect_robust"), 3)

* P3: linear analysis in a linear world, panel
powerape dgp p3, did model(logit) baseline(0.20) gap(0.20) trend(0.10) rho(0.5) trends(linear) analysis(linear) nodisplay
powerape setatt p3, target(0.10) nodisplay
powerape power p3, n(900) claim(detect) nsim(`NS') seed(409) engine(mata) nodisplay
mata: recd("P3 linear analysis, panel, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P3", "detect"), 3)
mata: recd("P3 mean of the estimates", st_numscalar("r(est_mean)"), st_numscalar("r(est_mcse)"), refd("P3", "detect_estmean"), 3)
powerape power p3, n(900) claim(detect) nsim(`NSS') seed(410) engine(stata) nodisplay
mata: recd("P3 engine(stata) regress vce(cluster), nsim 1000", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P3", "detect"), 3)

* P4: linear analysis, cross-sections with a covariate, unequal shares
powerape dgp p4, did model(probit) baseline(0.25) gap(0.10) trend(0.08) design(cross-sections) ///
    ptreated(0.3) ppost(0.6) trends(linear) analysis(linear) covariates((z, normal)) signal(0.2) nodisplay
powerape setatt p4, target(0.08) nodisplay
powerape power p4, n(2500) claim(detect) nsim(`NS') seed(411) engine(mata) nodisplay
mata: recd("P4 linear analysis, cross-sections, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P4", "detect"), 3)
powerape power p4, n(2500) claim(detect) nsim(`NSS') seed(412) engine(stata) nodisplay
mata: recd("P4 engine(stata) regress vce(robust), nsim 1000", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P4", "detect"), 3)

* P5: analyses that assume another scale than the world has (rejection rates)
powerape dgp p5a, did model(logit) baseline(0.20) gap(0.20) trend(0.10) rho(0.4) trends(linear) analysis(index) nodisplay
powerape setatt p5a, target(0.10) nodisplay
powerape power p5a, n(1000) claim(detect) nsim(`NS') seed(413) engine(mata) nodisplay
mata: recd("P5a logit analysis in a linear world: rejection rate", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P5a", "reject"), 3)
mata: recd("P5a mean of the estimates (truth .100)", st_numscalar("r(est_mean)"), st_numscalar("r(est_mcse)"), refd("P5a", "reject_estmean"), 3)
mata: recr("P5a flagged as not valid (r(valid) = 0)", st_numscalar("r(valid)"), 0, 0)
powerape dgp p5b, did model(logit) baseline(0.20) gap(0.20) trend(0.10) rho(0.4) trends(index) analysis(linear) nodisplay
powerape setatt p5b, target(0.0667) nodisplay
powerape power p5b, n(1000) claim(detect) nsim(`NS') seed(414) engine(mata) nodisplay
mata: recd("P5b linear analysis in an index world: rejection rate", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P5b", "reject"), 3)
mata: recd("P5b mean of the estimates (truth .0667)", st_numscalar("r(est_mean)"), st_numscalar("r(est_mcse)"), refd("P5b", "reject_estmean"), 3)

* P6: several periods
powerape dgp p6, did model(logit) baseline(0.20) gap(0.20) trend(0.10) rho(0.5) npre(2) npost(2) nodisplay
powerape setatt p6, target(0.067) nodisplay
powerape power p6, n(1500) claim(detect) nsim(`NS') seed(415) engine(mata) nodisplay
mata: recd("P6 two periods a side, logit panel, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P6", "detect"), 3)
powerape power p6, n(1500) claim(detect) nsim(`NSS') seed(416) engine(stata) nodisplay
mata: recd("P6 engine(stata) logit vce(cluster), nsim 1000", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P6", "detect"), 3)
powerape dgp p6x, did model(probit) baseline(0.20) gap(0.20) trend(0.10) design(cross-sections) npre(2) npost(3) ///
    trends(linear) analysis(linear) nodisplay
powerape setatt p6x, target(0.08) nodisplay
powerape power p6x, n(2000) claim(detect) nsim(`NS') seed(417) engine(mata) nodisplay
mata: recd("P6x 2 + 3 periods, cross-sections, linear analysis", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P6x", "detect"), 3)

* P7: a binary subgroup (difference between the two subgroups' effects)
powerape dgp p7, did model(logit) baseline(0.30 0.20) gap(0.10 0.15) trend(0.10 0.05) rho(0.5) subgroup(female, p(0.4)) nodisplay
powerape setatt p7, target(0.10 0.03) nodisplay
powerape power p7, n(4000) claim(detect) nsim(`NS') seed(418) engine(mata) nodisplay
mata: recd("P7 subgroup difference, logit, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P7", "detect"), 3)
mata: recd("P7 mean of the estimates (truth -.070)", st_numscalar("r(est_mean)"), st_numscalar("r(est_mcse)"), refd("P7", "detect_estmean"), 3)
powerape power p7, n(4000) claim(detect) nsim(600) seed(419) engine(stata) nodisplay
mata: recd("P7 engine(stata) logit vce(cluster), nsim 600", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P7", "detect"), 3)
powerape power p7, n(4000) claim(detect) nsim(200) seed(420) engine(stata) effects(margins) nodisplay
mata: recd("P7 engine(stata) margins over(), nsim 200", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P7", "detect"), 3)
powerape dgp p7l, did model(logit) baseline(0.30 0.20) gap(0.10 0.15) trend(0.10 0.05) rho(0.5) subgroup(female, p(0.4)) ///
    trends(linear) analysis(linear) nodisplay
powerape setatt p7l, target(0.10 0.03) nodisplay
powerape power p7l, n(4000) claim(detect) nsim(`NS') seed(421) engine(mata) nodisplay
mata: recd("P7l subgroup difference, linear analysis, engine(mata)", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P7l", "detect"), 3)
powerape power p7l, n(4000) claim(detect) nsim(`NSS') seed(422) engine(stata) nodisplay
mata: recd("P7l engine(stata) regress vce(cluster), nsim 1000", st_numscalar("r(power)"), st_numscalar("r(mcse)"), refd("P7l", "detect"), 3)

* P8: required n (answers differ by search noise: within 6% of R's, and each
*     must satisfy the confirmation rule)
powerape n p1, power(0.8) claim(detect) nsim(1500) seed(423) engine(mata) nodisplay
mata: recr("P8 required n, P1 detect 80%", st_numscalar("r(n)"), refd("P8", "n_p1_detect80")[1], 0.06)
mata: st_local("okc", strofreal(st_numscalar("r(power)") >= 0.8 - max((0.005, 1.5 * st_numscalar("r(mcse)")))))
mata: recr("P8 confirmation rule satisfied (1 = yes)", `okc', 1, 0)
powerape n p3, power(0.8) claim(detect) nsim(1500) seed(424) engine(mata) nodisplay
mata: recr("P8 required n, P3 (linear analysis) detect 80%", st_numscalar("r(n)"), refd("P8", "n_p3_detect80")[1], 0.06)
powerape n p6, power(0.8) claim(detect) nsim(1500) seed(425) engine(mata) nodisplay
mata: recr("P8 required n, P6 (two periods a side) detect 80%", st_numscalar("r(n)"), refd("P8", "n_p6_detect80")[1], 0.06)

* P9: minimum detectable effect
powerape mde p1, n(1500) claim(detect) nsim(1000) seed(426) engine(mata) nodisplay
mata: recr("P9 MDE, P1 at 1,500 units", st_numscalar("r(mde)"), refd("P9", "mde_p1")[1], 0.06)
powerape mde p3, n(900) claim(detect) nsim(1000) seed(427) engine(mata) nodisplay
mata: recr("P9 MDE, P3 (linear analysis) at 900 units", st_numscalar("r(mde)"), refd("P9", "mde_p3")[1], 0.06)
powerape mde p7, n(4000) claim(detect) nsim(1000) seed(428) engine(mata) nodisplay
mata: recr("P9 MDE of the subgroup difference, P7 at 4,000 units", st_numscalar("r(mde)"), refd("P9", "mde_p7")[1], 0.06)

* P10: robustness sweep (rho; the scale on which trends are parallel)
powerape robust p1, n(1500) claim(detect) vary(rho(0.2 0.6)) nsim(2000) seed(429) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)")
forvalues i = 1/3 {
    mata: recd("P10 sweep over rho, scenario `i' power", S[`i', 3], S[`i', 4], refd("P10", "rho_power_`i'"), 3)
}
powerape dgp pg, did model(logit) baseline(0.25) gap(0.20) trend(0.08) rho(0.4) nodisplay
powerape setatt pg, target(0.07) nodisplay
powerape robust pg, n(1500) claim(detect) vary(trends(index linear)) nsim(2000) seed(430) nonmax engine(mata) nodisplay
mata: S = st_matrix("r(scenarios)")
forvalues i = 1/2 {
    mata: recd("P10 sweep over trends, scenario `i' power / rejection rate", S[`i', 3], S[`i', 4], refd("P10", "trends_power_`i'"), 3)
    mata: recd("P10 sweep over trends, scenario `i' mean estimate", S[`i', 7], 0.0015, (refd("P10", "trends_estmean_`i'")[1], 0.0015), 3)
    mata: recr("P10 sweep over trends, scenario `i' valid flag", S[`i', 8], refd("P10", "trends_valid_`i'")[1], 0)
}

* P11: justified alpha
powerape power p1, n(1500) claim(minimum) sesoi(0.02) nsim(`NS') seed(431) engine(mata) nodisplay
powerape powerat, alpha(0.01 0.025 0.05 0.10)
mata: T = st_matrix("r(table)")
mata: recd("P11 power at alpha .01", T[1, 3], T[1, 4], refd("P11", "power_at_0.01"), 3)
mata: recd("P11 power at alpha .025", T[2, 3], T[2, 4], refd("P11", "power_at_0.025"), 3)
mata: recd("P11 power at alpha .05", T[3, 3], T[3, 4], refd("P11", "power_at_0.05"), 3)
mata: recd("P11 power at alpha .10", T[4, 3], T[4, 4], refd("P11", "power_at_0.1"), 3)
powerape alpha, cost(4) prior(1) nodisplay
mata: recr("P11 alpha* (flat objective)", st_numscalar("r(alpha)"), refd("P11", "alpha_star")[1], 0.5)
mata: recr("P11 weighted combined error at alpha*", st_numscalar("r(wcer)"), refd("P11", "wcer")[1], 0.05)

mata:
printf("\n")
for (j = 1; j <= rows(RES); j++) printf("{txt}%s\n", RES[j])
printf("\n{txt}Tier D (difference-in-differences): %f PASS, %f FAIL\n", NP, NF)
unlink("did-tierD-results.txt")
fh = fopen("did-tierD-results.txt", "w")
for (j = 1; j <= rows(RES); j++) fput(fh, RES[j])
fput(fh, sprintf("Tier D (difference-in-differences): %f PASS, %f FAIL", NP, NF))
fclose(fh)
end
log close tier
cd "`root'"
exit, clear
