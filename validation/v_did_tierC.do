* v_did_tierC.do -- difference-in-differences route, calibration constants
* against the R package: the index coefficients (a, b, l, th, de), the block
* rates, the two readings of the rates, the true effect, the coefficient
* vector in design-matrix order, and the pilot-data estimates. Designs
* without covariates are closed forms and must agree to root-finding
* tolerance; designs calibrated on a Monte Carlo integration draw differ by
* integration noise (the two packages draw their own sets). Reference values
* come from validation/did_export_fixtures.R, run on the development version
* of the R package.
* Run from the repository root (the folder that holds ado/ and validation/):
*     do validation/v_did_tierC.do
* Writes validation/did-tierC-results.txt and the log validation/v_did_tierC.log.
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "v_did_tierC.do: run it from the repository root:  do validation/v_did_tierC.do"
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
log using v_did_tierC.log, text replace nomsg name(tier)
di as txt "powerape `pav'; Stata `c(stata_version)' (`c(born_date)')"
import delimited data/did-reference.csv, clear varnames(1) stringcols(1 2) asdouble
mata:
REFC = st_sdata(., "case"); REFI = st_sdata(., "item"); REFV = st_data(., "value")
real scalar refv(string scalar c, string scalar it)
{
    external string colvector REFC, REFI
    external real colvector REFV
    real scalar j
    for (j = 1; j <= rows(REFC); j++) if (REFC[j] == c & REFI[j] == it) return(REFV[j])
    return(.)
}
RES = J(0, 1, ""); NP = 0; NF = 0
// absolute difference against tol x max(1, |reference|)
void recc(string scalar cs, string scalar item, real scalar stata, real scalar r, real scalar tol, string scalar note)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar d, ok
    d = abs(stata - r)
    ok = (stata != . & r != . & d <= tol * max((1, abs(r))))
    RES = RES \ sprintf("%-6s %-15s Stata %13.8f  R %13.8f  |diff| %9.2e  tol %7.1e  %s  %s", cs, item, stata, r, d, tol, (ok ? "PASS" : "FAIL"), note)
    if (ok) NP++
    else NF++
}
end

* compare everything a design stores; M = number of subgroups (1 or 2);
* tolc for coefficients, tolr for rates, readings and the true effect
program define calcheck
    args cs M tolc tolr note
    powerape describe `cs'
    tempname R C
    matrix `R' = r(rates)
    matrix `C' = r(didcoef)
    local rl = r(reading_linear)
    local ri = r(reading_index)
    local sig = r(signal)
    forvalues m = 1/`M' {
        local sfx = cond(`M' == 1, "", "_" + string(`m' - 1))
        local j 0
        foreach nm in a b l th de {
            local ++j
            mata: recc("`cs'", "`nm'`sfx'", st_matrix("`C'")[`m', `j'], refv("`cs'", "`nm'`sfx'"), `tolc', "`note'")
        }
        local j 0
        foreach nm in c00 c01 c10 cf c11 {
            local ++j
            mata: recc("`cs'", "`nm'`sfx'", st_matrix("`R'")[`m', `j'], refv("`cs'", "`nm'`sfx'"), `tolr', "`note'")
        }
    }
    mata: recc("`cs'", "reading_linear", `rl', refv("`cs'", "reading_linear"), `tolr', "`note'")
    mata: recc("`cs'", "reading_index", `ri', refv("`cs'", "reading_index"), `tolr', "`note'")
    mata: recc("`cs'", "signal_implied", `sig', refv("`cs'", "signal_implied"), `tolr', "`note'")
    powerape trueatt `cs'
    mata: recc("`cs'", "true effect", st_numscalar("r(att)"), refv("`cs'", "true"), `tolr', "`note'")
    mata: bt = pape_did_beta(pape_get_dgp("`cs'"))
    mata: for (j = 1; j <= rows(bt); j++) recc("`cs'", sprintf("beta%f", j), bt[j], refv("`cs'", sprintf("beta%f", j)), `tolc', "`note'")
    mata: recc("`cs'", "length of beta", rows(bt), sum(REFC :== "`cs'" :& substr(REFI, 1, 4) :== "beta"), 0, "`note'")
end

* closed forms (no covariates)
powerape dgp c1, did model(logit) baseline(0.20) gap(0.20) trend(0.10) rho(0.4) nodisplay
powerape setatt c1, target(0.07) nodisplay
calcheck c1 1 1e-8 1e-9 "closed form"
powerape dgp c2, did model(probit) baseline(0.20) gap(0.20) trend(0.10) trends(linear) analysis(linear) nodisplay
powerape setatt c2, target(0.10) nodisplay
calcheck c2 1 1e-8 1e-9 "closed form"
powerape dgp c3, did model(probit) baseline(0.25) gap(0.10) trend(0.08) design(cross-sections) npre(2) npost(3) nodisplay
powerape setatt c3, target(-0.05) nodisplay
calcheck c3 1 1e-8 1e-9 "closed form"
powerape dgp c4, did model(logit) baseline(0.30 0.20) gap(0.10 0.15) trend(0.10 0.05) rho(0.5) subgroup(female, p(0.4)) nodisplay
powerape setatt c4, target(0.10 0.03) nodisplay
calcheck c4 2 1e-8 1e-9 "closed form"
powerape dgp c5, did model(probit) baseline(0.30 0.05) gap(0.10 0.05) trend(0.10 0.10) trends(linear) subgroup(female, p(0.5)) nodisplay
powerape setatt c5, target(0.10 0.03) nodisplay
calcheck c5 2 1e-8 1e-9 "closed form"

* integration draws (nint 4e5 by default; each package draws its own set)
powerape dgp c6, did model(logit) baseline(0.25) gap(0.10) trend(0.08) rho(0.4) ///
    covariates((z, normal) (w, binary p(0.4))) correlation(0.3) signal(0.2) nodisplay
powerape setatt c6, target(0.07) nodisplay
calcheck c6 1 1e-2 2e-3 "integration draw, nint 4e5"
powerape describe c6
mata: recc("c6", "gamma1", st_matrix("r(gamma)")[1], refv("c6", "gamma1"), 5e-3, "integration draw, nint 4e5")
mata: recc("c6", "gamma2", st_matrix("r(gamma)")[2], refv("c6", "gamma2"), 5e-3, "integration draw, nint 4e5")
powerape dgp c7, did model(probit) baseline(0.25) gap(0.10) trend(0.08) trends(linear) ///
    covariates((z, normal) (w, binary p(0.4))) correlation(0.3) signal(0.2) nodisplay
powerape setatt c7, target(0.07) nodisplay
calcheck c7 1 1e-2 2e-3 "integration draw, nint 4e5"
powerape dgp c6big, did model(logit) baseline(0.25) gap(0.10) trend(0.08) rho(0.4) ///
    covariates((z, normal) (w, binary p(0.4))) correlation(0.3) signal(0.2) nint(4e6) nodisplay
powerape setatt c6big, target(0.07) nodisplay
calcheck c6big 1 3e-3 6e-4 "integration draw, nint 4e6"
powerape dgp c8, did model(logit) baseline(0.30 0.20) gap(0.10 0.15) trend(0.10 0.05) rho(0.5) ///
    subgroup(female, p(0.4)) covariates((z, normal)) signal(0.15) nodisplay
powerape setatt c8, target(0.10 0.03) nodisplay
calcheck c8 2 1e-2 2e-3 "integration draw, nint 4e5"

* the copula's joint rate, the pilot's method-of-moments solver, and the
* estimates from a pilot panel that both packages read
mata: recc("pilot", "joint rate", pape_did_joint_rate(0.30, 0.30, 0.5), refv("pilot", "joint_30_05"), 1e-9, "bivariate normal")
mata: recc("pilot", "joint rate", pape_did_joint_rate(0.08, 0.08, 0.2), refv("pilot", "joint_08_02"), 1e-9, "bivariate normal")
mata: recc("pilot", "rho solver", pape_did_pilot_rho((0.25 \ 0.35), (0.11 \ 0.19), (800 \ 760)), refv("pilot", "rho_solve"), 1e-7, "method of moments")
import delimited data/didfix_pilot.csv, clear varnames(1) asdouble
powerape dgp dp, didpilot outcome(y) group(g) id(id) period(period) trend(0.05) model(logit) nodisplay
matrix PL = r(pilot)
matrix RT = r(rates)
mata: PL = st_matrix("PL")
mata: recc("pilot", "baseline", PL[5], refv("pilot", "baseline"), 1e-12, "same pilot file")
mata: recc("pilot", "gap", PL[6], refv("pilot", "gap"), 1e-12, "same pilot file")
mata: recc("pilot", "p_treated", PL[7], refv("pilot", "p_treated"), 1e-12, "same pilot file")
mata: recc("pilot", "rho", PL[8], refv("pilot", "rho"), 1e-7, "same pilot file")
mata: recc("pilot", "control_change", PL[9], refv("pilot", "control_change"), 1e-12, "same pilot file")
mata: recc("pilot", "n_pairs", PL[4], refv("pilot", "n_pairs"), 0, "same pilot file")
mata: recc("pilot", "cf", st_matrix("RT")[1, 4], refv("pilot", "cf"), 1e-9, "same pilot file")
clear

mata:
printf("\n")
for (j = 1; j <= rows(RES); j++) printf("{txt}%s\n", RES[j])
printf("\n{txt}Tier C (difference-in-differences): %f PASS, %f FAIL\n", NP, NF)
unlink("did-tierC-results.txt")
fh = fopen("did-tierC-results.txt", "w")
for (j = 1; j <= rows(RES); j++) fput(fh, RES[j])
fput(fh, sprintf("Tier C (difference-in-differences): %f PASS, %f FAIL", NP, NF))
fclose(fh)
end
log close tier
cd "`root'"
exit, clear
