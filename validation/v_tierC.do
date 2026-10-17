* v_tierC.do -- calibration and inversion constants against the R package.
* Closed-form worlds (no covariates) must agree to root-finding tolerance;
* worlds calibrated on a Monte Carlo integration draw differ by integration
* noise (the two packages draw their own 400,000-row sets), which shrinks
* with the draw size (checked at n_int = 4e6 too).
* Run from the repository root (the folder that holds ado/ and validation/):
*     do validation/v_tierC.do
* Writes validation/tierC-results.txt and the log validation/v_tierC.log.
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "v_tierC.do: run it from the repository root:  do validation/v_tierC.do"
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
log using v_tierC.log, text replace nomsg name(tier)
di as txt "powerape `pav'; Stata `c(stata_version)' (`c(born_date)')"
import delimited data/reference.csv, clear varnames(1) stringcols(1 2) asdouble
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
void recc(string scalar cs, string scalar item, real scalar stata, real scalar r, real scalar reltol, string scalar note)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar d, rel, ok
    d = abs(stata - r)
    rel = d / max((abs(r), 1e-8))
    ok = (rel <= reltol)
    RES = RES \ sprintf("%-6s %-11s Stata %13.7f  R %13.7f  rel.diff %9.2e  tol %7.1e  %s  %s", cs, item, stata, r, rel, reltol, (ok ? "PASS" : "FAIL"), note)
    if (ok) NP++
    else NF++
}
end
* cal0: closed form (no covariates)
powerape dgp c0, focal(t, binary p(0.5)) baseline(0.30) nodisplay
powerape setape c0, target(0.10) nodisplay
mata: recc("cal0", "beta0", st_numscalar("r(beta0)"), refv("cal0", "beta0"), 1e-9, "closed form")
mata: recc("cal0", "beta_focal", st_numscalar("r(beta_focal)"), refv("cal0", "beta_focal"), 1e-8, "closed form")
mata: recc("cal0", "p1", st_numscalar("r(p1)"), refv("cal0", "p1"), 1e-8, "closed form")
* cal1: one normal covariate, signal .15 (integration draw, n_int 4e5)
powerape dgp c1, focal(treat, binary p(0.5)) covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.15) nodisplay
powerape setape c1, target(0.05) nodisplay
mata: recc("cal1", "beta0", st_numscalar("r(beta0)"), refv("cal1", "beta0"), 1e-2, "n_int 4e5, own draws")
mata: recc("cal1", "gamma1", st_matrix("r(gamma)")[1], refv("cal1", "gamma1"), 1e-2, "n_int 4e5, own draws")
mata: recc("cal1", "beta_focal", st_numscalar("r(beta_focal)"), refv("cal1", "beta_focal"), 1e-2, "n_int 4e5, own draws")
mata: recc("cal1", "p1", st_numscalar("r(p1)"), refv("cal1", "p1"), 1e-3, "implied cell rate")
powerape dgp c1h, focal(treat, binary p(0.5)) covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.15) nint(4000000) nodisplay
powerape setape c1h, target(0.05) nodisplay
mata: recc("cal1h", "beta0", st_numscalar("r(beta0)"), refv("cal1h", "beta0"), 3e-3, "n_int 4e6")
mata: recc("cal1h", "gamma1", st_matrix("r(gamma)")[1], refv("cal1h", "gamma1"), 3e-3, "n_int 4e6")
mata: recc("cal1h", "beta_focal", st_numscalar("r(beta_focal)"), refv("cal1h", "beta_focal"), 3e-3, "n_int 4e6")
* cal2: AIE, no covariates (exact)
powerape dgp c2, focal(treat, binary p(0.5)) moderator(female, binary p(0.55)) baseline(0.30) nodisplay
powerape setaie c2, target(0.04) mainfocal(0.10) mainmoderator(0.05) nodisplay
mata: recc("cal2", "beta_focal", st_numscalar("r(beta_focal)"), refv("cal2", "beta_focal"), 1e-7, "closed form")
mata: recc("cal2", "beta_mod", st_numscalar("r(beta_mod)"), refv("cal2", "beta_mod"), 1e-7, "closed form")
mata: recc("cal2", "beta_int", st_numscalar("r(beta_int)"), refv("cal2", "beta_int"), 1e-6, "closed form (grid+Brent)")
* cal3: logit, continuous focal, two correlated covariates
powerape dgp c3, model(logit) focal(dose, normal mean(2) sd(1)) covariates((z1, normal) (z2, binary p(0.4))) ///
    correlation(0.2) baseline(0.25) signal(0.10) nodisplay
powerape setape c3, target(0.03) nodisplay
mata: recc("cal3", "beta0", st_numscalar("r(beta0)"), refv("cal3", "beta0"), 1e-2, "n_int 4e5, own draws")
mata: recc("cal3", "gamma1", st_matrix("r(gamma)")[1], refv("cal3", "gamma1"), 1e-2, "n_int 4e5, own draws")
mata: recc("cal3", "gamma2", st_matrix("r(gamma)")[2], refv("cal3", "gamma2"), 1e-2, "n_int 4e5, own draws")
mata: recc("cal3", "beta_focal", st_numscalar("r(beta_focal)"), refv("cal3", "beta_focal"), 1e-2, "n_int 4e5, own draws")
powerape dgp c3h, model(logit) focal(dose, normal mean(2) sd(1)) covariates((z1, normal) (z2, binary p(0.4))) ///
    correlation(0.2) baseline(0.25) signal(0.10) nint(4000000) nodisplay
powerape setape c3h, target(0.03) nodisplay
mata: recc("cal3h", "beta0", st_numscalar("r(beta0)"), refv("cal3h", "beta0"), 3e-3, "n_int 4e6")
mata: recc("cal3h", "beta_focal", st_numscalar("r(beta_focal)"), refv("cal3h", "beta_focal"), 3e-3, "n_int 4e6")
* calp: panel (n_int 1e5 units)
powerape dgp cp, panel focal(treat, binary p(0.5) icc(1)) covariates((size, normal icc(0.6))) ///
    nperiods(4) rho(0.3) creshare(0.5) baseline(0.30) signal(0.10) nodisplay
powerape setape cp, target(0.08) nodisplay
mata: recc("calp", "beta0", st_numscalar("r(beta0)"), refv("calp", "beta0"), 2e-2, "n_int 1e5 units")
mata: recc("calp", "gamma1", st_matrix("r(gamma)")[1], refv("calp", "gamma1"), 2e-2, "n_int 1e5 units")
mata: recc("calp", "beta_focal", st_numscalar("r(beta_focal)"), refv("calp", "beta_focal"), 2e-2, "n_int 1e5 units")
mata: d = pape_get_dgp("cp"); recc("calp", "xi1", d.xi[1], refv("calp", "xi1"), 2e-2, "Mundlak loading"); recc("calp", "var_a", d.var_a, refv("calp", "var_a"), 1e-9, "closed form")
powerape trueape cp
mata: recc("calp", "true_ape", st_numscalar("r(ape)"), refv("calp", "true_ape"), 1e-6, "pinned target")
* cali: IV (n_int 1e5)
powerape dgp ci, iv focal(training, binary p(0.4)) covariates((size, normal)) instruments((subsidy, binary p(0.5))) ///
    endogeneity(0.4) ivstrength(0.2) baseline(0.30) signal(0.10) nodisplay
powerape setape ci, target(0.10) nodisplay
mata: d = pape_get_dgp("ci")
mata: recc("cali", "beta0", d.beta0, refv("cali", "beta0"), 2e-2, "n_int 1e5")
mata: recc("cali", "gamma1", d.gamma[1], refv("cali", "gamma1"), 2e-2, "n_int 1e5")
mata: recc("cali", "pi_z1", d.pi_z[1], refv("cali", "pi_z1"), 2e-2, "n_int 1e5")
mata: recc("cali", "kappa0", d.kappa0, refv("cali", "kappa0"), 2e-2, "n_int 1e5")
mata: recc("cali", "beta_focal", d.beta_focal, refv("cali", "beta_focal"), 2e-2, "n_int 1e5")
mata:
printf("\n")
for (j = 1; j <= rows(RES); j++) printf("{txt}%s\n", RES[j])
printf("\n{txt}Tier C: %f PASS, %f FAIL\n", NP, NF)
unlink("tierC-results.txt")
fh = fopen("tierC-results.txt", "w")
for (j = 1; j <= rows(RES); j++) fput(fh, RES[j])
fput(fh, sprintf("Tier C: %f PASS, %f FAIL", NP, NF))
fclose(fh)
end
log close tier
cd "`root'"
exit, clear
