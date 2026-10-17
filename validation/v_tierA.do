* v_tierA.do -- exact finite-sample power by triple binomial enumeration
* (powerape's V1 / helper-exact.R) against the Stata port's simulated power.
* In the saturated design (binary focal, no covariates) the whole decision
* rule is a function of four binomial counts, so its power has an exact value
* that neither implementation can tune to. Both engines are checked.
* Run from the repository root (the folder that holds ado/ and validation/):
*     do validation/v_tierA.do
* Writes validation/tierA-results.txt and the log validation/v_tierA.log.
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "v_tierA.do: run it from the repository root:  do validation/v_tierA.do"
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
log using v_tierA.log, text replace nomsg name(tier)
di as txt "powerape `pav'; Stata `c(stata_version)' (`c(born_date)')"

mata:
// port of exact_power_sat (helper-exact.R)
real scalar exact_power_sat(real scalar n, real scalar pd, real scalar p0, real scalar p1,
    real scalar conf, string scalar claim, real scalar sesoi, string scalar separation)
{
    real scalar z, total, n1, n0, pn1, c1, c0
    real colvector x1, x0, p1h, p0h
    real matrix est, se, dec, w, bad
    z = invnormal(1 - (1 - conf) / 2)
    total = 0
    for (n1 = 0; n1 <= n; n1++) {
        pn1 = binomialp(n, n1, pd)
        if (pn1 < 1e-14) continue
        n0 = n - n1
        if (n1 == 0 | n0 == 0) continue
        if (separation == "fail") {
            if (n1 < 2 | n0 < 2) continue
            x1 = (1::(n1 - 1)); x0 = (1::(n0 - 1))
        }
        else {
            x1 = (0::n1); x0 = (0::n0)
        }
        p1h = x1 :/ n1; p0h = x0 :/ n0
        c1 = rows(p1h); c0 = rows(p0h)
        est = p1h * J(1, c0, 1) - J(c1, 1, 1) * p0h'
        se = sqrt((p1h :* (1 :- p1h) :/ n1) * J(1, c0, 1) + J(c1, 1, 1) * (p0h :* (1 :- p0h) :/ n0)')
        if (claim == "detect") dec = (est :- z :* se :> 0)
        else if (claim == "minimum") dec = (est :- z :* se :> sesoi)
        else dec = (est :+ z :* se :< sesoi) :& (est :- z :* se :> -sesoi)
        if (separation == "keep") {
            bad = ((x1 :== 0) * (x0 :== 0)') :| ((x1 :== n1) * (x0 :== n0)')
            dec = dec :& !bad
        }
        w = binomialp(n1, x1, p1) * binomialp(n0, x0, p0)'
        total = total + pn1 * sum(w :* dec)
    }
    return(total)
}
RES = J(0, 1, "")
NP = 0; NF = 0
void reca(string scalar label, real scalar sim, real scalar ex, real scalar nsim)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar mcse, ok
    mcse = sqrt(max((ex * (1 - ex), 1e-6)) / nsim)
    ok = (abs(sim - ex) <= 4 * mcse + 1e-4)
    RES = RES \ sprintf("%-52s exact %7.4f  simulated %7.4f  MCSE %6.4f  |diff| %6.4f  %s", label, ex, sim, mcse, abs(sim - ex), (ok ? "PASS" : "FAIL"))
    if (ok) NP++
    else NF++
}
end

* the V1 design: probit, p(focal) .5, baseline .30, APE .10 (p1 .40), n = 712
powerape dgp d10, focal(t, binary p(0.5)) baseline(0.30) nodisplay
powerape setape d10, target(0.10) nodisplay
powerape dgp d00, focal(t, binary p(0.5)) baseline(0.30) nodisplay
powerape setape d00, target(0) nodisplay

* Mata engine at nsim = 40,000 (as the R battery)
mata: ex = exact_power_sat(712, .5, .30, .40, .95, "detect", ., "fail")
powerape power d10, n(712) claim(detect) nsim(40000) seed(1) engine(mata) nodisplay nodraws
mata: reca("V1 detect, n = 712, engine(mata), nsim 40000", st_numscalar("r(power)"), ex, 40000)
mata: ex = exact_power_sat(712, .5, .30, .40, .90, "minimum", .04, "fail")
powerape power d10, n(712) claim(minimum) sesoi(0.04) nsim(40000) seed(2) engine(mata) nodisplay nodraws
mata: reca("V1 minimum (sesoi .04), engine(mata), nsim 40000", st_numscalar("r(power)"), ex, 40000)
mata: ex = exact_power_sat(712, .5, .30, .30, .90, "equivalence", .05, "fail")
powerape power d00, n(712) claim(equivalence) sesoi(0.05) nsim(40000) seed(3) engine(mata) nodisplay nodraws
mata: reca("V1 equivalence (sesoi .05, truth 0), engine(mata)", st_numscalar("r(power)"), ex, 40000)
mata: ex = exact_power_sat(712, .5, .30, .30, .90, "equivalence", .07, "fail")
powerape power d00, n(712) claim(equivalence) sesoi(0.07) nsim(40000) seed(4) engine(mata) nodisplay nodraws
mata: reca("equivalence (sesoi .07, truth 0), engine(mata)", st_numscalar("r(power)"), ex, 40000)

* Stata engine (probit in the loop) at nsim = 4,000
mata: ex = exact_power_sat(712, .5, .30, .40, .95, "detect", ., "fail")
powerape power d10, n(712) claim(detect) nsim(4000) seed(5) engine(stata) nodisplay nodraws
mata: reca("V1 detect, n = 712, engine(stata), nsim 4000", st_numscalar("r(power)"), ex, 4000)
mata: ex = exact_power_sat(712, .5, .30, .40, .90, "minimum", .04, "fail")
powerape power d10, n(712) claim(minimum) sesoi(0.04) nsim(4000) seed(6) engine(stata) nodisplay nodraws
mata: reca("V1 minimum (sesoi .04), engine(stata), nsim 4000", st_numscalar("r(power)"), ex, 4000)
* margins in the loop at nsim = 1,000
mata: ex = exact_power_sat(712, .5, .30, .40, .95, "detect", ., "fail")
powerape power d10, n(712) claim(detect) nsim(1000) seed(7) engine(stata) effects(margins) nodisplay nodraws
mata: reca("V1 detect, engine(stata) effects(margins), nsim 1000", st_numscalar("r(power)"), ex, 1000)

* rare outcome, protective effect: .05 -> .01 at n = 300 (separation matters)
powerape dgp drare, focal(t, binary p(0.5)) baseline(0.05) nodisplay
powerape setape drare, target(-0.04) nodisplay
mata: ex = exact_power_sat(300, .5, .01, .05, .95, "detect", ., "fail")
powerape power drare, n(300) claim(detect) nsim(20000) seed(8) engine(mata) nodisplay nodraws separation(fail)
mata: reca("rare .05->.01, n = 300, separation(fail), engine(mata)", st_numscalar("r(power)"), ex, 20000)
mata: st_local("sepd", strofreal(st_numscalar("r(separated)")))
di as txt "share of simulated studies with separation: `sepd' (R: > .15)"
mata: ex = exact_power_sat(300, .5, .01, .05, .95, "detect", ., "keep")
powerape power drare, n(300) claim(detect) nsim(20000) seed(9) engine(mata) nodisplay nodraws separation(keep)
mata: reca("rare .05->.01, n = 300, separation(keep), engine(mata)", st_numscalar("r(power)"), ex, 20000)
mata: ex = exact_power_sat(300, .5, .01, .05, .95, "detect", ., "fail")
powerape power drare, n(300) claim(detect) nsim(4000) seed(10) engine(stata) nodisplay nodraws separation(fail)
mata: reca("rare .05->.01, n = 300, separation(fail), engine(stata)", st_numscalar("r(power)"), ex, 4000)
* 0.1.1: the Stata engine keeps the asis fit's separated coefficient (zero variance
* row), so separation(keep) now scores the degenerate interval as R's glm does
mata: ex = exact_power_sat(300, .5, .01, .05, .95, "detect", ., "keep")
powerape power drare, n(300) claim(detect) nsim(4000) seed(10) engine(stata) nodisplay nodraws separation(keep)
mata: reca("rare .05->.01, n = 300, separation(keep), engine(stata)", st_numscalar("r(power)"), ex, 4000)
mata: st_local("nfail", strofreal(st_numscalar("r(n_failed)")))
di as txt "separation(keep), engine(stata): replications failed to fit: `nfail' (asis fits not converged or omitted terms)"

* logit gives the identical saturated-case decisions on the same seed
powerape dgp dl, model(logit) focal(t, binary p(0.5)) baseline(0.30) nodisplay
powerape setape dl, target(0.10) nodisplay
powerape power dl, n(500) claim(detect) nsim(2000) seed(11) engine(mata) nodisplay nodraws
local pl = r(power)
powerape power d10, n(500) claim(detect) nsim(2000) seed(11) engine(mata) nodisplay nodraws
mata: reca("probit = logit, saturated case, same seed", `pl', st_numscalar("r(power)"), 1e12)

mata:
printf("\n")
for (j = 1; j <= rows(RES); j++) printf("{txt}%s\n", RES[j])
printf("\n{txt}Tier A: %f PASS, %f FAIL\n", NP, NF)
unlink("tierA-results.txt")
fh = fopen("tierA-results.txt", "w")
for (j = 1; j <= rows(RES); j++) fput(fh, RES[j])
fput(fh, sprintf("Tier A: %f PASS, %f FAIL", NP, NF))
fclose(fh)
end
log close tier
cd "`root'"
exit, clear
