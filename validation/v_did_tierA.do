* v_did_tierA.do -- difference-in-differences route, exact power by
* enumeration. Repeated cross-sections without covariates are saturated in
* the four group-by-period cells: given the cell sizes every study is one of
* prod(n_c + 1) outcomes with a known binomial probability, and the estimate
* and its standard error are closed forms in the four cell means (index
* analysis: p11 - G(G^-1(p10) + G^-1(p01) - G^-1(p00)) with its delta-method
* SE; linear analysis: the cross difference with the HC1 SE). The exact
* conditional power is averaged over 1,000 draws of the multinomial cell
* sizes (the draws the R reference used, read from a file), so the value has
* a small standard error of its own. Neither engine can tune to it.
* Checks: (1) the Mata oracle reproduces the R oracle on the same cell sizes;
* (2) both engines reproduce the exact power and the share of studies
* without a usable estimate. n = 120, detection claim, control rates .30 and
* .40, treated .40 before and an effect of .30.
* Run from the repository root (the folder that holds ado/ and validation/):
*     do validation/v_did_tierA.do
* Writes validation/did-tierA-results.txt and the log validation/v_did_tierA.log.
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "v_did_tierA.do: run it from the repository root:  do validation/v_did_tierA.do"
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
log using v_did_tierA.log, text replace nomsg name(tier)
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
end
import delimited data/did-tierA-sizes.csv, clear varnames(1) asdouble
mata:
SIZES = st_data(., ("n00", "n01", "n10", "n11"))
// exact power of the detection claim and exact share of studies without a
// usable estimate; pc and the columns of sizes are in the order c00, c01,
// c10, c11. Returns (power, its SE over the size draws, failed, its SE).
real rowvector exact_did(real rowvector pc, real matrix sizes, string scalar link,
    string scalar analysis, real scalar conf)
{
    real scalar z, n, draws, r, T, n4, n1, n2, n3
    real rowvector nc
    real colvector outp, fail, p1, p2, p3, p4, w1, w2, w3, w4, l1, l2, l3, v1, v2, v3
    real colvector A1, A2, A3, W, L, GL, AA, base, one, edge
    real matrix D
    z = invnormal(1 - (1 - conf) / 2)
    draws = rows(sizes)
    n = sum(sizes[1, .])
    outp = J(draws, 1, 0); fail = J(draws, 1, 0)
    for (r = 1; r <= draws; r++) {
        nc = sizes[r, .]
        if (any(nc :== 0)) {
            fail[r] = 1
            continue
        }
        p1 = (0::nc[1]) :/ nc[1]; w1 = binomialp(nc[1], 0::nc[1], pc[1])
        p2 = (0::nc[2]) :/ nc[2]; w2 = binomialp(nc[2], 0::nc[2], pc[2])
        p3 = (0::nc[3]) :/ nc[3]; w3 = binomialp(nc[3], 0::nc[3], pc[3])
        p4 = (0::nc[4]) :/ nc[4]; w4 = binomialp(nc[4], 0::nc[4], pc[4])
        if (analysis == "linear") {
            // every outcome; no usable estimate when the standard error is
            // zero (all four cell means at 0 or 1)
            fail[r] = (w1[1] + w1[rows(w1)]) * (w2[1] + w2[rows(w2)]) * (w3[1] + w3[rows(w3)]) * (w4[1] + w4[rows(w4)])
            n1 = rows(p1); n2 = rows(p2); n3 = rows(p3); n4 = rows(p4)
            A1 = p1 # J(n2 * n3, 1, 1)
            A2 = J(n1, 1, 1) # p2 # J(n3, 1, 1)
            A3 = J(n1 * n2, 1, 1) # p3
            W = (w1 # J(n2 * n3, 1, 1)) :* (J(n1, 1, 1) # w2 # J(n3, 1, 1)) :* (J(n1 * n2, 1, 1) # w3)
            base = A1 :- A2 :- A3
            AA = A1 :* (1 :- A1) :/ nc[1] :+ A2 :* (1 :- A2) :/ nc[2] :+ A3 :* (1 :- A3) :/ nc[3]
            T = rows(base)
            one = J(T, 1, 1)
            D = (one * p4' + base * J(1, n4, 1)) - z :* sqrt((n / (n - 4)) :* (AA * J(1, n4, 1) + one * (p4 :* (1 :- p4) :/ nc[4])'))
            D = (D :> 0)
            edge = (AA :== 0)
            D[., 1] = D[., 1] :* !edge
            D[., n4] = D[., n4] :* !edge
            outp[r] = quadcross(W, D * w4)
        }
        else {
            // the maximum-likelihood effect exists only when every cell mean
            // is strictly inside (0, 1)
            fail[r] = 1 - (1 - w1[1] - w1[rows(w1)]) * (1 - w2[1] - w2[rows(w2)]) * (1 - w3[1] - w3[rows(w3)]) * (1 - w4[1] - w4[rows(w4)])
            if (any(nc :< 2)) continue
            p1 = p1[2::nc[1]]; w1 = w1[2::nc[1]]
            p2 = p2[2::nc[2]]; w2 = w2[2::nc[2]]
            p3 = p3[2::nc[3]]; w3 = w3[2::nc[3]]
            p4 = p4[2::nc[4]]; w4 = w4[2::nc[4]]
            n1 = rows(p1); n2 = rows(p2); n3 = rows(p3); n4 = rows(p4)
            if (link == "probit") {
                l1 = invnormal(p1); l2 = invnormal(p2); l3 = invnormal(p3)
                v1 = p1 :* (1 :- p1) :/ (nc[1] :* normalden(l1):^2)
                v2 = p2 :* (1 :- p2) :/ (nc[2] :* normalden(l2):^2)
                v3 = p3 :* (1 :- p3) :/ (nc[3] :* normalden(l3):^2)
            }
            else {
                l1 = logit(p1); l2 = logit(p2); l3 = logit(p3)
                v1 = p1 :* (1 :- p1) :/ (nc[1] :* logisticden(l1):^2)
                v2 = p2 :* (1 :- p2) :/ (nc[2] :* logisticden(l2):^2)
                v3 = p3 :* (1 :- p3) :/ (nc[3] :* logisticden(l3):^2)
            }
            L = (J(n1 * n2, 1, 1) # l3) :+ (J(n1, 1, 1) # l2 # J(n3, 1, 1)) :- (l1 # J(n2 * n3, 1, 1))
            AA = (v1 # J(n2 * n3, 1, 1)) :+ (J(n1, 1, 1) # v2 # J(n3, 1, 1)) :+ (J(n1 * n2, 1, 1) # v3)
            W = (w1 # J(n2 * n3, 1, 1)) :* (J(n1, 1, 1) # w2 # J(n3, 1, 1)) :* (J(n1 * n2, 1, 1) # w3)
            if (link == "probit") {
                GL = normal(L)
                AA = normalden(L):^2 :* AA
            }
            else {
                GL = logistic(L)
                AA = logisticden(L):^2 :* AA
            }
            T = rows(L)
            one = J(T, 1, 1)
            D = (one * p4' - GL * J(1, n4, 1)) - z :* sqrt(one * (p4 :* (1 :- p4) :/ nc[4])' + AA * J(1, n4, 1))
            outp[r] = quadcross(W, (D :> 0) * w4)
        }
    }
    return((mean(outp), sqrt(variance(outp) / draws), mean(fail), sqrt(variance(fail) / draws)))
}
RES = J(0, 1, "")
NP = 0; NF = 0
void recx(string scalar label, real scalar a, real scalar b, real scalar tol)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar ok
    ok = (a != . & b != . & abs(a - b) <= tol)
    RES = RES \ sprintf("%-62s Stata %10.7f  R %10.7f  |diff| %8.1e  tol %7.1e  %s", label, a, b, abs(a - b), tol, (ok ? "PASS" : "FAIL"))
    if (ok) NP++
    else NF++
}
// simulated against exact: 4 x the two standard errors combined, + 1e-4
void reca(string scalar label, real scalar sim, real scalar ex, real scalar ex_se, real scalar nsim)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar mcse, tol, ok
    mcse = sqrt(max((ex * (1 - ex), 1e-4 * (1 - 1e-4))) / nsim)
    tol = 4 * sqrt(mcse^2 + ex_se^2) + 1e-4
    ok = (sim != . & abs(sim - ex) <= tol)
    RES = RES \ sprintf("%-62s exact %7.4f  simulated %7.4f  MCSE %6.4f  |diff| %6.4f  tol %6.4f  %s", label, ex, sim, mcse, abs(sim - ex), tol, (ok ? "PASS" : "FAIL"))
    if (ok) NP++
    else NF++
}
end

program define tiera
    args model scale seed
    powerape dgp a1, did model(`model') baseline(0.30) gap(0.10) trend(0.10) design(cross-sections) ///
        trends(`scale') analysis(`scale') nodisplay replace
    powerape setatt a1, target(0.30) nodisplay
    matrix RT = r(rates)
    local cs "A_`model'_`scale'"
    mata: pc = (st_matrix("RT")[1, 1], st_matrix("RT")[1, 2], st_matrix("RT")[1, 3], st_matrix("RT")[1, 5])
    mata: for (j = 1; j <= 4; j++) recx(sprintf("`model', `scale': cell rate %f equals R's", j), pc[j], refv("`cs'", sprintf("pc%f", j)), 1e-10)
    mata: ex = exact_did(pc, SIZES, "`model'", "`scale'", 0.95)
    mata: recx("`model', `scale': exact power, Mata oracle = R oracle", ex[1], refv("`cs'", "power"), 1e-9)
    mata: recx("`model', `scale': exact failed share, Mata oracle = R oracle", ex[3], refv("`cs'", "failed"), 1e-9)
    powerape power a1, n(120) claim(detect) nsim(40000) seed(`seed') engine(mata) nodisplay nodraws
    mata: o = st_matrix("r(outcomes)")
    mata: reca("`model', `scale': power, engine(mata), nsim 40000", st_numscalar("r(power)"), ex[1], ex[2], 40000)
    mata: reca("`model', `scale': failed share, engine(mata)", o[cols(o)], ex[3], ex[4], 40000)
    powerape power a1, n(120) claim(detect) nsim(4000) seed(`=`seed' + 1') engine(stata) nodisplay nodraws
    mata: o = st_matrix("r(outcomes)")
    mata: reca("`model', `scale': power, engine(stata), nsim 4000", st_numscalar("r(power)"), ex[1], ex[2], 4000)
    mata: reca("`model', `scale': failed share, engine(stata)", o[cols(o)], ex[3], ex[4], 4000)
end

tiera logit index 2701
tiera probit index 2703
tiera logit linear 2705
* margins in the loop (logit, index analysis)
powerape dgp a1, did model(logit) baseline(0.30) gap(0.10) trend(0.10) design(cross-sections) nodisplay replace
powerape setatt a1, target(0.30) nodisplay
mata: ex = (refv("A_logit_index", "power"), refv("A_logit_index", "power_se"))
powerape power a1, n(120) claim(detect) nsim(1500) seed(2707) engine(stata) effects(margins) nodisplay nodraws
mata: reca("logit, index: power, engine(stata) effects(margins), nsim 1500", st_numscalar("r(power)"), ex[1], ex[2], 1500)

mata:
printf("\n")
for (j = 1; j <= rows(RES); j++) printf("{txt}%s\n", RES[j])
printf("\n{txt}Tier A (difference-in-differences): %f PASS, %f FAIL\n", NP, NF)
unlink("did-tierA-results.txt")
fh = fopen("did-tierA-results.txt", "w")
for (j = 1; j <= rows(RES); j++) fput(fh, RES[j])
fput(fh, sprintf("Tier A (difference-in-differences): %f PASS, %f FAIL", NP, NF))
fclose(fh)
end
log close tier
cd "`root'"
exit, clear
