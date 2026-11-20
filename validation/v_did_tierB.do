* v_did_tierB.do -- difference-in-differences route, fixed-data parity: the
* port's estimators against the R package's internal estimators (same
* algorithms: expected-information ML, the effect on the treated with its
* delta-method Jacobian, HC0 / cluster sandwiches, OLS with HC1 / cluster
* variances) and against Stata's own commands (logit/probit + margins,
* regress, didregress, xtdidregress). Reference values come from
* validation/did_export_fixtures.R, run on version 1.12.1 of the R package.
* Run from the repository root (the folder that holds ado/ and validation/):
*     do validation/v_did_tierB.do
* Writes validation/did-tierB-results.txt and the log validation/v_did_tierB.log.
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "v_did_tierB.do: run it from the repository root:  do validation/v_did_tierB.do"
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
log using v_did_tierB.log, text replace nomsg name(tier)
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
RES = J(0, 1, "")
NP = 0; NF = 0
void rec(string scalar cs, string scalar item, string scalar against, real scalar stata, real scalar r, real scalar tol)
{
    external string colvector RES
    external real scalar NP, NF
    real scalar d, ok
    d = abs(stata - r)
    ok = (stata != . & r != . & d <= tol * max((1, abs(r))))
    RES = RES \ sprintf("%-5s %-15s %-34s %14.8f %14.8f %10.2e %8.0e %s", cs, item, against, stata, r, d, tol, (ok ? "PASS" : "FAIL"))
    if (ok) NP++
    else NF++
}
end

* One fixture through every path of the port. Arguments: case, link, number
* of periods, the standard-error kind of the index analysis in cross-sections
* (model or robust; clustered fixtures ignore it).
program define tierb
    args cs link T sekind
    import delimited data/didfix_`cs'.csv, clear varnames(1) asdouble
    mata: K = refv("`cs'", "k"); fc = refv("`cs'", "focal_col")
    mata: X = (J(st_nobs(), 1, 1), st_data(., "x" :+ strofreal(2..K))); y = st_data(., "y")
    mata: id = st_data(., "id"); msub = st_data(., "m")
    mata: id = (all(id :== 0) ? J(0, 1, .) : id)
    mata: msub = (refv("`cs'", "by0") == . ? J(0, 1, .) : msub)
    mata: sen = (rows(id) > 0 ? "se_cluster" : "se_`sekind'")
    mata: lsen = (rows(id) > 0 ? "lin_se_cluster" : "lin_se_hc1")
    mata: est = .; se = .; b = .; jac = .; by = .
    mata: rec("`cs'", "sum_y", "the data R wrote", sum(y), refv("`cs'", "sum_y"), 0)
    * engine(mata), index analysis: coefficients, estimate, Jacobian, SE
    mata: ok = pape_did_analyze(X, y, id, msub, "`link'", `T', "index", "`sekind'", "mata", "closed", est, se, b, jac, by)
    mata: for (j = 1; j <= K; j++) rec("`cs'", sprintf("b%f", j), "R fit_index_model", b[j], refv("`cs'", sprintf("b%f", j)), 1e-6)
    mata: for (j = 1; j <= K; j++) rec("`cs'", sprintf("jac%f", j), "R att_est / datt_est", jac[j], refv("`cs'", sprintf("jac%f", j)), 1e-6)
    mata: rec("`cs'", "est", "R, engine(mata)", est, refv("`cs'", "est"), 1e-7)
    mata: rec("`cs'", sen, "R, engine(mata)", se, refv("`cs'", sen), 1e-6)
    mata: for (j = 1; j <= 2 * (cols(by) == 2); j++) rec("`cs'", sprintf("by%f", j - 1), "R datt_est", by[j], refv("`cs'", sprintf("by%f", j - 1)), 1e-7)
    mata: est_m = est; se_m = se
    * engine(stata), effects(closed): Stata's fit, closed-form effect. Stata's
    * probit variance uses the observed information and R's the expected, so
    * probit SEs are compared at 1e-2; logit SEs coincide
    * (Stata's vce(robust) also carries the factor N/(N-1), which the R sandwich
    * does not: 1e-3 for those rows)
    mata: tolse = ("`link'" == "probit" ? 1e-2 : (rows(id) == 0 & "`sekind'" == "robust" ? 1e-3 : 1e-5))
    mata: ok = pape_did_analyze(X, y, id, msub, "`link'", `T', "index", "`sekind'", "stata", "closed", est, se, b, jac, by)
    mata: rec("`cs'", "est", "R, engine(stata) closed", est, refv("`cs'", "est"), 1e-5)
    mata: rec("`cs'", sen, "R, engine(stata) closed", se, refv("`cs'", sen), tolse)
    mata: est_c = est; se_c = se
    * engine(stata), effects(margins): margins, dydx() subpop() [over()]
    mata: ok = pape_did_analyze(X, y, id, msub, "`link'", `T', "index", "`sekind'", "stata", "margins", est, se, b, jac, by)
    mata: rec("`cs'", "est", "R, engine(stata) margins", est, refv("`cs'", "est"), 1e-5)
    mata: rec("`cs'", sen, "R, engine(stata) margins", se, refv("`cs'", sen), tolse)
    mata: rec("`cs'", "est", "margins = closed form", est, est_c, 1e-7)
    mata: rec("`cs'", sen, "margins = closed form", se, se_c, 1e-6)
    * linear analysis: OLS in Mata and Stata's regress
    mata: ok = pape_did_analyze(X, y, id, msub, "`link'", `T', "linear", "model", "mata", "closed", est, se, b, jac, by)
    mata: rec("`cs'", "lin_b", "R fit_ols, engine(mata)", est, refv("`cs'", "lin_b"), 1e-9)
    mata: rec("`cs'", lsen, "R vcov_ols, engine(mata)", se, refv("`cs'", lsen), 1e-9)
    mata: ok = pape_did_analyze(X, y, id, msub, "`link'", `T', "linear", "model", "stata", "closed", est, se, b, jac, by)
    mata: rec("`cs'", "lin_b", "R fit_ols, engine(stata) regress", est, refv("`cs'", "lin_b"), 1e-9)
    mata: rec("`cs'", lsen, "R vcov_ols, engine(stata) regress", se, refv("`cs'", lsen), 1e-9)
end

tierb did logit 2 model
tierb didx logit 4 model
tierb didx logit 4 robust
tierb sub logit 2 model
tierb rp probit 4 model
tierb rs logit 3 model
tierb rs logit 3 robust

* ---- Stata's own commands on the fixtures with named variables -------------
* did: panel, 400 units, one period before and after, covariate z
import delimited data/didfix_did.csv, clear varnames(1) asdouble
quietly logit y i.d i.period i.w z, vce(cluster id)
quietly margins, dydx(w) subpop(if w == 1)
mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
mata: rec("did", "est", "logit + margins subpop()", b[cols(b)], refv("did", "est"), 1e-5)
mata: rec("did", "se_cluster", "logit + margins subpop()", sqrt(V[cols(b), cols(b)]), refv("did", "se_cluster"), 1e-5)
quietly regress y i.d i.period i.w z, vce(cluster id)
scalar sb = _b[1.w]
mata: rec("did", "lin_b", "regress, vce(cluster)", st_numscalar("sb"), refv("did", "lin_b"), 1e-9)
scalar sse = _se[1.w]
mata: rec("did", "lin_se_cluster", "regress, vce(cluster)", st_numscalar("sse"), refv("did", "lin_se_cluster"), 1e-9)
xtset id period
quietly xtdidregress (y z) (w), group(id) time(period)
mata: b = st_matrix("e(b)"); V = st_matrix("e(V)"); N = refv("did", "n")
mata: rec("did", "lin_b", "xtdidregress", b[1], refv("did", "lin_b"), 1e-9)
* xtdidregress counts the within regressors (period, w, constant), the pooled
* regression all five columns: its SE differs by that degrees-of-freedom factor
mata: rec("did", "lin_se_cluster", "xtdidregress x sqrt((N-3)/(N-5))", sqrt(V[1, 1]) * sqrt((N - 3) / (N - 5)), refv("did", "lin_se_cluster"), 1e-9)

* didx: repeated cross-sections, four periods, treatment from period 3
import delimited data/didfix_didx.csv, clear varnames(1) asdouble
quietly logit y i.d i.period i.w
quietly margins, dydx(w) subpop(if w == 1)
mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
mata: rec("didx", "est", "logit + margins subpop()", b[cols(b)], refv("didx", "est"), 1e-5)
mata: rec("didx", "se_model", "logit + margins subpop()", sqrt(V[cols(b), cols(b)]), refv("didx", "se_model"), 1e-5)
quietly logit y i.d i.period i.w, vce(robust)
quietly margins, dydx(w) subpop(if w == 1)
mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
* Stata's vce(robust) carries the factor N/(N-1) on the variance; the R
* sandwich (HC0) does not
mata: N = refv("didx", "n")
mata: rec("didx", "se_robust", "logit vce(robust) x sqrt((N-1)/N)", sqrt(V[cols(b), cols(b)]) * sqrt((N - 1) / N), refv("didx", "se_robust"), 1e-5)
quietly regress y i.d i.period i.w, vce(robust)
scalar sb = _b[1.w]
mata: rec("didx", "lin_b", "regress, vce(robust)", st_numscalar("sb"), refv("didx", "lin_b"), 1e-9)
scalar sse = _se[1.w]
mata: rec("didx", "lin_se_hc1", "regress, vce(robust)", st_numscalar("sse"), refv("didx", "lin_se_hc1"), 1e-9)
quietly didregress (y) (w), group(d) time(period) vce(robust)
mata: b = st_matrix("e(b)"); V = st_matrix("e(V)")
mata: rec("didx", "lin_b", "didregress, vce(robust)", b[1], refv("didx", "lin_b"), 1e-9)
mata: rec("didx", "lin_se_hc1", "didregress, vce(robust)", sqrt(V[1, 1]), refv("didx", "lin_se_hc1"), 1e-9)

* sub: panel with a binary subgroup m (one period before and after)
import delimited data/didfix_sub.csv, clear varnames(1) asdouble
generate double md = m * d
generate double mpost = m * post
generate double mw = m * w
quietly logit y i.w##i.m d post md mpost, vce(cluster id)
quietly margins, dydx(w) subpop(if w == 1) over(m) post
matrix bb = e(b)
local nm : colfullnames bb
local k = colsof(bb)
local k0 = `k' - 1
local n0 : word `k0' of `nm'
local n1 : word `k' of `nm'
mata: rec("sub", "by0", "margins over(m)", st_matrix("e(b)")[`k0'], refv("sub", "by0"), 1e-5)
mata: rec("sub", "by1", "margins over(m)", st_matrix("e(b)")[`k'], refv("sub", "by1"), 1e-5)
quietly lincom _b[`n1'] - _b[`n0']
mata: rec("sub", "est", "margins over(m) + lincom", st_numscalar("r(estimate)"), refv("sub", "est"), 1e-5)
mata: rec("sub", "se_cluster", "margins over(m) + lincom", st_numscalar("r(se)"), refv("sub", "se_cluster"), 1e-5)
quietly regress y d post w m md mpost mw, vce(cluster id)
scalar sb = _b[mw]
mata: rec("sub", "lin_b", "regress, vce(cluster)", st_numscalar("sb"), refv("sub", "lin_b"), 1e-9)
scalar sse = _se[mw]
mata: rec("sub", "lin_se_cluster", "regress, vce(cluster)", st_numscalar("sse"), refv("sub", "lin_se_cluster"), 1e-9)
* the official triple difference (exact with two periods)
quietly didregress (y) (mw), group(d m) time(period)
mata: b = st_matrix("e(b)")
mata: rec("sub", "lin_b", "didregress, group(d m) (DDD)", b[1], refv("sub", "lin_b"), 1e-9)

mata:
printf("\n")
for (j = 1; j <= rows(RES); j++) printf("{txt}%s\n", RES[j])
printf("\n{txt}Tier B (difference-in-differences): %f PASS, %f FAIL\n", NP, NF)
unlink("did-tierB-results.txt")
fh = fopen("did-tierB-results.txt", "w")
fput(fh, sprintf("%-5s %-15s %-34s %14s %14s %10s %8s", "case", "item", "against", "Stata", "reference", "|diff|", "tol"))
for (j = 1; j <= rows(RES); j++) fput(fh, RES[j])
fput(fh, sprintf("Tier B (difference-in-differences): %f PASS, %f FAIL", NP, NF))
fclose(fh)
end
log close tier
cd "`root'"
exit, clear
