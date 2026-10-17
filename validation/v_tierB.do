* v_tierB.do -- fixed-data parity: the Stata port's estimators against the R
* package's internal estimators (identical algorithms: expected-information
* ML, closed-form delta-method Jacobians, HC0 / cluster sandwiches, the
* control-function stacked vcov) and against Stata's own commands
* (probit/logit + margins, ginteff, cfprobit + margins). Reference values
* come from validation/export_fixtures.R (powerape 1.11.3).
* Run from the repository root (the folder that holds ado/ and validation/):
*     do validation/v_tierB.do
* Writes validation/tierB-results.txt and the log validation/v_tierB.log.
version 19
clear all
local root "`c(pwd)'"
capture confirm file "`root'/ado/powerape.ado"
if _rc {
    di as err "v_tierB.do: run it from the repository root:  do validation/v_tierB.do"
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
log using v_tierB.log, text replace nomsg name(tier)
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
RES = J(0, 1, "")
RESN = J(0, 5, .)
void rec(string scalar cs, string scalar item, string scalar against, real scalar stata, real scalar r, real scalar tol)
{
    external string colvector RES
    external real matrix RESN
    real scalar d, rel
    d = abs(stata - r)
    rel = (r != 0 ? d / abs(r) : d)
    RES = RES \ sprintf("%-8s %-14s %-24s %14.8f %14.8f %10.2e %8.0e %s", cs, item, against, stata, r, d, tol, (d <= tol * max((1, abs(r))) ? "PASS" : "FAIL"))
    RESN = RESN \ (stata, r, d, tol, d <= tol * max((1, abs(r))))
}
end

* ---- A: probit, binary focal, two covariates -------------------------------
import delimited data/fixA.csv, clear varnames(1)
mata:
X = (J(st_nobs(), 1, 1), st_data(., ("d", "z1", "z2"))); y = st_data(., "y")
ft = pape_fit_index(X, y, "probit", J(0, 1, .))
e = pape_ape_est(ft.b, X, 2, "probit", "binary")
for (j = 1; j <= 4; j++) rec("A", sprintf("b%f", j), "R fit_index_model", ft.b[j], refv("A", sprintf("b%f", j)), 1e-6)
for (j = 1; j <= 4; j++) rec("A", sprintf("jac%f", j), "R ape_binary_est", e.jac[j], refv("A", sprintf("jac%f", j)), 1e-6)
rec("A", "ape", "R ape_binary_est", e.ape, refv("A", "ape"), 1e-7)
rec("A", "se_model", "R vcov_from_glmfit", sqrt(e.jac' * ft.V * e.jac), refv("A", "se_model"), 1e-6)
rec("A", "se_robust", "R vcov_sandwich HC0", sqrt(e.jac' * pape_vcov_sandwich(ft, X, y, "probit") * e.jac), refv("A", "se_robust"), 1e-6)
end
quietly probit y i.d z1 z2
quietly margins, dydx(d)
mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
mata: rec("A", "ape", "Stata probit+margins", b[2], refv("A", "ape"), 1e-6)
mata: rec("A", "se_model", "Stata margins (OIM)", sqrt(V[2, 2]), refv("A", "se_model"), 1e-2)
mata: Vm_oim = V
quietly probit y i.d z1 z2, vce(robust)
quietly margins, dydx(d)
mata: V = st_matrix("r(V)")
mata: rec("A", "se_robust", "Stata vce(robust) (OIM)", sqrt(V[2, 2]), refv("A", "se_robust"), 1e-2)
* the Stata engine's fit function (OIM) versus the direct command
mata: pape_frame_open(2000, 3, 0)
mata: fs = pape_fit_stata(X, y, J(0, 1, .), "probit", "oim", 0)
mata: rec("A", "b2 (engine stata)", "R fit_index_model", fs.b[2], refv("A", "b2"), 1e-5)
mata: es = pape_ape_est(fs.b, X, 2, "probit", "binary")
mata: rec("A", "ape (engine stata)", "R ape_binary_est", es.ape, refv("A", "ape"), 1e-6)
mata: rec("A", "se_model (engine stata)", "Stata margins (OIM)", sqrt(es.jac' * fs.V * es.jac), sqrt(Vm_oim[2, 2]), 1e-2)
mata: pape_frame_close()

* ---- B: logit, continuous focal --------------------------------------------
import delimited data/fixB.csv, clear varnames(1)
mata:
X = (J(st_nobs(), 1, 1), st_data(., ("x", "z1", "z2"))); y = st_data(., "y")
ft = pape_fit_index(X, y, "logit", J(0, 1, .))
e = pape_ape_est(ft.b, X, 2, "logit", "normal")
for (j = 1; j <= 4; j++) rec("B", sprintf("b%f", j), "R fit_index_model", ft.b[j], refv("B", sprintf("b%f", j)), 1e-6)
for (j = 1; j <= 4; j++) rec("B", sprintf("jac%f", j), "R ape_cont_est", e.jac[j], refv("B", sprintf("jac%f", j)), 1e-6)
rec("B", "ape", "R ape_cont_est", e.ape, refv("B", "ape"), 1e-7)
rec("B", "se_model", "R vcov_from_glmfit", sqrt(e.jac' * ft.V * e.jac), refv("B", "se_model"), 1e-6)
rec("B", "se_robust", "R vcov_sandwich HC0", sqrt(e.jac' * pape_vcov_sandwich(ft, X, y, "logit") * e.jac), refv("B", "se_robust"), 1e-6)
end
quietly logit y x z1 z2
quietly margins, dydx(x)
mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
mata: rec("B", "ape", "Stata logit+margins", b[1], refv("B", "ape"), 1e-6)
mata: rec("B", "se_model", "Stata margins (logit OIM=EIM)", sqrt(V[1, 1]), refv("B", "se_model"), 1e-5)
quietly logit y x z1 z2, vce(robust)
quietly margins, dydx(x)
mata: V = st_matrix("r(V)")
mata: rec("B", "se_robust", "Stata vce(robust) (N/(N-1))", sqrt(V[1, 1]), refv("B", "se_robust"), 1e-3)

* ---- C-F: AIE, four type pairs --------------------------------------------
foreach c in C D E F {
    import delimited data/fix`c'.csv, clear varnames(1)
    local link = cond("`c'" == "E", "logit", "probit")
    local ft = cond(inlist("`c'", "C", "E"), "binary", "normal")
    local mt = cond(inlist("`c'", "C", "D"), "binary", "normal")
    mata: X = (J(st_nobs(), 1, 1), st_data(., ("d", "m")), st_data(., "d") :* st_data(., "m"), st_data(., "z")); y = st_data(., "y")
    mata: ft = pape_fit_index(X, y, "`link'", J(0, 1, .))
    mata: e = pape_aie_est(ft.b, X, "`link'", "`ft'", "`mt'")
    mata: for (j = 1; j <= 5; j++) rec("`c'", sprintf("b%f", j), "R fit_index_model", ft.b[j], refv("`c'", sprintf("b%f", j)), 1e-6)
    mata: for (j = 1; j <= 5; j++) rec("`c'", sprintf("jac%f", j), "R aie_est", e.jac[j], refv("`c'", sprintf("jac%f", j)), 1e-6)
    mata: rec("`c'", "aie", "R aie_est", e.ape, refv("`c'", "aie"), 1e-7)
    mata: rec("`c'", "se_model", "R vcov_from_glmfit", sqrt(e.jac' * ft.V * e.jac), refv("`c'", "se_model"), 1e-6)
    local dterm = cond("`ft'" == "binary", "i.d", "c.d")
    local mterm = cond("`mt'" == "binary", "i.m", "c.m")
    quietly `link' y `dterm'##`mterm' z
    quietly ginteff, dydxs(d m)
    mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
    mata: rec("`c'", "aie", "Stata ginteff", b[1, 1], refv("`c'", "aie"), 1e-6)
    local tolse = cond("`link'" == "logit", "1e-5", "1e-2")
    mata: rec("`c'", "se_model", "Stata ginteff (`=cond("`link'"=="logit","OIM=EIM","OIM")')", sqrt(V[1, 1]), refv("`c'", "se_model"), `tolse')
}

* ---- G: clustered probit with a Mundlak mean --------------------------------
import delimited data/fixG.csv, clear varnames(1)
mata:
X = (J(st_nobs(), 1, 1), st_data(., ("d", "size", "sbar"))); y = st_data(., "y"); id = st_data(., "id")
ft = pape_fit_index(X, y, "probit", J(0, 1, .))
e = pape_ape_est(ft.b, X, 2, "probit", "binary")
for (j = 1; j <= 4; j++) rec("G", sprintf("b%f", j), "R fit_index_model", ft.b[j], refv("G", sprintf("b%f", j)), 1e-6)
rec("G", "ape", "R ape_binary_est", e.ape, refv("G", "ape"), 1e-7)
rec("G", "se_model", "R vcov_from_glmfit", sqrt(e.jac' * ft.V * e.jac), refv("G", "se_model"), 1e-6)
rec("G", "se_cluster", "R vcov_cluster G/(G-1)", sqrt(e.jac' * pape_vcov_cluster(ft, X, y, id, "probit") * e.jac), refv("G", "se_cluster"), 1e-6)
end
quietly probit y i.d size sbar, vce(cluster id)
quietly margins, dydx(d)
mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
mata: rec("G", "ape", "Stata probit+margins", b[2], refv("G", "ape"), 1e-6)
mata: rec("G", "se_cluster", "Stata vce(cluster) (OIM)", sqrt(V[2, 2]), refv("G", "se_cluster"), 2e-2)

* ---- H: control-function probit ----------------------------------------------
import delimited data/fixH.csv, clear varnames(1)
mata:
x1 = st_data(., "x1"); z1 = st_data(., "z1"); z2 = st_data(., "z2"); cid = st_data(., "cid")
y2 = st_data(., "y2"); y0 = st_data(., "y0"); y2b = st_data(., "y2b"); y0b = st_data(., "y0b")
m = st_data(., "m"); y0i = st_data(., "y0i")
X1 = (J(rows(x1), 1, 1), x1, z1, z2)
f = pape_cf_fit(y0, y2, J(0, 1, .), 0, x1, X1, "linear", J(0, 1, .), J(0, 1, .))
e = pape_cf_ape_est(f.theta, f.W, "normal", f.cf_col, ., 0)
for (j = 1; j <= 4; j++) rec("H_lin", sprintf("theta%f", j), "R cf_fit", f.theta[j], refv("H_lin", sprintf("theta%f", j)), 1e-6)
for (j = 1; j <= 4; j++) rec("H_lin", sprintf("se_theta%f", j), "R cf_fit stacked vcov", sqrt(f.V_main[j, j]), refv("H_lin", sprintf("se_theta%f", j)), 1e-6)
rec("H_lin", "ape", "R cf_ape_est", e.ape, refv("H_lin", "ape"), 1e-7)
rec("H_lin", "se", "R cf_ape_est", sqrt(e.jac' * f.V_main * e.jac), refv("H_lin", "se"), 1e-6)
fc = pape_cf_fit(y0, y2, J(0, 1, .), 0, x1, X1, "linear", cid, J(0, 1, .))
ec = pape_cf_ape_est(fc.theta, fc.W, "normal", fc.cf_col, ., 0)
rec("H_lin_cl", "se", "R cf_fit cluster", sqrt(ec.jac' * fc.V_main * ec.jac), refv("H_lin_cl", "se"), 1e-6)
fb = pape_cf_fit(y0b, y2b, J(0, 1, .), 0, x1, X1, "probit", J(0, 1, .), J(0, 1, .))
eb = pape_cf_ape_est(fb.theta, fb.W, "binary", fb.cf_col, ., 0)
for (j = 1; j <= 4; j++) rec("H_prb", sprintf("theta%f", j), "R cf_fit", fb.theta[j], refv("H_prb", sprintf("theta%f", j)), 1e-6)
for (j = 1; j <= 4; j++) rec("H_prb", sprintf("se_theta%f", j), "R cf_fit stacked vcov", sqrt(fb.V_main[j, j]), refv("H_prb", sprintf("se_theta%f", j)), 1e-6)
rec("H_prb", "ape", "R cf_ape_est", eb.ape, refv("H_prb", "ape"), 1e-7)
rec("H_prb", "se", "R cf_ape_est", sqrt(eb.jac' * fb.V_main * eb.jac), refv("H_prb", "se"), 1e-6)
X1m = (J(rows(x1), 1, 1), x1, m, z1, z2, z1 :* m, z2 :* m)
fi = pape_cf_fit(y0i, y2b, m, 1, x1, X1m, "probit", J(0, 1, .), J(0, 1, .))
ei = pape_cf_aie_est(fi.theta, fi.W, "binary", fi.cf_col, fi.cfm_col)
for (j = 1; j <= 7; j++) rec("H_aie", sprintf("theta%f", j), "R cf_fit", fi.theta[j], refv("H_aie", sprintf("theta%f", j)), 1e-6)
for (j = 1; j <= 7; j++) rec("H_aie", sprintf("se_theta%f", j), "R cf_fit stacked vcov", sqrt(fi.V_main[j, j]), refv("H_aie", sprintf("se_theta%f", j)), 1e-6)
rec("H_aie", "aie", "R cf_aie_est", ei.ape, refv("H_aie", "aie"), 1e-7)
rec("H_aie", "se", "R cf_aie_est", sqrt(ei.jac' * fi.V_main * ei.jac), refv("H_aie", "se"), 1e-6)
end
* Stata's own cfprobit + margins
quietly cfprobit y0 x1 (y2 = z1 z2, linear), vce(robust)
mata: b = st_matrix("e(b)"); V = st_matrix("e(V)")
mata: rec("H_lin", "theta2 (y2)", "Stata cfprobit", b[1], refv("H_lin", "theta2"), 1e-5)
mata: rec("H_lin", "theta4 (cf)", "Stata cfprobit", b[4], refv("H_lin", "theta4"), 1e-5)
mata: rec("H_lin", "se_theta2", "Stata cfprobit vce(robust)", sqrt(V[1, 1]), refv("H_lin", "se_theta2"), 1e-5)
mata: rec("H_lin", "se_theta4", "Stata cfprobit vce(robust)", sqrt(V[4, 4]), refv("H_lin", "se_theta4"), 1e-5)
quietly margins, dydx(y2)
mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
mata: rec("H_lin", "ape", "Stata cfprobit+margins", b[1], refv("H_lin", "ape"), 1e-6)
mata: rec("H_lin", "se", "Stata cfprobit+margins", sqrt(V[1, 1]), refv("H_lin", "se"), 1e-5)
quietly cfprobit y0 x1 (y2 = z1 z2, linear), vce(cluster cid)
quietly margins, dydx(y2)
mata: V = st_matrix("r(V)")
mata: rec("H_lin_cl", "se", "Stata cfprobit vce(cluster)+margins", sqrt(V[1, 1]), refv("H_lin_cl", "se"), 1e-5)
quietly cfprobit y0b x1 (y2b = z1 z2, probit), vce(robust)
mata: b = st_matrix("e(b)")
mata: rec("H_prb", "theta2 (1.y2b)", "Stata cfprobit", b[2], refv("H_prb", "theta2"), 1e-4)
quietly margins, dydx(y2b)
mata: b = st_matrix("r(b)"); V = st_matrix("r(V)")
mata: rec("H_prb", "ape", "Stata cfprobit+margins", b[2], refv("H_prb", "ape"), 1e-5)
mata: rec("H_prb", "se", "Stata cfprobit+margins", sqrt(V[2, 2]), refv("H_prb", "se"), 1e-5)
quietly cfprobit y0i x1 m (y2b = z1 z2 c.z1#c.m c.z2#c.m, probit interact(m)), mainonly(c.y2b#c.m) vce(robust)
mata: b = st_matrix("e(b)")
mata: rec("H_aie", "theta2 (1.y2b)", "Stata cfprobit interact()", b[2], refv("H_aie", "theta2"), 1e-4)
mata: rec("H_aie", "theta4 (y2b#m)", "Stata cfprobit interact()", b[3], refv("H_aie", "theta4"), 1e-4)
* the engine's cfprobit mapping on the same data reproduces the closed-form AIE
mata:
pape_frame_open(rows(x1), 2, 3)
d = pape_dgp_blank()
d.k = 1; d.has_m = 1; d.n_inst = 2; d.first_stage = "probit"; d.focal = pape_var_new("f", "binary", .4, ., ., .)
d.moderator = pape_var_new("m", "binary", .5, ., ., .)
xx = pape_xdraw()
xx.y = y0i; xx.d = y2b; xx.m = m; xx.Zx = x1; xx.X1 = X1m
cfs = pape_cf_fit_stata(d, xx, 0, 0, apem = ., sem = .)
es = pape_cf_aie_est(cfs.theta, cfs.W, "binary", cfs.cf_col, cfs.cfm_col)
rec("H_aie", "aie (engine stata)", "R cf_aie_est", es.ape, refv("H_aie", "aie"), 1e-5)
rec("H_aie", "se (engine stata)", "R cf_aie_est", sqrt(es.jac' * cfs.V_main * es.jac), refv("H_aie", "se"), 1e-5)
// 0.1.1: effects(margins) for the IV AIE = margins, dydx(d) at(m = 0 1) after cfprobit
// with the factor-form interaction (the control function held fixed: cfprobit's
// predict(pr) is the ASF). Against the R closed form on the same data.
cfm = pape_cf_fit_stata(d, xx, 0, 1, apem = ., sem = .)
rec("H_aie", "aie (margins dydx at)", "R cf_aie_est", apem, refv("H_aie", "aie"), 1e-6)
rec("H_aie", "se (margins dydx at)", "R cf_aie_est", sem, refv("H_aie", "se"), 1e-6)
// continuous focal (linear first stage) with interact(m): the R-port closed form on
// the same data is the reference (its APE counterparts are certified above, H_lin);
// y2, y0 were read from the fixture before the engine's temporary frame was opened
fl = pape_cf_fit(y0, y2, m, 1, x1, X1m, "linear", J(0, 1, .), J(0, 1, .))
el = pape_cf_aie_est(fl.theta, fl.W, "normal", fl.cf_col, fl.cfm_col)
d2 = pape_dgp_blank()
d2.k = 1; d2.has_m = 1; d2.n_inst = 2; d2.first_stage = "linear"
d2.focal = pape_var_new("f", "normal", ., 0, 1, .)
d2.moderator = pape_var_new("m", "binary", .5, ., ., .)
xx2 = pape_xdraw()
xx2.y = y0; xx2.d = y2; xx2.m = m; xx2.Zx = x1; xx2.X1 = X1m
cfs2 = pape_cf_fit_stata(d2, xx2, 0, 0, apem2 = ., sem2 = .)
es2 = pape_cf_aie_est(cfs2.theta, cfs2.W, "normal", cfs2.cf_col, cfs2.cfm_col)
rec("H_aie_c", "aie (engine stata)", "R port cf_aie_est", es2.ape, el.ape, 1e-5)
rec("H_aie_c", "se (engine stata)", "R port cf_aie_est", sqrt(es2.jac' * cfs2.V_main * es2.jac), sqrt(el.jac' * fl.V_main * el.jac), 1e-5)
cfm2 = pape_cf_fit_stata(d2, xx2, 0, 1, apem2 = ., sem2 = .)
rec("H_aie_c", "aie (margins dydx at)", "R port cf_aie_est", apem2, el.ape, 1e-6)
rec("H_aie_c", "se (margins dydx at)", "R port cf_aie_est", sem2, sqrt(el.jac' * fl.V_main * el.jac), 1e-6)
pape_frame_close()
end

* ---- S: separation(keep) under the Stata engine (0.1.1) ------------------------
* A saturated design whose focal arm has no events: Stata's asis fit returns the
* separated coefficient with a zero variance row; the delta method then collapses
* to the identified cell, the limit of R's glm (engine(mata)). Reference = the Mata
* engine on the same data (expected information; OIM differs at the 1e-3 level).
clear
mata:
rseed(20261016)
n = 300
t = (1::n) :> 150
y = (runiform(n, 1) :< 0.06) :& !t
z = rnormal(n, 1, 0, 1)
X = (J(n, 1, 1), t)
fm = pape_fit_index(X, y, "probit", J(0, 1, .))
em = pape_ape_est(fm.b, X, 2, "probit", "binary")
se_m = sqrt(em.jac' * fm.V * em.jac)
p0 = mean(select(y, !t)); n0 = sum(!t)
se_bin = sqrt(p0 * (1 - p0) / n0)
pape_frame_open(n, 2, 0)
fs = pape_fit_stata(X, y, J(0, 1, .), "probit", "oim", 1)
rec("sep_sat", "fit ok (asis)", "engine mata ok", fs.ok, fm.ok, 0)
rec("sep_sat", "V[t,t] zero kept", "rank k-1", fs.rank, cols(X) - 1, 0)
es = pape_ape_est(fs.b, X, 2, "probit", "binary")
se_s = sqrt(es.jac' * fs.V * es.jac)
// Stata's asis fit stops on its ll-change criterion with the intercept about 1e-4 off
// the closed form (the likelihood is flat along the separated coefficient): APE to 1e-5,
// negligible against the SE of .019; R's IRLS (engine mata) converges the intercept exactly
rec("sep_sat", "ape (asis)", "engine mata (R glm)", es.ape, em.ape, 1e-4)
rec("sep_sat", "ape (asis)", "-p0hat", es.ape, -p0, 1e-4)
rec("sep_sat", "se (asis)", "engine mata (R glm)", se_s, se_m, 1e-2)
rec("sep_sat", "se (asis)", "other-cell binomial", se_s, se_bin, 1e-2)
// the same without asis fails (the predictor and its observations are dropped)
fs0 = pape_fit_stata(X, y, J(0, 1, .), "probit", "oim", 0)
rec("sep_sat", "fit ok (no asis)", "fails", fs0.ok, 0, 0)
// with a covariate the asis fit has a finite, huge variance: kept before 0.1.1 too
Xz = (J(n, 1, 1), t, z)
fmz = pape_fit_index(Xz, y, "probit", J(0, 1, .))
emz = pape_ape_est(fmz.b, Xz, 2, "probit", "binary")
pape_frame_resize(n)
fsz = pape_fit_stata(Xz, y, J(0, 1, .), "probit", "oim", 1)
esz = pape_ape_est(fsz.b, Xz, 2, "probit", "binary")
rec("sep_cov", "fit ok (asis)", "engine mata ok", fsz.ok, fmz.ok, 0)
rec("sep_cov", "ape (asis)", "engine mata (R glm)", esz.ape, emz.ape, 1e-4)
rec("sep_cov", "se (asis)", "engine mata (R glm)", sqrt(esz.jac' * fsz.V * esz.jac), sqrt(emz.jac' * fmz.V * emz.jac), 1e-2)
// logit: same convention
fml = pape_fit_index(X, y, "logit", J(0, 1, .))
eml = pape_ape_est(fml.b, X, 2, "logit", "binary")
pape_frame_resize(n)
fsl = pape_fit_stata(X, y, J(0, 1, .), "logit", "oim", 1)
esl = pape_ape_est(fsl.b, X, 2, "logit", "binary")
rec("sep_sat", "ape logit (asis)", "engine mata (R glm)", esl.ape, eml.ape, 1e-6)
rec("sep_sat", "se logit (asis)", "engine mata (R glm)", sqrt(esl.jac' * fsl.V * esl.jac), sqrt(eml.jac' * fml.V * eml.jac), 1e-2)
pape_frame_close()
end

* ---- closed-form helpers ------------------------------------------------------
mata:
rec("util", "req_delta_1", "R req_delta", pape_req_delta(0.8, invnormal(0.95), 1), refv("util", "req_delta_1"), 1e-9)
rec("util", "req_delta_2", "R req_delta", pape_req_delta(0.9, invnormal(0.95), 3), refv("util", "req_delta_2"), 1e-9)
rec("util", "req_margin", "R req_margin", pape_req_margin(0.8, 0.02, 0.01, invnormal(0.95)), refv("util", "req_margin"), 1e-9)
rec("util", "frontier", "R alpha_frontier", pape_alpha_frontier(4, 0.8, 1, "detect"), refv("util", "frontier"), 1e-9)
rec("util", "frontier_min", "R alpha_frontier", pape_alpha_frontier(1, 0.8, 1, "minimum"), refv("util", "frontier_min"), 1e-9)
rec("util", "mz_r2", "R mz_r2_logit", pape_mz_r2_logit(0.1 / 0.9), refv("util", "mz_r2"), 1e-12)
end

* ---- report -------------------------------------------------------------------
mata:
printf("\n{txt}%-8s %-14s %-24s %14s %14s %10s %8s %s\n", "case", "item", "against", "Stata", "reference", "|diff|", "tol", "")
for (j = 1; j <= rows(RES); j++) printf("%s\n", RES[j])
printf("\n{txt}Tier B: %f comparisons, %f PASS, %f FAIL\n", rows(RESN), sum(RESN[., 5]), rows(RESN) - sum(RESN[., 5]))
unlink("tierB-results.txt")
fh = fopen("tierB-results.txt", "w")
for (j = 1; j <= rows(RES); j++) fput(fh, RES[j])
fput(fh, sprintf("Tier B: %f comparisons, %f PASS, %f FAIL", rows(RESN), sum(RESN[., 5]), rows(RESN) - sum(RESN[., 5])))
fclose(fh)
end
log close tier
cd "`root'"
exit, clear
