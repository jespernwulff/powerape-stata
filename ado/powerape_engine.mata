*! powerape_engine.mata -- simulation engine: ML fit, APE/AIE estimators,
*! delta-method covariances, separation checks, claims and outcomes.
*! Mirrors R/ape-math.R and R/engine.R of powerape 1.11.3.
version 16
mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// Structures
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Claim levels (engine.R: claim_levels). A sample-size rule is a Stata
// expression in @n evaluated at n.
// ---------------------------------------------------------------------------

real scalar pape_eval_rule(string scalar rule, real scalar n)
{
    string scalar ex
    real scalar a, rc
    ex = subinstr(rule, "@n", strofreal(n, "%21.0g"), .)
    rc = _stata("capture scalar PAPE_RULE_TMP = " + ex, 1)
    if (rc) pape_err(sprintf("The alpha rule '%s' could not be evaluated at n = %f.", rule, n))
    a = st_numscalar("PAPE_RULE_TMP")
    stata("capture scalar drop PAPE_RULE_TMP", 1)
    if (a == . | a <= 0 | a >= 0.5)
        pape_err(sprintf("The alpha rule must return a single number in (0, 0.5); at n = %f it returned %f.", n, a))
    return(a)
}

// returns (alpha, conf_det, conf_one, conf_claim)
real rowvector pape_claim_levels(string scalar claim, real scalar alpha,
    real scalar conf, real scalar n, string scalar rule)
{
    real scalar a
    if (rule != "") {
        if (conf != .) pape_err("Supply either an alpha rule or conf(), not both.")
        if (n == .) pape_err("A sample-size rule for alpha needs n to be evaluated.")
        a = pape_eval_rule(rule, n)
    }
    else if (conf != .) {
        if (conf <= 0.5 | conf >= 1) pape_err("conf() must lie in (0.5, 1).")
        a = (claim == "detect" ? 1 - conf : (1 - conf) / 2)
    }
    else {
        a = (alpha == . ? 0.05 : alpha)
        if (a <= 0 | a >= 0.5) pape_err("alpha() must lie in (0, 0.5).")
    }
    return((a, 1 - a, 1 - 2 * a, (claim == "detect" ? 1 - a : 1 - 2 * a)))
}

// ---------------------------------------------------------------------------
// Maximum-likelihood fit by Fisher scoring (IRLS), as R's glm.fit
// ---------------------------------------------------------------------------

struct pape_fit scalar pape_fit_index(real matrix X, real colvector y,
    string scalar link, real colvector start)
{
    struct pape_fit scalar f
    real colvector eta, mu, gd, w, z, b, bold, mustart, varmu
    real matrix XWX, V
    real scalar dev, devold, it, k, n, j, sdj, runaway, converged, eps
    n = rows(X); k = cols(X)
    eps = 1e-10
    f.ok = 0; f.converged = 0; f.rank = k
    if (rows(start) == k & !hasmissing(start)) {
        b = start
        eta = X * b
    }
    else {
        mustart = (y :+ 0.5) :/ 2
        eta = (link == "probit" ? invnormal(mustart) : logit(mustart))
        b = J(k, 1, 0)
    }
    devold = .
    converged = 0
    for (it = 1; it <= 100; it++) {
        mu = pape_G(link, eta)
        mu = rowmax((rowmin((mu, J(n, 1, 1 - eps))), J(n, 1, eps)))
        gd = pape_g(link, eta)
        varmu = mu :* (1 :- mu)
        w = gd:^2 :/ varmu
        z = eta :+ (y :- mu) :/ gd
        XWX = cross(X, w, X)
        V = invsym(XWX)
        if (any(diagonal(V) :== 0)) {
            f.rank = k - sum(diagonal(V) :== 0)
            f.b = b
            return(f)
        }
        bold = b
        b = V * cross(X, w, z)
        eta = X * b
        mu = pape_G(link, eta)
        mu = rowmax((rowmin((mu, J(n, 1, 1 - eps))), J(n, 1, eps)))
        dev = -2 * sum(y :* ln(mu) :+ (1 :- y) :* ln(1 :- mu))
        if (dev == .) {
            // step-halve back toward the previous coefficients
            for (j = 1; j <= 20; j++) {
                b = (b + bold) / 2
                eta = X * b
                mu = pape_G(link, eta)
                mu = rowmax((rowmin((mu, J(n, 1, 1 - eps))), J(n, 1, eps)))
                dev = -2 * sum(y :* ln(mu) :+ (1 :- y) :* ln(1 :- mu))
                if (dev != .) break
            }
            if (dev == .) break
        }
        if (devold != . & abs(dev - devold) / (abs(dev) + 0.1) < 1e-8) {
            converged = 1
            break
        }
        devold = dev
    }
    f.converged = converged
    f.b = b
    f.eta = eta
    // covariance at the converged coefficients (expected information)
    mu = pape_G(link, eta)
    mu = rowmax((rowmin((mu, J(n, 1, 1 - eps))), J(n, 1, eps)))
    gd = pape_g(link, eta)
    w = gd:^2 :/ (mu :* (1 :- mu))
    f.V = invsym(cross(X, w, X))
    f.mu = mu
    if (any(diagonal(f.V) :== 0)) f.rank = k - sum(diagonal(f.V) :== 0)
    runaway = 0
    if (k > 1 & !hasmissing(b)) {
        for (j = 2; j <= k; j++) {
            sdj = sqrt(variance(X[., j]))
            if (abs(b[j]) * max((sdj, 1e-12)) > 15) runaway = 1
        }
    }
    f.ok = (converged & f.rank == k & !hasmissing(b) & !runaway)
    return(f)
}

// Heteroskedasticity-robust (HC0) sandwich, expected-information bread.
real matrix pape_vcov_sandwich(struct pape_fit scalar f, real matrix X,
    real colvector y, string scalar link)
{
    real colvector mu, w
    real matrix B
    mu = pape_G(link, f.eta)
    w = pape_g(link, f.eta) :* (y :- mu) :/ rowmax((mu :* (1 :- mu), J(rows(mu), 1, 1e-12)))
    B = f.V
    return(B * cross(X :* w, X :* w) * B)
}

// Cluster-robust sandwich with the G/(G-1) multiplier (Stata's vce(cluster)).
real matrix pape_vcov_cluster(struct pape_fit scalar f, real matrix X,
    real colvector y, real colvector id, string scalar link)
{
    real colvector mu, w, ord
    real matrix B, U, info, S
    real scalar Gn
    mu = pape_G(link, f.eta)
    w = pape_g(link, f.eta) :* (y :- mu) :/ rowmax((mu :* (1 :- mu), J(rows(mu), 1, 1e-12)))
    S = X :* w
    ord = order(id, 1)
    info = panelsetup(id[ord], 1)
    U = panelsum(S[ord, .], info)
    B = f.V
    Gn = rows(U)
    return((Gn / (Gn - 1)) :* (B * cross(U, U) * B))
}

// ---------------------------------------------------------------------------
// APE estimators and Jacobians (ape-math.R)
// ---------------------------------------------------------------------------

struct pape_est scalar pape_ape_binary_est(real colvector beta, real matrix X,
    real scalar focal_col, string scalar link)
{
    struct pape_est scalar e
    real matrix X1, X0
    real colvector a1, a0
    X1 = X; X1[., focal_col] = J(rows(X), 1, 1)
    X0 = X; X0[., focal_col] = J(rows(X), 1, 0)
    a1 = X1 * beta; a0 = X0 * beta
    e.ape = mean(pape_G(link, a1) :- pape_G(link, a0))
    e.jac = (mean(pape_g(link, a1) :* X1) :- mean(pape_g(link, a0) :* X0))'
    return(e)
}

struct pape_est scalar pape_ape_cont_est(real colvector beta, real matrix X,
    real scalar focal_col, string scalar link)
{
    struct pape_est scalar e
    real colvector a, ga
    real scalar b_d
    a = X * beta
    ga = pape_g(link, a)
    b_d = beta[focal_col]
    e.ape = b_d * mean(ga)
    e.jac = (b_d :* mean(pape_gprime(link, a) :* X))'
    e.jac[focal_col] = e.jac[focal_col] + mean(ga)
    return(e)
}

struct pape_est scalar pape_ape_est(real colvector beta, real matrix X,
    real scalar focal_col, string scalar link, string scalar focal_type)
{
    if (focal_type == "binary") return(pape_ape_binary_est(beta, X, focal_col, link))
    return(pape_ape_cont_est(beta, X, focal_col, link))
}

// AIE: columns 2 = focal, 3 = moderator, 4 = interaction
struct pape_est scalar pape_aie_est(real colvector beta, real matrix X,
    string scalar link, string scalar focal_type, string scalar mod_type)
{
    struct pape_est scalar e
    real matrix X11, X01, X10, X00, X1, X0
    real colvector a11, a01, a10, a00, a1, a0, a, dv, mv, t1, t2, gpa, ga
    real scalar bd, bm, bi, n
    n = rows(X)
    if (focal_type == "binary" & mod_type == "binary") {
        X11 = X; X11[., 2] = J(n, 1, 1); X11[., 3] = J(n, 1, 1); X11[., 4] = J(n, 1, 1)
        X01 = X; X01[., 2] = J(n, 1, 0); X01[., 3] = J(n, 1, 1); X01[., 4] = J(n, 1, 0)
        X10 = X; X10[., 2] = J(n, 1, 1); X10[., 3] = J(n, 1, 0); X10[., 4] = J(n, 1, 0)
        X00 = X; X00[., 2] = J(n, 1, 0); X00[., 3] = J(n, 1, 0); X00[., 4] = J(n, 1, 0)
        a11 = X11 * beta; a01 = X01 * beta; a10 = X10 * beta; a00 = X00 * beta
        e.ape = mean(pape_G(link, a11) :- pape_G(link, a01) :- pape_G(link, a10) :+ pape_G(link, a00))
        e.jac = (mean(pape_g(link, a11) :* X11) :- mean(pape_g(link, a01) :* X01) :-
                 mean(pape_g(link, a10) :* X10) :+ mean(pape_g(link, a00) :* X00))'
        return(e)
    }
    if (focal_type != "binary" & mod_type == "binary") {
        bd = beta[2]; bi = beta[4]
        X1 = X; X1[., 3] = J(n, 1, 1); X1[., 4] = X[., 2]
        X0 = X; X0[., 3] = J(n, 1, 0); X0[., 4] = J(n, 1, 0)
        a1 = X1 * beta; a0 = X0 * beta
        e.ape = (bd + bi) * mean(pape_g(link, a1)) - bd * mean(pape_g(link, a0))
        e.jac = ((bd + bi) :* mean(pape_gprime(link, a1) :* X1) :- bd :* mean(pape_gprime(link, a0) :* X0))'
        e.jac[2] = e.jac[2] + mean(pape_g(link, a1)) - mean(pape_g(link, a0))
        e.jac[4] = e.jac[4] + mean(pape_g(link, a1))
        return(e)
    }
    if (focal_type == "binary" & mod_type != "binary") {
        bm = beta[3]; bi = beta[4]
        X1 = X; X1[., 2] = J(n, 1, 1); X1[., 4] = X[., 3]
        X0 = X; X0[., 2] = J(n, 1, 0); X0[., 4] = J(n, 1, 0)
        a1 = X1 * beta; a0 = X0 * beta
        e.ape = (bm + bi) * mean(pape_g(link, a1)) - bm * mean(pape_g(link, a0))
        e.jac = ((bm + bi) :* mean(pape_gprime(link, a1) :* X1) :- bm :* mean(pape_gprime(link, a0) :* X0))'
        e.jac[3] = e.jac[3] + mean(pape_g(link, a1)) - mean(pape_g(link, a0))
        e.jac[4] = e.jac[4] + mean(pape_g(link, a1))
        return(e)
    }
    bd = beta[2]; bm = beta[3]; bi = beta[4]
    a = X * beta
    dv = X[., 2]; mv = X[., 3]
    t1 = bd :+ bi :* mv
    t2 = bm :+ bi :* dv
    gpa = pape_gprime(link, a); ga = pape_g(link, a)
    e.ape = mean(gpa :* t1 :* t2 :+ ga :* bi)
    e.jac = (mean(pape_gpp(link, a) :* t1 :* t2 :* X) :+ bi :* mean(gpa :* X))'
    e.jac[2] = e.jac[2] + mean(gpa :* t2)
    e.jac[3] = e.jac[3] + mean(gpa :* t1)
    e.jac[4] = e.jac[4] + mean(gpa :* (mv :* t2 :+ dv :* t1)) + mean(ga)
    return(e)
}

// ---------------------------------------------------------------------------
// Separation checks (engine.R: cell_check, hauck_donner)
// ---------------------------------------------------------------------------

// returns (sep, min_cell); b1/b2 empty (rows 0) when absent
real rowvector pape_cell_check(real colvector y, real colvector b1, real colvector b2)
{
    real colvector g, nn, ev, keep
    real scalar j, n, has1, has2
    has1 = rows(b1) > 0; has2 = rows(b2) > 0
    if (has1 & has2) g = 2 :* b1 :+ b2
    else if (has1) g = b1
    else if (has2) g = b2
    else return((0, .))
    n = rows(g)
    nn = J(4, 1, 0); ev = J(4, 1, 0)
    for (j = 1; j <= n; j++) {
        nn[g[j] + 1] = nn[g[j] + 1] + 1
        if (y[j] == 1) ev[g[j] + 1] = ev[g[j] + 1] + 1
    }
    keep = selectindex(nn :> 0)
    nn = nn[keep]; ev = ev[keep]
    return((any(ev :== 0 :| ev :== nn), min(rowmin((ev, nn :- ev)))))
}

real scalar pape_hauck_donner(struct pape_fit scalar f, real matrix X, real rowvector cols)
{
    real scalar j, c, sdj
    for (j = 1; j <= cols(cols); j++) {
        c = cols[j]
        sdj = sqrt(variance(X[., c]))
        if (sqrt(max((f.V[c, c], 0))) * sdj > 50) return(1)
    }
    return(0)
}

// ---------------------------------------------------------------------------
// One simulated design matrix (dgp.R: draw_x) -- cross-sectional routes
// ---------------------------------------------------------------------------

struct pape_xdraw scalar pape_draw_x(struct pape_dgp scalar d, real scalar n)
{
    struct pape_xdraw scalar xx
    real matrix U, Z
    real colvector dd, m, idx
    real scalar j, zoff
    xx.has_id = 0; xx.has_start = 0
    if (d.route == "panel") return(pape_draw_x_panel(d, n))
    if (d.route == "iv") return(pape_draw_x_iv(d, n))
    if (d.route == "empirical") {
        idx = ceil(runiform(n, 1) :* d.n_emp)
        if (d.emp_focal_joint) dd = d.emp_xd[idx]
        else dd = pape_transform_u(rnormal(n, 1, 0, 1), d.focal)
        if (d.k > 0) xx.X = (J(n, 1, 1), dd, d.emp_x[idx, .])
        else xx.X = (J(n, 1, 1), dd)
        xx.focal_col = 2
        return(xx)
    }
    U = pape_draw_std_mvn(n, d.Rchol)
    dd = pape_transform_u(U[., 1], d.focal)
    if (d.has_m) m = pape_transform_u(U[., 2], d.moderator)
    zoff = 1 + d.has_m
    if (d.k > 0) {
        Z = U[., (zoff + 1)..(zoff + d.k)]
        for (j = 1; j <= d.k; j++) Z[., j] = pape_transform_u(Z[., j], d.covariates[j])
    }
    else Z = J(n, 0, .)
    if (d.has_m) xx.X = (J(n, 1, 1), dd, m, dd :* m, Z)
    else xx.X = (J(n, 1, 1), dd, Z)
    xx.focal_col = 2
    return(xx)
}

// ---------------------------------------------------------------------------
// Core loop (engine.R: sim_ci)
// ---------------------------------------------------------------------------

struct pape_sim scalar pape_sim_ci(struct pape_dgp scalar d, real scalar n,
    real scalar nsim, real scalar conf, real scalar seed, string scalar se_type,
    string scalar separation)
{
    struct pape_sim scalar s
    struct pape_xdraw scalar xx
    struct pape_fit scalar ft
    struct pape_est scalar est
    struct pape_cffit scalar cf
    real colvector bt, y, pr, b1, b2, start, idv
    real rowvector cc, hd_cols
    real matrix V
    real scalar z, is_aie, is_panel, is_iv, sep_fail, fb, mb, r, s_, ape_v
    real scalar use_stata, use_margins, frame_open, ape_m, se_m, closed
    string scalar st, vce
    z = pape_zcrit(conf)
    is_aie = (d.estimand == "aie")
    is_panel = (d.route == "panel")
    is_iv = (d.route == "iv")
    sep_fail = (separation != "keep")
    fb = (d.focal.type == "binary")
    mb = (d.has_m & d.moderator.type == "binary")
    hd_cols = (d.has_m ? (2, 3, 4) : 2)
    use_stata = (pape_engine() == "stata")
    use_margins = (use_stata & pape_effects() == "margins")
    frame_open = 0
    if (!is_panel & !is_iv) bt = pape_beta_true(d)
    if (seed != .) {
        st = rngstate()
        rseed(seed)
    }
    s.est = J(nsim, 1, .); s.l = J(nsim, 1, .); s.u = J(nsim, 1, .)
    s.se = J(nsim, 1, .); s.min_cell = J(nsim, 1, .)
    s.ok = J(nsim, 1, 0); s.sep = J(nsim, 1, 0)
    for (r = 1; r <= nsim; r++) {
        xx = pape_draw_x(d, n)
        closed = 1
        if (is_iv) {
            y = xx.y
            b1 = (fb ? xx.d : J(0, 1, .))
            b2 = (mb ? xx.m : J(0, 1, .))
            cc = pape_cell_check(y, b1, b2)
            s.min_cell[r] = cc[2]
            if (all(y :== y[1])) continue
            if (cc[1]) {
                s.sep[r] = 1
                if (sep_fail) continue
            }
            if (use_stata) {
                if (!frame_open) {
                    pape_frame_open(rows(y), 1 + d.k, d.has_m + d.n_inst)
                    frame_open = 1
                }
                else pape_frame_resize(rows(y))
                cf = pape_cf_fit_stata(d, xx, !sep_fail, use_margins, ape_m, se_m)
                if (!cf.ok) continue
                if (use_margins) {
                    ape_v = ape_m; s_ = se_m
                    closed = 0
                }
            }
            else {
                cf = pape_cf_fit(y, xx.d, xx.m, d.has_m, xx.Zx, xx.X1, d.first_stage, J(0, 1, .), xx.start)
                if (!cf.ok) continue
            }
            if (closed) {
                if (is_aie) est = pape_cf_aie_est(cf.theta, cf.W, d.focal.type, cf.cf_col, cf.cfm_col)
                else est = pape_cf_ape_est(cf.theta, cf.W, d.focal.type, cf.cf_col, cf.cfm_col, cf.has_m)
                V = cf.V_main
            }
        }
        else {
            if (is_panel) pr = xx.pr
            else pr = pape_G(d.link, xx.X * bt)
            y = runiform(rows(pr), 1) :< pr
            b1 = (fb ? xx.X[., 2] : J(0, 1, .))
            b2 = (mb ? xx.X[., 3] : J(0, 1, .))
            cc = pape_cell_check(y, b1, b2)
            s.min_cell[r] = cc[2]
            if (all(y :== y[1])) continue
            if (cc[1]) {
                s.sep[r] = 1
                if (sep_fail) continue
            }
            idv = (is_panel ? xx.id : J(0, 1, .))
            if (use_stata) {
                if (!frame_open) {
                    pape_frame_open(rows(y), cols(xx.X) - 1, 0)
                    frame_open = 1
                }
                else pape_frame_resize(rows(y))
                vce = (is_panel ? "cluster" : (se_type == "robust" ? "robust" : "oim"))
                if (use_margins) {
                    if (!pape_effects_margins(d, xx.X, y, idv, vce, !sep_fail, ape_m, se_m)) continue
                    ape_v = ape_m; s_ = se_m
                    closed = 0
                }
                else {
                    ft = pape_fit_stata(xx.X, y, idv, d.link, vce, !sep_fail)
                    if (!ft.ok) continue
                }
            }
            else {
                start = (is_panel ? xx.start : bt)
                ft = pape_fit_index(xx.X, y, d.link, start)
                if (!ft.ok) continue
            }
            if (closed) {
                // Hauck-Donner backup for separation the cell check cannot see;
                // a kept asis fit with a zero variance row (ft.rank < k, Stata
                // engine under separation(keep)) counts as separated too
                if (!s.sep[r] & (ft.rank < cols(xx.X) | pape_hauck_donner(ft, xx.X, hd_cols))) {
                    s.sep[r] = 1
                    if (sep_fail) continue
                }
                if (is_aie) est = pape_aie_est(ft.b, xx.X, d.link, d.focal.type, d.moderator.type)
                else est = pape_ape_est(ft.b, xx.X, xx.focal_col, d.link, d.focal.type)
                if (use_stata) V = ft.V
                else if (is_panel) V = pape_vcov_cluster(ft, xx.X, y, xx.id, d.link)
                else if (se_type == "robust") V = pape_vcov_sandwich(ft, xx.X, y, d.link)
                else V = ft.V
            }
        }
        if (closed) {
            s_ = sqrt(max((0, est.jac' * V * est.jac)))
            ape_v = est.ape
        }
        if (s_ == . | s_ <= 0 | ape_v == .) continue
        s.ok[r] = 1
        s.se[r] = s_
        s.est[r] = ape_v
        s.l[r] = ape_v - z * s_
        s.u[r] = ape_v + z * s_
    }
    if (frame_open) pape_frame_close()
    if (seed != .) rngstate(st)
    return(s)
}

// ---------------------------------------------------------------------------
// Coherence guards (engine.R: check_coherence)
// ---------------------------------------------------------------------------

void pape_check_coherence(struct pape_dgp scalar d, string scalar claim,
    real scalar sesoi, real scalar conf)
{
    string scalar lab
    real scalar t_abs
    if (!d.has_effect)
        pape_err("This DGP has no target effect yet; run powerape setape (or setaie) first.")
    lab = strupper(d.estimand == "" ? "ape" : d.estimand)
    t_abs = abs(d.target_est)
    if (claim == "minimum" | claim == "equivalence") {
        if (sesoi == . | sesoi <= 0)
            pape_err(sprintf("claim(%s) needs a single positive sesoi() (in APE units).", claim))
    }
    else if (sesoi != . & sesoi <= 0)
        pape_err("sesoi(), when supplied, must be a single positive number (in APE units).")
    if (claim == "minimum" & t_abs <= sesoi)
        pape_err(sprintf("Assumed true %s (%9.4f in magnitude) does not exceed the SESOI (%9.4f). On or below the null boundary of the minimum-effect test the claim succeeds with probability ~%3.1f%% (the test size) at ANY sample size, so no n delivers meaningful power. Either assume a planning value above the SESOI, or run the classic SESOI analysis instead: claim(detect) with the target set to the SESOI.",
            lab, t_abs, sesoi, 100 * (1 - conf) / 2))
    if (claim == "equivalence" & t_abs >= sesoi)
        pape_err(sprintf("Equivalence power requires |true %s| (%9.4f) strictly below the SESOI (%9.4f); typically the planning value is target = 0. On or beyond the boundary the equivalence claim cannot exceed the test size.",
            lab, t_abs, sesoi))
}

// ---------------------------------------------------------------------------
// Summaries (engine.R: summarize_sim)
// ---------------------------------------------------------------------------

void pape_summarize_sim(struct pape_sim scalar sim, real scalar target,
    string scalar claim, real scalar sesoi, real scalar nsim, real scalar alpha,
    struct pape_res scalar r)
{
    real colvector est, se, l_det, l, u, ok, det, mn, eq, d_only, inc
    real scalar z1
    est = (target < 0 ? -sim.est : sim.est)
    se = sim.se
    l_det = est :- invnormal(1 - alpha / 2) :* se
    z1 = invnormal(1 - alpha)
    l = est :- z1 :* se
    u = est :+ z1 :* se
    ok = sim.ok
    det = ok :& (l_det :> 0)
    if (sesoi != .) {
        mn = ok :& (l :> sesoi)
        eq = ok :& (l :> -sesoi) :& (u :< sesoi)
        d_only = det :& !mn :& !eq
        inc = ok :& !mn :& !eq :& !d_only
        r.outcomes = (mean(mn), mean(d_only), mean(inc), mean(eq), mean(!ok))
        r.outnames = ("minimum", "detect_only", "inconclusive", "equivalence", "failed")
        if (claim == "minimum") r.power = mean(mn)
        else if (claim == "detect") r.power = mean(det)
        else r.power = mean(eq)
        r.has_sesoi = 1
    }
    else {
        r.outcomes = (mean(det), mean(ok :& !det), mean(!ok))
        r.outnames = ("detect", "inconclusive", "failed")
        r.power = mean(det)
        r.has_sesoi = 0
    }
    r.mcse = sqrt(r.power * (1 - r.power) / nsim)
    r.n_failed = sum(!ok)
    r.separated = mean(sim.sep)
    r.n_separated = sum(sim.sep)
    if (all(sim.min_cell :== .)) r.min_cell = .
    else r.min_cell = mean(select(sim.min_cell, sim.min_cell :!= .))
}

// ---------------------------------------------------------------------------
// Warnings for separation / sparse cells (engine.R: sparse_warning)
// ---------------------------------------------------------------------------

real scalar pape_sparse_warning(real scalar separated, real scalar min_cell,
    string scalar separation, string scalar where, real scalar aie)
{
    string scalar cell, msg
    cell = (aie ? "focal-by-moderator cell" : "focal cell")
    msg = ""
    if (separated != . & separated > 0.01)
        msg = sprintf("%3.1f%% of simulated studies had separation (a %s with no events or no non-events, where the maximum-likelihood effect does not exist)%s",
            100 * separated, cell,
            (separation == "keep" ? "; they were kept with their degenerate Wald intervals (separation(keep))" :
             "; they count as failed (separation(fail), the default; Stata's probit also refuses such fits)"))
    if (min_cell != . & min_cell < 10)
        msg = msg + (msg != "" ? "; " : "") +
            sprintf("the smallest %s averages %3.1f events (or non-events) per simulated study; delta-method Wald inference is unreliable at this sparsity (one-sided error rates can run up to about twice nominal, in the direction in which the sparse cell has the lower rate)",
                cell, min_cell)
    if (msg != "") {
        printf("{err}warning: %s%s. See help powerape, section 'Separation and sparse cells'.\n",
            (where != "" ? where + ": " : ""), msg)
        return(1)
    }
    return(0)
}

// ---------------------------------------------------------------------------
// Shared power computation (engine.R: power_once)
// ---------------------------------------------------------------------------

struct pape_res scalar pape_power_once(struct pape_dgp scalar d, real scalar n,
    string scalar claim, real scalar sesoi, real scalar conf, real scalar nsim,
    real scalar seed, real scalar enforce, string scalar se, real scalar alpha,
    string scalar rule, real scalar keep_draws, string scalar separation)
{
    struct pape_res scalar r
    struct pape_sim scalar sim
    real rowvector lv
    string scalar se_use
    if (n == . | n < 20) pape_err("n() must be at least 20.")
    if (d.route == "panel" & n < 30)
        printf("{err}warning: fewer than 30 units (clusters): cluster-robust inference is unreliable at this size.\n")
    lv = pape_claim_levels(claim, alpha, conf, n, rule)
    if (nsim == . | nsim < 20) pape_err("nsim() must be at least 20.")
    se_use = se
    if (se != "model" & (d.route == "panel" | d.route == "iv")) {
        printf("{err}warning: se() is fixed by the route: panel designs use unit-clustered SEs, IV designs the stacked robust sandwich; se(%s) is ignored.\n", se)
        se_use = "model"
    }
    if (enforce) pape_check_coherence(d, claim, sesoi, lv[4])
    sim = pape_sim_ci(d, n, nsim, lv[4], seed, se_use, separation)
    pape_summarize_sim(sim, d.target_est, claim, sesoi, nsim, lv[1], r)
    r.kind = "power"
    r.claim = claim; r.sesoi = sesoi; r.conf = lv[4]; r.alpha = lv[1]
    r.alpha_rule = rule; r.has_rule = (rule != "")
    r.n = n; r.nsim = nsim; r.target = d.target_est
    r.estimand = (d.estimand == "" ? "ape" : d.estimand)
    r.se = se_use; r.separation = separation; r.model = d.model
    r.dgp = d
    pape_stamp_fit(r)
    if (keep_draws) {
        r.draws = sim
        r.has_draws = 1
    }
    else r.has_draws = 0
    r.direction = ""
    r.goal = .; r.confirmed = 0; r.confirm_requested = 0; r.nsim_confirm = .
    r.mde = .; r.is_margin = 0
    r.history = J(0, 0, .); r.results = J(0, 0, .)
    return(r)
}

// ---------------------------------------------------------------------------
// Printing (print.R)
// ---------------------------------------------------------------------------

void pape_print_sep_note(struct pape_res scalar x)
{
    real scalar sepd, keep, other, mc
    string scalar cell
    sepd = (x.separated == . ? 0 : x.separated)
    keep = (x.separation == "keep")
    cell = (x.estimand == "aie" ? "focal-by-moderator cell" : "focal cell")
    if (x.kind == "power") {
        other = x.outcomes[cols(x.outcomes)] - (keep ? 0 : sepd)
        if (other > 1e-12)
            printf("  note: %3.1f%% of replications failed to fit; they count against power.\n", 100 * other)
    }
    if (sepd > 0)
        printf("  note: %3.1f%% of replications had separation (a %s with no events or no non-events)%s.\n",
            100 * sepd, cell,
            (keep ? "; kept with their degenerate Wald intervals (separation(keep))" : "; they count as failed (separation(fail))"))
    mc = x.min_cell
    if (mc != . & mc < 10)
        printf("  note: the smallest %s averages %3.1f events (or non-events) per study; Wald inference is unreliable at this sparsity.\n", cell, mc)
}

string scalar pape_claim_lab(string scalar claim, real scalar sesoi)
{
    if (claim == "minimum") return(sprintf("minimum-effect claim (CI lower bound > %5.3f)", sesoi))
    if (claim == "detect") return("detection claim (CI excludes 0, directional)")
    return(sprintf("equivalence claim (CI within +/-%5.3f)", sesoi))
}

string scalar pape_n_lab(struct pape_res scalar x)
{
    if (x.dgp.route == "panel")
        return(sprintf("n = %f units x %s periods (%s obs)", x.n, pape_panel_len_label(x.dgp), pape_panel_obs_label(x.dgp, x.n)))
    return(sprintf("n = %f", x.n))
}

string scalar pape_se_lab(struct pape_res scalar x)
{
    if (x.dgp.route == "iv") return(", stacked robust SEs")
    if (x.se == "robust") return(", robust SEs")
    return("")
}

void pape_print_power(struct pape_res scalar x)
{
    string scalar al_lab, o
    real scalar j
    printf("{txt}powerape -- %s\n", pape_claim_lab(x.claim, x.sesoi))
    if (x.has_rule) al_lab = sprintf("alpha = %6.4g by sample-size rule", x.alpha)
    else al_lab = sprintf("alpha = %g", x.alpha)
    printf("  %s, %s, assumed true %s %s%7.4f, %6.4g%% CI (%s), nsim = %f%s\n",
        x.model, pape_n_lab(x), strupper(x.estimand), pape_sgn(x.target), x.target, 100 * x.conf, al_lab, x.nsim, pape_se_lab(x))
    printf("  fit: %s\n", x.fit)
    printf("  power = %5.3f (MCSE %5.3f)\n", x.power, x.mcse)
    o = ""
    for (j = 1; j <= cols(x.outcomes); j++)
        o = o + (j > 1 ? " | " : "") + sprintf("%s %5.3f", subinstr(x.outnames[j], "_", "-", .), x.outcomes[j])
    printf("  outcomes: %s\n", o)
    pape_print_sep_note(x)
}

end
