*! powerape_panel.mata -- panel designs: correlated random effects probit
*! (Mundlak means, unit-clustered SEs, unbalanced panels). Port of R/panel.R
*! of powerape 1.11.3 (DESIGN.md sections 11 and 13). Probit-only.
version 16
mata:
mata set matastrict on

// Draw panel lengths from the support/probability vectors.
real colvector pape_draw_T(real scalar n_units, real rowvector t_support,
    real rowvector t_prob)
{
    real colvector u, T, idx
    real rowvector cp
    real scalar j
    cp = runningsum(t_prob)
    u = runiform(n_units, 1)
    T = J(n_units, 1, t_support[cols(t_support)])
    for (j = cols(t_support) - 1; j >= 1; j--) {
        idx = selectindex(u :< cp[j])
        if (rows(idx) > 0 & cols(idx) > 0) T[idx] = J(rows(idx), 1, t_support[j])
    }
    return(T)
}

// Unit index of a unit-major panel (R: rep(seq_len(n_units), times = T_i)):
// a 1 on each unit's first row, then a running sum (every T_i >= 1). Replaces
// an append loop that was quadratic in the number of units (0.1.3); same
// integers and no random draws, so every result is bit-identical.
real colvector pape_unit_index(real colvector T_i)
{
    real colvector first
    if (any(T_i :< 1)) pape_err("Internal error: every unit needs at least one period.")
    first = J(sum(T_i), 1, 0)
    first[runningsum(T_i) :- T_i :+ 1] = J(rows(T_i), 1, 1)
    return(runningsum(first))
}

struct pape_pdraw scalar pape_panel_reg_draw(struct pape_var colvector vars,
    real colvector iccs, real matrix R, real scalar n_units, real scalar T_p,
    real scalar seed, real rowvector t_support, real rowvector t_prob)
{
    struct pape_pdraw scalar rd
    real matrix Rc, V, E, info
    real colvector u, xj
    real scalar k, N, j, unbal
    string scalar st
    k = rows(vars)
    if (seed != .) {
        st = rngstate()
        rseed(seed)
    }
    unbal = (cols(t_support) > 1)
    if (unbal) rd.T_i = pape_draw_T(n_units, t_support, t_prob)
    else rd.T_i = J(n_units, 1, T_p)
    N = sum(rd.T_i)
    Rc = cholesky(R)'
    V = rnormal(n_units, k, 0, 1) * Rc
    E = rnormal(N, k, 0, 1) * Rc
    rd.unit = pape_unit_index(rd.T_i)
    info = panelsetup(rd.unit, 1)
    rd.x = J(N, k, .)
    rd.xbar = J(n_units, k, .)
    for (j = 1; j <= k; j++) {
        u = sqrt(iccs[j]) :* V[rd.unit, j] :+ sqrt(1 - iccs[j]) :* E[., j]
        xj = pape_transform_u(u, vars[j])
        rd.x[., j] = xj
        rd.xbar[., j] = panelsum(xj, info) :/ rd.T_i
    }
    if (seed != .) rngstate(st)
    return(rd)
}

struct pape_var colvector pape_panel_vars(struct pape_dgp scalar d)
{
    struct pape_var colvector vars
    real scalar j, off
    vars = pape_var(1 + d.has_m + d.k, 1)
    vars[1] = d.focal
    if (d.has_m) vars[2] = d.moderator
    off = 1 + d.has_m
    for (j = 1; j <= d.k; j++) vars[off + j] = d.covariates[j]
    return(vars)
}

struct pape_idraw scalar pape_panel_integration(struct pape_dgp scalar d)
{
    struct pape_idraw scalar id
    struct pape_pdraw scalar rd
    real colvector idxz, m_i
    real matrix Z, Mbar
    real scalar off
    rd = pape_panel_reg_draw(pape_panel_vars(d), d.iccs, d.R, d.n_int, d.n_periods,
        d.seed_int, d.t_support, d.t_prob)
    off = 1 + d.has_m
    if (d.k > 0) {
        Z = rd.x[., (off + 1)..(off + d.k)]
        idxz = Z * d.gamma
    }
    else idxz = J(rows(rd.unit), 1, 0)
    if (rows(d.tv_idx) > 0 & any(d.xi :!= 0)) {
        Mbar = rd.xbar[., d.tv_idx']
        m_i = (Mbar :- d.xbar_center) * d.xi
    }
    else m_i = J(d.n_int, 1, 0)
    id.xd_c = rd.x[., 1] :- d.focal_ref
    id.has_xd = 1
    if (d.has_m) {
        id.md_c = rd.x[., 2] :- d.mod_ref
        id.has_md = 1
    }
    else id.has_md = 0
    id.q = idxz :+ m_i[rd.unit]
    id.idxz = id.q
    return(id)
}

string scalar pape_panel_len_label(struct pape_dgp scalar d)
{
    if (d.unbalanced)
        return(sprintf("%f-%f (unbalanced, mean %4.2f)", min(d.t_support), max(d.t_support), d.expected_T))
    return(sprintf("%f", d.n_periods))
}

string scalar pape_panel_obs_label(struct pape_dgp scalar d, real scalar n)
{
    if (d.unbalanced) return(sprintf("~%f", round(n * d.expected_T)))
    return(sprintf("%f", n * d.n_periods))
}

// E_a[Phi(x + a)], a ~ N(0, var_a): closed form (probit)
real matrix pape_mixP(struct pape_dgp scalar d, real matrix x) return(normal(x :/ d.scale_a))
real matrix pape_mixp(struct pape_dgp scalar d, real matrix x) return(normalden(x :/ d.scale_a) :/ d.scale_a)

real scalar pape_f_b0_panel(real scalar b0, struct pape_env scalar e)
    return(mean(pape_mixP(e.dgp, b0 :+ e.idxz)) - e.target)

// ape_dgp_panel
struct pape_dgp scalar pape_dgp_panel(string scalar name, string scalar model,
    struct pape_var colvector vars, real scalar has_m,
    real rowvector n_periods, real rowvector p_periods, real scalar retention,
    real scalar rho, real scalar cre_share,
    real scalar corr_given, real scalar corr_scalar, real matrix corr_matrix,
    real scalar baseline, real scalar signal, real scalar n_int, real scalar seed_int)
{
    struct pape_dgp scalar d
    struct pape_env scalar e
    struct pape_pdraw scalar rd
    real rowvector t_support, t_prob, T_in, o
    real colvector g_unit, idxz, m_unit, m_i, w_unit, sds, keep
    real matrix Z, Mbar
    real scalar k, j, v_unit, var_c, T_p, off, unbalanced, expected_T
    d = pape_dgp_blank()
    d.name = name
    if (model == "logit")
        pape_err("Panel designs are probit-only. The device that lets a pooled CRE estimator recover the average structural function is specific to the probit link: a normal unit effect convolved with a normal error just rescales the index by sqrt(1 + var_a). The logistic link has no such stability, so a pooled CRE logit estimates only a quasi-ML approximation of the ASF, and powerape does not price designs with an approximate estimand. Use model(probit).")
    if (model != "probit") pape_err("model() must be probit for panel DGPs.")
    d.model = "probit"; d.link = "probit"
    d.route = "panel"; d.builder = "panel"
    d.focal = vars[1]
    d.has_m = has_m
    if (has_m) {
        d.moderator = vars[2]
        if (d.moderator.name == d.focal.name) pape_err("moderator() must be a different variable than focal().")
        if (d.focal.type != "binary" | d.moderator.type != "binary")
            pape_err("Panel AIE designs currently support a binary focal x binary moderator; continuous pairs are on the roadmap.")
    }
    k = rows(vars) - 1 - has_m
    d.k = k
    if (k > 0) d.covariates = vars[(2 + has_m)..rows(vars)]
    // panel-length distribution (DESIGN.md section 13)
    if (cols(p_periods) > 0 & retention != .) pape_err("Give either pperiods() or retention(), not both.")
    T_in = n_periods
    if (cols(T_in) < 1 | any(T_in :< 1) | any(T_in :!= floor(T_in))) pape_err("nperiods() must be positive integers.")
    if (rows(uniqrows(T_in')) != cols(T_in)) pape_err("nperiods() must not contain duplicates.")
    if (retention != .) {
        if (cols(T_in) != 1 | retention <= 0 | retention > 1)
            pape_err("retention() needs a scalar nperiods() and a rate in (0, 1].")
        if (retention < 1) {
            t_support = 1..T_in
            t_prob = retention:^(t_support :- 1) :* (1 - retention)
            t_prob[T_in] = retention^(T_in - 1)
        }
        else {
            t_support = T_in
            t_prob = 1
        }
    }
    else if (cols(T_in) > 1) {
        o = order(T_in', 1)'
        t_support = T_in[o]
        if (cols(p_periods) == 0) t_prob = J(1, cols(T_in), 1 / cols(T_in))
        else {
            if (cols(p_periods) != cols(T_in) | any(p_periods :< 0) | sum(p_periods) <= 0)
                pape_err("pperiods() must match nperiods() in length, be nonnegative, and not all zero.")
            t_prob = p_periods[o] :/ sum(p_periods)
        }
        keep = selectindex(t_prob :> 0)
        t_support = t_support[keep]
        t_prob = t_prob[keep]
    }
    else {
        if (cols(p_periods) > 0) pape_err("pperiods() needs a vector nperiods() of the same length.")
        t_support = T_in
        t_prob = 1
    }
    if (max(t_support) < 2) pape_err("Panel designs need at least two periods for some units.")
    unbalanced = (cols(t_support) > 1)
    T_p = (unbalanced ? . : t_support[1])
    expected_T = sum(t_support :* t_prob)
    if (rho == . | rho < 0 | rho >= 1) pape_err("rho() must lie in [0, 1).")
    if (cre_share == . | cre_share < 0 | cre_share > 1) pape_err("creshare() must lie in [0, 1].")
    if (baseline == . | baseline <= 0 | baseline >= 1) pape_err("baseline() must lie strictly between 0 and 1.")
    if (signal == . | signal < 0 | signal >= 1) pape_err("signal() must lie in [0, 1).")
    if (n_int < 2e4) pape_err("nint() must be at least 20000 for panel DGPs.")
    d.iccs = J(rows(vars), 1, 0)
    for (j = 1; j <= rows(vars); j++) d.iccs[j] = (vars[j].icc == . ? 0 : vars[j].icc)
    d.tv_idx = selectindex(d.iccs :< 1)
    if (rows(d.tv_idx) == 0 & cols(d.tv_idx) == 0) d.tv_idx = J(0, 1, .)
    d.tv_idx = vec(d.tv_idx)
    if (k == 0 & signal > 0) pape_err("signal() > 0 requires at least one covariate.")
    if (cre_share > 0 & rho > 0 & rows(d.tv_idx) == 0)
        pape_err("creshare() > 0 needs at least one time-varying regressor (icc < 1) for the heterogeneity to load on.")
    d.R = pape_build_corr(k + 1 + has_m, corr_given, corr_scalar, corr_matrix)
    d.Rchol = cholesky(d.R)'
    d.corr_given = corr_given; d.corr_scalar = corr_scalar; d.corr_matrix = corr_matrix
    d.n_periods = T_p
    d.rho = rho; d.cre_share = cre_share
    d.unbalanced = unbalanced
    d.t_support = (unbalanced ? t_support : J(1, 0, .))
    d.t_prob = (unbalanced ? t_prob : J(1, 0, .))
    d.expected_T = expected_T
    d.retention = retention
    d.n_int = n_int; d.seed_int = seed_int
    d.baseline = baseline; d.signal = signal
    d.focal_ref = (d.focal.type == "binary" ? 0 : d.focal.mean)
    d.mod_ref = (has_m ? 0 : .)

    rd = pape_panel_reg_draw(vars, d.iccs, d.R, n_int, T_p, seed_int, d.t_support, d.t_prob)
    off = 1 + has_m
    if (k > 0) {
        Z = rd.x[., (off + 1)..(off + k)]
        if (signal == 0) d.gamma = J(k, 1, 0)
        else {
            g_unit = J(k, 1, .)
            for (j = 1; j <= k; j++) g_unit[j] = 1 / pape_var_sd(d.covariates[j])
            v_unit = variance(Z * g_unit)
            d.gamma = sqrt((signal / (1 - signal)) / v_unit) :* g_unit
        }
        idxz = Z * d.gamma
    }
    else {
        d.gamma = J(0, 1, .)
        idxz = J(rows(rd.unit), 1, 0)
    }
    d.s2 = (k > 0 ? variance(idxz) : 0)
    d.signal_implied = d.s2 / (1 + d.s2)
    var_c = (rho > 0 ? rho / (1 - rho) : 0)
    if (var_c > 0 & cre_share > 0) {
        Mbar = rd.xbar[., d.tv_idx']
        d.xbar_center = mean(Mbar)
        sds = sqrt(diagonal(variance(Mbar)))
        if (any(sds :< 1e-10)) pape_err("A time-varying regressor has (numerically) constant unit means.")
        w_unit = 1 :/ sds
        m_unit = (Mbar :- d.xbar_center) * w_unit
        d.xi = sqrt(cre_share * var_c / variance(m_unit)) :* w_unit
        m_i = (Mbar :- d.xbar_center) * d.xi
    }
    else {
        if (rows(d.tv_idx) > 0) d.xbar_center = mean(rd.xbar[., d.tv_idx'])
        else d.xbar_center = J(1, 0, .)
        d.xi = J(rows(d.tv_idx), 1, 0)
        m_i = J(n_int, 1, 0)
    }
    d.var_a = (1 - cre_share) * var_c
    d.scale_a = sqrt(1 + d.var_a)
    e.dgp = d
    e.idxz = idxz :+ m_i[rd.unit]
    e.target = baseline
    d.beta0 = pape_uniroot(&pape_f_b0_panel(), e, -15 * d.scale_a, 15 * d.scale_a, 1, 1e-10)
    return(d)
}

// ---- inversion and truth (panel) -------------------------------------------

real scalar pape_f_bshift_panel(real scalar b, struct pape_env scalar e)
    return(mean(pape_mixP(e.dgp, e.b0 :+ b :+ e.idxz)) - e.p_ref - e.target)

real scalar pape_solve_binary_panel(struct pape_dgp scalar d, real scalar target,
    real colvector q, string scalar label, real scalar root, string scalar msg)
{
    struct pape_env scalar e
    real scalar p0, eps
    p0 = mean(pape_mixP(d, d.beta0 :+ q))
    eps = 1e-4
    if (target <= -p0 + eps | target >= (1 - p0) - eps) {
        msg = sprintf("Target %s %9.4f is infeasible for this panel DGP: with baseline P(Y=1 | reference) = %5.3f the attainable range is (%6.3f, %5.3f).",
            label, target, p0, -p0, 1 - p0)
        return(1)
    }
    e.dgp = d; e.b0 = d.beta0; e.idxz = q; e.p_ref = p0; e.target = target
    root = pape_uniroot(&pape_f_bshift_panel(), e, -8 * d.scale_a, 8 * d.scale_a, 1, 1e-9)
    msg = ""
    return(0)
}

real scalar pape_f_cont_panel(real scalar b, struct pape_env scalar e)
    return(b * mean(pape_mixp(e.dgp, e.b0 :+ b :* e.x_c :+ e.idxz)) - e.target)

real scalar pape_solve_cont_panel(struct pape_dgp scalar d, real scalar target,
    real colvector xd_c, real colvector q, real scalar root, string scalar msg)
{
    struct pape_env scalar e
    real colvector bgrid, vals
    real scalar s, sd_d, j, amax, i, lo, hi, tmp
    s = d.scale_a
    e.dgp = d; e.b0 = d.beta0; e.idxz = q; e.x_c = xd_c; e.target = 0
    sd_d = pape_var_sd(d.focal)
    bgrid = (target < 0 ? -1 : 1) :* exp(rangen(ln(1e-3 / sd_d), ln(200 * s / sd_d), 100))
    vals = J(100, 1, .)
    for (j = 1; j <= 100; j++) vals[j] = pape_f_cont_panel(bgrid[j], e)
    amax = max(abs(vals))
    if (abs(target) >= amax * (1 - 1e-3)) {
        msg = sprintf("Target APE %9.4f is infeasible for this panel DGP: the attainable range is about (%9.4f, %9.4f). The average derivative plateaus as the coefficient grows.",
            target, -amax, amax)
        return(1)
    }
    i = min(selectindex(abs(vals) :>= abs(target)))
    lo = (i == 1 ? 0 : bgrid[i - 1]); hi = bgrid[i]
    if (lo > hi) {
        tmp = lo; lo = hi; hi = tmp
    }
    e.target = target
    root = pape_uniroot(&pape_f_cont_panel(), e, lo, hi, 0, 1e-9)
    msg = ""
    return(0)
}

real scalar pape_set_ape_panel(struct pape_dgp scalar d, real scalar target, string scalar msg)
{
    struct pape_idraw scalar id
    real scalar root, rc
    id = pape_panel_integration(d)
    if (d.focal.type == "binary") rc = pape_solve_binary_panel(d, target, id.q, "APE", root, msg)
    else rc = pape_solve_cont_panel(d, target, id.xd_c, id.q, root, msg)
    if (rc) return(rc)
    d.p1 = (d.focal.type == "binary" ? mean(pape_mixP(d, d.beta0 :+ root :+ id.q)) : .)
    d.beta_focal = root
    d.estimand = "ape"
    d.target_est = target
    d.has_effect = 1
    return(0)
}

real scalar pape_true_ape_panel(struct pape_dgp scalar d)
{
    struct pape_idraw scalar id
    id = pape_panel_integration(d)
    if (d.focal.type == "binary")
        return(mean(pape_mixP(d, d.beta0 :+ d.beta_focal :+ id.q) :- pape_mixP(d, d.beta0 :+ id.q)))
    return(d.beta_focal * mean(pape_mixp(d, d.beta0 :+ d.beta_focal :* id.xd_c :+ id.q)))
}

real scalar pape_f_p11_panel(real scalar b, struct pape_env scalar e)
    return(mean(pape_mixP(e.dgp, e.b0 :+ e.bf :+ e.bm :+ b :+ e.idxz)) - e.target)

real scalar pape_set_aie_panel(struct pape_dgp scalar d, real scalar target,
    real scalar main_focal, real scalar main_moderator, string scalar msg)
{
    struct pape_idraw scalar id
    struct pape_env scalar e
    real scalar p00, bf, bm, bi, rc, p11_target, eps
    id = pape_panel_integration(d)
    p00 = mean(pape_mixP(d, d.beta0 :+ id.q))
    rc = pape_solve_binary_panel(d, main_focal, id.q, "main-effect APE (focal)", bf, msg)
    if (rc) return(rc)
    rc = pape_solve_binary_panel(d, main_moderator, id.q, "main-effect APE (moderator)", bm, msg)
    if (rc) return(rc)
    p11_target = p00 + main_focal + main_moderator + target
    eps = 1e-4
    if (p11_target <= eps | p11_target >= 1 - eps) {
        msg = sprintf("Target AIE %9.4f is infeasible for this panel DGP: it implies the counterfactual cell rate P(Y=1 | focal=1, moderator=1) = %9.4f, outside (0,1). With baseline %5.3f and main-effect anchors %9.4f and %9.4f, the attainable AIE range is (%9.4f, %9.4f).",
            target, p11_target, p00, main_focal, main_moderator,
            eps - p00 - main_focal - main_moderator, 1 - eps - p00 - main_focal - main_moderator)
        return(1)
    }
    e.dgp = d; e.b0 = d.beta0; e.bf = bf; e.bm = bm; e.idxz = id.q; e.target = p11_target
    bi = pape_uniroot(&pape_f_p11_panel(), e, -8 * d.scale_a, 8 * d.scale_a, 1, 1e-9)
    d.beta_focal = bf; d.beta_mod = bm; d.beta_int = bi
    d.main_focal = main_focal; d.main_moderator = main_moderator
    d.estimand = "aie"
    d.target_est = target
    d.p1 = .
    d.has_effect = 1
    msg = ""
    return(0)
}

real scalar pape_true_aie_panel(struct pape_dgp scalar d)
{
    struct pape_idraw scalar id
    real scalar b0, bf, bm, bi
    id = pape_panel_integration(d)
    b0 = d.beta0; bf = d.beta_focal; bm = d.beta_mod; bi = d.beta_int
    return(mean(pape_mixP(d, b0 + bf + bm + bi :+ id.q)) - mean(pape_mixP(d, b0 + bm :+ id.q)) -
           mean(pape_mixP(d, b0 + bf :+ id.q)) + mean(pape_mixP(d, b0 :+ id.q)))
}

// ---- engine draw (panel) ------------------------------------------------------

struct pape_xdraw scalar pape_draw_x_panel(struct pape_dgp scalar d, real scalar n_units)
{
    struct pape_xdraw scalar xx
    struct pape_pdraw scalar rd
    real colvector dd, m, dm, idxz, m_i, a_i, eta, d_c, lev, dmbar, T_i, start
    real matrix Z, Mbar_cols, Tdum, info
    real scalar off, n_dum, fr, mr, s, b0_raw, j, has_dmbar
    if (!d.has_effect) pape_err("This DGP has no target effect yet; run powerape setape (or setaie) first.")
    rd = pape_panel_reg_draw(pape_panel_vars(d), d.iccs, d.R, n_units, d.n_periods, .,
        d.t_support, d.t_prob)
    off = 1 + d.has_m
    dd = rd.x[., 1]
    if (d.has_m) m = rd.x[., 2]
    if (d.k > 0) Z = rd.x[., (off + 1)..(off + d.k)]
    else Z = J(rows(dd), 0, .)
    if (rows(d.tv_idx) > 0) Mbar_cols = rd.xbar[rd.unit, d.tv_idx']
    else Mbar_cols = J(rows(dd), 0, .)
    Tdum = J(rows(dd), 0, .)
    if (d.unbalanced) {
        lev = uniqrows(rd.T_i)
        if (rows(lev) > 1) {
            Tdum = J(rows(dd), rows(lev) - 1, .)
            for (j = 2; j <= rows(lev); j++) Tdum[., j - 1] = (rd.T_i[rd.unit] :== lev[j])
        }
    }
    n_dum = cols(Tdum)
    has_dmbar = 0
    if (d.has_m) {
        dm = dd :* m
        if (d.iccs[1] < 1 | d.iccs[2] < 1) {
            info = panelsetup(rd.unit, 1)
            dmbar = (panelsum(dm, info) :/ rd.T_i)[rd.unit]
            has_dmbar = 1
            xx.X = (J(rows(dd), 1, 1), dd, m, dm, Z, Mbar_cols, dmbar, Tdum)
        }
        else xx.X = (J(rows(dd), 1, 1), dd, m, dm, Z, Mbar_cols, Tdum)
    }
    else xx.X = (J(rows(dd), 1, 1), dd, Z, Mbar_cols, Tdum)
    // structural probabilities
    if (d.k > 0) idxz = Z * d.gamma
    else idxz = J(rows(dd), 1, 0)
    if (rows(d.tv_idx) > 0 & any(d.xi :!= 0)) m_i = (rd.xbar[., d.tv_idx'] :- d.xbar_center) * d.xi
    else m_i = J(n_units, 1, 0)
    if (d.var_a > 0) a_i = rnormal(n_units, 1, 0, sqrt(d.var_a))
    else a_i = J(n_units, 1, 0)
    fr = d.focal_ref
    mr = (d.has_m ? d.mod_ref : 0)
    d_c = dd :- fr
    eta = d.beta0 :+ d.beta_focal :* d_c :+ idxz :+ (m_i :+ a_i)[rd.unit]
    if (d.has_m) eta = eta :+ d.beta_mod :* (m :- mr) :+ d.beta_int :* d_c :* (m :- mr)
    xx.pr = normal(eta)
    // pseudo-true (attenuated) start values on the raw columns
    s = d.scale_a
    b0_raw = d.beta0 - d.beta_focal * fr
    if (d.has_m) b0_raw = b0_raw - d.beta_mod * mr + d.beta_int * fr * mr
    start = (b0_raw - sum(d.xi :* d.xbar_center'))
    if (d.has_m) start = start \ (d.beta_focal - d.beta_int * mr) \ (d.beta_mod - d.beta_int * fr) \ d.beta_int
    else start = start \ d.beta_focal
    start = start \ d.gamma \ d.xi
    if (has_dmbar) start = start \ 0
    start = start :/ s
    if (n_dum > 0) start = start \ J(n_dum, 1, 0)
    xx.start = start
    xx.has_start = 1
    xx.id = rd.unit
    xx.has_id = 1
    xx.focal_col = 2
    return(xx)
}

end
