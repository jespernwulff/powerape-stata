*! powerape_iv.mata -- IV designs: endogenous focal, control-function probit
*! (the estimator of Stata's cfprobit, with stacked method-of-moments SEs).
*! Port of R/iv.R and R/cf.R of powerape 1.11.3 (DESIGN.md section 12).
version 16
mata:
mata set matastrict on

// Probit score factor s(q, y) and its derivative in q (cf.R: probit_score_bits)
void pape_score_bits(real colvector q, real colvector y, real colvector s, real colvector ds)
{
    real colvector g, G, B, A
    g = normalden(q)
    G = normal(q)
    B = rowmax((G :* (1 :- G), J(rows(q), 1, 1e-12)))
    A = y :- G
    s = g :* A :/ B
    ds = (-q :* g :* A :- g:^2) :/ B :- g :* A :* g :* (1 :- 2 :* G) :/ B:^2
}

// Poisson IRLS fit (first stage only). Returns the coefficient vector; ok flag via rc.
real colvector pape_fit_poisson(real matrix X, real colvector y, real scalar ok)
{
    real colvector eta, mu, z, b, w
    real matrix V
    real scalar dev, devold, it, k, converged
    k = cols(X)
    eta = ln(y :+ 0.1)
    devold = .
    converged = 0
    b = J(k, 1, .)
    for (it = 1; it <= 100; it++) {
        mu = exp(eta)
        w = mu
        z = eta :+ (y :- mu) :/ mu
        V = invsym(cross(X, w, X))
        if (any(diagonal(V) :== 0)) {
            ok = 0
            return(b)
        }
        b = V * cross(X, w, z)
        eta = X * b
        mu = exp(eta)
        dev = 2 * sum((y :> 0) :* y :* ln(rowmax((y, J(rows(y), 1, 1e-300))) :/ mu) :- (y :- mu))
        if (devold != . & abs(dev - devold) / (abs(dev) + 0.1) < 1e-8) {
            converged = 1
            break
        }
        devold = dev
    }
    ok = (converged & !hasmissing(b))
    return(b)
}

// First-stage fit for one endogenous variable (cf.R: cf_first)
real scalar pape_cf_first(real colvector y2, real matrix X1, string scalar type,
    real colvector pi, real colvector q, real colvector c, real colvector dcdq)
{
    struct pape_fit scalar ft
    real matrix V
    real colvector s, ds
    real scalar ok
    if (type == "linear") {
        V = invsym(cross(X1, X1))
        if (any(diagonal(V) :== 0)) return(0)
        pi = V * cross(X1, y2)
        q = X1 * pi
        c = y2 :- q
        dcdq = J(rows(q), 1, -1)
        return(1)
    }
    if (type == "poisson") {
        pi = pape_fit_poisson(X1, y2, ok)
        if (!ok) return(0)
        q = X1 * pi
        c = y2 :- exp(q)
        dcdq = -exp(q)
        return(1)
    }
    ft = pape_fit_index(X1, y2, "probit", J(0, 1, .))
    if (!ft.converged | ft.rank != cols(X1) | hasmissing(ft.b)) return(0)
    pi = ft.b
    q = X1 * pi
    pape_score_bits(q, y2, s, ds)
    c = s
    dcdq = ds
    return(1)
}

// Two-step control-function probit with stacked no-bootstrap vcov (cf.R: cf_fit)
struct pape_cffit scalar pape_cf_fit(real colvector y0, real colvector y2,
    real colvector m, real scalar has_m, real matrix Z, real matrix X1,
    string scalar first, real colvector id, real colvector start)
{
    struct pape_cffit scalar f
    struct pape_fit scalar ft
    real colvector pi, q, c, dcdq, idx, s, ds, didx_dc, ord
    real matrix W, A_mm, A_m1, A_11, A, S, B, Ainv, V, info
    real scalar k_m, k_1, n
    f.ok = 0
    n = rows(y0)
    if (!pape_cf_first(y2, X1, first, pi, q, c, dcdq)) return(f)
    if (has_m) W = (J(n, 1, 1), y2, m, y2 :* m, Z, c, c :* m)
    else W = (J(n, 1, 1), y2, Z, c)
    f.cf_col = cols(W) - (has_m ? 1 : 0)
    f.cfm_col = (has_m ? cols(W) : .)
    f.has_m = has_m
    ft = pape_fit_index(W, y0, "probit", start)
    if (!ft.ok) return(f)
    f.theta = ft.b
    idx = W * f.theta
    pape_score_bits(idx, y0, s, ds)
    k_m = cols(W); k_1 = cols(X1)
    A_mm = cross(W, ds, W)
    if (has_m) didx_dc = f.theta[f.cf_col] :+ f.theta[f.cfm_col] :* m
    else didx_dc = J(n, 1, f.theta[f.cf_col])
    A_m1 = cross((ds :* didx_dc) :* W, dcdq :* X1)
    A_11 = cross(X1, dcdq, X1)
    A = (A_mm, A_m1) \ (J(k_1, k_m, 0), A_11)
    S = (s :* W, c :* X1)
    if (rows(id) > 0) {
        ord = order(id, 1)
        info = panelsetup(id[ord], 1)
        B = cross(panelsum(S[ord, .], info), panelsum(S[ord, .], info))
    }
    else B = cross(S, S)
    Ainv = luinv(A)
    if (hasmissing(Ainv)) return(f)
    V = Ainv * B * Ainv'
    f.V = V
    f.V_main = V[1..k_m, 1..k_m]
    f.W = W
    f.cf = c
    f.pi = pi
    f.ok = 1
    return(f)
}

// ASF-based APE of the endogenous focal with the control function held fixed
struct pape_est scalar pape_cf_ape_est(real colvector theta, real matrix W,
    string scalar focal_type, real scalar cf_col, real scalar cfm_col, real scalar has_m)
{
    struct pape_est scalar e
    real matrix W1, W0
    real colvector a1, a0, a, slope, ga, gp
    real scalar n
    n = rows(W)
    if (focal_type == "binary") {
        W1 = W; W1[., 2] = J(n, 1, 1)
        W0 = W; W0[., 2] = J(n, 1, 0)
        if (has_m) {
            W1[., 4] = W[., 3]
            W0[., 4] = J(n, 1, 0)
        }
        a1 = W1 * theta; a0 = W0 * theta
        e.ape = mean(normal(a1) :- normal(a0))
        e.jac = (mean(normalden(a1) :* W1) :- mean(normalden(a0) :* W0))'
        return(e)
    }
    a = W * theta
    if (has_m) slope = theta[2] :+ theta[4] :* W[., 3]
    else slope = J(n, 1, theta[2])
    ga = normalden(a)
    e.ape = mean(ga :* slope)
    gp = -a :* ga
    e.jac = (mean(gp :* slope :* W))'
    e.jac[2] = e.jac[2] + mean(ga)
    if (has_m) e.jac[4] = e.jac[4] + mean(ga :* W[., 3])
    return(e)
}

real matrix pape_cf_cell(real matrix W, real colvector cf, real scalar dval,
    real scalar mval, real scalar cfm_col)
{
    real matrix Wc
    real scalar n
    n = rows(W)
    Wc = W
    Wc[., 2] = J(n, 1, dval)
    Wc[., 3] = J(n, 1, mval)
    Wc[., 4] = J(n, 1, dval * mval)
    Wc[., cfm_col] = cf :* mval
    return(Wc)
}

// ASF-based AIE (endogenous focal x exogenous binary moderator)
struct pape_est scalar pape_cf_aie_est(real colvector theta, real matrix W,
    string scalar focal_type, real scalar cf_col, real scalar cfm_col)
{
    struct pape_est scalar e
    real matrix W11, W01, W10, W00, Wc
    real colvector cf, a11, a01, a10, a00, d_obs, a, ga, gp
    real colvector jac1, jac0
    real scalar v1, v0, mval, slope, j, n
    cf = W[., cf_col]
    n = rows(W)
    if (focal_type == "binary") {
        W11 = pape_cf_cell(W, cf, 1, 1, cfm_col)
        W01 = pape_cf_cell(W, cf, 0, 1, cfm_col)
        W10 = pape_cf_cell(W, cf, 1, 0, cfm_col)
        W00 = pape_cf_cell(W, cf, 0, 0, cfm_col)
        a11 = W11 * theta; a01 = W01 * theta; a10 = W10 * theta; a00 = W00 * theta
        e.ape = mean(normal(a11) :- normal(a01) :- normal(a10) :+ normal(a00))
        e.jac = (mean(normalden(a11) :* W11) :- mean(normalden(a01) :* W01) :-
                 mean(normalden(a10) :* W10) :+ mean(normalden(a00) :* W00))'
        return(e)
    }
    d_obs = W[., 2]
    for (j = 1; j <= 2; j++) {
        mval = (j == 1 ? 1 : 0)
        Wc = W
        Wc[., 3] = J(n, 1, mval)
        Wc[., 4] = d_obs :* mval
        Wc[., cfm_col] = cf :* mval
        a = Wc * theta
        slope = theta[2] + theta[4] * mval
        ga = normalden(a)
        gp = -a :* ga
        if (j == 1) {
            jac1 = (mean(gp :* slope :* Wc))'
            jac1[2] = jac1[2] + mean(ga)
            jac1[4] = jac1[4] + mean(ga) * mval
            v1 = mean(ga :* slope)
        }
        else {
            jac0 = (mean(gp :* slope :* Wc))'
            jac0[2] = jac0[2] + mean(ga)
            jac0[4] = jac0[4] + mean(ga) * mval
            v0 = mean(ga :* slope)
        }
    }
    e.ape = v1 - v0
    e.jac = jac1 - jac0
    return(e)
}

// ---------------------------------------------------------------------------
// IV DGP: draws, integration, construction (iv.R)
// ---------------------------------------------------------------------------

struct pape_ivdraw scalar pape_iv_draw(struct pape_dgp scalar d, real scalar n, real scalar seed)
{
    struct pape_ivdraw scalar rd
    real matrix U
    real colvector qz
    real scalar j, k, nz
    string scalar st
    if (seed != .) {
        st = rngstate()
        rseed(seed)
    }
    U = pape_draw_std_mvn(n, d.Rchol)
    k = d.k; nz = d.n_inst
    if (d.has_m) rd.m = pape_transform_u(U[., 1], d.moderator)
    else rd.m = J(0, 1, .)
    if (k > 0) {
        rd.X = U[., (d.has_m + 1)..(d.has_m + k)]
        for (j = 1; j <= k; j++) rd.X[., j] = pape_transform_u(rd.X[., j], d.covariates[j])
    }
    else rd.X = J(n, 0, .)
    rd.Zi = U[., (d.has_m + k + 1)..(d.has_m + k + nz)]
    for (j = 1; j <= nz; j++) rd.Zi[., j] = pape_transform_u(rd.Zi[., j], d.instruments[j])
    qz = rd.Zi * d.pi_z
    if (d.focal.type == "binary") {
        rd.v = rnormal(n, 1, 0, 1)
        rd.d = (d.kappa0 :+ qz :+ rd.v :> 0)
    }
    else {
        rd.v = rnormal(n, 1, 0, d.sigma_v)
        rd.d = d.kappa0 :+ qz :+ rd.v
    }
    if (seed != .) rngstate(st)
    return(rd)
}

struct pape_idraw scalar pape_iv_integration(struct pape_dgp scalar d)
{
    struct pape_idraw scalar id
    struct pape_ivdraw scalar rd
    rd = pape_iv_draw(d, d.n_int, d.seed_int)
    if (d.k > 0) id.idxz = rd.X * d.gamma
    else id.idxz = J(d.n_int, 1, 0)
    id.xd_c = rd.d :- d.focal_ref
    id.has_xd = 1
    if (d.has_m) {
        id.md_c = rd.m :- d.mod_ref
        id.has_md = 1
    }
    else id.has_md = 0
    return(id)
}

real scalar pape_f_kappa0(real scalar b, struct pape_env scalar e)
    return(mean(normal(b :+ e.idxz)) - e.target)

// vars: focal, [moderator], covariates; instruments separately
struct pape_dgp scalar pape_dgp_iv(string scalar name, struct pape_var colvector vars,
    real scalar has_m, struct pape_var colvector instruments,
    real scalar endogeneity, real scalar iv_strength,
    real scalar corr_given, real scalar corr_scalar, real matrix corr_matrix,
    real scalar baseline, real scalar signal, real scalar n_int, real scalar seed_int)
{
    struct pape_dgp scalar d
    struct pape_env scalar e
    real matrix U, Xc, Zi
    real colvector w_unit, g_unit, qz, idxz
    real scalar k, nz, j, v_unit, var_qz, sd_d, v_g
    string scalar st
    d = pape_dgp_blank()
    d.name = name
    d.model = "probit"; d.link = "probit"
    d.route = "iv"; d.builder = "iv"
    d.focal = vars[1]
    d.has_m = has_m
    if (has_m) {
        d.moderator = vars[2]
        if (d.moderator.type != "binary") pape_err("IV AIE designs currently require a binary moderator.")
        if (d.moderator.name == d.focal.name) pape_err("moderator() must be a different variable than focal().")
    }
    k = rows(vars) - 1 - has_m
    d.k = k
    if (k > 0) d.covariates = vars[(2 + has_m)..rows(vars)]
    nz = rows(instruments)
    if (nz < 1) pape_err("instruments() needs at least one instrument.")
    d.instruments = instruments
    d.n_inst = nz
    if (endogeneity == . | abs(endogeneity) >= 1) pape_err("endogeneity() must lie in (-1, 1).")
    if (iv_strength == . | iv_strength <= 0 | iv_strength >= 1) pape_err("ivstrength() must lie in (0, 1).")
    if (baseline == . | baseline <= 0 | baseline >= 1) pape_err("baseline() must lie strictly between 0 and 1.")
    if (signal == . | signal < 0 | signal >= 1) pape_err("signal() must lie in [0, 1).")
    if (k == 0 & signal > 0) pape_err("signal() > 0 requires at least one covariate.")
    d.R = pape_build_corr(has_m + k + nz, corr_given, corr_scalar, corr_matrix)
    d.Rchol = cholesky(d.R)'
    d.corr_given = corr_given; d.corr_scalar = corr_scalar; d.corr_matrix = corr_matrix
    d.n_int = n_int; d.seed_int = seed_int
    d.endogeneity = endogeneity; d.iv_strength = iv_strength
    d.baseline = baseline; d.signal = signal
    d.first_stage = (d.focal.type == "binary" ? "probit" : "linear")
    d.focal_ref = (d.focal.type == "binary" ? 0 : d.focal.mean)
    d.mod_ref = (has_m ? 0 : .)
    st = rngstate()
    rseed(seed_int)
    U = pape_draw_std_mvn(n_int, d.Rchol)
    rngstate(st)
    if (k > 0) {
        Xc = U[., (has_m + 1)..(has_m + k)]
        for (j = 1; j <= k; j++) Xc[., j] = pape_transform_u(Xc[., j], d.covariates[j])
    }
    else Xc = J(n_int, 0, .)
    Zi = U[., (has_m + k + 1)..(has_m + k + nz)]
    for (j = 1; j <= nz; j++) Zi[., j] = pape_transform_u(Zi[., j], instruments[j])
    w_unit = J(nz, 1, .)
    for (j = 1; j <= nz; j++) w_unit[j] = 1 / pape_var_sd(instruments[j])
    v_unit = variance(Zi * w_unit)
    if (d.focal.type == "binary") {
        var_qz = iv_strength / (1 - iv_strength)
        d.pi_z = sqrt(var_qz / v_unit) :* w_unit
        d.sigma_v = 1
        qz = Zi * d.pi_z
        e.idxz = qz; e.target = d.focal.p
        d.kappa0 = pape_uniroot(&pape_f_kappa0(), e, -10, 10, 1, 1e-10)
    }
    else {
        sd_d = pape_var_sd(d.focal)
        d.pi_z = sqrt(iv_strength * sd_d^2 / v_unit) :* w_unit
        d.sigma_v = sqrt((1 - iv_strength) * sd_d^2)
        qz = Zi * d.pi_z
        d.kappa0 = d.focal.mean - mean(qz)
    }
    if (k > 0) {
        if (signal == 0) d.gamma = J(k, 1, 0)
        else {
            g_unit = J(k, 1, .)
            for (j = 1; j <= k; j++) g_unit[j] = 1 / pape_var_sd(d.covariates[j])
            v_g = variance(Xc * g_unit)
            d.gamma = sqrt((signal / (1 - signal)) / v_g) :* g_unit
        }
        idxz = Xc * d.gamma
    }
    else {
        d.gamma = J(0, 1, .)
        idxz = J(n_int, 1, 0)
    }
    d.s2 = (k > 0 ? variance(idxz) : 0)
    d.signal_implied = d.s2 / (1 + d.s2)
    e.link = "probit"; e.idxz = idxz; e.target = baseline
    d.beta0 = pape_uniroot(&pape_f_b0(), e, -15, 15, 1, 1e-10)
    return(d)
}

// One simulated study: exogenous block + first stage + structural outcome
struct pape_xdraw scalar pape_draw_x_iv(struct pape_dgp scalar d, real scalar n)
{
    struct pape_xdraw scalar xx
    struct pape_ivdraw scalar rd
    real colvector v_std, u, d_c, eta, start
    real scalar rho, fr, mr, s2s, lam
    if (!d.has_effect) pape_err("This DGP has no target effect yet; run powerape setape (or setaie) first.")
    rd = pape_iv_draw(d, n, .)
    rho = d.endogeneity
    v_std = (d.focal.type == "binary" ? rd.v : rd.v :/ d.sigma_v)
    u = rho :* v_std :+ sqrt(1 - rho^2) :* rnormal(n, 1, 0, 1)
    fr = d.focal_ref
    mr = (d.has_m ? d.mod_ref : 0)
    d_c = rd.d :- fr
    eta = d.beta0 :+ d.beta_focal :* d_c
    if (d.k > 0) eta = eta :+ rd.X * d.gamma
    if (d.has_m) eta = eta :+ d.beta_mod :* (rd.m :- mr) :+ d.beta_int :* d_c :* (rd.m :- mr)
    xx.y = (eta :+ u :> 0)
    s2s = sqrt(1 - rho^2)
    lam = (d.focal.type == "binary" ? rho / s2s : rho / (d.sigma_v * s2s))
    if (!d.has_m) start = ((d.beta0 - d.beta_focal * fr) \ d.beta_focal \ d.gamma) :/ s2s \ lam
    else start = ((d.beta0 - d.beta_focal * fr - d.beta_mod * mr + d.beta_int * fr * mr) \
                  (d.beta_focal - d.beta_int * mr) \ (d.beta_mod - d.beta_int * fr) \
                  d.beta_int \ d.gamma) :/ s2s \ lam \ 0
    xx.start = start
    xx.has_start = 1
    xx.d = rd.d
    xx.m = rd.m
    xx.Zx = rd.X
    if (d.has_m) xx.X1 = (J(n, 1, 1), rd.X, rd.m, rd.Zi, rd.Zi :* rd.m)
    else xx.X1 = (J(n, 1, 1), rd.X, rd.Zi)
    xx.focal_col = 2
    xx.has_id = 0
    return(xx)
}

end
