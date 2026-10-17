*! powerape_core.mata -- Stata/Mata port of the powerape R package (core)
*! Structures, link functions, copula draws, root-finding, DGP calibration,
*! and the APE/AIE inversion. Mirrors R/dgp.R, R/var.R, R/utils.R,
*! R/invert.R, R/aie.R, R/empirical.R, R/from-fit.R of powerape 1.11.3.
version 16
mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// Structures
// ---------------------------------------------------------------------------

struct pape_var {
    string scalar name
    string scalar type      // "binary" | "normal"
    real scalar p
    real scalar mean
    real scalar sd
    real scalar icc         // . when not set
}

struct pape_dgp {
    string scalar name
    string scalar model, link, route, builder, estimand
    struct pape_var scalar focal, moderator
    real scalar has_m, k
    struct pape_var colvector covariates     // parametric / panel / iv routes
    string colvector covnames                // empirical / from_fit routes
    real matrix R, Rchol                     // Rchol upper: R = Rchol' Rchol
    real scalar corr_given                   // 0 none, 1 scalar, 2 matrix
    real scalar corr_scalar
    real matrix corr_matrix
    real scalar baseline, signal, signal_implied, s2
    real scalar gamma_given
    real colvector gamma
    real scalar beta0, beta_focal, beta_mod, beta_int
    real scalar has_effect                   // 1 after set_ape / set_aie
    real scalar target_est, p1, main_focal, main_moderator
    real scalar focal_ref, mod_ref
    real scalar n_int, seed_int
    // empirical / from_fit
    real matrix emp_x
    real colvector emp_xd
    real scalar n_emp, emp_focal_joint
    string scalar formula_str
    // panel
    real scalar n_periods, rho, cre_share, unbalanced, expected_T, retention
    real rowvector t_support, t_prob
    real scalar var_a, scale_a
    real colvector iccs
    real colvector tv_idx
    real colvector xi
    real rowvector xbar_center
    // iv
    struct pape_var colvector instruments
    real scalar n_inst, endogeneity, iv_strength
    string scalar first_stage
    real colvector pi_z
    real scalar kappa0, sigma_v
    // spec (for rebuilds in the robustness sweep)
    string scalar spec_opts
}

struct pape_idraw {
    real colvector xd_c, md_c, idxz
    real scalar has_xd, has_md
    real colvector q                        // panel: idxz + Mundlak part
}

// generic environment for root-finding objective functions
struct pape_env {
    string scalar link
    real scalar b0, target, p_ref, scale, bf, bm, bi, z, p, r, s, t
    real colvector idxz, x_c, m_c
    struct pape_dgp scalar dgp
    struct pape_idraw scalar id
}

// ---------------------------------------------------------------------------
// Errors and globals
// ---------------------------------------------------------------------------

void pape_err(string scalar msg)
{
    pointer(transmorphic) scalar p
    real scalar quiet
    p = findexternal("PAPE_LASTERR")
    if (p == NULL) p = crexternal("PAPE_LASTERR")
    *p = msg
    quiet = 0
    p = findexternal("PAPE_QUIET")
    if (p != NULL) quiet = (*p == 1)
    if (!quiet) errprintf("%s\n", msg)
    exit(498)
}

// Storage of DGPs and results: associative arrays held as Mata externals.
transmorphic pape_store(string scalar which)
{
    pointer(transmorphic) scalar p
    string scalar nm
    nm = "PAPE_STORE_" + which
    p = findexternal(nm)
    if (p == NULL) {
        p = crexternal(nm)
        *p = asarray_create()
    }
    return(*p)
}

void pape_store_put(string scalar which, string scalar key, transmorphic val)
{
    pointer(transmorphic) scalar p
    string scalar nm
    nm = "PAPE_STORE_" + which
    p = findexternal(nm)
    if (p == NULL) {
        p = crexternal(nm)
        *p = asarray_create()
    }
    asarray(*p, key, val)
}

real scalar pape_store_has(string scalar which, string scalar key)
{
    pointer(transmorphic) scalar p
    p = findexternal("PAPE_STORE_" + which)
    if (p == NULL) return(0)
    return(asarray_contains(*p, key))
}

void pape_store_drop(string scalar which, string scalar key)
{
    pointer(transmorphic) scalar p
    p = findexternal("PAPE_STORE_" + which)
    if (p == NULL) return
    if (asarray_contains(*p, key)) asarray_remove(*p, key)
}

string colvector pape_store_keys(string scalar which)
{
    pointer(transmorphic) scalar p
    string colvector k
    p = findexternal("PAPE_STORE_" + which)
    if (p == NULL) return(J(0, 1, ""))
    k = asarray_keys(*p)
    if (rows(k) == 0) return(J(0, 1, ""))
    return(sort(k, 1))
}

struct pape_dgp scalar pape_get_dgp(string scalar name)
{
    struct pape_dgp scalar d
    if (!pape_store_has("DGP", name))
        pape_err(sprintf("DGP '%s' not found; create it with powerape dgp.", name))
    d = asarray(pape_store("DGP"), name)
    return(d)
}

void pape_put_dgp(struct pape_dgp scalar d)
{
    pape_store_put("DGP", d.name, d)
}

// ---------------------------------------------------------------------------
// Links
// ---------------------------------------------------------------------------

real matrix pape_G(string scalar link, real matrix x)
{
    if (link == "probit") return(normal(x))
    return(logistic(x))
}

real matrix pape_g(string scalar link, real matrix x)
{
    if (link == "probit") return(normalden(x))
    return(logisticden(x))
}

real scalar pape_Ginv(string scalar link, real scalar p)
{
    if (link == "probit") return(invnormal(p))
    return(logit(p))
}

// g'(z): probit -z phi(z); logit g(1 - 2G)
real matrix pape_gprime(string scalar link, real matrix z)
{
    if (link == "probit") return(-z :* normalden(z))
    return(logisticden(z) :* (1 :- 2 :* logistic(z)))
}

// g''(z): probit (z^2 - 1) phi(z); logit g(1-2G)^2 - 2 g^2
real matrix pape_gpp(string scalar link, real matrix z)
{
    real matrix g
    if (link == "probit") return((z:^2 :- 1) :* normalden(z))
    g = logisticden(z)
    return(g :* (1 :- 2 :* logistic(z)):^2 :- 2 :* g:^2)
}

real scalar pape_zcrit(real scalar conf) return(invnormal(1 - (1 - conf) / 2))

real scalar pape_mz_r2_logit(real scalar s2) return(s2 / (s2 + pi()^2 / 3))

// ---------------------------------------------------------------------------
// Variables and the Gaussian copula
// ---------------------------------------------------------------------------

struct pape_var scalar pape_var_new(string scalar name, string scalar type,
    real scalar p, real scalar mean, real scalar sd, real scalar icc)
{
    struct pape_var scalar v
    v.name = name
    v.type = type
    v.p = p
    v.mean = mean
    v.sd = sd
    v.icc = icc
    if (type == "binary") {
        if (p == . | p <= 0 | p >= 1)
            pape_err(sprintf("A binary variable (%s) needs p strictly between 0 and 1.", name))
        v.mean = .; v.sd = .
    }
    else if (type == "normal") {
        if (mean == .) v.mean = 0
        if (sd == .) v.sd = 1
        if (v.sd <= 0) pape_err(sprintf("A normal variable (%s) needs sd > 0.", name))
        v.p = .
    }
    else pape_err(sprintf("Unknown variable type '%s' for %s.", type, name))
    if (icc != . & (icc < 0 | icc > 1))
        pape_err(sprintf("icc for %s must lie in [0, 1].", name))
    return(v)
}

real scalar pape_var_sd(struct pape_var scalar v)
{
    if (v.type == "binary") return(sqrt(v.p * (1 - v.p)))
    return(v.sd)
}

// Map standard-normal draws to the marginal (Gaussian copula).
real colvector pape_transform_u(real colvector u, struct pape_var scalar v)
{
    if (v.type == "binary") return(u :> invnormal(1 - v.p))
    return(v.mean :+ v.sd :* u)
}

string scalar pape_var_label(struct pape_var scalar v)
{
    string scalar s
    if (v.type == "binary") s = sprintf("%s (binary, p = %4.2f", v.name, v.p)
    else s = sprintf("%s (normal, mean %4.2f, sd %4.2f", v.name, v.mean, v.sd)
    if (v.icc != .) s = sprintf("%s, icc = %4.2f", s, v.icc)
    return(s + ")")
}

real matrix pape_build_corr(real scalar K, real scalar corr_given,
    real scalar corr_scalar, real matrix corr_matrix)
{
    real matrix R, L
    real scalar rc
    if (corr_given == 0) R = I(K)
    else if (corr_given == 1) {
        R = J(K, K, corr_scalar)
        _diag(R, 1)
    }
    else R = corr_matrix
    if (rows(R) != K | cols(R) != K)
        pape_err(sprintf("correlation() must be a %f x %f matrix over (focal%s, covariates).",
            K, K, (K > 1 ? ", moderator" : "")))
    if (max(abs(R - R')) > 1e-8 | max(abs(diagonal(R) :- 1)) > 1e-8)
        pape_err("correlation() must be symmetric with unit diagonal.")
    L = cholesky(R)
    if (hasmissing(L)) pape_err("correlation() is not positive definite.")
    return(R)
}

// n x K standard-normal rows with correlation R (Rchol upper triangular).
real matrix pape_draw_std_mvn(real scalar n, real matrix Rchol)
{
    return(rnormal(n, cols(Rchol), 0, 1) * Rchol)
}

// Deterministic parametric regressor draw under seed (seed = . : ambient).
void pape_reg_draw_par(struct pape_dgp scalar d, real scalar n, real scalar seed,
    real colvector xd, real colvector md, real matrix Z)
{
    real matrix U
    real scalar j, zoff
    string scalar st
    if (seed != .) {
        st = rngstate()
        rseed(seed)
    }
    U = pape_draw_std_mvn(n, d.Rchol)
    xd = pape_transform_u(U[., 1], d.focal)
    if (d.has_m) md = pape_transform_u(U[., 2], d.moderator)
    else md = J(0, 1, .)
    zoff = 1 + d.has_m
    if (d.k > 0) {
        Z = U[., (zoff + 1)..(zoff + d.k)]
        for (j = 1; j <= d.k; j++) Z[., j] = pape_transform_u(Z[., j], d.covariates[j])
    }
    else Z = J(n, 0, .)
    if (seed != .) rngstate(st)
}

// Integration set: centered focal/moderator draws + nuisance index.
struct pape_idraw scalar pape_integration_draw(struct pape_dgp scalar d)
{
    struct pape_idraw scalar id
    real colvector xd, md
    real matrix Z
    real scalar m
    string scalar st
    if (d.route == "panel") return(pape_panel_integration(d))
    if (d.route == "iv") return(pape_iv_integration(d))
    if (d.route == "empirical") {
        if (d.k > 0) id.idxz = d.emp_x * d.gamma
        else id.idxz = J(d.n_emp, 1, 0)
        id.has_md = 0
        if (d.emp_focal_joint) {
            id.xd_c = d.emp_xd :- d.focal_ref
            id.has_xd = 1
        }
        else if (d.focal.type != "binary") {
            m = max((1, ceil(d.n_int / rows(id.idxz))))
            id.idxz = J(m, 1, id.idxz)
            st = rngstate()
            rseed(d.seed_int)
            id.xd_c = pape_transform_u(rnormal(rows(id.idxz), 1, 0, 1), d.focal) :- d.focal_ref
            rngstate(st)
            id.has_xd = 1
        }
        else id.has_xd = 0
        return(id)
    }
    pape_reg_draw_par(d, d.n_int, d.seed_int, xd, md, Z)
    if (d.k > 0) id.idxz = Z * d.gamma
    else id.idxz = J(d.n_int, 1, 0)
    id.xd_c = xd :- d.focal_ref
    id.has_xd = 1
    if (d.has_m) {
        id.md_c = md :- d.mod_ref
        id.has_md = 1
    }
    else id.has_md = 0
    return(id)
}

// ---------------------------------------------------------------------------
// Root finding (Brent's zeroin, as in R's uniroot) with bracket extension
// ---------------------------------------------------------------------------

// extend: 0 none, 1 "upX" (f increasing in x), 2 "yes" (either direction)
real scalar pape_uniroot(pointer(real scalar function) scalar f,
    struct pape_env scalar e, real scalar lo, real scalar hi,
    real scalar extend, real scalar tol)
{
    real scalar a, b, c, fa, fb, fc, prev_step, tol_act, new_step, cb
    real scalar t1, t2, p, q, EPS, it, maxit, delta
    EPS = epsilon(1)
    a = lo; b = hi
    fa = (*f)(a, e); fb = (*f)(b, e)
    if (extend == 1) {
        delta = 0.01 * max((1e-4, abs(a)))
        it = 0
        while (fa > 0 & it < 100) {
            a = a - delta; fa = (*f)(a, e); delta = 2 * delta; it++
        }
        delta = 0.01 * max((1e-4, abs(b)))
        it = 0
        while (fb < 0 & it < 100) {
            b = b + delta; fb = (*f)(b, e); delta = 2 * delta; it++
        }
    }
    else if (extend == 2) {
        delta = 0.01 * max((1e-4, abs(b - a)))
        it = 0
        while (sign(fa) == sign(fb) & fa != 0 & fb != 0 & it < 100) {
            a = a - delta; fa = (*f)(a, e)
            if (sign(fa) != sign(fb)) break
            b = b + delta; fb = (*f)(b, e)
            delta = 2 * delta; it++
        }
    }
    if (fa == 0) return(a)
    if (fb == 0) return(b)
    if (sign(fa) == sign(fb))
        pape_err("Root-finding failed: f() has the same sign at both bracket ends.")
    c = a; fc = fa
    maxit = 1000
    while (maxit-- > 0) {
        prev_step = b - a
        if (abs(fc) < abs(fb)) {
            a = b; b = c; c = a
            fa = fb; fb = fc; fc = fa
        }
        tol_act = 2 * EPS * abs(b) + tol / 2
        new_step = (c - b) / 2
        if (abs(new_step) <= tol_act | fb == 0) return(b)
        if (abs(prev_step) >= tol_act & abs(fa) > abs(fb)) {
            cb = c - b
            if (a == c) {
                t1 = fb / fa
                p = cb * t1
                q = 1 - t1
            }
            else {
                q = fa / fc; t1 = fb / fc; t2 = fb / fa
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1))
                q = (q - 1) * (t1 - 1) * (t2 - 1)
            }
            if (p > 0) q = -q
            else p = -p
            if (p < (0.75 * cb * q - abs(tol_act * q) / 2) & p < abs(prev_step * q / 2))
                new_step = p / q
        }
        if (abs(new_step) < tol_act) new_step = (new_step > 0 ? tol_act : -tol_act)
        a = b; fa = fb
        b = b + new_step
        fb = (*f)(b, e)
        if ((fb > 0 & fc > 0) | (fb < 0 & fc < 0)) {
            c = a; fc = fa
        }
    }
    return(b)
}

// ---------------------------------------------------------------------------
// Inverse power maps (utils.R: req_delta, req_margin)
// ---------------------------------------------------------------------------

real scalar pape_f_reqdelta(real scalar dd, struct pape_env scalar e)
    return(normal(dd - e.z) + normal(e.r * dd - e.z) - 1 - e.p)

real scalar pape_req_delta(real scalar p, real scalar z, real scalar r)
{
    struct pape_env scalar e
    real scalar hi
    if (r == .) return(z + invnormal(p))
    hi = z + invnormal((1 + p) / 2)
    e.z = z; e.r = r; e.p = p
    if (pape_f_reqdelta(hi, e) <= 1e-12) return(hi)
    return(pape_uniroot(&pape_f_reqdelta(), e, max((0, z + invnormal(p))), hi, 0, 1e-12))
}

real scalar pape_f_reqmargin(real scalar m, struct pape_env scalar e)
    return(normal((m - e.t) / e.s - e.z) + normal((m + e.t) / e.s - e.z) - 1 - e.p)

real scalar pape_req_margin(real scalar p, real scalar s, real scalar t, real scalar z)
{
    struct pape_env scalar e
    real scalar hi
    e.p = p; e.s = s; e.t = t; e.z = z
    hi = t + s * (z + invnormal((1 + p) / 2))
    if (pape_f_reqmargin(hi, e) <= 1e-12) return(hi)
    return(pape_uniroot(&pape_f_reqmargin(), e, t + s * max((0, z + invnormal(p))), hi, 0, 1e-12))
}

// ---------------------------------------------------------------------------
// DGP construction: parametric route (ape_dgp)
// ---------------------------------------------------------------------------

real scalar pape_f_b0(real scalar b0, struct pape_env scalar e)
    return(mean(pape_G(e.link, b0 :+ e.idxz)) - e.target)

struct pape_dgp scalar pape_dgp_blank()
{
    struct pape_dgp scalar d
    d.has_m = 0; d.k = 0
    d.has_effect = 0
    d.estimand = ""
    d.beta_focal = .; d.beta_mod = .; d.beta_int = .
    d.target_est = .; d.p1 = .; d.main_focal = .; d.main_moderator = .
    d.gamma = J(0, 1, .)
    d.gamma_given = 0
    d.covariates = pape_var(0, 1)
    d.covnames = J(0, 1, "")
    d.corr_given = 0; d.corr_scalar = .; d.corr_matrix = J(0, 0, .)
    d.n_emp = 0; d.emp_focal_joint = 0
    d.emp_x = J(0, 0, .); d.emp_xd = J(0, 1, .)
    d.formula_str = ""
    d.n_periods = .; d.rho = 0; d.cre_share = 0; d.unbalanced = 0
    d.expected_T = .; d.retention = .
    d.t_support = J(1, 0, .); d.t_prob = J(1, 0, .)
    d.var_a = 0; d.scale_a = 1
    d.iccs = J(0, 1, .); d.tv_idx = J(0, 1, .); d.xi = J(0, 1, .)
    d.xbar_center = J(1, 0, .)
    d.instruments = pape_var(0, 1)
    d.n_inst = 0; d.endogeneity = .; d.iv_strength = .
    d.first_stage = ""
    d.pi_z = J(0, 1, .); d.kappa0 = .; d.sigma_v = .
    d.spec_opts = ""
    d.mod_ref = .
    return(d)
}

// vars: focal first, moderator second if has_m, then covariates
struct pape_dgp scalar pape_dgp_par(string scalar name, string scalar model,
    struct pape_var colvector vars, real scalar has_m,
    real scalar corr_given, real scalar corr_scalar, real matrix corr_matrix,
    real scalar baseline, real scalar signal, real scalar gamma_given,
    real colvector gamma, real scalar n_int, real scalar seed_int)
{
    struct pape_dgp scalar d
    struct pape_env scalar e
    real colvector xd, md, g_unit, idx
    real matrix Z
    real scalar K, j, v_unit
    d = pape_dgp_blank()
    d.name = name
    if (model != "probit" & model != "logit") pape_err("model() must be probit or logit.")
    d.model = model; d.link = model
    d.route = "parametric"; d.builder = "parametric"
    d.focal = vars[1]
    d.has_m = has_m
    if (has_m) {
        d.moderator = vars[2]
        if (d.moderator.name == d.focal.name)
            pape_err("moderator() must be a different variable than focal().")
    }
    d.k = rows(vars) - 1 - has_m
    if (d.k > 0) d.covariates = vars[(2 + has_m)..rows(vars)]
    if (baseline == . | baseline <= 0 | baseline >= 1)
        pape_err("baseline() must lie strictly between 0 and 1.")
    if (n_int < 1e4) pape_err("nint() must be at least 10000.")
    K = 1 + has_m + d.k
    d.R = pape_build_corr(K, corr_given, corr_scalar, corr_matrix)
    d.Rchol = cholesky(d.R)'
    d.corr_given = corr_given; d.corr_scalar = corr_scalar; d.corr_matrix = corr_matrix
    d.baseline = baseline
    d.n_int = n_int; d.seed_int = seed_int
    d.gamma_given = gamma_given
    if (!gamma_given) {
        if (signal == . | signal < 0 | signal >= 1)
            pape_err("signal() must lie in [0, 1).")
        if (d.k == 0 & signal > 0)
            pape_err("signal() > 0 requires at least one covariate.")
        d.signal = signal
    }
    else {
        if (rows(gamma) != d.k)
            pape_err(sprintf("gamma() must have %f element(s), one per covariate.", d.k))
        d.signal = .
    }
    if (d.k == 0) {
        d.gamma = J(0, 1, .)
        d.s2 = 0
        d.signal_implied = 0
        d.beta0 = pape_Ginv(model, baseline)
    }
    else {
        pape_reg_draw_par(d, n_int, seed_int, xd, md, Z)
        if (!gamma_given) {
            if (signal == 0) d.gamma = J(d.k, 1, 0)
            else {
                g_unit = J(d.k, 1, .)
                for (j = 1; j <= d.k; j++) g_unit[j] = 1 / pape_var_sd(d.covariates[j])
                v_unit = variance(Z * g_unit)
                d.gamma = sqrt((signal / (1 - signal)) / v_unit) :* g_unit
            }
        }
        else d.gamma = gamma
        idx = Z * d.gamma
        d.s2 = variance(idx)
        d.signal_implied = d.s2 / (1 + d.s2)
        e.link = model; e.idxz = idx; e.target = baseline
        d.beta0 = pape_uniroot(&pape_f_b0(), e, -15, 15, 1, 1e-10)
    }
    d.focal_ref = (d.focal.type == "binary" ? 0 : d.focal.mean)
    if (has_m) d.mod_ref = (d.moderator.type == "binary" ? 0 : d.moderator.mean)
    return(d)
}

// ---------------------------------------------------------------------------
// DGP construction: empirical route (ape_dgp_empirical)
// ---------------------------------------------------------------------------

// M: pilot matrix; focal_col: column index of the focal in M (0 = none,
// focal drawn independently from `fvar`); covariates = the other columns.
struct pape_dgp scalar pape_dgp_emp(string scalar name, string scalar model,
    real matrix M, string colvector colnames, real scalar focal_col,
    struct pape_var scalar fvar, real scalar baseline, real scalar signal,
    real scalar gamma_given, real colvector gamma, real scalar seed_int)
{
    struct pape_dgp scalar d
    struct pape_env scalar e
    real colvector xd, g_unit, sds, idxz, keep
    real matrix covs
    real scalar k, j, v_unit, nr
    d = pape_dgp_blank()
    d.name = name
    if (model != "probit" & model != "logit") pape_err("model() must be probit or logit.")
    d.model = model; d.link = model
    d.route = "empirical"; d.builder = "empirical"
    if (baseline == . | baseline <= 0 | baseline >= 1)
        pape_err("baseline() must lie strictly between 0 and 1.")
    nr = rows(M)
    if (nr == 0) pape_err("The pilot data have no rows.")
    if (hasmissing(M)) pape_err("The pilot data contain missing values; supply complete cases.")
    if (nr < 30) pape_err("The empirical route needs at least 30 pilot rows.")
    if (nr < 500)
        printf("{txt}note: fewer than 500 pilot rows: pilot sampling noise becomes part of the design (cf. Hancock & Feng, 2025).\n")
    if (focal_col > 0) {
        xd = M[., focal_col]
        keep = J(cols(M), 1, 1); keep[focal_col] = 0
        covs = select(M, keep')
        d.covnames = select(colnames, keep)
        if (rows(uniqrows(xd)) < 2) pape_err("The focal column is constant.")
        if (all(xd :== 0 :| xd :== 1))
            d.focal = pape_var_new(colnames[focal_col], "binary", mean(xd), ., ., .)
        else
            d.focal = pape_var_new(colnames[focal_col], "normal", ., mean(xd), sqrt(variance(xd)), .)
        d.emp_xd = xd
        d.emp_focal_joint = 1
    }
    else {
        d.focal = fvar
        covs = M
        d.covnames = colnames
        d.emp_xd = J(0, 1, .)
        d.emp_focal_joint = 0
    }
    k = cols(covs)
    d.k = k
    d.gamma_given = gamma_given
    if (!gamma_given) {
        if (signal == . | signal < 0 | signal >= 1) pape_err("signal() must lie in [0, 1).")
        if (k == 0 & signal > 0) pape_err("signal() > 0 requires at least one covariate column.")
        d.signal = signal
        if (k == 0 | signal == 0) d.gamma = J(k, 1, 0)
        else {
            sds = sqrt(diagonal(variance(covs)))
            if (any(sds :< 1e-12)) pape_err("Constant covariate column(s) in the pilot data.")
            g_unit = 1 :/ sds
            v_unit = variance(covs * g_unit)
            d.gamma = sqrt((signal / (1 - signal)) / v_unit) :* g_unit
        }
    }
    else {
        if (rows(gamma) != k) pape_err(sprintf("gamma() must have %f element(s).", k))
        d.gamma = gamma
        d.signal = .
    }
    if (k > 0) idxz = covs * d.gamma
    else idxz = J(nr, 1, 0)
    d.s2 = (k > 0 ? variance(idxz) : 0)
    d.signal_implied = d.s2 / (1 + d.s2)
    e.link = model; e.idxz = idxz; e.target = baseline
    d.beta0 = pape_uniroot(&pape_f_b0(), e, -15, 15, 1, 1e-10)
    d.baseline = baseline
    d.emp_x = covs
    d.n_emp = nr
    d.focal_ref = (d.focal.type == "binary" ? 0 : (d.emp_focal_joint ? mean(xd) : d.focal.mean))
    d.n_int = 4e5; d.seed_int = seed_int
    return(d)
}

// ---------------------------------------------------------------------------
// DGP construction: pilot-model route (ape_dgp_from_fit)
// ---------------------------------------------------------------------------

// X: pilot regressors WITHOUT the intercept (focal + covariates), b: fitted
// coefficients in the same order followed by the intercept last.
struct pape_dgp scalar pape_dgp_fromfit(string scalar name, string scalar link,
    real matrix X, string colvector colnames, real scalar focal_col,
    real colvector b, real scalar b0_pilot, real scalar baseline_override,
    string scalar formula_str)
{
    struct pape_dgp scalar d
    struct pape_env scalar e
    real colvector xd, keep, idxz
    real matrix covs
    real scalar k, nr
    d = pape_dgp_blank()
    d.name = name
    d.model = link; d.link = link
    d.route = "empirical"; d.builder = "from_fit"
    nr = rows(X)
    if (nr < 30) pape_err("The pilot-model route needs at least 30 pilot rows.")
    if (nr < 500)
        printf("{txt}note: fewer than 500 pilot rows: pilot sampling noise becomes part of the design (cf. Hancock & Feng, 2025).\n")
    xd = X[., focal_col]
    if (rows(uniqrows(xd)) < 2) pape_err("The focal column is constant in the pilot data.")
    keep = J(cols(X), 1, 1); keep[focal_col] = 0
    covs = select(X, keep')
    d.covnames = select(colnames, keep)
    k = cols(covs)
    d.k = k
    if (all(xd :== 0 :| xd :== 1))
        d.focal = pape_var_new(colnames[focal_col], "binary", mean(xd), ., ., .)
    else
        d.focal = pape_var_new(colnames[focal_col], "normal", ., mean(xd), sqrt(variance(xd)), .)
    d.focal_ref = (d.focal.type == "binary" ? 0 : mean(xd))
    d.gamma = (k > 0 ? select(b, keep) : J(0, 1, .))
    d.gamma_given = 1
    d.signal = .
    if (k > 0) idxz = covs * d.gamma
    else idxz = J(nr, 1, 0)
    d.s2 = (k > 0 ? variance(idxz) : 0)
    d.signal_implied = d.s2 / (1 + d.s2)
    if (baseline_override == .) {
        d.beta0 = b0_pilot + b[focal_col] * d.focal_ref
        d.baseline = mean(pape_G(link, d.beta0 :+ idxz))
    }
    else {
        if (baseline_override <= 0 | baseline_override >= 1)
            pape_err("baseline() must lie strictly between 0 and 1.")
        e.link = link; e.idxz = idxz; e.target = baseline_override
        d.beta0 = pape_uniroot(&pape_f_b0(), e, -15, 15, 1, 1e-10)
        d.baseline = baseline_override
    }
    d.emp_x = covs
    d.emp_xd = xd
    d.emp_focal_joint = 1
    d.n_emp = nr
    d.formula_str = formula_str
    d.n_int = 4e5; d.seed_int = 20260812
    return(d)
}

// ---------------------------------------------------------------------------
// Inversion helpers (invert.R)
// ---------------------------------------------------------------------------

real scalar pape_f_bshift(real scalar b, struct pape_env scalar e)
    return(mean(pape_G(e.link, e.b0 :+ b :+ e.idxz)) - e.p_ref - e.target)

// returns 0 on success (root in `root`), else sets msg and returns 1
real scalar pape_solve_binary_shift(struct pape_dgp scalar d, real scalar target,
    real colvector idxz, real scalar p_ref, string scalar label,
    real scalar root, string scalar msg)
{
    struct pape_env scalar e
    real scalar eps
    eps = 1e-4
    if (target <= -p_ref + eps | target >= (1 - p_ref) - eps) {
        msg = sprintf("Target %s %9.4f is infeasible for this DGP: with baseline P(Y=1 | reference) = %5.3f the attainable range is (%6.3f, %5.3f). Adjust the target or the baseline.",
            label, target, p_ref, -p_ref, 1 - p_ref)
        return(1)
    }
    e.link = d.link; e.b0 = d.beta0; e.idxz = idxz; e.p_ref = p_ref; e.target = target
    root = pape_uniroot(&pape_f_bshift(), e, -8, 8, 1, 1e-9)
    msg = ""
    return(0)
}

real scalar pape_f_cont(real scalar b, struct pape_env scalar e)
    return(b * mean(pape_g(e.link, e.b0 :+ b :* e.x_c :+ e.idxz)) - e.target)

real scalar pape_solve_cont_main(struct pape_dgp scalar d, real scalar target,
    real colvector x_c, real colvector idxz, real scalar scale, string scalar label,
    real scalar root, string scalar msg)
{
    struct pape_env scalar e
    real colvector bgrid, vals
    real scalar s, i, amax, lo, hi, j
    e.link = d.link; e.b0 = d.beta0; e.idxz = idxz; e.x_c = x_c; e.target = 0
    s = (target < 0 ? -1 : 1)
    bgrid = s :* exp(rangen(ln(1e-3 / scale), ln(200 / scale), 100))
    vals = J(100, 1, .)
    for (j = 1; j <= 100; j++) vals[j] = pape_f_cont(bgrid[j], e)
    amax = max(abs(vals))
    if (abs(target) >= amax * (1 - 1e-3)) {
        msg = sprintf("Target %s %9.4f is infeasible for this DGP: the attainable range is about (%9.4f, %9.4f). The average derivative plateaus as the coefficient grows because the index disperses; adjust the target, the baseline, or the focal variable's spread.",
            label, target, -amax, amax)
        return(1)
    }
    i = min(selectindex(abs(vals) :>= abs(target)))
    lo = (i == 1 ? 0 : bgrid[i - 1])
    hi = bgrid[i]
    e.target = target
    if (lo > hi) {
        j = lo; lo = hi; hi = j
    }
    root = pape_uniroot(&pape_f_cont(), e, lo, hi, 0, 1e-9)
    msg = ""
    return(0)
}

// Interaction coefficient: f(bi) = truth_aie(bi) - target over a signed grid;
// pick the crossing nearest bi = 0.
real scalar pape_solve_interaction(pointer(real scalar function) scalar f,
    struct pape_env scalar e, real scalar target, real scalar scale,
    string scalar label, real scalar root, string scalar msg)
{
    real colvector gr, bs, vals, h, cross, mid
    real scalar j, nb, jj
    gr = exp(rangen(ln(1e-3 / scale), ln(50 / scale), 60))
    bs = -gr[rows(gr)..1] \ 0 \ gr
    nb = rows(bs)
    vals = J(nb, 1, .)
    for (j = 1; j <= nb; j++) vals[j] = (*f)(bs[j], e)
    h = vals :- target
    cross = selectindex(h[2..nb] :* h[1..(nb - 1)] :<= 0)
    if (rows(cross) == 0 | cols(cross) == 0) {
        msg = sprintf("Target %s %9.4f is infeasible for this DGP: the attainable range is about (%9.4f, %9.4f). Adjust the target, the main effects, or the baseline.",
            label, target, min(vals), max(vals))
        return(1)
    }
    cross = vec(cross)
    mid = abs((bs[cross] :+ bs[cross :+ 1]) :/ 2)
    jj = cross[selectindex(mid :== min(mid))[1]]
    root = pape_uniroot_target(f, e, bs[jj], bs[jj + 1], target)
    msg = ""
    return(0)
}

// Brent root of f(b) - target on [lo, hi]
real scalar pape_uniroot_target(pointer(real scalar function) scalar f,
    struct pape_env scalar e, real scalar lo, real scalar hi, real scalar target)
{
    real scalar a, b, c, fa, fb, fc, prev_step, tol_act, new_step, cb
    real scalar t1, t2, p, q, EPS, maxit, tol
    tol = 1e-9
    EPS = epsilon(1)
    a = lo; b = hi
    fa = (*f)(a, e) - target; fb = (*f)(b, e) - target
    if (fa == 0) return(a)
    if (fb == 0) return(b)
    c = a; fc = fa
    maxit = 1000
    while (maxit-- > 0) {
        prev_step = b - a
        if (abs(fc) < abs(fb)) {
            a = b; b = c; c = a
            fa = fb; fb = fc; fc = fa
        }
        tol_act = 2 * EPS * abs(b) + tol / 2
        new_step = (c - b) / 2
        if (abs(new_step) <= tol_act | fb == 0) return(b)
        if (abs(prev_step) >= tol_act & abs(fa) > abs(fb)) {
            cb = c - b
            if (a == c) {
                t1 = fb / fa
                p = cb * t1
                q = 1 - t1
            }
            else {
                q = fa / fc; t1 = fb / fc; t2 = fb / fa
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1))
                q = (q - 1) * (t1 - 1) * (t2 - 1)
            }
            if (p > 0) q = -q
            else p = -p
            if (p < (0.75 * cb * q - abs(tol_act * q) / 2) & p < abs(prev_step * q / 2))
                new_step = p / q
        }
        if (abs(new_step) < tol_act) new_step = (new_step > 0 ? tol_act : -tol_act)
        a = b; fa = fb
        b = b + new_step
        fb = (*f)(b, e) - target
        if ((fb > 0 & fc > 0) | (fb < 0 & fc < 0)) {
            c = a; fc = fa
        }
    }
    return(b)
}

// ---------------------------------------------------------------------------
// set_ape / true_ape
// ---------------------------------------------------------------------------

real scalar pape_set_ape(struct pape_dgp scalar d, real scalar target, string scalar msg)
{
    struct pape_idraw scalar id
    real colvector idxz
    real scalar p0, root, rc
    msg = ""
    if (target == .) {
        msg = "target() must be a number."
        return(1)
    }
    if (d.has_m) {
        msg = "This DGP has a moderator: pin the estimand with powerape setaie."
        return(1)
    }
    if (d.route == "panel") return(pape_set_ape_panel(d, target, msg))
    if (d.focal.type == "binary") {
        if (d.route == "parametric" & d.k == 0) {
            idxz = J(1, 1, 0)
            p0 = pape_G(d.link, d.beta0)
        }
        else {
            id = pape_integration_draw(d)
            idxz = id.idxz
            p0 = mean(pape_G(d.link, d.beta0 :+ idxz))
        }
        rc = pape_solve_binary_shift(d, target, idxz, p0, "APE", root, msg)
        if (rc) return(rc)
        d.p1 = mean(pape_G(d.link, d.beta0 :+ root :+ idxz))
    }
    else {
        id = pape_integration_draw(d)
        rc = pape_solve_cont_main(d, target, id.xd_c, id.idxz, pape_var_sd(d.focal), "APE", root, msg)
        if (rc) return(rc)
        d.p1 = .
    }
    d.beta_focal = root
    d.estimand = "ape"
    d.target_est = target
    d.has_effect = 1
    return(0)
}

real scalar pape_true_ape(struct pape_dgp scalar d)
{
    struct pape_idraw scalar id
    real colvector idxz
    if (d.estimand != "ape") pape_err("This DGP has no pinned APE; run powerape setape first.")
    if (d.route == "panel") return(pape_true_ape_panel(d))
    if (d.focal.type == "binary") {
        if (d.route == "parametric" & d.k == 0)
            return(pape_G(d.link, d.beta0 + d.beta_focal) - pape_G(d.link, d.beta0))
        id = pape_integration_draw(d)
        idxz = id.idxz
        return(mean(pape_G(d.link, d.beta0 :+ d.beta_focal :+ idxz) :- pape_G(d.link, d.beta0 :+ idxz)))
    }
    id = pape_integration_draw(d)
    return(d.beta_focal * mean(pape_g(d.link, d.beta0 :+ d.beta_focal :* id.xd_c :+ id.idxz)))
}

// ---------------------------------------------------------------------------
// AIE: truth and set_aie (aie.R)
// ---------------------------------------------------------------------------

real scalar pape_truth_aie(struct pape_dgp scalar d, real scalar bf, real scalar bm,
    real scalar bi, struct pape_idraw scalar id)
{
    string scalar link, td, tm
    real scalar b0
    real colvector xz, xd, md, idx
    link = d.link
    b0 = d.beta0
    xz = id.idxz
    td = d.focal.type; tm = d.moderator.type
    xd = id.xd_c; md = id.md_c
    if (td == "binary" & tm == "binary") {
        return(mean(pape_G(link, b0 + bf + bm + bi :+ xz)) - mean(pape_G(link, b0 + bm :+ xz)) -
               mean(pape_G(link, b0 + bf :+ xz)) + mean(pape_G(link, b0 :+ xz)))
    }
    else if (td != "binary" & tm == "binary") {
        return((bf + bi) * mean(pape_g(link, b0 :+ bf :* xd :+ bm :+ bi :* xd :+ xz)) -
               bf * mean(pape_g(link, b0 :+ bf :* xd :+ xz)))
    }
    else if (td == "binary" & tm != "binary") {
        return((bm + bi) * mean(pape_g(link, b0 :+ bf :+ bm :* md :+ bi :* md :+ xz)) -
               bm * mean(pape_g(link, b0 :+ bm :* md :+ xz)))
    }
    idx = b0 :+ bf :* xd :+ bm :* md :+ bi :* xd :* md :+ xz
    return(mean(pape_gprime(link, idx) :* (bf :+ bi :* md) :* (bm :+ bi :* xd) :+ pape_g(link, idx) :* bi))
}

real scalar pape_f_aie(real scalar b, struct pape_env scalar e)
    return(pape_truth_aie(e.dgp, e.bf, e.bm, b, e.id))

real scalar pape_set_aie(struct pape_dgp scalar d, real scalar target,
    real scalar main_focal, real scalar main_moderator, string scalar msg)
{
    struct pape_idraw scalar id
    struct pape_env scalar e
    real scalar p0, bf, bm, bi, rc, s_d, s_m
    msg = ""
    if (!d.has_m) {
        msg = "This DGP has no moderator: add one in powerape dgp, or pin an APE with powerape setape."
        return(1)
    }
    if (main_focal == . | main_moderator == .) {
        msg = "powerape setaie needs explicit main-effect anchors: mainfocal() and mainmoderator() are the conditional-at-reference APEs (the APE of each variable with the other held at its reference value)."
        return(1)
    }
    if (target == .) {
        msg = "target() must be a number."
        return(1)
    }
    if (d.route == "panel") return(pape_set_aie_panel(d, target, main_focal, main_moderator, msg))
    id = pape_integration_draw(d)
    p0 = mean(pape_G(d.link, d.beta0 :+ id.idxz))
    if (d.focal.type == "binary")
        rc = pape_solve_binary_shift(d, main_focal, id.idxz, p0, "main-effect APE (focal)", bf, msg)
    else
        rc = pape_solve_cont_main(d, main_focal, id.xd_c, id.idxz, pape_var_sd(d.focal), "main-effect APE (focal)", bf, msg)
    if (rc) return(rc)
    if (d.moderator.type == "binary")
        rc = pape_solve_binary_shift(d, main_moderator, id.idxz, p0, "main-effect APE (moderator)", bm, msg)
    else
        rc = pape_solve_cont_main(d, main_moderator, id.md_c, id.idxz, pape_var_sd(d.moderator), "main-effect APE (moderator)", bm, msg)
    if (rc) return(rc)
    s_d = (d.focal.type == "binary" ? 1 : pape_var_sd(d.focal))
    s_m = (d.moderator.type == "binary" ? 1 : pape_var_sd(d.moderator))
    e.dgp = d; e.id = id; e.bf = bf; e.bm = bm
    rc = pape_solve_interaction(&pape_f_aie(), e, target, s_d * s_m, "AIE", bi, msg)
    if (rc) return(rc)
    d.beta_focal = bf
    d.beta_mod = bm
    d.beta_int = bi
    d.main_focal = main_focal
    d.main_moderator = main_moderator
    d.estimand = "aie"
    d.target_est = target
    d.p1 = .
    d.has_effect = 1
    return(0)
}

real scalar pape_true_aie(struct pape_dgp scalar d)
{
    if (d.estimand != "aie") pape_err("This DGP has no pinned AIE; run powerape setaie first.")
    if (d.route == "panel") return(pape_true_aie_panel(d))
    return(pape_truth_aie(d, d.beta_focal, d.beta_mod, d.beta_int, pape_integration_draw(d)))
}

// Re-pin at a candidate target, whatever the estimand (mde.R: repin)
real scalar pape_repin(struct pape_dgp scalar d, real scalar target,
    real scalar mf, real scalar mm, string scalar msg)
{
    if (d.has_m) return(pape_set_aie(d, target, mf, mm, msg))
    return(pape_set_ape(d, target, msg))
}

// True coefficient vector in raw design-matrix order (dgp.R: beta_true)
real colvector pape_beta_true(struct pape_dgp scalar d)
{
    real scalar bf, bm, bi, rd, rm
    if (!d.has_effect) pape_err("This DGP has no target effect yet; run powerape setape or setaie first.")
    if (!d.has_m) return((d.beta0 - d.beta_focal * d.focal_ref) \ d.beta_focal \ d.gamma)
    bf = d.beta_focal; bm = d.beta_mod; bi = d.beta_int
    rd = d.focal_ref; rm = d.mod_ref
    return((d.beta0 - bf * rd - bm * rm + bi * rd * rm) \ (bf - bi * rm) \ (bm - bi * rd) \ bi \ d.gamma)
}

// ---------------------------------------------------------------------------
// Printing a DGP (print.R: print.powerape_dgp)
// ---------------------------------------------------------------------------

void pape_print_dgp(struct pape_dgp scalar x)
{
    string scalar labs, sig, note
    real scalar j
    printf("{txt}powerape DGP -- %s%s\n", x.model, (x.route == "empirical" ? " (empirical covariates)" : ""))
    note = ""
    if (x.route == "empirical")
        note = (x.emp_focal_joint ? " [resampled jointly with covariates]" : " [drawn independently of covariates]")
    printf("  focal: %s%s\n", pape_var_label(x.focal), note)
    if (x.has_m)
        printf("  moderator: %s (focal-by-moderator interaction included)\n", pape_var_label(x.moderator))
    if (x.k > 0) {
        if (x.route == "empirical") labs = invtokens(x.covnames', ", ")
        else {
            labs = ""
            for (j = 1; j <= x.k; j++)
                labs = labs + (j > 1 ? ", " : "") + pape_var_label(x.covariates[j])
        }
        printf("  covariates: %s\n", labs)
    }
    else printf("  covariates: none\n")
    if (x.route == "panel")
        printf("  panel: %s periods per unit; rho = %4.2f (latent unit share), cre_share = %4.2f; CRE %s w/ Mundlak means%s, unit-clustered SEs\n",
            pape_panel_len_label(x), x.rho, x.cre_share, x.model,
            (x.unbalanced ? " over observed periods + T-cohort dummies (Wooldridge 2019)" : ""))
    if (x.route == "iv") {
        labs = ""
        for (j = 1; j <= x.n_inst; j++) labs = labs + (j > 1 ? ", " : "") + x.instruments[j].name
        printf("  iv: endogenous focal (%s first stage; instruments: %s); endogeneity rho = %4.2f, first-stage strength = %4.2f; control-function probit, stacked robust SEs\n",
            x.first_stage, labs, x.endogeneity, x.iv_strength)
    }
    if (x.route == "empirical")
        printf("  pilot rows: %f (resampled with replacement)\n", x.n_emp)
    if (x.builder == "from_fit")
        printf("  pilot model: %s (nuisance coefficients lifted; focal coefficient re-solved)\n", x.formula_str)
    if (x.link == "logit" & x.s2 > 0)
        sig = sprintf("nuisance signal = %5.3f (latent index at unit error variance; McKelvey-Zavoina R2 on the logistic scale = %5.3f)",
            x.signal_implied, pape_mz_r2_logit(x.s2))
    else sig = sprintf("nuisance signal (latent pseudo-R2) = %5.3f", x.signal_implied)
    printf("  baseline P(Y=1 | reference) = %5.3f, %s\n", x.baseline, sig)
    if (x.has_effect) {
        if (x.estimand == "aie") {
            printf("  true AIE pinned at %s%7.4f (beta_int = %9.4f)\n", pape_sgn(x.target_est), x.target_est, x.beta_int)
            printf("  main APEs at reference: %s %6.3f (beta = %9.4f), %s %6.3f (beta = %9.4f)\n",
                x.focal.name, x.main_focal, x.beta_focal, x.moderator.name, x.main_moderator, x.beta_mod)
        }
        else if (x.focal.type == "binary")
            printf("  true APE pinned at %s%7.4f (beta_focal = %9.4f, implied P(Y=1 | focal=1) = %5.3f)\n",
                pape_sgn(x.target_est), x.target_est, x.beta_focal, x.p1)
        else
            printf("  true APE pinned at %s%7.4f per unit of %s (beta_focal = %9.4f, reference = %4.2f)\n",
                pape_sgn(x.target_est), x.target_est, x.focal.name, x.beta_focal, x.focal_ref)
    }
    else printf("  target effect not set yet - run powerape setape or powerape setaie.\n")
}

end

mata:
mata set matastrict on
string scalar pape_sgn(real scalar x) return(x >= 0 ? "+" : "")
end

mata:
mata set matastrict on
// Recalibrate an empirical / pilot-model DGP at a new baseline and (unless
// signal is missing) a new nuisance signal, from the stored pilot rows.
// The effect coefficients are cleared (the sweep re-pins them).
struct pape_dgp scalar pape_recal_emp(struct pape_dgp scalar d0, real scalar baseline, real scalar signal)
{
    struct pape_dgp scalar d
    struct pape_env scalar e
    real colvector sds, g_unit, idxz
    real scalar v_unit
    d = d0
    if (baseline == . | baseline <= 0 | baseline >= 1) pape_err("baseline() must lie strictly between 0 and 1.")
    if (signal != .) {
        if (d.gamma_given) pape_err("Cannot vary signal when the DGP was built with explicit gamma().")
        if (signal < 0 | signal >= 1) pape_err("signal() must lie in [0, 1).")
        if (d.k == 0 & signal > 0) pape_err("signal() > 0 requires at least one covariate column.")
        if (d.k == 0 | signal == 0) d.gamma = J(d.k, 1, 0)
        else {
            sds = sqrt(diagonal(variance(d.emp_x)))
            g_unit = 1 :/ sds
            v_unit = variance(d.emp_x * g_unit)
            d.gamma = sqrt((signal / (1 - signal)) / v_unit) :* g_unit
        }
        d.signal = signal
    }
    if (d.k > 0) idxz = d.emp_x * d.gamma
    else idxz = J(d.n_emp, 1, 0)
    d.s2 = (d.k > 0 ? variance(idxz) : 0)
    d.signal_implied = d.s2 / (1 + d.s2)
    e.link = d.link; e.idxz = idxz; e.target = baseline
    d.beta0 = pape_uniroot(&pape_f_b0(), e, -15, 15, 1, 1e-10)
    d.baseline = baseline
    d.beta_focal = .; d.beta_mod = .; d.beta_int = .
    d.has_effect = 0; d.estimand = ""; d.target_est = .; d.p1 = .
    return(d)
}
end
