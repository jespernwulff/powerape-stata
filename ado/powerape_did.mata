*! powerape_did.mata -- difference-in-differences designs with a binary outcome.
*! A treated and a control group observed in n_pre periods before and n_post
*! periods after one treatment date. The target is the average treatment
*! effect on the treated (ATT) after treatment, in probability points; with a
*! binary subgroup, the difference between the two subgroups' ATTs. The
*! simulated world has untreated trends that are parallel on the probit/logit
*! index or in probabilities; each study is analysed by a pooled probit/logit
*! (ATT = average partial effect of the treatment indicator over the treated
*! post-period rows) or by a linear probability model (OLS). Mirrors R/did.R,
*! R/did-subgroup.R, R/did-statement.R, R/did-pilot.R and the OLS part of
*! R/ape-math.R in the development version of the powerape R package.
version 16
mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// Predicates and labels
// ---------------------------------------------------------------------------

real scalar pape_is_did(struct pape_dgp scalar d) return(d.route == "did")

real scalar pape_did_panel(struct pape_dgp scalar d)
    return(d.route == "did" & d.design == "panel")

real scalar pape_did_linear(struct pape_dgp scalar d)
    return(d.route == "did" & d.analysis == "linear")

// an analysis is valid in the simulated world exactly when it assumes the
// scale on which the world's untreated trends are parallel
real scalar pape_did_valid(struct pape_dgp scalar d) return(d.trends == d.analysis)

// designs whose n counts units observed over several periods
real scalar pape_unit_design(struct pape_dgp scalar d)
    return(d.route == "panel" | (d.route == "did" & d.design == "panel"))

// routes whose standard errors are fixed by the design (se() is ignored)
real scalar pape_did_se_fixed(struct pape_dgp scalar d)
    return(d.route == "did" & (d.design == "panel" | d.analysis == "linear"))

string scalar pape_model_label(struct pape_dgp scalar d)
{
    if (pape_did_linear(d)) return("linear probability (OLS)")
    return(d.model)
}

string scalar pape_did_se_label(struct pape_dgp scalar d)
{
    if (d.design == "panel") return(", unit-clustered SEs")
    return(", robust (HC1) SEs")
}

string scalar pape_est_label(string scalar est) return(strupper(est == "" ? "ape" : est))

// the cells that identify the estimand, as named in warnings and notes
string scalar pape_cell_label(string scalar est)
{
    if (est == "aie") return("focal-by-moderator cell")
    if (est == "att") return("group-by-before/after cell")
    if (est == "datt") return("group-by-before/after cell within a subgroup")
    return("focal cell")
}

// code for pape_sparse_warning(): the caller's 0/1 for the routes that
// existed before, 2 (ATT) or 3 (subgroup difference) for these designs
real scalar pape_cell_code(struct pape_dgp scalar d, real scalar old)
{
    if (d.route == "did") return(d.has_sub ? 3 : 2)
    return(old)
}

string scalar pape_did_mismatch_text(struct pape_dgp scalar d)
    return(sprintf("the simulated world has %s-parallel trends and the %s analysis assumes %s-parallel trends",
        d.trends, pape_model_label(d), d.analysis))

// A search needs an analysis that is valid in the simulated world.
void pape_did_need_valid(struct pape_dgp scalar d, string scalar what)
{
    if (d.route == "did" & !pape_did_valid(d))
        pape_err(sprintf("%s needs an analysis that is valid in the simulated world. Here %s, so the analysis is biased for the effect on the treated: powerape power and powerape curve report its rejection rate and the mean of its estimates. Search with trends() and analysis() on the same scale (both index or both linear).",
            what, pape_did_mismatch_text(d)))
}

// A result whose analysis assumes another scale than its world is a
// rejection rate: the commands that re-read or report power stop on it.
void pape_did_power_only(struct pape_dgp scalar d, string scalar what)
{
    if (d.route == "did" & !pape_did_valid(d))
        pape_err(sprintf("%s needs power, and this result is a rejection rate: %s, so the analysis is biased for the effect on the treated. Use a design whose trends() and analysis() are on the same scale.",
            what, pape_did_mismatch_text(d)))
}

// last focal-side design column (group, period dummies, treatment indicator
// and their subgroup counterparts): the range of the Hauck-Donner backup
real scalar pape_did_last_col(struct pape_dgp scalar d)
    return(d.has_sub ? 2 * d.n_periods + 4 : d.n_periods + 2)

// ---------------------------------------------------------------------------
// Calibration
// ---------------------------------------------------------------------------

// deterministic draw of (group, subgroup, covariates) under the design's seed
void pape_did_reg_draw(struct pape_dgp scalar d, real scalar n, real scalar seed,
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
    if (d.has_sub) md = pape_transform_u(U[., 2], d.subgroup)
    else md = J(n, 1, 0)
    zoff = 1 + d.has_sub
    if (d.k > 0) {
        Z = U[., (zoff + 1)..(zoff + d.k)]
        for (j = 1; j <= d.k; j++) Z[., j] = pape_transform_u(Z[., j], d.covariates[j])
    }
    else Z = J(n, 0, .)
    if (seed != .) rngstate(st)
}

// nuisance-index draws by group (first digit) and subgroup (second digit);
// a single 0 per cell without covariates
void pape_did_qsplit(struct pape_dgp scalar d, real colvector q00, real colvector q01,
    real colvector q10, real colvector q11)
{
    real colvector xd, md, q
    real matrix Z
    if (d.k == 0) {
        q00 = 0; q01 = 0; q10 = 0; q11 = 0
        return
    }
    pape_did_reg_draw(d, d.n_int, d.seed_int, xd, md, Z)
    q = Z * d.gamma
    q00 = select(q, (xd :== 0) :& (md :== 0))
    q01 = select(q, (xd :== 0) :& (md :== 1))
    q10 = select(q, (xd :== 1) :& (md :== 0))
    q11 = select(q, (xd :== 1) :& (md :== 1))
}

// solve c in mean(G(offset + c + q)) = rate (closed form without covariates)
real scalar pape_did_solve(string scalar link, real colvector q, real scalar rate,
    real scalar offset)
{
    struct pape_env scalar e
    if (rows(q) == 1) return(pape_Ginv(link, rate) - offset - q[1])
    e.link = link; e.idxz = offset :+ q; e.target = rate
    return(pape_uniroot(&pape_f_b0(), e, -15, 15, 1, 1e-10))
}

// one number (the same in both subgroups) or two
real rowvector pape_did_two(real rowvector x, real scalar M, string scalar nm, string scalar sname)
{
    if (cols(x) < 1 | cols(x) > M | hasmissing(x)) {
        if (M == 1) pape_err(sprintf("%s() must be one number (two numbers are for designs with subgroup()).", nm))
        pape_err(sprintf("%s() must be one number (the same in both subgroups) or two: the reference subgroup (%s = 0), then %s = 1.", nm, sname, sname))
    }
    if (cols(x) == M) return(x)
    return(J(1, M, x[1]))
}

struct pape_dgp scalar pape_dgp_did(string scalar name, string scalar model,
    real rowvector baseline_in, real rowvector gap_in, real rowvector trend_in,
    real scalar p_treated, string scalar design, real scalar rho, real scalar p_post,
    string scalar trends, string scalar analysis, real scalar n_pre, real scalar n_post,
    real scalar has_sub, struct pape_var scalar subgroup, struct pape_var colvector covariates,
    real scalar corr_given, real scalar corr_scalar, real matrix corr_matrix,
    real scalar signal, real scalar gamma_given, real colvector gamma,
    real scalar n_int, real scalar seed_int)
{
    struct pape_dgp scalar d
    real rowvector baseline, gap, trend, a, b, l, th, cf
    real colvector xd, md, g_unit, q, q0, q1
    real matrix Z
    real scalar M, m, K, j, v_unit, eps, x
    string scalar lab
    d = pape_dgp_blank()
    d.name = name
    if (model != "probit" & model != "logit") pape_err("model() must be probit or logit.")
    if (design != "panel" & design != "cross-sections")
        pape_err("design() must be panel or cross-sections.")
    if (trends != "index" & trends != "linear") pape_err("trends() must be index or linear.")
    if (analysis != "index" & analysis != "linear") pape_err("analysis() must be index or linear.")
    d.model = model; d.link = model
    d.route = "did"; d.builder = "did"
    d.design = design; d.trends = trends; d.analysis = analysis
    d.has_sub = has_sub
    M = (has_sub ? 2 : 1)
    if (has_sub) {
        if (subgroup.type != "binary")
            pape_err("subgroup() must describe a binary variable, e.g. subgroup(female, p(0.4)).")
        d.subgroup = subgroup
    }
    baseline = pape_did_two(baseline_in, M, "baseline", subgroup.name)
    gap = pape_did_two((cols(gap_in) == 0 ? 0 : gap_in), M, "gap", subgroup.name)
    trend = pape_did_two((cols(trend_in) == 0 ? 0 : trend_in), M, "trend", subgroup.name)
    if (any(baseline :<= 0) | any(baseline :>= 1)) pape_err("baseline() must lie strictly between 0 and 1.")
    if (p_treated == . | p_treated <= 0 | p_treated >= 1) pape_err("ptreated() must lie strictly between 0 and 1.")
    if (rho == . | rho < 0 | rho >= 1) pape_err("rho() must lie in [0, 1).")
    if (n_pre == . | n_pre < 1 | n_pre != round(n_pre) | n_post == . | n_post < 1 | n_post != round(n_post))
        pape_err("npre() and npost() must be whole numbers of at least 1.")
    if (n_int == . | n_int < 1e4) pape_err("nint() must be at least 10000.")
    d.p_post_given = (p_post != .)
    if (p_post == .) p_post = n_post / (n_pre + n_post)
    if (p_post <= 0 | p_post >= 1) pape_err("ppost() must lie strictly between 0 and 1.")
    eps = 1e-4
    for (m = 1; m <= M; m++) {
        lab = (has_sub ? sprintf("Subgroup %s = %f: ", subgroup.name, m - 1) : "")
        x = baseline[m] + gap[m]
        if (!(x > eps & x < 1 - eps))
            pape_err(sprintf("%sbaseline() + gap() = %6.4f is not a rate strictly inside (0, 1): the treated group's pre-period rate must be.", lab, x))
        x = baseline[m] + trend[m]
        if (!(x > eps & x < 1 - eps))
            pape_err(sprintf("%sbaseline() + trend() = %6.4f is not a rate strictly inside (0, 1): the control group's post-period rate must be.", lab, x))
        x = baseline[m] + gap[m] + trend[m]
        if (trends == "linear" & !(x > eps & x < 1 - eps))
            pape_err(sprintf("%swith trends(linear) the treated group's untreated post-period rate is baseline + gap + trend = %6.4f, outside (0, 1). Linear parallel trends cannot hold for these rates (a bounded outcome cannot follow the control group's change in probability points); use trends(index) or other rates.", lab, x))
    }
    if (design == "cross-sections" & rho != 0)
        printf("{err}warning: rho() is ignored for design(cross-sections) (independent observations).\n")
    d.k = rows(covariates)
    if (d.k > 0) d.covariates = covariates
    d.focal = pape_var_new("treated group", "binary", p_treated, ., ., .)
    d.has_m = 0
    K = 1 + has_sub + d.k
    if (corr_given == 2 & (rows(corr_matrix) != K | cols(corr_matrix) != K))
        pape_err(sprintf("correlation() must be a %f x %f matrix over (group%s, covariates).", K, K, (has_sub ? ", subgroup" : "")))
    d.R = pape_build_corr(K, corr_given, corr_scalar, corr_matrix)
    d.Rchol = cholesky(d.R)'
    d.corr_given = corr_given; d.corr_scalar = corr_scalar; d.corr_matrix = corr_matrix
    d.n_int = n_int; d.seed_int = seed_int
    d.gamma_given = gamma_given
    if (!gamma_given) {
        if (signal == . | signal < 0 | signal >= 1) pape_err("signal() must lie in [0, 1).")
        if (d.k == 0 & signal > 0) pape_err("signal() > 0 requires at least one covariate.")
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
        q = 0; xd = 0; md = 0
    }
    else {
        pape_did_reg_draw(d, n_int, seed_int, xd, md, Z)
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
        q = Z * d.gamma
        d.s2 = variance(q)
    }
    d.signal_implied = d.s2 / (1 + d.s2)
    // calibration within each subgroup, on its own covariate distributions
    a = J(1, M, .); b = J(1, M, .); l = J(1, M, .); th = J(1, M, .); cf = J(1, M, .)
    for (m = 1; m <= M; m++) {
        if (d.k == 0) {
            q0 = 0; q1 = 0
        }
        else {
            q0 = select(q, (xd :== 0) :& (md :== m - 1))
            q1 = select(q, (xd :== 1) :& (md :== m - 1))
            if (has_sub & (rows(q0) < 100 | rows(q1) < 100))
                pape_err("A group-by-subgroup cell is (nearly) empty in the integration draw; check the shares and the correlation.")
        }
        a[m] = pape_did_solve(model, q0, baseline[m], 0)
        b[m] = pape_did_solve(model, q1, baseline[m] + gap[m], a[m])
        l[m] = pape_did_solve(model, q0, baseline[m] + trend[m], a[m])
        if (trends == "index") {
            th[m] = 0
            cf[m] = mean(pape_G(model, a[m] + b[m] + l[m] :+ q1))
        }
        else {
            cf[m] = baseline[m] + gap[m] + trend[m]
            th[m] = pape_did_solve(model, q1, cf[m], a[m] + b[m] + l[m])
        }
    }
    d.baseline = baseline[1]
    d.beta0 = a[1]
    d.did_base = baseline; d.did_gap = gap; d.did_trend = trend
    d.did_a = a; d.did_b = b; d.did_l = l; d.did_th = th; d.did_de = J(1, M, .)
    d.c00 = baseline; d.c10 = baseline :+ gap; d.c01 = baseline :+ trend
    d.cf = cf; d.c11 = J(1, M, .)
    d.p_treated = p_treated; d.p_post = p_post
    d.rho = (design == "panel" ? rho : 0)
    d.n_pre = n_pre; d.n_post = n_post
    d.n_periods = n_pre + n_post
    d.unbalanced = 0
    d.focal_ref = 0
    return(d)
}

// Pin the effect on the treated: one target, or the two subgroups' effects
// (reference first; the estimand is then their difference).
real scalar pape_set_att(struct pape_dgp scalar d, real rowvector target, string scalar msg)
{
    real colvector q00, q01, q10, q11
    real rowvector de
    real scalar M, m, cf, eps
    msg = ""
    if (d.route != "did") {
        msg = "powerape setatt pins the effect on the treated of a difference-in-differences design; build one with powerape dgp NAME, did. For other designs use powerape setape or powerape setaie."
        return(1)
    }
    M = (d.has_sub ? 2 : 1)
    if (cols(target) != M | hasmissing(target)) {
        if (d.has_sub) msg = sprintf("This design has a subgroup: give target() the two effects on the treated, the effect for %s = 0 and then the effect for %s = 1. The estimand is their difference.", d.subgroup.name, d.subgroup.name)
        else msg = "target() must be one number (the effect on the treated in probability points)."
        return(1)
    }
    eps = 1e-4
    for (m = 1; m <= M; m++) {
        cf = d.cf[m]
        if (target[m] <= -cf + eps | target[m] >= (1 - cf) - eps) {
            if (d.has_sub) msg = sprintf("Target ATT %6.4f is infeasible for subgroup %s = %f: without treatment its treated group's post-period rate is %5.3f, so the attainable range is (%6.3f, %5.3f).",
                target[m], d.subgroup.name, m - 1, cf, -cf, 1 - cf)
            else msg = sprintf("Target ATT %6.4f is infeasible for this design: without treatment the treated group's post-period rate is %5.3f, so the attainable range is (%6.3f, %5.3f). Adjust the target, the rates, or trends().",
                target[m], cf, -cf, 1 - cf)
            return(1)
        }
    }
    pape_did_qsplit(d, q00, q01, q10, q11)
    de = J(1, M, .)
    for (m = 1; m <= M; m++)
        de[m] = pape_did_solve(d.link, (m == 1 ? q10 : q11), d.cf[m] + target[m],
            d.did_a[m] + d.did_b[m] + d.did_l[m] + d.did_th[m])
    d.did_de = de
    d.c11 = d.cf :+ target
    if (d.has_sub) {
        d.att_sub = target
        d.estimand = "datt"
        d.target_est = target[2] - target[1]
        d.p1 = .
    }
    else {
        d.estimand = "att"
        d.target_est = target[1]
        d.p1 = d.cf[1] + target[1]
    }
    d.has_effect = 1
    return(0)
}

// the true effect(s) recomputed on the integration set: (ATT), or
// (reference, subgroup, difference)
real rowvector pape_true_att(struct pape_dgp scalar d)
{
    real colvector q00, q01, q10, q11, q1
    real rowvector e
    real scalar M, m, off
    if (d.route != "did" | !d.has_effect)
        pape_err("This DGP has no pinned ATT; build it with powerape dgp NAME, did and run powerape setatt.")
    M = (d.has_sub ? 2 : 1)
    pape_did_qsplit(d, q00, q01, q10, q11)
    e = J(1, M, .)
    for (m = 1; m <= M; m++) {
        q1 = (m == 1 ? q10 : q11)
        off = d.did_a[m] + d.did_b[m] + d.did_l[m] + d.did_th[m]
        e[m] = mean(pape_G(d.link, off + d.did_de[m] :+ q1) :- pape_G(d.link, off :+ q1))
    }
    if (d.has_sub) return((e, e[2] - e[1]))
    return(e)
}

// What the design's observable block rates are read as under each
// parallel-trends scale, by subgroup: the cross difference of probabilities
// and c11 - G(G^-1(c10) + G^-1(c01) - G^-1(c00)). Marginal rates.
void pape_did_readings(struct pape_dgp scalar d, real rowvector lin, real rowvector idx)
{
    real scalar M, m
    M = (d.has_sub ? 2 : 1)
    lin = J(1, M, .); idx = J(1, M, .)
    if (!d.has_effect) return
    for (m = 1; m <= M; m++) {
        lin[m] = d.c11[m] - (d.c10[m] + d.c01[m] - d.c00[m])
        idx[m] = d.c11[m] - pape_G(d.link, pape_Ginv(d.link, d.c10[m]) + pape_Ginv(d.link, d.c01[m]) - pape_Ginv(d.link, d.c00[m]))
    }
}

// True coefficients in design-matrix order. Without a subgroup:
// [1, D, f_2..f_T, W, Z] (0 on the remaining pre-period dummies, l on every
// post-period dummy, th + de on W). With a subgroup: the reference
// subgroup's block, the other subgroup's differences from it, then Z.
real colvector pape_did_blk(struct pape_dgp scalar d, real scalar m)
    return(d.did_a[m] \ d.did_b[m] \ J(d.n_pre - 1, 1, 0) \ J(d.n_post, 1, d.did_l[m]) \ (d.did_th[m] + d.did_de[m]))

real colvector pape_did_beta(struct pape_dgp scalar d)
{
    real colvector b1
    if (!d.has_effect) pape_err("This DGP has no target effect yet; run powerape setatt first.")
    b1 = pape_did_blk(d, 1)
    if (!d.has_sub) return(b1 \ d.gamma)
    return(b1 \ (pape_did_blk(d, 2) - b1) \ d.gamma)
}

// ---------------------------------------------------------------------------
// One simulated study
// ---------------------------------------------------------------------------

// Panel: n units x T rows (unit-major, periods in order), outcome rates
// conditional on the unit draw of a Gaussian copula, so that the marginal
// rates stay exact under either link. Cross-sections: n independent rows,
// after-block membership Bernoulli(p_post) and, with more than two periods,
// a uniform period within the block.
struct pape_xdraw scalar pape_draw_x_did(struct pape_dgp scalar d, real scalar n)
{
    struct pape_xdraw scalar xx
    real matrix U, Z, Fm, base
    real colvector bt, dd, mm, id, D, Mv, period, post, p, a_i, pp, pq
    real scalar j, zoff, Tn, n_pre, N
    bt = pape_did_beta(d)
    U = pape_draw_std_mvn(n, d.Rchol)
    dd = pape_transform_u(U[., 1], d.focal)
    if (d.has_sub) mm = pape_transform_u(U[., 2], d.subgroup)
    zoff = 1 + d.has_sub
    if (d.k > 0) {
        Z = U[., (zoff + 1)..(zoff + d.k)]
        for (j = 1; j <= d.k; j++) Z[., j] = pape_transform_u(Z[., j], d.covariates[j])
    }
    n_pre = d.n_pre
    Tn = d.n_periods
    if (d.design == "panel") {
        id = (1::n) # J(Tn, 1, 1)
        D = dd[id]
        if (d.has_sub) Mv = mm[id]
        period = J(n, 1, 1) # (1::Tn)
        post = (period :> n_pre)
        if (d.k > 0) Z = Z[id, .]
        xx.id = id
        xx.has_id = 1
    }
    else {
        D = dd
        if (d.has_sub) Mv = mm
        post = (runiform(n, 1) :< d.p_post)
        if (Tn == 2) period = 1 :+ post
        else {
            pp = n_pre :+ ceil(runiform(n, 1) :* d.n_post)
            pq = ceil(runiform(n, 1) :* n_pre)
            period = post :* pp :+ (1 :- post) :* pq
        }
        xx.id = J(0, 1, .)
        xx.has_id = 0
    }
    N = rows(D)
    if (Tn == 2) Fm = post
    else {
        Fm = J(N, Tn - 1, 0)
        for (j = 2; j <= Tn; j++) Fm[., j - 1] = (period :== j)
    }
    base = (J(N, 1, 1), D, Fm, D :* post)
    if (d.has_sub) xx.X = (base, Mv :* base)
    else xx.X = base
    if (d.k > 0) xx.X = (xx.X, Z)
    p = pape_G(d.link, xx.X * bt)
    // keep the rates strictly inside (0, 1): invnormal() is missing at 0 and 1
    p = rowmax((rowmin((p, J(N, 1, 1 - 1e-15))), J(N, 1, 1e-15)))
    if (xx.has_id & d.rho > 0) {
        a_i = rnormal(n, 1, 0, 1)
        xx.pr = normal((invnormal(p) :- sqrt(d.rho) :* a_i[id]) :/ sqrt(1 - d.rho))
    }
    else xx.pr = p
    xx.w_col = Tn + 2
    xx.focal_col = (d.has_sub ? 2 * Tn + 4 : Tn + 2)
    xx.start = bt
    xx.has_start = 1
    xx.d = D
    xx.post = post
    if (d.has_sub) xx.sub = Mv
    else xx.sub = J(0, 1, .)
    return(xx)
}

// ---------------------------------------------------------------------------
// Estimators
// ---------------------------------------------------------------------------

// mean over the selected rows of G(x'b) - G(x'b with the columns `off` at 0)
// and its gradient; rc 0 when no row is selected
real scalar pape_did_part(real colvector beta, real matrix X, real colvector sel,
    real rowvector off, string scalar link, real scalar ape, real colvector jac)
{
    real matrix X1, X0
    real colvector a1, a0, idx
    idx = selectindex(sel)
    if (rows(idx) == 0 | cols(idx) == 0) return(0)
    X1 = X[idx, .]
    X0 = X1
    X0[., off] = J(rows(X1), cols(off), 0)
    a1 = X1 * beta; a0 = X0 * beta
    ape = mean(pape_G(link, a1) :- pape_G(link, a0))
    jac = (mean(pape_g(link, a1) :* X1) :- mean(pape_g(link, a0) :* X0))'
    return(1)
}

// ATT: average partial effect of the treatment indicator (column w_col)
// over the rows where it is 1. Returns 0 without a treated post-period row.
real scalar pape_att_est(real colvector beta, real matrix X, string scalar link,
    real scalar w_col, struct pape_est scalar e)
{
    real scalar ape
    real colvector jac
    if (!pape_did_part(beta, X, X[., w_col] :== 1, w_col, link, ape, jac)) return(0)
    e.ape = ape
    e.jac = jac
    return(1)
}

// Difference between the two subgroups' ATTs and its gradient. Subgroup 1:
// its treated post-period rows, with the treatment indicator and its
// subgroup product switched off; subgroup 0: its treated post-period rows,
// with the treatment indicator switched off.
real scalar pape_datt_est(real colvector beta, real matrix X, real colvector sub,
    real scalar w_col, real scalar mw_col, string scalar link, struct pape_est scalar e,
    real rowvector by)
{
    real scalar a1, a0
    real colvector j1, j0
    if (!pape_did_part(beta, X, X[., mw_col] :== 1, (w_col, mw_col), link, a1, j1)) return(0)
    if (!pape_did_part(beta, X, (X[., w_col] :== 1) :& (sub :== 0), w_col, link, a0, j0)) return(0)
    e.ape = a1 - a0
    e.jac = j1 - j0
    by = (a0, a1)
    return(1)
}

real scalar pape_did_est(real colvector beta, struct pape_xdraw scalar xx,
    string scalar link, struct pape_est scalar e)
{
    real rowvector by
    if (rows(xx.sub) == 0) return(pape_att_est(beta, xx.X, link, xx.focal_col, e))
    return(pape_datt_est(beta, xx.X, xx.sub, xx.w_col, xx.focal_col, link, e, by))
}

// identifying cells: group x before/after, within each subgroup if any;
// returns (sep, min_cell)
real rowvector pape_did_cell_check(real colvector y, struct pape_xdraw scalar xx)
{
    real colvector s0, s1
    real rowvector c0, c1
    if (rows(xx.sub) == 0) return(pape_cell_check(y, xx.d, xx.post))
    s0 = selectindex(xx.sub :== 0)
    s1 = selectindex(xx.sub :== 1)
    if (rows(s0) == 0 | cols(s0) == 0 | rows(s1) == 0 | cols(s1) == 0) return((0, .))
    c0 = pape_cell_check(y[s0], xx.d[s0], xx.post[s0])
    c1 = pape_cell_check(y[s1], xx.d[s1], xx.post[s1])
    return(((c0[1] | c1[1]), min((c0[2], c1[2]))))
}

// ---------------------------------------------------------------------------
// Ordinary least squares for the linear probability analysis
// ---------------------------------------------------------------------------

// rc 0 for a rank-deficient design (e.g. an empty group-by-period cell)
real scalar pape_fit_ols(real matrix X, real colvector y, real colvector b,
    real colvector resid, real matrix XtXi)
{
    if (rows(X) <= cols(X)) return(0)
    XtXi = invsym(quadcross(X, X))
    if (any(diagonal(XtXi) :== 0)) return(0)
    b = XtXi * quadcross(X, y)
    if (hasmissing(b)) return(0)
    resid = y - X * b
    return(1)
}

// HC1: the HC0 sandwich times N / (N - K), as Stata's regress, vce(robust)
real matrix pape_vcov_ols_hc1(real matrix X, real colvector resid, real matrix XtXi)
{
    real scalar n, k
    real matrix S
    n = rows(X); k = cols(X)
    S = X :* resid
    return((n / (n - k)) :* (XtXi * quadcross(S, S) * XtXi))
}

// cluster-robust with the factor G / (G - 1) x (N - 1) / (N - K) of Stata's
// regress, vce(cluster)
real matrix pape_vcov_ols_cluster(real matrix X, real colvector resid, real matrix XtXi,
    real colvector id)
{
    real scalar n, k, Gn
    real matrix S, U, info
    real colvector ord
    n = rows(X); k = cols(X)
    S = X :* resid
    if (n > 1 & any(id[2::n] :< id[1::(n - 1)])) {
        ord = order(id, 1)
        info = panelsetup(id[ord], 1)
        U = panelsum(S[ord, .], info)
    }
    else {
        info = panelsetup(id, 1)
        U = panelsum(S, info)
    }
    Gn = rows(U)
    return((Gn / (Gn - 1)) * ((n - 1) / (n - k)) :* (XtXi * quadcross(U, U) * XtXi))
}

// regress by Stata: coefficient and robust / cluster standard error of
// design column lc. X includes the constant in column 1. rc 0 when the
// study has no usable estimate (a column omitted, observations dropped).
real scalar pape_ols_stata(real matrix X, real colvector y, real colvector id,
    string scalar vce, real scalar lc, real scalar b_out, real scalar se_out)
{
    real scalar n, k, rc, j, i
    real rowvector bs
    real matrix Vs
    string colvector names
    string scalar cmd
    n = rows(X); k = cols(X) - 1
    st_store(., "pape_y", y)
    st_store(., tokens(pape_xnames(k)), X[., 2..(k + 1)])
    if (vce == "cluster") st_store(., "pape_id", id)
    cmd = "quietly regress pape_y" + pape_xnames(k) + (vce == "cluster" ? ", vce(cluster pape_id)" : ", vce(robust)")
    rc = _stata(cmd, 1)
    if (rc) return(0)
    if (st_numscalar("e(N)") != n) return(0)
    bs = st_matrix("e(b)")
    Vs = st_matrix("e(V)")
    if (cols(bs) != k + 1) return(0)
    names = st_matrixcolstripe("e(b)")[., 2]
    for (j = 1; j <= rows(names); j++) {
        if (substr(names[j], 1, 2) == "o.") return(0)
    }
    i = lc - 1                       // e(b) holds the constant last
    b_out = bs[i]
    se_out = sqrt(Vs[i, i])
    return(b_out != . & se_out != .)
}

// Effect on the treated by Stata's margins after probit/logit with the
// treatment indicator as a factor variable:
//   margins, dydx(w) subpop(if w == 1)                 (one effect)
//   margins, dydx(w) subpop(if w == 1) over(m)         (subgroup difference)
// Stored columns: pape_x1 = group, pape_x2..x_T period dummies, pape_x(T+1)
// = W, and with a subgroup pape_x(T+2) = M, then M x (group, periods), with
// M x W supplied by the factor interaction. rc 0 when not estimable.
real scalar pape_did_margins(struct pape_dgp scalar d, struct pape_xdraw scalar xx,
    real colvector y, string scalar vce, real scalar asis, real scalar ape_out,
    real scalar se_out)
{
    real scalar n, k, rc, j, Tn, wj, mj, kk
    real rowvector b
    real matrix V
    string scalar cmd, vars
    n = rows(xx.X); k = cols(xx.X) - 1
    Tn = d.n_periods
    wj = Tn + 1
    st_store(., "pape_y", y)
    st_store(., tokens(pape_xnames(k)), xx.X[., 2..(k + 1)])
    if (vce == "cluster") st_store(., "pape_id", xx.id)
    vars = ""
    for (j = 1; j <= Tn; j++) vars = vars + sprintf(" pape_x%f", j)
    if (d.has_sub) {
        mj = Tn + 2
        vars = vars + sprintf(" i.pape_x%f##i.pape_x%f", wj, mj)
        for (j = Tn + 3; j <= 2 * Tn + 2; j++) vars = vars + sprintf(" pape_x%f", j)
        for (j = 2 * Tn + 4; j <= k; j++) vars = vars + sprintf(" pape_x%f", j)
    }
    else {
        vars = vars + sprintf(" i.pape_x%f", wj)
        for (j = Tn + 2; j <= k; j++) vars = vars + sprintf(" pape_x%f", j)
    }
    cmd = "quietly " + d.link + " pape_y" + vars + ", iterate(100)"
    if (vce == "robust") cmd = cmd + " vce(robust)"
    else if (vce == "cluster") cmd = cmd + " vce(cluster pape_id)"
    if (asis) cmd = cmd + " asis"
    rc = _stata(cmd, 1)
    if (rc) return(0)
    if (st_numscalar("e(converged)") != 1 | st_numscalar("e(N)") != n) return(0)
    cmd = sprintf("quietly margins, dydx(pape_x%f) subpop(if pape_x%f == 1)", wj, wj)
    if (d.has_sub) cmd = cmd + sprintf(" over(pape_x%f)", mj)
    rc = _stata(cmd, 1)
    if (rc) return(0)
    b = st_matrix("r(b)")
    V = st_matrix("r(V)")
    kk = cols(b)
    if (kk == 0) return(0)
    if (d.has_sub) {
        // the last two columns are the effect at subgroup 0 and at subgroup 1
        if (kk < 2) return(0)
        ape_out = b[kk] - b[kk - 1]
        se_out = sqrt(V[kk, kk] + V[kk - 1, kk - 1] - 2 * V[kk, kk - 1])
    }
    else {
        ape_out = b[kk]
        se_out = sqrt(V[kk, kk])
    }
    if (ape_out == . | se_out == . | se_out <= 0) return(0)
    return(1)
}

// who fitted the simulated studies of a difference-in-differences design
string scalar pape_did_fit_lab(struct pape_dgp scalar d, string scalar se,
    string scalar separation, string scalar engine, string scalar effects)
{
    string scalar cmd, opts
    if (d.analysis == "linear") {
        if (engine == "mata") return("Mata engine (OLS)")
        return(sprintf("regress, vce(%s) via engine(stata)", (d.design == "panel" ? "cluster" : "robust")))
    }
    if (engine == "mata") return("Mata engine (Fisher scoring)")
    cmd = d.link
    opts = ""
    if (d.design == "panel") opts = " vce(cluster)"
    else if (se == "robust") opts = " vce(robust)"
    if (separation == "keep") opts = opts + " asis"
    if (opts != "") cmd = cmd + "," + opts
    return(sprintf("%s via engine(%s), effects(%s)", cmd, engine, effects))
}

// ---------------------------------------------------------------------------
// One given study, analysed as the simulation loop analyses a simulated one
// (no separation rules). X is the design matrix in the route's layout with
// the constant in column 1; id and sub are 0 x 1 without clustering / without
// a subgroup. Returns 0 when there is no usable estimate. b and jac are
// filled for the closed-form paths, by for subgroup designs. Used by the
// fixed-data tier of the validation.
// ---------------------------------------------------------------------------
real scalar pape_did_analyze(real matrix X, real colvector y, real colvector id,
    real colvector sub, string scalar link, real scalar n_periods, string scalar analysis,
    string scalar se_type, string scalar engine, string scalar effects,
    real scalar est_out, real scalar se_out, real colvector b, real colvector jac,
    real rowvector by)
{
    struct pape_dgp scalar d
    struct pape_xdraw scalar xx
    struct pape_fit scalar ft
    struct pape_est scalar est
    real colvector resid
    real matrix V, XtXi
    real scalar panel, ok, lc
    string scalar vce
    panel = (rows(id) > 0)
    d = pape_dgp_blank()
    d.route = "did"; d.link = link; d.model = link
    d.n_periods = n_periods
    d.has_sub = (rows(sub) > 0)
    d.design = (panel ? "panel" : "cross-sections")
    d.analysis = analysis
    xx.X = X
    xx.id = id
    xx.has_id = panel
    xx.sub = sub
    xx.w_col = n_periods + 2
    xx.focal_col = (d.has_sub ? 2 * n_periods + 4 : n_periods + 2)
    lc = xx.focal_col
    b = J(0, 1, .); jac = J(0, 1, .); by = J(1, 0, .)
    est_out = .; se_out = .
    if (engine == "stata") pape_frame_open(rows(y), cols(X) - 1, 0)
    ok = 0
    if (analysis == "linear") {
        if (engine == "stata") ok = pape_ols_stata(X, y, id, (panel ? "cluster" : "robust"), lc, est_out, se_out)
        else if (pape_fit_ols(X, y, b, resid, XtXi)) {
            if (panel) V = pape_vcov_ols_cluster(X, resid, XtXi, id)
            else V = pape_vcov_ols_hc1(X, resid, XtXi)
            est_out = b[lc]
            se_out = sqrt(V[lc, lc])
            ok = 1
        }
    }
    else {
        vce = (panel ? "cluster" : (se_type == "robust" ? "robust" : "oim"))
        if (engine == "stata" & effects == "margins") ok = pape_did_margins(d, xx, y, vce, 0, est_out, se_out)
        else {
            if (engine == "stata") ft = pape_fit_stata(X, y, id, link, vce, 0)
            else ft = pape_fit_index(X, y, link, J(0, 1, .))
            if (ft.ok) {
                if (d.has_sub) ok = pape_datt_est(ft.b, X, sub, xx.w_col, xx.focal_col, link, est, by)
                else ok = pape_att_est(ft.b, X, link, xx.focal_col, est)
                if (ok) {
                    if (engine == "stata") V = ft.V
                    else if (panel) V = pape_vcov_cluster(ft, X, y, id, link)
                    else if (se_type == "robust") V = pape_vcov_sandwich(ft, X, y, link)
                    else V = ft.V
                    b = ft.b
                    jac = est.jac
                    est_out = est.ape
                    se_out = sqrt(est.jac' * V * est.jac)
                }
            }
        }
    }
    if (engine == "stata") pape_frame_close()
    return(ok)
}

// ---------------------------------------------------------------------------
// Core loop for difference-in-differences designs (engine.R: sim_ci)
// ---------------------------------------------------------------------------

struct pape_sim scalar pape_sim_ci_did(struct pape_dgp scalar d, real scalar n,
    real scalar nsim, real scalar conf, real scalar seed, string scalar se_type,
    string scalar separation)
{
    struct pape_sim scalar s
    struct pape_xdraw scalar xx
    struct pape_fit scalar ft
    struct pape_est scalar est
    real colvector y, b, resid
    real rowvector cc, hd_cols
    real matrix V, XtXi
    real scalar z, lin, panel, sep_fail, r, s_, ape_v, lc, closed
    real scalar use_stata, use_margins, frame_open
    string scalar st, vce
    z = pape_zcrit(conf)
    lin = (d.analysis == "linear")
    panel = (d.design == "panel")
    sep_fail = (separation != "keep")
    hd_cols = (2..pape_did_last_col(d))
    use_stata = (pape_engine() == "stata")
    use_margins = (use_stata & pape_effects() == "margins")
    frame_open = 0
    (void) pape_did_beta(d)          // stops when the effect is not pinned
    if (seed != .) {
        st = rngstate()
        rseed(seed)
    }
    s.est = J(nsim, 1, .); s.l = J(nsim, 1, .); s.u = J(nsim, 1, .)
    s.se = J(nsim, 1, .); s.min_cell = J(nsim, 1, .)
    s.ok = J(nsim, 1, 0); s.sep = J(nsim, 1, 0)
    for (r = 1; r <= nsim; r++) {
        xx = pape_draw_x_did(d, n)
        y = runiform(rows(xx.pr), 1) :< xx.pr
        cc = pape_did_cell_check(y, xx)
        s.min_cell[r] = cc[2]
        if (all(y :== y[1])) continue
        if (use_stata) {
            if (!frame_open) {
                pape_frame_open(rows(y), cols(xx.X) - 1, 0)
                frame_open = 1
            }
            else pape_frame_resize(rows(y))
        }
        if (lin) {
            // linear probability analysis: no separation, one coefficient,
            // always robust (HC1 for independent rows, unit-clustered in panels)
            lc = xx.focal_col
            if (use_stata) {
                if (!pape_ols_stata(xx.X, y, xx.id, (panel ? "cluster" : "robust"), lc, ape_v, s_)) continue
            }
            else {
                if (!pape_fit_ols(xx.X, y, b, resid, XtXi)) continue
                if (panel) V = pape_vcov_ols_cluster(xx.X, resid, XtXi, xx.id)
                else V = pape_vcov_ols_hc1(xx.X, resid, XtXi)
                ape_v = b[lc]
                s_ = sqrt(max((0, V[lc, lc])))
            }
        }
        else {
            if (cc[1]) {
                s.sep[r] = 1
                if (sep_fail) continue
            }
            closed = 1
            if (use_stata) {
                vce = (panel ? "cluster" : (se_type == "robust" ? "robust" : "oim"))
                if (use_margins) {
                    if (!pape_did_margins(d, xx, y, vce, !sep_fail, ape_v, s_)) continue
                    closed = 0
                }
                else {
                    ft = pape_fit_stata(xx.X, y, xx.id, d.link, vce, !sep_fail)
                    if (!ft.ok) continue
                }
            }
            else {
                ft = pape_fit_index(xx.X, y, d.link, xx.start)
                if (!ft.ok) continue
            }
            if (closed) {
                if (!s.sep[r] & (ft.rank < cols(xx.X) | pape_hauck_donner(ft, xx.X, hd_cols))) {
                    s.sep[r] = 1
                    if (sep_fail) continue
                }
                if (!pape_did_est(ft.b, xx, d.link, est)) continue
                if (use_stata) V = ft.V
                else if (panel) V = pape_vcov_cluster(ft, xx.X, y, xx.id, d.link)
                else if (se_type == "robust") V = pape_vcov_sandwich(ft, xx.X, y, d.link)
                else V = ft.V
                s_ = sqrt(max((0, est.jac' * V * est.jac)))
                ape_v = est.ape
            }
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

end
