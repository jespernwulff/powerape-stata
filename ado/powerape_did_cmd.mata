*! powerape_did_cmd.mata -- command layer for difference-in-differences
*! designs: the design print, setatt / trueatt, the robustness sweep, the
*! methods paragraph, and the design built from pre-treatment pilot data.
version 16
mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// Small formatting helpers
// ---------------------------------------------------------------------------

string scalar pape_fsgn(real scalar x)
    return((x >= 0 ? "+" : "") + strtrim(sprintf("%12.4f", x)))

string scalar pape_g3(real scalar x) return(strtrim(sprintf("%12.3g", x)))

string scalar pape_bigmark(real scalar x) return(strtrim(sprintf("%18.0fc", x)))

string scalar pape_numlist(real rowvector x)
{
    string scalar s
    real scalar j
    s = ""
    for (j = 1; j <= cols(x); j++) s = s + (j > 1 ? " " : "") + strofreal(x[j], "%21.0g")
    return(s)
}

// ---------------------------------------------------------------------------
// Printing a design
// ---------------------------------------------------------------------------

// 1 when the pilot's rho was solved given the design's covariates
real scalar pape_did_pilot_converted(struct pape_dgp scalar x)
{
    if (cols(x.pilot) < 10) return(0)
    if (x.pilot[10] == . | x.pilot[8] == .) return(0)
    if (strpos(" " + x.pilot_over + " ", " rho ")) return(0)
    return(abs(x.pilot[8] - x.pilot[10]) > 1e-12)
}

// what rho is, as named in the design print and the methods paragraph: with
// covariates that move the outcome, the latent correlation given them
real scalar pape_did_cov_dep(struct pape_dgp scalar x)
{
    if (x.k == 0) return(0)
    return(x.s2 > 0)
}

string scalar pape_did_pilot_line(struct pape_dgp scalar x)
{
    string rowvector est
    string scalar s
    real scalar j
    est = J(1, 0, "")
    if (!strpos(" " + x.pilot_over + " ", " baseline ")) est = est, "baseline"
    if (!strpos(" " + x.pilot_over + " ", " gap ")) est = est, "gap"
    if (!strpos(" " + x.pilot_over + " ", " ptreated ")) est = est, "treated share"
    if (x.pilot[8] != . & !strpos(" " + x.pilot_over + " ", " rho ")) est = est, "rho"
    s = ""
    for (j = 1; j <= cols(est); j++) s = s + (j > 1 ? ", " : "") + est[j]
    return(sprintf("  pilot: %s estimated from %f pre-treatment rows%s%s%s", s, x.pilot[1],
        (x.pilot[2] == . ? "" : sprintf(" (%f units)", x.pilot[2])),
        (x.pilot[9] == . ? "" : sprintf("; control group's change between pilot periods %s",
            (x.pilot[9] >= 0 ? "+" : "") + strtrim(sprintf("%12.3f", x.pilot[9])))),
        (pape_did_pilot_converted(x) ? sprintf("; rho given the covariates %5.3f (raw within-unit correlation %5.3f)",
            x.pilot[8], x.pilot[10]) : "")))
}

void pape_print_did_dgp(struct pape_dgp scalar x)
{
    string scalar labs, sig, nm, per, dat, an, tr
    real rowvector lin, idx
    real scalar j, m, M, panel, lnr, two
    panel = (x.design == "panel")
    lnr = (x.analysis == "linear")
    two = (x.n_periods == 2)
    M = (x.has_sub ? 2 : 1)
    nm = (x.has_sub ? x.subgroup.name : "")
    per = (two ? "two periods" : sprintf("%f periods: %f before, %f after", x.n_periods, x.n_pre, x.n_post))
    printf("{txt}powerape DGP -- %s, difference-in-differences%s (%s)\n", x.model,
        (x.has_sub ? " with a subgroup" : ""), per)
    if (panel)
        dat = sprintf("panel, %s per unit; rho = %4.2f (latent within-unit correlation%s); unit-clustered SEs",
            (two ? "both periods" : sprintf("all %f periods", x.n_periods)), x.rho,
            (pape_did_cov_dep(x) ? " given the covariates" : ""))
    else
        dat = sprintf("repeated cross-sections, %1.0f%% of observations in the post %s",
            100 * x.p_post, (two ? "period" : "periods"))
    printf("  data: %s\n", dat)
    if (x.has_sub) {
        if (lnr) an = sprintf("linear probability model (OLS), interacted with %s; estimate = coefficient of %s x treatment;%s",
            nm, nm, substr(pape_did_se_label(x), 2, .))
        else an = sprintf("pooled %s, fully interacted with %s; estimate = difference of the two subgroups' average partial effects of treatment",
            x.model, nm)
    }
    else {
        if (lnr) an = sprintf("linear probability model (OLS); ATT = coefficient of the treatment indicator;%s",
            substr(pape_did_se_label(x), 2, .))
        else an = sprintf("pooled %s; ATT = average partial effect of the treatment indicator over treated post-period rows", x.model)
    }
    printf("  analysis: %s\n", an)
    if (x.has_sub)
        printf("  treated group: %1.0f%% of %s; subgroup %s = 1: %1.0f%%\n", 100 * x.p_treated,
            (panel ? "units" : "observations"), nm, 100 * x.subgroup.p)
    else
        printf("  treated group: %1.0f%% of %s\n", 100 * x.p_treated, (panel ? "units" : "observations"))
    if (x.k > 0) {
        labs = ""
        for (j = 1; j <= x.k; j++) labs = labs + (j > 1 ? ", " : "") + pape_var_label(x.covariates[j])
        printf("  covariates: %s\n", labs)
        if (x.link == "logit" & x.s2 > 0)
            sig = sprintf("nuisance signal = %5.3f (latent index at unit error variance; McKelvey-Zavoina R2 on the logistic scale = %5.3f)",
                x.signal_implied, pape_mz_r2_logit(x.s2))
        else sig = sprintf("nuisance signal (latent pseudo-R2) = %5.3f", x.signal_implied)
        printf("  %s\n", sig)
    }
    else printf("  covariates: none\n")
    if (x.has_sub)
        tr = sprintf("parallel within each subgroup, %s",
            (x.trends == "index" ? sprintf("on the %s index", x.model) : "in probabilities (linear)"))
    else if (x.trends == "index") {
        if (lnr) tr = sprintf("parallel on the %s index; the linear probability analysis assumes parallel trends in probabilities", x.model)
        else tr = sprintf("parallel on the %s index (the scale the %s analysis assumes)", x.model, x.model)
    }
    else {
        if (lnr) tr = "parallel in probabilities (linear; the scale the linear probability analysis assumes)"
        else tr = sprintf("parallel in probabilities (linear); the %s analysis assumes the index scale", x.model)
    }
    printf("  untreated trends: %s\n", tr)
    for (m = 1; m <= M; m++)
        printf("  rates before -> after%s: control %5.3f -> %5.3f; treated %5.3f -> %s\n",
            (x.has_sub ? sprintf(", %s = %f", nm, m - 1) : ""), x.c00[m], x.c01[m], x.c10[m],
            (x.has_effect ? sprintf("%5.3f (%5.3f without treatment)", x.c11[m], x.cf[m]) :
                sprintf("%5.3f without treatment", x.cf[m])))
    if (cols(x.pilot) > 0) printf("%s\n", pape_did_pilot_line(x))
    if (!x.has_effect) {
        if (x.has_sub) printf("  effects not set yet - run powerape setatt %s, target(# #): the effect for %s = 0, then for %s = 1.\n", x.name, nm, nm)
        else printf("  target effect not set yet - run powerape setatt.\n")
        return
    }
    pape_did_readings(x, lin, idx)
    if (x.has_sub) {
        printf("  effects on the treated: %s (%s = 0), %s (%s = 1); true difference (DATT) pinned at %s\n",
            pape_fsgn(x.att_sub[1]), nm, pape_fsgn(x.att_sub[2]), nm, pape_fsgn(x.target_est))
        if (x.trends == "index")
            printf("  the same eight rates read under linear-parallel trends: difference %s (their triple cross difference)\n", pape_fsgn(lin[2] - lin[1]))
        else
            printf("  the same eight rates read on the %s index: difference %s (covariate-free reading)\n", x.model, pape_fsgn(idx[2] - idx[1]))
        pape_did_scale_note(x, "the true difference")
    }
    else {
        printf("  true ATT pinned at %s (treatment coefficient = %6.4f)\n", pape_fsgn(x.target_est), x.did_de[1])
        if (x.trends == "index")
            printf("  the same four rates read under linear-parallel trends: %s (their cross difference)\n", pape_fsgn(lin[1]))
        else
            printf("  the same four rates read on the %s index: %s (covariate-free reading)\n", x.model, pape_fsgn(idx[1]))
        pape_did_scale_note(x, "the ATT")
    }
}

// ---------------------------------------------------------------------------
// r() results of a design
// ---------------------------------------------------------------------------

void pape_did_rreturn(struct pape_dgp scalar d)
{
    real rowvector lin, idx, tv
    real scalar M
    M = (d.has_sub ? 2 : 1)
    pape_rstr("design", d.design)
    pape_rstr("trends", d.trends)
    pape_rstr("analysis", d.analysis)
    pape_rstr("subgroup", (d.has_sub ? d.subgroup.name : ""))
    pape_rnum("p_treated", d.p_treated)
    pape_rnum("p_post", d.p_post)
    pape_rnum("rho", d.rho)
    pape_rnum("n_pre", d.n_pre)
    pape_rnum("n_post", d.n_post)
    st_matrix("pape_r_rates", (d.c00', d.c01', d.c10', d.cf', d.c11'))
    st_matrixcolstripe("pape_r_rates", (J(5, 1, ""), ("c00" \ "c01" \ "c10" \ "cf" \ "c11")))
    st_matrix("pape_r_didcoef", (d.did_a', d.did_b', d.did_l', d.did_th', d.did_de'))
    st_matrixcolstripe("pape_r_didcoef", (J(5, 1, ""), ("a" \ "b" \ "l" \ "th" \ "de")))
    if (d.has_effect) {
        pape_did_readings(d, lin, idx)
        pape_rnum("reading_linear", (d.has_sub ? lin[2] - lin[1] : lin[1]))
        pape_rnum("reading_index", (d.has_sub ? idx[2] - idx[1] : idx[1]))
        if (d.has_sub) {
            pape_rnum("att_ref", d.att_sub[1])
            pape_rnum("att_sub", d.att_sub[2])
        }
    }
    if (cols(d.pilot) > 0) {
        tv = d.pilot
        st_matrix("pape_r_pilot", tv)
        st_matrixcolstripe("pape_r_pilot", (J(10, 1, ""), ("n_rows" \ "n_units" \ "n_periods" \ "n_pairs" \ "baseline" \ "gap" \ "p_treated" \ "rho" \ "control_change" \ "rho_raw")))
    }
}

// ---------------------------------------------------------------------------
// powerape dgp NAME, did ...
// ---------------------------------------------------------------------------

// the canonical option string of a design, from which the sweep rebuilds it
string scalar pape_did_spec(string scalar model, real rowvector baseline, real rowvector gap,
    real rowvector trend, real scalar p_treated, string scalar design, real scalar rho,
    real scalar p_post, string scalar trends, string scalar analysis, real scalar n_pre,
    real scalar n_post, string scalar sub_name, real scalar sub_p, string scalar covraw,
    string scalar corrraw, real scalar signal, string scalar gamma_str, real scalar n_int,
    real scalar seed_int)
{
    string scalar s
    s = "did model(" + model + ") baseline(" + pape_numlist(baseline) + ")"
    if (cols(gap) > 0) s = s + " gap(" + pape_numlist(gap) + ")"
    if (cols(trend) > 0) s = s + " trend(" + pape_numlist(trend) + ")"
    s = s + " ptreated(" + strofreal(p_treated, "%21.0g") + ") design(" + design + ")"
    s = s + " rho(" + strofreal(rho, "%21.0g") + ")"
    if (p_post != . & design != "panel") s = s + " ppost(" + strofreal(p_post, "%21.0g") + ")"
    s = s + " trends(" + trends + ") analysis(" + analysis + ")"
    s = s + " npre(" + strofreal(n_pre, "%21.0g") + ") npost(" + strofreal(n_post, "%21.0g") + ")"
    if (sub_name != "") s = s + " subgroup(" + sub_name + ", binary p(" + strofreal(sub_p, "%21.0g") + "))"
    if (strtrim(covraw) != "") s = s + " covariates(" + covraw + ")"
    if (strtrim(corrraw) != "") s = s + " correlation(" + corrraw + ")"
    if (strtrim(gamma_str) != "") s = s + " gamma(" + gamma_str + ")"
    else s = s + " signal(" + strofreal(signal, "%21.0g") + ")"
    s = s + " nint(" + strofreal(n_int, "%21.0g") + ") seedint(" + strofreal(seed_int, "%21.0g") + ")"
    return(s)
}

void pape_cmd_dgp_did(string scalar name, string scalar model, string scalar baseline_str,
    string scalar gap_str, string scalar trend_str, real scalar p_treated,
    string scalar design, real scalar rho, real scalar p_post, string scalar trends,
    string scalar analysis, real scalar n_pre, real scalar n_post, string scalar sub_name,
    real scalar sub_p, string scalar cnames, string scalar ctypes, string scalar cps,
    string scalar cmeans, string scalar csds, string scalar ciccs, real scalar corr_kind,
    real scalar corr_scalar, string scalar corr_matname, string scalar covraw,
    string scalar corrraw, real scalar signal, string scalar gamma_str, real scalar n_int,
    real scalar seed_int, real scalar display)
{
    struct pape_dgp scalar d
    struct pape_var scalar sv
    real matrix Cm
    real colvector gamma
    real rowvector baseline, gap, trend
    Cm = (corr_kind == 2 ? st_matrix(corr_matname) : J(0, 0, .))
    gamma = pape_num_col(gamma_str)
    baseline = pape_num_col(baseline_str)'
    gap = pape_num_col(gap_str)'
    trend = pape_num_col(trend_str)'
    if (cols(baseline) == 0) pape_err("baseline() is required.")
    if (sub_name != "") sv = pape_var_new(sub_name, "binary", sub_p, ., ., .)
    d = pape_dgp_did(name, model, baseline, gap, trend, p_treated, design, rho, p_post,
        trends, analysis, n_pre, n_post, sub_name != "", sv,
        pape_vars_parse(cnames, ctypes, cps, cmeans, csds, ciccs), corr_kind, corr_scalar, Cm,
        signal, rows(gamma) > 0, gamma, n_int, seed_int)
    d.spec_opts = pape_did_spec(model, baseline, gap, trend, p_treated, design, rho, p_post,
        trends, analysis, n_pre, n_post, sub_name, sub_p, covraw, corrraw, signal, gamma_str,
        n_int, seed_int)
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

void pape_cmd_setatt(string scalar name, string scalar targets, real scalar display)
{
    struct pape_dgp scalar d
    string scalar msg
    d = pape_get_dgp(name)
    if (pape_set_att(d, pape_num_col(targets)', msg)) pape_err(msg)
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

void pape_cmd_trueatt(string scalar name)
{
    struct pape_dgp scalar d
    real rowvector v
    d = pape_get_dgp(name)
    v = pape_true_att(d)
    if (d.has_sub) {
        printf("{txt}true effects on the treated of DGP %s on the integration set: {res}%12.0g{txt} (%s = 0), {res}%12.0g{txt} (%s = 1); difference {res}%12.0g\n",
            name, v[1], d.subgroup.name, v[2], d.subgroup.name, v[3])
        pape_rnum("att_ref", v[1])
        pape_rnum("att_sub", v[2])
        pape_rnum("value", v[3])
    }
    else {
        printf("{txt}true ATT of DGP %s on the integration set: {res}%12.0g\n", name, v[1])
        pape_rnum("value", v[1])
    }
}

// ---------------------------------------------------------------------------
// A design from pre-treatment pilot data
// ---------------------------------------------------------------------------

// P(Y_s = 1, Y_t = 1) for two periods of one unit with marginal rates p and
// q and latent correlation rho: the bivariate normal probability
real scalar pape_did_joint_rate(real scalar p, real scalar q, real scalar rho)
{
    if (rho <= 0) return(p * q)
    return(binormal(invnormal(p), invnormal(q), rho))
}

real scalar pape_f_pilot_rho(real scalar rho, struct pape_env scalar e)
{
    real scalar k, s
    s = 0
    for (k = 1; k <= rows(e.idxz); k++)
        s = s + e.m_c[k] * (pape_did_joint_rate(e.idxz[k], e.idxz[k], rho) - e.x_c[k])
    return(s)
}

// latent within-unit correlation by the method of moments; p, m, w: per
// group, the marginal rate, the observed joint rate and the number of pairs
real scalar pape_did_pilot_rho(real colvector p, real colvector m, real colvector w)
{
    struct pape_env scalar e
    e.idxz = p; e.x_c = m; e.m_c = w
    if (pape_f_pilot_rho(0, e) >= 0) return(0)
    if (pape_f_pilot_rho(0.995, e) <= 0) return(0.995)
    return(pape_uniroot(&pape_f_pilot_rho(), e, 0, 0.995, 0, 1e-9))
}

// The same moment condition for a design with covariates. Covariates are
// constant over time, so they are a second source of dependence between a
// unit's outcomes, and the design's rho is the latent correlation that is
// left GIVEN the covariates. The raw pairs of the pilot show the whole
// dependence, so the rho to use is the one at which the design's joint rate,
// averaged over the covariate index q of each group, equals the pilot's.

// an equally weighted grid standing in for the draws of q: quantile
// midpoints (R's quantile(), type 7)
real colvector pape_did_q_grid(real colvector q, real scalar K)
{
    real colvector s, h, lo, hi
    real scalar n
    n = rows(q)
    if (n <= K) return(q)
    s = sort(q, 1)
    h = (n - 1) :* (((1::K) :- 0.5) :/ K) :+ 1
    lo = floor(h); hi = ceil(h)
    return(s[lo] :+ (h :- lo) :* (s[hi] :- s[lo]))
}

// P(Y_s = 1, Y_t = 1) for two periods of one unit whose rate is G(cc + q),
// averaged over the values q, at latent correlation rho given q
real scalar pape_did_joint_rate_cov(string scalar link, real scalar cc, real colvector q,
    real scalar rho)
{
    real colvector p, thr
    p = pape_G(link, cc :+ q)
    if (rho <= 0) return(mean(p:^2))
    p = (p :< 1e-15) :* 1e-15 :+ (p :>= 1e-15) :* p
    p = (p :> 1 - 1e-15) :* (1 - 1e-15) :+ (p :<= 1 - 1e-15) :* p
    thr = invnormal(p)
    return(mean(binormal(thr, thr, J(rows(thr), 1, rho))))
}

// e.idxz, e.x_c: the covariate index of the control and of the treated
// group; e.b0, e.bf: their intercepts; e.p, e.r: their numbers of pairs;
// e.s, e.t: their observed joint rates
real scalar pape_f_pilot_rho_cov(real scalar rho, struct pape_env scalar e)
{
    real scalar v
    v = 0
    if (e.p > 0) v = v + e.p * (pape_did_joint_rate_cov(e.link, e.b0, e.idxz, rho) - e.s)
    if (e.r > 0) v = v + e.r * (pape_did_joint_rate_cov(e.link, e.bf, e.x_c, rho) - e.t)
    return(v)
}

// rho given the covariates; the group intercepts are solved at the pilot's
// own rates. exceeds = 1 when the covariates alone already imply at least
// the dependence the pilot shows (then rho = 0)
real scalar pape_did_pilot_rho_cov(string scalar link, real colvector p, real colvector m,
    real colvector w, real colvector q0, real colvector q1, real scalar exceeds)
{
    struct pape_env scalar e
    e.link = link
    e.idxz = pape_did_q_grid(q0, 2000)
    e.x_c = pape_did_q_grid(q1, 2000)
    e.b0 = pape_did_solve(link, e.idxz, p[1], 0)
    e.bf = pape_did_solve(link, e.x_c, p[2], 0)
    e.p = w[1]; e.r = w[2]; e.s = m[1]; e.t = m[2]
    exceeds = 0
    if (pape_f_pilot_rho_cov(0, e) >= 0) {
        exceeds = 1
        return(0)
    }
    if (pape_f_pilot_rho_cov(0.995, e) <= 0) return(0.995)
    return(pape_uniroot(&pape_f_pilot_rho_cov(), e, 0, 0.995, 0, 1e-9))
}

// estimates baseline, gap, the treated share and (panel) rho from the data
// in memory, then builds the design; *_given flags mark user overrides
void pape_cmd_dgp_didpilot(string scalar name, string scalar yvar, string scalar gvar,
    string scalar idvar, string scalar pervar, string scalar touse, string scalar model,
    string scalar baseline_str, string scalar gap_str, string scalar trend_str,
    real scalar p_treated_in, string scalar design_in, real scalar rho_in, real scalar p_post,
    string scalar trends, string scalar analysis, real scalar n_pre, real scalar n_post,
    string scalar cnames, string scalar ctypes, string scalar cps, string scalar cmeans,
    string scalar csds, string scalar ciccs, real scalar corr_kind, real scalar corr_scalar,
    string scalar corr_matname, string scalar covraw, string scalar corrraw,
    real scalar signal, string scalar gamma_str, real scalar n_int, real scalar seed_int,
    real scalar display)
{
    struct pape_dgp scalar d
    struct pape_var scalar sv
    real colvector y, g, u, ord, gi, Ti, si, per, pu, r0, gamma, w, m
    real matrix info, Cm
    real rowvector baseline, gap, trend
    real scalar panel, p0, p1, n_units, n_pairs, rho_hat, pt_hat, ctrl_change, n_per, j, k
    real scalar p_treated, rho, rho_raw, exceeds
    real colvector q00, q01, q10, q11
    string scalar design, over
    y = st_data(., yvar, touse)
    g = st_data(., gvar, touse)
    if (rows(y) == 0) pape_err("The pilot data have no rows.")
    if (!all(y :== 0 :| y :== 1) | !all(g :== 0 :| g :== 1))
        pape_err("outcome() and group() must be binary (0/1).")
    if (!any(g :== 0) | !any(g :== 1)) pape_err("The pilot needs observations of both groups.")
    panel = (idvar != "")
    p0 = mean(select(y, g :== 0))
    p1 = mean(select(y, g :== 1))
    if (p0 <= 0 | p0 >= 1 | p1 <= 0 | p1 >= 1)
        pape_err("A group has no events or no non-events in the pilot: its rate cannot be planned from these data.")
    rho_hat = .; n_pairs = 0; n_units = .
    if (panel) {
        u = st_data(., idvar, touse)
        ord = order(u, 1)
        info = panelsetup(u[ord], 1)
        n_units = rows(info)
        gi = J(n_units, 1, .); Ti = J(n_units, 1, .); si = J(n_units, 1, .)
        for (j = 1; j <= n_units; j++) {
            gi[j] = g[ord[info[j, 1]]]
            if (any(g[ord[info[j, 1]::info[j, 2]]] :!= gi[j])) pape_err("group() must be constant within a unit.")
            Ti[j] = info[j, 2] - info[j, 1] + 1
            si[j] = sum(y[ord[info[j, 1]::info[j, 2]]])
        }
        pt_hat = mean(gi)
        w = J(2, 1, 0); m = J(2, 1, 0)
        for (k = 0; k <= 1; k++) {
            w[k + 1] = sum(select(Ti :* (Ti :- 1), gi :== k))
            if (w[k + 1] > 0) m[k + 1] = sum(select(si :* (si :- 1), gi :== k)) / w[k + 1]
        }
        n_pairs = sum(w) / 2
        if (n_pairs > 0) rho_hat = pape_did_pilot_rho((p0 \ p1), m, w)
        else if (rho_in == .)
            pape_err("Every unit has one pre-treatment period, so the dependence between a unit's periods cannot be estimated: supply rho().")
    }
    else pt_hat = mean(g)
    // the control group's change between pre-treatment periods, as a guide
    ctrl_change = .; n_per = .
    if (pervar != "") {
        per = st_data(., pervar, touse)
        pu = uniqrows(per)
        n_per = rows(pu)
        pu = uniqrows(select(per, g :== 0))
        if (rows(pu) >= 2) {
            r0 = J(rows(pu), 1, .)
            for (j = 1; j <= rows(pu); j++) r0[j] = mean(select(y, (g :== 0) :& (per :== pu[j])))
            ctrl_change = mean(r0[2::rows(pu)] :- r0[1::(rows(pu) - 1)])
        }
    }
    trend = pape_num_col(trend_str)'
    if (cols(trend) == 0)
        pape_err("trend() is required: pre-treatment data do not show how the control group will change from before to after treatment." +
            (ctrl_change == . ? "" : sprintf(" (Between the pilot's pre-treatment periods the control group's rate changed by %s%5.3f on average.)", (ctrl_change >= 0 ? "+" : ""), ctrl_change)))
    over = ""
    baseline = pape_num_col(baseline_str)'
    if (cols(baseline) > 0) over = over + " baseline"
    else baseline = p0
    gap = pape_num_col(gap_str)'
    if (cols(gap) > 0) over = over + " gap"
    else gap = p1 - p0
    if (p_treated_in != .) {
        over = over + " ptreated"
        p_treated = p_treated_in
    }
    else p_treated = pt_hat
    if (design_in != "") design = design_in
    else design = (panel ? "panel" : "cross-sections")
    if (rho_in != .) {
        rho = rho_in
        if (panel) over = over + " rho"
    }
    else rho = (rho_hat == . | design == "cross-sections" ? 0 : rho_hat)
    if (rows(y) < 500)
        printf("{err}warning: the pilot has fewer than 500 rows: the estimated rates and the within-unit correlation are imprecise. Vary them with powerape robust before relying on a sample size.\n")
    Cm = (corr_kind == 2 ? st_matrix(corr_matname) : J(0, 0, .))
    gamma = pape_num_col(gamma_str)
    d = pape_dgp_did(name, model, baseline, gap, trend, p_treated, design, rho, p_post,
        trends, analysis, n_pre, n_post, 0, sv,
        pape_vars_parse(cnames, ctypes, cps, cmeans, csds, ciccs), corr_kind, corr_scalar, Cm,
        signal, rows(gamma) > 0, gamma, n_int, seed_int)
    // covariates add within-unit dependence of their own: the design's rho
    // is the part left given them. The calibration does not depend on rho,
    // so the design is updated in place.
    rho_raw = rho_hat
    if (rho_hat != . & rho_in == . & d.design == "panel" & d.k > 0 & d.s2 > 0) {
        pape_did_qsplit(d, q00, q01, q10, q11)
        rho_hat = pape_did_pilot_rho_cov(d.link, (p0 \ p1), m, w, q00, q10, exceeds)
        if (exceeds)
            printf("{err}warning: the covariates as specified imply more within-unit dependence than the pilot shows (its raw latent correlation is %5.3f): rho given the covariates is set to 0, and the design still has more dependence than the pilot. The stated signal looks too high for these data.\n", rho_raw)
        rho = rho_hat
        d.rho = rho
    }
    d.spec_opts = pape_did_spec(model, baseline, gap, trend, p_treated, design, rho, p_post,
        trends, analysis, n_pre, n_post, "", ., covraw, corrraw, signal, gamma_str, n_int, seed_int)
    d.pilot = (rows(y), n_units, n_per, n_pairs, p0, p1 - p0, pt_hat, rho_hat, ctrl_change, rho_raw)
    d.pilot_over = strtrim(over)
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

// ---------------------------------------------------------------------------
// Methods paragraph
// ---------------------------------------------------------------------------

string scalar pape_did_nper(real scalar k)
    return(k == 1 ? "one period" : sprintf("%f periods", k))

string scalar pape_did_rate_txt(struct pape_dgp scalar d, real scalar m)
    return(sprintf("%5.3f and %5.3f in the control group and %5.3f and %5.3f in the treated group (%5.3f without treatment)",
        d.c00[m], d.c01[m], d.c10[m], d.c11[m], d.cf[m]))

string scalar pape_did_pilot_clause(struct pape_dgp scalar d)
{
    string rowvector lab
    string scalar s
    real scalar j
    if (cols(d.pilot) == 0) return("")
    lab = J(1, 0, "")
    if (!strpos(" " + d.pilot_over + " ", " baseline ")) lab = lab, "the control group's rate before treatment"
    if (!strpos(" " + d.pilot_over + " ", " gap ")) lab = lab, "the gap between the groups"
    if (!strpos(" " + d.pilot_over + " ", " ptreated ")) lab = lab, "the treated share"
    if (d.pilot[8] != . & !strpos(" " + d.pilot_over + " ", " rho "))
        lab = lab, (pape_did_pilot_converted(d) ? "the within-unit correlation given the covariates" : "the within-unit correlation")
    if (cols(lab) == 0) return("")
    if (cols(lab) == 1) s = lab[1]
    else {
        s = ""
        for (j = 1; j < cols(lab); j++) s = s + (j > 1 ? ", " : "") + lab[j]
        s = s + " and " + lab[cols(lab)]
    }
    return(sprintf("%s estimated from %f pre-treatment pilot observations", s, d.pilot[1]))
}

// the assumed world of a difference-in-differences design, in words
string scalar pape_describe_did(struct pape_dgp scalar d)
{
    string scalar s_design, s_data, s_rates, s_trend, s_cov, sig, dep, nm, vs, pc
    real matrix R
    real colvector off
    real scalar panel, i, j
    panel = (d.design == "panel")
    s_design = sprintf("a difference-in-differences design with a treated group (%1.0f%% of the %s) and a control group, observed in %s before and %s after one common treatment date",
        100 * d.p_treated, (panel ? "units" : "observations"), pape_did_nper(d.n_pre), pape_did_nper(d.n_post))
    if (panel) s_data = sprintf("as a panel of the same units, with a latent within-unit correlation of %4.2f%s shared by all periods (Gaussian copula)", d.rho,
        (pape_did_cov_dep(d) ? " given the covariates," : ""))
    else s_data = sprintf("as repeated cross-sections, with %1.0f%% of the observations after the treatment date", 100 * d.p_post)
    if (d.has_sub) {
        nm = d.subgroup.name
        s_rates = sprintf("a binary subgroup variable %s (prevalence %4.2f); outcome rates before and after treatment of %s for %s = 0, and of %s for %s = 1",
            nm, d.subgroup.p, pape_did_rate_txt(d, 1), nm, pape_did_rate_txt(d, 2), nm)
    }
    else s_rates = sprintf("outcome rates before and after treatment of %s", pape_did_rate_txt(d, 1))
    s_trend = sprintf("untreated trends parallel %s%s",
        (d.trends == "index" ? sprintf("on the %s index", d.model) : "in probabilities"),
        (d.has_sub ? " within each subgroup" : ""))
    if (d.k > 0) {
        if (d.link == "logit" & d.s2 > 0)
            sig = sprintf("contributing %4.2f of the latent-index variance at unit error variance (a McKelvey-Zavoina R-squared of %5.3f on the logistic scale)", d.signal_implied, pape_mz_r2_logit(d.s2))
        else sig = sprintf("contributing a latent pseudo-R-squared of %4.2f", d.signal_implied)
        off = J(0, 1, .)
        if (d.corr_given == 2) {
            R = d.corr_matrix
            for (i = 1; i <= rows(R); i++) for (j = i + 1; j <= cols(R); j++) off = off \ R[i, j]
        }
        else if (d.corr_given == 1) off = d.corr_scalar
        if (rows(off) == 0) dep = "independent of group membership"
        else if (all(abs(off) :< 1e-12)) dep = "independent of group membership"
        else if (d.corr_given == 2)
            dep = sprintf("with a latent (Gaussian-copula) correlation matrix over group membership and the covariates (off-diagonal entries from %4.2f to %4.2f)", min(off), max(off))
        else dep = sprintf("with an exchangeable latent (Gaussian-copula) correlation of %4.2f among group membership and the covariates", d.corr_scalar)
        vs = ""
        for (j = 1; j <= d.k; j++) {
            if (d.covariates[j].type == "binary")
                vs = vs + (j > 1 ? ", " : "") + sprintf("%s (binary, prevalence %4.2f)", d.covariates[j].name, d.covariates[j].p)
            else
                vs = vs + (j > 1 ? ", " : "") + sprintf("%s (continuous, mean %s, SD %s)", d.covariates[j].name, pape_g3(d.covariates[j].mean), pape_g3(d.covariates[j].sd))
        }
        s_cov = sprintf("%scovariates %s, %s, %s", (panel ? "time-constant " : ""), vs, dep, sig)
    }
    else s_cov = "no covariates"
    pc = pape_did_pilot_clause(d)
    return(s_design + " " + s_data + "; " + s_rates + "; " + s_trend + "; " + s_cov + (pc != "" ? "; " + pc : ""))
}

// how the simulated studies were analysed, in words (one clause)
string scalar pape_did_analysis_text(struct pape_dgp scalar d, string scalar se)
{
    string scalar regs, se_txt
    regs = "group, period indicators and the treatment indicator" + (d.k > 0 ? ", and the covariates" : "") +
        (d.has_sub ? sprintf(", fully interacted with %s", d.subgroup.name) : "")
    if (d.analysis == "linear")
        return(sprintf("ordinary least squares of the outcome on %s (a linear probability model), with the %s as the coefficient of %s and %s standard errors",
            regs, (d.has_sub ? "difference between the subgroups' effects" : "effect on the treated"),
            (d.has_sub ? sprintf("the treatment indicator's product with %s", d.subgroup.name) : "the treatment indicator"),
            (d.design == "panel" ? "unit-clustered" : "heteroskedasticity-robust (HC1)")))
    if (d.design == "panel") se_txt = "unit-clustered standard errors"
    else if (se == "robust") se_txt = "heteroskedasticity-robust (sandwich) standard errors"
    else se_txt = "model-based standard errors"
    return(sprintf("a pooled %s regression of the outcome on %s, with %s computed as %s the average partial effect of the treatment indicator over the treated post-treatment observations (Wooldridge, 2023), and delta-method Wald confidence intervals from %s",
        d.model, regs, (d.has_sub ? "the estimate" : "the effect on the treated"),
        (d.has_sub ? "the difference between the two subgroups in" : ""), se_txt))
}

// the result sentence shared by power, n and mde results
string scalar pape_stmt_result(struct pape_res scalar x, string scalar n_text,
    string scalar est_short, string scalar size_noun)
{
    string scalar conf_frag
    real scalar rule
    rule = x.has_rule
    if (x.kind == "mde") {
        if (x.confirmed) conf_frag = sprintf(", confirmed by a high-precision verification run (simulated power %5.3f, Monte Carlo SE %5.3f, %f replications)", x.power, x.mcse, x.nsim_confirm)
        else conf_frag = sprintf(" (simulated power %5.3f, Monte Carlo SE %5.3f)", x.power, x.mcse)
        if (x.is_margin)
            return(sprintf("At a sample size of %s, the tightest equivalence bounds establishable with %2.0f%% power were +/-%5.3f (%3.1f percentage points)%s.", n_text, 100 * x.goal, x.mde, 100 * x.mde, conf_frag))
        return(sprintf("At a sample size of %s, the minimum detectable %s at %2.0f%% power for %s was %5.3f (%s of %3.1f percentage points)%s.",
            n_text, est_short, 100 * x.goal, pape_claim_text(x.claim, x.sesoi, x.conf, rule), x.mde,
            (x.mde < 0 ? "a decrease" : "an increase"), 100 * abs(x.mde), conf_frag))
    }
    if (x.kind == "n") {
        if (x.confirmed)
            return(sprintf("The required %s for %2.0f%% power for %s was %s, confirmed by a high-precision verification run (simulated power %5.3f, Monte Carlo SE %5.3f, %f replications).",
                size_noun, 100 * x.goal, pape_claim_text(x.claim, x.sesoi, x.conf, rule), n_text, x.power, x.mcse, x.nsim_confirm))
        return(sprintf("The required %s for %2.0f%% power for %s was %s (simulated power %5.3f, Monte Carlo SE %5.3f, %f replications per search step).",
            size_noun, 100 * x.goal, pape_claim_text(x.claim, x.sesoi, x.conf, rule), n_text, x.power, x.mcse, x.nsim))
    }
    return(sprintf("At a sample size of %s, simulated power for %s was %5.3f (Monte Carlo SE %5.3f; %f replications).",
        n_text, pape_claim_text(x.claim, x.sesoi, x.conf, rule), x.power, x.mcse, x.nsim))
}

string scalar pape_did_statement(struct pape_res scalar x)
{
    struct pape_dgp scalar d
    string scalar est_short, est_long, s1, s2, s3, s4, s5, s6, n_text, an, cell
    real scalar keep, sepd, other, failed
    d = x.dgp
    pape_did_power_only(d, "powerape statement")
    est_short = (d.has_sub ? "difference in ATTs" : "ATT")
    if (d.has_sub) est_long = sprintf("difference between two subgroups' average treatment effects on the treated (by %s)", d.subgroup.name)
    else est_long = "average treatment effect on the treated (ATT)"
    s1 = sprintf("We conducted a simulation-based power analysis for the %s in a difference-in-differences design, in probability points, using the powerape package for Stata (version %s), following the confidence-interval approach of Riesthuis (2024).",
        est_long, pape_version())
    s2 = sprintf("The assumed data-generating process was %s.", pape_describe_did(d))
    if (x.kind == "mde") s3 = ""
    else if (d.has_sub)
        s3 = sprintf("The assumed effects on the treated were %5.3f (%s = 0) and %5.3f (%s = 1), a true difference of %5.3f (%3.1f percentage points).",
            d.att_sub[1], d.subgroup.name, d.att_sub[2], d.subgroup.name, x.target, 100 * x.target)
    else s3 = sprintf("The assumed true ATT was %5.3f (%3.1f percentage points).", x.target, 100 * x.target)
    if (d.design == "panel")
        n_text = sprintf("n = %f units observed in %f periods (%s unit-period observations)", x.n, d.n_periods, pape_bigmark(x.n * d.n_periods))
    else n_text = sprintf("n = %f observations", x.n)
    s4 = pape_stmt_result(x, n_text, est_short, (d.design == "panel" ? "number of units" : "total sample size"))
    if (x.kind == "power" & x.has_sesoi)
        s5 = sprintf("Across replications, the probability of concluding a meaningful effect was %5.3f, of detection without meaningfulness %5.3f, of an inconclusive result %5.3f, and of equivalence %5.3f.",
            x.outcomes[1], x.outcomes[2], x.outcomes[3], x.outcomes[4])
    else s5 = ""
    keep = (x.separation == "keep")
    sepd = (x.separated == . ? 0 : x.separated)
    cell = pape_cell_label(d.estimand)
    an = pape_did_analysis_text(d, x.se)
    failed = (x.kind == "power" ? x.outcomes[cols(x.outcomes)] : .)
    if (d.analysis == "linear") {
        // ordinary least squares has no separation
        if (x.kind == "power")
            s6 = sprintf("Each simulated study was analysed by %s; replications whose fit failed (%3.1f%%) counted against the claim.", an, 100 * failed)
        else s6 = sprintf("Each simulated study was analysed by %s; replications whose fit failed counted against the claim.", an)
    }
    else if (x.kind == "power") {
        other = failed - (keep ? 0 : sepd)
        if (keep)
            s6 = sprintf("Each simulated study was analysed by %s; replications whose model fit failed (%3.1f%%) counted against the claim, and replications with separation (a %s with no events or no non-events; %3.1f%%) were retained with their Wald intervals.",
                an, 100 * other, cell, 100 * sepd)
        else
            s6 = sprintf("Each simulated study was analysed by %s; replications without a maximum-likelihood estimate, from a failed fit or from separation (a %s with no events or no non-events), counted against the claim (%3.1f%% of replications%s).",
                an, cell, 100 * failed, (sepd > 0 ? sprintf(", %3.1f%% with separation", 100 * sepd) : ""))
    }
    else
        s6 = sprintf("Each simulated study was analysed by %s; replications without a maximum-likelihood estimate (failed fits%s) counted against the claim%s.",
            an, (keep ? "" : sprintf(" and separation, a %s with no events or no non-events", cell)),
            (keep ? ", while replications with separation were retained with their Wald intervals" : ""))
    return(s1 + " " + s2 + " " + (s3 != "" ? s3 + " " : "") + s4 + " " + (s5 != "" ? s5 + " " : "") + s6)
}

// ---------------------------------------------------------------------------
// Robustness sweep (robust.R: ape_robust, difference-in-differences branch)
// ---------------------------------------------------------------------------

// replace name(...) in an option string by name(text)
string scalar pape_opt_replace_str(string scalar opts, string scalar nm, string scalar val)
{
    real scalar i, j, depth
    string scalar s
    s = opts
    i = strpos(s, nm + "(")
    while (i > 0) {
        if (i == 1) break
        if (strpos(" ,", substr(s, i - 1, 1))) break
        j = strpos(substr(s, i + 1, .), nm + "(")
        i = (j > 0 ? i + j : 0)
    }
    if (i == 0) return(s + " " + nm + "(" + val + ")")
    depth = 0
    for (j = i + strlen(nm); j <= strlen(s); j++) {
        if (substr(s, j, 1) == "(") depth++
        else if (substr(s, j, 1) == ")") {
            depth--
            if (depth == 0) break
        }
    }
    return(substr(s, 1, i - 1) + nm + "(" + val + ")" + substr(s, j + 1, .))
}

void pape_cmd_robust_did(string scalar dname, real scalar n, string scalar claim,
    real scalar sesoi, real scalar alpha, real scalar conf, string scalar rule,
    real scalar nsim, real scalar seed, string scalar varynames, string scalar varyvals,
    string scalar pin, string scalar mode, real scalar power, real scalar grid_points,
    real scalar nmax, real scalar nmax_power, string scalar se, string scalar direction,
    string scalar separation, string scalar resname, real scalar display)
{
    struct pape_dgp scalar dgp, d_i
    struct pape_res scalar pw, m_i, nr, best_nr
    real rowvector lv, tw, target0
    string rowvector vn, vv, allowed, notes, names
    string scalar opts, msg, dir, se_use, fitlab, eng, eff, lab
    real matrix grid, scen, cand_tab, L, out
    real colvector key, crit, ok_idx, cr, o, cand, nn, nl, valid, est_mean
    real scalar nv, j, i, G, rc, rule_on, worst_i, nmax_i, k, sk, dist, best_j, world, jt
    struct pape_dgp scalar few_d
    real scalar few_min, fcnt
    string scalar fcell
    dgp = pape_get_dgp(dname)
    // the smallest group over the base design and the scenarios
    few_min = pape_did_few_units(dgp, n, fcell)
    few_d = dgp
    // the base design must be analysed on the scale its trends are parallel
    // on; the other scale is a scenario (vary(trends(index linear)))
    pape_did_need_valid(dgp, "powerape robust")
    if (pin == "coefficients")
        pape_err("pin(coefficients) is not available for difference-in-differences designs: the sweep re-pins the effect on the treated in every scenario (pin(ape)).")
    if (!dgp.has_effect) pape_err("This DGP has no target effect yet; run powerape setatt first.")
    se_use = se
    if (se != "model" & pape_did_se_fixed(dgp)) {
        if (dgp.analysis != "linear")
            printf("{err}warning: se() is fixed by the route; se(%s) is ignored.\n", se)
        se_use = "model"
    }
    eng = pape_engine()
    eff = (eng == "stata" ? pape_effects() : "closed")
    fitlab = pape_fit_lab(dgp, se_use, separation, eng, eff)
    if (n == . | n < 20) pape_err("n() must be at least 20.")
    lv = pape_claim_levels(claim, alpha, conf, n, rule)
    conf = lv[4]
    rule_on = (rule != "")
    if (mode == "mde" & claim == "equivalence")
        pape_err("MDE mode searches the effect for detect/minimum claims; for the equivalence analog run powerape mde, claim(equivalence) per scenario directly.")
    if (mode == "power") pape_check_coherence(dgp, claim, sesoi, conf)
    dir = direction
    if (dir == "") dir = (dgp.target_est < 0 ? "negative" : "positive")
    vn = tokens(varynames)
    nv = cols(vn)
    if (nv == 0) pape_err("vary() must name at least one contextual input.")
    allowed = J(1, 0, "")
    if (!dgp.has_sub) allowed = ("baseline", "gap", "trend")
    if (dgp.design == "panel") allowed = allowed, "rho"
    allowed = allowed, "ptreated", "signal", "trends"
    for (j = 1; j <= nv; j++) {
        if (!anyof(allowed, vn[j]))
            pape_err(sprintf("Cannot vary %s for this DGP route; varyable inputs: %s.", vn[j], invtokens(allowed, ", ")))
    }
    if (anyof(vn, "signal") & dgp.gamma_given) pape_err("Cannot vary signal when the DGP was built with explicit gamma().")
    world = anyof(vn, "trends")
    jt = (world ? selectindex(vn :== "trends")[1] : 0)
    vv = tokens(varyvals, "|")
    vv = select(vv, vv :!= "|")
    if (cols(vv) != nv) pape_err("vary(): internal parse error.")
    // level lists in a padded matrix (one row per input; nl = length); the
    // trends scale is coded 1 = index, 2 = linear
    nl = J(nv, 1, .)
    L = J(nv, 0, .)
    G = 1
    for (j = 1; j <= nv; j++) {
        if (j == jt) {
            tw = J(1, 0, .)
            names = tokens(strlower(vv[j]))
            if (cols(names) != 2) pape_err("vary(trends()) must be trends(index linear).")
            for (k = 1; k <= 2; k++) {
                if (names[k] == "index") tw = tw, 1
                else if (names[k] == "linear") tw = tw, 2
                else pape_err("vary(trends()) must be trends(index linear).")
            }
            if (tw[1] == tw[2]) pape_err("vary(trends()) must be trends(index linear).")
            key = tw'
        }
        else {
            key = pape_num_col(vv[j])
            if (rows(key) < 2 | hasmissing(key)) pape_err(sprintf("vary(%s()) needs at least two values.", vn[j]))
            if (rows(key) == 2) key = rangen(key[1], key[2], grid_points)
        }
        nl[j] = rows(key)
        if (rows(key) > cols(L)) L = (L, J(nv, rows(key) - cols(L), .))
        L[j, 1..rows(key)] = key'
        G = G * rows(key)
    }
    grid = J(G, nv, .)
    for (i = 1; i <= G; i++) {
        k = i - 1
        for (j = 1; j <= nv; j++) {
            grid[i, j] = L[j, mod(k, nl[j]) + 1]
            k = floor(k / nl[j])
        }
    }
    target0 = (dgp.has_sub ? dgp.att_sub : dgp.target_est)
    // scenario table: inputs, implied effect or mde, power, mcse, separated, feasible
    scen = J(G, nv + 5, .)
    notes = J(1, G, "")
    crit = J(G, 1, .)
    valid = J(G, 1, 1)
    est_mean = J(G, 1, .)
    pape_set_quiet(1)
    for (i = 1; i <= G; i++) {
        opts = dgp.spec_opts
        for (j = 1; j <= nv; j++) {
            if (j == jt) opts = pape_opt_replace_str(opts, "trends", (grid[i, j] == 1 ? "index" : "linear"))
            else opts = pape_opt_replace(opts, vn[j], grid[i, j])
        }
        scen[i, 1..nv] = grid[i, .]
        rc = _stata("capture quietly powerape dgp _pape_scn, " + opts + " nodisplay replace", 1)
        if (rc) {
            notes[i] = pape_last_err()
            scen[i, nv + 5] = 0
            continue
        }
        d_i = pape_get_dgp("_pape_scn")
        valid[i] = pape_did_valid(d_i)
        fcnt = pape_did_few_units(d_i, n, fcell)
        if (fcnt < few_min) {
            few_min = fcnt
            few_d = d_i
        }
        if (mode == "mde") {
            // a subgroup design must be pinned before its effect is searched
            if (dgp.has_sub) {
                if (pape_set_att(d_i, dgp.att_sub, msg)) {
                    notes[i] = msg
                    scen[i, nv + 5] = 0
                    continue
                }
            }
            m_i = pape_mde_safe(d_i, n, power, claim, sesoi, alpha, conf, rule, nsim,
                (seed == . ? . : seed + i), ., ., 6, 0, 4 * nsim, se_use, dir, separation, rc, msg)
            if (rc) {
                notes[i] = msg
                scen[i, nv + 5] = 0
                continue
            }
            scen[i, nv + 1] = m_i.mde; scen[i, nv + 2] = m_i.power
            scen[i, nv + 3] = m_i.mcse; scen[i, nv + 4] = m_i.separated
            scen[i, nv + 5] = 1
            continue
        }
        if (pape_set_att(d_i, target0, msg)) {
            notes[i] = msg
            scen[i, nv + 5] = 0
            continue
        }
        pw = pape_power_once(d_i, n, claim, sesoi, conf, nsim, (seed == . ? . : seed + i), 0, se_use, ., "", 1, separation)
        scen[i, nv + 1] = d_i.target_est; scen[i, nv + 2] = pw.power
        scen[i, nv + 3] = pw.mcse; scen[i, nv + 4] = pw.separated
        scen[i, nv + 5] = 1
        est_mean[i] = pw.est_mean
        d_i.name = sprintf("_pape_scn%f", i)
        pape_put_dgp(d_i)
        if (!valid[i]) continue           // a rejection rate: shown, not ranked
        dist = pape_claim_distance(claim, d_i.target_est, sesoi)
        if (any(pw.draws.ok)) crit[i] = (dist > 0 ? mean(select(pw.draws.se, pw.draws.ok)) * sqrt(n) / dist : .)
        if (dist <= 0 & any(pw.draws.ok)) crit[i] = 1e300
    }
    pape_set_quiet(0)
    // few independent units behind one group, in the base design or in a
    // scenario (the treated share may be varied): one warning, for the
    // smallest such group
    if (few_min != .)
        pape_did_few_warning(few_d, n, sprintf("At n = %f, in the scenario with the smallest group", n))
    key = (mode == "mde" ? abs(scen[., nv + 1]) : scen[., nv + 2])
    ok_idx = selectindex((key :!= .) :& valid)
    if (rows(ok_idx) == 0 | cols(ok_idx) == 0) pape_err("No scenario was feasible; widen or shift the vary() ranges.")
    ok_idx = vec(ok_idx)
    if (mode == "mde") worst_i = ok_idx[selectindex(key[ok_idx] :== max(key[ok_idx]))[1]]
    else {
        cr = crit[ok_idx]
        cr = editmissing(cr, -1e300)
        o = order((scen[ok_idx, nv + 2], -cr), (1, 2))
        worst_i = ok_idx[o[1]]
    }
    nmax_i = .
    cand_tab = J(0, nv + 2, .)
    if (nmax & mode == "power") {
        cr = crit[ok_idx]
        if (any(cr :== 1e300)) {
            j = ok_idx[selectindex(cr :== 1e300)[1]]
            printf("{err}warning: n_max not reported: in %f scenario(s) the implied effect sits on or beyond the claim boundary (e.g. row %f), so no sample size reaches the target power there.\n",
                sum(cr :== 1e300), j)
        }
        else if (all(cr :== .)) printf("{err}warning: n_max not reported: no scenario produced usable estimates at this n to rank the scenarios by.\n")
        else {
            o = ok_idx[order(-cr, 1)]
            o = select(o, crit[o] :!= .)
            cand = select(o, crit[o] :>= 0.95 * crit[o[1]])
            cand = cand[1..min((3, rows(cand)))]
            nn = J(rows(cand), 1, .)
            best_j = 0
            for (k = 1; k <= rows(cand); k++) {
                sk = (seed == . ? . : seed + G + 1 + 1000 * (k - 1))
                d_i = pape_get_dgp(sprintf("_pape_scn%f", cand[k]))
                nr = pape_n(d_i, nmax_power, claim, sesoi, (rule_on ? alpha : lv[1]), (rule_on ? . : conf), rule,
                    nsim, sk, 30, 2e6, 5, 1, 4 * nsim, se_use, separation)
                nn[k] = nr.n
                if (best_j == 0) {
                    best_j = k
                    best_nr = nr
                }
                else if (nr.n > nn[best_j]) {
                    best_j = k
                    best_nr = nr
                }
            }
            nmax_i = cand[best_j]
            cand_tab = (grid[cand, .], crit[cand], nn)
        }
    }
    for (i = 1; i <= G; i++) {
        if (pape_store_has("DGP", sprintf("_pape_scn%f", i))) pape_store_drop("DGP", sprintf("_pape_scn%f", i))
    }
    if (pape_store_has("DGP", "_pape_scn")) pape_store_drop("DGP", "_pape_scn")
    lab = pape_est_label(dgp.estimand)
    if (display) {
        if (dgp.analysis == "linear")
            printf("{txt}powerape -- linear probability analysis (OLS)%s\n", pape_did_se_label(dgp))
        if (mode == "mde") {
            printf("{txt}powerape robustness sweep -- minimum detectable %s, %s claim%s\n", lab, claim, (dir == "negative" ? ", decreases" : ""))
            printf("  n = %f, %2.0f%% target power, nsim = %f per scenario, %f scenario(s) over: %s\n",
                n, 100 * power, nsim, G, invtokens(vn, ", "))
            if (any(!valid))
                printf("  note: in %f scenario(s) (valid = 0) the analysis assumes another parallel-trends scale than the world has; no effect is searched there.\n", sum(!valid))
            printf("  MDE range [%9.4f, %9.4f]; worst (largest-MDE) scenario: row %f\n",
                min(scen[ok_idx, nv + 1]), max(scen[ok_idx, nv + 1]), worst_i)
        }
        else {
            printf("{txt}powerape robustness sweep -- %s claim, %s, pin = %s\n", claim, lab, pin)
            printf("  n = %f, nsim = %f per scenario, %f scenario(s) over: %s\n", n, nsim, G, invtokens(vn, ", "))
            if (any(!valid)) {
                printf("  note: in %f scenario(s) (valid = 0) the analysis assumes another parallel-trends scale than the world has: power is a rejection rate there\n", sum(!valid))
                printf("    and est_mean shows where the estimates centre (the true %s is %s); they are left out of the range, the worst case and n_max.\n", lab, pape_fsgn(dgp.target_est))
            }
            printf("  power range [%5.3f, %5.3f]; worst scenario: row %f\n",
                min(scen[ok_idx, nv + 2]), max(scen[ok_idx, nv + 2]), worst_i)
        }
        printf("  fit: %s\n", fitlab)
        printf("{txt}%4s", "row")
        for (j = 1; j <= nv; j++) printf(" %12s", vn[j])
        printf(" %14s %8s %8s %9s", (mode == "mde" ? "mde" : "implied_eff"), "power", "mcse", "separated")
        if (world) printf("%s %6s", (mode == "power" ? sprintf(" %9s", "est_mean") : ""), "valid")
        printf("\n")
        for (i = 1; i <= G; i++) {
            printf("{res}%4.0f", i)
            for (j = 1; j <= nv; j++) {
                if (j == jt) printf(" %12s", (scen[i, j] == 1 ? "index" : "linear"))
                else printf(" %12.4g", scen[i, j])
            }
            if (scen[i, nv + 5] == 1) {
                printf(" %14.4f %8.4f %8.4f %9.4f", scen[i, nv + 1], scen[i, nv + 2], scen[i, nv + 3], scen[i, nv + 4])
                if (world) {
                    if (mode == "power") printf(" %9.4f", est_mean[i])
                    printf(" %6.0f", valid[i])
                }
                printf("\n")
            }
            else printf("   infeasible: %s\n", substr(notes[i], 1, 70))
        }
        if (nmax_i != .)
            printf("{txt}  n_max: n = %f reaches %2.0f%% power in every scenario; set by the least favorable one (row %f), achieved %5.3f (MCSE %5.3f).\n",
                best_nr.n, 100 * nmax_power, nmax_i, best_nr.power, best_nr.mcse)
        if (any(scen[ok_idx, nv + 4] :> 0.01))
            printf("{txt}  note: separation in more than 1%% of simulated studies in %f scenario(s) (column separated; separation(%s)).\n",
                sum(scen[ok_idx, nv + 4] :> 0.01), separation)
    }
    out = scen
    names = vn, (mode == "mde" ? "mde" : "implied_effect"), "power", "mcse", "separated", "feasible"
    if (world) {
        // which scenarios the analysis is valid in, and where its estimates
        // centre (power mode); the columns appear only when trends is varied
        if (mode == "power") {
            out = out, est_mean
            names = names, "est_mean"
        }
        out = out, valid
        names = names, "valid"
    }
    st_matrix("pape_r_scenarios", out)
    st_matrixcolstripe("pape_r_scenarios", (J(cols(names), 1, ""), names'))
    pape_rstr("fit", fitlab)
    pape_rstr("engine", eng)
    pape_rstr("effects", eff)
    pape_rstr("analysis", dgp.analysis)
    pape_rnum("worst", worst_i)
    pape_rnum("nmax", (nmax_i != . ? best_nr.n : .))
    pape_rnum("nmax_power", (nmax_i != . ? best_nr.power : .))
    pape_rnum("nmax_mcse", (nmax_i != . ? best_nr.mcse : .))
    pape_rnum("nmax_scenario", nmax_i)
    if (rows(cand_tab) > 0) {
        st_matrix("pape_r_nmaxcand", cand_tab)
        st_matrixcolstripe("pape_r_nmaxcand", (J(nv + 2, 1, ""), (vn' \ "criterion" \ "n")))
        pape_rstr("has_cand", "1")
    }
    else pape_rstr("has_cand", "0")
    if (nmax_i != .) pape_store_put("RES", resname + "_nmax", best_nr)
    for (i = 1; i <= G; i++) pape_rstr(sprintf("note%f", i), notes[i])
}

end
