*! powerape_cmd.mata -- command layer: parses the ado's option strings into
*! structures, runs the engine, stores results, prints, and posts r().
version 16
mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// Variable specs from parallel space-separated strings
// ---------------------------------------------------------------------------

struct pape_var colvector pape_vars_parse(string scalar names, string scalar types,
    string scalar ps, string scalar means, string scalar sds, string scalar iccs)
{
    struct pape_var colvector v
    string rowvector nm, ty
    real rowvector p, mu, sd, icc
    real scalar j, K
    nm = tokens(names); ty = tokens(types)
    K = cols(nm)
    if (K == 0) return(pape_var(0, 1))
    p = strtoreal(tokens(ps)); mu = strtoreal(tokens(means))
    sd = strtoreal(tokens(sds)); icc = strtoreal(tokens(iccs))
    v = pape_var(K, 1)
    for (j = 1; j <= K; j++) v[j] = pape_var_new(nm[j], ty[j], p[j], mu[j], sd[j], icc[j])
    return(v)
}

real colvector pape_num_col(string scalar s)
{
    if (strtrim(s) == "") return(J(0, 1, .))
    return(strtoreal(tokens(s))')
}

// ---------------------------------------------------------------------------
// powerape dgp (parametric / panel / iv / empirical / fromfit)
// ---------------------------------------------------------------------------

void pape_cmd_dgp_par(string scalar name, string scalar model, string scalar vnames,
    string scalar vtypes, string scalar vp, string scalar vmean, string scalar vsd,
    string scalar vicc, real scalar has_m, real scalar corr_kind, real scalar corr_scalar,
    string scalar corr_matname, real scalar baseline, real scalar signal,
    string scalar gamma_str, real scalar n_int, real scalar seed_int,
    string scalar spec_opts, real scalar display)
{
    struct pape_dgp scalar d
    real matrix Cm
    real colvector gamma
    real scalar gamma_given
    Cm = (corr_kind == 2 ? st_matrix(corr_matname) : J(0, 0, .))
    gamma = pape_num_col(gamma_str)
    gamma_given = (rows(gamma) > 0)
    d = pape_dgp_par(name, model, pape_vars_parse(vnames, vtypes, vp, vmean, vsd, vicc), has_m,
        corr_kind, corr_scalar, Cm, baseline, signal, gamma_given, gamma, n_int, seed_int)
    d.spec_opts = spec_opts
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

void pape_cmd_dgp_panel(string scalar name, string scalar model, string scalar vnames,
    string scalar vtypes, string scalar vp, string scalar vmean, string scalar vsd,
    string scalar vicc, real scalar has_m, string scalar nperiods, string scalar pperiods,
    real scalar retention, real scalar rho, real scalar cre_share, real scalar corr_kind,
    real scalar corr_scalar, string scalar corr_matname, real scalar baseline,
    real scalar signal, real scalar n_int, real scalar seed_int, string scalar spec_opts,
    real scalar display)
{
    struct pape_dgp scalar d
    real matrix Cm
    Cm = (corr_kind == 2 ? st_matrix(corr_matname) : J(0, 0, .))
    d = pape_dgp_panel(name, model, pape_vars_parse(vnames, vtypes, vp, vmean, vsd, vicc), has_m,
        pape_num_col(nperiods)', pape_num_col(pperiods)', retention, rho, cre_share,
        corr_kind, corr_scalar, Cm, baseline, signal, n_int, seed_int)
    d.spec_opts = spec_opts
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

void pape_cmd_dgp_iv(string scalar name, string scalar vnames, string scalar vtypes,
    string scalar vp, string scalar vmean, string scalar vsd, string scalar vicc,
    real scalar has_m, string scalar inames, string scalar itypes, string scalar ip,
    string scalar imean, string scalar isd, real scalar endogeneity, real scalar iv_strength,
    real scalar corr_kind, real scalar corr_scalar, string scalar corr_matname,
    real scalar baseline, real scalar signal, real scalar n_int, real scalar seed_int,
    string scalar spec_opts, real scalar display)
{
    struct pape_dgp scalar d
    real matrix Cm
    string scalar iicc
    real scalar j
    Cm = (corr_kind == 2 ? st_matrix(corr_matname) : J(0, 0, .))
    iicc = ""
    for (j = 1; j <= cols(tokens(inames)); j++) iicc = iicc + " ."
    d = pape_dgp_iv(name, pape_vars_parse(vnames, vtypes, vp, vmean, vsd, vicc), has_m,
        pape_vars_parse(inames, itypes, ip, imean, isd, iicc), endogeneity, iv_strength,
        corr_kind, corr_scalar, Cm, baseline, signal, n_int, seed_int)
    d.spec_opts = spec_opts
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

// empirical: varlist = pilot columns (focal first when joint); focal_joint
// 1 => first column is the focal; else fspec describes an independent focal
void pape_cmd_dgp_emp(string scalar name, string scalar model, string scalar varlist,
    string scalar touse, real scalar focal_joint, string scalar fname, string scalar ftype,
    real scalar fp, real scalar fmean, real scalar fsd, real scalar baseline,
    real scalar signal, string scalar gamma_str, real scalar seed_int,
    string scalar spec_opts, real scalar display)
{
    struct pape_dgp scalar d
    struct pape_var scalar fv
    real matrix M
    real colvector gamma
    M = st_data(., tokens(varlist), touse)
    gamma = pape_num_col(gamma_str)
    if (focal_joint) fv = pape_var_new("none", "normal", ., 0, 1, .)
    else fv = pape_var_new(fname, ftype, fp, fmean, fsd, .)
    d = pape_dgp_emp(name, model, M, tokens(varlist)', (focal_joint ? 1 : 0), fv,
        baseline, signal, rows(gamma) > 0, gamma, seed_int)
    d.spec_opts = spec_opts
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

void pape_cmd_dgp_fromfit(string scalar name, string scalar focal, real scalar baseline,
    string scalar touse, string scalar spec_opts, real scalar display)
{
    struct pape_dgp scalar d
    string matrix stripe
    string scalar cmd, link, depvar
    string colvector vnames
    real rowvector b
    real matrix X
    real scalar j, k, fcol, b0
    cmd = st_global("e(cmd)")
    if (cmd != "probit" & cmd != "logit")
        pape_err("powerape dgp, fromfit needs the results of a probit or logit fit in e().")
    link = cmd
    b = st_matrix("e(b)")
    stripe = st_matrixcolstripe("e(b)")
    k = cols(b) - 1
    if (stripe[cols(b), 2] != "_cons") pape_err("The pilot model must include an intercept.")
    vnames = J(k, 1, "")
    fcol = 0
    for (j = 1; j <= k; j++) {
        vnames[j] = stripe[j, 2]
        if (strpos(vnames[j], ".") | strpos(vnames[j], "#") | strpos(vnames[j], "("))
            pape_err(sprintf("The pilot-model route needs plain numeric regressors; refit without factor-variable or operator notation (found %s).", vnames[j]))
        if (vnames[j] == focal) fcol = j
    }
    if (fcol == 0) pape_err(sprintf("Focal variable %s is not a regressor of the pilot model.", focal))
    if (any(b :== .)) pape_err("The pilot fit has missing coefficients; fix the pilot model first.")
    if (st_global("e(wtype)") != "")
        printf("{err}warning: weights in the pilot fit are ignored; rows are resampled unweighted.\n")
    X = st_data(., vnames', touse)
    if (hasmissing(X)) pape_err("The pilot sample contains missing regressor values.")
    b0 = b[k + 1]
    depvar = st_global("e(depvar)")
    d = pape_dgp_fromfit(name, link, X, vnames, fcol, b[1..k]', b0, baseline,
        sprintf("%s %s %s", cmd, depvar, invtokens(vnames')))
    d.spec_opts = spec_opts
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

void pape_dgp_rreturn(struct pape_dgp scalar d)
{
    pape_rstr("model", d.model)
    pape_rstr("route", d.route)
    pape_rstr("estimand", d.estimand)
    pape_rstr("focal", d.focal.name)
    pape_rnum("baseline", d.baseline)
    pape_rnum("signal", d.signal_implied)
    pape_rnum("beta0", d.beta0)
    pape_rnum("beta_focal", d.beta_focal)
    pape_rnum("beta_mod", d.beta_mod)
    pape_rnum("beta_int", d.beta_int)
    pape_rnum("target", d.target_est)
    pape_rnum("p1", d.p1)
    pape_rnum("k", d.k)
    st_matrix("pape_r_gamma", d.gamma')
}

// ---------------------------------------------------------------------------
// setape / setaie / true values / describe / list / drop
// ---------------------------------------------------------------------------

void pape_cmd_setape(string scalar name, real scalar target, real scalar display)
{
    struct pape_dgp scalar d
    string scalar msg
    d = pape_get_dgp(name)
    if (pape_set_ape(d, target, msg)) pape_err(msg)
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

void pape_cmd_setaie(string scalar name, real scalar target, real scalar mf,
    real scalar mm, real scalar display)
{
    struct pape_dgp scalar d
    string scalar msg
    d = pape_get_dgp(name)
    if (pape_set_aie(d, target, mf, mm, msg)) pape_err(msg)
    pape_put_dgp(d)
    if (display) pape_print_dgp(d)
    pape_dgp_rreturn(d)
}

void pape_cmd_true(string scalar name, string scalar which)
{
    struct pape_dgp scalar d
    real scalar v
    d = pape_get_dgp(name)
    if (which == "ape") v = pape_true_ape(d)
    else v = pape_true_aie(d)
    printf("{txt}true %s of DGP %s on the integration set: {res}%12.0g\n", strupper(which), name, v)
    pape_rnum("value", v)
}

void pape_cmd_describe(string scalar name)
{
    pape_print_dgp(pape_get_dgp(name))
    pape_dgp_rreturn(pape_get_dgp(name))
}

void pape_cmd_list()
{
    string colvector k
    real scalar j
    k = pape_store_keys("DGP")
    printf("{txt}DGPs: ")
    if (rows(k) == 0) printf("(none)\n")
    else {
        for (j = 1; j <= rows(k); j++) printf("%s%s", (j > 1 ? ", " : ""), k[j])
        printf("\n")
    }
    pape_rstr("dgps", invtokens(k'))
    k = pape_store_keys("RES")
    printf("{txt}results: ")
    if (rows(k) == 0) printf("(none)\n")
    else {
        for (j = 1; j <= rows(k); j++) printf("%s%s", (j > 1 ? ", " : ""), k[j])
        printf("\n")
    }
    pape_rstr("results", invtokens(k'))
}

void pape_cmd_drop(string scalar name)
{
    string colvector k
    real scalar j
    if (name == "_all") {
        k = pape_store_keys("DGP")
        for (j = 1; j <= rows(k); j++) pape_store_drop("DGP", k[j])
        k = pape_store_keys("RES")
        for (j = 1; j <= rows(k); j++) pape_store_drop("RES", k[j])
        k = pape_store_keys("ALPHA")
        for (j = 1; j <= rows(k); j++) pape_store_drop("ALPHA", k[j])
        return
    }
    if (pape_store_has("DGP", name)) pape_store_drop("DGP", name)
    else if (pape_store_has("RES", name)) pape_store_drop("RES", name)
    else if (pape_store_has("ALPHA", name)) pape_store_drop("ALPHA", name)
    else pape_err(sprintf("Nothing named %s is stored.", name))
}

// ---------------------------------------------------------------------------
// Results: store, fetch, post to r()
// ---------------------------------------------------------------------------

struct pape_res scalar pape_get_res(string scalar name)
{
    struct pape_res scalar r
    if (!pape_store_has("RES", name))
        pape_err(sprintf("Result '%s' not found; run powerape power / n / mde first.", name))
    r = asarray(pape_store("RES"), name)
    return(r)
}

void pape_res_rreturn(struct pape_res scalar r, string scalar resname)
{
    real scalar j
    string scalar s
    pape_rstr("name", resname)
    pape_rstr("kind", r.kind)
    pape_rstr("claim", r.claim)
    pape_rstr("estimand", r.estimand)
    pape_rstr("model", r.model)
    pape_rstr("se", r.se)
    pape_rstr("separation", r.separation)
    pape_rstr("direction", r.direction)
    pape_rstr("dgp", r.dgp.name)
    pape_rstr("alpharule", r.alpha_rule)
    pape_rstr("engine", r.engine)
    pape_rstr("effects", r.effects)
    pape_rstr("fit", r.fit)
    pape_rnum("power", r.power)
    pape_rnum("mcse", r.mcse)
    pape_rnum("n", r.n)
    pape_rnum("nsim", r.nsim)
    pape_rnum("alpha", r.alpha)
    pape_rnum("conf", r.conf)
    pape_rnum("sesoi", r.sesoi)
    pape_rnum("target", r.target)
    pape_rnum("separated", r.separated)
    pape_rnum("min_cell", r.min_cell)
    pape_rnum("n_failed", r.n_failed)
    pape_rnum("goal", r.goal)
    pape_rnum("confirmed", r.confirmed)
    pape_rnum("nsim_confirm", r.nsim_confirm)
    pape_rnum("mde", r.mde)
    pape_rnum("is_margin", r.is_margin)
    if (cols(r.outcomes) > 0) {
        st_matrix("pape_r_outcomes", r.outcomes)
        st_matrixcolstripe("pape_r_outcomes", (J(cols(r.outcomes), 1, ""), r.outnames'))
        pape_rstr("has_outcomes", "1")
    }
    else pape_rstr("has_outcomes", "0")
    if (rows(r.history) > 0) {
        st_matrix("pape_r_history", r.history)
        s = (r.kind == "n" ? "n" : "value")
        st_matrixcolstripe("pape_r_history", (J(6, 1, ""), ("stage" \ "iter" \ s \ "alpha" \ "power" \ "mcse")))
        pape_rstr("has_history", "1")
    }
    else pape_rstr("has_history", "0")
    if (rows(r.results) > 0) {
        st_matrix("pape_r_results", r.results)
        st_matrixcolstripe("pape_r_results", (J(7, 1, ""), ("n" \ "alpha" \ "power" \ "mcse" \ "failed" \ "separated" \ "min_cell")))
        pape_rstr("has_results", "1")
    }
    else pape_rstr("has_results", "0")
}

// ---------------------------------------------------------------------------
// power / curve / n / mde
// ---------------------------------------------------------------------------

void pape_cmd_power(string scalar dname, real scalar n, string scalar claim,
    real scalar sesoi, real scalar alpha, real scalar conf, string scalar rule,
    real scalar nsim, real scalar seed, string scalar se, string scalar separation,
    real scalar keep_draws, string scalar resname, real scalar display)
{
    struct pape_dgp scalar d
    struct pape_res scalar r
    d = pape_get_dgp(dname)
    r = pape_power_once(d, n, claim, sesoi, conf, nsim, seed, 1, se, alpha, rule, keep_draws, separation)
    (void) pape_sparse_warning(r.separated, r.min_cell, separation, "", r.estimand == "aie")
    pape_store_put("RES", resname, r)
    if (display) pape_print_power(r)
    pape_res_rreturn(r, resname)
}

void pape_cmd_curve(string scalar dname, string scalar nlist, string scalar claim,
    real scalar sesoi, real scalar alpha, real scalar conf, string scalar rule,
    real scalar nsim, real scalar seed, string scalar se, string scalar separation,
    string scalar resname, real scalar display)
{
    struct pape_dgp scalar d
    struct pape_res scalar r
    d = pape_get_dgp(dname)
    r = pape_curve(d, pape_num_col(nlist), claim, sesoi, alpha, conf, rule, nsim, seed, se, separation)
    pape_store_put("RES", resname, r)
    if (display) pape_print_curve(r)
    pape_res_rreturn(r, resname)
}

void pape_cmd_n(string scalar dname, real scalar power, string scalar claim,
    real scalar sesoi, real scalar alpha, real scalar conf, string scalar rule,
    real scalar nsim, real scalar seed, real scalar n_lo, real scalar n_hi,
    real scalar max_iter, real scalar confirm, real scalar nsim_confirm,
    string scalar se, string scalar separation, string scalar resname, real scalar display)
{
    struct pape_dgp scalar d
    struct pape_res scalar r
    d = pape_get_dgp(dname)
    r = pape_n(d, power, claim, sesoi, alpha, conf, rule, nsim, seed, n_lo, n_hi, max_iter,
        confirm, nsim_confirm, se, separation)
    pape_store_put("RES", resname, r)
    if (display) pape_print_n(r)
    pape_res_rreturn(r, resname)
}

void pape_cmd_mde(string scalar dname, real scalar n, real scalar power, string scalar claim,
    real scalar sesoi, real scalar alpha, real scalar conf, string scalar rule,
    real scalar nsim, real scalar seed, real scalar mf, real scalar mm, real scalar max_iter,
    real scalar confirm, real scalar nsim_confirm, string scalar se, string scalar direction,
    string scalar separation, string scalar resname, real scalar display)
{
    struct pape_dgp scalar d
    struct pape_res scalar r
    d = pape_get_dgp(dname)
    r = pape_mde(d, n, power, claim, sesoi, alpha, conf, rule, nsim, seed, mf, mm, max_iter,
        confirm, nsim_confirm, se, direction, separation)
    pape_store_put("RES", resname, r)
    if (display) pape_print_mde(r)
    pape_res_rreturn(r, resname)
}

// ---------------------------------------------------------------------------
// statement / alpha / powerat / frontier
// ---------------------------------------------------------------------------

void pape_cmd_statement(string scalar resname)
{
    struct pape_res scalar r
    struct pape_alpha scalar a
    string scalar s
    if (pape_store_has("RES", resname)) {
        r = pape_get_res(resname)
        if (r.kind == "curve") pape_err("powerape statement needs a power, n, mde, or alpha result.")
        s = pape_statement(r)
    }
    else if (pape_store_has("ALPHA", resname)) {
        a = asarray(pape_store("ALPHA"), resname)
        s = pape_alpha_statement(a)
    }
    else pape_err(sprintf("Result '%s' not found.", resname))
    pape_print_wrapped(s)
    pape_rstr("statement", s)
}

void pape_cmd_alpha(string scalar resname, real scalar cost, real scalar prior,
    string scalar error, real scalar cap, real scalar size, real scalar nsim_size,
    real scalar seed, string scalar outname, real scalar display)
{
    struct pape_res scalar x
    struct pape_alpha scalar a
    x = pape_get_res(resname)
    a = pape_alpha_fn(x, cost, prior, error, cap, size, nsim_size, seed)
    pape_store_put("ALPHA", outname, a)
    if (display) pape_print_alpha(a)
    pape_rnum("alpha", a.alpha)
    pape_rnum("power", a.power)
    pape_rnum("mcse", a.mcse)
    pape_rnum("beta", a.beta)
    pape_rnum("wcer", a.wcer)
    pape_rnum("power_default", a.power_default)
    pape_rnum("wcer_default", a.wcer_default)
    pape_rnum("range_lo", a.range[1])
    pape_rnum("range_hi", a.range[2])
    pape_rnum("capped", a.capped)
    pape_rnum("cap", a.cap)
    pape_rnum("alpha_uncapped", a.alpha_uncapped)
    pape_rnum("wcer_uncapped", a.wcer_uncapped)
    pape_rnum("size", (a.has_size ? a.size : .))
    pape_rnum("size_mcse", (a.has_size ? a.size_mcse : .))
    pape_rnum("size_truth", (a.has_size ? a.size_truth : .))
    st_matrix("pape_r_grid", a.grid)
    st_matrixcolstripe("pape_r_grid", (J(4, 1, ""), ("alpha" \ "power" \ "beta" \ "wcer")))
}

void pape_cmd_powerat(string scalar resname, string scalar alphas)
{
    struct pape_res scalar x
    real matrix out
    real colvector a
    real scalar j
    string colvector nm
    x = pape_get_res(resname)
    if (x.kind != "power") pape_err("powerape powerat needs a powerape power result.")
    if (!x.has_draws) pape_err("The result carries no stored draws; rerun powerape power without nodraws.")
    a = pape_num_col(alphas)
    if (rows(a) == 0 | any(a :<= 0) | any(a :>= 0.5)) pape_err("alpha() values must lie in (0, 0.5).")
    out = pape_power_at_draws(x.draws, x.target, x.claim, x.sesoi, x.nsim, a)
    nm = "alpha" \ "conf" \ "power" \ "mcse"
    for (j = 1; j <= cols(x.outnames); j++) nm = nm \ x.outnames[j]
    st_matrix("pape_r_table", out)
    st_matrixcolstripe("pape_r_table", (J(rows(nm), 1, ""), nm))
    printf("{txt}power at other error rates (re-thresholding the stored draws of %s)\n", resname)
    printf("{txt}%9s %9s %8s %8s\n", "alpha", "conf", "power", "mcse")
    for (j = 1; j <= rows(out); j++)
        printf("{res}%9.4f %9.4f %8.4f %8.4f\n", out[j, 1], out[j, 2], out[j, 3], out[j, 4])
}

void pape_cmd_frontier(real scalar cost, real scalar power, real scalar prior, string scalar claim)
{
    real scalar a
    a = pape_alpha_frontier(cost, power, prior, claim)
    printf("{txt}frontier alpha (%s claim, cost %g, power %4.2f, prior %g) = {res}%8.4f\n", claim, cost, power, prior, a)
    pape_rnum("alpha", a)
}

// ---------------------------------------------------------------------------
// Robustness sweep (robust.R: ape_robust); DGP rebuilds go through the ado
// ---------------------------------------------------------------------------

string scalar pape_last_err()
{
    pointer(transmorphic) scalar p
    p = findexternal("PAPE_LASTERR")
    if (p == NULL) return("")
    return(*p)
}

void pape_set_quiet(real scalar q)
{
    pointer(transmorphic) scalar p
    p = findexternal("PAPE_QUIET")
    if (p == NULL) p = crexternal("PAPE_QUIET")
    *p = q
}

// replace name(...) in an option string (numeric scalar options only)
string scalar pape_opt_replace(string scalar opts, string scalar nm, real scalar val)
{
    real scalar i, j, depth
    string scalar out, s
    s = opts
    i = strpos(s, nm + "(")
    // find a token boundary before nm
    while (i > 0) {
        if (i == 1) break
        if (strpos(" ,", substr(s, i - 1, 1))) break
        j = strpos(substr(s, i + 1, .), nm + "(")
        i = (j > 0 ? i + j : 0)
    }
    if (i == 0) return(s + " " + nm + "(" + strofreal(val, "%21.0g") + ")")
    depth = 0
    for (j = i + strlen(nm); j <= strlen(s); j++) {
        if (substr(s, j, 1) == "(") depth++
        else if (substr(s, j, 1) == ")") {
            depth--
            if (depth == 0) break
        }
    }
    out = substr(s, 1, i - 1) + nm + "(" + strofreal(val, "%21.0g") + ")" + substr(s, j + 1, .)
    return(out)
}

real scalar pape_claim_distance(string scalar claim, real scalar t, real scalar sesoi)
{
    if (claim == "detect") return(abs(t))
    if (claim == "minimum") return(abs(t) - sesoi)
    return(sesoi - abs(t))
}

void pape_cmd_robust(string scalar dname, real scalar n, string scalar claim,
    real scalar sesoi, real scalar alpha, real scalar conf, string scalar rule,
    real scalar nsim, real scalar seed, string scalar varynames, string scalar varyvals,
    string scalar pin, string scalar mode, real scalar power, real scalar grid_points,
    real scalar nmax, real scalar nmax_power, string scalar se, string scalar direction,
    string scalar separation, string scalar resname, real scalar display)
{
    struct pape_dgp scalar dgp, d_i
    struct pape_res scalar pw, m_i, nr, best_nr
    struct pape_sim scalar sim
    real rowvector lv
    string rowvector vn, vv, allowed, notes
    string scalar opts, msg, dir, se_use, where, fitlab, eng, eff
    real matrix levels, grid, scen, cand_tab
    real colvector key, crit, ok_idx, cr, o, cand, nn
    real scalar nv, j, i, G, rc, is_aie, rule_on, worst_i, nmax_i, k, sk, dist, best_j
    real colvector nl
    real matrix L
    dgp = pape_get_dgp(dname)
    se_use = se
    if (se != "model" & (dgp.route == "panel" | dgp.route == "iv")) {
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
    if (mode == "mde" & dgp.has_m & (dgp.main_focal == . | dgp.main_moderator == .))
        pape_err("MDE mode for an AIE design needs the DGP pinned with powerape setaie first.")
    dir = direction
    if (dir == "") dir = (dgp.has_effect & dgp.target_est < 0 ? "negative" : "positive")
    vn = tokens(varynames)
    nv = cols(vn)
    if (nv == 0) pape_err("vary() must name at least one contextual input.")
    if (dgp.builder == "from_fit") allowed = ("baseline")
    else if (dgp.builder == "empirical") allowed = ("baseline", "signal")
    else if (dgp.builder == "panel") allowed = ("baseline", "signal", "rho", "retention")
    else if (dgp.builder == "iv") allowed = ("baseline", "signal", "endogeneity", "ivstrength")
    else allowed = ("baseline", "signal", "correlation")
    for (j = 1; j <= nv; j++) {
        if (!anyof(allowed, vn[j]))
            pape_err(sprintf("Cannot vary %s for this DGP route; varyable inputs: %s.", vn[j], invtokens(allowed, ", ")))
    }
    if (anyof(vn, "signal") & dgp.gamma_given) pape_err("Cannot vary signal when the DGP was built with explicit gamma().")
    if (anyof(vn, "correlation") & dgp.corr_given == 2) pape_err("Varying correlation needs a scalar correlation in the base DGP.")
    // expand levels
    vv = tokens(varyvals, "|")
    vv = select(vv, vv :!= "|")
    if (cols(vv) != nv) pape_err("vary(): internal parse error.")
    // level lists kept in a padded matrix (one row per input; nl = length)
    nl = J(nv, 1, .)
    L = J(nv, 0, .)
    G = 1
    for (j = 1; j <= nv; j++) {
        key = pape_num_col(vv[j])
        if (rows(key) < 2) pape_err(sprintf("vary(%s()) needs at least two values.", vn[j]))
        if (rows(key) == 2) key = rangen(key[1], key[2], grid_points)
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
    is_aie = (dgp.estimand == "aie")
    // scenario table: inputs, implied effect or mde, power, mcse, separated, feasible
    scen = J(G, nv + 5, .)
    notes = J(1, G, "")
    crit = J(G, 1, .)
    pape_set_quiet(1)
    for (i = 1; i <= G; i++) {
        opts = dgp.spec_opts
        for (j = 1; j <= nv; j++) opts = pape_opt_replace(opts, vn[j], grid[i, j])
        scen[i, 1..nv] = grid[i, .]
        if (dgp.route == "empirical") {
            // pilot-row routes recalibrate from the stored rows (no data, no e() needed)
            d_i = pape_recal_emp(dgp, (anyof(vn, "baseline") ? grid[i, selectindex(vn :== "baseline")[1]] : dgp.baseline),
                (anyof(vn, "signal") ? grid[i, selectindex(vn :== "signal")[1]] : .))
        }
        else {
            rc = _stata("capture quietly powerape dgp _pape_scn, " + opts + " nodisplay replace", 1)
            if (rc) {
                notes[i] = pape_last_err()
                scen[i, nv + 5] = 0
                continue
            }
            d_i = pape_get_dgp("_pape_scn")
        }
        if (mode == "mde") {
            // MDE per scenario; errors are reported as infeasible rows
            m_i = pape_mde_safe(d_i, n, power, claim, sesoi, alpha, conf, rule, nsim,
                (seed == . ? . : seed + i), dgp.main_focal, dgp.main_moderator, 6, 0, 4 * nsim,
                se_use, dir, separation, rc, msg)
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
        if (pin == "ape") {
            if (is_aie) rc = pape_set_aie(d_i, dgp.target_est, dgp.main_focal, dgp.main_moderator, msg)
            else rc = pape_set_ape(d_i, dgp.target_est, msg)
            if (rc) {
                notes[i] = msg
                scen[i, nv + 5] = 0
                continue
            }
        }
        else {
            d_i.beta_focal = dgp.beta_focal; d_i.beta_mod = dgp.beta_mod; d_i.beta_int = dgp.beta_int
            d_i.estimand = dgp.estimand; d_i.target_est = dgp.target_est
            d_i.main_focal = dgp.main_focal; d_i.main_moderator = dgp.main_moderator
            d_i.has_effect = 1
            d_i.target_est = (is_aie ? pape_true_aie(d_i) : pape_true_ape(d_i))
        }
        pw = pape_power_once(d_i, n, claim, sesoi, conf, nsim, (seed == . ? . : seed + i), 0, se_use, ., "", 1, separation)
        scen[i, nv + 1] = d_i.target_est; scen[i, nv + 2] = pw.power
        scen[i, nv + 3] = pw.mcse; scen[i, nv + 4] = pw.separated
        scen[i, nv + 5] = 1
        d_i.name = sprintf("_pape_scn%f", i)
        pape_put_dgp(d_i)
        dist = pape_claim_distance(claim, d_i.target_est, sesoi)
        if (any(pw.draws.ok)) crit[i] = (dist > 0 ? mean(select(pw.draws.se, pw.draws.ok)) * sqrt(n) / dist : .)
        if (dist <= 0 & any(pw.draws.ok)) crit[i] = 1e300
    }
    pape_set_quiet(0)
    key = (mode == "mde" ? abs(scen[., nv + 1]) : scen[., nv + 2])
    ok_idx = selectindex(key :!= .)
    if (rows(ok_idx) == 0 & cols(ok_idx) == 0) pape_err("No scenario was feasible; widen or shift the vary() ranges.")
    ok_idx = vec(ok_idx)
    if (mode == "mde") worst_i = ok_idx[selectindex(key[ok_idx] :== max(key[ok_idx]))[1]]
    else {
        cr = crit[ok_idx]
        cr = editmissing(cr, -1e300)
        o = order((scen[ok_idx, nv + 2], -cr), (1, 2))
        worst_i = ok_idx[o[1]]
    }
    // n_max
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
    // print
    if (display) {
        if (mode == "mde") {
            printf("{txt}powerape robustness sweep -- minimum detectable %s, %s claim%s\n",
                strupper(is_aie ? "aie" : "ape"), claim, (dir == "negative" ? ", decreases" : ""))
            printf("  n = %f, %2.0f%% target power, nsim = %f per scenario, %f scenario(s) over: %s\n",
                n, 100 * power, nsim, G, invtokens(vn, ", "))
            printf("  MDE range [%9.4f, %9.4f]; worst (largest-MDE) scenario: row %f\n",
                min(scen[ok_idx, nv + 1]), max(scen[ok_idx, nv + 1]), worst_i)
        }
        else {
            printf("{txt}powerape robustness sweep -- %s claim, %s, pin = %s\n", claim, strupper(is_aie ? "aie" : "ape"), pin)
            printf("  n = %f, nsim = %f per scenario, %f scenario(s) over: %s\n", n, nsim, G, invtokens(vn, ", "))
            printf("  power range [%5.3f, %5.3f]; worst scenario: row %f\n",
                min(scen[ok_idx, nv + 2]), max(scen[ok_idx, nv + 2]), worst_i)
        }
        printf("  fit: %s\n", fitlab)
        printf("{txt}%4s", "row")
        for (j = 1; j <= nv; j++) printf(" %12s", vn[j])
        printf(" %14s %8s %8s %9s\n", (mode == "mde" ? "mde" : "implied_eff"), "power", "mcse", "separated")
        for (i = 1; i <= G; i++) {
            printf("{res}%4.0f", i)
            for (j = 1; j <= nv; j++) printf(" %12.4g", scen[i, j])
            if (scen[i, nv + 5] == 1)
                printf(" %14.4f %8.4f %8.4f %9.4f\n", scen[i, nv + 1], scen[i, nv + 2], scen[i, nv + 3], scen[i, nv + 4])
            else printf("   infeasible: %s\n", substr(notes[i], 1, 70))
        }
        if (nmax_i != .)
            printf("{txt}  n_max: n = %f reaches %2.0f%% power in every scenario; set by the least favorable one (row %f), achieved %5.3f (MCSE %5.3f).\n",
                best_nr.n, 100 * nmax_power, nmax_i, best_nr.power, best_nr.mcse)
        if (any(scen[ok_idx, nv + 4] :> 0.01))
            printf("{txt}  note: separation in more than 1%% of simulated studies in %f scenario(s) (column separated; separation(%s)).\n",
                sum(scen[ok_idx, nv + 4] :> 0.01), separation)
    }
    st_matrix("pape_r_scenarios", scen)
    st_matrixcolstripe("pape_r_scenarios", (J(nv + 5, 1, ""), (vn' \ (mode == "mde" ? "mde" : "implied_effect") \ "power" \ "mcse" \ "separated" \ "feasible")))
    pape_rstr("fit", fitlab)
    pape_rstr("engine", eng)
    pape_rstr("effects", eff)
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

// MDE with errors returned instead of raised (for the sweep)
struct pape_res scalar pape_mde_safe(struct pape_dgp scalar d, real scalar n, real scalar power,
    string scalar claim, real scalar sesoi, real scalar alpha, real scalar conf, string scalar rule,
    real scalar nsim, real scalar seed, real scalar mf, real scalar mm, real scalar max_iter,
    real scalar confirm, real scalar nsim_confirm, string scalar se, string scalar direction,
    string scalar separation, real scalar rc, string scalar msg)
{
    struct pape_res scalar r
    pointer(transmorphic) scalar p
    rc = 0; msg = ""
    // The searches raise errors through pape_err(); in quiet mode they are
    // recorded and the Stata-level capture around the sweep's scenario
    // rebuild does not apply here, so run under _stata with a wrapper.
    p = findexternal("PAPE_MDE_TMP")
    if (p == NULL) p = crexternal("PAPE_MDE_TMP")
    *p = J(0, 0, .)
    pape_store_put("DGP", "_pape_mde_tmp", d)
    rc = _stata("capture mata: pape_mde_wrap(" + pape_ns(n) + ", " + pape_ns(power) + ", " +
        pape_qs(claim) + ", " + pape_ns(sesoi) + ", " + pape_ns(alpha) + ", " + pape_ns(conf) + ", " +
        pape_qs(rule) + ", " + pape_ns(nsim) + ", " + pape_ns(seed) + ", " + pape_ns(mf) + ", " +
        pape_ns(mm) + ", " + pape_ns(max_iter) + ", " + pape_ns(confirm) + ", " + pape_ns(nsim_confirm) + ", " +
        pape_qs(se) + ", " + pape_qs(direction) + ", " + pape_qs(separation) + ")", 1)
    if (rc) {
        msg = pape_last_err()
        pape_store_drop("DGP", "_pape_mde_tmp")
        return(r)
    }
    r = asarray(pape_store("RES"), "_pape_mde_tmp")
    pape_store_drop("RES", "_pape_mde_tmp")
    pape_store_drop("DGP", "_pape_mde_tmp")
    return(r)
}

void pape_mde_wrap(real scalar n, real scalar power, string scalar claim, real scalar sesoi,
    real scalar alpha, real scalar conf, string scalar rule, real scalar nsim, real scalar seed,
    real scalar mf, real scalar mm, real scalar max_iter, real scalar confirm,
    real scalar nsim_confirm, string scalar se, string scalar direction, string scalar separation)
{
    struct pape_dgp scalar d
    struct pape_res scalar r
    d = pape_get_dgp("_pape_mde_tmp")
    r = pape_mde(d, n, power, claim, sesoi, alpha, conf, rule, nsim, seed, mf, mm, max_iter,
        confirm, nsim_confirm, se, direction, separation)
    pape_store_put("RES", "_pape_mde_tmp", r)
}

end

mata:
mata set matastrict on
// post results as Stata scalars / globals with the pape_r_ prefix
void pape_rnum(string scalar name, real scalar value) st_numscalar("pape_r_" + name, value)
void pape_rstr(string scalar name, string scalar value) st_global("pape_r_" + name, value)
end

mata:
mata set matastrict on
string scalar pape_qs(string scalar s) return(char(34) + s + char(34))
string scalar pape_ns(real scalar x) return(strofreal(x, "%21.0g"))
end
