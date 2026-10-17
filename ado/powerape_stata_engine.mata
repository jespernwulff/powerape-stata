*! powerape_stata_engine.mata -- fits by Stata's own estimation commands.
*! Each simulated study is written to a temporary frame and fit with
*! probit/logit (vce(), asis) or cfprobit; effects come from closed-form
*! delta-method formulas applied to e(b)/e(V) (equal to margins/ginteff to
*! numerical precision, see validation) or, with effects(margins), from
*! margins / ginteff themselves.
version 16
mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// Engine settings (Mata externals set by the ado)
// ---------------------------------------------------------------------------

void pape_set_engine(string scalar engine, string scalar effects)
{
    pointer(transmorphic) scalar p
    if (engine != "stata" & engine != "mata") pape_err("engine() must be stata or mata.")
    if (effects != "closed" & effects != "margins") pape_err("effects() must be closed or margins.")
    p = findexternal("PAPE_ENGINE")
    if (p == NULL) p = crexternal("PAPE_ENGINE")
    *p = engine
    p = findexternal("PAPE_EFFECTS")
    if (p == NULL) p = crexternal("PAPE_EFFECTS")
    *p = effects
}

string scalar pape_engine()
{
    pointer(transmorphic) scalar p
    p = findexternal("PAPE_ENGINE")
    if (p == NULL) return("stata")
    return(*p)
}

string scalar pape_effects()
{
    pointer(transmorphic) scalar p
    p = findexternal("PAPE_EFFECTS")
    if (p == NULL) return("closed")
    return(*p)
}

// Who fitted the simulated studies (0.1.2). pape_fit_lab() builds the label
// that the simulation commands print as "fit: ..." and store in r(fit);
// pape_stamp_fit() records it on a result together with r(engine) and
// r(effects). Both run while the result is built, inside the mata call of
// the command, so a per-command engine()/effects() is still in force (the
// ado restores the session default only after the call returns); the print
// routines read the stored fields and never consult the externals.
string scalar pape_fit_lab(struct pape_dgp scalar d, string scalar se,
    string scalar separation, string scalar engine, string scalar effects)
{
    string scalar cmd, opts
    if (engine == "mata") {
        if (d.route == "iv") return("Mata engine (Fisher scoring; control-function probit)")
        return("Mata engine (Fisher scoring)")
    }
    cmd = (d.route == "iv" ? "cfprobit" : d.link)
    opts = ""
    if (d.route == "panel") opts = " vce(cluster)"
    else if (d.route != "iv" & se == "robust") opts = " vce(robust)"
    if (separation == "keep") opts = opts + " asis"
    if (opts != "") cmd = cmd + "," + opts
    return(sprintf("%s via engine(%s), effects(%s)", cmd, engine, effects))
}

void pape_stamp_fit(struct pape_res scalar r)
{
    r.engine = pape_engine()
    r.effects = (r.engine == "stata" ? pape_effects() : "closed")
    r.fit = pape_fit_lab(r.dgp, r.se, r.separation, r.engine, r.effects)
}

// ---------------------------------------------------------------------------
// Temporary frame for the simulated studies
// ---------------------------------------------------------------------------

void pape_frame_open(real scalar n, real scalar k, real scalar extra)
{
    pointer(transmorphic) scalar p
    string rowvector names
    real scalar j
    p = findexternal("PAPE_USERFRAME")
    if (p == NULL) p = crexternal("PAPE_USERFRAME")
    *p = st_framecurrent()
    if (!st_frameexists("_pape_sim")) st_framecreate("_pape_sim")
    st_framecurrent("_pape_sim")
    stata("clear", 1)
    st_addobs(n)
    names = "pape_y"
    for (j = 1; j <= k; j++) names = names, sprintf("pape_x%f", j)
    names = names, "pape_id"
    for (j = 1; j <= extra; j++) names = names, sprintf("pape_w%f", j)
    (void) st_addvar("double", names)
}

void pape_frame_resize(real scalar n)
{
    real scalar cur
    cur = st_nobs()
    if (n > cur) st_addobs(n - cur)
    else if (n < cur) st_dropobsin((n + 1, cur))
}

void pape_frame_close()
{
    pointer(transmorphic) scalar p
    p = findexternal("PAPE_USERFRAME")
    if (p != NULL) {
        if (st_frameexists(*p)) st_framecurrent(*p)
    }
    if (st_frameexists("_pape_sim")) st_framedrop("_pape_sim")
}

// ---------------------------------------------------------------------------
// probit / logit by Stata (OIM, vce(robust), vce(cluster)); X includes the
// constant in column 1. Returns the fit in X's column order.
// ---------------------------------------------------------------------------

string scalar pape_xnames(real scalar k)
{
    string scalar s
    real scalar j
    s = ""
    for (j = 1; j <= k; j++) s = s + sprintf(" pape_x%f", j)
    return(s)
}

struct pape_fit scalar pape_fit_stata(real matrix X, real colvector y, real colvector id,
    string scalar link, string scalar vce, real scalar asis)
{
    struct pape_fit scalar f
    real colvector bs, perm, b, zero
    real matrix Vs
    real scalar k, n, rc, j, sdj, runaway, kept_sep, omitted
    string scalar cmd
    string colvector names
    n = rows(X); k = cols(X) - 1
    f.ok = 0; f.converged = 0; f.rank = .
    st_store(., "pape_y", y)
    if (k > 0) st_store(., tokens(pape_xnames(k)), X[., 2..(k + 1)])
    if (vce == "cluster") st_store(., "pape_id", id)
    cmd = "quietly " + link + " pape_y" + pape_xnames(k) + ", iterate(100)"
    if (vce == "robust") cmd = cmd + " vce(robust)"
    else if (vce == "cluster") cmd = cmd + " vce(cluster pape_id)"
    if (asis) cmd = cmd + " asis"
    rc = _stata(cmd, 1)
    if (rc) return(f)
    bs = st_matrix("e(b)")'
    Vs = st_matrix("e(V)")
    if (rows(bs) != k + 1) return(f)
    perm = (k + 1) \ (1::k)
    if (k == 0) perm = 1
    b = bs[perm]
    f.b = b
    f.V = Vs[perm, perm]
    f.converged = (st_numscalar("e(converged)") == 1)
    zero = (diagonal(f.V) :== 0)
    f.rank = (k + 1) - sum(zero)
    // separation(keep): Stata's asis fit returns a separated coefficient with
    // a zero variance row (its SE prints as missing). Keep it: the delta
    // method then collapses to the identified cells, which is the limit of
    // R's glm fit (the degenerate Wald interval of engine(mata)). Omitted
    // (collinear, o.) terms still fail, as they do without asis.
    kept_sep = 0
    if (asis & f.rank < k + 1) {
        names = st_matrixcolstripe("e(b)")[., 2]
        names = names[perm]
        omitted = 0
        for (j = 1; j <= k + 1; j++) {
            if (zero[j] & (strpos(names[j], "o.") > 0 | b[j] == 0)) omitted = 1
        }
        kept_sep = !omitted
    }
    f.eta = X * b
    f.mu = pape_G(link, f.eta)
    runaway = 0
    if (k > 0 & !hasmissing(b)) {
        for (j = 2; j <= k + 1; j++) {
            sdj = sqrt(variance(X[., j]))
            if (abs(b[j]) * max((sdj, 1e-12)) > 15) runaway = 1
        }
    }
    f.ok = (f.converged & (f.rank == k + 1 | kept_sep) & !hasmissing(b) & st_numscalar("e(N)") == n & !runaway)
    return(f)
}

// leading integer of a stripe label such as "2._at" or "1bn._at" (0 if none)
real scalar pape_leading_int(string scalar s)
{
    real scalar j
    j = 1
    while (j <= strlen(s) & strpos("0123456789", substr(s, j, 1)) > 0) j++
    if (j == 1) return(0)
    return(strtoreal(substr(s, 1, j - 1)))
}

// Factor-variable term for a stored column: i.name for binary, c.name otherwise
string scalar pape_fvterm(real scalar j, string scalar type)
    return((type == "binary" ? "i." : "c.") + sprintf("pape_x%f", j))

// Effects by margins / ginteff after a Stata fit with factor-variable
// notation. layout: "ape" (X = 1, d, Z...) or "aie" (X = 1, d, m, dm, Z...).
// Columns after the first block are plain continuous regressors. Returns 0
// when the effect is not estimable.
real scalar pape_effects_margins(struct pape_dgp scalar d, real matrix X,
    real colvector y, real colvector id, string scalar vce, real scalar asis,
    real scalar ape_out, real scalar se_out)
{
    real scalar k, n, rc, j, kk
    real rowvector b
    real matrix V
    string scalar cmd, rest, first
    n = rows(X); k = cols(X) - 1
    st_store(., "pape_y", y)
    st_store(., tokens(pape_xnames(k)), X[., 2..(k + 1)])
    if (vce == "cluster") st_store(., "pape_id", id)
    if (d.estimand == "aie") {
        first = pape_fvterm(1, d.focal.type) + "##" + pape_fvterm(2, d.moderator.type)
        rest = ""
        for (j = 4; j <= k; j++) rest = rest + sprintf(" pape_x%f", j)
    }
    else {
        first = pape_fvterm(1, d.focal.type)
        rest = ""
        for (j = 2; j <= k; j++) rest = rest + sprintf(" pape_x%f", j)
    }
    cmd = "quietly " + d.link + " pape_y " + first + rest + ", iterate(100)"
    if (vce == "robust") cmd = cmd + " vce(robust)"
    else if (vce == "cluster") cmd = cmd + " vce(cluster pape_id)"
    if (asis) cmd = cmd + " asis"
    rc = _stata(cmd, 1)
    if (rc) return(0)
    if (st_numscalar("e(converged)") != 1 | st_numscalar("e(N)") != n) return(0)
    if (d.estimand == "aie") rc = _stata("quietly ginteff, dydxs(pape_x1 pape_x2)", 1)
    else rc = _stata("quietly margins, dydx(pape_x1)", 1)
    if (rc) return(0)
    b = st_matrix("r(b)")
    V = st_matrix("r(V)")
    if (cols(b) == 0) return(0)
    kk = cols(b)
    ape_out = b[kk]
    se_out = sqrt(V[kk, kk])
    if (ape_out == . | se_out == . | se_out <= 0) return(0)
    return(1)
}

// ---------------------------------------------------------------------------
// cfprobit by Stata for the IV route; maps e(b)/e(V) onto the control-function
// design W = [1, d, m, d*m, Z, c, c*m] (powerape's layout).
// ---------------------------------------------------------------------------

real scalar pape_stripe_index(string matrix stripe, string scalar eq, string scalar name)
{
    real scalar j
    for (j = 1; j <= rows(stripe); j++) {
        if (stripe[j, 2] == name & (eq == "" | stripe[j, 1] == eq)) return(j)
    }
    return(0)
}

struct pape_cffit scalar pape_cf_fit_stata(struct pape_dgp scalar d, struct pape_xdraw scalar xx,
    real scalar asis, real scalar use_margins, real scalar ape_out, real scalar se_out)
{
    struct pape_cffit scalar f
    real colvector pi, q, c, dcdq, idxv, cvec
    real matrix Vs, W
    real rowvector bs
    string matrix stripe
    string scalar cmd, zs, inst, dn, dmn, mn, eqm, eqe, cfn, eqn
    real scalar n, k, nz, j, rc, has_m, kw, kk, i1, i2, at
    f.ok = 0
    n = rows(xx.y); k = d.k; nz = d.n_inst; has_m = d.has_m
    // stored columns: pape_x1 = d, pape_x2.. = Z (k), then m (if any) as pape_w1,
    // instruments as pape_w2.. (or pape_w1.. without a moderator)
    st_store(., "pape_y", xx.y)
    st_store(., "pape_x1", xx.d)
    if (k > 0) st_store(., tokens(pape_xnames(k + 1))[2..(k + 1)], xx.Zx)
    zs = ""
    for (j = 2; j <= k + 1; j++) zs = zs + sprintf(" pape_x%f", j)
    if (has_m) {
        st_store(., "pape_w1", xx.m)
        st_store(., tokens(pape_wnames(1 + nz))[2..(1 + nz)], xx.X1[., (k + 3)..(k + 2 + nz)])
        inst = ""
        for (j = 2; j <= 1 + nz; j++) inst = inst + sprintf(" pape_w%f c.pape_w%f#c.pape_w1", j, j)
        // the focal-by-moderator term enters the main equation only. A binary
        // focal (probit first stage) is a factor variable in cfprobit, so the
        // interaction is written in factor form too (1.d#c.m): the fit is the
        // same as with c.d#c.m, but margins can then differentiate d (a variable
        // may not be both factor and continuous).
        dmn = (d.first_stage == "probit" ? "1.pape_x1#c.pape_w1" : "c.pape_x1#c.pape_w1")
        cmd = "quietly cfprobit pape_y" + zs + " pape_w1 (pape_x1 =" + inst + ", " + d.first_stage +
            " interact(pape_w1)), mainonly(" + dmn + ") vce(robust) iterate(100)"
    }
    else {
        st_store(., tokens(pape_wnames(nz)), xx.X1[., (2 + k)..(1 + k + nz)])
        inst = ""
        for (j = 1; j <= nz; j++) inst = inst + sprintf(" pape_w%f", j)
        cmd = "quietly cfprobit pape_y" + zs + " (pape_x1 =" + inst + ", " + d.first_stage + "), vce(robust) iterate(100)"
    }
    if (asis) cmd = cmd + " asis"
    rc = _stata(cmd, 1)
    if (rc) return(f)
    if (st_numscalar("e(converged)") != 1 | st_numscalar("e(N)") != n) return(f)
    if (use_margins & !has_m) {
        rc = _stata("quietly margins, dydx(pape_x1)", 1)
        if (rc) return(f)
        bs = st_matrix("r(b)")
        Vs = st_matrix("r(V)")
        if (cols(bs) == 0) return(f)
        kk = cols(bs)
        ape_out = bs[kk]
        se_out = sqrt(Vs[kk, kk])
        f.ok = (ape_out != . & se_out != . & se_out > 0)
        f.has_m = 0
        return(f)
    }
    if (use_margins & has_m) {
        // AIE = [dE(y)/dd at m = 1] - [dE(y)/dd at m = 0] with the control
        // function held at its estimated value (cfprobit's predict(pr) is the
        // average structural function). dydx() of the focal at the two
        // moderator values gives both effects and their joint covariance.
        // (margins, at(d=(0 1) m=(0 1)) is NOT equivalent: setting the
        // endogenous variable through at() lets the control function move.)
        rc = _stata("quietly margins, dydx(pape_x1) at(pape_w1=(0 1))", 1)
        if (rc) return(f)
        bs = st_matrix("r(b)")
        Vs = st_matrix("r(V)")
        if (cols(bs) == 0) return(f)
        stripe = st_matrixcolstripe("r(b)")
        eqn = (d.first_stage == "probit" ? "1.pape_x1" : "pape_x1")
        i1 = 0; i2 = 0
        for (j = 1; j <= rows(stripe); j++) {
            if (stripe[j, 1] == eqn) {
                at = pape_leading_int(stripe[j, 2])
                if (at == 1) i1 = j
                else if (at == 2) i2 = j
            }
        }
        if (i1 == 0 | i2 == 0) return(f)
        cvec = J(cols(bs), 1, 0)
        cvec[i1] = -1; cvec[i2] = 1
        ape_out = cvec' * bs'
        se_out = sqrt(cvec' * Vs * cvec)
        f.ok = (ape_out != . & se_out != . & se_out > 0)
        f.has_m = 1
        return(f)
    }
    // map onto W's order
    if (!pape_cf_first(xx.d, xx.X1, d.first_stage, pi, q, c, dcdq)) return(f)
    if (has_m) W = (J(n, 1, 1), xx.d, xx.m, xx.d :* xx.m, xx.Zx, c, c :* xx.m)
    else W = (J(n, 1, 1), xx.d, xx.Zx, c)
    kw = cols(W)
    bs = st_matrix("e(b)")
    Vs = st_matrix("e(V)")
    stripe = st_matrixcolstripe("e(b)")
    eqm = stripe[1, 1]
    eqe = "e." + eqm
    dn = (d.first_stage == "probit" ? "1.pape_x1" : "pape_x1")
    cfn = (d.first_stage == "probit" ? "cf(1.pape_x1)" : "cf(pape_x1)")
    idxv = J(kw, 1, 0)
    idxv[1] = pape_stripe_index(stripe, eqm, "_cons")
    idxv[2] = pape_stripe_index(stripe, eqm, dn)
    if (has_m) {
        idxv[3] = pape_stripe_index(stripe, eqm, "pape_w1")
        idxv[4] = pape_stripe_index(stripe, eqm, dmn)
        for (j = 1; j <= k; j++) idxv[4 + j] = pape_stripe_index(stripe, eqm, sprintf("pape_x%f", j + 1))
        idxv[kw - 1] = pape_stripe_index(stripe, eqe, cfn)
        idxv[kw] = pape_stripe_index(stripe, eqe, cfn + "#c.pape_w1")
    }
    else {
        for (j = 1; j <= k; j++) idxv[2 + j] = pape_stripe_index(stripe, eqm, sprintf("pape_x%f", j + 1))
        idxv[kw] = pape_stripe_index(stripe, eqe, cfn)
    }
    if (any(idxv :== 0)) {
        printf("{err}powerape: could not map cfprobit's coefficients (names: %s)\n", invtokens(stripe[., 2]'))
        return(f)
    }
    f.theta = bs[idxv]'
    f.V_main = Vs[idxv, idxv]
    f.V = f.V_main
    f.W = W
    f.cf = c
    f.pi = pi
    f.has_m = has_m
    f.cf_col = kw - (has_m ? 1 : 0)
    f.cfm_col = (has_m ? kw : .)
    // under separation(keep) (asis) a separated coefficient comes back with a
    // zero variance row; it is kept (see pape_fit_stata). Omitted (o.) terms
    // never reach this point: their names are not found in the stripe above.
    f.ok = (!hasmissing(f.theta) & (asis ? 1 : all(diagonal(f.V_main) :> 0)))
    return(f)
}

string scalar pape_wnames(real scalar k)
{
    string scalar s
    real scalar j
    s = ""
    for (j = 1; j <= k; j++) s = s + sprintf(" pape_w%f", j)
    return(s)
}

end
