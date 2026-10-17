*! powerape_search.mata -- power curves, required n, minimum detectable
*! effect, justified alpha, and the citable power statement.
*! Port of R/curve.R, R/mde.R, R/alpha.R, R/statement.R of powerape 1.11.3.
version 16
mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// Power curve (curve.R: ape_curve)
// ---------------------------------------------------------------------------

struct pape_res scalar pape_curve(struct pape_dgp scalar d, real colvector n,
    string scalar claim, real scalar sesoi, real scalar alpha, real scalar conf,
    string scalar rule, real scalar nsim, real scalar seed, string scalar se,
    string scalar separation)
{
    struct pape_res scalar r, p
    real rowvector lv1
    real colvector sparse_at
    real scalar i, si, first
    string scalar where
    n = sort(uniqrows(n), 1)
    if (rows(n) < 2 | any(n :< 20)) pape_err("n() needs at least two sample sizes, all >= 20.")
    lv1 = pape_claim_levels(claim, alpha, conf, n[1], rule)
    pape_check_coherence(d, claim, sesoi, lv1[4])
    r.results = J(rows(n), 7, .)
    sparse_at = J(0, 1, .)
    first = 0
    for (i = 1; i <= rows(n); i++) {
        si = (seed == . ? . : seed + i - 1)
        p = pape_power_once(d, n[i], claim, sesoi, conf, nsim, si, 0, se, alpha, rule, 0, separation)
        r.results[i, .] = (n[i], p.alpha, p.power, p.mcse, p.outcomes[cols(p.outcomes)], p.separated, p.min_cell)
        if ((p.separated != . & p.separated > 0.01) | (p.min_cell != . & p.min_cell < 10)) {
            sparse_at = sparse_at \ n[i]
            if (first == 0) first = i
        }
    }
    if (rows(sparse_at) > 0) {
        where = "At n = " + invtokens(strofreal(sparse_at', "%12.0g"), ", ")
        (void) pape_sparse_warning(r.results[first, 6], r.results[first, 7], separation, where, d.estimand == "aie")
    }
    r.kind = "curve"
    r.claim = claim; r.sesoi = sesoi
    r.has_rule = (rule != "")
    r.alpha_rule = rule
    r.conf = (r.has_rule ? . : lv1[4])
    r.alpha = (r.has_rule ? . : lv1[1])
    r.nsim = nsim; r.target = d.target_est
    r.estimand = (d.estimand == "" ? "ape" : d.estimand)
    r.model = d.model; r.separation = separation; r.se = se
    r.dgp = d
    pape_stamp_fit(r)
    r.has_draws = 0; r.has_sesoi = (sesoi != .)
    r.n = .; r.power = .; r.mcse = .
    r.outcomes = J(1, 0, .); r.outnames = J(1, 0, "")
    r.separated = .; r.min_cell = .; r.n_failed = .; r.n_separated = .
    r.goal = .; r.confirmed = 0; r.confirm_requested = 0; r.nsim_confirm = .
    r.mde = .; r.is_margin = 0; r.direction = ""
    r.history = J(0, 0, .)
    return(r)
}

void pape_print_curve(struct pape_res scalar x)
{
    real scalar i
    printf("{txt}powerape power curve -- %s claim, %s (%s), nsim = %f per point%s\n",
        x.claim, x.model, strupper(x.estimand), x.nsim,
        (x.has_rule ? "; alpha by sample-size rule (column alpha)" : ""))
    printf("  fit: %s\n", x.fit)
    printf("\n{txt}%10s %9s %8s %8s %8s %10s\n", "n", "alpha", "power", "mcse", "failed", "separated")
    for (i = 1; i <= rows(x.results); i++)
        printf("{res}%10.0f %9.4f %8.4f %8.4f %8.4f %10.4f\n",
            x.results[i, 1], x.results[i, 2], x.results[i, 3], x.results[i, 4], x.results[i, 5], x.results[i, 6])
}

// ---------------------------------------------------------------------------
// Required sample size (curve.R: ape_n)
// ---------------------------------------------------------------------------

real scalar pape_z_at(string scalar claim, real scalar alpha, real scalar conf,
    string scalar rule, real scalar n_eval)
{
    return(pape_zcrit(pape_claim_levels(claim, alpha, conf, n_eval, rule)[4]))
}

real scalar pape_newton_n(real scalar n_from, real scalar p_hat, real scalar power,
    string scalar claim, real scalar alpha, real scalar conf, string scalar rule,
    real scalar r_far)
{
    real scalar p_cl, denom, n_new, n_fp, k
    p_cl = min((max((p_hat, 0.02)), 0.998))
    denom = max((pape_req_delta(p_cl, pape_z_at(claim, alpha, conf, rule, n_from), r_far), 0.1))
    n_new = ceil(n_from * (pape_req_delta(power, pape_z_at(claim, alpha, conf, rule, n_from), r_far) / denom)^2)
    if (rule != "") {
        for (k = 1; k <= 8; k++) {
            n_fp = ceil(n_from * (pape_req_delta(power, pape_z_at(claim, alpha, conf, rule, max((n_new, 20))), r_far) / denom)^2)
            if (n_fp == n_new) break
            n_new = n_fp
        }
    }
    return(n_new)
}

struct pape_res scalar pape_n(struct pape_dgp scalar d, real scalar power,
    string scalar claim, real scalar sesoi, real scalar alpha, real scalar conf,
    string scalar rule, real scalar nsim, real scalar seed, real scalar n_lo,
    real scalar n_hi, real scalar max_iter, real scalar confirm,
    real scalar nsim_confirm, string scalar se, string scalar separation)
{
    struct pape_res scalar r, pw, pwc
    struct pape_sim scalar ps
    real rowvector lv_p, lv_b
    real scalar n_pilot, t_abs, dist, r_far, c_hat, n_i, n_fp, k, it, si, cr, sc, n_c
    real scalar best_n, best_power, best_mcse, best_sep, best_mc, confirmed, found
    real matrix hist
    string scalar se_use
    se_use = se
    if (se != "model" & (d.route == "panel" | d.route == "iv")) {
        printf("{err}warning: se() is fixed by the route; se(%s) is ignored.\n", se)
        se_use = "model"
    }
    n_pilot = 2000
    lv_p = pape_claim_levels(claim, alpha, conf, n_pilot, rule)
    if (power == . | power <= 0.5 | power >= 0.999) pape_err("power() must lie in (0.5, 0.999).")
    pape_check_coherence(d, claim, sesoi, lv_p[4])
    t_abs = abs(d.target_est)
    if (claim == "detect") dist = t_abs
    else if (claim == "minimum") dist = t_abs - sesoi
    else dist = sesoi - t_abs
    if (dist <= 0) pape_err("No solvable sample size: the assumed truth sits on the claim boundary.")
    r_far = (claim == "equivalence" ? (sesoi + t_abs) / dist : .)
    ps = pape_sim_ci(d, n_pilot, 300, lv_p[4], (seed == . ? . : seed - 1), se_use, separation)
    if (mean(ps.ok) < 0.5)
        pape_err("Most pilot replications failed to fit or had separation; the design looks degenerate (e.g. a very rare outcome or a tiny focal group).")
    c_hat = mean(select(ps.se, ps.ok)) * sqrt(n_pilot)
    n_i = ceil((c_hat * pape_req_delta(power, pape_z_at(claim, alpha, conf, rule, n_pilot), r_far) / dist)^2)
    if (rule != "") {
        for (k = 1; k <= 10; k++) {
            n_fp = ceil((c_hat * pape_req_delta(power, pape_z_at(claim, alpha, conf, rule, max((n_i, 20))), r_far) / dist)^2)
            if (n_fp == n_i) break
            n_i = n_fp
        }
    }
    hist = J(0, 6, .)     // stage (1 search, 2 confirm), iter, n, alpha, power, mcse
    best_n = .; best_power = .; best_mcse = .; best_sep = .; best_mc = .
    for (it = 1; it <= max_iter; it++) {
        n_i = max((n_lo, min((n_hi, n_i))))
        si = (seed == . ? . : seed + it)
        pw = pape_power_once(d, n_i, claim, sesoi, conf, nsim, si, 0, se_use, alpha, rule, 0, separation)
        hist = hist \ (1, it, n_i, pw.alpha, pw.power, pw.mcse)
        if (best_n == . | abs(pw.power - power) < abs(best_power - power)) {
            best_n = n_i; best_power = pw.power; best_mcse = pw.mcse
            best_sep = pw.separated; best_mc = pw.min_cell
        }
        if (abs(pw.power - power) <= max((0.01, 2 * pw.mcse))) break
        n_i = pape_newton_n(n_i, pw.power, power, claim, alpha, conf, rule, r_far)
        found = 0
        for (k = 1; k <= rows(hist); k++) if (hist[k, 3] == n_i) found = 1
        if (found) break
    }
    confirmed = 0
    if (confirm) {
        if (nsim_confirm == . | nsim_confirm < nsim) pape_err("nsimconfirm() must be at least nsim().")
        n_c = best_n
        for (cr = 1; cr <= 3; cr++) {
            sc = (seed == . ? . : seed + 100 + cr)
            pwc = pape_power_once(d, n_c, claim, sesoi, conf, nsim_confirm, sc, 0, se_use, alpha, rule, 0, separation)
            hist = hist \ (2, cr, n_c, pwc.alpha, pwc.power, pwc.mcse)
            best_n = n_c; best_power = pwc.power; best_mcse = pwc.mcse
            best_sep = pwc.separated; best_mc = pwc.min_cell
            if (pwc.power >= power - max((0.005, 1.5 * pwc.mcse))) {
                confirmed = 1
                break
            }
            if (cr == 3) break
            n_c = min((n_hi, max((pape_newton_n(n_c, pwc.power, power, claim, alpha, conf, rule, r_far), ceil(n_c * 1.02), n_c + 1))))
        }
        if (!confirmed)
            printf("{err}warning: the confirmation stage did not reach the target within tolerance; the returned n = %f has high-precision power %5.3f. Increase nsimconfirm() or treat the answer as a lower bound.\n", best_n, best_power)
    }
    if (best_n >= n_hi)
        printf("{err}warning: required n reached the nrange() ceiling; treat the result as a lower bound.\n")
    (void) pape_sparse_warning(best_sep, best_mc, separation, sprintf("At the returned n = %f", best_n), d.estimand == "aie")
    lv_b = pape_claim_levels(claim, alpha, conf, best_n, rule)
    r.kind = "n"
    r.n = best_n; r.power = best_power; r.mcse = best_mcse
    r.goal = power; r.claim = claim; r.sesoi = sesoi
    r.conf = lv_b[4]; r.alpha = lv_b[1]
    r.alpha_rule = rule; r.has_rule = (rule != "")
    r.nsim = nsim
    r.nsim_confirm = (confirm ? nsim_confirm : .)
    r.confirm_requested = confirm
    r.confirmed = confirmed
    r.history = hist
    r.se = se_use; r.separation = separation
    r.separated = best_sep; r.min_cell = best_mc
    r.target = d.target_est
    r.estimand = (d.estimand == "" ? "ape" : d.estimand)
    r.model = d.model
    r.dgp = d
    pape_stamp_fit(r)
    r.has_draws = 0; r.has_sesoi = (sesoi != .)
    r.outcomes = J(1, 0, .); r.outnames = J(1, 0, "")
    r.n_failed = .; r.n_separated = .
    r.mde = .; r.is_margin = 0; r.direction = ""
    r.results = J(0, 0, .)
    return(r)
}

void pape_print_n(struct pape_res scalar x)
{
    string scalar lab, unit_lab
    real scalar n_search, n_conf
    printf("{txt}powerape required sample size -- %s claim\n", x.claim)
    lab = (x.confirmed ? "confirmed" : "achieved")
    if (x.dgp.route == "panel") unit_lab = sprintf("n = %f units (x %s periods)", x.n, pape_panel_len_label(x.dgp))
    else unit_lab = sprintf("n = %f", x.n)
    printf("  %s for %2.0f%% target power (%s %5.3f, MCSE %5.3f)\n", unit_lab, 100 * x.goal, lab, x.power, x.mcse)
    printf("  assumed true %s %s%7.4f, sesoi %s, %6.4g%% CI (alpha = %6.4g%s), %s%s\n",
        strupper(x.estimand), pape_sgn(x.target), x.target, (x.sesoi == . ? "-" : sprintf("%5.3f", x.sesoi)),
        100 * x.conf, x.alpha, (x.has_rule ? " by sample-size rule" : ""), x.model, pape_se_lab(x))
    printf("  fit: %s\n", x.fit)
    n_search = sum(x.history[., 1] :== 1)
    n_conf = sum(x.history[., 1] :== 2)
    if (x.confirm_requested) {
        if (x.confirmed)
            printf("  search: %f step(s); confirmed in %f round(s) at nsim = %f.\n", n_search, n_conf, x.nsim_confirm)
        else
            printf("  search: %f step(s); confirmation (%f round(s), nsim = %f) fell short -- treat as a lower bound.\n", n_search, n_conf, x.nsim_confirm)
    }
    else printf("  search: %f step(s); confirm with powerape power at a larger nsim.\n", n_search)
    pape_print_sep_note(x)
}

// ---------------------------------------------------------------------------
// Minimum detectable effect (mde.R: ape_mde, finish_mde)
// ---------------------------------------------------------------------------

// mean SE of the estimator in world d_at at the design's n; rc 1 on degenerate
real scalar pape_pilot_se(struct pape_dgp scalar d_at, real scalar n, real scalar conf,
    real scalar seed_use, string scalar st, string scalar separation, real scalar se_out)
{
    struct pape_sim scalar ps
    ps = pape_sim_ci(d_at, n, 300, conf, seed_use, st, separation)
    if (mean(ps.ok) < 0.5) return(1)
    se_out = mean(select(ps.se, ps.ok))
    return(0)
}

struct pape_res scalar pape_mde(struct pape_dgp scalar dgp, real scalar n,
    real scalar power, string scalar claim, real scalar sesoi, real scalar alpha,
    real scalar conf, string scalar rule, real scalar nsim, real scalar seed,
    real scalar mf_in, real scalar mm_in, real scalar max_iter, real scalar confirm,
    real scalar nsim_confirm, string scalar se, string scalar direction,
    string scalar separation)
{
    struct pape_res scalar pw, pc
    struct pape_dgp scalar d_ref, d_k, d_t
    real rowvector lv
    real scalar sgn, is_aie, mf, mm, z, zg, t_abs, se_hat, d_i, it, p_cl, se_impl
    real scalar delta0, t_i, t_ref, k, se_k, rc, ceilv, lo, hi, mid, p_ceil, rc2
    real scalar best_value, best_power, best_mcse, best_sep, best_mc
    real matrix hist
    real colvector t_try
    string scalar msg, st, dir, se_use, est_lab
    se_use = se
    if (se != "model" & (dgp.route == "panel" | dgp.route == "iv")) {
        printf("{err}warning: se() is fixed by the route; se(%s) is ignored.\n", se)
        se_use = "model"
    }
    if (n == . | n < 20) pape_err("n() must be at least 20.")
    if (power == . | power <= 0.5 | power >= 0.999) pape_err("power() must lie in (0.5, 0.999).")
    if (nsim == . | nsim < 20) pape_err("nsim() must be at least 20.")
    dir = direction
    if (dir == "") dir = (dgp.has_effect & dgp.target_est < 0 ? "negative" : "positive")
    if (dir != "positive" & dir != "negative") pape_err("direction() must be positive or negative.")
    sgn = (dir == "negative" ? -1 : 1)
    lv = pape_claim_levels(claim, alpha, conf, n, rule)
    conf = lv[4]
    if (claim == "minimum" & (sesoi == . | sesoi <= 0))
        pape_err("claim(minimum) needs a single positive sesoi() (in APE units).")
    is_aie = dgp.has_m
    mf = (mf_in == . ? dgp.main_focal : mf_in)
    mm = (mm_in == . ? dgp.main_moderator : mm_in)
    if (is_aie & (mf == . | mm == .))
        pape_err("An AIE MDE needs the main-effect anchors: pin the DGP with powerape setaie first, or pass mainfocal() and mainmoderator().")
    est_lab = (is_aie ? "AIE" : "APE")
    z = pape_zcrit(conf)
    zg = invnormal(power)
    st = ((dgp.route == "panel" | dgp.route == "iv") ? "model" : se_use)
    hist = J(0, 6, .)   // stage, iter, value, alpha, power, mcse
    best_value = .; best_power = .; best_mcse = .; best_sep = .; best_mc = .

    // --- equivalence: search the margin at the pinned target ----------------
    if (claim == "equivalence") {
        if (!dgp.has_effect)
            pape_err("An equivalence MDE holds the truth fixed: pin the DGP first (typically at target = 0) with powerape setape or setaie.")
        t_abs = abs(dgp.target_est)
        if (pape_pilot_se(dgp, n, conf, (seed == . ? . : seed - 1), st, separation, se_hat))
            pape_err("Most pilot replications failed to fit or had separation; the design looks degenerate.")
        d_i = pape_req_margin(power, se_hat, t_abs, z)
        for (it = 1; it <= max_iter; it++) {
            pw = pape_power_once(dgp, n, "equivalence", d_i, conf, nsim, (seed == . ? . : seed + it), 0, se_use, ., "", 0, separation)
            hist = hist \ (1, it, d_i, pw.alpha, pw.power, pw.mcse)
            if (best_value == . | abs(pw.power - power) < abs(best_power - power)) {
                best_value = d_i; best_power = pw.power; best_mcse = pw.mcse
                best_sep = pw.separated; best_mc = pw.min_cell
            }
            if (abs(pw.power - power) <= max((0.01, 2 * pw.mcse))) break
            p_cl = min((max((pw.power, 0.02)), 0.998))
            se_impl = (d_i - t_abs) / max((pape_req_delta(p_cl, z, (d_i + t_abs) / (d_i - t_abs)), 0.1))
            d_i = pape_req_margin(power, se_impl, t_abs, z)
        }
        return(pape_finish_mde(dgp, n, best_value, best_power, best_mcse, best_sep, best_mc, hist,
            power, "equivalence", best_value, conf, nsim, seed, confirm, nsim_confirm, se_use,
            mf, mm, 1, 0, t_abs, lv[1], rule, separation, 1))
    }

    // --- detect / minimum: search the target ---------------------------------
    delta0 = (claim == "minimum" ? sesoi : 0)
    if (dgp.has_effect) d_ref = dgp
    else {
        t_try = sgn :* (0.5 * min((dgp.baseline, 1 - dgp.baseline)) \ 0.05 \ 0.02)
        rc = 1
        for (k = 1; k <= 3; k++) {
            d_ref = dgp
            rc = pape_repin(d_ref, t_try[k], mf, mm, msg)
            if (!rc) break
        }
        if (rc) pape_err("Could not pin a reference effect for the pilot stage.")
    }
    if (pape_pilot_se(d_ref, n, conf, (seed == . ? . : seed - 1), st, separation, se_hat))
        pape_err("Most pilot replications failed to fit or had separation; the design looks degenerate.")
    t_i = sgn * (delta0 + se_hat * (z + zg))
    t_ref = d_ref.target_est
    for (k = 1; k <= 2; k++) {
        if (abs(t_i - t_ref) <= 0.02 * abs(t_i)) break
        d_k = dgp
        if (pape_repin(d_k, t_i, mf, mm, msg)) break
        if (pape_pilot_se(d_k, n, conf, (seed == . ? . : seed - 1 - k), st, separation, se_k)) break
        t_ref = t_i
        t_i = sgn * (delta0 + se_k * (z + zg))
    }
    for (it = 1; it <= max_iter; it++) {
        d_t = dgp
        rc = pape_repin(d_t, t_i, mf, mm, msg)
        if (rc) {
            // candidate beyond the feasible range: report the ceiling honestly
            lo = delta0 + 1e-4; hi = abs(t_i)
            for (k = 1; k <= 30; k++) {
                mid = (lo + hi) / 2
                d_t = dgp
                if (pape_repin(d_t, sgn * mid, mf, mm, msg)) hi = mid
                else lo = mid
            }
            ceilv = lo
            d_t = dgp
            rc2 = pape_repin(d_t, sgn * 0.98 * ceilv, mf, mm, msg)
            p_ceil = .
            if (!rc2) {
                pc = pape_power_once(d_t, n, claim, sesoi, conf, nsim, (seed == . ? . : seed + 50), 0, se_use, ., "", 0, separation)
                p_ceil = pc.power
            }
            pape_err(sprintf("No detectable %s at %2.0f%% power within the feasible range: the largest attainable effect is about %6.3f, where simulated power reaches only %4.2f. Increase n (or loosen the goal).",
                est_lab, 100 * power, sgn * ceilv, p_ceil))
        }
        pw = pape_power_once(d_t, n, claim, sesoi, conf, nsim, (seed == . ? . : seed + it), 0, se_use, ., "", 0, separation)
        hist = hist \ (1, it, t_i, pw.alpha, pw.power, pw.mcse)
        if (best_value == . | abs(pw.power - power) < abs(best_power - power)) {
            best_value = t_i; best_power = pw.power; best_mcse = pw.mcse
            best_sep = pw.separated; best_mc = pw.min_cell
        }
        if (abs(pw.power - power) <= max((0.01, 2 * pw.mcse))) break
        p_cl = min((max((pw.power, 0.02)), 0.998))
        se_impl = (abs(t_i) - delta0) / max((z + invnormal(p_cl), 0.1))
        t_i = sgn * (delta0 + se_impl * (z + zg))
    }
    return(pape_finish_mde(dgp, n, best_value, best_power, best_mcse, best_sep, best_mc, hist,
        power, claim, sesoi, conf, nsim, seed, confirm, nsim_confirm, se_use,
        mf, mm, 0, delta0, 0, lv[1], rule, separation, sgn))
}

struct pape_res scalar pape_finish_mde(struct pape_dgp scalar dgp, real scalar n,
    real scalar best_value, real scalar best_power, real scalar best_mcse,
    real scalar best_sep, real scalar best_mc, real matrix hist, real scalar power,
    string scalar claim, real scalar sesoi, real scalar conf, real scalar nsim,
    real scalar seed, real scalar confirm, real scalar nsim_confirm, string scalar se,
    real scalar mf, real scalar mm, real scalar is_margin, real scalar delta0,
    real scalar t_abs, real scalar alpha_val, string scalar rule, string scalar separation,
    real scalar sgn)
{
    struct pape_res scalar r, pwc
    struct pape_dgp scalar d_t, dgp_at
    real scalar z, zg, anchor, confirmed, v_c, cr, sc, p_cl, se_impl, rc
    string scalar msg
    z = pape_zcrit(conf)
    zg = invnormal(power)
    anchor = (is_margin ? t_abs : delta0)
    confirmed = 0
    if (confirm) {
        if (nsim_confirm == . | nsim_confirm < nsim) pape_err("nsimconfirm() must be at least nsim().")
        v_c = best_value
        for (cr = 1; cr <= 3; cr++) {
            sc = (seed == . ? . : seed + 100 + cr)
            if (is_margin)
                pwc = pape_power_once(dgp, n, "equivalence", v_c, conf, nsim_confirm, sc, 0, se, ., "", 0, separation)
            else {
                d_t = dgp
                rc = pape_repin(d_t, v_c, mf, mm, msg)
                if (rc) break
                pwc = pape_power_once(d_t, n, claim, sesoi, conf, nsim_confirm, sc, 0, se, ., "", 0, separation)
            }
            hist = hist \ (2, cr, v_c, pwc.alpha, pwc.power, pwc.mcse)
            best_value = v_c; best_power = pwc.power; best_mcse = pwc.mcse
            best_sep = pwc.separated; best_mc = pwc.min_cell
            if (pwc.power >= power - max((0.005, 1.5 * pwc.mcse))) {
                confirmed = 1
                break
            }
            if (cr == 3) break
            p_cl = min((max((pwc.power, 0.02)), 0.998))
            if (is_margin) {
                se_impl = (v_c - anchor) / max((pape_req_delta(p_cl, z, (v_c + anchor) / (v_c - anchor)), 0.1))
                v_c = pape_req_margin(power, se_impl, anchor, z)
            }
            else {
                se_impl = (abs(v_c) - anchor) / max((z + invnormal(p_cl), 0.1))
                v_c = sgn * (anchor + se_impl * (z + zg))
            }
        }
    }
    (void) pape_sparse_warning(best_sep, best_mc, separation, "At the returned " + (is_margin ? "margin" : "MDE"), dgp.has_m)
    dgp_at = dgp
    if (!is_margin) {
        rc = pape_repin(dgp_at, best_value, mf, mm, msg)
        if (rc) pape_err(msg)
    }
    r.kind = "mde"
    r.mde = best_value; r.power = best_power; r.mcse = best_mcse
    r.goal = power; r.claim = claim
    r.sesoi = (is_margin ? best_value : sesoi)
    r.conf = conf; r.alpha = alpha_val
    r.alpha_rule = rule; r.has_rule = (rule != "")
    r.n = n; r.nsim = nsim
    r.nsim_confirm = nsim_confirm
    r.confirmed = confirmed; r.confirm_requested = confirm
    r.is_margin = is_margin; r.se = se
    r.direction = (is_margin ? "" : (sgn < 0 ? "negative" : "positive"))
    r.separation = separation; r.separated = best_sep; r.min_cell = best_mc
    r.estimand = (dgp.has_m ? "aie" : "ape")
    r.model = dgp.model
    r.history = hist
    r.dgp = dgp_at
    r.target = dgp_at.target_est
    pape_stamp_fit(r)
    r.has_draws = 0; r.has_sesoi = (r.sesoi != .)
    r.outcomes = J(1, 0, .); r.outnames = J(1, 0, "")
    r.n_failed = .; r.n_separated = .
    r.results = J(0, 0, .)
    return(r)
}

void pape_print_mde(struct pape_res scalar x)
{
    string scalar unit_lab
    real scalar n_search, n_conf
    if (x.dgp.route == "panel") unit_lab = sprintf("n = %f units (x %s periods)", x.n, pape_panel_len_label(x.dgp))
    else unit_lab = sprintf("n = %f", x.n)
    if (x.is_margin) {
        printf("{txt}powerape minimum establishable equivalence margin -- %s\n", strupper(x.estimand))
        printf("  bounds +/-%9.4f at %s for %2.0f%% target power (%s %5.3f, MCSE %5.3f)\n",
            x.mde, unit_lab, 100 * x.goal, (x.confirmed ? "confirmed" : "achieved"), x.power, x.mcse)
    }
    else {
        printf("{txt}powerape minimum detectable %s -- %s claim%s\n", strupper(x.estimand), x.claim,
            (x.direction == "negative" ? ", decreases" : ""))
        printf("  MDE = %9.4f at %s for %2.0f%% target power (%s %5.3f, MCSE %5.3f)\n",
            x.mde, unit_lab, 100 * x.goal, (x.confirmed ? "confirmed" : "achieved"), x.power, x.mcse)
        if (x.claim == "minimum") {
            if (x.direction == "negative") printf("  smallest decrease demonstrably beyond -sesoi = %6.3f\n", -x.sesoi)
            else printf("  smallest effect demonstrably above sesoi %5.3f\n", x.sesoi)
        }
    }
    if (x.has_rule) printf("  alpha = %6.4g by sample-size rule (%6.4g%% CI)\n", x.alpha, 100 * x.conf)
    printf("  fit: %s\n", x.fit)
    n_search = sum(x.history[., 1] :== 1)
    n_conf = sum(x.history[., 1] :== 2)
    if (x.confirm_requested) {
        if (x.confirmed)
            printf("  search: %f step(s); confirmed in %f round(s) at nsim = %f.\n", n_search, n_conf, x.nsim_confirm)
        else
            printf("  search: %f step(s); confirmation (%f round(s), nsim = %f) fell short -- treat as a lower bound.\n", n_search, n_conf, x.nsim_confirm)
    }
    else printf("  search: %f step(s); confirm with powerape power at a larger nsim.\n", n_search)
    pape_print_sep_note(x)
}

// ---------------------------------------------------------------------------
// Justified alpha (alpha.R)
// ---------------------------------------------------------------------------

// rows: alpha, conf, power, mcse, outcomes...
real matrix pape_power_at_draws(struct pape_sim scalar draws, real scalar target,
    string scalar claim, real scalar sesoi, real scalar nsim, real colvector alpha)
{
    struct pape_res scalar r
    real matrix out
    real scalar i
    out = J(rows(alpha), 4 + (sesoi != . ? 5 : 3), .)
    for (i = 1; i <= rows(alpha); i++) {
        pape_summarize_sim(draws, target, claim, sesoi, nsim, alpha[i], r)
        out[i, .] = (alpha[i], pape_claim_levels(claim, alpha[i], ., ., "")[4], r.power, r.mcse, r.outcomes)
    }
    return(out)
}

struct pape_dgp scalar pape_boundary_dgp(struct pape_dgp scalar dgp, string scalar claim,
    real scalar sesoi, real scalar target)
{
    struct pape_dgp scalar d
    real scalar sgn, b, rc
    string scalar msg
    sgn = (target < 0 ? -1 : 1)
    b = (claim == "detect" ? 0 : sgn * sesoi)
    d = dgp
    if (d.estimand == "aie") rc = pape_set_aie(d, b, d.main_focal, d.main_moderator, msg)
    else rc = pape_set_ape(d, b, msg)
    if (rc) pape_err(msg)
    return(d)
}

struct pape_alpha scalar pape_alpha_fn(struct pape_res scalar x, real scalar cost,
    real scalar prior, string scalar error, real scalar cap, real scalar size,
    real scalar nsim_size, real scalar seed)
{
    struct pape_alpha scalar a
    struct pape_res scalar r
    struct pape_sim scalar sim
    struct pape_dgp scalar d_b
    real colvector grid, pw, beta, wcer, obj, in_cap, flat
    real matrix pa
    real scalar i_free, i_cap, tol, i_def, ns, sgn, j
    string scalar se_type
    if (x.kind != "power") pape_err("powerape alpha needs a powerape power result.")
    if (!x.has_draws) pape_err("The result carries no stored draws; rerun powerape power without nodraws.")
    if (cost == . | cost <= 0 | prior == . | prior <= 0) pape_err("cost() and prior() must be positive.")
    if (cap == . | cap <= 0 | cap >= 0.5) pape_err("cap() must lie in (0, 0.5).")
    if (cap > 0.05 + 1e-9)
        printf("{err}warning: cap() exceeds .05: an error rate above the conventional level needs a justification of its own (Maier & Lakens, 2022).\n")
    grid = rangen(1e-4, 0.01, 100) \ rangen(0.0105, 0.4995, 979)
    pa = pape_power_at_draws(x.draws, x.target, x.claim, x.sesoi, x.nsim, grid)
    pw = pa[., 3]
    beta = 1 :- pw
    wcer = (cost :* grid :+ prior :* beta) :/ (cost + prior)
    obj = (error == "minimize" ? wcer : abs(cost :* grid :- prior :* beta))
    i_free = selectindex(obj :== min(obj))[1]
    in_cap = (grid :<= cap + 1e-9)
    i_cap = 0
    for (j = 1; j <= rows(grid); j++) {
        if (in_cap[j]) {
            if (i_cap == 0) i_cap = j
            else if (obj[j] < obj[i_cap]) i_cap = j
        }
    }
    a.capped = (grid[i_free] > cap + 1e-9)
    a.alpha = grid[i_cap]
    a.power = pw[i_cap]
    a.mcse = sqrt(a.power * (1 - a.power) / x.nsim)
    tol = (error == "minimize" ? prior * a.mcse / (cost + prior) : prior * a.mcse)
    flat = in_cap :& (obj :<= obj[i_cap] + tol)
    a.range = (min(select(grid, flat)), max(select(grid, flat)))
    i_def = selectindex(abs(grid :- 0.05) :== min(abs(grid :- 0.05)))[1]
    a.beta = beta[i_cap]; a.wcer = wcer[i_cap]
    a.power_default = pw[i_def]; a.wcer_default = wcer[i_def]
    a.cap = cap
    a.alpha_uncapped = grid[i_free]; a.beta_uncapped = beta[i_free]; a.wcer_uncapped = wcer[i_free]
    a.cost = cost; a.prior = prior; a.error = error
    a.claim = x.claim; a.sesoi = x.sesoi; a.n = x.n; a.nsim = x.nsim
    a.target = x.target; a.estimand = x.estimand; a.model = x.model; a.se = x.se
    a.grid = (grid, pw, beta, wcer)
    a.dgp = x.dgp
    a.has_size = 0
    if (size) {
        ns = (nsim_size == . ? x.nsim : nsim_size)
        if (ns < 20) pape_err("nsimsize() must be at least 20.")
        d_b = pape_boundary_dgp(x.dgp, x.claim, x.sesoi, x.target)
        sgn = (x.target < 0 ? -1 : 1)
        se_type = ((x.dgp.route == "panel" | x.dgp.route == "iv") ? "model" : x.se)
        sim = pape_sim_ci(d_b, x.n, ns, pape_claim_levels(x.claim, a.alpha, ., ., "")[4], seed, se_type, x.separation)
        pape_summarize_sim(sim, sgn, x.claim, x.sesoi, ns, a.alpha, r)
        a.has_size = 1
        a.size = r.power
        a.size_mcse = sqrt(r.power * (1 - r.power) / ns)
        a.size_nsim = ns
        a.size_truth = d_b.target_est
        a.size_failed = r.n_failed
        a.size_separated = r.separated
    }
    return(a)
}

void pape_print_alpha(struct pape_alpha scalar x)
{
    printf("{txt}powerape justified alpha -- %s, n = %f, nsim = %f\n", pape_claim_lab(x.claim, x.sesoi), x.n, x.nsim)
    printf("  objective: %s (Type I cost x%g, prior odds H1:H0 = %g)\n",
        (x.error == "minimize" ? "minimize the weighted combined error rate" : "balance the cost-weighted Type I and Type II error rates"),
        x.cost, x.prior)
    printf("  alpha* = %6.4f -> power %5.3f (MCSE %5.3f), beta %5.3f, weighted combined error %6.4f\n",
        x.alpha, x.power, x.mcse, x.beta, x.wcer)
    printf("  at alpha = 0.05: power %5.3f, weighted combined error %6.4f\n", x.power_default, x.wcer_default)
    printf("  flat range: alpha in [%6.4f, %6.4f] is indistinguishable from the optimum at this nsim\n", x.range[1], x.range[2])
    if (x.capped) {
        printf("  cap %5.3f BINDS: unconstrained optimum alpha = %6.4f (weighted combined error %6.4f);\n", x.cap, x.alpha_uncapped, x.wcer_uncapped)
        printf("    raising alpha above the conventional level needs its own justification (Maier & Lakens, 2022).\n")
    }
    else printf("  cap %5.3f not binding (unconstrained optimum %6.4f)\n", x.cap, x.alpha_uncapped)
    if (x.has_size)
        printf("  realized size at the claim boundary (true %s = %5.3f): %6.4f (MCSE %6.4f, nsim %f)\n",
            strupper(x.estimand), x.size_truth, x.size, x.size_mcse, x.size_nsim)
}

real scalar pape_alpha_frontier(real scalar cost, real scalar power, real scalar prior, string scalar claim)
{
    real scalar m, rhs, a
    if (cost == . | cost <= 0 | prior == . | prior <= 0) pape_err("cost() and prior() must be positive.")
    if (power == . | power <= 0 | power >= 1) pape_err("power() must lie in (0, 1).")
    m = (claim == "detect" ? 2 : 1)
    rhs = prior * normalden(invnormal(power)) / (m * cost)
    if (rhs >= normalden(0))
        pape_err(sprintf("No interior optimum at the frontier: with cost = %g and prior = %g, the weighted combined error rate keeps falling as alpha rises toward 1 for the just-detectable effect, so no level below 1 is optimal. Weight Type I errors more heavily (larger cost / smaller prior), or run powerape alpha on a stored simulation, which reports the capped optimum for a specific design.", cost, prior))
    a = m * normal(-sqrt(-2 * ln(sqrt(2 * pi()) * rhs)))
    if (a >= 0.5)
        printf("{err}warning: the frontier level (%5.3f) is not below 0.5 and cannot be used as an alpha in powerape's power functions; this cost/prior weighting does not support a conventional test.\n", a)
    return(a)
}

// ---------------------------------------------------------------------------
// Citable power statement (statement.R)
// ---------------------------------------------------------------------------

string scalar pape_var_desc(struct pape_dgp scalar d, struct pape_var scalar v)
{
    string scalar base
    real scalar icc, i, j
    struct pape_var colvector pvars
    if (v.type == "binary") base = sprintf("binary, prevalence %4.2f", v.p)
    else base = sprintf("continuous, mean %6.3g, SD %6.3g", v.mean, v.sd)
    if (d.route == "panel") {
        pvars = pape_panel_vars(d)
        i = 0
        for (j = 1; j <= rows(pvars); j++) if (pvars[j].name == v.name) i = j
        if (i > 0) {
            icc = d.iccs[i]
            if (icc >= 1) base = base + ", time-constant"
            else if (icc <= 0) base = base + ", redrawn each period"
            else base = base + sprintf(", within-unit persistence (icc) %4.2f", icc)
        }
    }
    return(sprintf("%s (%s)", v.name, base))
}

string scalar pape_dep_desc(struct pape_dgp scalar d, string scalar over)
{
    real matrix R
    real colvector off
    real scalar i, j
    string scalar lvl
    if (d.corr_given == 0) return("")
    if (d.corr_given == 2) {
        R = d.corr_matrix
        off = J(0, 1, .)
        for (i = 1; i <= rows(R); i++) for (j = i + 1; j <= cols(R); j++) off = off \ R[i, j]
    }
    else off = d.corr_scalar
    if (rows(off) == 0 | all(abs(off) :< 1e-12)) return("")
    lvl = (d.route == "panel" ? ", at both the unit and the within-unit level" : "")
    if (d.corr_given == 2)
        return(sprintf("a latent (Gaussian-copula) correlation matrix over %s with off-diagonal entries from %4.2f to %4.2f%s", over, min(off), max(off), lvl))
    return(sprintf("an exchangeable latent (Gaussian-copula) correlation of %4.2f among %s%s", d.corr_scalar, over, lvl))
}

string scalar pape_describe_dgp(struct pape_dgp scalar d)
{
    string scalar parts, covs, foc, cv, dep, sig, over, s, insts
    string rowvector blk
    real scalar is_panel, is_iv, n_vars, j
    is_panel = (d.route == "panel"); is_iv = (d.route == "iv")
    if (is_panel)
        parts = sprintf("a correlated random effects %s model over %s periods per unit (latent unit-effect share rho = %4.2f, %2.0f%% of it loaded on regressor means), estimated by pooled %s with Mundlak means%s and unit-clustered standard errors, with focal variable %s",
            d.model, pape_panel_len_label(d), d.rho, 100 * d.cre_share, d.model,
            (d.unbalanced ? " over observed periods plus period-count cohort indicators (Wooldridge, 2019)" : ""),
            pape_var_desc(d, d.focal))
    else if (is_iv) {
        insts = ""
        for (j = 1; j <= d.n_inst; j++) insts = insts + (j > 1 ? ", " : "") + pape_var_desc(d, d.instruments[j])
        parts = sprintf("a probit model with an endogenous focal variable %s (error correlation %4.2f with the first stage), instrumented by %s (%s first stage, instrument strength %4.2f as a share of the focal's %svariance) and estimated by control-function probit with stacked method-of-moments robust standard errors (the estimator of Stata's cfprobit)",
            pape_var_desc(d, d.focal), d.endogeneity, insts, d.first_stage, d.iv_strength,
            (d.focal.type == "binary" ? "latent " : ""))
    }
    else parts = sprintf("a %s model with focal variable %s", d.model, pape_var_desc(d, d.focal))
    if (d.has_m) parts = parts + sprintf(", moderator %s, and their interaction", pape_var_desc(d, d.moderator))
    if (d.builder == "from_fit")
        covs = sprintf("covariate rows and nuisance coefficients taken from a fitted pilot model (%s; %f pilot observations, resampled with replacement, the focal jointly with them)", d.formula_str, d.n_emp)
    else if (d.route == "empirical") {
        foc = (d.emp_focal_joint ? "the focal resampled jointly with the pilot rows" : "the focal drawn independently of the resampled rows")
        if (d.k > 0) covs = sprintf("covariates resampled with replacement from %f pilot observations (%s), %s", d.n_emp, invtokens(d.covnames', ", "), foc)
        else covs = sprintf("no covariates, %s (%f pilot observations)", foc, d.n_emp)
    }
    else {
        blk = J(1, 0, "")
        if (!is_iv) blk = blk, "the focal variable"
        if (d.has_m) blk = blk, "the moderator"
        if (d.k > 0) blk = blk, "the covariates"
        if (is_iv) blk = blk, "the instruments"
        n_vars = (!is_iv) + d.has_m + d.k + (is_iv ? d.n_inst : 0)
        if (cols(blk) > 1) over = invtokens(blk[1..(cols(blk) - 1)], ", ") + " and " + blk[cols(blk)]
        else over = blk[1]
        dep = (n_vars > 1 ? pape_dep_desc(d, over) : "")
        if (d.k > 0) {
            s = ""
            for (j = 1; j <= d.k; j++) s = s + (j > 1 ? ", " : "") + pape_var_desc(d, d.covariates[j])
            cv = "parametric covariates " + s
        }
        else cv = "no additional covariates"
        covs = cv + (dep != "" ? ", with " + dep : (n_vars > 1 ? ", with " + over + " mutually independent" : ""))
    }
    if (d.link == "logit" & d.s2 > 0)
        sig = sprintf("nuisance covariates contribute %4.2f of the latent-index variance at unit error variance (a McKelvey-Zavoina R-squared of %5.3f on the logistic scale)", d.signal_implied, pape_mz_r2_logit(d.s2))
    else sig = sprintf("nuisance covariates contribute a latent pseudo-R-squared of %4.2f", d.signal_implied)
    return(sprintf("%s; %s; baseline outcome rate %5.3f with the focal%s at reference; %s",
        parts, covs, d.baseline, (d.has_m ? " and moderator" : ""), sig))
}

string scalar pape_claim_text(string scalar claim, real scalar sesoi, real scalar conf, real scalar rule)
{
    string scalar cl, al_txt, by_n
    real scalar al
    cl = sprintf("%g%%", round(100 * conf, 0.01))
    al = (claim == "detect" ? 1 - conf : (1 - conf) / 2)
    al_txt = sprintf("%g%%", round(100 * al, 0.001))
    by_n = (rule ? ", set as a function of the sample size" : "")
    if (claim == "minimum")
        return(sprintf("the minimum-effect claim (the %s confidence interval's lower bound exceeding the smallest effect size of interest, %5.3f; a one-sided test at alpha = %s%s)", cl, sesoi, al_txt, by_n))
    if (claim == "detect")
        return(sprintf("the detection claim (the %s confidence interval excluding zero in the hypothesized direction; a two-sided test at alpha = %s%s)", cl, al_txt, by_n))
    return(sprintf("the equivalence claim (the %s confidence interval lying within +/-%5.3f; two one-sided tests at alpha = %s%s)", cl, sesoi, al_txt, by_n))
}

string scalar pape_version() return("0.1.3 (Stata port of powerape 1.11.3)")

string scalar pape_statement(struct pape_res scalar x)
{
    struct pape_dgp scalar d
    string scalar est_long, est_short, subject, s1, s2, s3, s4, s5, s6, n_text, conf_frag, est_desc, cell, size_noun
    real scalar is_aie, is_panel, rule, keep, sepd, other
    d = x.dgp
    rule = x.has_rule
    is_aie = (x.estimand == "aie")
    est_long = (is_aie ? "average interaction effect (AIE)" : "average partial effect (APE)")
    est_short = (is_aie ? "AIE" : "APE")
    subject = (is_aie ? sprintf("%s-by-%s", d.focal.name, d.moderator.name) : d.focal.name)
    s1 = sprintf("We conducted a simulation-based power analysis for the %s of %s using the powerape package for Stata (version %s), following the confidence-interval approach of Riesthuis (2024).",
        est_long, subject, pape_version())
    s2 = sprintf("The assumed data-generating process was %s.", pape_describe_dgp(d))
    if (x.kind == "mde") s3 = ""
    else s3 = sprintf("The assumed true %s was %5.3f (%3.1f percentage points)%s.", est_short, x.target, 100 * x.target,
        (is_aie ? sprintf(", with conditional-at-reference main-effect APEs of %5.3f (%s) and %5.3f (%s)", d.main_focal, d.focal.name, d.main_moderator, d.moderator.name) : ""))
    is_panel = (d.route == "panel")
    if (is_panel) n_text = sprintf("n = %f units observed over %s periods (%s unit-period observations)", x.n, pape_panel_len_label(d), pape_panel_obs_label(d, x.n))
    else n_text = sprintf("n = %f", x.n)
    if (x.kind == "mde") {
        if (x.confirmed) conf_frag = sprintf(", confirmed by a high-precision verification run (simulated power %5.3f, Monte Carlo SE %5.3f, %f replications)", x.power, x.mcse, x.nsim_confirm)
        else conf_frag = sprintf(" (simulated power %5.3f, Monte Carlo SE %5.3f)", x.power, x.mcse)
        if (x.is_margin)
            s4 = sprintf("At a sample size of %s, the tightest equivalence bounds establishable with %2.0f%% power were +/-%5.3f (%3.1f percentage points)%s.", n_text, 100 * x.goal, x.mde, 100 * x.mde, conf_frag)
        else
            s4 = sprintf("At a sample size of %s, the minimum detectable %s at %2.0f%% power for %s was %5.3f (%s of %3.1f percentage points)%s.",
                n_text, est_short, 100 * x.goal, pape_claim_text(x.claim, x.sesoi, x.conf, rule), x.mde,
                (x.mde < 0 ? "a decrease" : "an increase"), 100 * abs(x.mde), conf_frag)
    }
    else if (x.kind == "n") {
        size_noun = (is_panel ? "number of units" : "total sample size")
        if (x.confirmed)
            s4 = sprintf("The required %s for %2.0f%% power for %s was %s, confirmed by a high-precision verification run (simulated power %5.3f, Monte Carlo SE %5.3f, %f replications).",
                size_noun, 100 * x.goal, pape_claim_text(x.claim, x.sesoi, x.conf, rule), n_text, x.power, x.mcse, x.nsim_confirm)
        else
            s4 = sprintf("The required %s for %2.0f%% power for %s was %s (simulated power %5.3f, Monte Carlo SE %5.3f, %f replications per search step).",
                size_noun, 100 * x.goal, pape_claim_text(x.claim, x.sesoi, x.conf, rule), n_text, x.power, x.mcse, x.nsim)
    }
    else
        s4 = sprintf("At a sample size of %s, simulated power for %s was %5.3f (Monte Carlo SE %5.3f; %f replications).",
            n_text, pape_claim_text(x.claim, x.sesoi, x.conf, rule), x.power, x.mcse, x.nsim)
    if (x.kind == "power" & x.has_sesoi)
        s5 = sprintf("Across replications, the probability of concluding a meaningful effect was %5.3f, of detection without meaningfulness %5.3f, of an inconclusive result %5.3f, and of equivalence %5.3f.",
            x.outcomes[1], x.outcomes[2], x.outcomes[3], x.outcomes[4])
    else s5 = ""
    if (d.route == "iv") est_desc = "two-step control-function estimation with stacked method-of-moments robust standard errors and delta-method Wald confidence intervals"
    else if (x.se == "robust") est_desc = "maximum likelihood with heteroskedasticity-robust (sandwich) standard errors and delta-method Wald confidence intervals"
    else est_desc = "maximum likelihood with delta-method Wald confidence intervals"
    keep = (x.separation == "keep")
    sepd = (x.separated == . ? 0 : x.separated)
    cell = (is_aie ? "focal-by-moderator cell" : "focal cell")
    if (x.kind == "power") {
        other = x.outcomes[cols(x.outcomes)] - (keep ? 0 : sepd)
        if (keep)
            s6 = sprintf("Estimation used %s; replications whose model fit failed (%3.1f%%) counted against the claim, and replications with separation (a %s with no events or no non-events; %3.1f%%) were retained with their Wald intervals.",
                est_desc, 100 * other, cell, 100 * sepd)
        else
            s6 = sprintf("Estimation used %s; replications without a maximum-likelihood estimate, from a failed fit or from separation (a %s with no events or no non-events), counted against the claim (%3.1f%% of replications%s).",
                est_desc, cell, 100 * x.outcomes[cols(x.outcomes)],
                (sepd > 0 ? sprintf(", %3.1f%% with separation", 100 * sepd) : ""))
    }
    else
        s6 = sprintf("Estimation used %s; replications without a maximum-likelihood estimate (failed fits%s) counted against the claim%s.",
            est_desc, (keep ? "" : sprintf(" and separation, a %s with no events or no non-events", cell)),
            (keep ? ", while replications with separation were retained with their Wald intervals" : ""))
    return(s1 + " " + s2 + " " + (s3 != "" ? s3 + " " : "") + s4 + " " + (s5 != "" ? s5 + " " : "") + s6)
}

string scalar pape_alpha_statement(struct pape_alpha scalar x)
{
    string scalar est_short, how, n_text, s1, s2, s3, s4, s5
    est_short = (x.estimand == "aie" ? "AIE" : "APE")
    how = (x.error == "minimize" ? "minimizing the weighted combined Type I and Type II error rate" : "balancing the cost-weighted Type I and Type II error rates")
    n_text = (x.dgp.route == "panel" ? sprintf("n = %f units", x.n) : sprintf("n = %f", x.n))
    s1 = sprintf("We chose the error rate for %s by %s (Maier & Lakens, 2022; Mudge et al., 2012), weighting a Type I error %g times a Type II error and taking prior odds of %g for the planning value (%s = %5.3f) against the claim boundary, using the powerape package for Stata (version %s).",
        pape_claim_text(x.claim, x.sesoi, pape_claim_levels(x.claim, x.alpha, ., ., "")[4], 0),
        how, x.cost, x.prior, est_short, x.target, pape_version())
    s2 = sprintf("The assumed data-generating process was %s.", pape_describe_dgp(x.dgp))
    s3 = sprintf("At %s the optimum within the conventional cap of alpha = %g was alpha = %6.4f: simulated power %5.3f (Monte Carlo SE %5.3f; %f replications), a Type II error rate of %5.3f, and a weighted combined error rate of %6.4f versus %6.4f at alpha = 0.05; error rates between %6.4f and %6.4f are indistinguishable from the optimum at this precision.",
        n_text, x.cap, x.alpha, x.power, x.mcse, x.nsim, x.beta, x.wcer, x.wcer_default, x.range[1], x.range[2])
    s4 = (x.capped ? sprintf("The unconstrained optimum, alpha = %6.4f (weighted combined error rate %6.4f), lies above the cap; because raising the error rate above the conventional level needs a justification of its own (Maier & Lakens, 2022), the capped value is reported.", x.alpha_uncapped, x.wcer_uncapped) : "")
    s5 = (x.has_size ? sprintf("A verification simulation at the claim boundary (true %s = %5.3f; %f replications) gave a realized error rate of %6.4f (Monte Carlo SE %6.4f) at the chosen alpha.", est_short, x.size_truth, x.size_nsim, x.size, x.size_mcse) : "")
    return(s1 + " " + s2 + " " + s3 + (s4 != "" ? " " + s4 : "") + (s5 != "" ? " " + s5 : ""))
}

// Print a long string wrapped at ~78 characters.
void pape_print_wrapped(string scalar s)
{
    string rowvector w
    string scalar line
    real scalar j
    w = tokens(s)
    line = ""
    for (j = 1; j <= cols(w); j++) {
        if (strlen(line) + strlen(w[j]) + 1 > 78 & line != "") {
            printf("{txt}%s\n", line)
            line = w[j]
        }
        else line = (line == "" ? w[j] : line + " " + w[j])
    }
    if (line != "") printf("{txt}%s\n", line)
}

end
