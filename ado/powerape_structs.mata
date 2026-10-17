*! powerape_structs.mata -- shared structures (compiled after powerape_core)
version 16
mata:
mata set matastrict on

struct pape_fit {
    real colvector b
    real scalar ok, converged, rank
    real matrix V            // expected-information covariance
    real colvector eta, mu
}

struct pape_est {
    real scalar ape
    real colvector jac
}

struct pape_xdraw {
    real matrix X
    real scalar focal_col
    real colvector id, pr, start, y, d, m
    real matrix Zx, X1
    real scalar has_id, has_start
}

struct pape_sim {
    real colvector est, l, u, se, ok, sep, min_cell
}

struct pape_res {
    string scalar kind          // "power" | "n" | "mde" | "curve" | "robust" | "alpha"
    real scalar power, mcse
    real rowvector outcomes
    string rowvector outnames
    real scalar n_failed, separated, n_separated, min_cell
    string scalar claim, se, separation, model, estimand, direction
    string scalar engine, effects, fit   // who fitted the simulated studies (0.1.2): r(engine), r(effects), the "fit:" line
    real scalar sesoi, conf, alpha, n, nsim, target
    string scalar alpha_rule
    real scalar has_rule
    real scalar has_sesoi
    struct pape_sim scalar draws
    real scalar has_draws
    struct pape_dgp scalar dgp
    // ape_n / ape_mde
    real scalar goal, confirmed, confirm_requested, nsim_confirm, mde, is_margin
    real matrix history
    // curve
    real matrix results
}

struct pape_pdraw {
    real matrix x          // N x K regressor paths (focal, moderator, covariates)
    real matrix xbar       // n_units x K unit means over observed periods
    real colvector unit    // N x 1 unit index (unit-major)
    real colvector T_i     // n_units x 1 panel lengths
}

struct pape_cffit {
    real scalar ok, has_m, cf_col, cfm_col
    real colvector theta, cf, pi
    real matrix W, V_main, V
}

struct pape_ivdraw {
    real colvector d, m, v
    real matrix X, Zi
}

struct pape_alpha {
    real scalar alpha, beta, power, mcse, wcer, power_default, wcer_default
    real rowvector range
    real scalar capped, cap, alpha_uncapped, beta_uncapped, wcer_uncapped
    real scalar cost, prior
    string scalar error, claim, estimand, model, se
    real scalar sesoi, n, nsim, target
    real scalar has_size, size, size_mcse, size_nsim, size_truth, size_failed, size_separated
    real matrix grid
    struct pape_dgp scalar dgp
}

end
