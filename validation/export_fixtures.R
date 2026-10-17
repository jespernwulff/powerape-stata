## export_fixtures.R -- fixed datasets and R reference values for the
## fixed-data parity tier of the Stata port's validation (Tier B).
## Run from the repository root, with the R package powerape 1.11.3 installed:
##     Rscript validation/export_fixtures.R   (writes to validation/data/)
suppressMessages(library(powerape))
if (!dir.exists("validation")) stop("run export_fixtures.R from the repository root")
out <- file.path("validation", "data")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ref <- data.frame(case = character(), item = character(), value = numeric(),
                  stringsAsFactors = FALSE)
add <- function(case, item, value) {
  ref <<- rbind(ref, data.frame(case = case, item = item, value = as.numeric(value),
                                stringsAsFactors = FALSE))
}
pa_se <- function(e, V) sqrt(drop(t(e$jac) %*% V %*% e$jac))
wcsv <- function(df, name) write.csv(df, file.path(out, name), row.names = FALSE)

## ---- A: probit, binary focal, two covariates ---------------------------------
set.seed(1)
n <- 2000
d <- rbinom(n, 1, 0.5); z1 <- rnorm(n); z2 <- rbinom(n, 1, 0.4)
y <- rbinom(n, 1, pnorm(-0.6 + 0.3 * d + 0.4 * z1 + 0.2 * z2))
X <- cbind(1, d, z1, z2)
ft <- powerape:::fit_index_model(X, y, "probit")
e <- powerape:::ape_est(ft$fit$coefficients, X, 2L, "probit", "binary")
add("A", "ape", e$ape)
add("A", "se_model", pa_se(e, powerape:::vcov_from_glmfit(ft$fit)))
add("A", "se_robust", pa_se(e, powerape:::vcov_sandwich(ft$fit, X, "probit")))
for (j in 1:4) add("A", paste0("b", j), ft$fit$coefficients[j])
for (j in 1:4) add("A", paste0("jac", j), e$jac[j])
wcsv(data.frame(y, d, z1, z2), "fixA.csv")

## ---- B: logit, continuous focal, two covariates ------------------------------
set.seed(2)
x <- rnorm(n, 2, 1); z1 <- rnorm(n); z2 <- rbinom(n, 1, 0.4)
y <- rbinom(n, 1, plogis(-1.2 + 0.35 * x + 0.5 * z1 + 0.3 * z2))
X <- cbind(1, x, z1, z2)
ft <- powerape:::fit_index_model(X, y, "logit")
e <- powerape:::ape_est(ft$fit$coefficients, X, 2L, "logit", "normal")
add("B", "ape", e$ape)
add("B", "se_model", pa_se(e, powerape:::vcov_from_glmfit(ft$fit)))
add("B", "se_robust", pa_se(e, powerape:::vcov_sandwich(ft$fit, X, "logit")))
for (j in 1:4) add("B", paste0("b", j), ft$fit$coefficients[j])
for (j in 1:4) add("B", paste0("jac", j), e$jac[j])
wcsv(data.frame(y, x, z1, z2), "fixB.csv")

## ---- C-F: AIE, four type pairs ------------------------------------------------
aie_case <- function(case, seed, link, ftype, mtype) {
  set.seed(seed)
  G <- if (link == "probit") pnorm else plogis
  dv <- if (ftype == "binary") rbinom(n, 1, 0.5) else rnorm(n)
  mv <- if (mtype == "binary") rbinom(n, 1, 0.45) else rnorm(n)
  z <- rnorm(n)
  y <- rbinom(n, 1, G(-0.4 + 0.35 * dv + 0.25 * mv + 0.3 * dv * mv + 0.4 * z))
  X <- cbind(1, dv, mv, dv * mv, z)
  ft <- powerape:::fit_index_model(X, y, link)
  e <- powerape:::aie_est(ft$fit$coefficients, X, link, ftype, mtype)
  add(case, "aie", e$ape)
  add(case, "se_model", pa_se(e, powerape:::vcov_from_glmfit(ft$fit)))
  for (j in 1:5) add(case, paste0("b", j), ft$fit$coefficients[j])
  for (j in 1:5) add(case, paste0("jac", j), e$jac[j])
  wcsv(data.frame(y, d = dv, m = mv, z), paste0("fix", case, ".csv"))
}
aie_case("C", 3, "probit", "binary", "binary")
aie_case("D", 4, "probit", "normal", "binary")
aie_case("E", 5, "logit", "binary", "normal")
aie_case("F", 6, "probit", "normal", "normal")

## ---- G: clustered (panel-style) probit with a Mundlak mean ----------------------
set.seed(7)
nu <- 100; Tp <- 5
id <- rep(seq_len(nu), each = Tp)
d <- rbinom(nu * Tp, 1, 0.5)
size <- rnorm(nu * Tp)
sbar <- ave(size, id)
a <- rnorm(nu, 0, 0.6)[id]
y <- rbinom(nu * Tp, 1, pnorm(-0.5 + 0.4 * d + 0.3 * size + 0.2 * sbar + a))
X <- cbind(1, d, size, sbar)
ft <- powerape:::fit_index_model(X, y, "probit")
e <- powerape:::ape_est(ft$fit$coefficients, X, 2L, "probit", "binary")
add("G", "ape", e$ape)
add("G", "se_model", pa_se(e, powerape:::vcov_from_glmfit(ft$fit)))
add("G", "se_cluster", pa_se(e, powerape:::vcov_cluster(ft$fit, X, id, "probit")))
for (j in 1:4) add("G", paste0("b", j), ft$fit$coefficients[j])
wcsv(data.frame(y, d, size, sbar, id), "fixG.csv")

## ---- H: control-function fixture (test-cf.R) ------------------------------------
set.seed(20260814)
n <- 1500
x1 <- rnorm(n); z1 <- rnorm(n); z2 <- rnorm(n)
cid <- rep(1:75, each = 20)
rho <- 0.5
v <- rnorm(n); u <- rho * v + sqrt(1 - rho^2) * rnorm(n)
y2 <- 0.3 + 0.4 * x1 + 0.5 * z1 + 0.3 * z2 + v
y0 <- as.integer(-0.3 + 0.5 * y2 + 0.4 * x1 + u > 0)
y2b <- as.integer(0.2 + 0.4 * x1 + 0.6 * z1 + 0.4 * z2 + v > 0)
y0b <- as.integer(-0.3 + 0.7 * y2b + 0.4 * x1 + u > 0)
m <- rbinom(n, 1, 0.4)
y0i <- as.integer(-0.4 + 0.6 * y2b + 0.3 * m + 0.4 * y2b * m + 0.4 * x1 + u > 0)
X1 <- cbind(1, x1, z1, z2)
f <- powerape:::cf_fit(y0, y2, Z = cbind(x1), X1 = X1, first = "linear")
e <- powerape:::cf_ape_est(f$theta, f$W, "normal", f$cf_col)
add("H_lin", "ape", e$ape); add("H_lin", "se", pa_se(e, f$V_main))
for (j in seq_along(f$theta)) add("H_lin", paste0("theta", j), f$theta[j])
for (j in seq_along(f$theta)) add("H_lin", paste0("se_theta", j), sqrt(diag(f$V_main))[j])
fc <- powerape:::cf_fit(y0, y2, Z = cbind(x1), X1 = X1, first = "linear", id = cid)
ec <- powerape:::cf_ape_est(fc$theta, fc$W, "normal", fc$cf_col)
add("H_lin_cl", "ape", ec$ape); add("H_lin_cl", "se", pa_se(ec, fc$V_main))
fb <- powerape:::cf_fit(y0b, y2b, Z = cbind(x1), X1 = X1, first = "probit")
eb <- powerape:::cf_ape_est(fb$theta, fb$W, "binary", fb$cf_col)
add("H_prb", "ape", eb$ape); add("H_prb", "se", pa_se(eb, fb$V_main))
for (j in seq_along(fb$theta)) add("H_prb", paste0("theta", j), fb$theta[j])
for (j in seq_along(fb$theta)) add("H_prb", paste0("se_theta", j), sqrt(diag(fb$V_main))[j])
## AIE with interact(m) and instrument-by-moderator first stage
X1m <- cbind(1, x1, m, z1, z2, z1 * m, z2 * m)
fi <- powerape:::cf_fit(y0i, y2b, m = m, Z = cbind(x1), X1 = X1m, first = "probit")
ei <- powerape:::cf_aie_est(fi$theta, fi$W, "binary", fi$cf_col, fi$cfm_col)
add("H_aie", "aie", ei$ape); add("H_aie", "se", pa_se(ei, fi$V_main))
for (j in seq_along(fi$theta)) add("H_aie", paste0("theta", j), fi$theta[j])
for (j in seq_along(fi$theta)) add("H_aie", paste0("se_theta", j), sqrt(diag(fi$V_main))[j])
wcsv(data.frame(x1, z1, z2, cid, y2, y0, y2b, y0b, m, y0i), "fixH.csv")

## ---- calibration constants (Tier C): closed-form and integration-draw worlds -----
d0 <- set_ape(ape_dgp(focal = pa_var("t", "binary", p = 0.5), baseline = 0.30), 0.10)
add("cal0", "beta0", d0$beta0); add("cal0", "beta_focal", d0$beta_focal); add("cal0", "p1", d0$p1)
d1 <- set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                      covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
                      baseline = 0.30, signal = 0.15), 0.05)
add("cal1", "beta0", d1$beta0); add("cal1", "gamma1", d1$gamma[1]); add("cal1", "beta_focal", d1$beta_focal)
add("cal1", "p1", d1$p1)
d2 <- set_aie(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                      moderator = pa_var("female", "binary", p = 0.55), baseline = 0.30),
              0.04, main_focal = 0.10, main_moderator = 0.05)
add("cal2", "beta_focal", d2$beta_focal); add("cal2", "beta_mod", d2$beta_mod); add("cal2", "beta_int", d2$beta_int)
d3 <- set_ape(ape_dgp("logit", focal = pa_var("dose", "normal", mean = 2, sd = 1),
                      covariates = list(pa_var("z1", "normal"), pa_var("z2", "binary", p = 0.4)),
                      correlation = 0.2, baseline = 0.25, signal = 0.10), 0.03)
add("cal3", "beta0", d3$beta0); add("cal3", "gamma1", d3$gamma[1]); add("cal3", "gamma2", d3$gamma[2])
add("cal3", "beta_focal", d3$beta_focal)
## high-precision versions (n_int = 4e6) to show convergence of the two draws
d1h <- set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                       covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
                       baseline = 0.30, signal = 0.15, n_int = 4e6), 0.05)
add("cal1h", "beta0", d1h$beta0); add("cal1h", "gamma1", d1h$gamma[1]); add("cal1h", "beta_focal", d1h$beta_focal)
d3h <- set_ape(ape_dgp("logit", focal = pa_var("dose", "normal", mean = 2, sd = 1),
                       covariates = list(pa_var("z1", "normal"), pa_var("z2", "binary", p = 0.4)),
                       correlation = 0.2, baseline = 0.25, signal = 0.10, n_int = 4e6), 0.03)
add("cal3h", "beta0", d3h$beta0); add("cal3h", "beta_focal", d3h$beta_focal)
## panel and IV calibration
dp <- set_ape(ape_dgp_panel(focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                            covariates = list(pa_var("size", "normal", icc = 0.6)),
                            n_periods = 4, rho = 0.3, cre_share = 0.5,
                            baseline = 0.30, signal = 0.10), 0.08)
add("calp", "beta0", dp$beta0); add("calp", "gamma1", dp$gamma[1]); add("calp", "xi1", dp$xi[1])
add("calp", "beta_focal", dp$beta_focal); add("calp", "var_a", dp$var_a); add("calp", "true_ape", true_ape(dp))
di <- set_ape(ape_dgp_iv(focal = pa_var("training", "binary", p = 0.4),
                         covariates = list(pa_var("size", "normal")),
                         instruments = pa_var("subsidy", "binary", p = 0.5),
                         endogeneity = 0.4, iv_strength = 0.2,
                         baseline = 0.30, signal = 0.10), 0.10)
add("cali", "beta0", di$beta0); add("cali", "gamma1", di$gamma[1]); add("cali", "pi_z1", di$pi_z[1])
add("cali", "kappa0", di$kappa0); add("cali", "beta_focal", di$beta_focal)
## closed-form helpers
add("util", "req_delta_1", powerape:::req_delta(0.8, qnorm(0.95), 1))
add("util", "req_delta_2", powerape:::req_delta(0.9, qnorm(0.95), 3))
add("util", "req_margin", powerape:::req_margin(0.8, 0.02, 0.01, qnorm(0.95)))
add("util", "frontier", alpha_frontier(cost = 4, power = 0.8))
add("util", "frontier_min", alpha_frontier(cost = 1, power = 0.8, claim = "minimum"))
add("util", "mz_r2", powerape:::mz_r2_logit(0.1 / 0.9))

write.csv(ref, file.path(out, "reference.csv"), row.names = FALSE)
cat("wrote", nrow(ref), "reference values and fixtures to", out, "\n")
