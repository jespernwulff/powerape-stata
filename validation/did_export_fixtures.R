## did_export_fixtures.R -- R reference values for the difference-in-differences
## route of the Stata port: fixed datasets with R's estimates (Tier B),
## calibration constants (Tier C), and the inputs and exact values of the
## enumeration oracle (Tier A). Written to validation/data/ as didfix_*.csv,
## did-reference.csv (case, item, value) and did-tierA-sizes.csv.
##
## The route is in the development version of the R package powerape (after
## 1.11.4). Point the environment variable POWERAPE_R_DEV at that source tree
## (it is loaded with pkgload), or install a powerape that has ape_dgp_did().
## Run from the repository root:
##     Rscript validation/did_export_fixtures.R
dev <- Sys.getenv("POWERAPE_R_DEV")
if (nzchar(dev)) {
  suppressMessages(pkgload::load_all(dev, quiet = TRUE))
} else {
  suppressMessages(library(powerape))
}
if (!exists("ape_dgp_did")) stop("this powerape has no difference-in-differences route")
if (!dir.exists("validation")) stop("run did_export_fixtures.R from the repository root")
out <- file.path("validation", "data")
ns <- asNamespace("powerape")
g <- function(nm) get(nm, envir = ns)
ref <- data.frame(case = character(), item = character(), value = numeric(),
                  stringsAsFactors = FALSE)
add <- function(case, item, value) {
  ref <<- rbind(ref, data.frame(case = case, item = item, value = as.numeric(value),
                                stringsAsFactors = FALSE))
}
se_of <- function(e, V) sqrt(drop(t(e$jac) %*% V %*% e$jac))
wcsv <- function(df, name) {
  write.csv(as.data.frame(lapply(df, function(v) sprintf("%.17g", v))),
            file.path(out, name), row.names = FALSE, quote = FALSE)
}

## ---- Tier B: fixed data -------------------------------------------------------
## Every fixture is written with the design matrix columns x2..xK (the
## constant is column 1), the outcome, the unit identifier (0 without one),
## the group, after-block and subgroup indicators, and the estimand's column.
one_fixture <- function(case, X, y, id, d, post, m, link, focal_col, w_col,
                        se_kinds, extra = NULL) {
  sub <- !is.null(m)
  ft <- g("fit_index_model")(X, y, link)
  stopifnot(ft$ok)
  b <- ft$fit$coefficients
  e <- if (sub) {
    g("datt_est")(b, list(X = X, focal_col = focal_col, w_col = w_col, m = m), link)
  } else {
    g("att_est")(b, X, link, focal_col)
  }
  for (j in seq_along(b)) add(case, paste0("b", j), b[j])
  for (j in seq_along(b)) add(case, paste0("jac", j), e$jac[j])
  add(case, "est", e$ape)
  if (sub) { add(case, "by0", e$by[1]); add(case, "by1", e$by[2]) }
  if ("cluster" %in% se_kinds) add(case, "se_cluster", se_of(e, g("vcov_cluster")(ft$fit, X, id, link)))
  if ("model" %in% se_kinds) add(case, "se_model", se_of(e, g("vcov_from_glmfit")(ft$fit)))
  if ("robust" %in% se_kinds) add(case, "se_robust", se_of(e, g("vcov_sandwich")(ft$fit, X, link)))
  fo <- g("fit_ols")(X, y)
  stopifnot(fo$ok)
  add(case, "lin_b", fo$coefficients[focal_col])
  if ("cluster" %in% se_kinds) add(case, "lin_se_cluster", sqrt(g("vcov_ols_cluster")(fo, X, id)[focal_col, focal_col]))
  if (any(c("model", "robust") %in% se_kinds)) add(case, "lin_se_hc1", sqrt(g("vcov_ols_hc1")(fo, X)[focal_col, focal_col]))
  add(case, "n", nrow(X)); add(case, "k", ncol(X))
  add(case, "focal_col", focal_col); add(case, "w_col", w_col)
  add(case, "sum_y", sum(y))
  df <- data.frame(y = y, id = if (is.null(id)) 0 else id, d = d, post = post,
                   m = if (sub) m else 0)
  Xc <- as.data.frame(X[, -1, drop = FALSE])
  names(Xc) <- paste0("x", seq_len(ncol(X))[-1])
  df <- cbind(df, Xc)
  if (!is.null(extra)) df <- cbind(df, extra)
  wcsv(df, paste0("didfix_", case, ".csv"))
  cat(sprintf("%-6s %5d rows, %2d design columns, estimate %.6f\n", case, nrow(X), ncol(X), e$ape))
}

## (1) the fixtures of the R package's regression tests against Stata
hf <- if (nzchar(dev)) file.path(dev, "tests", "testthat", "helper-stata-fixtures.R") else ""
if (!file.exists(hf)) stop("helper-stata-fixtures.R not found: set POWERAPE_R_DEV to the R source tree")
source(hf)
f <- stata_fixture("did")              # panel, 400 units, one period a side, a covariate
X <- cbind(1, f$d, as.numeric(f$period == 2), f$w, f$z)
one_fixture("did", X, f$y, f$id, f$d, as.numeric(f$period == 2), NULL, "logit", 4L, 4L,
            "cluster", extra = data.frame(period = f$period, w = f$w, z = f$z))
f <- stata_fixture("didx")             # cross-sections, four periods, treatment from period 3
X <- cbind(1, f$d, as.numeric(f$period == 2), as.numeric(f$period == 3),
           as.numeric(f$period == 4), f$w)
one_fixture("didx", X, f$y, NULL, f$d, as.numeric(f$period >= 3), NULL, "logit", 6L, 6L,
            c("model", "robust"), extra = data.frame(period = f$period, w = f$w))
f <- stata_fixture("sub")              # panel with a binary subgroup
post <- as.numeric(f$period == 2)
X <- cbind(1, f$d, post, f$w, f$m, f$m * f$d, f$m * post, f$m * f$w)
one_fixture("sub", X, f$y, f$id, f$d, post, f$m, "logit", 8L, 4L, "cluster",
            extra = data.frame(period = f$period, w = f$w))

## (2) studies drawn from the route itself
draw_fix <- function(d, n, seed) {
  g("with_seed")(seed, {
    xx <- g("draw_x")(d, n)
    y <- rbinom(length(xx$pr), 1, xx$pr)
    list(xx = xx, y = y)
  })
}
## probit panel, two periods before and two after, a covariate correlated with the group
d_rp <- set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                            rho = 0.4, n_pre = 2, n_post = 2,
                            covariates = list(pa_var("z", "normal")), correlation = 0.3,
                            signal = 0.2, n_int = 1e5), 0.08)
s <- draw_fix(d_rp, 500, 20261301L)
one_fixture("rp", s$xx$X, s$y, s$xx$id, s$xx$d, s$xx$post, NULL, "probit",
            s$xx$focal_col, s$xx$focal_col, "cluster")
## logit cross-sections with a subgroup, one period before and two after, a covariate
d_rs <- set_att(ape_dgp_did(model = "logit", baseline = c(0.30, 0.20), gap = c(0.10, 0.15),
                            trend = c(0.10, 0.05), design = "cross-sections",
                            n_pre = 1, n_post = 2,
                            subgroup = pa_var("female", "binary", p = 0.4),
                            covariates = list(pa_var("z", "normal")), signal = 0.15,
                            n_int = 1e5), c(0.10, 0.03))
s <- draw_fix(d_rs, 3000, 20261302L)
one_fixture("rs", s$xx$X, s$y, NULL, s$xx$d, s$xx$post, s$xx$m, "logit",
            s$xx$focal_col, s$xx$w_col, c("model", "robust"))

## ---- Tier C: calibration constants ---------------------------------------------
cal <- function(case, d) {
  cc <- d$did; r <- d$rates
  M <- length(cc$a)
  for (m in seq_len(M)) {
    sfx <- if (M == 1L) "" else paste0("_", m - 1L)
    for (nm in c("a", "b", "l", "th", "de")) add(case, paste0(nm, sfx), cc[[nm]][m])
    for (nm in c("c00", "c01", "c10", "cf", "c11")) add(case, paste0(nm, sfx), r[[nm]][m])
  }
  rd <- g("did_readings")(d)
  add(case, "reading_linear", rd$linear); add(case, "reading_index", rd$index)
  ta <- true_att(d)
  add(case, "true", ta[length(ta)])
  add(case, "signal_implied", d$signal_implied)
  if (d$k > 0L) for (j in seq_len(d$k)) add(case, paste0("gamma", j), d$gamma[j])
  bt <- g("did_beta")(d)
  for (j in seq_along(bt)) add(case, paste0("beta", j), bt[j])
}
## closed form: no covariates
cal("c1", set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20, trend = 0.10, rho = 0.4), 0.07))
cal("c2", set_att(ape_dgp_did(model = "probit", baseline = 0.20, gap = 0.20, trend = 0.10,
                              trends = "linear", analysis = "linear"), 0.10))
cal("c3", set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                              design = "cross-sections", n_pre = 2, n_post = 3), -0.05))
cal("c4", set_att(ape_dgp_did(model = "logit", baseline = c(0.30, 0.20), gap = c(0.10, 0.15),
                              trend = c(0.10, 0.05), rho = 0.5,
                              subgroup = pa_var("female", "binary", p = 0.4)), c(0.10, 0.03)))
cal("c5", set_att(ape_dgp_did(model = "probit", baseline = c(0.30, 0.05), gap = c(0.10, 0.05),
                              trend = c(0.10, 0.10), trends = "linear",
                              subgroup = pa_var("female", "binary", p = 0.5)), c(0.10, 0.03)))
## integration draw: covariates (the two packages draw their own sets)
cov2 <- list(pa_var("z", "normal"), pa_var("w", "binary", p = 0.4))
cal("c6", set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.10, trend = 0.08, rho = 0.4,
                              covariates = cov2, correlation = 0.3, signal = 0.2), 0.07))
cal("c7", set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                              trends = "linear", covariates = cov2, correlation = 0.3,
                              signal = 0.2), 0.07))
cal("c6big", set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.10, trend = 0.08, rho = 0.4,
                                 covariates = cov2, correlation = 0.3, signal = 0.2, n_int = 4e6), 0.07))
cal("c8", set_att(ape_dgp_did(model = "logit", baseline = c(0.30, 0.20), gap = c(0.10, 0.15),
                              trend = c(0.10, 0.05), rho = 0.5,
                              subgroup = pa_var("female", "binary", p = 0.4),
                              covariates = list(pa_var("z", "normal")), signal = 0.15), c(0.10, 0.03)))

## the copula's joint rate and the pilot's method-of-moments rho
add("pilot", "joint_30_05", g("did_joint_rate")(0.30, 0.5))
add("pilot", "joint_08_02", g("did_joint_rate")(0.08, 0.2))
add("pilot", "rho_solve", g("did_pilot_rho")(c(0.25, 0.35), c(0.11, 0.19), c(800, 760)))
## a pilot panel drawn from the route: 3,000 units, three pre-treatment periods
set.seed(20261303L)
G <- 3000L
grp <- rbinom(G, 1, 0.45)
a_i <- rnorm(G)
p_g <- c(0.22, 0.34)[grp + 1L]
pil <- do.call(rbind, lapply(1:3, function(t) {
  pr <- pnorm((qnorm(p_g) - sqrt(0.35) * a_i) / sqrt(1 - 0.35))
  data.frame(id = seq_len(G), period = t, g = grp, y = rbinom(G, 1, pr))
}))
dp <- suppressWarnings(ape_dgp_did_pilot(pil, outcome = "y", group = "g", id = "id",
                                         period = "period", trend = 0.05, model = "logit"))
for (nm in c("baseline", "gap", "p_treated", "rho", "control_change", "n_pairs"))
  add("pilot", nm, dp$pilot[[nm]])
add("pilot", "cf", dp$rates$cf)
wcsv(pil, "didfix_pilot.csv")

## ---- Tier A: exact power by enumeration -------------------------------------------
## Repeated cross-sections without covariates are saturated in the four
## group-by-period cells: given the cell sizes every study is one of
## prod(n_c + 1) outcomes with a known binomial probability, and the estimate
## and its standard error are closed forms in the four cell means. The exact
## conditional power is averaged over 1,000 draws of the multinomial cell
## sizes, which are written out so that the Stata oracle uses the same ones.
exact_did <- function(pc, sizes, link, analysis, conf = 0.95) {
  G <- if (link == "probit") pnorm else plogis
  Gi <- if (link == "probit") qnorm else qlogis
  gd <- if (link == "probit") dnorm else dlogis
  z <- qnorm(1 - (1 - conf) / 2)
  n <- sum(sizes[1, ])
  draws <- nrow(sizes)
  outp <- fail <- numeric(draws)
  for (r in seq_len(draws)) {
    nc <- sizes[r, ]
    if (any(nc == 0)) { fail[r] <- 1; next }
    dims <- nc + 1L
    E <- function(x, i) aperm(array(x, dim = c(dims[i], dims[-i])), order(c(i, setdiff(1:4, i))))
    ph <- lapply(nc, function(m) (0:m) / m)
    W <- E(dbinom(0:nc[1], nc[1], pc[1]), 1) * E(dbinom(0:nc[2], nc[2], pc[2]), 2) *
      E(dbinom(0:nc[3], nc[3], pc[3]), 3) * E(dbinom(0:nc[4], nc[4], pc[4]), 4)
    P <- lapply(1:4, function(i) E(ph[[i]], i))
    if (analysis == "linear") {
      est <- P[[4]] - P[[3]] - P[[2]] + P[[1]]
      v <- lapply(1:4, function(i) E(ph[[i]] * (1 - ph[[i]]) / nc[i], i))
      se <- sqrt(n / (n - 4) * (v[[1]] + v[[2]] + v[[3]] + v[[4]]))
      allsame <- (P[[1]] == P[[2]]) & (P[[2]] == P[[3]]) & (P[[3]] == P[[4]]) &
        (P[[1]] == 0 | P[[1]] == 1)
      ok <- se > 0 & !allsame
    } else {
      inner <- lapply(ph, function(x) x > 0 & x < 1)
      ok <- E(inner[[1]], 1) & E(inner[[2]], 2) & E(inner[[3]], 3) & E(inner[[4]], 4)
      li <- lapply(seq_along(ph), function(i) ifelse(inner[[i]], Gi(pmin(pmax(ph[[i]], 1e-12), 1 - 1e-12)), NA))
      L <- E(li[[3]], 3) + E(li[[2]], 2) - E(li[[1]], 1)
      est <- P[[4]] - G(L)
      vl <- lapply(1:3, function(i) E(ph[[i]] * (1 - ph[[i]]) / (nc[i] * gd(li[[i]])^2), i))
      se <- sqrt(E(ph[[4]] * (1 - ph[[4]]) / nc[4], 4) + gd(L)^2 * (vl[[1]] + vl[[2]] + vl[[3]]))
    }
    det <- ok & (est - z * se > 0)
    det[is.na(det)] <- FALSE
    outp[r] <- sum(W[det])
    fail[r] <- sum(W[!ok])
  }
  list(power = mean(outp), power_se = sd(outp) / sqrt(draws), failed = mean(fail),
       failed_se = sd(fail) / sqrt(draws))
}
set.seed(99L)
sizes <- t(rmultinom(1000L, 120L, rep(0.25, 4)))
write.csv(data.frame(n00 = sizes[, 1], n01 = sizes[, 2], n10 = sizes[, 3], n11 = sizes[, 4]),
          file.path(out, "did-tierA-sizes.csv"), row.names = FALSE, quote = FALSE)
for (cs in list(c("logit", "index"), c("probit", "index"), c("logit", "linear"))) {
  d <- set_att(ape_dgp_did(model = cs[1], baseline = 0.30, gap = 0.10, trend = 0.10,
                           design = "cross-sections", trends = cs[2], analysis = cs[2]), 0.30)
  pc <- unname(unlist(d$rates[c("c00", "c01", "c10", "c11")]))
  ex <- exact_did(pc, sizes, cs[1], cs[2])
  case <- paste0("A_", cs[1], "_", cs[2])
  add(case, "power", ex$power); add(case, "power_se", ex$power_se)
  add(case, "failed", ex$failed); add(case, "failed_se", ex$failed_se)
  for (j in 1:4) add(case, paste0("pc", j), pc[j])
  cat(sprintf("%-16s exact power %.4f (se %.4f), failed %.5f\n", case, ex$power, ex$power_se, ex$failed))
}

write.csv(transform(ref, value = sprintf("%.17g", value)), file.path(out, "did-reference.csv"),
          row.names = FALSE, quote = FALSE)
cat(sprintf("did-reference.csv: %d rows; powerape %s\n", nrow(ref),
            as.character(utils::packageVersion("powerape"))))
