## tierD_r.R -- R reference runs for the stochastic-parity tier (Tier D):
## matched designs simulated by powerape 1.11.3 at high replication counts.
## Writes validation/data/tierD-reference.csv (design, quantity, value, mcse)
## and the pilot data used by the empirical-route design
## (validation/data/pilot600.csv). Run from the repository root (about 16 min):
##     Rscript validation/tierD_r.R > validation/tierD_r.log
suppressMessages(library(powerape))
out <- file.path("validation", "data")
if (!dir.exists(out)) stop("run tierD_r.R from the repository root")
ref <- data.frame(design = character(), quantity = character(), value = numeric(),
                  mcse = numeric(), stringsAsFactors = FALSE)
add <- function(design, quantity, value, mcse = NA) {
  ref <<- rbind(ref, data.frame(design = design, quantity = quantity,
                                value = as.numeric(value), mcse = as.numeric(mcse),
                                stringsAsFactors = FALSE))
}
t0 <- Sys.time()
NS <- 4000

## D1: Design A world, minimum-effect claim
d1 <- set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                      covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
                      baseline = 0.30, signal = 0.15), 0.05)
p <- ape_power(d1, n = 2000, claim = "minimum", sesoi = 0.02, nsim = NS, seed = 101)
add("D1", "power", p$power, p$mcse)
add("D1", "detect_only", p$outcomes[["detect_only"]]); add("D1", "inconclusive", p$outcomes[["inconclusive"]])
p <- ape_power(d1, n = 2000, claim = "detect", nsim = NS, seed = 102)
add("D1", "power_detect", p$power, p$mcse)

## D2: logit, continuous focal, robust SEs
d3 <- set_ape(ape_dgp("logit", focal = pa_var("dose", "normal", mean = 2, sd = 1),
                      covariates = list(pa_var("z1", "normal"), pa_var("z2", "binary", p = 0.4)),
                      correlation = 0.2, baseline = 0.25, signal = 0.10), 0.03)
p <- ape_power(d3, n = 800, claim = "detect", nsim = NS, seed = 103, se = "robust")
add("D2", "power", p$power, p$mcse)
p <- ape_power(d3, n = 800, claim = "detect", nsim = NS, seed = 104)
add("D2", "power_model", p$power, p$mcse)

## D3: AIE binary x binary
d2 <- set_aie(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                      moderator = pa_var("female", "binary", p = 0.55), baseline = 0.30),
              0.04, main_focal = 0.10, main_moderator = 0.05)
p <- ape_power(d2, n = 4000, claim = "detect", nsim = NS, seed = 105)
add("D3", "power", p$power, p$mcse)
## D3b: AIE continuous x binary (logit)
d2b <- set_aie(ape_dgp("logit", focal = pa_var("x", "normal"),
                       moderator = pa_var("m", "binary", p = 0.5),
                       covariates = list(pa_var("z", "normal")), baseline = 0.35, signal = 0.1),
               0.03, main_focal = 0.05, main_moderator = 0.08)
p <- ape_power(d2b, n = 3000, claim = "detect", nsim = NS, seed = 106)
add("D3b", "power", p$power, p$mcse)

## D4: equivalence at truth 0
d4 <- set_ape(ape_dgp(focal = pa_var("t", "binary", p = 0.5), baseline = 0.30), 0)
p <- ape_power(d4, n = 2000, claim = "equivalence", sesoi = 0.05, nsim = NS, seed = 107)
add("D4", "power", p$power, p$mcse)

## D5: panel CRE probit
dp <- set_ape(ape_dgp_panel(focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                            covariates = list(pa_var("size", "normal", icc = 0.6)),
                            n_periods = 4, rho = 0.3, cre_share = 0.5,
                            baseline = 0.30, signal = 0.10), 0.08)
p <- ape_power(dp, n = 150, claim = "detect", nsim = NS, seed = 108)
add("D5", "power", p$power, p$mcse)
dpu <- set_ape(ape_dgp_panel(focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                             covariates = list(pa_var("size", "normal", icc = 0.6)),
                             n_periods = 6, retention = 0.8, rho = 0.3, cre_share = 0.5,
                             baseline = 0.30, signal = 0.10), 0.08)
p <- ape_power(dpu, n = 150, claim = "detect", nsim = NS, seed = 109)
add("D5u", "power", p$power, p$mcse)

## D6: IV, binary endogenous focal (probit first stage) and continuous (linear)
di <- set_ape(ape_dgp_iv(focal = pa_var("training", "binary", p = 0.4),
                         covariates = list(pa_var("size", "normal")),
                         instruments = pa_var("subsidy", "binary", p = 0.5),
                         endogeneity = 0.4, iv_strength = 0.2,
                         baseline = 0.30, signal = 0.10), 0.10)
p <- ape_power(di, n = 1500, claim = "detect", nsim = NS, seed = 110)
add("D6", "power", p$power, p$mcse)
di2 <- set_ape(ape_dgp_iv(focal = pa_var("rd", "normal", mean = 0, sd = 1),
                          instruments = pa_var("z1", "normal"),
                          endogeneity = 0.5, iv_strength = 0.3, baseline = 0.30), 0.08)
p <- ape_power(di2, n = 900, claim = "detect", nsim = NS, seed = 111)
add("D6c", "power", p$power, p$mcse)

## D7: empirical route on a shared pilot file
set.seed(600)
pilot <- data.frame(treat = rbinom(600, 1, 0.5), age = rnorm(600, 45, 12),
                    female = rbinom(600, 1, 0.55))
write.csv(pilot, file.path(out, "pilot600.csv"), row.names = FALSE)
de <- set_ape(ape_dgp_empirical(data = pilot, focal = "treat", baseline = 0.30, signal = 0.15), 0.05)
add("D7", "beta0", de$beta0); add("D7", "gamma1", de$gamma[1]); add("D7", "gamma2", de$gamma[2])
add("D7", "beta_focal", de$beta_focal)
p <- ape_power(de, n = 1500, claim = "minimum", sesoi = 0.02, nsim = NS, seed = 112)
add("D7", "power", p$power, p$mcse)

## D8: required n at the article's Design A/B anchors (seeded; confirmation at 4 x nsim)
nA <- ape_n(d1, power = 0.80, claim = "detect", nsim = 1500, seed = 113)
add("D8", "n_detect80", nA$n); add("D8", "power_at_n", nA$power, nA$mcse)
nB <- ape_n(d1, power = 0.80, claim = "minimum", sesoi = 0.02, nsim = 1500, seed = 114)
add("D8", "n_minimum80", nB$n); add("D8", "power_at_nmin", nB$power, nB$mcse)
d0 <- set_ape(ape_dgp(focal = pa_var("t", "binary", p = 0.5), baseline = 0.30), 0.10)
n0 <- ape_n(d0, power = 0.80, claim = "detect", nsim = 1500, seed = 115)
add("D8", "n_sat80", n0$n); add("D8", "power_at_nsat", n0$power, n0$mcse)

## D9: MDE at n = 712 (the four-tool anchor)
m <- ape_mde(d0, n = 712, claim = "detect", nsim = 1000, seed = 116)
add("D9", "mde", m$mde); add("D9", "power_at_mde", m$power, m$mcse)
m2 <- ape_mde(d4, n = 2000, claim = "equivalence", nsim = 1000, seed = 117)
add("D9", "margin", m2$mde); add("D9", "power_at_margin", m2$power, m2$mcse)

## D10: alpha rule and robustness sweep (power mode, pin ape)
p <- ape_power(d0, n = 712, claim = "detect", alpha = function(n) 0.05 * sqrt(500 / n),
               nsim = NS, seed = 118)
add("D10", "power_rule", p$power, p$mcse); add("D10", "alpha_rule", p$alpha)
rb <- ape_robust(d1, n = 700, claim = "detect", vary = list(baseline = c(0.20, 0.40)),
                 nsim = 2000, seed = 119, nmax = FALSE)
for (i in seq_len(nrow(rb$scenarios)))
  add("D10", paste0("robust_power_", i), rb$scenarios$power[i], rb$scenarios$mcse[i])

## D11: justified alpha on stored draws
p <- ape_power(d1, n = 2000, claim = "minimum", sesoi = 0.02, nsim = NS, seed = 120)
ja <- ape_alpha(p, cost = 4, prior = 1)
add("D11", "alpha_star", ja$alpha); add("D11", "wcer", ja$wcer); add("D11", "power_at_alpha", ja$power, ja$mcse)
pa <- power_at(p, alpha = c(0.01, 0.025, 0.05, 0.10))
for (i in 1:4) add("D11", paste0("power_at_", pa$alpha[i]), pa$power[i], pa$mcse[i])

write.csv(ref, file.path(out, "tierD-reference.csv"), row.names = FALSE)
cat("done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min;",
    nrow(ref), "rows\n")
