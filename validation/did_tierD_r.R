## did_tierD_r.R -- R reference runs for the stochastic-parity tier of the
## difference-in-differences route: matched designs simulated by the R
## package at 4,000 replications. Writes validation/data/did-tierD-reference.csv
## (design, quantity, value, mcse).
##
## The route is in the development version of the R package powerape (after
## 1.11.4). Point the environment variable POWERAPE_R_DEV at that source tree
## (it is loaded with pkgload), or install a powerape that has ape_dgp_did().
## Run from the repository root:
##     Rscript validation/did_tierD_r.R > validation/did_tierD_r.log
dev <- Sys.getenv("POWERAPE_R_DEV")
if (nzchar(dev)) {
  suppressMessages(pkgload::load_all(dev, quiet = TRUE))
} else {
  suppressMessages(library(powerape))
}
if (!exists("ape_dgp_did")) stop("this powerape has no difference-in-differences route")
out <- file.path("validation", "data")
if (!dir.exists(out)) stop("run did_tierD_r.R from the repository root")
ref <- data.frame(design = character(), quantity = character(), value = numeric(),
                  mcse = numeric(), stringsAsFactors = FALSE)
add <- function(design, quantity, value, mcse = NA) {
  ref <<- rbind(ref, data.frame(design = design, quantity = quantity,
                                value = as.numeric(value), mcse = as.numeric(mcse),
                                stringsAsFactors = FALSE))
}
addp <- function(design, quantity, p) {
  add(design, quantity, p$power, p$mcse)
  if (!is.null(p$est_mean)) add(design, paste0(quantity, "_estmean"), p$est_mean, p$est_mcse)
  add(design, paste0(quantity, "_failed"), p$outcomes[["failed"]])
}
sw <- suppressWarnings
t0 <- Sys.time()
NS <- 4000

## P1: logit panel, index world and analysis
p1 <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.10, trend = 0.08, rho = 0.4), 0.07)
addp("P1", "detect", sw(ape_power(p1, 1500, "detect", nsim = NS, seed = 301)))
p <- sw(ape_power(p1, 1500, "minimum", sesoi = 0.02, nsim = NS, seed = 302))
addp("P1", "minimum", p)
add("P1", "detect_only", p$outcomes[["detect_only"]]); add("P1", "inconclusive", p$outcomes[["inconclusive"]])

## P2: probit cross-sections with a covariate, model and robust SEs
p2 <- set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                          design = "cross-sections", covariates = list(pa_var("z", "normal")),
                          signal = 0.2), 0.08)
addp("P2", "detect", sw(ape_power(p2, 2500, "detect", nsim = NS, seed = 303)))
addp("P2", "detect_robust", sw(ape_power(p2, 2500, "detect", nsim = NS, seed = 304, se = "robust")))

## P3: linear analysis in a linear world, panel
p3 <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20, trend = 0.10, rho = 0.5,
                          trends = "linear", analysis = "linear"), 0.10)
addp("P3", "detect", sw(ape_power(p3, 900, "detect", nsim = NS, seed = 305)))

## P4: linear analysis in a linear world, cross-sections with a covariate, unequal shares
p4 <- set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                          design = "cross-sections", p_treated = 0.3, p_post = 0.6,
                          trends = "linear", analysis = "linear",
                          covariates = list(pa_var("z", "normal")), signal = 0.2), 0.08)
addp("P4", "detect", sw(ape_power(p4, 2500, "detect", nsim = NS, seed = 306)))

## P5: analyses that assume another scale than the world has (rejection rates)
p5a <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20, trend = 0.10, rho = 0.4,
                           trends = "linear", analysis = "index"), 0.10)
addp("P5a", "reject", sw(ape_power(p5a, 1000, "detect", nsim = NS, seed = 307)))
p5b <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20, trend = 0.10, rho = 0.4,
                           trends = "index", analysis = "linear"), 0.0667)
addp("P5b", "reject", sw(ape_power(p5b, 1000, "detect", nsim = NS, seed = 308)))

## P6: several periods
p6 <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20, trend = 0.10, rho = 0.5,
                          n_pre = 2, n_post = 2), 0.067)
addp("P6", "detect", sw(ape_power(p6, 1500, "detect", nsim = NS, seed = 309)))
p6x <- set_att(ape_dgp_did(model = "probit", baseline = 0.20, gap = 0.20, trend = 0.10,
                           design = "cross-sections", n_pre = 2, n_post = 3,
                           trends = "linear", analysis = "linear"), 0.08)
addp("P6x", "detect", sw(ape_power(p6x, 2000, "detect", nsim = NS, seed = 310)))

## P7: a binary subgroup (difference between the two subgroups' effects)
sub_args <- list(model = "logit", baseline = c(0.30, 0.20), gap = c(0.10, 0.15),
                 trend = c(0.10, 0.05), rho = 0.5,
                 subgroup = pa_var("female", "binary", p = 0.4))
p7 <- set_att(do.call(ape_dgp_did, sub_args), c(0.10, 0.03))
addp("P7", "detect", sw(ape_power(p7, 4000, "detect", nsim = NS, seed = 311)))
p7l <- set_att(do.call(ape_dgp_did, c(sub_args, list(trends = "linear", analysis = "linear"))),
               c(0.10, 0.03))
addp("P7l", "detect", sw(ape_power(p7l, 4000, "detect", nsim = NS, seed = 312)))

## P8: required n
nn <- sw(ape_n(p1, power = 0.8, claim = "detect", nsim = 1500, seed = 313))
add("P8", "n_p1_detect80", nn$n); add("P8", "n_p1_power", nn$power, nn$mcse)
nn <- sw(ape_n(p3, power = 0.8, claim = "detect", nsim = 1500, seed = 314))
add("P8", "n_p3_detect80", nn$n); add("P8", "n_p3_power", nn$power, nn$mcse)
nn <- sw(ape_n(p6, power = 0.8, claim = "detect", nsim = 1500, seed = 315))
add("P8", "n_p6_detect80", nn$n)

## P9: minimum detectable effect
m <- sw(ape_mde(p1, n = 1500, claim = "detect", nsim = 1000, seed = 316))
add("P9", "mde_p1", m$mde); add("P9", "mde_p1_power", m$power, m$mcse)
m <- sw(ape_mde(p3, n = 900, claim = "detect", nsim = 1000, seed = 317))
add("P9", "mde_p3", m$mde)
m <- sw(ape_mde(p7, n = 4000, claim = "detect", nsim = 1000, seed = 318))
add("P9", "mde_p7", m$mde)

## P10: robustness sweep (rho; the scale on which trends are parallel)
r <- sw(ape_robust(p1, n = 1500, claim = "detect", nsim = 2000, seed = 319,
                   vary = list(rho = c(0.2, 0.6)), nmax = FALSE))
for (i in 1:3) add("P10", paste0("rho_power_", i), r$scenarios$power[i], r$scenarios$mcse[i])
pg <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.20, trend = 0.08, rho = 0.4), 0.07)
r <- sw(ape_robust(pg, n = 1500, claim = "detect", nsim = 2000, seed = 320,
                   vary = list(trends = c("index", "linear")), nmax = FALSE))
for (i in 1:2) {
  add("P10", paste0("trends_power_", i), r$scenarios$power[i], r$scenarios$mcse[i])
  add("P10", paste0("trends_estmean_", i), r$scenarios$est_mean[i])
  add("P10", paste0("trends_valid_", i), as.numeric(r$scenarios$valid[i]))
}

## P11: justified alpha
p <- sw(ape_power(p1, 1500, "minimum", sesoi = 0.02, nsim = NS, seed = 321))
pa <- power_at(p, c(0.01, 0.025, 0.05, 0.10))
for (i in seq_len(nrow(pa))) add("P11", paste0("power_at_", pa$alpha[i]), pa$power[i], pa$mcse[i])
ja <- ape_alpha(p, cost = 4, prior = 1)
add("P11", "alpha_star", ja$alpha); add("P11", "wcer", ja$wcer)

write.csv(transform(ref, value = sprintf("%.10g", value), mcse = sprintf("%.10g", mcse)),
          file.path(out, "did-tierD-reference.csv"), row.names = FALSE, quote = FALSE)
cat(sprintf("done in %.1f min; %d rows; powerape %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(ref),
            as.character(utils::packageVersion("powerape"))))
