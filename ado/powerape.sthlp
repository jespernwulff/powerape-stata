{smcl}
{* *! version 0.2.0  13nov2026}{...}
{vieweralsosee "[R] probit" "help probit"}{...}
{vieweralsosee "[R] margins" "help margins"}{...}
{vieweralsosee "[R] cfprobit" "help cfprobit"}{...}
{vieweralsosee "[PSS] power" "help power"}{...}
{viewerjumpto "Syntax" "powerape##syntax"}{...}
{viewerjumpto "Description" "powerape##description"}{...}
{viewerjumpto "Variable specifications" "powerape##varspec"}{...}
{viewerjumpto "Subcommands" "powerape##subcommands"}{...}
{viewerjumpto "Options" "powerape##options"}{...}
{viewerjumpto "Difference-in-differences designs" "powerape##did"}{...}
{viewerjumpto "Engines" "powerape##engines"}{...}
{viewerjumpto "Separation and sparse cells" "powerape##separation"}{...}
{viewerjumpto "Stored results" "powerape##results"}{...}
{viewerjumpto "Examples" "powerape##examples"}{...}
{viewerjumpto "Installation and requirements" "powerape##install"}{...}
{viewerjumpto "Validation" "powerape##validation"}{...}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col:{bf:powerape} {hline 2}}Power analysis for average partial effects (APE) and average interaction effects (AIE){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}Describe the world (a data-generating process, DGP) and store it under a name:

{p 8 17 2}
{cmd:powerape dgp} {it:name}{cmd:,} {opt f:ocal(varspec)} {opt base:line(#)} [{opt mod:erator(varspec)} {opt cov:ariates(varspecs)} {opt corr:elation(#|matname)} {opt sig:nal(#)} {opt gam:ma(numlist)} {opt m:odel(probit|logit)} {opt nint(#)} {opt seed:int(#)} {opt replace} {opt nodis:play}]

{p 8 17 2}
{cmd:powerape dgp} {it:name} {ifin}{cmd:,} {opt emp:irical} {opt f:ocal(varname|varspec)} {opt base:line(#)} [{opt cov:ariates(varlist)} {opt sig:nal(#)} {opt gam:ma(numlist)} {opt m:odel(probit|logit)} ...]

{p 8 17 2}
{cmd:powerape dgp} {it:name}{cmd:,} {opt from:fit} {opt f:ocal(varname)} [{opt base:line(#)} ...]{space 8}(after {cmd:probit} or {cmd:logit})

{p 8 17 2}
{cmd:powerape dgp} {it:name}{cmd:,} {opt pan:el} {opt f:ocal(varspec)} {opt nper:iods(numlist)} {opt base:line(#)} [{opt pper:iods(numlist)} {opt ret:ention(#)} {opt rho(#)} {opt cre:share(#)} {opt mod:erator(varspec)} {opt cov:ariates(varspecs)} {opt corr:elation()} {opt sig:nal(#)} ...]

{p 8 17 2}
{cmd:powerape dgp} {it:name}{cmd:,} {opt iv} {opt f:ocal(varspec)} {opt inst:ruments(varspecs)} {opt endog:eneity(#)} {opt ivs:trength(#)} {opt base:line(#)} [{opt mod:erator(varspec)} {opt cov:ariates(varspecs)} {opt corr:elation()} {opt sig:nal(#)} ...]

{p 8 17 2}
{cmd:powerape dgp} {it:name}{cmd:,} {opt did} {opt base:line(# [#])} [{opt gap(# [#])} {opt trend(# [#])} {opt ptr:eated(#)} {opt des:ign(panel|cross-sections)} {opt rho(#)} {opt ppo:st(#)} {opt trends(index|linear)} {opt anal:ysis(index|linear)} {opt npre(#)} {opt npost(#)} {opt sub:group(varspec)} {opt cov:ariates(varspecs)} {opt corr:elation()} {opt sig:nal(#)} {opt gam:ma(numlist)} {opt m:odel(probit|logit)} ...]

{p 8 17 2}
{cmd:powerape dgp} {it:name} {ifin}{cmd:,} {opt didpilot} {opt out:come(varname)} {opt group(varname)} {opt trend(#)} [{opt id(varname)} {opt per:iod(varname)} {it:options of} {opt did}]

{pstd}Pin the effect in effect-size units (the inversion step):

{p 8 17 2}
{cmd:powerape setape} {it:name}{cmd:,} {opt t:arget(#)}

{p 8 17 2}
{cmd:powerape setaie} {it:name}{cmd:,} {opt t:arget(#)} {opt mainf:ocal(#)} {opt mainm:oderator(#)}

{p 8 17 2}
{cmd:powerape setatt} {it:name}{cmd:,} {opt t:arget(# [#])}{space 8}(difference-in-differences designs)

{pstd}Power, curves, required sample size, minimum detectable effect, robustness:

{p 8 17 2}
{cmd:powerape power} {it:name}{cmd:,} {opt n(#)} [{it:sim_options}]

{p 8 17 2}
{cmd:powerape curve} {it:name}{cmd:,} {opt n(numlist)} [{it:sim_options}]

{p 8 17 2}
{cmd:powerape n} {it:name} [{cmd:,} {opt pow:er(#)} {opt nrange(# #)} {opt max:iter(#)} {opt noconf:irm} {opt nsimc:onfirm(#)} {it:sim_options}]

{p 8 17 2}
{cmd:powerape mde} {it:name}{cmd:,} {opt n(#)} [{opt pow:er(#)} {opt dir:ection(positive|negative)} {opt mainf:ocal(#)} {opt mainm:oderator(#)} {opt max:iter(#)} {opt noconf:irm} {opt nsimc:onfirm(#)} {it:sim_options}]

{p 8 17 2}
{cmd:powerape robust} {it:name}{cmd:,} {opt n(#)} {opt vary(input(values) ...)} [{opt pin(ape|coefficients)} {opt mode(power|mde)} {opt pow:er(#)} {opt grid:points(#)} {opt nonmax} {opt nmaxp:ower(#)} {opt dir:ection()} {it:sim_options}]

{pstd}Reporting and justified error rates:

{p 8 17 2}
{cmd:powerape statement} [{it:resname}]

{p 8 17 2}
{cmd:powerape alpha} [{it:resname}] [{cmd:,} {opt cost(#)} {opt prior(#)} {opt err:or(minimize|balance)} {opt cap(#)} {opt size} {opt nsims:ize(#)} {opt seed(#)} {opt n:ame(resname)}]

{p 8 17 2}
{cmd:powerape powerat} [{it:resname}]{cmd:,} {opt alp:ha(numlist)}

{p 8 17 2}
{cmd:powerape frontier} [{cmd:,} {opt cost(#)} {opt pow:er(#)} {opt prior(#)} {opt cl:aim(detect|minimum)}]

{pstd}Housekeeping:

{p 8 17 2}
{cmd:powerape} {cmd:describe} {it:name} | {cmd:trueape} {it:name} | {cmd:trueaie} {it:name} | {cmd:trueatt} {it:name} | {cmd:list} | {cmd:drop} {it:name}|{cmd:_all} | {cmd:engine} [{cmd:,} {opt eng:ine(stata|mata)} {opt eff:ects(closed|margins)}] | {cmd:version}

{synoptset 26 tabbed}{...}
{synopthdr:sim_options}
{synoptline}
{synopt:{opt cl:aim(minimum|detect|equivalence)}}claim to power; default {cmd:minimum} ({cmd:detect} for {cmd:mde}){p_end}
{synopt:{opt ses:oi(#)}}smallest effect size of interest, in APE units{p_end}
{synopt:{opt alp:ha(#)}}the claim's error rate; default 0.05{p_end}
{synopt:{opt alphar:ule(exp)}}sample-size rule: a Stata expression in {cmd:@n}, e.g. {cmd:alpharule(0.05*sqrt(1000/@n))}{p_end}
{synopt:{opt conf(#)}}interval level override for the claim{p_end}
{synopt:{opt nsim(#)}}replications (power 1000, curve 1000, n 1500, mde 1000, robust 800){p_end}
{synopt:{opt seed(#)}}seed; the ambient RNG state is restored afterwards{p_end}
{synopt:{opt se(model|robust)}}model-based or heteroskedasticity-robust SEs (cross-sectional routes; ignored where the design fixes the SEs: panels, IV, and the linear analysis of a difference-in-differences design){p_end}
{synopt:{opt sep:aration(fail|keep)}}scoring of replications with a separated identifying cell{p_end}
{synopt:{opt eng:ine(stata|mata)}}who fits the simulated studies (see {help powerape##engines:Engines}) for this command only; default {cmd:stata}, or the session default set by {cmd:powerape engine}{p_end}
{synopt:{opt eff:ects(closed|margins)}}closed-form delta-method effects or {cmd:margins}/{cmd:ginteff} per replication, for this command only; default {cmd:closed}, or the session default{p_end}
{synopt:{opt n:ame(resname)}}store the result under {it:resname} (default {cmd:_last}){p_end}
{synopt:{opt nodr:aws}}do not keep the per-replication draws ({cmd:power} only){p_end}
{synopt:{opt nodis:play}}suppress output{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:powerape} is a Stata port of the R package {bf:powerape} (version 1.12.0):
simulation-based power analysis in which the effect size of interest is the
average partial effect (APE) or the average interaction effect (AIE) of a
probit or logit model. Targets are stated in probability points; the package
calibrates the data-generating process (baseline rate plus nuisance signal),
inverts the target into coefficients, checks feasibility, and computes power
for three confidence-interval claims (Riesthuis, 2024): {cmd:detect} (the
interval excludes zero in the hypothesized direction), {cmd:minimum} (the
lower bound exceeds the SESOI) and {cmd:equivalence} (the interval lies
within +/- SESOI). It provides power curves, required sample sizes with a
high-precision confirmation stage, minimum detectable effects, robustness
sweeps with the insurance-premium n_max (Hancock and Feng, 2025), citable
power statements, and justified error rates (a sample-size rule, or the
error-cost optimum of Maier and Lakens, 2022).

{pstd}
Routes: parametric covariates with Gaussian-copula dependence; covariates
resampled from pilot data in memory ({cmd:empirical}); a fitted {cmd:probit}
or {cmd:logit} pilot model ({cmd:fromfit}); correlated random effects panel
probit with Mundlak means, unit-clustered SEs and unbalanced panels
({cmd:panel}); and an endogenous focal variable estimated by control-function
probit, the estimator of Stata's {cmd:cfprobit} ({cmd:iv}). A sixth route,
{cmd:did}, plans a difference-in-differences study with a binary outcome: its
target is the effect on the treated in probability points (see
{help powerape##did:Difference-in-differences designs}).

{pstd}
Every simulated study is, by default, fit by Stata's own estimation commands
({cmd:probit}, {cmd:logit}, {cmd:cfprobit}; see {help powerape##engines:Engines}).


{marker varspec}{...}
{title:Variable specifications}

{pstd}
A {it:varspec} is {it:name}{cmd:,} [{cmd:binary} {cmd:p(#)} | {cmd:normal} {cmd:mean(#)} {cmd:sd(#)}] [{cmd:icc(#)}].
The type is inferred when omitted: {cmd:p()} means binary, otherwise normal
(default mean 0, sd 1). {cmd:icc()} is the within-unit persistence on the
latent scale and matters only for {cmd:panel} DGPs ({cmd:icc(1)} =
time-constant). Several variables are given as parenthesized specs:
{cmd:covariates((age, normal mean(45) sd(12)) (female, binary p(0.55)))}.


{marker subcommands}{...}
{title:Subcommands}

{phang}{cmd:dgp} describes the world. {opt baseline()} is the average
P(Y=1) with the focal (and moderator) at their reference values (0 for
binary, the mean for continuous); {opt signal()} is the nuisance covariates'
share of the latent-index variance at unit error variance (the McKelvey-Zavoina
R-squared under probit; a weaker world under logit, see the R documentation of
{cmd:ape_dgp}); {opt gamma()} gives explicit nuisance coefficients instead;
{opt correlation()} is a scalar exchangeable Gaussian-copula correlation or
a matrix over (focal, moderator, covariates). The panel route is
probit-only; {opt nperiods()} is a scalar (balanced) or a vector of lengths
with {opt pperiods()} probabilities, or a scalar with {opt retention()} for
monotone attrition; {opt rho()} is the latent unit-effect share and
{opt creshare()} the fraction loaded on the Mundlak means. The IV route
takes {opt endogeneity()} (the error correlation) and {opt ivstrength()}
(the instruments' share of the focal's (latent) variance); a binary focal
gets a probit first stage, a continuous one a linear first stage.

{phang}{cmd:setape} / {cmd:setaie} solve the focal (and interaction)
coefficients so the true APE / AIE equals {opt target()}; {cmd:setaie}
anchors the two main effects by their conditional-at-reference APEs.
Infeasible targets are refused with the attainable range.

{phang}{cmd:setatt} pins the effect on the treated of a
difference-in-differences design ({cmd:dgp ..., did}): {opt target()} is the
treated group's rate after treatment minus the rate it would have had
without, in probability points. A design with {opt subgroup()} takes two
numbers, the effect in the reference subgroup and then in the other one; its
estimand is their difference. {cmd:trueatt} recomputes the effect(s) from the
stored design.

{phang}{cmd:power} simulates power at {opt n()} and prints the outcome
distribution (minimum / detect-only / inconclusive / equivalence / failed).
{cmd:curve} does so over a grid. {cmd:n} searches the required sample size
for {opt power()} (default 0.90) and confirms it at {opt nsimconfirm()}
(default 4 x nsim), never trimming n on overshoot. {cmd:mde} searches the
smallest effect that reaches {opt power()} (default 0.80) at a fixed n (for
{cmd:equivalence}: the tightest margin at the pinned truth). {cmd:robust}
re-runs the analysis over a factorial grid of contextual inputs
({cmd:baseline}, {cmd:signal}, {cmd:correlation}; {cmd:rho}, {cmd:retention};
{cmd:endogeneity}, {cmd:ivstrength}) and reports the worst scenario and
n_max. Difference-in-differences designs can vary {cmd:baseline}, {cmd:gap}
and {cmd:trend} (not with a subgroup), {cmd:rho} (panels), {cmd:ptreated},
{cmd:signal}, and the scale on which untreated trends are parallel,
{cmd:vary(trends(index linear))}. For panel DGPs {opt n()} counts units. Sweeps over {cmd:empirical} and {cmd:fromfit} DGPs recalibrate from the pilot rows stored in the DGP (the data need not be in memory); parametric, panel and IV DGPs are rebuilt from their options. The engine subcommands hold and restore the user's e() around the simulated fits.

{phang}{cmd:statement} prints a citable methods paragraph for the last (or
named) {cmd:power}, {cmd:n}, {cmd:mde} or {cmd:alpha} result.
{cmd:alpha} chooses the error rate on the stored draws of a {cmd:power}
result (weighted combined error rate; cap 0.05); {cmd:powerat} re-thresholds
those draws at other error rates; {cmd:frontier} is the closed-form level at
the detectability frontier.

{phang}{cmd:version} prints the version of the ado-file and of the Mata
library {cmd:lpowerape.mlib} with their locations and warns when they differ
(a stale library after an update); it stores {cmd:r(version)},
{cmd:r(ado_version)}, {cmd:r(mlib_version)} and {cmd:r(match)}. {cmd:engine}
without options shows the session defaults; with options it sets them.
{opt engine()} or {opt effects()} given on a simulation command apply to that
command only and leave the session defaults unchanged.


{marker did}{...}
{title:Difference-in-differences designs}

{pstd}
{cmd:powerape dgp} {it:name}{cmd:, did} describes a difference-in-differences
study with one treatment date: a treated and a control group observed before
and after the treated group receives the intervention, in one period each by
default or in {opt npre()} periods before and {opt npost()} periods after. The
target is the {bf:average treatment effect on the treated (ATT)} after
treatment, in probability points, pinned with {cmd:powerape setatt}. The
outcome is regressed on the group indicator, period indicators, the treatment
indicator (treated group after the treatment date) and any covariates.

{pstd}
{bf:The world} is described on the probability scale by three rates:
{opt baseline()} (control group, before), {opt gap()} (treated minus control,
before; default 0) and {opt trend()} (the control group's change from before
to after; default 0). With several periods the rates are the same in every
period before and in every period after the treatment date, and so is the
effect. What the treated group's rate would have been after the treatment
date {it:without} treatment depends on the scale on which the two groups'
untreated trends are parallel, set by {opt trends()}: {cmd:index} (default)
means parallel on the probit/logit index, the assumption under which the
probit/logit analysis identifies the ATT (Wooldridge, 2023; Puhani, 2012);
{cmd:linear} means parallel in probabilities, the assumption of the linear
analysis (with covariates: parallel on average over each group's covariates).
{opt ptreated()} is the treated share (default 0.5). Covariates are constant
over time; {opt correlation()} runs over (group, covariates), or over (group,
subgroup, covariates).

{pstd}
{bf:The analysis} each simulated study receives is set by {opt analysis()}:
{cmd:index} (default) is a pooled probit or logit, with the ATT as the average
partial effect of the treatment indicator over the treated post-period
observations and a delta-method standard error (in Stata:
{cmd:margins, dydx(}{it:w}{cmd:) subpop(if }{it:w}{cmd: == 1)}); {cmd:linear} is a linear
probability model fitted by ordinary least squares, with the ATT as the
coefficient of the treatment indicator (the standard difference-in-differences
estimate). Each analysis is unbiased only in the world whose scale it assumes.
When {opt trends()} and {opt analysis()} name the same scale, {cmd:powerape power}
reports power. When they differ, the analysis centres on another number than
the ATT: {cmd:power} and {cmd:curve} then report a {bf:rejection rate} together
with the mean of the estimates, and {cmd:n}, {cmd:mde}, {cmd:powerat},
{cmd:alpha} and {cmd:statement} refuse the design or the result. With two
periods and no covariates the data cannot tell the two scales apart, and the
same four rates imply different effects under each (control .20 to .30 and
treated .40 to .60 is an effect of .100 under linear trends, .074 on the
probit index and .067 on the logit index): the printed design shows both
readings. The linear analysis always centres on the cross difference of the
four rates.

{pstd}
{bf:Data structure.} {cmd:design(panel)} (default) follows the same units over
all periods (a balanced panel); {opt n()} counts units and standard errors are
clustered by unit. Outcomes of one unit are dependent through a Gaussian
copula with latent correlation {opt rho()} that leaves the marginal rates
exact under either link (under probit this is the random-effects probit with
rho = Var(c)/(Var(c) + 1)). {cmd:design(cross-sections)} draws independent
observations, each after the treatment date with probability {opt ppost()}
(default {it:npost}/({it:npre} + {it:npost})) and then in one of that block's
periods with equal probability; {opt n()} counts all observations and
{opt se()} chooses model-based or robust standard errors for the probit/logit
analysis.

{pstd}
{bf:More than two periods.} The analysis includes an indicator for every
period but the first and one treatment indicator, so it estimates one effect,
not one per period. In a panel the dependence between two outcomes of a unit
is the same however far apart the periods are (one shared unit component).
That is the favourable case for adding periods: the gain is larger than under
dependence that fades with the distance between periods, which this version
does not offer. (With rates of .20 and .30 in the control group, .40 and .60
in the treated group and {cmd:rho(0.5)}, two periods on each side instead of
one cut the variance of the linear panel estimate to 51%; if the latent
correlation halved with every period of distance it would be 71%.) Treat the
precision bought by extra periods as an upper bound.

{pstd}
{bf:Comparing two subgroups.} With {opt subgroup(name, p(#))} both subgroups
are treated in the treated group and the target is the {bf:difference between the two subgroups' effects on the treated},
in probability points (the estimand of a triple difference that compares, say,
the effect on women with the effect on men; Caron, 2025). {opt baseline()},
{opt gap()} and {opt trend()} may then hold two numbers, for the reference
subgroup ({it:name} = 0) and for {it:name} = 1, and {cmd:setatt} takes the two
effects. The probit/logit analysis is fully interacted with the subgroup and
the estimate is the difference of the two subgroups' average partial effects
of treatment ({cmd:margins, dydx() subpop() over()}); the linear analysis
estimates the coefficient of the subgroup-by-treatment product. Trends are
parallel within each subgroup, on the scale {opt trends()} names. The same
eight rates imply different subgroup differences on the two scales (the design
print shows both); a difference of zero in probability points is not a zero
interaction coefficient when the subgroups have different rates; and with
equal subgroups a difference as large as the effect itself needs about four
times the units the effect needs. {cmd:mde} searches the difference, holding
the reference subgroup's effect at its pinned value.

{pstd}
{bf:The linear analysis.} A linear probability model is heteroskedastic by
construction, so its standard errors are always robust and {opt se()} is
ignored: HC1 in cross-sections ({cmd:regress, vce(robust)}) and clustered by
unit in panels ({cmd:regress, vce(cluster)}, with its small-sample factor).
Ordinary least squares has no separation, so no simulated study is discarded
for an all-zero or all-one cell; the sparse-cell warning still applies. In a
balanced panel a time-constant covariate changes neither the estimate nor its
clustered standard error beyond the degrees of freedom, so {opt signal()} buys
the linear panel analysis no precision; in cross-sections it does.

{pstd}
{bf:Interval coverage in small samples.} Under the index analysis the effect
on the treated is a nonlinear function of the estimated rates, and its
delta-method interval (the one {cmd:margins} reports) covers slightly less
than its nominal level when the group-by-period cells are small: for a 95%
interval at a true effect of zero, 94.0% with 125 observations per cell,
94.5% with 250, 94.7% with 500 and 94.9% with 1,000 under logit (exact values
in the design without covariates). The simulated studies use this interval,
so the power reported is the power of the analysis a study will run.

{pstd}
{bf:A design from pre-treatment pilot data.} {cmd:powerape dgp} {it:name}{cmd:, didpilot}
estimates, from data in memory observed {it:before} the treatment date, the
inputs such data can inform and builds the design: {opt baseline()} (the
control group's rate; {opt group()} is the 0/1 indicator of the group that
will be treated), {opt gap()}, {opt ptreated()} and, with {opt id()} and two
or more periods per unit, {opt rho()} by the method of moments (the value at
which the model's probability of the outcome in two periods of one unit equals
the observed share of such pairs). All pre-treatment periods are pooled.
Pre-treatment data cannot show how the control group will change from before
to after treatment, so {opt trend()} is required; with {opt period()} the
control group's average change between pilot periods is reported as a guide.
Without {opt id()} the design is repeated cross-sections. Naming
{opt baseline()}, {opt gap()}, {opt ptreated()}, {opt rho()} or {opt design()}
overrides the estimate; covariates are not taken from the pilot, and subgroup
designs are not supported. The estimates are used as if they were known: vary
them with {cmd:powerape robust}, in particular for a small pilot (a note is
printed below 500 rows).

{pstd}
{bf:Limits of this version.} One treated group, one treatment date, treatment
assigned at the unit level, balanced panels, effects that are the same in
every period after treatment, covariates that are constant over time and not
interacted with period or treatment. Staggered adoption, period-specific
effects, pre-trend tests and treatment assigned to whole clusters are not
available (unit clustering would overstate power under cluster-level
assignment).


{marker engines}{...}
{title:Engines}

{pstd}
{cmd:engine(stata)} (the default) writes every simulated study to a
temporary frame and fits it with {cmd:probit} or {cmd:logit}
({cmd:vce(robust)} for {cmd:se(robust)}, {cmd:vce(cluster)} on the panel
route, {cmd:asis} for {cmd:separation(keep)}) or with {cmd:cfprobit} on the
IV route. With {cmd:effects(closed)} the APE/AIE and its delta-method
standard error are then computed in Mata from e(b) and e(V) by closed-form
Jacobians; these equal {cmd:margins, dydx()} and {cmd:ginteff} to numerical
precision (validation tier B). With {cmd:effects(margins)} the effects come
from {cmd:margins} / {cmd:ginteff} themselves, which is several times slower.
On the IV route the APE is {cmd:margins, dydx()} after {cmd:cfprobit}, and
the AIE is the difference of {cmd:margins, dydx(}{it:focal}{cmd:) at(}{it:moderator}{cmd:=(0 1))}:
the focal-by-moderator term is entered in factor form ({cmd:1.}{it:d}{cmd:#c.}{it:m},
the same fit as {cmd:c.}{it:d}{cmd:#c.}{it:m}) so that {cmd:margins} can
differentiate the focal, and the control function is held at its estimated
value, as in {cmd:cfprobit}'s average structural function. (Setting the
endogenous focal through {cmd:at()} instead would let the control function
move and is not the AIE.) Both agree with the closed forms to 1e-7
(validation tier B). Inference follows Stata's conventions (observed
information, Stata's small-sample factors).

{pstd}
For difference-in-differences designs {cmd:engine(stata)} fits {cmd:logit} or
{cmd:probit} with {cmd:vce(cluster)} in panels (model-based or
{cmd:vce(robust)} in cross-sections) and computes the effect on the treated
in closed form, or with {cmd:effects(margins)} by
{cmd:margins, dydx(}{it:w}{cmd:) subpop(if }{it:w}{cmd: == 1)} (adding {cmd:over(}{it:m}{cmd:)} and
taking the difference for a subgroup design); the two agree to 1e-16 on fixed
data. The linear analysis is {cmd:regress} with {cmd:vce(robust)} or
{cmd:vce(cluster)}; {opt effects()} does not apply to it. {cmd:engine(mata)}
uses Fisher scoring with the cluster sandwich, and ordinary least squares with
the same HC1 and cluster factors as {cmd:regress}.

{pstd}
{cmd:engine(mata)} is a self-contained Mata engine that reproduces the R
package's conventions (Fisher scoring with expected information, HC0 and
G/(G-1) cluster sandwiches, the stacked control-function covariance). It is
about 10-30 times faster and is the right choice for long searches and
sweeps; the two engines agree within Monte Carlo error (tiers A and D).

{pstd}
Every simulation command ({cmd:power}, {cmd:curve}, {cmd:n}, {cmd:mde},
{cmd:robust}) prints a line {cmd:fit:} that names who fitted the simulated
studies, for example {cmd:fit: probit via engine(stata), effects(closed)},
{cmd:fit: probit, vce(cluster) via engine(stata), effects(closed)} on the
panel route, {cmd:fit: cfprobit via engine(stata), effects(margins)} on the
IV route, or {cmd:fit: Mata engine (Fisher scoring)}; {cmd:asis} appears
under {cmd:separation(keep)}. The same text is stored in {cmd:r(fit)}, with
its parts in {cmd:r(engine)} and {cmd:r(effects)}. {cmd:powerape engine}
without options shows the session default.


{marker separation}{...}
{title:Separation and sparse cells}

{pstd}
A simulated study may contain an identifying cell (a level of a binary focal
variable, a focal-by-moderator cell, or for difference-in-differences designs a
group-by-before/after cell, within each subgroup if there is one) with no events or no non-events, where
the maximum-likelihood effect does not exist. By default such replications
count as failed ({cmd:separation(fail)}), as Stata's {cmd:probit} refuses the
fit; {cmd:separation(keep)} retains them with their degenerate Wald intervals
(R's glm behaviour). {cmd:engine(mata)} keeps the Fisher-scoring fit;
{cmd:engine(stata)} fits with {cmd:asis} and keeps the separated coefficient,
whose variance row Stata sets to zero (its standard error prints as missing),
so the delta-method interval collapses onto the identified cells, which is the
same limit: for a binary focal with an empty arm the estimate is minus the
other arm's event rate, with that arm's binomial standard error. For
difference-in-differences designs this holds under {cmd:engine(mata)} only:
under {cmd:engine(stata)} the {cmd:asis} fit is kept when the separated cell
is the treated group's post-period cell (one coefficient without a variance),
but when another group-by-before/after cell is separated several coefficients
diverge together, Stata's fit does not converge and the study still counts as
failed, so use {cmd:engine(mata)} with {cmd:separation(keep)} there. Terms omitted
for collinearity still count as failures under either engine. The share with
separation and the average smallest cell count are reported and warned about
(more than 1% separated, or fewer than 10 events in the smallest cell).


{marker results}{...}
{title:Stored results}

{pstd}{cmd:powerape power} stores in {cmd:r()}: scalars {cmd:r(power)},
{cmd:r(mcse)}, {cmd:r(n)}, {cmd:r(nsim)}, {cmd:r(alpha)}, {cmd:r(conf)},
{cmd:r(sesoi)}, {cmd:r(target)}, {cmd:r(separated)}, {cmd:r(min_cell)},
{cmd:r(n_failed)}; macros {cmd:r(claim)}, {cmd:r(estimand)}, {cmd:r(model)},
{cmd:r(se)}, {cmd:r(separation)}, {cmd:r(dgp)}, {cmd:r(name)}, and {cmd:r(fit)},
{cmd:r(engine)}, {cmd:r(effects)} (who fitted the simulated studies, as on the
{cmd:fit:} line of the output); matrix
{cmd:r(outcomes)}. {cmd:curve} adds {cmd:r(results)}; {cmd:n} adds
{cmd:r(goal)}, {cmd:r(confirmed)}, {cmd:r(nsim_confirm)}, {cmd:r(history)};
{cmd:mde} adds {cmd:r(mde)}, {cmd:r(is_margin)}, {cmd:r(direction)};
{cmd:robust} stores {cmd:r(scenarios)}, {cmd:r(worst)}, {cmd:r(nmax)},
{cmd:r(nmax_power)}, {cmd:r(nmax_scenario)}, {cmd:r(nmax_candidates)},
{cmd:r(fit)}, {cmd:r(engine)}, {cmd:r(effects)};
{cmd:alpha} stores {cmd:r(alpha)}, {cmd:r(power)}, {cmd:r(beta)},
{cmd:r(wcer)}, {cmd:r(range_lo)}, {cmd:r(range_hi)}, {cmd:r(capped)},
{cmd:r(alpha_uncapped)}, {cmd:r(size)} and {cmd:r(grid)}; {cmd:dgp},
{cmd:setape}, {cmd:setaie} and {cmd:describe} store {cmd:r(beta0)},
{cmd:r(beta_focal)}, {cmd:r(beta_mod)}, {cmd:r(beta_int)}, {cmd:r(target)},
{cmd:r(p1)}, {cmd:r(gamma)}; {cmd:statement} stores {cmd:r(statement)}.

{pstd}Difference-in-differences designs add: after {cmd:dgp}, {cmd:setatt}
and {cmd:describe}, macros {cmd:r(design)}, {cmd:r(trends)},
{cmd:r(analysis)}, {cmd:r(subgroup)}, scalars {cmd:r(p_treated)},
{cmd:r(p_post)}, {cmd:r(rho)}, {cmd:r(n_pre)}, {cmd:r(n_post)},
{cmd:r(reading_linear)}, {cmd:r(reading_index)} (the two readings of the
rates; differences between the subgroups in a subgroup design),
{cmd:r(att_ref)}, {cmd:r(att_sub)}, and matrices {cmd:r(rates)} (one row per
subgroup: control before, control after, treated before, treated after
without treatment, treated after) and {cmd:r(didcoef)} (the index
coefficients a, b, l, th, de); {cmd:didpilot} adds {cmd:r(pilot)} (rows,
units, periods, pairs, and the estimates). {cmd:trueatt} stores {cmd:r(att)}
(for a subgroup design the difference, with {cmd:r(att_ref)} and
{cmd:r(att_sub)}). {cmd:power}, {cmd:n} and {cmd:mde} add {cmd:r(analysis)},
{cmd:r(trends)} and {cmd:r(valid)} (1 when the analysis assumes the scale the
world has), and {cmd:power} adds {cmd:r(est_mean)} and {cmd:r(est_mcse)}, the
mean of the usable estimates and its Monte Carlo standard error;
{cmd:r(results)} of {cmd:curve} gains a column {cmd:est_mean}. When
{cmd:robust} varies {cmd:trends}, {cmd:r(scenarios)} codes it 1 = index,
2 = linear and gains the columns {cmd:est_mean} (power mode) and {cmd:valid};
scenarios with {cmd:valid} = 0 are shown but left out of the range, the worst
case and n_max.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. powerape dgp d1, focal(treat, binary p(0.5)) covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.15)}{p_end}
{phang2}{cmd:. powerape setape d1, target(0.05)}{p_end}
{phang2}{cmd:. powerape power d1, n(2000) claim(minimum) sesoi(0.02) nsim(1000) seed(1)}{p_end}
{phang2}{cmd:. powerape statement}{p_end}
{phang2}{cmd:. powerape n d1, power(0.8) claim(detect) seed(2) engine(mata)}{p_end}
{phang2}{cmd:. powerape mde d1, n(1500) claim(detect) seed(3) engine(mata)}{p_end}
{phang2}{cmd:. powerape robust d1, n(1500) claim(detect) vary(baseline(0.2 0.4) signal(0 0.3)) seed(4) engine(mata)}{p_end}

{phang2}{cmd:. powerape dgp d2, focal(treat, binary p(0.5)) moderator(female, binary p(0.55)) baseline(0.30)}{p_end}
{phang2}{cmd:. powerape setaie d2, target(0.04) mainfocal(0.10) mainmoderator(0.05)}{p_end}
{phang2}{cmd:. powerape power d2, n(4000) claim(detect) effects(margins)}{p_end}

{phang2}{cmd:. probit y treat age female}{p_end}
{phang2}{cmd:. powerape dgp f1, fromfit focal(treat)}{p_end}
{phang2}{cmd:. powerape setape f1, target(0.10)}{p_end}

{phang2}{cmd:. powerape dgp p1, panel focal(treat, binary p(0.5) icc(1)) covariates((size, normal icc(0.6))) nperiods(4) rho(0.3) creshare(0.5) baseline(0.30) signal(0.10)}{p_end}
{phang2}{cmd:. powerape setape p1, target(0.08)}{p_end}
{phang2}{cmd:. powerape power p1, n(150) claim(detect)}{p_end}

{phang2}{cmd:. powerape dgp i1, iv focal(training, binary p(0.4)) covariates((size, normal)) instruments((subsidy, binary p(0.5))) endogeneity(0.4) ivstrength(0.2) baseline(0.30) signal(0.10)}{p_end}
{phang2}{cmd:. powerape setape i1, target(0.10)}{p_end}
{phang2}{cmd:. powerape power i1, n(1500) claim(detect)}{p_end}

{pstd}Difference-in-differences: a logit panel, the linear analysis, a subgroup, a sweep over the scale of the trends{p_end}
{phang2}{cmd:. powerape dgp dd, did model(logit) baseline(0.20) gap(0.20) trend(0.10) rho(0.4)}{p_end}
{phang2}{cmd:. powerape setatt dd, target(0.07)}{p_end}
{phang2}{cmd:. powerape power dd, n(1500) claim(detect)}{p_end}
{phang2}{cmd:. powerape n dd, power(0.8) claim(detect) engine(mata)}{p_end}
{phang2}{cmd:. powerape robust dd, n(1500) claim(detect) vary(trends(index linear)) engine(mata)}{p_end}
{phang2}{cmd:. powerape dgp dl, did baseline(0.20) gap(0.20) trend(0.10) rho(0.5) trends(linear) analysis(linear) npre(2) npost(2)}{p_end}
{phang2}{cmd:. powerape setatt dl, target(0.10)}{p_end}
{phang2}{cmd:. powerape power dl, n(500) claim(detect)}{p_end}
{phang2}{cmd:. powerape dgp ds, did model(logit) baseline(0.30 0.20) gap(0.10 0.15) trend(0.10 0.05) rho(0.5) subgroup(female, p(0.4))}{p_end}
{phang2}{cmd:. powerape setatt ds, target(0.10 0.03)}{p_end}
{phang2}{cmd:. powerape mde ds, n(4000) claim(detect) engine(mata)}{p_end}
{phang2}{cmd:. powerape dgp dp, didpilot outcome(y) group(treated) id(firm) period(year) trend(0.05) model(logit)}{p_end}


{marker install}{...}
{title:Installation and requirements}

{pstd}
Install from the distribution directory (or a URL that serves it) with
{cmd:net install powerape, from(}{it:path}{cmd:)}; then run
{cmd:mata: mata mlib index} once, so that the Mata library
{cmd:lpowerape.mlib} is found, and {cmd:powerape version} to check that the
ado-file and the library match. {cmd:net get powerape, from(}{it:path}{cmd:)}
copies {cmd:build.do} (the {cmd:.mata} sources are installed with the
package). In this project the alternative is
{cmd:adopath ++ "}{it:...}{cmd:/powerape-stata/ado"} followed by
{cmd:mata: mata mlib index}.

{pstd}
Requirements: Stata 16 or later (frames). The default {cmd:engine(stata)} on
the IV route calls {cmd:cfprobit} (Stata 18 or later); {cmd:engine(mata)}
covers every route on Stata 16 or later; {cmd:effects(margins)} for the AIE
calls the community-contributed {cmd:ginteff} (Radean). The package was
tested on StataNow 19.5 (update 14 Oct 2026) only. The shipped Mata library
was compiled by StataNow 19.5; a library compiled by a newer Stata may not
load in an older one, so on Stata 16-18 rebuild it by running {cmd:build.do}
(it finds the installed sources on the adopath and replaces
{cmd:lpowerape.mlib} in place) and then {cmd:mata: mata mlib index}.


{marker validation}{...}
{title:Validation}

{pstd}
The port is validated against the R package in four tiers (see
{bf:powerape-stata/validation/RESULTS.md}): exact finite-sample power by
binomial enumeration; fixed-data parity of every estimator with R's internals
and with {cmd:probit}/{cmd:logit} + {cmd:margins}, {cmd:ginteff} and
{cmd:cfprobit}; calibration and inversion constants; and stochastic agreement
of power, required n, MDE, robustness sweeps and justified alpha on matched
designs. Version 0.1.1 adds the {cmd:effects(margins)} IV AIE and the Stata
engine's {cmd:separation(keep)} to tiers A and B. Version 0.1.2 (also
2026-10-16) adds the {cmd:fit:} line and {cmd:r(fit)}, {cmd:r(engine)},
{cmd:r(effects)}, a print-only change that leaves every validated number
unchanged, and includes the fix by which a per-command {opt engine()} or
{opt effects()} no longer alters the session default. Version 0.1.3
(2026-10-17) builds the unit index of a panel draw in one step instead of a
loop, which makes panel calibration about as fast as in R (seconds instead of
about a minute at the default {opt nint()}); it changes no number: the panel
constants of tiers C and D are bit-identical to 0.1.2.

{pstd}
Version 0.2.0 (2026-11-13) adds the difference-in-differences route, which has
its own four tiers, against version 1.12.0 of the R package: exact power
by enumeration of the four group-by-period cells under both engines and both
analyses; fixed-data parity with R's internals and with {cmd:logit}/{cmd:probit}
+ {cmd:margins}, {cmd:regress}, {cmd:didregress} and {cmd:xtdidregress};
calibration constants; and stochastic agreement of power, rejection rates,
required n, MDE, sweeps and justified alpha. The results of the routes that
existed before are unchanged. Not ported from version 1.12.0 of the R package:
the three-way interaction design, and the linear probability analysis of the
designs other than difference-in-differences.


{title:References}

{phang}Caron, L. 2025. Triple difference designs with heterogeneous treatment
effects. arXiv:2502.19620.{p_end}
{phang}Hancock, G. R., and Feng, Y. 2025. Power analysis under uncertainty
about contextual parameters. {it:Psychological Methods}.{p_end}
{phang}Maier, M., and Lakens, D. 2022. Justify your alpha: A primer on two
practical approaches. {it:Advances in Methods and Practices in Psychological Science} 5(2).{p_end}
{phang}Puhani, P. A. 2012. The treatment effect, the cross difference, and the
interaction term in nonlinear "difference-in-differences" models.
{it:Economics Letters} 115(1): 85-87.{p_end}
{phang}Riesthuis, P. 2024. Simulation-based power analyses for the smallest
effect size of interest. {it:Advances in Methods and Practices in Psychological Science} 7.{p_end}
{phang}Wooldridge, J. M. 2023. Simple approaches to nonlinear
difference-in-differences with panel data. {it:The Econometrics Journal}
26(3): C31-C66.{p_end}
{phang}Wooldridge, J. M. 2019. Correlated random effects models with
unbalanced panels. {it:Journal of Econometrics} 211: 137-150.{p_end}
{phang}Wulff, J. N. powerape: Power analysis for average partial effects. R
package version 1.12.0. https://github.com/jespernwulff/powerape{p_end}


{title:Author}

{pstd}Stata port of Jesper Wulff's R package powerape; Aarhus University.{p_end}
