{smcl}
{* *! version 1.1.0  22oct2026  Sánchez-Aragón, Acharya, Menon, and Voia}{...}
{vieweralsosee "[ST] streg" "help streg"}{...}
{vieweralsosee "[R] ivregress" "help ivregress"}{...}
{viewerjumpto "Syntax" "aaaft##syntax"}{...}
{viewerjumpto "Description" "aaaft##description"}{...}
{viewerjumpto "Options" "aaaft##options"}{...}
{viewerjumpto "Remarks" "aaaft##remarks"}{...}
{viewerjumpto "Exporting results" "aaaft##export"}{...}
{viewerjumpto "Stored results" "aaaft##results"}{...}
{viewerjumpto "Examples" "aaaft##examples"}{...}
{viewerjumpto "References" "aaaft##references"}{...}
{viewerjumpto "Authors" "aaaft##authors"}{...}

{title:Title}

{phang}
{bf:aaaft} {hline 2} Rank-based instrumental variables inference for accelerated failure time models


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:aaaft}
{depvar}
{it:treatvar}
{it:covariates}
{ifin}
[{cmd:,}
{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opth instr:ument(varname)}}instrumental variable for {it:treatvar}{p_end}

{syntab:Grid}
{synopt:{opth gridmin(#)}}lower bound of the beta grid; default is {cmd:gridmin(0.01)}{p_end}
{synopt:{opth gridmax(#)}}upper bound of the beta grid; default is {cmd:gridmax(0.99)}{p_end}
{synopt:{opth gridstep(#)}}step of the beta grid; default is {cmd:gridstep(0.01)}{p_end}

{syntab:Simulation}
{synopt:{opth reps(#)}}number of simulation draws for the confidence sets; default is {cmd:reps(1000)}{p_end}
{synopt:{opth level(#)}}confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt seed(#)}}random-number seed for the simulations{p_end}

{syntab:Empirical quantile score}
{synopt:{opt fo(filename)}}file with draws of the error used to build an empirical quantile score{p_end}
{synopt:{opt fovar(name)}}variable in {it:filename} that holds the draws; default is {cmd:eps}{p_end}

{syntab:Instrument}
{synopt:{opt zcenter}}use the centered instrument instead of the instrument residualized on the covariates{p_end}

{syntab:Reporting}
{synopt:{opt norank}}skip the rank estimates and the rank confidence sets{p_end}
{synopt:{opt nols}}skip the least-squares IV estimate and confidence sets{p_end}

{syntab:Output}
{synopt:{opt export(filename)}}write the results to {it:filename} as a CSV file{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:depvar} is the log duration outcome. {it:treatvar} is the endogenous
treatment variable. {it:covariates} is a list of one or more exogenous
control variables. All variables must be numeric. At least one covariate is
required.{p_end}

{p 4 6 2}
{cmd:aaaft} requires the {cmd:moremata} package:
{stata ssc install moremata}. The command is made of two files,
{cmd:aaaft.ado} and {cmd:aaaft_export.ado}, which must sit in the same
folder on the ado path.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:aaaft} implements rank-based instrumental variables (IV) inference for
accelerated failure time (AFT) duration models. It follows the framework of
Acharya, Khalaf, Menon, Voia, Yazbeck, and Wensley (2026). The rank
confidence sets are exact and robust to weak instruments.

{pstd}
The AFT model specifies the log duration outcome as

{phang2}
ln(T{sub:i}) = {it:beta}*{it:treat}{sub:i} + {it:x}{sub:i}'{it:gamma} + {it:epsilon}{sub:i}

{pstd}
where {it:beta} is the scalar coefficient on the endogenous treatment,
{it:x}{sub:i} holds the exogenous covariates, and {it:epsilon}{sub:i} is an
error term with unspecified distribution.

{pstd}
With {cmd:instrument()}, {cmd:aaaft} reports four sets of results:

{phang2}
1. Hodges-Lehmann-Sen (HLS) point estimates. For each candidate {it:beta},
the command partials the covariates out of ln(T) - {it:beta}*{it:treat},
ranks these residuals, and applies a score function to the ranks. The rank
statistic G({it:beta}) is the squared inner product of the scores with the
instrument residualized on the covariates, divided by the sum of squares of
that instrument. The HLS estimate minimizes G({it:beta}) over a grid. Three
score functions are always used: logistic, Gumbel (extreme value), and
exponential. With {cmd:fo()}, a fourth empirical quantile score is added.

{phang2}
2. Rank confidence sets. The null distribution of G does not depend on
{it:beta}, so it is simulated once from random rank permutations. The set
keeps every grid value with G({it:beta}) below the simulated critical value.
These sets are exact and robust to weak instruments.

{phang2}
3. The least-squares (LS) IV point estimate, the value of {it:beta} where the
generalized Anderson-Rubin (GAR) statistic is zero. With one instrument it
equals the usual IV estimate.

{phang2}
4. LS confidence sets, obtained by inverting the GAR statistic. The critical
value is simulated with errors drawn from the distribution of each score,
so these sets are exact only if the errors follow that distribution.

{pstd}
Without {cmd:instrument()}, {cmd:aaaft} runs in rank-only mode. The scores
are projected on the treatment residualized on the covariates instead of on
the instrument. The command reports HLS point estimates with no IV correction
and no confidence sets. With {cmd:fo()}, it also reports the HLS estimate
for the empirical quantile score. This lets you compare naive rank
estimates with IV estimates under the same scores.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth instrument(varname)} specifies the instrument for {it:treatvar}. The
instrument must be numeric, exogenous, and correlated with {it:treatvar}.
Without it, {cmd:aaaft} computes rank point estimates only.

{dlgtab:Grid}

{phang}
{opth gridmin(#)} sets the lower bound of the grid for {it:beta}. It must be
strictly less than {cmd:gridmax()}. Default is {cmd:0.01}.

{phang}
{opth gridmax(#)} sets the upper bound of the grid for {it:beta}. Default is
{cmd:0.99}.

{phang}
{opth gridstep(#)} sets the step of the grid. Point estimates and rank
confidence bounds are only as precise as this step. Default is {cmd:0.01};
{cmd:0.005} or finer is recommended for published results.

{dlgtab:Simulation}

{phang}
{opth reps(#)} sets the number of simulation draws used for the critical
values of both the rank and the LS confidence sets. It must be at least 100.
Default is {cmd:1000}, which is fine for exploration. For published results
{cmd:reps(5000)} or more is recommended; see {help aaaft##mc:Monte Carlo error}.

{phang}
{opth level(#)} sets the confidence level of the rank and LS confidence sets.
Default is {cmd:95}.

{phang}
{opt seed(#)} sets the random-number seed before the simulations, so the
confidence sets can be reproduced. With the same seed, results do not depend
on the sort order of the data.

{dlgtab:Empirical quantile score}

{phang}
{opt fo(filename)} names a Stata dataset with draws of the error
{it:epsilon}, for example simulated draws of a baseline and frailty
convolution. {cmd:aaaft} builds the empirical distribution Fo of these draws
and adds the score c(u) = Fo{sup:-1}(u). The score is used for the HLS
estimate, with or without {cmd:instrument()}, and, with {cmd:instrument()},
for the rank confidence set. At
least 1,000 nonmissing draws are needed; 100,000 or more is recommended.

{phang}
{opt fovar(name)} names the variable in the {cmd:fo()} file that holds the
draws. Default is {cmd:eps}.

{dlgtab:Instrument}

{phang}
{opt zcenter} makes the rank statistic use the instrument only centered at
its mean. By default the instrument is
residualized on the covariates, as in the least-squares routine. Both
choices keep the rank confidence set exact. The residualized instrument
removes the part of the instrument explained by the covariates, which
matters when the two are correlated.

{dlgtab:Reporting}

{phang}
{opt norank} skips the HLS point estimates and the rank confidence sets.

{phang}
{opt nols} skips the LS point estimate and the LS confidence sets. It saves
time when only the rank results are needed.

{dlgtab:Output}

{phang}
{opt export(filename)} writes the results to {it:filename} as a CSV file,
one row per score; see {help aaaft##export:Exporting results}. Any existing
file is overwritten. A name without an extension gets {cmd:.csv}, and any
other extension is replaced by {cmd:.csv}.


{marker remarks}{...}
{title:Remarks}

{pstd}
{ul:Score functions}

{pstd}
The rank statistic is computed for three score functions. Each one matches
an error distribution in the AFT model.

{phang2}
Logistic: c(u) = ln(u/(1-u)), for a logistic error distribution.

{phang2}
Gumbel (extreme value): c(u) = -ln(-ln(u)), for a maximum extreme-value
(Gumbel) error distribution. Note that a Weibull duration model has a
minimum extreme-value error for ln(T), which is the reflection of this
distribution.

{phang2}
Exponential: c(u) = -ln(1-u), the quantile function of an exponential
distribution (a Savage-type score). It is not the empirical quantile score of
Acharya et al. (2026).

{phang2}
Empirical quantile, with {cmd:fo()}: c(u) = Fo{sup:-1}(u), where Fo is the
empirical distribution of the supplied draws, with linear interpolation
between order statistics.

{pstd}
Adding a constant to a score or multiplying it by a positive number does not
change the results, so only the shape of Fo matters, not its location or
scale.

{pstd}
{ul:Empirical quantile score and exactness}

{pstd}
The rank confidence set is exact for any score fixed in advance. When Fo is
calibrated with parameters estimated from the same sample, the score depends
on the data and the set is exact only conditionally on that calibration. A
split sample, with Fo calibrated on one part and {cmd:aaaft} run on the
other, keeps exactness.

{pstd}
{ul:Weak instruments}

{pstd}
The rank confidence sets keep their coverage whatever the strength of the
instrument. With a weak instrument the set can be wide, reach the edge of the
grid, or have gaps. The output marks a set that reaches the grid edge with
{cmd:*}; widen {cmd:gridmin()} or {cmd:gridmax()} in that case. It marks a set
with gaps on the grid with {cmd:+}; the reported bounds are then the smallest
and largest values kept. An empty set means that no grid value is kept.

{pstd}
{ul:Discrete treatment variables}

{pstd}
When {it:treatvar} is discrete, such as an integer severity score, the
residuals ln(T) - {it:beta}*{it:treatvar} can have ties. {cmd:aaaft} adds a
tiny deterministic jitter, 1e-12 times the observation index, to break them.
In the pediatric ICU application the point estimates and rank confidence
sets do not change when the data are sorted in different random orders.

{marker mc}{...}
{pstd}
{ul:Monte Carlo error}

{pstd}
Both sets of confidence bounds use simulated critical values, so they carry
Monte Carlo error. The critical value is the floor({it:reps} x {it:level}/100)-th
ordered simulated value. In the pediatric ICU application, with
{cmd:reps(5000)} the logistic rank bounds move at most one grid step across
seeds and the LS bounds less than 0.003. Where G is flat the bounds move
more: the upper exponential bound for PRISM moves up to six steps. A larger
{cmd:reps()} reduces this error. Report the seed and the number
of draws with any published result.

{pstd}
{ul:Grid selection}

{pstd}
If an HLS estimate or a confidence bound falls at the edge of the grid, the
grid is too narrow. A good practice is to run first with a wide grid, such as
{cmd:gridmin(0.01) gridmax(0.99)}, and then refine {cmd:gridstep()} around the
region of interest.

{pstd}
{ul:Censored durations}

{pstd}
{cmd:aaaft} treats every duration as observed and does not adjust for
censoring. It suits data where all durations are complete. With censored
durations the rank test keeps its size, but the point estimate is biased.


{marker export}{...}
{title:Exporting results}

{pstd}
{cmd:aaaft_export} writes the results of the last {cmd:aaaft} estimation to a
CSV file. It reads from {cmd:e()}, so it must be run after {cmd:aaaft}. The
{cmd:export()} option calls it automatically.

{p 8 17 2}
{cmd:aaaft_export}
{cmd:using} {it:filename}
[{cmd:,} {opt replace}]

{pstd}
{opt replace} overwrites {it:filename} if it already exists. The file has one
row per score and these columns: {cmd:score}, {cmd:hls}, {cmd:rank_lo},
{cmd:rank_hi}, {cmd:ls_point}, {cmd:ls_lo}, {cmd:ls_hi}, {cmd:depvar},
{cmd:treatvar}, {cmd:instrument}, {cmd:N}, and {cmd:level}. Cells are empty
for results that were not computed. The file can be read back with
{cmd:import delimited}.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:aaaft} stores the following in {cmd:r()}, where {it:s} is
{cmd:logistic}, {cmd:gumbel}, or {cmd:exponential}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(beta_}{it:s}{cmd:)}}HLS point estimate{p_end}
{synopt:{cmd:r(cir_l_}{it:s}{cmd:)}}lower bound of rank confidence set{p_end}
{synopt:{cmd:r(cir_u_}{it:s}{cmd:)}}upper bound of rank confidence set{p_end}
{synopt:{cmd:r(cir_edge_}{it:s}{cmd:)}}1 if rank set reaches the grid edge{p_end}
{synopt:{cmd:r(cir_gap_}{it:s}{cmd:)}}1 if rank set has gaps on the grid{p_end}
{synopt:{cmd:r(gcrit_}{it:s}{cmd:)}}simulated critical value of G{p_end}
{synopt:{cmd:r(beta_ls)}}LS IV point estimate{p_end}
{synopt:{cmd:r(ci_l_}{it:s}{cmd:)}}lower bound of LS confidence set{p_end}
{synopt:{cmd:r(ci_u_}{it:s}{cmd:)}}upper bound of LS confidence set{p_end}
{synopt:{cmd:r(fo_draws)}}number of draws used for Fo, with {cmd:fo()}{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(fo)}}name of the {cmd:fo()} file{p_end}
{p2colreset}{...}

{pstd}
With {cmd:fo()}, {it:s} can also be {cmd:empq}, the empirical quantile score,
for the HLS and rank set results. The LS results do not use {cmd:empq}.

{pstd}
The rank set results need {cmd:instrument()} and no {cmd:norank}. The LS
results need {cmd:instrument()} and no {cmd:nols}.

{pstd}
{cmd:aaaft} also stores the following in {cmd:e()}, so results work with
{cmd:estimates table}, {cmd:estimates store}, {cmd:outreg2}, and {cmd:esttab}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{synopt:{cmd:e(HLS_}{it:s}{cmd:)}}HLS point estimate{p_end}
{synopt:{cmd:e(CIR_lo_}{it:s}{cmd:)}}lower bound of rank confidence set{p_end}
{synopt:{cmd:e(CIR_hi_}{it:s}{cmd:)}}upper bound of rank confidence set{p_end}
{synopt:{cmd:e(beta_LS)}}LS IV point estimate{p_end}
{synopt:{cmd:e(CI_lo_}{it:s}{cmd:)}}lower bound of LS confidence set{p_end}
{synopt:{cmd:e(CI_hi_}{it:s}{cmd:)}}upper bound of LS confidence set{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:aaaft}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(treatvar)}}name of treatment variable{p_end}
{synopt:{cmd:e(instrument)}}name of instrument, if specified{p_end}
{synopt:{cmd:e(fo)}}name of the {cmd:fo()} file, if specified{p_end}
{synopt:{cmd:e(title)}}title of the estimation{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}point estimates: {cmd:HLS_logistic}, {cmd:HLS_gumbel}, {cmd:HLS_exponential}, {cmd:HLS_empq} (with {cmd:fo()}), and {cmd:LS_IV}{p_end}
{synopt:{cmd:e(V)}}zero matrix (see below){p_end}

{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{pstd}
{cmd:e(b)} holds only the point estimates that were computed. {cmd:e(V)} is a
zero matrix because {cmd:aaaft} reports confidence sets, not standard errors;
it is posted only because {cmd:estimates table} and similar commands require
it. Do not read it as a measure of precision.


{marker examples}{...}
{title:Examples}

{pstd}
{ul:IV estimation with reproducible confidence sets}

{phang2}{cmd:. aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, instrument(Trauma) gridmin(0.05) gridmax(0.50) gridstep(0.005) reps(5000) seed(12345)}{p_end}

{pstd}
{ul:Rank estimates only, no instrument}

{phang2}{cmd:. aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, gridmin(0.05) gridmax(0.50) gridstep(0.005)}{p_end}

{pstd}
{ul:Rank results only, faster}

{phang2}{cmd:. aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, instrument(Trauma) gridmin(0.05) gridmax(0.50) gridstep(0.005) seed(12345) nols}{p_end}

{pstd}
{ul:Compare naive and IV rank estimates}

{phang2}{cmd:. aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, gridmin(0.05) gridmax(0.50) gridstep(0.005)}{p_end}
{phang2}{cmd:. scalar beta_noiv = r(beta_logistic)}{p_end}
{phang2}{cmd:. aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, instrument(Trauma) gridmin(0.05) gridmax(0.50) gridstep(0.005) seed(12345)}{p_end}
{phang2}{cmd:. scalar beta_iv = r(beta_logistic)}{p_end}
{phang2}{cmd:. display "Endogeneity correction: " beta_iv - beta_noiv}{p_end}

{pstd}
{ul:Empirical quantile score from simulated error draws in eps_draws.dta}

{phang2}{cmd:. aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, instrument(Trauma) gridmin(0.05) gridmax(0.50) gridstep(0.005) reps(5000) seed(12345) fo("eps_draws.dta")}{p_end}

{pstd}
{ul:Tabulate point estimates}

{phang2}{cmd:. estimates store prism}{p_end}
{phang2}{cmd:. estimates table prism}{p_end}

{pstd}
{ul:Estimate and export in one step}

{phang2}{cmd:. aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, instrument(Trauma) gridmin(0.05) gridmax(0.50) gridstep(0.005) seed(12345) export("prism.csv")}{p_end}

{pstd}
{ul:Export the last estimation without rerunning it}

{phang2}{cmd:. aaaft_export using "prism.csv", replace}{p_end}
{phang2}{cmd:. import delimited using "prism.csv", clear}{p_end}


{marker references}{...}
{title:References}

{phang}
Acharya, A., L. Khalaf, K. Menon, M. Voia, M. Yazbeck, and D. Wensley. 2026.
Robust Duration Inference and Bias in Allocating Critical Hospital Resources.
{it:Canadian Journal of Economics}. Forthcoming.

{phang}
Anderson, T. W., and H. Rubin. 1949.
Estimation of the parameters of a single equation in a complete system of
stochastic equations.
{it:Annals of Mathematical Statistics} 20: 46-63.

{phang}
Hodges, J. L., and E. L. Lehmann. 1963.
Estimates of location based on rank tests.
{it:Annals of Mathematical Statistics} 34: 598-611.


{marker authors}{...}
{title:Authors}

{pstd}
Leonardo Sánchez-Aragón{break}
ESPOL Polytechnic University{break}
Department of Economics{break}
Guayaquil, Ecuador{break}
lfsanche@espol.edu.ec

{pstd}
Anand Acharya (deceased){break}
Carleton University{break}
Department of Economics{break}
Ottawa, ON, Canada

{pstd}
Kusum Menon{break}
University of Ottawa{break}
Pediatric Critical Care Research{break}
Ottawa, ON, Canada{break}
KMenon@cheo.on.ca

{pstd}
Marcel Voia{break}
Université d'Orléans and University of Bucharest{break}
Laboratoire d'Économie d'Orléans{break}
Orléans, France{break}
marcel.voia@univ-orleans.fr

{pstd}
Bug reports and suggestions are welcome.


{title:Also see}

{psee}
Manual: {bf:[ST] streg}, {bf:[R] ivregress}

{psee}
{stata "help streg":help streg},
{stata "help ivregress":help ivregress},
{stata "help stcox":help stcox}
{p_end}
