{smcl}
{* *! version 0.95 12Aug2017}{...}
{viewerjumpto "Syntax" "rdbounds##syntax"}{...}
{viewerjumpto "Description" "rdbounds##description"}{...}
{viewerjumpto "Options" "rdbounds##options"}{...}
{viewerjumpto "Examples" "rdbounds##examples"}{...}
{viewerjumpto "Saved results" "rdbounds##saved_results"}{...}

{title:Title}

{p 4 8}{cmd:rdbounds} {hline 2} Manipulation Robust Regression Discontinuity Bounds Estimation.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdbounds } {it:depvar} {it:runningvar} {ifin} 
[{cmd:,} 
{cmd:covs(}{it:covsvars}{cmd:)}
{cmd:treatment(}{it:treatmentvar}{cmd:)}
{cmd:c(}{it:#}{cmd:)}
{cmd:discretex}
{cmd:discretey}
]
{cmd:bwsx(}{it:#}{cmd:)}
[
{cmd:bwy(}{it:#}{cmd:)}
{cmd:kernel(}{it:kernelname}{cmd:)}
{cmd:orders(}{it:#}{cmd:)}
{cmd:evaluation_ys(}{it:yvalues}{cmd:)}
{cmd:ymin(}{it:#}{cmd:)}
{cmd:ymax(}{it:#}{cmd:)}
{cmd:type(}{it:typename}{cmd:)}
{cmd:percentiles(}{it:percentiles}{cmd:)}
{cmd:num_tau_pairs(}{it:#}{cmd:)}
{cmd:num_bootstraps(}{it:#}{cmd:)}
{cmd:kn(}{it:#}{cmd:)}
{cmd:alpha(}{it:#}{cmd:)}
{cmd:potential_taus(}{it:taulist}{cmd:)}
{cmd:righteffects}
{cmd:yextremes(}{it:# #}{cmd:)}
{cmd:num_lambdas(}{it:#}{cmd:)}
{cmd:num_rs(}{it:#}{cmd:)}
{cmd:kernel_y(}{it:kernelname}{cmd:)}
{cmd:outputfile(}{it:stringname}{cmd:)}
{cmd:outputfileInputs(}{it:stringname}{cmd:)}
{cmd:replace}
{cmd:verbose}
]
{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}This function implements the estimation procedure developed in {browse "http://www.nber.org/papers/w22892":Gerard, Rokkanen, and Rothe (2018)}, to estimate bounds on treatment effects under potential manipulation of the running varible.{p_end}

{title:Dependencies}

{p 4 8}This function requires the package {it:moremata}, which is used to speed up matrix computations. It can be installed via {cmd: ssc install moremata}.{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{it:depvar} specifies the dependent/outcome variable for the RDD.{p_end}

{p 4 8}{it:runningvar} specifies the running variable for the RDD.{p_end}

{p 4 8}{cmd:covs(}{it:covsvars}{cmd:)} specifies covariate(s) to implement the covariate-based refinement.{p_end}

{p 4 8}{cmd:treatment(}{it:treatmentvar}{cmd:)} specifies the treatment status variable if implementing a Fuzzy RDD. Defaults to computation of Sharp RDD results only.{p_end}

{p 4 8}{cmd:c(}{it:#}{cmd:)} specifies the threshold for assignment to treatment (assigned iff {it: runningvar} >= {it:#}}). Defaults to 0.{p_end}

{p 4 8}{cmd:discretex} if set, treat each value of x as a mass-point for density estimation.{p_end}

{p 4 8}{cmd:discretey} if set, treat each value of y as a mass-point for density estimation.{p_end}

{p 4 8}{cmd:bwsx(}{it:#}{cmd:)} is a {it: numlist} of bandwidths in x, respectively for 1) estimation of the discontinuity in the density of x at the cutoff; and 2)local polynomial estimation of conditional means. Expects either a single bandwidth to be used for both or a {it: numlist} of two. Required.{p_end}

{p 4 8}{cmd:bwy(}{it:#}{cmd:)} is a bandwidth for density estimation of y, implemented if {cmd:discretey} not set. Required if {cmd:discretey} not set.{p_end}

{p 4 8}{cmd:kernel(}{it:kernelname}{cmd:)} specificies a kernel function to be used throughout estimation for {it:runningvar}. Choices are {cmd: triangular}, {cmd: rectangular}, {cmd: gaussian} and {cmd: epanechnikov}. Defaults to {cmd: triangular}. {p_end}

{p 4 8}{cmd:orders(}{it:#}{cmd:)} specifies the order of polynomial regression, for: 1) estimation of the discontinuity in the density at the cutoff (tau in paper), and 2) local polynomial regressions. Expects either a single integer to be used for both or a vector of two values. Defaults to 1 (local linear regression) for all. Estimation of tau can only be implemented up to quadratic order if {cmd: discretex} not set.{p_end}

{p 4 8}{cmd:evaluation_ys(}{it:yvalues}{cmd:)} an explicit {it: numlist} of values of {it: depvar} (i.e. "y") to evaluate CDF's at (and PDF's if {cmd: discretey} not set). Values must be in ascending order with no repeats. If {cmd: evaluation_ys()} is not set, the set of unique values of y in the sample will be used. Caution is required if {cmd: discretey} is set, because computation will assume a probability mass function can be estimated from differences in estimated CDF's at subsequent values of {cmd: evaluation_ys()}. This can bias FRD estimates if {cmd: evaluation_ys()} does not contain all values in the support of dependent variable. {cmd: evaluation_ys()} can take a {it: numlist} in the form e.g.{cmd: evaluation_ys("0 0.1 to 10")}, for instance to create a list of values from 0 to 10 by increments of 0.1. {p_end}

{p 4 8}{cmd:ymin(}{it:#}{cmd:)} left/lower bound on {it: depvar} at which to implement a boundary kernel correction if {cmd: discretey} is not set and {it: depvar} is a variable with bounded support (e.g. after censoring). If omitted, no boundary kernel correction is implemented on the left side of the support of {it: depvar}.{p_end}

{p 4 8}{cmd:ymax(}{it:#}{cmd:)} right/upper bound on {it: depvar} at which to implement a boundary kernel correction if {cmd: discretey} is not set and {it: depvar} is a variable with bounded support (e.g. after censoring). If omitted, no boundary kernel correction is implemented on the right side of the support of {it: depvar}.{p_end}

{p 4 8}{cmd:type(}{it:typename}{cmd:)} {it:typename} = {cmd: ate} for average treatment effects (default) or {cmd: qte} for quantile treatment effects at the percentiles given by parameter {cmd: percentiles}. Defaults to {cmd: ate}.{p_end}

{p 4 8}{cmd:percentiles(}{it:percentiles}{cmd:)} {it: numlist} of percentiles at which to asses quantile treatment effects. Defaults to median (50). User may add -1 as a percentile, in order to estimate average treatment effects along with QTE's. Must be integers between 1 and 99 (aside from -1 if included, which must come first). For example, {cmd: percentiles(-1,30,50)} will compute ATEs as well as the 30 percent and 50 percent QTEs.{p_end}

{p 4 8}{cmd:num_tau_pairs(}{it:#}{cmd:)} integer number of points to search over in the set of possible values for (tau0, tau1) in notation of paper, for fuzzy RD estimation. Defaults to 50. If set to 1, the single tau is set to the "rightmost" (t=1) extreme of the set T (Section 3.1 of paper), such that user can enforce the assumption that always-assigned units always receive treatment ("Refinement B", or Theorem 4 in paper), if this is consistent with data. If only sharp estimands are desired (i.e. {cmd:treatment} is not set) it is most efficient to set this to zero.{p_end}

{p 4 8}{cmd:refinementA} if set, additionally calculate refined bounds with the restriction that always assigned units are at least as likely to be treated as potentially assigned units (i.e. tau1 ge tau; see Corollary 1 in paper){p_end}

{p 4 8}{cmd:refinementB} if set, additionally calculate refined bounds with the restriction that always assigned units are always treated (i.e. tau0 = 0; see Corollary 2 in paper){p_end}

{p 4 8}{cmd:num_bootstraps(}{it:#}{cmd:)} number of bootstrap samples for estimating confidence intervals (and for diagnostic testing of the estimated discontinuity in the density at the cutoff). To avoid bootstrap testing altogether, set {cmd: num_bootstraps(0)}.{p_end}

{p 4 8}{cmd:kn(}{it:#}{cmd:)} a hardcoded constant for kappa_n (see Section 5.2 on inference in paper). Defaults to log(n)^{1/2}, where n is the number of observations{p_end}

{p 4 8}{cmd:alpha(}{it:#}{cmd:)} sets the level for confidence intervals. Defaults to alpha=.05 for 95 percent confidence intervals{p_end}

{p 4 8}{cmd:potential_taus(}{it:taulist}{cmd:)} {it: numlist} of different values of tau to use for the confidence intervals estimating the potential impact of manipulation, e.g. {cmd: potential_taus(.025 .05 .1 .2)}.{p_end}

{p 4 8}{cmd:righteffects} if set, additionally estimate causal effects for units just to the right of the cutoff (Section 4.3 of paper).{p_end}

{p 4 8}{cmd:yextremes(}{it:# #}{cmd:)} extreme values (lower and upper) to assume if {cmd: righteffects} is set, e.g. {cmd: yextremes(0 10)}. Defaults to the sample range of{p_end}

{p 4 8}{cmd:num_lambdas(}{it:#}{cmd:)} integer number of points to search over for the causal effect of units just to the right of the cutoff (lambda in paper). Defaults to 50.{p_end}

{p 4 8}{cmd:num_rs(}{it:#}{cmd:)} integer number of points to search over to find r_alpha when computing confidence intervals (see Section 5.2 of paper). {p_end}

{p 4 8}{cmd:kernel_y(}{it:kernelname}{cmd:)} allows a separate kernel for density estimation of {it: depvar}. Same choices as kernel for {it: runningvar}. Defaults to kernel specified for use with {it: runningvar}. {p_end}

{p 4 8}{cmd:bwsxcov(}{it:#}{cmd:)} an optional separate {cmd:bwsx} to use for quantities that are computed on a subsample conditioned on a value of {cmd:covs} (e.g. covariate-conditional CDFs).

{p 4 8}{cmd:bwycov(}{it:#}{cmd:)} an optional separate {cmd: bwy} to use for quantities that are computed on a subsample conditioned on a value of {cmd:covs} (e.g. covariate-conditional CDFs).

{p 4 8}{cmd:outputfile(}{it:filename}{cmd:)} saves as {it:filename} a Stata dataset containing intermediate estimands, in particular averages (and/or quantiles) of potential outcome distributions corresponding to all of the estimands. {p_end}

{p 4 8}{cmd:outputfileInputs(}{it:filename}{cmd:)} saves as {it:filename} a Stata dataset containing fundamental estimated inputs (densities, CDFs, and basic quantities like G and s(y)(1-tau0)) (see paper for notation) for the main sample (not for bootstrap samples). {p_end}

{p 4 8}{cmd:replace} allows {cmd: outputfile} and {cmd: outputfileInputs} to overwrite an existing file. {p_end}

{p 4 8}{cmd:verbose} allows more detailed reporting of warnings in estimation. {p_end}

    {hline}


{marker examples}{...}
{title:Example: simulated sample data}

    
{p 4 8}Generate sample data{p_end}
{p 8 8}{cmd:. rdbounds_sampledata, samplesize(50000) covs clear}{p_end}

{p 4 8}Run estimation{p_end}
{p 8 8}{cmd:. rdbounds y x, c(0) covs(cov) righteffects bwsx(.2,1) bwy(.1) evaluation_ys("0 0.1 to 23") orders(1) yextremes(0 23) ymin(0) ymax(23) kernel(triangular) num_bootstraps(0)	}{p_end}

{p 4 8}Display results{p_end}
{p 8 8}{cmd:. disp "tau_hat: `e(tau_hat)'"}{p_end}
{p 8 8}{cmd:. matrix list e(treatment_effects_ATE)}{p_end}


{marker saved_results}{...}
{title:Saved results}

{p 4 8}{cmd:rdbounds} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(tau_hat)}}estimated discontinuity in the density of running variable at cutoff{p_end}
{synopt:{cmd:e(takeup_increase)}}estimated increase in probability of treatment at the cutoff (if fuzzy RDD){p_end}
{synopt:{cmd:e(seconds_taken)}}length of time taken for computation in seoonds{p_end}

{synopt:{cmd:e(tau_hat_CI_lower)}} lower limit of confidence interavl for {cmd: tau_hat}{p_end}
{synopt:{cmd:e(tau_hat_CI_upper)}} upper limit of confidence interavl for {cmd: tau_hat}{p_end}
{synopt:{cmd:e(takeup_CI_lower)}} lower limit of confidence interavl for {cmd: takeup}{p_end}
{synopt:{cmd:e(takeup_CI_upper)}} upper limit of confidence interavl for {cmd: takeup}{p_end}

{synopt:{cmd:e(treatment_effects_{it:effect)}} table of treatment effect estimates and confidence intervals for treatment effect {it:effect}, e.g. ATE, QTE10, QTE50, etc.}{p_end}
{synopt:{cmd:e(fixedtau_{it:effect)}} table of treatment effect estimates and confidence intervals for treatment effect {it:effect} and fixed values of tau specified through the {cmd:potential_taus} option}{p_end}

{title:References}

{p 4 8}F. Gerard, M. Rokkanen, and C. Rothe. 2016. Bounds on Treatment Effects in Regression Discontinuity Designs under Manipulation of the Running Variable, with an Application to Unemployment Insurance in Brazil.
Working Paper.
{browse "http://www.nber.org/papers/w22892"}.{p_end}

{title:Authors}

{p 4 8}Francois Gerard, Columbia University.
{browse "mailto:fgerard@columbia.edu":fgerard@columbia.edu}.{p_end}

{p 4 8}Leonard Goff, Columbia University.
{browse "mailto:leonard.goff@columbia.edu":leonard.goff@columbia.edu}.{p_end}

{p 4 8}Miikka Rokkanen, Amazon.
{browse "mailto:mr3454@columbia.edu":mr3454@columbia.edu}.{p_end}

{p 4 8}Christoph Rothe, University of Mannheim.
{browse "mailto:christoph.s.rothe@gmail.com ":christoph.s.rothe@gmail.com}.{p_end}

{title:See Also}

{p 4 8}{help rdbounds_sampledata:rdbounds_sampledata}{p_end}
