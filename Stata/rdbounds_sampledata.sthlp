{smcl}
{* *! version 0.95 12Aug2017}{...}
{viewerjumpto "Syntax" "rdbounds##syntax"}{...}
{viewerjumpto "Description" "rdbounds##description"}{...}
{viewerjumpto "Options" "rdbounds##options"}{...}
{viewerjumpto "Examples" "rdbounds##examples"}{...}
{viewerjumpto "Saved results" "rdbounds##saved_results"}{...}

{title:Title}

{p 4 8}{cmd:rdbounds_sampledata} {hline 2} Generate a simulated dataset for use with the {help rdbounds:rdbounds} command {p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdbounds_sampledata }
{cmd:,}[
{cmd:samplesize(}{it:#}{cmd:)}
{cmd:covs}
]
{cmd:clear}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8} This function replaces the current dataset in memory with a simulated sample dataset for testing of {help rdbounds:rdbounds}.{p_end}

{p 4 8}The x-values of potentially-assigned units (95% of sample) are normally distributed (variance of 5, censored at -10 and 10) around zero (the cutoff) and always-assigned units (5% of sample) follow a triangular distribution from the cutoff (x=0) to x=5. A random 5% of all units are never-takers, another 25% are always-takers, and the remaining units are compliers. Outcome values are then generated with a treatment effect of 2 for potentially-assigned units, 5 for always-assigned units, and following a linear trend with some normally distributed noise.
Specifically: y=(x+10)/2*treatment*(alwaysassigned=0)+5*treatment*(alwaysassigned=1)+normal(0,1)) and y is censored at 0 and 23.
{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{cmd:samplesize({it: #})} number {it: #} of observations to include in the dataset. Defaults to 30,000.{p_end}
{p 4 8}{cmd:covs} if set, a covariate will be included in the dataset.{p_end}
{p 4 8}{cmd:covs} required to allow {cmd:rdbounds_sampledata} to replace the current data in memory. {p_end}

    {hline}


{marker examples}{...}
{title:Example: simulated sample data}

    
{p 4 8}Generate sample data of 50,000 observations{p_end}
{p 8 8}{cmd:. rdbounds_sampledata, samplesize(50000) covs clear}{p_end}

{title:References}

{p 4 8}F. Gerard, M. Rokkanen, and C. Rothe. 2016. Bounds on Treatment Effects in Regression Discontinuity Designs under Manipulation of the Running Variable, with an Application to Unemployment Insurance in Brazil. NBER Working Paper 22892. {browse "http://www.nber.org/papers/w22892"}.{p_end}
{title:Authors}

{p 4 8}Francois Gerard, Columbia University.
{browse "mailto:fgerard@columbia.edu":fgerard@columbia.edu}.{p_end}

{p 4 8}Leonard Goff, Columbia University.
{browse "mailto:leonard.goff@columbia.edu":leonard.goff@columbia.edu}.{p_end}

{p 4 8}Miikka Rokkanen, Amazon.
{browse "mailtomiikka.rokkanen@gmail.com":mr3454@columbia.edu}.{p_end}

{p 4 8}Christoph Rothe, University of Mannheim.
{browse "mailto:rothe@vwl.uni-mannheim.de":christoph.s.rothe@gmail.com}.{p_end}

{title:See Also}

{p 4 8}{help rdbounds:rdbounds}{p_end}

