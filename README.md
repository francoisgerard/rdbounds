## Manipulation Robust Regression Discontinuity Bounds Estimation in Stata and R

This is a public repository for the package ```rdbounds``` for Stata and R, which implements the estimation procedure developed in the paper [Bounds on Treatment Effects in Regression Discontinuity Designs under Manipulation of the Running Variable, with an Application to Unemployment Insurance in Brazil](http://www.nber.org/papers/w22892 "NBER Working Paper"), by François Gerard, Miikka Rokkanen, and Christoph Rothe.

This is a preliminary version of the code and is offered without warranty. We appreciate any feedback or issues noted. Please direct your comments to leonard.goff at columbia dot edu. The current version is 1.00, dated July 1st, 2018.

The Stata version of the code generally runs a bit faster than the R version with many bootstrap resamples, but this may depend on your version of Stata and the number of cores on your computer, among other things.

## Stata Version

The file [Stata/rdbounds.ado](Stata/rdbounds.ado) is the main code file, and [Stata/rdbounds.sthlp](Stata/rdbounds_sampledata.sthlp) is a Stata help files for the ```rdbounds``` function. The second .ado file [Stata/rdbounds_sampledata.ado](Stata/rdbounds_sampledata.ado) generates a sample dataset to experiment with, and its associated help file is  [Stata/rdbounds_sampledata.sthlp](Stata/rdbounds.sthlp).

The package can be installed directly from within Stata by running
```
net from https://raw.githubusercontent.com/leonardgoff/rdbounds/master/Stata/
net describe rdbounds
net install rdbounds
```

Alternatively, you can install the package by downloading ```rdbounds.ado``` and ```rdbounds_sampledata.ado``` to your ado directory for Stata, e.g. C:\ado\personal. The help files should also be dropped into your local ado directory or can be viewed directly in the Stata help file viewer. 

The ```rdbounds``` function requires the Stata package ```moremata```, which can be installed by running:
```
ssc install moremata
```

Here's some test code to get you going, once ```rdbounds``` is installed:

```
set seed 1
rdbounds_sampledata, samplesize(50000) clear covs
rdbounds y x, c(0) treatment(treatment) covs(cov) bwsx(.2, .5) bwy(.1) evaluation_ys("0 .2 to 23") orders(1) kernel(epanechnikov) type(ate) refinementA refinementB righteffects yextremes(0 23) num_bootstraps(0)
disp "tau_hat: `e(tau_hat)'"
disp "takeup increase: `e(takeup_increase)'"
matrix list e(treatment_effects_ATE)
```

## R Version

The main code file is [R/R/rdbounds.R](R/R/rdbounds.R), and the documentation can be found here: [R/rdbounds.pdf](R/rdbounds.pdf).

You may install the package directly from this repository using the ```devtools``` package, by running the following commands in R (the first three lines may not be necessary if you already have the corresponding packages installed):

```{r}
install.packages("formattable")
install.packages("data.table")
install.packages("devtools")
library(devtools)
install_github("leonardgoff/rdbounds/R")
```

Alternatively, you may download [R/rdbounds_1.0.tar.gz](R/rdbounds_1.0.tar.gz) and install the package from source code. You will need to also install the packages ```formattable``` and ```data.table``` if you do not have them already.

Here's some test code to get you going, once ```rdbounds``` is installed:

```{r}
library(formattable)
library(data.table)
library(rdbounds)
set.seed(1)
df<-rdbounds_sampledata(50000, covs=TRUE)
rdbounds_est<-rdbounds(y=df$y,x=df$x, covs=as.factor(df$cov), treatment=df$treatment, c=0,
                       discrete_x=FALSE, discrete_y=FALSE, bwsx=c(.2,.5), bwy = .1,
                       kernel="epanechnikov", orders=1,
                       evaluation_ys = seq(from = 0, to=23, by=.2), ymin=0, ymax=23,
                       num_bootstraps=0, refinement_A=TRUE, refinement_B=TRUE, 
                       right_effects=TRUE, yextremes = c(0,23))
rdbounds_summary(rdbounds_est, title_prefix="Sample Data Results")
```