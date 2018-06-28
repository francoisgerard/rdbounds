## Manipulation Robust Regression Discontinuity Bounds Estimation in Stata and R

This is a public repository for the package ```rdbounds``` for Stata and R, implementing the estimation procedure developed in the paper: [Bounds on Treatment Effects in Regression Discontinuity Designs under Manipulation of the Running Variable, with an Application to Unemployment Insurance in Brazil](http://www.nber.org/papers/w22892 "NBER Working Paper"), by François Gerard, Miikka Rokkanen, and Christoph Rothe.

This is a preliminary version of the code and is still undergoing testing. We appreciate any feedback or issues noted. Please direct your comments to leonard.goff at columbia dot edu. The current version is 1.00, dated July 1st, 2018.

The Stata version of the code generally runs a bit faster, but this may depend on your version of Stata and the number of cores on your computer, among other things.

## Stata Version

The file [Stata/rdbounds.ado](Stata/rdbounds.ado) is the main code file. You can install ```rdbounds``` by downloading this to your ado directory for Stata, e.g. C:\ado\personal.

The files [Stata/rdbounds.sthlp](Stata/rdbounds_sampledata.sthlp) and [Stata/rdbounds_sampledata.sthlp](Stata/rdbounds.sthlp) are Stata help files for the main ``rdbounds''' function and for the function ``rdbounds_sampledata''' (which generates a sample dataset to experiment with) respectively. These should also be dropped into your local ado directory or can be viewed directly in the Stata help file viewer.
 
Alternatively, the package can be installed directly from Stata by running
```
net describe rdbounds from https://github.com/leonardgoff/rdbounds/Stata
net install rdbounds
```
 
The ```rdbounds``` function also requires the Stata package ```moremata```, which can be installed by running:
```
ssc install moremata
```

## R Version

See [R/rdbounds.pdf](R/rdbounds.pdf) for documentation. The code is viewable in [R/rdbounds.R](R/rdbounds.R).

You may install the ```rdbounds``` package directly from this repository using the ```devtools``` package, by running the following commands in R (the first three lines may not be necessary if you already have the corresponding packages installed):

```{r}
install.packages("formattable")
install.packages("data.table")
install.packages("devtools")
library(devtools)
install_github("leonardgoff/rdbounds/R")
```

Alternatively, you may download [R/rdbounds_1.0.tar.gz](R/rdbounds_1.0.tar.gz) and install the package from source code. You will need to also install the packages ```formattable``` and ```data.table``` if you do not have them already.

Here's some test code to get you going, once installed:

```{r}
library(formattable)
library(data.table)
library(rdbounds)
df<-rdbounds_sampledata(30000, covs=TRUE)
rdbounds_est<-rdbounds(y=df$y,x=df$x, covs=df$cov, treatment=df$treatment, c=0,
                       discrete_x=FALSE, discrete_y=FALSE, bwsx=c(.1,1), bwy = .1,
                       kernel="triangular", orders=c(1,1),
                       evaluation_ys = seq(from = 0, to=23, by=.1), ymin=0, ymax=23,
                       right_effects=TRUE, yextremes = c(0,23),
                       num_bootstraps=c(1,1))
rdbounds_summary(rdbounds_est, title_prefix="Sample Data Results")
```