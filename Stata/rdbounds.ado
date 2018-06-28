/*
********************************************************************************
VERSION 1.00 (JULY 1, 2018)
DISCLAIMER: THIS CODE SHOULD BE CONSIDERED PRELIMINARY AND IS OFFERED WITHOUT WARRANTY.
WE APPRECIATE YOUR COMMENTS AND ANY ISSUES NOTICED, AT LEONARD.GOFF@COLUMBIA.EDU
Copyright (C) 2018 Francois Gerard, Leonard Goff, Miikka Rokkanen, & Christoph Rothe.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
********************************************************************************
*/

capture program drop rdbounds 
program define rdbounds, eclass

timer clear 51
timer on 51

set more off

qui set type float, permanently

*Preparation and checking input
********************************************************************************

	*Check syntax
	syntax varlist [if] [in],[covs(varlist fv) treatment(string) c(real 0) discretex discretey] bwsx(numlist) [bwy(numlist) kernel(string) orders(numlist) evaluation_ys(numlist ascending) ymin(numlist min=1 max=1) ymax(numlist min=1 max=1) type(string) percentiles(numlist ascending) num_tau_pairs(integer 50) refinementA refinementB num_bootstraps(numlist min = 0) kn(string) alpha(real 0.05) potential_taus(numlist) righteffects yextremes(numlist ascending) num_lambdas(real 50) num_rs(real 100) kernel_y(string) bwsxcov(numlist) bwycov(numlist) outputfile(string) outputfileInputs(string) replace verbose]
	
	capt mata mata which mm_collapse()
	if _rc {
		di as result in smcl `"The function {cmd:_mm_collapse()} from the package -moremata- is required for rdbounds; You can obtain it by clicking this link: {stata "ssc install moremata":ssc install moremata}"'
		exit 499
	}
	
	*Limit sample if requested
	marksample touse
	preserve
	tempfile tmpfilesample tmpfileinputs tmpfileoutcomes_b tmpfileTEs tmpfile tmpfile2
	qui keep if `touse'
	
	gen temporder=_n
	
	*Define main variables
	tokenize "`varlist'"
	local w : word count `varlist'
	if `w' < 2 {
		di as error  "{err}{cmd:varlist()} requires an outcome variable name followed by a single running variable name"  
		exit 102
	}
	if `w' == 2 {
		qui gen rdbounds_y=`1'
		qui gen rdbounds_dist_cut=`2'-`c'
	}
	if `w' > 2 {
		di as error  "{err}{cmd:varlist()} too many variables detected (expects one outcome variable and one running variable). To use higher orders of the running variable in local polynomial regression, set the orders() option"  
		exit 103
	}
	
	*Check treatment variable
	local w : word count `treatment'
	if `w' == 0 {
		local fuzzy = 0
		qui gen rdbounds_treatment=(rdbounds_dist_cut>=0)
		di as result "Note: treatment() variable omitted, only computing sharp RDD"
	}
	if `w' == 1 {
		local fuzzy = 1
		qui gen rdbounds_treatment=`treatment'
		qui gen tmp=abs(rdbounds_treatment-(rdbounds_dist_cut>=0))
		qui summ tmp
		if `r(max)'==0 {
			di as result "Warning: treatment variable passed by argument treatment appears to be identical to the condition (x>=c), computing Sharp RDD estimators only"
			local fuzzy = 0
		}
	}
	if `w' >= 2 {
		di as error  "{err}{cmd:treatment()} only accepts a single variable name"  
		exit 103
	}

	*Read in covariates
	local w : word count `covs'
	if `w' == 0 {
		local gotcovs = 0
		qui gen rdbounds_cov = .
		local ncovlevels=0
	}
	else if `w' >= 1{
		local gotcovs = 1
		by `covs', sort: gen nvals = _n == 1 
		qui count if nvals
		local ncovlevels=`r(N)'
		drop nvals
		qui egen rdbounds_cov = group(`covs')
		qui save "`tmpfile'", replace
				keep rdbounds_cov `covs'
				collapse `covs', by(rdbounds_cov)
				sort rdbounds_cov `covs'
				order rdbounds_cov `covs'
				mkmat rdbounds_cov `covs', matrix(wvals)
		qui use "`tmpfile'", clear
		di as result "Covariate list `covs' interpreted as a factor variable with `ncovlevels' levels"
		if `num_tau_pairs'^`ncovlevels' > 2500 {
			di as result "Warning: computational time and storage will scale as {cmd:num_tau_pairs}^{cmd:k}, where k is the number of levels of the varlist {cmd:covs}"
		}
		if `ncovlevels'==1 {
			di as error "{err}{cmd:covs()} has only one value, please drop"  
			exit
		}
	}

	*Check for missing values
	
	*Finalize dataset consisting only of y, dist_cut, treatment (and covs if set)
	keep rdbounds_* temporder
	rename (rdbounds_y rdbounds_dist_cut rdbounds_treatment rdbounds_cov) (y dist_cut treatment cov)	
	
	*Check for x or y flagged as having a discrete distribution
	if ("`discretex'"=="discretex") {
		local discretex=1
	}
	else{
		local discretex=0
	}
	
	if ("`discretey'"=="discretey") {
		local discretey=1
	}
	else{
		local discretey=0
	}
		
	*Check bandwidths, etc.
		*X bandwidth
		tokenize `bwsx'	
		local w : word count `bwsx'
		if `w' == 0 {
			di as error  "{err}{cmd:bwsx()} requires at least one number"  
			exit 122
		}
		if `w' == 1 {
			scalar bwx_tau=`bwsx'
			scalar bwx=`bwsx'
		}
		if `w' == 2 {
			scalar bwx_tau=`1'
			scalar bwx=`2'
		}
		if `w' > 2 {
			di as error  "{err}{cmd:bwsx()} accepts at most two numbers"  
			exit 123
		}
		
		*Y bandwidth
		scalar bwy=0
		
		local w : word count `bwy'
		if `w' == 0 & `discretey'==0 {
			di as error  "{err}{cmd:bwy()} required unless discretey is set"  
			exit
		}
		if `w'==1 {
			scalar bwy = `bwy'
		}
		if `w' > 1 {
			di as error  "{err}{cmd:bwy()} accepts only one number"  
			exit 123
		}
		
		*covariate conditional bandwidths
			tokenize `bwsxcov'	
			local w : word count `bwsxcov'
			if `w' == 0 {
				scalar bwx_tau_cov=scalar(bwx_tau)
				scalar bwx_cov=scalar(bwx)
			}
			if `w' == 1 {
				scalar bwx_tau_cov=`bwsxcov'
				scalar bwx_cov=`bwsxcov'
			}
			if `w' == 2 {
				scalar bwx_tau_cov=`1'
				scalar bwx_cov=`2'
			}
			if `w' > 2 {
				di as error  "{err}{cmd:bwsxcov()} accepts at most two numbers"  
				exit 123
			}
			
			local w : word count `bwycov'
			if `w' == 0 & `discretey'==0 {
				scalar bwy_cov=scalar(bwy)
			}
			if `w'==1 {
				scalar bwy_cov = `bwycov'
			}
			if `w' > 1 {
				di as error  "{err}{cmd:bwycov()} accepts only one number"  
				exit 123
			}
			
		
	*Kernel
	local kernel   = lower("`kernel'")
	local w : word count `kernel'
	if `w' == 0 {
		local kernel="triangular"
	}
	
	if ("`kernel'"=="triangular"){
		qui gen W = (1/scalar(bwx))*(1-abs(dist_cut/scalar(bwx)))*(abs(dist_cut/scalar(bwx))<=1)
		qui gen Wtau = (1/scalar(bwx_tau))*(1-abs(dist_cut/scalar(bwx_tau)))*(abs(dist_cut/scalar(bwx_tau))<=1)
		if `gotcovs'==1{
			qui gen W_cov = (1/scalar(bwx_cov))*(1-abs(dist_cut/scalar(bwx_cov)))*(abs(dist_cut/scalar(bwx_cov))<=1)
			qui gen Wtau_cov = (1/scalar(bwx_tau_cov))*(1-abs(dist_cut/scalar(bwx_tau_cov)))*(abs(dist_cut/scalar(bwx_tau_cov))<=1)
		}
		matrix kernelconstants = [1/2, 1/6, 1/12, 1/20, 1/30 \ 1/2, -1/6, 1/12, -1/20, 1/30]
	}
	else if("`kernel'"=="rectangular") {
		qui gen W = (1/scalar(bwx))*(1/2)*(abs(dist_cut/scalar(bwx))<=1)
		qui gen Wtau = (1/scalar(bwx_tau))*(1/2)*(abs(dist_cut/scalar(bwx_tau))<=1)
		if `gotcovs'==1{
			qui gen W_cov = (1/scalar(bwx_cov))*(1/2)*(abs(dist_cut/scalar(bwx_cov))<=1)
			qui gen Wtau_cov = (1/scalar(bwx_tau_cov))*(1/2)*(abs(dist_cut/scalar(bwx_tau_cov))<=1)
		}
		matrix kernelconstants = [1/2, 1/4, 1/6, 1/8, 1/10 \ 1/2, -1/4, 1/6, -1/8, 1/10]
	}
	else if("`kernel'"=="gaussian") {
		qui gen W = (1/scalar(bwx))*(1/sqrt(2*c(pi)))*exp(-(dist_cut/scalar(bwx))^2/2)
		qui gen Wtau = (1/scalar(bwx_tau))*(1/sqrt(2*c(pi)))*exp(-(dist_cut/scalar(bwx_tau))^2/2)
		if `gotcovs'==1{
			qui gen W_cov = (1/scalar(bwx_cov))*(1/sqrt(2*c(pi)))*exp(-(dist_cut/scalar(bwx_cov))^2/2)
			qui gen Wtau_cov = (1/scalar(bwx_tau_cov))*(1/sqrt(2*c(pi)))*exp(-(dist_cut/scalar(bwx_tau_cov))^2/2)
		}
		matrix kernelconstants = [1/2, 1/sqrt(2*c(pi)), 1/2, sqrt(2/c(pi)), 3/2 \ 1/2, -1/sqrt(2*c(pi)), 1/2, -sqrt(2/c(pi)), 3/2]
	}
	else if("`kernel'"=="epanechnikov") {
		qui gen W = (1/scalar(bwx))*(1-abs(dist_cut/scalar(bwx))^2)*(abs(dist_cut/scalar(bwx))<=1)
		qui gen Wtau = (1/scalar(bwx_tau))*(1-abs(dist_cut/scalar(bwx_tau))^2)*(abs(dist_cut/scalar(bwx_tau))<=1)
		if `gotcovs'==1{
			qui gen W_cov = (1/scalar(bwx_cov))*(1-abs(dist_cut/scalar(bwx_cov))^2)*(abs(dist_cut/scalar(bwx_cov))<=1)
			qui gen Wtau_cov = (1/scalar(bwx_tau_cov))*(1-abs(dist_cut/scalar(bwx_tau_cov))^2)*(abs(dist_cut/scalar(bwx_tau_cov))<=1)
		}
		matrix kernelconstants = [1/2, 3/16, 1/10, 1/16, 3/70 \ 1/2, -3/16, 1/10, -1/16, 3/70]
	}
	else{
		di as error  "{err}{cmd:kernel()} must be set to triangular, rectangular, gaussian, or epanechnikov"  
		exit
	}
	
	if `discretey'==0 {
		if "`kernel_y'" == "" {
			local kernel_y = "`kernel'"
		}
		else if "`kernel_y'"!="triangular" & "`kernel_y'"!="rectangular" & "`kernel_y'"!="gaussian" & "`kernel_y'"!="epanechnikov"{
			di as error  "{err}{cmd:kernel_y()} must be set to triangular, rectangular, gaussian, or epanechnikov"  
			exit
		}
	}
	
	*Polynomial orders
		tokenize `orders'	
		local w : word count `orders'
		if `w' == 0 {
			*Default to local linear:
			local order_tau=1
			local order_main=1
		}
		if `w' == 1 {
			local order_tau=`orders'
			local order_main=`orders'
		}
		if `w' == 2 {
			local order_tau=`1'
			local order_main=`2'
		}
		if `w' > 2 {
			di as error  "{err}{cmd:orders()} accepts at most two numbers"  
			exit 123
		}
		
		local order_max = max(`order_tau', `order_main')
		
		local polyvars
		local polyvars_tau
		
		qui gen polyvars0 = 1
		
		if `order_main' >= 1 {
			local polyvars dist_cut
			qui gen polyvars1 = dist_cut
		}
		if `order_tau' >= 1 {
			local polyvars_tau dist_cut
		}
				
		if `order_max' > 1 {
			local polyvars dist_cut
			foreach p of numlist 2/`order_max'{
				qui gen dist_cut`p' = dist_cut^`p'
				if `p' <= `order_main' {
					local polyvars `polyvars' dist_cut`p'
					qui gen polyvars`p' = dist_cut`p'
				}
				if `p' <= `order_tau' {
					local polyvars_tau `polyvars_tau' dist_cut`p'
				}
			}
		}
		*Note: the variables polyvars* and the varlist polyvars have the same content. The former is for use in mata, where we can't send a varlist.
		
		*Check ymin/ymax and yextremes
		qui summ y
		if "`ymin'"!="" {
			scalar ymin = `ymin'
			if scalar(ymin) > r(min) {
				disp "{err}{cmd:ymin()} There exist values of y smaller than ymin"
				exit
			}
			if "`ymax'"!="" {
				scalar ymax = `ymax'
				if scalar(ymax) < scalar(ymin) {
					disp "{err}{cmd:ymin()} ymin cannot be larger than ymax"
					exit
				}
			}
		}
		if "`ymax'"!="" {
			scalar ymax = `ymax'
			if scalar(ymax) < r(max) {
				disp "{err}{cmd:ymin()} There exist values of y larger than ymax"
				exit
			}
		}
		
		*Check for desire to compute refinements
		if ("`refinementA'"=="refinementA") {
			local dorefA=1
		}
		else{
			local dorefA=0
		}
		if ("`refinementB'"=="refinementB") {
			local dorefB=1
		}
		else{
			local dorefB=0
		}
		*Easter egg: "refinement C" sets tau0 = lambda=0 for treatment effects among units just to the right of the cutoff. 
		*It's no longer in paper and depreciated here, but can still be computed by setting dorefC to 1
		local dorefC=0

		if ("`righteffects'"=="righteffects") {
			local righteffects=1
			local w : word count `yextremes'
			if `w'==0{
				qui summ y
				di as result "Warning: {cmd:yextremes()} not specified while {cmd:right_effects} set, setting as range of y values in sample: `r(min)' to `r(max)'"
				scalar yextreme_min = r(min)
				scalar yextreme_max = r(max)
			}
			else if `w'==2 {
				tokenize `yextremes'
				scalar yextreme_min = `1'
				scalar yextreme_max = `2'
				qui summ y
				if scalar(yextreme_min) < r(min) {
					di as result "Warning: the lower bound in yextremes is smaller than the smallest y value in sample"
				}
				if scalar(yextreme_max) > r(max) {
					di as result "Warning: the upper bound in yextremes is larger than the largest y value in sample"
				}
				
				if scalar(yextreme_min) > r(min) {
					disp "{err}{cmd:yextremes()} there exist values of y smaller than the lower bound in yextremes"
					exit
				}
				if scalar(yextreme_max) < r(max) {
					disp "{err}{cmd:yextremes()} there exist values of y larger than the upper bound in yextremes"
					exit
				}
			}
			else{
				di as error "{err}{cmd:yxtremes()} expects two numbers"
				exit
			}
		}
		else{
			local righteffects=0
		}
		
		*Read in/determine "evaluation_ys", the values of Y at which cdfs/pdfs will be computed
		*First, get a list of the unique values of y (either to be used as evaluation_ys or to check against those given if discretey is set, to issue a possible warning)
		by y, sort: gen nvals = _n == 1 
		qui count if nvals
		local num_unique_ys=`r(N)'
		drop nvals
		
		qui query mem
		local matsizetoosmall = `num_unique_ys'>`r(matsize)'
		
		if `matsizetoosmall'==0{
			qui save "`tmpfile'", replace
					keep y dist_cut
					collapse dist_cut, by(y)
					keep y
					sort y
					mkmat y, matrix(unique_ys)
			qui use "`tmpfile'", clear
		}
			
		local w : word count `evaluation_ys'
		if `w' == 0 {
			if `matsizetoosmall'==1{
				di as error  "{err}{cmd:evaluation_ys()} is not set, so {cmd:rdbounds} is attempting to use all of the unique values of y to evaluate its distribution. However {cmd:matsize}=`r(matsize)' is too small for the `num_unique_ys' unique values of y. You must either increase {cmd:matsize} or specify a list of y values in {cmd: evaluation_ys}"  
				exit
			}
			local num_ys = `num_unique_ys'
			matrix evaluation_ys = unique_ys
		}
		else {
			qui query mem
			local num_ys = `w'
			if `num_ys'>`r(matsize)' {
				di as error  "{err}{cmd:evaluation_ys()} has more values than the maximum matrix size:`r(matsize)'. Please increase the maximum matrix size using {cmd: set matsize} or reduce the number of values in {cmd: evaluation_ys}"  
				exit
			
			}
		
			matrix evaluation_ys = J(`num_ys',1,0)
			
			forval i=1/`num_ys' {
				local el : word `i' of `evaluation_ys'
				mat evaluation_ys[`i',1]=`el'
			}
			
			*Warning #1					
			if `discretey' == 1 {
				if `matsizetoosmall'==1{
					di as result  "Warning: since {cmd:discretey} is set, {cmd:rdbounds} is attempting to compare the unique values of y to those specified in {cmd: evaluation_ys} to ensure that no values of y were skipped (if they are the calcualted probablity mass function over y will be innacurate). However {cmd:matsize}=`r(matsize)' is too small for the `num_unique_ys' unique values of y, so this diagnostic check is being skipped. If you want to run it, try increaseing {cmd:matsize}."  
				}
				else {
					matrix appended_ys = unique_ys
					matrix appended_ys = appended_ys \ evaluation_ys
					qui save "`tmpfile'", replace
						clear
						cap svmat appended_ys
						rename appended_ys y
						by y, sort: gen nvals = _n == 1 
						qui count if nvals
						local num_unique_appended_ys=r(N)
						drop nvals
					qui use "`tmpfile'", clear
					matrix drop appended_ys
					
					if `num_unique_appended_ys' > `num_ys' {
						di as result "Warning: when discrete_y is set to TRUE, a PMF function for y is be estimated by differencing the estimated CDF of y for consecutive values of evaluation_ys. However, the vector evaluation_ys does not contain all values of y detected in the sample, so this PMF will be innacurate. Re-running with evaluation_ys omitted recommended."
					}
				}
			}
			
			*Warning #2
			qui summ y
			if evaluation_ys[1,1] > r(min) {
				di as result "Warning: smallest value of evaluation_ys is greater than smallest value of y. CDFs will assume there is no mass to the left of this value. Suggested to re-run with evaluation_ys containing a lower bound for y."
			}
		}
	
		if `matsizetoosmall'==0{
			matrix drop unique_ys
		}
	
		if `righteffects'==1{
			if scalar(yextreme_min) < evaluation_ys[1,1] {
				di as result "Warning: the lower bound in yextremes is smaller than the smallest y value in evaluation_ys, so adding this value"
				local num_ys = `num_ys'+1
				matrix evaluation_ys = scalar(yextreme_min) \ evaluation_ys
			}
			if scalar(yextreme_max) > evaluation_ys[`num_ys',1] {
				di as result "Warning: the upper bound in yextremes is larger than the largest value of value in evaluation_ys, so adding this value"
				local num_ys = `num_ys'+1
				matrix evaluation_ys = evaluation_ys \ scalar(yextreme_max)
			}
		}
	
	*Number of bootstrap replications
	tokenize `num_bootstraps'	
		local w : word count `num_bootstraps'
		if `w' == 0 {
			local num_bs = 0
		}
		if `w' == 1 {
			local num_bs=`num_bootstraps'
		}
		if `w' > 1 {
			di as error  "{err}{cmd:num_bootstraps()} accepts at most one number"  
			exit 123
		}
	*Don't allow just one bootstrap, since we can't compute a standard deviation from this
	if `num_bs'==1{
		di as error "{err}{cmd:num_bootstraps()} num_bootstraps cannot be equal to one"
		exit
	}
		
	if "`kn'" == "" {
		scalar kn = sqrt(log(_N))
	}
	else{
		scalar kn = `kn'*1
		if !(scalar(kn) > 0) {
			di as error  "{err}{cmd:kn()} parameter {cmd:kn} must be a positive number"  
			exit
		}
	}
	
	scalar alpha = `alpha'
	if !(scalar(alpha) > 0 & scalar(alpha)<1 ) {
		di as error  "{err}{cmd:alpha()} parameter {cmd:alpha} must be between zero and one"  
		exit
	}
	
	if (scalar(alpha) >=0.5) {
		di as result "Warning: {cmd: alpha} >=0.5. For confidence intervals r_alpha calculated assuming that alpha<0.5. Since alpha>=0.5, true r_alpha may be very far from value used."
	}	
	
	*Initialize a matrix to store probabilities of each covariate value at the cutoff, for each bootstrap
	if `gotcovs'==1 {
		qui query mem
		if `num_bs'+1>`r(matsize)' {
			di as error "With covariate refinement turned on, {cmd:rdbounds} needs to create a matrix with {cmd: num_bs}+1 rows, but this is not possible due to {cmd: matsize}=`r(matsize)'. Please increase the maximum matrix size using {cmd: set matsize}, reduce the number of bootstrap replications, or turn off the covariate refinement."
			exit
		}			
		matrix probw=J(`ncovlevels', `num_bs'+1, .)
	}
	
	
	*Check type and percentiles parameters
	if "`type'" == "" | "`type'"=="ate" {
		if "`percentiles'"=="" {
			local percentiles = "-1"
		}
		else {
			di as result  "Warning: {cmd:type} is set to ate or is blank, but {cmd:percentiles} is set. Treatment effects will be set based on {cmd:percentiles}"  
		}
	}
	if "`type'" == "qte" {
		if "`percentiles'" == "" {
			di as error  "{err}{cmd:percentiles()} parameter {cmd:type} is set to qte, but percentiles at which to compute quantile treatment effects not set"  
			exit
		}
	}
	
	if "`percentiles'" != "" {
		foreach p of numlist `percentiles' {
			if (floor(`p')!=`p')|(`p'==0)|(`p'<-1)|(`p'>99) {
				di as error  "{err}{cmd:percentiles()} values must be integers between 1 and 99 or -1 (for ate)"  
				exit 125
			}
		}
	}
	
	*Check potential taus
	if "`potential_taus'" != "" {
		foreach t of numlist `potential_taus' {
			if (`t'<0)|(`t'>=1) {
				di as error  "{err}{cmd:potential_taus()} values must be greater than or equal to zero and less than one"  
				exit 125
			}
		}
	}
	
	*Get some functions ready for doing collapses in mata:
	local hasate=1
	local hasqte=1
	if `hasate'==1{
		cap mata: mata drop my_sum()
		cap mata: function my_sum(x,ignoreme) return(sum(x))
	}
	if `hasqte'==1{
		cap mata: mata drop my_min()
		cap mata: function my_min(x,ignoreme) return(min(x))
	}
	
	*Get tau values provided by user for CI's based on the potential impact of manipulation:
	local num_taus : word count `potential_taus'
	if `num_taus'>0 {
		matrix potential_taus = J(`num_taus',1,0)
		forval i=1/`num_taus' {
			local el : word `i' of `potential_taus'
			mat potential_taus[`i',1]=`el'
		}
	}
	
	*Detect if user requested an outputfile
	local saveToFile = 0
	local w: word count `outputfile'
	if `w'>0 {
		local saveToFile = 1
		  capture confirm new file "`outputfile'"
		  if _rc==602 {
			if ("`replace'"!="replace") {
				disp as error "{err} {cmd: outputfile()}: the file `outputfile' already exists. To overwrite it, please set the {cmd: replace} option"
				exit
			}
		  }
		  	if _rc==7 {
				disp as error "{err} {cmd: outputfile()}: the file `outputfile' does not appear to be a valid filename. Please make appropriate changes"
				exit
		  }
	}
	
	local saveToFileInputs = 0
	local w: word count `outputfileInputs'
	if `w'>0 {
		local saveToFileInputs = 1
		  capture confirm new file "`outputfileInputs'"
		  if _rc==602 {
			if ("`replace'"!="replace") {
				disp as error "{err} {cmd: outputfile()}: the file `outputfileInputs' already exists. To overwrite it, please set the {cmd: replace} option"
				exit
			}
		  }
		  	if _rc==7 {
				disp as error "{err} {cmd: outputfile()}: the file `outputfileInputs' does not appear to be a valid filename. Please make appropriate changes"
				exit
		  }
	}
	
	if ("`verbose'"=="verbose") {
		local verbose = 1
	}
	else {
		local verbose = 0
	}

*First, loop through bootstrap samples to calculate scalars (tau, treatLeft, treatRight), 
*and estimated inputs to CDF computation. Further steps will no longer require
*access to the samples themselves
********************************************************************************
********************************************************************************
********************************************************************************
qui save "`tmpfilesample'", replace

*Begin bootstrap/covariate loop
********************************************************************************
********************************************************************************

local num_bs_plusone = `num_bs'+1

forvalues b = 0/`num_bs' {
  
  *Start probSum at zero (this is for normalizing probailities of covariate values)
  scalar probSum = 0
  
  forvalues w = 0/`ncovlevels' {
  
	*Load bandwidths for covariate-conditional quantities (if not set to be different these just default to the specified full-sample bandwidths)
	qui capture drop Wtau_current
	if `w'==0{
		scalar hy = scalar(bwy)
		scalar hxtau = scalar(bwx_tau)
		local Wvar = "W"
		local Wtauvar = "Wtau"
	}
	if `w'>0{
		if `discretey'==0 {
			scalar hy = scalar(bwy_cov)
		}
		scalar hxtau = scalar(bwx_tau_cov)
		local Wvar = "W_cov"
		local Wtauvar = "Wtau_cov"
	}

	qui use "`tmpfilesample'", clear
	
	*b=0 are point estimates with full sample, for b>0 resample:
	if `b' > 0 {
		bsample
	}
	
	*Compute covariate probabilities if user selected covariate refinement (need to do this while we have the sample in memory)
	********************************************************************************
	if `gotcovs'==1 & `w' > 0{
			qui gen isw = (cov==`w')
			qui reg isw `polyvars' [aw=W] if dist_cut<0
			scalar probLimit = _b[_cons]
			drop isw
			if `b'==0 & probLimit < -0.00000001 {
				di as result "Warning: estimated left limit of prob(covs=w|x) at the cutoff for w=`w' is less than zero, replacing with zero"
				scalar probLimit = 0
			}
			if `b'==0 & probLimit > 1.00000001 {
				di as result "Warning: estimated left limit of prob(covs=w|x) at the cutoff for w=`w' is greater than one, replacing with one"
				scalar probLimit = 1
			}
			matrix probw[`w',`b'+1] = scalar(probLimit)
			scalar probSum = scalar(probSum)+scalar(probLimit)
	}

	
	*w=0 is estimation with full sample, w>0 subsampling by covariates value
	if `w' > 0 {
		qui keep if cov==`w'
	}
	
	*Estimate treatment probabilities
	********************************************************************************
		qui reg treatment `polyvars' [aw=`Wvar'] if dist_cut<0
		scalar treatLeft = _b[_cons]
		qui reg treatment `polyvars' [aw=`Wvar'] if dist_cut>=0
		scalar treatRight = _b[_cons]
				
	*Estimate tau
	********************************************************************************
		if `discretex' == 1 {
			qui bysort dist_cut: gen freq = _N if Wtau>0
			qui reg freq `polyvars_tau' [aw=`Wtauvar'] if dist_cut>=0
			scalar Fplus = _b[_cons]
			qui reg freq `polyvars_tau' [aw=`Wtauvar'] if dist_cut<0
			scalar Fminus = _b[_cons]
			drop freq
		} 
		else{
			if `order_tau' == 0 {
				qui gen tmp=(dist_cut >=0)*`Wtauvar'/kernelconstants[1,1]
				qui summ tmp, meanonly
				scalar Fplus = r(mean)
				drop tmp
				qui gen tmp=(dist_cut <0)*`Wtauvar'/kernelconstants[2,1]
				qui summ tmp, meanonly
				scalar Fminus = r(mean)
				drop tmp
			}
			if `order_tau' == 1 {
				qui gen tmp=(dist_cut >=0)*`Wtauvar'*(kernelconstants[1,3]-dist_cut/scalar(hxtau)*kernelconstants[1,2])/(kernelconstants[1,1]*kernelconstants[1,3]-kernelconstants[1,2]^2)
				qui summ tmp, meanonly
				gen tmp2 = (kernelconstants[1,3]-dist_cut/scalar(hxtau)*kernelconstants[1,2])
				scalar Fplus = r(mean)
				drop tmp
				qui gen tmp=(dist_cut<0)*`Wtauvar'*(kernelconstants[2,3]-dist_cut/scalar(hxtau)*kernelconstants[2,2])/(kernelconstants[2,1]*kernelconstants[2,3]-kernelconstants[2,2]^2)
				qui summ tmp, meanonly
				scalar Fminus = r(mean)
				drop tmp
			}
			if `order_tau' == 2 {
				local det = kernelconstants[1,1]*(kernelconstants[1,3]*kernelconstants[1,5]-kernelconstants[1,4]^2)-kernelconstants[1,2]*(kernelconstants[1,2]*kernelconstants[1,5]-kernelconstants[1,3]*kernelconstants[1,4])+kernelconstants[1,3]*(kernelconstants[1,2]*kernelconstants[1,4]-kernelconstants[1,3]^2)
				qui gen tmp=(dist_cut >=0)*`Wtauvar'/`det'*((kernelconstants[1,3]*kernelconstants[1,5]-kernelconstants[1,4]^2)+(kernelconstants[1,3]*kernelconstants[1,4]-kernelconstants[1,2]*kernelconstants[1,5])*(dist_cut/scalar(hxtau))+(kernelconstants[1,2]*kernelconstants[1,4]-kernelconstants[1,3]^2)*(dist_cut^2/(scalar(hxtau)^2)))
				qui summ tmp, meanonly
				scalar Fplus = r(mean)
				drop tmp
				local det = kernelconstants[2,1]*(kernelconstants[2,3]*kernelconstants[2,5]-kernelconstants[2,4]^2)-kernelconstants[2,2]*(kernelconstants[2,2]*kernelconstants[2,5]-kernelconstants[2,3]*kernelconstants[2,4])+kernelconstants[2,3]*(kernelconstants[2,2]*kernelconstants[2,4]-kernelconstants[2,3]^2)
				qui gen tmp=(dist_cut <0)*`Wtauvar'/`det'*((kernelconstants[2,3]*kernelconstants[2,5]-kernelconstants[2,4]^2)+(kernelconstants[2,3]*kernelconstants[2,4]-kernelconstants[2,2]*kernelconstants[2,5])*(dist_cut/scalar(hxtau))+(kernelconstants[2,2]*kernelconstants[2,4]-kernelconstants[2,3]^2)*(dist_cut^2/scalar(hxtau)^2))
				qui summ tmp, meanonly
				scalar Fminus = r(mean)
				drop tmp
			}
			if `order_tau' > 2 {
				di as error  "{err}{cmd:orders()} For estimation of discontinuity in density at cutoff (tau_hat) must be of quadratic order or lower if discretex is not set."  
				exit
			}
		}
		scalar tau_tilde = 1-scalar(Fminus)/scalar(Fplus)
		scalar tau_hat = max(0, scalar(tau_tilde))		
		
		local tau_hat = round(scalar(tau_hat), .00001)
		if(`b'==0 & `w'==0){
			di  "The proportion of always-assigned units just to the right of the cutoff is estimated to be `tau_hat'"
		}
		
		if scalar(tau_hat)==0 & `b'==0 & `w'==0{
			di as result  "Warning: No evidence of manipulation found, thus bounds based on estimated tau will not differ from estimation that assumes no manipulation"
		}
		
		if scalar(tau_hat)==0 & `b'==0 & `w'>0{
			di ""
			di as result  "<Note: no evidence of manipulation found for covariate subsample cov=`w'>"
		}
	
		*Now that we've displayed tau_hat for the main sample (b=0), let's start the counter which shows progress over bootstrap replications
		*We want this to display once per bootstrap, even if there is a covariate refinement, so we'll only do this when w=0 (full sample)
		
		if `w'==0{
		  if `b'==-1 {
		  
		  }
		  else if `b'==0 {
			di as result "Calculating CDF inputs. Dots indicate bootstrap repetitions completed (out of num_bootraps+1, including one iteration for the full sample)"
			di "CDF inputs progress (out of `num_bs_plusone'): " _cont
			di "*" _cont
		  }
		  else if floor(`b'/100) == (`b'/100) {
			di "|" _cont
		  }
		  else if floor(`b'/5) == (`b'/5) {
			di "+" _cont
		  }
		  else {
			di "." _cont
		  }
		}
	
		*Compute raw CDF inputs
		********************************************************************************
		*Note: regressions are implemented as explicit matrix equations in mata because it vastly improves performance for large `num_ys'
		mata: inputs_mata=J(`num_ys', 9, .)
		
		gen indicate_right = (dist_cut >= 0)
		gen indicate_left = (dist_cut < 0)
		gen indicate_right_treated = (dist_cut >= 0)*(treatment==1)
		gen indicate_right_untreated = (dist_cut >= 0)*(treatment==0)
		gen indicate_left_treated = (dist_cut < 0)*(treatment==1)
		gen indicate_left_untreated = (dist_cut < 0)*(treatment==0)		
		
		forvalues i = 1/`num_ys' {
			scalar little_y = evaluation_ys[`i',1]
			gen below_y = (y<=scalar(little_y))
			mata: inputs_mata[`i',1] = st_numscalar("little_y")
			
			mata: M=X=below_y=W=.; st_view(M, ., ("below_y", "`Wvar'", "polyvars*"), "indicate_right"); st_subview(below_y, M, ., 1); st_subview(W, M, ., 2); st_subview(X, M, ., (3\.))
			mata: XX = cross(X,W,X); Xbelowy = cross(X,W,below_y); b = cholsolve(XX,Xbelowy)
			mata: inputs_mata[`i',2] = b[1]
			
			mata: M=X=below_y=W=.; st_view(M, ., ("below_y", "`Wvar'", "polyvars*"), "indicate_left"); st_subview(below_y, M, ., 1); st_subview(W, M, ., 2); st_subview(X, M, ., (3\.))
			mata: XX = cross(X,W,X); Xbelowy = cross(X,W,below_y); b = cholsolve(XX,Xbelowy)
			mata: inputs_mata[`i',3] = b[1]
			
			if `fuzzy'==1 {
				mata: M=X=below_y=W=.; st_view(M, ., ("below_y", "`Wvar'", "polyvars*"), "indicate_right_treated"); st_subview(below_y, M, ., 1); st_subview(W, M, ., 2); st_subview(X, M, ., (3\.))
				mata: XX = cross(X,W,X); Xbelowy = cross(X,W,below_y); b = cholsolve(XX,Xbelowy)
				mata: inputs_mata[`i',4] = b[1]
				
				local norightuntreated = 1
				qui count if (dist_cut >= 0 & treatment==0 & W>0)
				if `r(N)' > 0 {
					local norightuntreated = 0
					mata: M=X=below_y=W=.; st_view(M, ., ("below_y", "`Wvar'", "polyvars*"), "indicate_right_untreated"); st_subview(below_y, M, ., 1); st_subview(W, M, ., 2); st_subview(X, M, ., (3\.))
					mata: XX = cross(X,W,X); Xbelowy = cross(X,W,below_y); b = cholsolve(XX,Xbelowy)
					mata: inputs_mata[`i',5] = b[1]
				}
				
				local nolefttreated = 1
				qui count if (dist_cut < 0 & treatment==1 & W>0)
				if `r(N)' > 0 {
					local nolefttreated = 0
					mata: M=X=below_y=W=.; st_view(M, ., ("below_y", "`Wvar'", "polyvars*"), "indicate_left_treated"); st_subview(below_y, M, ., 1); st_subview(W, M, ., 2); st_subview(X, M, ., (3\.))
					mata: XX = cross(X,W,X); Xbelowy = cross(X,W,below_y); b = cholsolve(XX,Xbelowy)
					mata: inputs_mata[`i',6] = b[1]
				}
				
				mata: M=X=below_y=W=.; st_view(M, ., ("below_y", "`Wvar'", "polyvars*"), "indicate_left_untreated"); st_subview(below_y, M, ., 1); st_subview(W, M, ., 2); st_subview(X, M, ., (3\.))
				mata: XX = cross(X,W,X); Xbelowy = cross(X,W,below_y); b = cholsolve(XX,Xbelowy)
				mata: inputs_mata[`i',7] = b[1]
			
				*Separately estimate untreated densities on either side of the cutoff, unless discretey is set (in which case PMF's will be estimated by differencing CDFs)
				if `discretey'==0 {
					scalar kernel_left = 1/2
					scalar kernel_right = 1/2
									
					if ("`kernel_y'"=="triangular"){
						qui gen Khy = (1/scalar(hy))*(1-abs((y-scalar(little_y))/scalar(hy)))*(abs((y-scalar(little_y))/scalar(hy))<=1)
						if "`ymin'"!=""  {
							scalar xx = (scalar(ymin)-scalar(little_y))/scalar(hy)
							if scalar(xx) <= -1 {
								scalar kernel_left = 1/2-0
							}
							if scalar(xx)>-1 & scalar(xx)<1 {
								scalar kernel_left = 1/2-(scalar(xx)-1/2*scalar(xx)^2*sign(scalar(xx))+1/2)
							}
							if scalar(xx) >= 1{
								scalar kernel_left = 1/2-1
							}	
						}
						if "`ymax'"!=""  {
							scalar xx = (scalar(ymax)-scalar(little_y))/scalar(hy)
							if scalar(xx) <= -1 {
								scalar kernel_right = 0-1/2
							}
							if scalar(xx)>-1 & scalar(xx)<1 {
								scalar kernel_right = (scalar(xx)-1/2*scalar(xx)^2*sign(scalar(xx))+1/2)-1/2
							}
							if scalar(xx) >= 1{
								scalar kernel_right = 1-1/2
							}
						}					
					}
					else if("`kernel_y'"=="rectangular") {
						qui gen Khy = (1/scalar(hy))*(1/2)*(abs((y-scalar(little_y))/scalar(hy))<=1)
						if "`ymin'"!=""  {
							scalar xx = (scalar(ymin)-scalar(little_y))/scalar(hy)
							if scalar(xx) <= -1 {
								scalar kernel_left = 1/2-0
							}
							if scalar(xx)>-1 & scalar(xx)<1 {
								scalar kernel_left = 1/2-(1/2*scalar(xx)+1/2)
							}
							if scalar(xx) >= 1{
								scalar kernel_left = 1/2-1
							}	
						}
						if "`ymax'"!=""  {
							scalar xx = (scalar(ymax)-scalar(little_y))/scalar(hy)
							if scalar(xx) <= -1 {
								scalar kernel_right = 0-1/2
							}
							if scalar(xx)>-1 & scalar(xx)<1 {
								scalar kernel_right = (1/2*scalar(xx)+1/2)-1/2
							}
							if scalar(xx) >= 1{
								scalar kernel_right = 1-1/2
							}
						}
					}
					else if("`kernel_y'"=="gaussian") {
						qui gen Khy = (1/scalar(hy))*(1/sqrt(2*c(pi)))*exp(-((y-scalar(little_y))/scalar(hy))^2/2)
						if "`ymin'"!=""  {
							scalar xx = (scalar(ymin)-scalar(little_y))/scalar(hy)
							scalar kernel_left = 1/2-normal(scalar(xx))
						}
						if "`ymax'"!=""  {
							scalar xx = (scalar(ymax)-scalar(little_y))/scalar(hy)
							scalar kernel_right = normal(scalar(xx))-1/2
						}
					}
					else if("`kernel_y'"=="epanechnikov") {
						qui gen Khy = (1/scalar(hy))*(1-abs((y-scalar(little_y))/scalar(hy))^2)*(abs((y-scalar(little_y))/scalar(hy))<=1)
						if "`ymin'"!=""  {
							scalar xx = (scalar(ymin)-scalar(little_y))/scalar(hy)
							if scalar(xx) <= -1 {
								scalar kernel_left = 1/2-0
							}
							if scalar(xx)>-1 & scalar(xx)<1 {
								scalar kernel_left = 1/2-(-3*(scalar(xx)^3/3-scalar(xx))/4+1/2)
							}
							if scalar(xx) >= 1{
								scalar kernel_left = 1/2-1
							}	
						}
						if "`ymax'"!=""  {
							scalar xx = (scalar(ymax)-scalar(little_y))/scalar(hy)
							if scalar(xx) <= -1 {
								scalar kernel_right = 0-1/2
							}
							if scalar(xx)>-1 & scalar(xx)<1 {
								scalar kernel_right = (-3*(scalar(xx)^3/3-scalar(xx))/4+1/2)-1/2
							}
							if scalar(xx) >= 1{
								scalar kernel_right = 1-1/2
							}
						}
					}
							
					scalar kernel_coverage = scalar(kernel_left)+scalar(kernel_right)
							
					qui replace Khy = Khy/scalar(kernel_coverage)
					
					if `norightuntreated' == 0 {
						mata: M=X=kernel_y=W=.; st_view(M, ., ("Khy", "`Wvar'", "polyvars*"), "indicate_right_untreated"); st_subview(kernel_y, M, ., 1); st_subview(W, M, ., 2); st_subview(X, M, ., (3\.))
						mata: XX = cross(X,W,X); Xkernely = cross(X,W,kernel_y); b = cholsolve(XX,Xkernely)
						mata: inputs_mata[`i',8] = b[1]
					}
					else {
						mata: inputs_mata[`i',8] = .
					}
					
					mata: M=X=kernel_y=W=.; st_view(M, ., ("Khy", "`Wvar'", "polyvars*"), "indicate_left_untreated"); st_subview(kernel_y, M, ., 1); st_subview(W, M, ., 2); st_subview(X, M, ., (3\.))
					mata: XX = cross(X,W,X); Xkernely = cross(X,W,kernel_y); b = cholsolve(XX,Xkernely)
					mata: inputs_mata[`i',9] = b[1]
					
					drop Khy
				}
				else {
					mata: inputs_mata[1,8] = .
					mata: inputs_mata[1,9] = .
				}				
			}
			drop below_y
		}
		
		drop indicate_*
		
		mata: st_matrix("inputs",inputs_mata)
		matrix colnames inputs="little_y" "F_right" "F_left" "F_right_treated" "F_right_untreated" "F_left_treated" "F_left_untreated" "dens_right_untreated" "dens_left_untreated"
		
		*Make estimated CDFs proper CDF functions
		clear
		qui svmat inputs, names(col)
		
			sort little_y
			gen rank_y = _n
		
			foreach F of varlist F_right F_left F_right_treated F_right_untreated F_left_treated F_left_untreated {
				*Censor CDFs to unit interval
					qui replace `F' = 1 if `F' > 1 & `F' != .
					qui replace `F' = 0 if `F' < 0 & `F' != .
				*Normalize CDFs
					qui summ `F'
					if r(N) > 0 {
						qui replace `F' = `F'/r(max) if `F' != .
					}
				*Monotonize all columns of inputs
					sort `F'
					gen rank = _n
					qui save "`tmpfile'", replace
					rename (`F' rank) (`F'_old rank_old)
					rename rank_y rank
					qui merge 1:1 rank using "`tmpfile'", keepusing(`F')
					
					rename rank rank_y
					drop _merge `F'_old rank_old
					
			}
			
			drop rank_y
			sort little_y			
			
			*Finalize density functions
			if `discretey'==1 {
				*if discretey is set, compute PMFs by differencing corresponding CDFs. Since this is done after making these proper CDF functions, PMFs will already be normalized
				qui replace dens_right_untreated = F_right_untreated[_n] - F_right_untreated[_n-1]
				qui replace dens_right_untreated = F_right_untreated if _n==1
				
				qui replace dens_left_untreated = F_left_untreated[_n] - F_left_untreated[_n-1]
				qui replace dens_left_untreated = F_left_untreated if _n==1
			}
			else {
				*If discretey is not set, PDF's have been estimated directly, so make sure PDFs are proper PDF functions
				foreach f of varlist dens_right_untreated dens_left_untreated {
				*Censor PDFs at zero
					qui replace `f' = 0 if `f' < 0 & `f' != .
				*Normalize PDFs
					sort little_y
					qui gen lengths = little_y[_n] - little_y[_n-1]
					qui replace lengths = 0 if _n==1
					qui gen temp = lengths*`f'
					qui summ temp, meanonly
					qui replace `f' = `f'/r(sum) if `f' != .
					drop lengths temp
				}
			}
			
		*W indicates a value of covs if this option is passed
		if `w' > 0 {
			qui gen w = `w'
		}
		else {
			qui gen w = 0
		}

		qui gen b = `b'
		
		*We'll store these scalars as variables(even though its inefficient), because it will be useful to work with them across bootstraps and covariate values later, and makes the code much cleaner
		qui gen treatRight = scalar(treatRight)
		qui gen treatLeft = scalar(treatLeft)
		qui gen tau_tilde = scalar(tau_tilde)
		qui gen tau_hat = scalar(tau_hat)
		
		if `b'==0 & `num_bs' > 0{
			*If b=0 (full sample), expand now to store nudged tau as b=-1
			qui expand 2
			qui bysort w little_y: replace b = -1 if _n==1
		}
		if `b'==0{
			scalar tau_hat_pt_`w' = scalar(tau_hat)
			scalar tau_tilde_pt_`w' = scalar(tau_tilde)
		}
		if `b' > 0 | `w' > 0 {
			qui append using "`tmpfileinputs'"
		}
		qui save "`tmpfileinputs'", replace
		
*End bootstrap/covariate loop
********************************************************************************
********************************************************************************
  }
  
  *Normalize vector of probabilities for this bootstrap
  if `gotcovs'==1{
	forvalues w = 1/`ncovlevels' {
		matrix probw[`w',`b'+1] = probw[`w',`b'+1]/scalar(probSum)
	}
  }
			
}
order b w little_y
sort b w little_y

disp _newline

*Store point estimates and standard deviations of tau and takeup increase for later
********************************************************************************
	qui gen takeup_increase = treatRight-treatLeft
	qui summ takeup_increase if b==0 & w==0 & little_y==evaluation_ys[1,1]
	scalar takeup_increase_pt = r(mean)
	qui summ takeup_increase if b > 0& w==0 & little_y==evaluation_ys[1,1]
	scalar sdtakeup = r(sd)
	drop takeup_increase

	qui summ tau_hat if b==0 & w==0 & little_y==evaluation_ys[1,1]
	scalar tau_hat_pt = r(mean)
	qui summ tau_hat if b > 0& w==0 & little_y==evaluation_ys[1,1]
	scalar sdtauhat = r(sd)

**Check for crazy tau values (>=1), which may happen with bad bandwidth choices
********************************************************************************
	qui summ tau_hat if b>-1
	if r(max) > 1 {
		di as result "{err} an estimated tau is greater than one, cannot proceed as is. Your bandwidths may not be set appropriately."
		di as result "The table below indicates tau_hat by bootstrap/covariate value. w=0 indicates the full sample (not conditional on covariate), and b=0 indicates the original sample (not a bootstrap sample):"
		di as result "lowest value of tau0U:"
		table b w, c(mean tau_hat)
		exit
	}	
	
*Compute "nudged" taus as necessary
********************************************************************************
	qui gen sdtautilde = .
	forvalues w = 0/`ncovlevels' {
		qui sum tau_tilde if b > 0 & w==`w' & little_y==evaluation_ys[1,1]
		qui replace sdtautilde = r(sd) if w==`w'
	}
	qui gen tau_star = max(tau_tilde, scalar(kn)*sdtautilde)
	
	qui gen tau = tau_hat
	qui replace tau = tau_star if b==-1
	
	forvalues w = 0/`ncovlevels' {
		qui replace tau = max(0,tau_tilde-scalar(tau_tilde_pt_`w')+max(scalar(tau_hat_pt_`w'), scalar(kn)*sdtautilde)) if b > 0 & w==`w'
		qui summ tau if b > 0 & w==`w'
		if (r(max)>1 & `num_bs' > 0) {
			di as result "Warning: the bootstrap nudge (see step 3 of inference section of paper) for w==`w' (w=0 indicates full sample) resulted in tau values greater than one, replacing with one"
			tab b if tau>1 & w==`w'
			qui replace tau = min(1, tau) if b > 0 & w==`w'
		}
	}
	qui drop sdtautilde
	
	*Save tau_star for w=0 for the tau confidence interval
	qui summ tau_star if b==-1 & w==0
	scalar tau_star=r(mean)
	
	*Use standard deviation of nudged taus for the tau confidence interval, so save this as a scalar now
	qui summ tau if b > 0& w==0 & little_y==evaluation_ys[1,1]
	scalar sdtauhatnudged = r(sd)
	
*Bring in any fixed taus provided by user
********************************************************************************
	qui gen fixedtau = 0
	order fixedtau
	qui expand `num_taus'+1
	qui bysort w b little_y: replace fixedtau = _n-1 if _n > 1
	*fixedtau=0 corresponds to tau estimated from data. Since we aren't nudging taus for fixedtau, we could drop either b=0 or b=-1 for fixedtau>0, since these are
	*redundant. But we'll keep them to keep the code the same as the other cases, since we'll compute the RAM-intensive refinements (righteffects and covariate) only for
	*estimated taus
	if `num_taus' > 0 {
		qui replace tau = potential_taus[fixedtau,1] if fixedtau > 0
	}
	*Don't need covariate CDFs for fixed tau
	qui drop if fixedtau > 0 & w > 0
	
*Do tau testing
********************************************************************************
	if `num_bs'>1{
		qui gen theta = tau_hat - (1-treatRight)
		qui egen sdtheta = sd(theta) if little_y==evaluation_ys[1,1], by(w fixedtau)
		qui count if abs(theta)<1.96*sdtheta & sdtheta!=.
		if r(N) > 0	{
			di as result "Warning: the constraint tau != P(D=0|c+) is close to binding (cannot reject the null hypothesis of equality) in the following cases (w=0 indicates full sample, not covariate conditional). CI's may not have correct coverage. "
			qui gen warninghere = (abs(theta)<1.96*sdtheta & sdtheta!=.)
				table fixedtau w, c(max warninghere)
			drop warninghere
		}
		drop theta sdtheta
	}
		
sort fixedtau w b little_y
qui save "`tmpfileinputs'", replace

*Second loop over bootstraps to get treatment effects
********************************************************************************
********************************************************************************
di ""
di as result "Calculating Treatment Effects"
di as result "-----------------------------"

local therewerewarnings=0
local therewereerrors=0

local firstb = 0
if `num_bs' > 1 {
	local firstb = -1
}
forvalues b = `firstb'/`num_bs' {

	*Check for errors from previous loop iterations
	if(`therewereerrors'==1){
		exit
	}

	qui use "`tmpfileinputs'" if b==`b', clear	
	
	*Compute CDFs for potential outcomes
	********************************************************************************		
		
		qui gen keepthisbootstrap = 1
			
		*First, do some prep work for the FRD before we expand the dataset:
		if `fuzzy'==1 {
								
				*Generate CDF G:
				*---------------
					qui gen kappa1 = (1-tau)*treatLeft/treatRight

					if `nolefttreated' == 0 {
						qui gen G = (F_right_treated-kappa1*F_left_treated)/(1-kappa1)
					}
					else {
						*If there are no treated units on the left, then kappa1 = 0 and the CDF G is just equal to F_right_treated
						qui gen G = F_right_treated
					}
					
					*G will already be normalized, but may not be monotonic or within the unit interval. Correct these problems if they occur and warn.
						qui count if G < 0 | G > 1
						if r(N) > (.02*_N) & "`F'"=="G"  & `b'==0 {
							local therewerewarnings=1
							di as result "Warning: : the function G(y) (see paper for definition) should be a CDF, but it contains more than 2% of values outside of the unit interval. Values have been censored to the unit interval."
							if `verbose'==1 {								
								di as result "--Positive numbers indicate values of G(y) outside the unit interval for that fixedtau/covariate value cell, for main sample (not bootstrap sample). w=0 indicates the full sample (not conditional on covariate), andfixedtau=0 indicates tau is estimated, not fixed via {cmd:potential_taus}:"
								qui gen warninghere = (G < 0 | G > 1)
									table fixedtau w, c(max warninghere)
								drop warninghere
							}	
						}
						qui replace G = 1 if G > 1 & G != .
						qui replace G = 0 if G < 0 & G != .
						
						sort fixedtau w b little_y
						qui by fixedtau w b: gen rank_y = _n
						sort fixedtau w b G little_y
						qui by fixedtau w b: gen rank = _n
						qui save "`tmpfile'", replace
						rename (G rank) (G_old rank_old)
						rename rank_y rank
						qui merge 1:1 fixedtau w b rank using "`tmpfile'", keepusing(G)
						qui count if rank!=rank_old 
						if r(N) > (.02*_N)  & `b'==0{{
							local therewerewarnings=1
							di as result "Warning: the function G(y) (see paper for definition) should be a CDF, but is not completely monotonic (at least 2% reordered). Values have been monotonized, but this could be evidence against the model."							
							if `verbose'==1 {
								di as result "--Positive numbers indicate non-monotonicities in that fixedtau/covariate value cell, for main sample (not bootstrap sample). w=0 indicates the full sample (not conditional on covariate), and fixedtau=0 indicates tau is estimated, not fixed via {cmd:potential_taus}:"
								qui gen warninghere = rank!=rank_old
									table fixedtau w, c(max warninghere)
								drop warninghere
								qui correlate rank rank_old if fixedtau==0 & w==0
								di as result "--Correlation between initial rank and rank after monotonization (for fixedtau=0 & w=0): `r(rho)'"
							}				
						}
						rename rank rank_y
						drop _merge G_old rank_old
						drop rank_y
						sort little_y
						
				*Generate density s:
				*-------------------
					qui gen kappa0 = 1/(1-tau)*(1-treatRight)/(1-treatLeft)
									
					if `norightuntreated' == 0 {
						qui gen s_notau = dens_right_untreated
						qui replace s_notau = dens_left_untreated/kappa0 if dens_left_untreated/kappa0 < dens_right_untreated
												
						*Calculate mass of s_notau function
						if `discretey'==0 {
							sort fixedtau w b little_y
							qui by fixedtau w b: gen lengths = little_y[_n] - little_y[_n-1]
							qui replace lengths = 0 if _n==1
							qui gen temp = lengths*s_notau
							drop lengths
						}
						else{
							qui gen temp = s_notau
						}
						qui egen Pi = sum(temp), by(fixedtau w b)
						drop temp
						
						*Warnings/errors for mass of s_notau (which we're calling Pi)						
						qui gen tau0L = min(1, tau/(1-treatRight))
						
						*Check if T set is null in main sample (not bootstrap)
						qui count if Pi < 1-tau0L & tau>0
						if r(N) > 0 & `b'==0 {
						
							local therewereerrors=1
							
							qui count if Pi < 1-tau0L & w==0 & fixedtau==0 & tau>0
							if r(N) > 0 {
								di as error "Error: The estimated identified set for tau0 is null. This occurs when s(y)*(1-tau0) integrates to less than 1-tau0L (see paper for definitions), making the lower bound for tau0 larger than the upper bound. Use the {cmd:verbose} command for more details."
							}
							else {
								qui count if Pi < 1-tau0L & w!=0 & tau>0
								if r(N) > 0 {
									di as error "Error: The estimated identified set for tau0 is null for some values of covariates. This occurs when, conditional on a value of {cmd:covs}, s(y)*(1-tau0) integrates to less than 1-tau0L (see paper for definitions), making the lower bound for tau0 larger than the upper bound. You'll need to re-run without the {cmd: covs} option specified. Use the {cmd:verbose} command for more details."
								}
							}
							else {
								qui count if Pi < 1-tau0L & fixedtau!=0 & tau>0
								if r(N) > 0 {
									di as error "Error: The estimated identified set for tau0 is null for some pre-specified values of {cmd: potential_taus}. You'll need to re-run with different values of {cmd: potential_taus}, or this option omitted. This occurs when, conditional on a value of {cmd:covs}, s(y)*(1-tau0) integrates to less than 1-tau0L (see paper for definitions), making the lower bound for tau0 larger than the upper bound. Use the {cmd:verbose} command for more details."
								}								
							}
							
							if `verbose'==1 { 		
								di as result "--------The tables below indicate the total mass of s(y)*(1-tau0), and the value of tau0L, respectively (for the main sample, not a bootstrap resample). w=0 indicates the full sample (not conditional on covariate), and fixedtau=0 indicates tau is estimated, not fixed via {cmd:potential_taus}:"
								di as result "mass of s(y)*(1-tau0):"
								table fixedtau w, c(mean Pi)
								di as result "value of tau0L:"
								table fixedtau w, c(mean tau0L)
							}	
							
							exit
						}
						*Check if T set is null in nudged-tau "sample"
						if r(N) > 0 & `b'==-1 {
						
							local therewereerrors=1
							
							di as error "Error: The estimated identified set for tau0 is null, for starred estimands in notation of the paper (see the use of bootstrap to nudge tau away from zero in section on inference). This occurs when s(y)*(1-tau0) integrates to less than 1-tau0L (see paper for definitions), making the lower bound for tau0 larger than the upper bound. You may wish to re-run to generate new bootstrap resamples. Use the {cmd:verbose} command for more details."							
							
							if `verbose'==1 { 		
								di as result "--------The tables below indicate the total mass of s(y)*(1-tau0), and the value of tau0L, respectively (for the 'starred' version of estimands). w=0 indicates the full sample (not conditional on covariate), and fixedtau=0 indicates tau is estimated, not fixed via {cmd:potential_taus}:"
								di as result "mass of s(y)*(1-tau0):"
								table fixedtau w, c(mean Pi)
								di as result "value of tau0L:"
								table fixedtau w, c(mean tau0L)
							}	
							
							exit
						}
						*Check if T set is null in bootstrap resamples (this won't trigger an error, but estimands will be flagged for ignoring later)
						if r(N) > 0 & `b' > 0 {
								local therewerewarnings=1
								di as result "-----Warning: The estimated identified set for tau0 is null in bootstrap sample `b'. This occurs when s(y)*(1-tau0) integrates to less than 1-tau0L (see paper for definitions), making the lower bound for tau0 larger than the upper bound. This bootstrap sample will be ignored for fuzzy bounds."
								qui replace keepthisbootstrap = 0 if Pi < 1-tau0L & tau > 0
									*Note keepthisbootstrap will later be switched to 1 for all observations with FRD=0 or naive=1
						}
						
						qui gen tau0U = max(0,tau-((1-tau)*(treatRight-treatLeft))/(1-treatRight))
						qui count if Pi < .99*(1-tau0U) & tau>0
						if r(N) > 0 & `b'==0 {
							local therewerewarnings=1
							di as result "Warning: The integral of stilde(y)=s(y)*(1-tau0) is small enough to shrink the identified set of (tau0, tau1) pairs (see paper for definitions)."							
							if `verbose'==1 { 		
								di as result "--------The tables below indicate the total mass of s(y)*(1-tau0), and the value of tau0U, respectively, (for the main sample, not a bootstrap resample). w=0 indicates the full sample (not conditional on covariate), and fixedtau=0 indicates tau is estimated, not fixed via {cmd:potential_taus}:"
								di as result "mass of s(y)*(1-tau0):"
								table fixedtau w, c(mean Pi)
								di as result "value of tau0U:"
								table fixedtau w, c(mean tau0U)
							}	
						}
						*Detect whether data is consistent with refinement B:
						qui count if Pi < 1 & tau>0 & tau <= ((1-tau)*(treatRight-treatLeft))/(1-treatRight) & fixedtau==0 & w==0
						if r(N) > 0 & `b'==0 & `dorefB'==1{
							local therewerewarnings=1
							di as result "Warning: Bounds for refinement B impose tau0=0 even though the integral of stilde(y)=s(y)*(1-tau0) is less than unity, suggesting that tau0>0."
						}
						
						qui count if tau>0 & tau > ((1-tau)*(treatRight-treatLeft))/(1-treatRight) & fixedtau==0 & w==0
						if r(N) > 0 & `b'==0 & `dorefB'==1{
							local therewerewarnings=1
							di as result "Warning: Refinement B will not be computed because tau/(1-tau) > (treatRight-treatLeft))/(1-treatRight), which implies that tau0 > 0."
							local dorefB = 0
						}
						
						qui count if Pi > 1.01*(1-tau0U) & tau>0
						if r(N) > 0 & `b'==0 {
							local therewerewarnings=1
							di as result "Warning: The integral of stilde(y)=s(y)*(1-tau0) is larger than one minus the smallest possible value of tau0. This could be evidence against the model."
							
							if `verbose'==1 { 		
								di as result "--------The tables below indicate the total mass of s(y)*(1-tau0), and the value of tau0U, respectively, (for the main sample, not a bootstrap resample). w=0 indicates the full sample (not conditional on covariate), andfixedtau=0 indicates tau is estimated, not fixed via {cmd:potential_taus}:"
								di as result "mass of s(y)*(1-tau0):"
								table fixedtau w, c(mean Pi)
								di as result "value of tau0U:"
								table fixedtau w, c(mean tau0U)
							}	
						}
						
						if `therewereerrors'==1{
							exit
						}
						
						drop tau0U tau0L
												
						if `discretey'==1 {
							sort fixedtau w b little_y
							qui by fixedtau w b: gen F = sum(s_notau)
						}
						else {
							sort fixedtau w b little_y
							qui by fixedtau w b: gen lengths = little_y[_n] - little_y[_n-1]
							qui replace lengths = 0 if _n==1
							qui gen temp = lengths*s_notau
							drop lengths
							qui bysort fixedtau w b: gen F = sum(temp)
							drop temp
						}
					}
		}
	
	if `b'==0 & `therewerewarnings'==1 & `verbose'==0 {
		di as result "There were warnings. For more details, please set the {cmd:verbose} command."
	}

	if `b'==0 & `saveToFileInputs'==1 {
		drop keepthisbootstrap
		qui save "`outputfileInputs'", replace
		gen keepthisbootstrap = 1
	}
	
	*Expand the dataset across treatment effects and values of t 
	*-----------------------------------------------------------------------
	sort fixedtau w b little_y
	
	qui expand `num_tau_pairs'+3
	qui bysort fixedtau w b little_y: gen t = _n-3 if _n>3
	qui by fixedtau w b little_y: gen naive = _n==1 | _n==3
	qui by fixedtau w b little_y: gen FRD = _n >= 3
	
	qui replace keepthisbootstrap = 1 if !(naive==0 & FRD==1)
	
	*Don't need naive CDFs for covariate refinement or fixed tau estimates
	qui drop if naive==1 & (w > 0 | fixedtau > 0)
	*Don't need covariate refinement for fixed tau
	qui drop if w > 0 & fixedtau > 0
	
	order fixedtau w b FRD naive t little_y
	gsort fixedtau w b FRD -naive t little_y
	
	*For point identified treatment effects (e.g. naive estimates, we'll store result in F_Y0_lower and leave F_Y1_upper missing
	qui gen F_Y0_lower = .
	qui gen F_Y0_upper = .
	qui gen F_Y1_lower = .
	qui gen F_Y1_upper = .
	
		*Sharp research design CDFs:
		*-----------------------------------------------------------------------
		*Naive effects ignoring manipulation
			qui replace F_Y0_lower  = F_left if naive==1 & FRD==0
			qui replace F_Y1_lower = F_right if naive==1 & FRD==0
	
		*Bounds accounting for manipulation
			qui replace F_Y0_lower = F_left if naive==0 & FRD==0
			qui replace F_Y0_upper = F_left if naive==0 & FRD==0
		
			*Y1 lower: trim tau from the top end of F_right, by covariate value (since we want to compute SRD treatment effect for each covariate value):
			qui gen qtiletemp = little_y if F_right >= 1 - tau & naive==0 & FRD==0
			qui egen qtile = min(qtile) if naive==0 & FRD==0, by(fixedtau w)
			qui replace F_Y1_lower = 1/(1-tau)*((little_y < qtile)*F_right+(little_y>=qtile)*(1-tau)) if naive==0 & FRD==0
			drop qtiletemp qtile
		
			*Y1 upper: trim tau_hat from the bottom end of F_left, by covariate value:
			qui gen qtiletemp = little_y if F_right >= tau & naive==0 & FRD==0
			qui egen qtile = min(qtile) if naive==0 & FRD==0, by(fixedtau w)
			qui replace F_Y1_upper = 1/(1-tau)*((little_y > qtile)*(F_right-tau)+(little_y==qtile)*(F_right-tau)) if naive==0 & FRD==0
			drop qtiletemp qtile
			
		*Fuzzy research design CDFs:
		*-----------------------------------------------------------------------
		if `fuzzy'==1{
		
			*Naive effects ignoring manipulation
			*-------------------------------------------------------------------	
				qui gen kappa0_naive = (1-treatRight)/(1-treatLeft)
				if `norightuntreated'==0 {
					qui replace F_Y0_lower = (F_left_untreated-kappa0_naive*F_right_untreated)/(1-kappa0_naive) if naive==1 & FRD==1
				}
				else{
					qui replace F_Y0_lower = F_left_untreated if naive==1 & FRD==1
				}
				
				qui gen kappa1_naive = treatLeft/treatRight
				if `nolefttreated'==0 {
					qui replace F_Y1_lower = (F_right_treated-kappa1_naive*F_left_treated)/(1-kappa1_naive) if naive==1 & FRD==1
				}
				else{
					qui replace F_Y1_lower = F_right_treated if naive==1 & FRD==1
				}
				
				drop kappa0_naive kappa1_naive
		
			*Bounds accounting for manipulation
			*-------------------------------------------------------------------	
				*Get grid of tau0s and tau1s corresponding to each t			
				qui gen tau1L = max(0, 1-(1-tau)/treatRight)
				qui gen tau1U = min(1-(1-tau)*treatLeft/treatRight, tau/treatRight)
				qui gen tau0L = min(1, tau/(1-treatRight))
				qui gen tau0U = max(0,tau-((1-tau)*(treatRight-treatLeft))/(1-treatRight))
				
				qui gen boundshrinks = (tau0U < 1-Pi) & (tau > 0)
					qui replace tau0U = 1-Pi if boundshrinks==1
					qui replace tau1U = (tau-tau0U*(1-treatRight))/treatRight if boundshrinks==1
				qui drop boundshrinks
				
				qui gen tau1 = tau1L+(tau1U-tau1L)/(`num_tau_pairs'-1)*(t-1) if t!=.
				qui gen tau0 = tau0L+(tau0U-tau0L)/(`num_tau_pairs'-1)*(t-1) if t!=.
				order fixedtau w b t tau1 tau0 little_y
				
				if `num_tau_pairs'==1 {
					*if num_tau_pairs = 1, make the single tau the "rightmost" extreme of the t-set, such that user can enforce the always treated only assumption if it is consistent with data
					qui replace tau1 = tau1U
					qui replace tau0 = tau0U
				}
				
				drop tau1L tau1U tau0L tau0U
				
				*Y1 CDFs
				*-----------------------------------------------------------
				qui gen toTrim=tau1/(1-kappa1)				
				qui summ toTrim
				if `r(max)'>=1 {
					di as result "Warning: the admissible value of tau1/(1-kappa1)=`r(max)' results in zero potentially assigned compliers on the right of the cutoff. This possible value of tau1 will be ignored."
					*disp "b: `b'"
					*table tau1 w if toTrim>=1
				}
				
				*Y1 lower: trim tau1/(1-kappa1) from the top end of G:
					qui gen qtiletemp = little_y if G >= 1 - toTrim & naive==0 & FRD==1
					qui egen qtile = min(qtiletemp) if naive==0 & FRD==1, by(fixedtau t w)
					qui replace F_Y1_lower = 1/(1-toTrim)*((little_y < qtile)*G+(little_y>=qtile)*(1-toTrim)) if naive==0 & FRD==1
					drop qtiletemp qtile 
							
				*Y1 upper: trim tau1/(1-kappa1) from the bottom end of G:
					qui gen qtiletemp = little_y if G >= toTrim & naive==0 & FRD==1
					qui egen qtile = min(qtiletemp) if naive==0 & FRD==1, by(fixedtau t w)
					qui replace F_Y1_upper = 1/(1-toTrim)*((little_y > qtile)*(G-toTrim)+(little_y==qtile)*(G-toTrim)) if naive==0 & FRD==1
					drop qtiletemp qtile 
					
				*replace F_Y1_lower = . if naive==0 & FRD==1 & toTrim >= 1
				*replace F_Y1_upper= . if naive==0 & FRD==1 & toTrim >= 1
				
				drop toTrim
				
				*Y0 CDFs:
				*-----------------------------------------------------------
				if `norightuntreated' == 1 {
					*If there are no untreated units on the right, then kappa0 = 0. There are no never-takers and the distribution of Y0 is point identified by untreated units on left
					qui replace F_Y0_lower = F_left_untreated if naive==0 & FRD==1
					qui replace F_Y0_upper = F_left_untreated if naive==0 & FRD==1
				}
				else {								
					qui gen Fs=F/(1-tau0)
					*Fs will have a mass of Pi/1-tau0 \ge 1
					
					*For F_N0_L, want to find point QL at which there is mass one to the left of QL: the one "quantile" of Fs
					qui gen qtiletemp = little_y if Fs >= 1 & naive==0 & FRD==1
					qui egen qtile = min(qtiletemp) if naive==0 & FRD==1, by(fixedtau t w)
					qui gen F_N0_FRD_lower = (little_y < qtile)*Fs+(little_y>=qtile)*1 if naive==0 & FRD==1
					drop qtiletemp qtile 
						
					*For F_N0_U, want to find point QU at which there is mass one to the right of QU: the Pi/1-tau0 - 1 "quantile" of Fs
					qui gen toTrim = Pi / (1-tau0) - 1 if naive==0 & FRD==1
					qui gen qtiletemp = little_y if Fs >= toTrim & naive==0 & FRD==1
					qui egen qtile = min(qtiletemp) if naive==0 & FRD==1, by(fixedtau t w)
					qui gen F_N0_FRD_upper = (little_y >= qtile)*(Fs-toTrim) if naive==0 & FRD==1
					drop qtiletemp qtile toTrim
					
					drop F Fs				
					
					qui replace F_N0_FRD_lower = F_right_untreated if naive==0 & FRD==1 & abs(tau0)<0.00000001
					qui replace F_N0_FRD_upper = F_right_untreated if naive==0 & FRD==1 & abs(tau0)<0.00000001
					
					qui gen N0_weight = kappa0*(1-tau0)
					qui replace F_Y0_lower = (F_left_untreated - N0_weight*F_N0_FRD_upper)/(1-N0_weight) if naive==0 & FRD==1 & tau0!=1
					qui replace F_Y0_upper = (F_left_untreated - N0_weight*F_N0_FRD_lower)/(1-N0_weight) if naive==0 & FRD==1 & tau0!=1
					
					drop F_N0_FRD_lower F_N0_FRD_upper N0_weight kappa0
					
					qui replace F_Y0_lower = F_left_untreated if naive==0 & FRD==1 & abs(1-tau0)<0.00000001
					qui replace F_Y0_upper = F_left_untreated if naive==0 & FRD==1 & abs(1-tau0)<0.00000001
					
			}
		}
		else {
			qui gen G = .
			qui gen s_notau = .
			qui gen tau0 = .
			qui gen tau1 = .
		}

		*CDFs for units on the right side of cutoff (doing this before covariate CDFs so that we don't need to save F_left, etc in the covariate collapse
		*-----------------------------------------------------------
		qui gen isrightside=0
		order fixedtau isrightside
		if `righteffects'==1 {
			*Get a copy of all standard bounds obervations (not naive, not covariate, and not fixed tau)
			qui expand 2 if fixedtau==0 & w==0 & naive==0
			qui bysort fixedtau w b little_y naive t: replace isrightside=1 if _n==2 & fixedtau==0 & w==0 & naive==0
			
			*SRD
			qui replace F_Y1_lower = F_right if isrightside==1 & FRD==0
			qui replace F_Y1_upper = F_right if isrightside==1 & FRD==0
			qui replace F_Y0_lower = (1-tau)*F_left+tau*(little_y >=scalar(yextreme_min)) if isrightside==1 & FRD==0
			qui replace F_Y0_upper = (1-tau)*F_left+tau*(little_y >=scalar(yextreme_max)) if isrightside==1 & FRD==0
			
			if `fuzzy'==1{
			
				*Expand by `num_lambdas' to get (t,lambda) pairs:
				qui expand `num_lambdas' if isrightside==1 & FRD==1
				qui bysort fixedtau w isrightside b little_y naive t: gen lambda=(_n-1)/(`num_lambdas'-1) if isrightside==1 & FRD==1
				order fixedtau w isrightside b t lambda little_y
				sort fixedtau w isrightside b t lambda little_y
			
				*Y1 CDFs
				*------------
				qui gen toTrim=tau1*lambda/(1-kappa1)
				
				*Y1 lower: trim lambda*tau1/(1-kappa1) from the top end of G:
					qui gen qtiletemp = little_y if G >= 1 - toTrim & isrightside==1 & FRD==1
					*(note fixedtau=0 and w=0 if isrightside==1)
					qui egen qtile = min(qtiletemp) if isrightside==1 & FRD==1, by(t lambda)
					qui replace F_Y1_lower = 1/(1-toTrim)*((little_y < qtile)*G+(little_y>=qtile)*(1-toTrim)) if isrightside==1 & FRD==1
					drop qtiletemp qtile 
				
				*Y1 upper: trim lambda*tau1/(1-kappa1) from the bottom end of G:
					qui gen qtiletemp = little_y if G >= toTrim & isrightside==1 & FRD==1
					qui egen qtile = min(qtiletemp) if isrightside==1 & FRD==1, by(t lambda)
					qui replace F_Y1_upper = 1/(1-toTrim)*((little_y > qtile)*(G-toTrim)+(little_y==qtile)*(G-toTrim)) if isrightside==1 & FRD==1
					drop qtiletemp qtile 
				
				drop toTrim
			
				*Y0 CDFs
				*------------
				qui gen staulambda = ((1-tau1)*treatRight-(1-tau)*treatLeft)/((1-lambda*tau1)*treatRight-(1-tau)*treatLeft)
				qui replace F_Y0_lower = staulambda*F_Y0_lower+(1-staulambda)*(little_y>=scalar(yextreme_min)) if isrightside==1 & FRD==1
				qui replace F_Y0_upper = staulambda*F_Y0_upper+(1-staulambda)*(little_y>=scalar(yextreme_max)) if isrightside==1 & FRD==1
				drop staulambda
				
				sort fixedtau w isrightside b little_y t lambda
				qui bysort fixedtau w isrightside b little_y: replace t = _n if isrightside==1 & FRD==1
				sort fixedtau b w isrightside little_y t
			}
			
			*Detect whether data is consistent with refinement C now, before we drop lambda:
			if `fuzzy'==1 {
				if `b'==0 {
					if `dorefC'==1 {
						qui summ tau0 if lambda==0 & FRD==1 & naive==0 & fixedtau==0 & w==0 & isrightside==1 & b==0
						if r(min) > 0.00000001 {
							di as result "Warning: The data are inconsistent with tau0=0 & lambda=0 so refinement C will not be computed"
							local dorefC = 0
						}
						else {
							qui summ t if abs(tau0)<0.00000001 & lambda==0 & FRD==1 & naive==0 & fixedtau==0 & w==0 & isrightside==1 & b==0
							if r(sd) == 0{
								local refCt = r(mean)
							}
							else{
								di as result "Warning: there was a problem detecting treatment effect estimates compatible with refinement C, so it will not be computed"
								local dorefC = 0
							}
						}
					}
				}
				drop lambda
			}
		}
	
		drop dens_right_untreated dens_left_untreated F_right F_left F_right_treated F_right_untreated F_left_treated F_left_untreated G s_notau

		*Covariate refinement
		*-----------------------------------------------------------		
		if `gotcovs'==1 {		
			qui gen isw = (w>0)
			order fixedtau isw isrightside
				
			*SRD
			qui gen pw=.
			forvalues w = 1/`ncovlevels' {
				*Full sample treatment probabilities are stored in the first column of matrix probw, bootstrap sample b is stored in column b+1
				if `b' < 1 {
					qui replace pw = probw[`w',1] if w==`w'
				}
				if `b' >= 1 {
					qui replace pw = probw[`w',`b'+1] if w==`w'
				}
					
			}
			
			*For Y0 ITT/SRD, the upper and lower bound under the covariate refinement is the same as the normal Y0 bounds, which is given by F_left 
			*on the full sample. This will generally have better finite sample performance than estimating things conditional on w and aggregating, 
			*so we will override the below later on (here we calculate it in the same manner as Y1 in case we detect an unexpected problem later)
				qui replace F_Y0_lower = F_Y0_lower*pw if isw==1 & FRD==0
				qui replace F_Y0_upper = F_Y0_upper*pw if isw==1 & FRD==0
			
			qui replace F_Y1_lower = F_Y1_lower*pw if isw==1 & FRD==0
			qui replace F_Y1_upper = F_Y1_upper*pw if isw==1 & FRD==0		
			
			*FRD			
			if `fuzzy'==1{
				qui gen pi_minus=(1-tau1)/(1-tau)*treatRight-treatLeft
				qui summ pi_minus
				if r(min) < 0 {
					di as result "Warning: estimated proportion of potentially-assigned compliers at the cutoff conditional on covariate value is negative for some values of the covariates, replacing with zero"
					qui replace pi_minus=0 if pi_minus < 0
				}
				if r(max) > 1 {
					di as result "Warning: estimated proportion of potentially-assigned compliers at the cutoff conditional on covariate value is greater than unity for some values of the covariates, replacing with one"
					qui replace pi_minus=1 if pi_minus > 1
				}
				
				qui gen t_w = t if isw==1 & FRD==1
				sort fixedtau w b t little_y
				
				*These next two loops are a fancy way of creating all combinations (t_1,t_2...t_k) (where k is ncovslevels) such that
				*for a given row with covariate value `w', the treatment effect is copmuted with the t-value t_`w'. Then we'll collapse over w.
				
				*First, expand the dataset by num_tau_pairs k-1 times and create indices ttemp_1..ttemp_(k-1)
				forvalues i = 2/`ncovlevels' {
					qui gen cur_unique_index = _n if isw==1 & FRD==1
					qui expand `num_tau_pairs' if isw==1 & FRD==1
					local w = `i'-1
					qui bysort cur_unique_index: gen ttemp_`w' = _n if isw==1 & FRD==1
					qui drop cur_unique_index
				}
				
				*Now generate t_1, t_2, such that t_`w' equals the original t (renamed as "t_w") and the other k-1 indices get ttemp1..ttemp_(k-1)
				forvalues w = 1/`ncovlevels' {
					qui gen t_`w' = t_w if w == `w'
					if `w' < `ncovlevels' {
						qui replace t_`w' = ttemp_`w' if `w' < w
					}
					if `w' > 1 {
						local wminus1 = `w'-1
						qui replace t_`w' = ttemp_`wminus1' if `w' > w
					}
				}
				drop t_w ttemp_*
				order t t_*
				
				*Now, with our num_tau_pairs^k copies of each observation by w, let's assign an overall group number "t" to the k rows sharing the same combination
				*Based on the presort below, the rownumber of an observation within a (little_y x w) cell (fixedtau=0 and b is fixed here) maps to the same (t_1..t_k) combo across w's
				sort fixedtau b little_y w t_*
				qui bysort fixedtau b little_y w: replace t = _n if isw==1 & FRD==1
				sort fixedtau b little_y t w
				
				qui replace F_Y0_lower = F_Y0_lower*pw*pi_minus if isw==1 & FRD==1
				qui replace F_Y0_upper = F_Y0_upper*pw*pi_minus if isw==1 & FRD==1
				qui replace F_Y1_lower = F_Y1_lower*pw*pi_minus if isw==1 & FRD==1
				qui replace F_Y1_upper = F_Y1_upper*pw*pi_minus if isw==1 & FRD==1
				qui gen omega_denom = pw*pi_minus if isw==1 & FRD==1
				
			}
			else {
				qui replace F_Y0_lower = . if isw==1 & FRD==1
				qui replace F_Y0_upper = . if isw==1 & FRD==1
				qui replace F_Y1_lower = . if isw==1 & FRD==1
				qui replace F_Y1_upper = . if isw==1 & FRD==1
				*making a dummy omega_denom just to share the below code whether or not fuzzy is set
				qui gen omega_denom=1
			}
				
			*Collapse over covariate values for both sharp and fuzzy covariate-refined bounds, then append back to the rest of the dataset
				qui save "`tmpfile'", replace
					qui keep if isw==0
				qui save "`tmpfile2'", replace
				
				qui use "`tmpfile'", clear
				qui keep if isw==1
			
				qui egen F_Y0_lower_sum=sum(F_Y0_lower),  by(fixedtau isw isrightside b FRD naive t little_y)
				qui egen F_Y0_upper_sum=sum(F_Y0_upper),  by(fixedtau isw isrightside b FRD naive t little_y)
				qui egen F_Y1_lower_sum=sum(F_Y1_lower),  by(fixedtau isw isrightside b FRD naive t little_y)
				qui egen F_Y1_upper_sum=sum(F_Y1_upper),  by(fixedtau isw isrightside b FRD naive t little_y)
				qui egen omega_denom_sum=sum(omega_denom),  by(fixedtau isw isrightside b FRD naive t little_y)	
				qui bysort fixedtau isw isrightside b FRD naive t little_y: keep if _n==1
				drop F_Y0_lower F_Y0_upper F_Y1_lower F_Y1_upper omega_denom
				rename (F_Y0_lower_sum F_Y0_upper_sum F_Y1_lower_sum F_Y1_upper_sum omega_denom_sum) (F_Y0_lower F_Y0_upper F_Y1_lower F_Y1_upper omega_denom)
			
				qui append using "`tmpfile2'"
			
			*Finish up fuzzy bounds
			if `fuzzy'==1{
				qui replace F_Y0_lower = F_Y0_lower/omega_denom if isw==1 & FRD==1
				qui replace F_Y0_upper = F_Y0_upper/omega_denom if isw==1 & FRD==1
				qui replace F_Y1_lower = F_Y1_lower/omega_denom if isw==1 & FRD==1
				qui replace F_Y1_upper = F_Y1_upper/omega_denom if isw==1 & FRD==1
				drop omega_denom
			}
		}
		else {
			drop w
			qui gen isw = 0
		}
	order fixedtau isw isrightside b FRD naive t little_y tau0 tau1
	gsort fixedtau isw isrightside b FRD -naive t little_y
		
	*Compute treatment effects from CDFs
	********************************************************************************
	
	foreach p of numlist `percentiles' {
		*Compute the functional corresponding to `p' for each of the four potential outcome CDFs:
		foreach F of varlist F_Y0_lower F_Y0_upper F_Y1_lower F_Y1_upper {
			if `p'==-1 {
				*Taking right Reimann sum to compute expectation, which is exact if y is discretely distributed. Treating cdfy[1] as the total mass at or to the left of little_ys[1].
				sort fixedtau isw isrightside b FRD naive t little_y
				qui bysort fixedtau isw isrightside b FRD naive t: gen mass_in_interval = `F'[_n] - `F'[_n-1]
				qui bysort fixedtau isw isrightside b FRD naive t: replace mass_in_interval = `F' if _n==1
				qui gen E_`F' = little_y*mass_in_interval
				drop mass_in_interval
			}
			else {
				qui gen Q`p'_`F' = little_y if `F' >= `p'/100 & `F' != . 
			}
		}			
	}
	
	qui egen groupid = group(fixedtau isw isrightside b FRD naive t),missing
	sort groupid
	qui save "`tmpfile'", replace
		collapse (first) fixedtau isw isrightside b FRD naive t tau1 tau0 tau treatLeft treatRight keepthisbootstrap, by(groupid)
		sort groupid
		qui save "`tmpfile2'", replace	
	qui use "`tmpfile'", clear
	
	mata: id = st_data(., ("groupid"))
		
	local num_p : word count `percentiles'
	forval i=1/`num_p' {
		local p : word `i' of `percentiles'
		
		qui use "`tmpfile'", clear	
		sort fixedtau isw isrightside b FRD naive t
	
		if `p'==-1{
			mata: M = st_data(., ("E_F_Y0_lower","E_F_Y0_upper","E_F_Y1_lower","E_F_Y1_upper"))
			mata: C=_mm_collapse(M, 1, id, &my_sum())
			mata: st_matrix("collapse_results", C[.,(2,3,4,5)])
			matrix colnames collapse_results="Y0_lower" "Y0_upper" "Y1_lower" "Y1_upper"
		} 
		else {
			mata: M = st_data(., ("Q`p'_F_Y0_lower","Q`p'_F_Y0_upper","Q`p'_F_Y1_lower","Q`p'_F_Y1_upper"))
			mata: C=_mm_collapse(M, 1, id, &my_min())
			mata: st_matrix("collapse_results", C[.,(2,3,4,5)])			
		}
		matrix colnames collapse_results="Y0_lower" "Y0_upper" "Y1_lower" "Y1_upper"	
		
		qui use "`tmpfile2'", clear
		sort groupid
		qui svmat collapse_results, names(col)
		qui gen p = `p'
		
		if `i'>1 {
			qui append using "`tmpfileoutcomes_b'"
		}
		qui save "`tmpfileoutcomes_b'", replace

	}
	drop groupid
		
	qui gen TE_lower = Y1_lower-Y0_upper if naive==0
	qui gen TE_upper = Y1_upper-Y0_lower if naive==0
	qui replace TE_lower = Y1_lower-Y0_lower if naive==1
	
	*When all are missing, my_sum returns 0. So if fuzzy=FALSE we should replace the 0 treatment effects by .
	if `fuzzy'==0 {
		qui replace TE_lower = . if FRD==1 & p==-1
		qui replace TE_upper = . if FRD==1 & p==-1
	}
	
	if `gotcovs'==1 {
		*Override the estimated Y0 bounds for ITT/SRD covariate refinement with the more precisely estimated bounds from the full sample,
		*which identify the same population parameter and are simply F_left. This line could be removed for simplicity with minimal effect
		foreach p of numlist `percentiles' {
			qui summ Y0_lower if isw==0 & FRD==0 & naive==0 & fixedtau==0 & isrightside==0 & p==`p'
			local Y0 = `r(mean)'
			if `r(N)'==1 {
				qui replace TE_upper = Y1_upper-`Y0' if isw==1 & FRD==0 & naive==0 & fixedtau==0 & isrightside==0 & p==`p'
				qui replace TE_lower = Y1_lower-`Y0' if isw==1 & FRD==0 & naive==0 & fixedtau==0 & isrightside==0 & p==`p'
			}
			else{
				di as result "Warning: expected to find a single observation with the below filters, in order to finalize Y0 estimate for the ITT/SRD covariate refinement. This could be a sign of an unanticipated problem."
				summ Y0_lower if isw==0 & FRD==0 & naive==0 & fixedtau==0 & isrightside==0 & p==`p'
			}
		}
	}
	
	sort fixedtau isw isrightside b FRD naive t p
	order fixedtau isw isrightside b FRD naive t tau0 tau1 p TE_lower TE_upper
	
	if `fuzzy'==1 {
		gen refA = tau1 >= tau
	}
	else {
		*Make a variable called refA later even if SRD, so that we can assume its existence
		gen refA = 0
	}
	
if `b'>-1 & `num_bs' > 0 {
	qui append using "`tmpfileTEs'"
}
qui save "`tmpfileTEs'", replace

*End bootstrap loop
********************************************************************************
********************************************************************************
}

*drop observations corresponding to fuzzy bounds in bootstraps for which the T set was null
qui drop if keepthisbootstrap==0
drop keepthisbootstrap

gsort fixedtau isw isrightside b p FRD -naive t

*Begin saving results
********************************************************************************
ereturn clear
ereturn scalar tau_hat = scalar(tau_hat_pt)
ereturn scalar takeup_increase = scalar(takeup_increase_pt)

*Confidence intervals for tau_hat and takeup_increase
ereturn scalar tau_hat_CI_lower = scalar(tau_star) - invnormal(1-scalar(alpha)/2)*scalar(sdtauhatnudged)
ereturn scalar tau_hat_CI_upper = scalar(tau_star) + invnormal(1-scalar(alpha)/2)*scalar(sdtauhatnudged)

ereturn scalar takeup_CI_lower = scalar(takeup_increase_pt) - invnormal(1-scalar(alpha)/2)*scalar(sdtakeup)
ereturn scalar takeup_CI_upper = scalar(takeup_increase_pt) + invnormal(1-scalar(alpha)/2)*scalar(sdtakeup)

*Save to file if requested
if `saveToFile'==1 {
	qui save "`outputfile'", replace
}

*Collect point estimates and compute confidence intervals for treatment effects
********************************************************************************

	di as result "Gathering point estimates and confidence intervals"
	
	*Collapse over bootstraps to extract standard deviations
	qui replace b = 1 if b > 1
	collapse (count) cnt=TE_lower (first) gammaL=TE_lower gammaU=TE_upper refA=refA (sd) sL=TE_lower sU=TE_upper, by(FRD naive fixedtau isw isrightside t p b)

	scalar minr = invnormal(1-scalar(alpha))
	scalar maxr = invnormal(1-scalar(alpha)/2)
	qui gen minr = .
	qui gen thisobjfn = .
	qui gen minobjfn = .
	local rwarning=0
	forvalues i = 1/`num_rs' {
		scalar r = scalar(minr)+(scalar(maxr)-scalar(minr))/(`num_rs'-1)*`i'
		qui replace thisobjfn = abs(normal(scalar(r)+(gammaU-gammaL)/max(sU, sL))-normal(-scalar(r))-(1-scalar(alpha))) if max(sU,sL) > 0 & b==1
		qui replace minr = scalar(r) if (thisobjfn <= minobjfn | `i'==1) & b==1
		qui replace minobjfn=thisobjfn if (thisobjfn <= minobjfn | `i'==1) & b==1
		qui count if (thisobjfn == minobjfn) & b==1
		if r(N)>0 & `rwarning==0' {
			di as result "Warning: r_alpha (see inference section of paper for definition) may not be unique: {cmd: rdbounds} will take the highest value minimizing objective function between invnormal(1-scalar(alpha)) and invnormal(1-scalar(alpha/2))"
			local rwarning=1
		}
	}
	drop thisobjfn minobjfn
	rename minr r
	qui replace r = invnormal(1-scalar(alpha)/2) if naive==1
	
	sort fixedtau isw isrightside naive FRD t p b
	qui by fixedtau isw isrightside naive FRD t p: gen CI_lower = gammaL[_n-2]-r[_n]*sL[_n] if b==1
	qui by fixedtau isw isrightside naive FRD t p: gen CI_upper = gammaU[_n-2]+r[_n]*sU[_n] if b==1
	qui by fixedtau isw isrightside naive FRD t p: replace CI_upper = gammaL[_n-2]+r[_n]*sL[_n] if b==1 & naive==1	
	qui drop if b==-1
		
	foreach p of numlist `percentiles' {
		
		if `p'==-1{
			local stub="ATE"
		}
		else{
			local stub = "QTE`p'"
		}
		matrix treatment_effects_`stub' = J(10, 4, .)
		matrix rownames treatment_effects_`stub'= "SRD naive" "SRD bounds" "FRD naive" "FRD bounds" "FRD refA" "FRD refB" "SRD cov" "FRD cov" "SRD right" "FRD right"
		matrix colnames treatment_effects_`stub'= "Lower/point" "Upper bound" "CI lower" "CI upper"
		
		*Get point estimates
		*-----------------------------------------------------------------------
			*SRD
			*---------------
			qui summ gammaL if FRD==0 & naive==1 & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[1,1] = r(mean)
			
			qui summ gammaL if FRD==0 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[2,1] = r(mean)
			qui summ gammaU if FRD==0 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[2,2] = r(mean)
			
			qui summ gammaL if FRD==0 & naive==0 & fixedtau==0 & isw==1 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[7,1] = r(mean)
			qui summ gammaU if FRD==0 & naive==0 & fixedtau==0 & isw==1 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[7,2] = r(mean)
			
			qui summ gammaL if FRD==0 & naive==0 & fixedtau==0 & isw==0 & isrightside==1 & b==0 & p==`p'
			matrix treatment_effects_`stub'[9,1] = r(mean)
			qui summ gammaU if FRD==0 & naive==0 & fixedtau==0 & isw==0 & isrightside==1 & b==0 & p==`p'
			matrix treatment_effects_`stub'[9,2] = r(mean)
			
			*FRD
			*---------------
			qui summ gammaL if FRD==1 & naive==1 & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[3,1] = r(mean)
			
			qui summ gammaL if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[4,1] = r(min)
			qui summ gammaU if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[4,2] = r(max)
			
			if `dorefA'==1 {
				qui summ gammaL if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p' & refA==1
				matrix treatment_effects_`stub'[5,1] = r(min)
				qui summ gammaU if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p' & refA==1
				matrix treatment_effects_`stub'[5,2] = r(max)
			}
			
			if `dorefB'==1 {
				qui summ gammaL if FRD==1 & naive==0 & t==`num_tau_pairs' & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p'
				matrix treatment_effects_`stub'[6,1] = r(min)			
				qui summ gammaU if FRD==1 & naive==0 & t==`num_tau_pairs' & fixedtau==0 & isw==0 & isrightside==0 & b==0 & p==`p'
				matrix treatment_effects_`stub'[6,2] = r(max)			
			}		
		
			qui summ gammaL if FRD==1 & naive==0 & fixedtau==0 & isw==1 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[8,1] = r(min)
			qui summ gammaU if FRD==1 & naive==0 & fixedtau==0 & isw==1 & isrightside==0 & b==0 & p==`p'
			matrix treatment_effects_`stub'[8,2] = r(max)
			
			qui summ gammaL if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==1 & b==0 & p==`p'
			matrix treatment_effects_`stub'[10,1] = r(min)
			qui summ gammaU if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==1 & b==0 & p==`p'
			matrix treatment_effects_`stub'[10,2] = r(max)
			
			*One could use this code to report refinement C
			*if `dorefC'==1 {
			*	qui summ gammaL if FRD==1 & naive==0 & t==`refCt' & fixedtau==0 & isw==0 & isrightside==1 & b==0 & p==`p'
			*	ereturn scalar refC_lower =  = r(min)			
			*	qui summ gammaU if FRD==1 & naive==0 & t==`refCt' & fixedtau==0 & isw==0 & isrightside==1 & b==0 & p==`p'
			*	ereturn scalar refC_upper =  = r(max)			
			*}
		
		
		*Compute confidence intervals for treatment effects
		*-----------------------------------------------------------------------
			*SRD
			*---------------
			qui summ CI_lower if FRD==0 & naive==1 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[1,3] = r(mean)
			qui summ CI_upper if FRD==0 & naive==1 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[1,4] = r(mean)
			
			qui summ CI_lower if FRD==0 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[2,3] = r(mean)
			qui summ CI_upper if FRD==0 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[2,4] = r(mean)
			
			qui summ CI_lower if FRD==0 & naive==0 & fixedtau==0 & isw==1 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[7,3] = r(mean)
			qui summ CI_upper if FRD==0 & naive==0 & fixedtau==0 & isw==1 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[7,4] = r(mean)
			
			qui summ CI_lower if FRD==0 & naive==0 & fixedtau==0 & isw==0 & isrightside==1 & b==1 & p==`p'
			matrix treatment_effects_`stub'[9,3] = r(mean)
			qui summ CI_upper if FRD==0 & naive==0 & fixedtau==0 & isw==0 & isrightside==1 & b==1 & p==`p'
			matrix treatment_effects_`stub'[9,4] = r(mean)
			
			*FRD
			*---------------
			qui summ CI_lower if FRD==1 & naive==1 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[3,3] = r(mean)
			qui summ CI_upper if FRD==1 & naive==1 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[3,4] = r(mean)
			
			qui summ CI_lower if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[4,3] = r(min)
			qui summ CI_upper if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[4,4] = r(max)
			
			if `dorefA'==1 {
				qui summ CI_lower if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p' & refA==1
				matrix treatment_effects_`stub'[5,3] = r(min)
				qui summ CI_upper if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p' & refA==1
				matrix treatment_effects_`stub'[5,4] = r(max)
			}
			
			if `dorefB'==1 {
				qui summ CI_lower if FRD==1 & naive==0 & t==`num_tau_pairs' & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
				matrix treatment_effects_`stub'[6,3] = r(min)
				qui summ CI_upper if FRD==1 & naive==0 & t==`num_tau_pairs' & fixedtau==0 & isw==0 & isrightside==0 & b==1 & p==`p'
				matrix treatment_effects_`stub'[6,4] = r(max)			
			}
			
			qui summ CI_lower if FRD==1 & naive==0 & fixedtau==0 & isw==1 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[8,3] = r(min)
			qui summ CI_upper if FRD==1 & naive==0 & fixedtau==0 & isw==1 & isrightside==0 & b==1 & p==`p'
			matrix treatment_effects_`stub'[8,4] = r(max)
			
			qui summ CI_lower if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==1 & b==1 & p==`p'
			matrix treatment_effects_`stub'[10,3] = r(min)
			qui summ CI_upper if FRD==1 & naive==0 & fixedtau==0 & isw==0 & isrightside==1 & b==1 & p==`p'
			matrix treatment_effects_`stub'[10,4] = r(max)
			
			*One could use this code to report refinement C
				*if `dorefC'==1 {
				*	qui summ CI_lower if FRD==1 & naive==0 & t==`refCt' & fixedtau==0 & isw==0 & isrightside==1 & b==1 & p==`p'
				*	ereturn scalar refC_CI_lower =  = r(min)	
				*	qui summ CI_upper if FRD==1 & naive==0 & t==`refCt' & fixedtau==0 & isw==0 & isrightside==1 & b==1 & p==`p'
				*	ereturn scalar refC_CI_upper =  = r(max)			
				*}
			
		ereturn matrix treatment_effects_`stub' = treatment_effects_`stub'
		
		*Get point estimates and confidence intervals for the potential impact of manipulation (fixed tau)
		*-----------------------------------------------------------------------

			if `num_taus' > 0 {
				matrix fixedtau_`stub' = J(`num_taus', 9, .)
				matrix colnames fixedtau_`stub'= "tau" "SRD lower" "SRD upper" "SRD CI lower" "SRD CI upper" "FRD lower" "FRD upper" "FRD CI lower" "FRD CI upper"
				
				forvalues i = 1/`num_taus' {
					
					matrix fixedtau_`stub'[`i',1] = potential_taus[`i',1]
				
					*Note: using b=-1 for point estimates since we've dropped b=0 to save space (there would be no difference since fixed taus are not nudged)
				
					*SRD
					qui summ gammaL if FRD==0 & naive==0 & fixedtau==`i' & isw==0 & isrightside==0 & b==0 & p==`p'
					matrix fixedtau_`stub'[`i',2] = r(mean)
					qui summ gammaU if FRD==0 & naive==0 & fixedtau==`i' & isw==0 & isrightside==0 & b==0 & p==`p'
					matrix fixedtau_`stub'[`i',3] = r(mean)
					
					qui summ CI_lower if FRD==0 & naive==0 & fixedtau==`i' & isw==0 & isrightside==0 & b==1 & p==`p'
					matrix fixedtau_`stub'[`i',4] = r(mean)
					qui summ CI_upper if FRD==0 & naive==0 & fixedtau==`i' & isw==0 & isrightside==0 & b==1 & p==`p'
					matrix fixedtau_`stub'[`i',5] = r(mean)
					
					*FRD
					qui summ gammaL if FRD==1 & naive==0 & fixedtau==`i' & isw==0 & isrightside==0 & b==0 & p==`p'
					matrix fixedtau_`stub'[`i',6] = r(min)
					qui summ gammaU if FRD==1 & naive==0 & fixedtau==`i' & isw==0 & isrightside==0 & b==0 & p==`p'
					matrix fixedtau_`stub'[`i',7] = r(max)
					
					qui summ CI_lower if FRD==1 & naive==0 & fixedtau==`i' & isw==0 & isrightside==0 & b==1 & p==`p'
					matrix fixedtau_`stub'[`i',8] = r(min)
					qui summ CI_upper if FRD==1 & naive==0 & fixedtau==`i' & isw==0 & isrightside==0 & b==1 & p==`p'
					matrix fixedtau_`stub'[`i',9] = r(max)
				}
				
				ereturn matrix fixedtau_`stub' = fixedtau_`stub'
			}
	}
	
mata mata clear

timer off 51
qui timer list 51
if r(t51) < 60 {
	disp "Time taken: `r(t51)' seconds"
}
else if r(t51) < 3600 {
	local minutes_taken = r(t51)/60
	disp "Time taken: `minutes_taken' minutes"
}
else {
	local hours_taken = r(t51)/3600
	disp "Time taken: `hours_taken' hours"
}
ereturn scalar seconds_taken = r(t51)

end


capture program drop rdbounds_sampledata
program define rdbounds_sampledata, eclass

	*Check syntax	
	syntax, [samplesize(real 50000) covs clear]

	
	if ("`clear'"=="clear") {
		clear
	}
	else{
		di as error  "{err} The function rdbounds_sampledata will clear the data in memory. Please set option{cmd: clear} to confirm."  
		exit
	}
	
	if ("`covs'"=="covs") {
		local covs=1
	}
	else{
		local covs=0
	}
	
	qui set obs `samplesize'
	
	if `covs'==0 {
		*x values for potentially-assigned, ~ N(0,5) censored at -10 and 10
		local n0 = floor(`samplesize'*.95)
		qui gen x = rnormal(0,1)*sqrt(5) if _n <= `n0'
		qui replace x = -10 if x<-10 & x !=.
		qui replace x = 10 if x>10 & x !=.
		
		*x values for always-assigned ~ triangular decline from 0 to 5. f(x) = 2/5*(1-x/5); F(x) = 2/5x-x^2/25; F^{-1}(u) = 5*(1-sqrt(1-u))
		qui replace x = 5*(1-sqrt(runiform())) if _n > `n0'
		
		qui gen alwaysassigned = (_n>`n0')
		
		*introduce taker groups
		qui gen temprand = runiform()
			qui gen alwaystaker = (temprand<0.05)
			qui gen nevertaker = (temprand>=0.05) & (temprand<0.3)
			qui gen complier = 1-alwaystaker-nevertaker
			qui gen treatment = (alwaystaker+(x>=0)*complier)*(1-nevertaker)
		drop temprand  
		
		*outcome variable. Treatment effect is 2 for potentially assigned and 5 for always assigned units
		qui gen y = (x+10)/2+2*treatment*(alwaysassigned==0)+5*treatment*(alwaysassigned==1)+rnormal(0,1)
		qui replace y = 23 if y>23
		qui replace y = 0 if y<0
	}	
	else {
		local ss0 = floor(`samplesize'/2)
		local ss1 = `samplesize'-`ss0'
		gen cov = _n > `ss0'
		
		qui gen x=.
		qui gen alwaysassigned =.
		
		*introduce taker groups
			qui gen temprand = runiform()
				qui gen alwaystaker = (temprand<0.05)
				qui gen nevertaker = (temprand>=0.05) & (temprand<0.3)
				qui gen complier = 1-alwaystaker-nevertaker
			drop temprand  
			
		*cov==0 ("Group 1")
			local n0 = floor(`ss0'*.975)
			qui replace alwaysassigned = (_n>`n0') if cov==0
			qui replace x = rnormal(0,1)*sqrt(5) if alwaysassigned==0 & cov==0
			qui replace x = 5*(1-sqrt(runiform())) if alwaysassigned==1 & cov==0
			
			local n1 = `ss0'-`n0'
			scalar tau = (.025/.975*2/5)/(normalden(0, 0, sqrt(5)))
			local tau = round(scalar(tau),.00001)
			local trueRHS = round((1-scalar(tau))*2+scalar(tau)*5,.00001)
			disp as result "---------------------------------------------"
			disp as result "Subsample: cov = 0 (group 1)"
			disp as result "---------------------------------------------"
			disp as result "True tau: `tau'"
			disp as result "True treatment effect on potentially-assigned: 2"
			disp as result "True treatment effect on right side of cutoff: `trueRHS'"

		*cov==1 ("Group 2")
			local n0 = floor(`ss1'*.925)
			qui replace alwaysassigned = (_n>`n0'+`ss0') if cov==1
			qui replace x = rnormal(0,1)*sqrt(5) if alwaysassigned==0 & cov==1
			qui replace x = 5*(1-sqrt(runiform())) if alwaysassigned==1 & cov==1
			
			local n1 = `ss1'-`n0'
			scalar tau = (.075/.925*2/5)/(normalden(0, 0, sqrt(5)))
			local tau = round(scalar(tau),.00001)
			local trueRHS = round((1-scalar(tau))*2+scalar(tau)*5,.00001)
			disp as result "---------------------------------------------"
			disp as result "Subsample: cov = 1 (group 2)"
			disp as result "---------------------------------------------"
			disp as result "True tau: `tau'"
			disp as result "True treatment effect on potentially-assigned: 2"
			disp as result "True treatment effect on right side of cutoff: `trueRHS'"
			
			disp as result "---------------------------------------------"
			disp as result "Full sample"
			disp as result "---------------------------------------------"
			
		*This part can be done across both groups simultaneously
		qui gen treatment = (alwaystaker+(x>=0)*complier)*(1-nevertaker)
		qui replace x = -10 if x<-10 & x !=.
		qui replace x = 10 if x>10 & x !=.
		qui gen y = (x+10)/2+2*treatment*(alwaysassigned==0)+5*treatment*(alwaysassigned==1)+rnormal(0,1)
		qui replace y = 23 if y>23
		qui replace y = 0 if y<0
		
	}
	
	local n0 = floor(`samplesize'*.95)
	local n1 = `samplesize'-`n0'
	scalar tau = (.05/.95*2/5)/(normalden(0, 0, sqrt(5)))
	local tau = round(scalar(tau),.00001)
	local trueRHS = round((1-scalar(tau))*2+scalar(tau)*5,.00001)
	disp as result "True tau: `tau'"
	disp as result "True treatment effect on potentially-assigned: 2"
	disp as result "True treatment effect on right side of cutoff: `trueRHS'"

end
