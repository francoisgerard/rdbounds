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
