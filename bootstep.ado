*! version 1.0 20211203 David Veenman

/* 
20211203: 1.0		First version
*/

program define bootstep, eclass sortpreserve
	syntax anything [in] [if], [cluster(varlist)] nboot(integer) [absorb(varlist)] [logit] [standardize]

	if `"`logit'"'!="" & `"`standardize'"'!="" {
		di as text "ERROR: options <logit> and <standardize> cannot be combined: standardization only works on OLS"			
		exit
	}		
	if `"`absorb'"'!=""{
		local feoption="absorb(`absorb')" 
		local est="areg"
	}
	else{
		local feoption=""
		local est="reg"
	}
	if `"`logit'"'!="" {
		local estlabel="logit"
	}	
	else{
		local estlabel="OLS"
	}
		
	marksample touse
	tokenize `anything', parse(|)
	local step1 "`1'"
	macro shift 1
	local step2 "`*'"
	tokenize `step2'
	macro shift 1
	local step2 "`*'"
	tokenize `step1'
	local depv1 `"`1'"'
	macro shift 1
	local indepv1 "`*'"
	tokenize `step2'
	local depv2 `"`1'"'
	macro shift 1
	local indepv2 "`*'"
	
	* Allow factor variables in second step:
	fvexpand `indepv2'
	local indepv2 `r(varlist)'

	* Ensure dv's are not factor variables:
	_fv_check_depvar `depv1'
	_fv_check_depvar `depv2'
	
	if `"`logit'"'!="" {
		qui sum `depv1'
		if r(min)==0 & r(max)==1{
			qui sum `depv1' if `depv1'>0 & `depv1'<1
			if r(N)!=0{
				di as text "ERROR: option <logit> should be combined with first-stage dummy variable as dependent"
				exit
			}
		}
	}
	
	local nc: word count `cluster'
	if (`nc'>2){
		di as text "ERROR: Maximum number of dimensions to cluster on is two"
		exit
	}
	if (`nc'>0){
		local clusterdim1: word 1 of `cluster'
	}
	if (`nc'>1){
		local clusterdim2: word 2 of `cluster'
	}
	
	tempvar n fit res f r res_depv1 res_fit cl1 cl2
	qui gen `n'=_n

	if (`nc'==2){	    
		tempvar intersection
		qui egen `intersection'=group(`clusterdim1' `clusterdim2') if `touse'
	}	
	if (`nc'>0){
		qui egen `cl1'=group(`clusterdim1') if `touse'
		qui sum `cl1'
		local ncluster1=r(max)
	}
	if (`nc'>1){
		qui egen `cl2'=group(`clusterdim2') if `touse'
		qui sum `cl2'
		local ncluster2=r(max)
	}
	
	di  ""
	di as text "Performing bootstrap resampling... (B=`nboot')"
	
	if(`nc'<=1){
		/* Bootstrap resampling in first cluster dimension*/
		forvalues i=1(1)`nboot'{
			preserve
			bsample, cluster(`clusterdim1')
			if `"`logit'"'!=""{
				qui logit `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui predict `res' if `touse', res			
				qui replace `fit'=exp(`fit')/(1+exp(`fit'))
				qui replace `res'=exp(`res')/(1+exp(`res'))
			}	
			else{
				qui reg `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui predict `res' if `touse', res
			}
			qui `est' `depv2' `indepv2' `fit'  if `touse' & `depv1'!=., `feoption'
			matrix b_0=e(b)'
			if (`i'==1){
				matrix b=b_0
			}
			else{
				matrix b=(b, b_0)
			}
			restore
			di ".", _continue
		}
		mata: b=st_matrix("b")
		mata: var=diagonal(variance(b'))
		mata: st_matrix("V", var)
		matrix rownames V=`indepv2' "`depv1'_hat" "_cons"
		matrix V=diag(V)
	}
	
	if(`nc'==2){
		forvalues i=1(1)`nboot'{
			/* Bootstrap resampling in first cluster dimension*/
			preserve
			bsample, cluster(`clusterdim1')
			if `"`logit'"'!=""{
				qui logit `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui predict `res' if `touse', res			
				qui replace `fit'=exp(`fit')/(1+exp(`fit'))
				qui replace `res'=exp(`res')/(1+exp(`res'))
			}	
			else{
				qui reg `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui predict `res' if `touse', res
			}
			qui `est' `depv2' `indepv2' `fit'  if `touse' & `depv1'!=., `feoption'
			matrix b1_0=e(b)'
			if (`i'==1){
				matrix b1=b1_0
			}
			else{
				matrix b1=(b1, b1_0)
			}
			restore

			mata: b1=st_matrix("b1")
			mata: var=diagonal(variance(b1'))
			mata: st_matrix("V1", var)
			matrix rownames V1=`indepv2' "`depv1'_hat" "_cons"
			matrix V1=diag(V1)
		
			/* Bootstrap resampling in second cluster dimension*/
			preserve
			bsample, cluster(`clusterdim2')
			if `"`logit'"'!=""{
				qui logit `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui predict `res' if `touse', res			
				qui replace `fit'=exp(`fit')/(1+exp(`fit'))
				qui replace `res'=exp(`res')/(1+exp(`res'))
			}	
			else{
				qui reg `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui predict `res' if `touse', res
			}
			qui `est' `depv2' `indepv2' `fit'  if `touse' & `depv1'!=., `feoption'
			matrix b2_0=e(b)'
			if (`i'==1){
				matrix b2=b2_0
			}
			else{
				matrix b2=(b2, b2_0)
			}
			restore

			mata: b2=st_matrix("b2")
			mata: var=diagonal(variance(b2'))
			mata: st_matrix("V2", var)
			matrix rownames V2=`indepv2' "`depv1'_hat" "_cons"
			matrix V2=diag(V2)
		
			/* Bootstrap resampling in intersection of cluster dimensions*/
			preserve
			bsample, cluster(`intersection')
			if `"`logit'"'!=""{
				qui logit `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui predict `res' if `touse', res			
				qui replace `fit'=exp(`fit')/(1+exp(`fit'))
				qui replace `res'=exp(`res')/(1+exp(`res'))
			}	
			else{
				qui reg `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui predict `res' if `touse', res
			}
			qui `est' `depv2' `indepv2' `fit'  if `touse' & `depv1'!=., `feoption'
			matrix b3_0=e(b)'
			if (`i'==1){
				matrix b3=b3_0
			}
			else{
				matrix b3=(b3, b3_0)
			}
			restore
			mata: b3=st_matrix("b3")
			mata: var=diagonal(variance(b3'))
			mata: st_matrix("V3", var)
			matrix rownames V3=`indepv2' "`depv1'_hat" "_cons"
			matrix V3=diag(V3)
			di ".", _continue
		}
		
		/* Construct twoway cluster-robust variance matrix */
		matrix V=V1+V2-V3
		
	}
	
	if `"`logit'"'!=""{
		qui logit `depv1' `indepv1' if `touse'
		qui predict `fit' if `touse' & `depv1'!=., xb
		qui predict `res' if `touse', res			
		qui replace `fit'=exp(`fit')/(1+exp(`fit'))
		qui replace `res'=exp(`res')/(1+exp(`res'))
	}	
	else{
		qui reg `depv1' `indepv1' if `touse'
		qui predict `fit' if `touse' & `depv1'!=., xb
		qui predict `res' if `touse', res
	}
				
	qui `est' `depv2' `indepv2' `fit'  if `touse' & `depv1'!=., `feoption'
	scalar e_N=e(N)
	scalar e_r2_a=e(r2_a)
	matrix coef=e(b)
	matrix colnames coef=`indepv2' "`depv1'_hat" "_cons"
		
	/* For standardization of fit coefficient */
	if `"`standardize'"'!=""{
		qui `est' `depv1' `indepv2', `feoption'
		qui predict `res_depv1', res
		qui `est' `fit' `indepv2', `feoption'
		qui predict `res_fit', res
		qui sum `res_depv1'
		qui scalar v1=r(Var)
		qui sum `res_fit'
		qui scalar v2=r(Var)
	}

	/* Post results in e() */
	ereturn clear
	tempname b V
	matrix `b' = coef
	matrix `V' = V
	ereturn post `b' `V'
	ereturn scalar N=e_N
	ereturn scalar r2_a=e_r2_a
	ereturn local depvar "`depv2'"
	eret local vcetype "Bootstrap"
	
	if `"`standardize'"'!=""{
		matrix b=e(b)'
		local dim=rowsof(b)
		matrix coef[1,`dim'-1]=coef[1,`dim'-1]*(v2/v1)
		matrix V[`dim'-1,`dim'-1]=V[`dim'-1,`dim'-1]*(v2/v1)^2
		ereturn clear
		tempname b V
		matrix `b' = coef
		matrix `V' = V
		ereturn post `b' `V'
		ereturn scalar N=e_N
		ereturn scalar r2_a=e_r2_a
		ereturn local depvar "`depv2'"
		eret local vcetype "Bootstrap"
	}
	
	di " "
	di " " 
	if (`nc'==0){
		di as text "Two-step estimator with bootstrapped standard errors" 
	}
	if (`nc'==1){
		di as text "Two-step estimator with bootstrapped standard errors" 
		di as text "Pairs-cluster bootstrap by cluster variable `clusterdim1'"
	}
	if (`nc'==2){
		di as text "Two-step estimator with bootstrapped standard errors" 
		di as text "Pairs-cluster bootstrap by cluster variables `clusterdim1' and `clusterdim2'"
	}
	di " "
	di as text "First-step estimation (`estlabel'): `depv1' = f(`indepv1')"
	di as text "Second-step estimation (OLS): `depv2' = f(`indepv2' `depv1'_hat)"
	
	di " "
	di as text "Number of bootstrap samples = " _column(31) %5.0f in yellow `nboot' ///
		_column(50) as text "Number of obs =       " %7.0f in yellow e(N)
	if (`nc'==0){
	    di _column(50) as text "Adj. R2 =             " %7.4f in yellow e(r2_a)
	}
	if (`nc'>0){
		di as text "Number of clusters (`clusterdim1') = " _column(31) %5.0f in yellow `ncluster1' ///
			_column(50) as text "Adj. R2 =             " %7.4f in yellow e(r2_a)
	}
	if (`nc'>1){
		di as text "Number of clusters (`clusterdim2') = " _column(31) %5.0f in yellow `ncluster2' 
	}

	ereturn display
	
	di ""
	if `"`standardize'"'!=""{
		loc url "https://doi.org/10.1016/j.jfineco.2016.02.013"
		local message "Note: coefficient on generated regressor is standardized following"
		di as text `"`message' {browse "`url'":Hou and Loh (2016, JFE).}"'
	}
	
	if (`nc'<=1){
		matrix drop coef b b_0 V
	}
	if (`nc'==2){
		matrix drop coef b1 b1_0 V1 b2 b2_0 V2 b3 b3_0 V3 V
	}
end
