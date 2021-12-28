*! version 1.0 20211203 David Veenman

/* 
20211223: 1.0.1		Added option "fast" for much faster execution based on mata
					Added small-sample correction to variance matrix and changed critical values based on 
					 standard normal to t(G-1) distribution following the advise of Cameron and Miller (2015)
20211203: 1.0.0		First version
*/

program define bootstep, eclass sortpreserve
	syntax anything [in] [if], nboot(integer) [cluster(varlist)] [absorb(varlist)] [logit] [standardize] [fast]

	if `"`absorb'"'!="" & `"`fast'"'!="" {
		di as text "ERROR: options <absorb> and <fast> cannot be combined"			
		exit
	}		
	if `"`logit'"'!="" & `"`fast'"'!="" {
		di as text "ERROR: options <logit> and <fast> cannot be combined: fast method only works on OLS"			
		exit
	}		
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
	
	* Expand factor variables in second step:
	fvexpand `indepv2'
	local indepv2 `r(varlist)'
	local fvcheck `r(fvops)'
	capture _fv_check_depvar `indepv1'
	if _rc>0 & `"`fast'"'!=""{
		di as text "ERROR: factor variables not allowed in first step with fast estimation; instead create separate dummies and include them."			
		exit		
	}

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
	qui gen `n'=_n if `touse'
	qui sum `n'
	local N=r(max)

	if (`nc'==2){	    
		tempvar intersection
		qui egen `intersection'=group(`clusterdim1' `clusterdim2') if `touse'
		qui sum `intersection'
		local ncluster3=r(max)
		if (`ncluster3'==`N'){
			local white="true"
		}
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

	if `"`fast'"'=="" {
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
		
			if "`fvcheck'"=="true"{
				local K=rowsof(V)-1 
			}
			else{
				local K=rowsof(V) 
			}
			
			if (`nc'==0){
				scalar factor1=`N'/(`N'-`K')
			}
			else{
				scalar factor1=(`ncluster1'/(`ncluster1'-1))*((`N'-1)/(`N'-`K'))
			}			
			matrix V=factor1*V
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

			// Obtain combined covariance matrix:
			matrix V=V1+V2-V3

			// Small-sample correction:
			if "`fvcheck'"=="true"{
				local K=rowsof(V)-1 
			}
			else{
				local K=rowsof(V) 
			}
			scalar factor1=(`ncluster1'/(`ncluster1'-1))*((`N'-1)/(`N'-`K'))
			scalar factor2=(`ncluster2'/(`ncluster2'-1))*((`N'-1)/(`N'-`K'))
				
			if `ncluster1'<`ncluster2'{
				scalar factormin=factor1
			}
			else{
				scalar factormin=factor2
			}
			matrix V=factormin*V
			
		}
		
		if (`nc'==0){
			scalar e_df_r=`N'-`K'
		}
		if (`nc'==1){
			scalar e_df_r=`ncluster1'-1
		}
		if (`nc'==2){
			if (`ncluster1'<`ncluster2'){
				scalar e_df_r=`ncluster1'-1
			}
			else{
				scalar e_df_r=`ncluster2'-1
			}
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
		ereturn scalar df_r=e_df_r
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
			ereturn scalar df_r=e_df_r
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
	}
	// Fast execution using mata optimized functions
	else{
		if `"`cluster'"'!=""{
			// One cluster dimension:
			if (`nc'==1){
				tempvar clusterid cvarn 
				egen `clusterid'=group(`clusterdim1') if `touse'

				// Mata's panelsetup requires sorting on the cluster variable
				sort `clusterid' 
				qui by `clusterid': gen _cvarn_temp=_n if `touse'
				qui replace _cvarn_temp=. if _cvarn_temp>1		
		
				gen _cons_temp=1 if `touse'
			
				local y "`depv1'"
				local xlist "`indepv1' _cons_temp"
				local cvar "`clusterid'"
				local cvarn "`clusterid' _cvarn_temp"
				local y2 "`depv2'"
				local xlist2 "`indepv2'"
				local xlistname "`indepv2' `depv1'_hat _cons"

				mata: _bootfast(`nboot')
				
				// Small-sample correction:
				if "`fvcheck'"=="true"{
					local K=rowsof(V)-1 
				}
				else{
					local K=rowsof(V) 
				}
				scalar factor1=(`ncluster1'/(`ncluster1'-1))*((`N'-1)/(`N'-`K'))
				matrix V=factor1*V
				
				drop _cons_temp _cvarn_temp
				matrix colnames b =`xlistname'
				matrix colnames V =`xlistname'
				matrix rownames V =`xlistname'	
			}
			// Two cluster dimensions:
			if (`nc'==2){
				tempvar clusterid1 cvarn1 clusterid2 cvarn2
				egen `clusterid1'=group(`clusterdim1') if `touse'
				egen `clusterid2'=group(`clusterdim2') if `touse'
		
				gen _cons_temp=1 if `touse'

				// Mata's panelsetup requires sorting on the cluster variable
				sort `clusterid1' 
				qui by `clusterid1': gen _cvarn_temp1=_n if `touse'
				qui replace _cvarn_temp1=. if _cvarn_temp1>1		
			
				local y "`depv1'"
				local xlist "`indepv1' _cons_temp"
				local cvar "`clusterid1'"
				local cvarn "`clusterid1' _cvarn_temp1"
				local y2 "`depv2'"
				local xlist2 "`indepv2'"
				local xlistname "`indepv2' `depv1'_hat _cons"

				// First dimension:
				mata: _bootfast(`nboot')
				matrix V1=V

				// Second dimension:
				sort `clusterid2' 
				qui by `clusterid2': gen _cvarn_temp2=_n if `touse'
				qui replace _cvarn_temp2=. if _cvarn_temp2>1		
				local cvar "`clusterid2'"
				local cvarn "`clusterid2' _cvarn_temp2"
				mata: _bootfast(`nboot')
				matrix V2=V
				
				// Third dimension:
				if "`white'"!="true"{
					sort `intersection' 
					qui by `intersection': gen _cvarn_temp3=_n if `touse'
					qui replace _cvarn_temp3=. if _cvarn_temp3>1		
					local cvar "`intersection'"
					local cvarn "`intersection' _cvarn_temp3"
					mata: _bootfast(`nboot')
					drop _cvarn_temp3
				}
				else{
					mata: _bootfast_nocl(`nboot')
				}
				matrix V3=V
				
				// Obtain combined covariance matrix:
				matrix V=V1+V2-V3

				// Small-sample correction:
				if "`fvcheck'"=="true"{
					local K=rowsof(V)-1 
				}
				else{
					local K=rowsof(V) 
				}
				scalar factor1=(`ncluster1'/(`ncluster1'-1))*((`N'-1)/(`N'-`K'))
				scalar factor2=(`ncluster2'/(`ncluster2'-1))*((`N'-1)/(`N'-`K'))
					
				if `ncluster1'<`ncluster2'{
					scalar factormin=factor1
				}
				else{
					scalar factormin=factor2
				}
				matrix V=factormin*V
				
				drop _cons_temp _cvarn_temp1 _cvarn_temp2 
				matrix colnames b =`xlistname'
				matrix colnames V =`xlistname'
				matrix rownames V =`xlistname'	
			}
				
			/* For standardization of fit coefficient */
			if `"`standardize'"'!=""{
				qui reg `depv1' `indepv1' if `touse'
				qui predict `fit' if `touse' & `depv1'!=., xb
				qui `est' `depv1' `indepv2', `feoption'
				qui predict `res_depv1', res
				qui `est' `fit' `indepv2', `feoption'
				qui predict `res_fit', res
				qui sum `res_depv1'
				qui scalar v1=r(Var)
				qui sum `res_fit'
				qui scalar v2=r(Var)
				qui `est' `depv2' `indepv2' `fit'  if `touse' & `depv1'!=., `feoption'
				matrix coef=e(b)
				matrix colnames coef=`indepv2' "`depv1'_hat" "_cons"
			}
		
			ereturn clear
			tempname b V
			matrix `b' = b
			matrix `V' = V
			ereturn post `b' `V'
			ereturn scalar N=e_N
			ereturn scalar r2_a=e_r2_a
			ereturn scalar df_r=e_df_r
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
				ereturn scalar df_r=e_df_r
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
		}
		else{
			local y "`depv1'"
			local xlist "`indepv1'"
			local y2 "`depv2'"
			local xlist2 "`indepv2'"
			local xlistname "`indepv2' `depv1'_hat _cons"
			
			mata: _bootfast_nocl(`nboot')
			matrix colnames b =`xlistname'
			matrix colnames V =`xlistname'
			matrix rownames V =`xlistname'	
			if "`fvcheck'"=="true"{
				local K=rowsof(V)-1 
			}
			else{
				local K=rowsof(V) 
			}
			
			matrix V=(`N'/(`N'-`K'))*V
			
			ereturn clear
			tempname b V
			matrix `b' = b
			matrix `V' = V
			ereturn post `b' `V'
			ereturn scalar N=e_N
			ereturn scalar r2_a=e_r2_a
			ereturn scalar df_r=`N'-`K'
			ereturn local depvar "`depv2'"
			eret local vcetype "Bootstrap"

			di " "
			di " " 
			if (`nc'==0){
				di as text "Two-step estimator with bootstrapped standard errors" 
			}
			di " "
			di as text "First-step estimation (`estlabel'): `depv1' = f(`indepv1')"
			di as text "Second-step estimation (OLS): `depv2' = f(`indepv2' `depv1'_hat)"
			
			di " "
			di as text "Number of bootstrap samples = " _column(31) %5.0f in yellow `nboot' ///
				_column(50) as text "Number of obs =       " %7.0f in yellow e(N)
			di _column(50) as text "Adj. R2 =             " %7.4f in yellow e(r2_a)
			ereturn display					
		}
	} 
end

mata:
	mata clear
	void _bootfast(real scalar B){
		st_view(y=., ., st_local("y"), st_local("touse"))
		st_view(X=., ., tokens(st_local("xlist")), st_local("touse"))
		st_view(cvar=., ., tokens(st_local("cvar")), st_local("touse"))
		st_view(cvarn=., ., tokens(st_local("cvarn")), st_local("touse"))
		st_view(y2=., ., st_local("y2"), st_local("touse"))
		st_view(X2=., ., tokens(st_local("xlist2")), st_local("touse"))
		
		n=rows(X)
		constant=J(n,1,1)
		XXinv=invsym(cross(X,X))
		b=XXinv*cross(X,y)
		coef=b'
		k=cols(X)
		k2=cols(X2)+1
		
		fit=X*coef'
		X2fit=(X2,fit,constant)
		XXinv2=invsym(cross(X2fit,X2fit))
		b2=XXinv2*cross(X2fit,y2)
		coef2=b2'

		e=y2-X2fit*b2
		my2=colsum(y2)/n
		m=y2 :- my2
		rss=(e'e)
		tss=(m'm)
		r2=1-rss/tss
		r2a=1-(rss/(n-(k2+1)))/(tss/(n-1))
		
        info=panelsetup(cvar, 1)
        nc=rows(info)		
		mm_panels(cvar, Cinfo=.)
		cl=1::rows(Cinfo)

		XXg=J(0,k,0)
        Xyg=J(0,1,0)
		for(i=1; i<=nc; i++) {
			xg=panelsubmatrix(X,i,info)
			xgxg=cross(xg,xg)
			XXg=(XXg\xgxg)
			yg=panelsubmatrix(y,i,info)
			xyg=cross(xg,yg)
			Xyg=(Xyg\xyg)
		}
		beta2=J(0,k2+1,0)
		for(b=1; b<=B; b++) {
			// STEP 1:
			// Create vectors indicating bootstrapped obs (n) and clusters (c)
			bsample_n=mm_sample(nc,.,Cinfo)
			bsample_c=cvarn[bsample_n,.]
			bsample_c=select(bsample_c,rowmissing(bsample_c):==0)
			bsample_c=bsample_c[.,1]
			
			// Obtain first-step coefficients:
			XXb=J(k,k,0)
			Xyb=J(k,1,0)
			submatx=J(k,k,0)
			submaty=J(n,1,0)
			a=bsample_c
			for(i=1; i<=nc; i++) {
				cluster=a[i]
				pos0=1+k*(cluster-1)
				pos=(pos0::pos0+k-1)
				submatx=XXg[pos,.]
				XXb=XXb+submatx
				submaty=Xyg[pos,.]
				Xyb=Xyb+submaty			
			}			
			XXbinv=invsym(XXb)
			beta1=cross(XXbinv,Xyb)

			// STEP 2:
			// get X matrix for same bootstrap sample --> NxK
			// multiply beta1 with X for bootstrap sample --> Nx1
			// add Nx1 vector to X2 matrix
			// get b2
			// get variance of b2			
			n2=rows(bsample_n)
			constant2=J(n2,1,1)
			Xb=X[bsample_n,.]
			X2b=X2[bsample_n,.]
			y2b=y2[bsample_n,.]
			fit1=Xb*beta1
			X2b=(X2b,fit1,constant2)
			bstep2=invsym(cross(X2b,X2b))*cross(X2b,y2b)
			beta2=(beta2 \ bstep2')
		}
		V=variance(beta2)
		dfr=nc-1
		st_matrix("b", coef2)
		st_matrix("V", V) 
		st_numscalar("e_N", n)
		st_numscalar("e_r2", r2)
		st_numscalar("e_r2_a", r2a)
		st_numscalar("e_df_r", dfr)
	}
	
	void _bootfast_nocl(real scalar B){
		st_view(y=., ., st_local("y"), st_local("touse"))
		st_view(X=., ., tokens(st_local("xlist")), st_local("touse"))
		st_view(y2=., ., st_local("y2"), st_local("touse"))
		st_view(X2=., ., tokens(st_local("xlist2")), st_local("touse"))
		
		n=rows(X)
		constant=J(n,1,1)
		X=(X,constant)
		XXinv=invsym(cross(X,X))
		b=XXinv*cross(X,y)
		coef=b'
		k=cols(X)
		k2=cols(X2)+1
		
		fit=X*coef'
		X2fit=(X2,fit,constant)
		XXinv2=invsym(cross(X2fit,X2fit))
		b2=XXinv2*cross(X2fit,y2)
		coef2=b2'

		e=y2-X2fit*b2
		my2=colsum(y2)/n
		m=y2 :- my2
		rss=(e'e)
		tss=(m'm)
		r2=1-rss/tss
		r2a=1-(rss/(n-(k2+1)))/(tss/(n-1))
		
		Xb=J(n,k,0)
		yb=J(n,1,0)
		X2b=J(n,k2+1,0)
		y2b=J(n,1,0)
		beta2=J(0,k2+1,0)
		for(b=1; b<=B; b++) {
			// Random sampling with replacement:
			p=mm_sample(rows(X),rows(X))
			Xb=X[p,.]
			yb=y[p,.]
			X2b=X2[p,.]
			y2b=y2[p,.]
			// First-step estimation:
			bstep1=invsym(cross(Xb,Xb))*cross(Xb,yb)
			fit1=Xb*bstep1
			// Second-step estimation:
			X2b=(X2b,fit1,constant)
			bstep2=invsym(cross(X2b,X2b))*cross(X2b,y2b)
			beta2=(beta2 \ bstep2')			
		}
		V=variance(beta2)
		st_matrix("b", coef2)
		st_matrix("V", V) 
		st_numscalar("e_N", n)
		st_numscalar("e_r2", r2)
		st_numscalar("e_r2_a", r2a)
	}
	
end
