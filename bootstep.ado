*! version 1.2.0 20241118 David Veenman

/*
20241118: 1.2.0     Added Cameron/Gelbach/Miller (2011) adjustment for non-positive-semidefinite VCE (see also Gu and Yoo 2019, DOI: 10.1177/1536867X19893637) 
20230106: 1.1.1     Minor update: changed reference for standardization option
20221216: 1.1.0     Major improvements:
                     - Changed default to fast Mata estimation
                     - Absorb option in step 2 now also works with Mata estimation
                     - Factor variables now allowed in step 1
                     - Added restriction to touse variable that observations much be available for all variables of both steps
					 - Added Mata timer to display status within bootstrap procedure
                     - Added seed option
20211229: 1.0.2     Fixed issue with fast option for twoway clustering 
20211223: 1.0.1     Added option "fast" for much faster execution based on mata
                    Added small-sample correction to variance matrix and changed critical values based on 
                     standard normal to t(G-1) distribution following the advise of Cameron and Miller (2015)
20211203: 1.0.0     First version
*/

program define bootstep, eclass sortpreserve
    syntax anything [in] [if], nboot(integer) [cluster(varlist)] [absorb(varlist)] [standardize] [logit] [seed(integer 0)]

    capture findfile lmoremata.mlib
    if _rc>0 {
        di as text "ERROR: Package moremata required (install via ssc inst moremata, replace)."            
        exit        
    }    
    if (`seed'==0){
        di ""
        di as text "Careful: make sure to use seed() option or the 'set seed' command in Stata to obtain reproducable results!"
    }
    else{
        set seed `seed'
    }    
    if `"`logit'"'!="" & `"`standardize'"'!="" {
        di as text "ERROR: options <logit> and <standardize> cannot be combined: standardization only works for OLS"            
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
    
    // Expand factor variables in first step:
    fvexpand `indepv1'
    local indepv1 `r(varlist)'

    // Expand factor variables in second step:
    fvexpand `indepv2'
    local indepv2 `r(varlist)'
    local fvcheck `r(fvops)'
    if "`fvcheck'"=="true"{
        scalar fvdum=1
    }
    else {
        scalar fvdum=0
    }
    
    // Ensure dv's are not factor variables:
    _fv_check_depvar `depv1'
    _fv_check_depvar `depv2'

    // Ensure observations with a missing value in any of the variables are not used:
    tempvar rowsum
    qui gen `rowsum'=`depv1'
    qui replace `rowsum'=`rowsum'+`depv2'
    foreach var of local indepv1 {
        qui replace `rowsum'=`rowsum'+`var'
    }
    foreach var of local indepv2 {
        qui replace `rowsum'=`rowsum'+`var'
    }
    qui replace `touse'=0 if `rowsum'==.
        
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

    // Logit in first step:
    if `"`logit'"'!="" {
        di  ""
        di as text "Performing bootstrap resampling... (B=`nboot')"
        
        if(`nc'<=1){
            /* Bootstrap resampling in first cluster dimension*/
            forvalues i=1(1)`nboot'{
                preserve
                bsample, cluster(`clusterdim1')
                qui logit `depv1' `indepv1' if `touse'
                qui predict `fit' if `touse' & `depv1'!=., xb
                qui predict `res' if `touse', res            
                qui replace `fit'=exp(`fit')/(1+exp(`fit'))
                qui replace `res'=exp(`res')/(1+exp(`res'))
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
                qui logit `depv1' `indepv1' if `touse'
                qui predict `fit' if `touse' & `depv1'!=., xb
                qui predict `res' if `touse', res            
                qui replace `fit'=exp(`fit')/(1+exp(`fit'))
                qui replace `res'=exp(`res')/(1+exp(`res'))
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
                qui logit `depv1' `indepv1' if `touse'
                qui predict `fit' if `touse' & `depv1'!=., xb
                qui predict `res' if `touse', res            
                qui replace `fit'=exp(`fit')/(1+exp(`fit'))
                qui replace `res'=exp(`res')/(1+exp(`res'))
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
                qui logit `depv1' `indepv1' if `touse'
                qui predict `fit' if `touse' & `depv1'!=., xb
                qui predict `res' if `touse', res            
                qui replace `fit'=exp(`fit')/(1+exp(`fit'))
                qui replace `res'=exp(`res')/(1+exp(`res'))
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
        
        qui logit `depv1' `indepv1' if `touse'
        qui predict `fit' if `touse' & `depv1'!=., xb
        qui predict `res' if `touse', res            
        qui replace `fit'=exp(`fit')/(1+exp(`fit'))
        qui replace `res'=exp(`res')/(1+exp(`res'))
                    
        qui `est' `depv2' `indepv2' `fit'  if `touse' & `depv1'!=., `feoption'
        scalar e_N=e(N)
        scalar e_r2_a=e(r2_a)
        matrix coef=e(b)
        matrix colnames coef=`indepv2' "`depv1'_hat" "_cons"
            
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
               _column(50) as text "Number of obs =         " %7.0f in yellow e(N)
        if (`nc'==0){
            di _column(50) as text "Adj. R2 (step 2) =      " %7.4f in yellow e(r2_a)
        }
        if (`nc'>0){
            di as text "Number of clusters (`clusterdim1') = " _column(31) %5.0f in yellow `ncluster1' ///
               _column(50) as text "Adj. R2 (step 2) =      " %7.4f in yellow e(r2_a)
        }
        if (`nc'>1){
            di as text "Number of clusters (`clusterdim2') = " _column(31) %5.0f in yellow `ncluster2' 
        }

        ereturn display
                
        if (`nc'<=1){
            matrix drop coef b b_0 V
        }
        if (`nc'==2){
            matrix drop coef b1 b1_0 V1 b2 b2_0 V2 b3 b3_0 V3 V
        }
    }
    // OLS: fast execution using mata optimized functions
    else{
        if `"`cluster'"'!=""{
            // One cluster dimension:
            if (`nc'==1){
                tempvar clusterid 
                qui egen `clusterid'=group(`clusterdim1') if `touse'

                // Mata's panelsetup requires sorting on the cluster variable
                sort `clusterid' 
                
                local y "`depv1'"
                local xlist "`indepv1'"
                local cvar "`clusterid'"
                local y2 "`depv2'"
                local xlist2 "`indepv2'"
                local xlistname "`indepv2' `depv1'_hat _cons"

                di  ""
                di as text "Performing bootstrap resampling... (B=`nboot')"
        
                if "`absorb'"==""{
                    mata: _bootfast(`nboot')
                }
                else {
                    local ivar "`absorb'"
                    mata: _bootfast_fe(`nboot')
                }                
                
                // Small-sample correction:
                if "`fvcheck'"=="true"{
                    local K=rowsof(V)-1 
                }
                else{
                    local K=rowsof(V) 
                }
                scalar factor1=(`ncluster1'/(`ncluster1'-1))*((`N'-1)/(`N'-`K'))
                matrix V=factor1*V
                
                matrix colnames b =`xlistname'
                matrix colnames V =`xlistname'
                matrix rownames V =`xlistname'    
            }
            // Two cluster dimensions:
            if (`nc'==2){
                tempvar clusterid1 clusterid2 
                qui egen `clusterid1'=group(`clusterdim1') if `touse'
                qui egen `clusterid2'=group(`clusterdim2') if `touse'
        
                // Mata's panelsetup requires sorting on the cluster variable
                sort `clusterid1' 
            
                local y "`depv1'"
                local xlist "`indepv1'"
                local cvar "`clusterid1'"
                local y2 "`depv2'"
                local xlist2 "`indepv2'"
                local xlistname "`indepv2' `depv1'_hat _cons"

                // First dimension:
                di  ""
                di as text "Performing bootstrap resampling in first clustering dimension (`clusterdim1')... (B=`nboot')"
                if "`absorb'"==""{
                    mata: _bootfast(`nboot')
                }
                else {
                    local ivar "`absorb'"
                    mata: _bootfast_fe(`nboot')
                }                
                matrix V1=V

                // Second dimension:
                sort `clusterid2' 
                local cvar "`clusterid2'"
                di  ""
                di as text "Performing bootstrap resampling in second clustering dimension (`clusterdim2')... (B=`nboot')"
                if "`absorb'"==""{
                    mata: _bootfast(`nboot')
                }
                else {
                    local ivar "`absorb'"
                    mata: _bootfast_fe(`nboot')
                }                
                matrix V2=V
                
                // Third dimension:
                di  ""
                di as text "Performing bootstrap resampling in intersection of clustering dimensions... (B=`nboot')"
                if "`white'"!="true"{
                    sort `intersection' 
                    local cvar "`intersection'"
                    if "`absorb'"==""{
                        mata: _bootfast(`nboot')
                    }
                    else {
                        local ivar "`absorb'"
                        mata: _bootfast_fe(`nboot')
                    }                
                }
                else{
                    sort `n'
                    if "`absorb'"==""{
                        mata: _bootfast_nocl(`nboot')
                    }
                    else {
                        local ivar "`absorb'"
                        mata: _bootfast_nocl_fe(`nboot')
                    }                
                }
                matrix V3=V
                
                // Obtain combined covariance matrix:
                matrix V=V1+V2-V3

				mata: _check_vce()
				if (negative==1) {
					matrix colnames Vcadj=`indepvnames'
					matrix rownames Vcadj=`indepvnames'	
					matrix Vc=Vcadj
					di ""
					di "Note: adjustment from Cameron, Gelbach, and Miller (2011) applied to non-positive semi-definite VCE"
				}
				
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
            ereturn scalar r2=e_r2
            ereturn scalar r2_a=e_r2_a
            if ("`absorb'"!=""){
                ereturn scalar r2_within=e_r2_w
            }
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
                ereturn scalar r2=e_r2
                ereturn scalar r2_a=e_r2_a
                if ("`absorb'"!=""){
                    ereturn scalar r2_within=e_r2_w
                }
                ereturn scalar df_r=e_df_r
                ereturn local depvar "`depv2'"
                eret local vcetype "Bootstrap"
            }

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
                di _column(50) as text "Adj. R2 (step 2) =    " %7.4f in yellow e(r2_a)
            }
            if (`nc'>0){
                di as text "Number of clusters (`clusterdim1') = " _column(31) %5.0f in yellow `ncluster1' ///
                   _column(50) as text "R2 (step 2) =         " %7.4f in yellow e(r2)
            }
            if (`nc'==1){
                di _column(50) as text "Adj. R2 (step 2) =    " %7.4f in yellow e(r2_a)
            }
            if (`nc'>1){
                di as text "Number of clusters (`clusterdim2') = " _column(31) %5.0f in yellow `ncluster2' ///
                    _column(50) as text "Adj. R2 (step 2) =    " %7.4f in yellow e(r2_a)
            }
            if ("`absorb'"!=""){
                di as text "Absorbed FE in second step: `absorb'" ///
                    _column(50) as text "Within R2 (step 2) =  " %7.4f in yellow e(r2_within)
            }
            ereturn display        
        
            di ""
            if `"`standardize'"'!=""{
                loc url "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4310933"
                local message "Note: coefficient on generated regressor is standardized following"
                di as text `"`message' {browse "`url'":Cascino, Szeles, and Veenman (2024).}"'
            }
        }
        else{
            local y "`depv1'"
            local xlist "`indepv1'"
            local y2 "`depv2'"
            local xlist2 "`indepv2'"
            local xlistname "`indepv2' `depv1'_hat _cons"
            
            di  ""
            di as text "Performing bootstrap resampling... (B=`nboot')"
            if "`absorb'"==""{
                mata: _bootfast_nocl(`nboot')
            }
            else {
                local ivar "`absorb'"
                mata: _bootfast_nocl_fe(`nboot')
            }                
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
            ereturn scalar r2=e_r2
            ereturn scalar r2_a=e_r2_a
            if ("`absorb'"!=""){
                ereturn scalar r2_within=e_r2_w
            }
            ereturn scalar df_r=`N'-`K'
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
                ereturn scalar r2=e_r2
                ereturn scalar r2_a=e_r2_a
                if ("`absorb'"!=""){
                    ereturn scalar r2_within=e_r2_w
                }
                ereturn scalar df_r=e_df_r
                ereturn local depvar "`depv2'"
                eret local vcetype "Bootstrap"
            }
            
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
            di _column(50) as text "R2 (step 2) =         " %7.4f in yellow e(r2)
            di _column(50) as text "Adj. R2 (step 2) =    " %7.4f in yellow e(r2_a)
            if ("`absorb'"!=""){
                di as text "Absorbed FE in second step: `absorb'" ///
                    _column(50) as text "Within R2 (step 2) =  " %7.4f in yellow e(r2_within)
            }
            ereturn display                    
                    
            di ""
            if `"`standardize'"'!=""{
                loc url "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4310933"
                local message "Note: coefficient on generated regressor is standardized following"
                di as text `"`message' {browse "`url'":Cascino, Szeles, and Veenman (2024).}"'
            }
        }
    } 
end

mata:
    mata clear
        void _bootfast(real scalar B){
        st_view(y=., ., st_local("y"), st_local("touse"))
        st_view(X=., ., tokens(st_local("xlist")), st_local("touse"))
        st_view(cvar=., ., tokens(st_local("cvar")), st_local("touse"))
        st_view(y2=., ., st_local("y2"), st_local("touse"))
        st_view(X2=., ., tokens(st_local("xlist2")), st_local("touse"))
        fvc=st_numscalar("fvdum")        

        n=rows(X)
        X=(X,J(n,1,1))
        
        XXinv=invsym(cross(X,X))
        b=XXinv*cross(X,y)
        coef=b'
        k=cols(X)
        k2=cols(X2)+1
        
        fit=X*coef'
        X2fit=(X2,fit,J(n,1,1))
        XXinv2=invsym(cross(X2fit,X2fit))
        b2=XXinv2*cross(X2fit,y2)
        coef2=b2'

        e=y2-X2fit*b2
        my2=colsum(y2)/n
        m=y2 :- my2
        rss=(e'e)
        tss=(m'm)
        r2=1-rss/tss
        r2a=1-(rss/(n-(k2+1-fvc)))/(tss/(n-1))
        
        info=panelsetup(cvar, 1)
        nc=rows(info)        
        mm_panels(cvar, Cinfo=.)

        beta2=J(0,k2+1,0)
        bcounter=0
        bcounter2=0
        pct=0        
        for(b=1; b<=B; b++) {
            // Bootstrap sample:
            bsample_n=mm_sample(nc,.,Cinfo)
            Xb=X[bsample_n,.]
            yb=y[bsample_n,.]
            X2b=X2[bsample_n,.]
            y2b=y2[bsample_n,.]
            n2=rows(bsample_n)
            // First-step estimation:
            bstep1=invsym(cross(Xb,Xb))*cross(Xb,yb)
            fit1=Xb*bstep1
            // Second-step estimation:
            X2b=(X2b,fit1,J(n2,1,1))
            bstep2=invsym(cross(X2b,X2b))*cross(X2b,y2b)
            beta2=(beta2 \ bstep2')            
            
            _bootcounter(b, B, bcounter, bcounter2, pct)

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

    void _bootfast_fe(real scalar B){
        st_view(y=., ., st_local("y"), st_local("touse"))
        st_view(X=., ., tokens(st_local("xlist")), st_local("touse"))
        st_view(ivar=., ., tokens(st_local("ivar")), st_local("touse"))
        st_view(cvar=., ., tokens(st_local("cvar")), st_local("touse"))
        st_view(y2=., ., st_local("y2"), st_local("touse"))
        st_view(X2=., ., tokens(st_local("xlist2")), st_local("touse"))
        fvc=st_numscalar("fvdum")        
            
        n=rows(X)
        X=(X,J(n,1,1))
        
        XXinv=invsym(cross(X,X))
        b=XXinv*cross(X,y)
        coef=b'
        k=cols(X)
        k2=cols(X2)+1        
        
        fit=X*coef'
        X2fit=(X2,fit)
        S=mm_areg(y2, ivar, X2fit, 1, 1)
        coef2=mm_areg_b(S)'
        r2=mm_areg_r2(S)
        k_levels=mm_areg_k_levels(S)
        r2a=1-(((1-r2)*(n-1))/(n-k2-k_levels+fvc))
        
        // Get within-R2:
        yd=mm_areg_yd(S)
        Xd=mm_areg_Xd(S)
        Xd=(Xd, J(n,1,1))
        ed=yd-Xd*coef2'
        r2within=1-variance(ed)/variance(yd)
        
        info=panelsetup(cvar, 1)
        nc=rows(info)        
        mm_panels(cvar, Cinfo=.)

        beta2=J(0,k2+1,0)
        bcounter=0
        bcounter2=0
        pct=0
        for(b=1; b<=B; b++) {
            // Bootstrap sample:
            bsample_n=mm_sample(nc,.,Cinfo)
            Xb=X[bsample_n,.]
            yb=y[bsample_n,.]
            X2b=X2[bsample_n,.]
            y2b=y2[bsample_n,.]
            ivarb=ivar[bsample_n,.]
            // First-step estimation:            
            beta1=invsym(cross(Xb,Xb))*cross(Xb,yb)
            fit1=Xb*beta1
            // Second-step estimation:    
            X2b=(X2b,fit1)
            S = mm_areg(y2b, ivarb, X2b, 1, 1)
            bstep2 = mm_areg_b(S)
            beta2=(beta2 \ bstep2')

            _bootcounter(b, B, bcounter, bcounter2, pct)
            
        }
        V=variance(beta2)
        dfr=nc-1
        st_matrix("b", coef2)
        st_matrix("V", V) 
        st_numscalar("e_N", n)
        st_numscalar("e_r2", r2)
        st_numscalar("e_r2_a", r2a)
        st_numscalar("e_r2_w", r2within)
        st_numscalar("e_df_r", dfr)
    }
        
    void _bootfast_nocl(real scalar B){
        st_view(y=., ., st_local("y"), st_local("touse"))
        st_view(X=., ., tokens(st_local("xlist")), st_local("touse"))
        st_view(y2=., ., st_local("y2"), st_local("touse"))
        st_view(X2=., ., tokens(st_local("xlist2")), st_local("touse"))
        fvc=st_numscalar("fvdum")        
        
        n=rows(X)
        X=(X,J(n,1,1))
        
        XXinv=invsym(cross(X,X))
        b=XXinv*cross(X,y)
        coef=b'
        k=cols(X)
        k2=cols(X2)+1
        
        fit=X*coef'
        X2fit=(X2,fit,J(n,1,1))
        XXinv2=invsym(cross(X2fit,X2fit))
        b2=XXinv2*cross(X2fit,y2)
        coef2=b2'

        e=y2-X2fit*b2
        my2=colsum(y2)/n
        m=y2 :- my2
        rss=(e'e)
        tss=(m'm)
        r2=1-rss/tss
        r2a=1-(rss/(n-(k2+1-fvc)))/(tss/(n-1))
        
        Xb=J(n,k,0)
        yb=J(n,1,0)
        X2b=J(n,k2+1,0)
        y2b=J(n,1,0)
        beta2=J(0,k2+1,0)
        bcounter=0
        bcounter2=0
        pct=0
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
            X2b=(X2b,fit1,J(n,1,1))
            bstep2=invsym(cross(X2b,X2b))*cross(X2b,y2b)
            beta2=(beta2 \ bstep2')            
            
            _bootcounter(b, B, bcounter, bcounter2, pct)

        }
        V=variance(beta2)
        st_matrix("b", coef2)
        st_matrix("V", V) 
        st_numscalar("e_N", n)
        st_numscalar("e_r2", r2)
        st_numscalar("e_r2_a", r2a)
    }

    void _bootfast_nocl_fe(real scalar B){
        st_view(y=., ., st_local("y"), st_local("touse"))
        st_view(X=., ., tokens(st_local("xlist")), st_local("touse"))
        st_view(y2=., ., st_local("y2"), st_local("touse"))
        st_view(X2=., ., tokens(st_local("xlist2")), st_local("touse"))
        st_view(ivar=., ., tokens(st_local("ivar")), st_local("touse"))
        fvc=st_numscalar("fvdum")        
                
        n=rows(X)
        X=(X,J(n,1,1))
        
        XXinv=invsym(cross(X,X))
        b=XXinv*cross(X,y)
        coef=b'
        k=cols(X)
        k2=cols(X2)+1
        
        fit=X*coef'
        X2fit=(X2,fit)
        
        S=mm_areg(y2, ivar, X2fit, 1, 1)
        coef2=mm_areg_b(S)'
        r2=mm_areg_r2(S)
        k_levels=mm_areg_k_levels(S)
        r2a=1-(((1-r2)*(n-1))/(n-k2-k_levels+fvc))

        // Get within-R2:
        yd=mm_areg_yd(S)
        Xd=mm_areg_Xd(S)
        Xd=(Xd, J(n,1,1))
        ed=yd-Xd*coef2'
        r2within=1-variance(ed)/variance(yd)
        
        Xb=J(n,k,0)
        yb=J(n,1,0)
        X2b=J(n,k2+1,0)
        y2b=J(n,1,0)
        beta2=J(0,k2+1,0)
        bcounter=0
        bcounter2=0
        pct=0
        for(b=1; b<=B; b++) {
            // Random sampling with replacement:
            p=mm_sample(rows(X),rows(X))
            Xb=X[p,.]
            yb=y[p,.]
            X2b=X2[p,.]
            y2b=y2[p,.]
            ivarb=ivar[p,.]
            // First-step estimation:
            bstep1=invsym(cross(Xb,Xb))*cross(Xb,yb)
            fit1=Xb*bstep1
            // Second-step estimation:
            X2b=(X2b,fit1)
            S = mm_areg(y2b, ivarb, X2b, 1, 1)
            bstep2 = mm_areg_b(S)
            beta2=(beta2 \ bstep2')
 
            _bootcounter(b, B, bcounter, bcounter2, pct)

        }
        V=variance(beta2)
        st_matrix("b", coef2)
        st_matrix("V", V) 
        st_numscalar("e_N", n)
        st_numscalar("e_r2", r2)
        st_numscalar("e_r2_a", r2a)
        st_numscalar("e_r2_w", r2within)
    }
    
    void _bootcounter(
        real scalar b,
        real scalar B,
        real scalar bcounter,
        real scalar bcounter2,
        real scalar pct
        ){
        if (b==1) {
            printf("0%%")
            displayflush()                
        }
        bcounter++
        bcounter2++
        if (floor(5*bcounter/B)==1) {
            pct=pct+20
            printf("%f", pct)
            printf("%%")
            displayflush()
            bcounter=0
            bcounter2=0
        }
        else{
            if (floor(25*bcounter2/B)==1) {
                printf(".")
                displayflush()                                    
                bcounter2=0
            }
        }
        if (b==B) {
            printf("\n")
            displayflush()                
        }    
    }

	void _check_vce() {
		real scalar neg
		real matrix diag, EVEC, eval
		
		neg=0
		diag=diagonal(st_matrix("Vc"))
		for (i=1; i<=rows(diag); i++) { 
			if (diag[i]<=0) {
				neg=1
			}
		}
		st_numscalar("negative", neg)
		if (neg==1) {
			// Cameron, Gelbach, and Miller (2011) adjustment from Gu and Yoo (2019, DOI: 10.1177/1536867X19893637) 
			symeigensystem(st_matrix("Vc"), EVEC = ., eval = .)
			eval = eval :* (eval :> 0)
			st_matrix("Vcadj", EVEC*diag(eval)*EVEC')
		}
	}
    
end
