// This code is written to test the performance of bootstep.ado using simulations
// Version: 17dec2022 (David Veenman)

	clear all
	set seed 1234
	set obs 10000
	local firms=1000
	local industries=100
	
	local iter=1000
	local nboot=1000 
	
	gen firm=ceil(_n/(_N/`firms'))
	bysort firm: gen year=_n
	gen industry=ceil(_n/(_N/`industries'))
	sum firm year industry
	
	bysort industry: gen industryn=_n
		
	capture program drop panelframe
	program panelframe
		gen double xg=rnormal() if industryn==1
		gegen max=max(xg), by(industry)
		replace xg=max
		drop max
		gen double eg=rnormal() if industryn==1
		gegen max=max(eg), by(industry)
		replace eg=max
		drop max
		gen double yg=rnormal() if industryn==1
		gegen max=max(yg), by(industry)
		replace yg=max
		drop max
		gen double x1=xg+rnormal()
		gen double e1=eg+rnormal()
		gen double s=x1+e1
		gen double x2=rnormal()
		gen double e2=yg+rnormal()
		gen double y=x1+x2+e2	
		gen double x3=rnormal()
		gen double x4=rnormal()
		gen double x5=rnormal()
		drop e1 e2 xg eg yg
	end

	gen n=_n	
	gen b_ols=.
	gen b_bootstep=.
	gen se_ols=.
	gen se_bootstep=.
	gen sig_ols=0 if n<=`iter'
	gen sig_bootstep=0 if n<=`iter'
	
	timer clear
	forvalues i=1(1)`iter'{
		qui panelframe
		qui timer on 1
		qui reg s x1
		qui predict s_hat, xb
		qui reg y x2 x3 x4 x5 s_hat, cluster(industry)	
		qui replace b_ols=_b[s_hat] if n==`i'
		qui replace se_ols=_se[s_hat] if n==`i'
		qui local t=(_b[s_hat]-1)/_se[s_hat]
		qui local p=2*ttail(e(df_r), abs(`t'))
		qui replace sig_ols=1 if `p'<=0.05 & n==`i'			
		qui timer off 1
		qui timer on 2
		qui bootstep s x1 | y x2 x3 x4 x5, nboot(`nboot') cluster(industry) 
		qui replace b_bootstep=_b[s_hat] if n==`i'
		qui replace se_bootstep=_se[s_hat] if n==`i'
		qui local t=(_b[s_hat]-1)/_se[s_hat]
		qui local p=2*ttail(e(df_r), abs(`t'))
		qui replace sig_bootstep=1 if `p'<=0.05 & n==`i'			
		qui timer off 2
		drop x1 s x2 y x3 x4 x5 s_hat
		di `i'
		tabstat b_* se_* sig_* if n<=`i', stats(mean sd)
	}
	
	tabstat b_* se_* sig_*, stats(mean sd)
	
	timer list
	
	/*
.         tabstat b_* se_* sig_*, stats(mean sd)

   Stats |     b_ols  b_boot~p    se_ols  se_boo~p   sig_ols  sig_bo~p
---------+------------------------------------------------------------
    Mean |  1.001036  1.001036  .0499406  .0711964      .165      .045
      SD |  .0712532  .0712532  .0064542  .0101848  .3713663  .2074079
----------------------------------------------------------------------

.         
.         timer list
   1:     30.42 /     1000 =       0.0304
   2:   3457.77 /     1000 =       3.4578
	
