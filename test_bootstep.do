// This code is written to test the performance of bootstep.ado
// Version: 16dec2022 (David Veenman)

	clear all
	set seed 1234
	set obs 10000
	local firms=1000
	local industries=100
	
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
	
	panelframe
	
	// Basic illustration of effects on standard errors and speed of the program:
	timer clear
	timer on 1
	reg s x1
	predict s_hat, xb
	reg y x2 x3 x4 x5 s_hat, cluster(industry)	
	timer off 1
	timer on 2
	bootstep s x1 | y x2 x3 x4 x5, nboot(1000) cluster(industry) seed(1234)
	timer off 2
	timer list
	
	// Twoway clustering:
	bootstep s x1 | y x2 x3 x4 x5, nboot(1000) cluster(industry year) seed(1234)
	
	// No clustering:
	bootstep s x1 | y x2 x3 x4 x5, nboot(1000) seed(1234)
	
	// Absorbing fixed effects in step 2:
	bootstep s x1 | y x2 x3 x4 x5, nboot(1000) cluster(industry) absorb(firm) seed(1234)
	
	// Absorbing fixed effects in step 2 and including factor variables in steps 2 (and/or 1):
	bootstep s x1 | y i.year x2 x3 x4 x5, nboot(1000) cluster(industry) absorb(firm) seed(1234)
	
	// Standardize coefficient on generated regressor:
	bootstep s x1 | y x2 x3 x4 x5, nboot(1000) cluster(industry) standardize seed(1234)
	
	