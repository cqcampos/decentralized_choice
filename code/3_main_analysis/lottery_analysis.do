/********************************************************************************
2_lottery_analysis.do
- Author: Chris Campos, Antonia Vazquez
- Description: Conducts lottery analysis for the decentralized choice project
- Date created: 07/10/2025 


* Installs: unique, ritest
*******************************************************************************/
clear all
set trace on 
set tracedepth 2
set maxvar 120000
*Packages needed
*ssc inst egenmore
*ssc install coefplot, replace
*ssc install seq
*ssc install unique
*net install grc1leg, from( http://www.stata.com/users/vwiggins/)

capture program drop main 
program define main 

	set_paths 

    * Table 2 
	*lottery_balance
    * Figure 5
	*distance_analysis
    * Figure 6
	*dist_hetero_robust

    *pref_hetero

	*peer_effects

	model_validate, eta("eta_out") homog("hetero") types(3) subject("math") blockfe("TRUE") app_2013(1)
	
	model_lates ,  eta("eta_out") homog("hetero") types(3) subject("math") blockfe("TRUE") app_2013(1)
end 



********************************************************************************
* set_paths 
* 
* Description: Establish paths 
********************************************************************************
capture program drop set_paths
program define set_paths

	if "`c(os)'"=="MacOSX" & "`c(username)'"=="cqcampos"{
		global ROOT "/Volumes/lausd/decentralized_choice"
		global BUILD "/Volumes/lausd/build/output/data"
    	global MAGNET "/Volumes/lausd/magnet"
		global RAWDATA "$BUILD/rawdata/Campos_Master_DUA" 
		global BUILDROOT "/Volumes/lausd/build/" 
	}
	else{
			global ROOT "/project/lausd/decentralized_choice"
			global BUILD "/project/lausd/build/output/data"
			global MAGNET "/project/lausd/magnet"
			global RAWDATA "$BUILD/rawdata/Campos_Master_DUA" 
			global BUILDROOT "/project/lausd/build/" 

	}

	* Data directories for Magnet project 
    global MAGNETDTA "$MAGNET/data/final"
    global MAGNETINT "$MAGNET/data/intermediate"
	global LOGS "$ROOT/logs"
	global OUT "$ROOT/output"

    * Data directories for decentralized choice folder 
	global DTARAW "$ROOT/data/raw" 
	global DTAVAM "$ROOT/rawdata/vam/data/intermediate"
	global DTAINT "$ROOT/data/intermediate" 
	global DTAFINAL "$ROOT/data/"

	* Where we save output 
	global tables "$ROOT/output/tables"
	global figures "$ROOT/output/figures/tempfigures"
	global slidefigures "$ROOT/output/slidefigures"

	* Old directories, to be updated 
	global LOGS "$ROOT/logs"

	* Raw source data
	global UE "$BUILDROOT/rawdata/Campos_Master_DUA/Unified_Enrollment"
end



********************************************************************************
* lottery_balance
*
* Description: Conducts lottery balance analysis 
********************************************************************************/
capture program drop lottery_balance
program define lottery_balance
	* ssc install ritest
	*seed for randomization inference
	set seed 2025 
	
	use "$DTAFINAL/lotteries.dta", clear
	keep if endyear <=2019

	ivreghdfe enrolled_app above_cutoff  ///
							if !missing(app_dist), absorb(lottery_ID ) cluster(lottery_ID)
	gen sample = 1 if e(sample)==1
	replace sample = 0 if e(sample)==0
	keep if sample==1
    gen rel_dist = app_nearest_mag_dist - app_nearest_schl_dist 
    local chars "female app_stu_black app_stu_latino english_learner special_edu poverty college_grad eng_home esp_home app_nearest_mag_dist L1numsuspensions L1mat L1ela anyScore_ms"

    gen ms_app_stu_black = missing(app_stu_black)
    gen ms_app_stu_latino = missing(app_stu_latino)
    gen ms_app_stu_asian = missing(app_stu_asian)


	texdoc init "$tables/lottery_balance.tex" , replace force 
	tex \begin{tabular}{lccc} \toprule \hline
	tex  & Control Mean & Above Cutoff Difference & P-Value \\ \hline \addlinespace
	
	* loop through each 
	foreach var of local chars{
		reghdfe `var' above_cutoff if ms_`var'==0, absorb(lottery_ID) cluster(lottery_ID)
		
		local coef_`var': dis %05.3f _b[above_cutoff]
		local se_`var': dis %05.3f _se[above_cutoff]
		local f_stat_`var': dis %05.3f e(F)
        sum `var' if ms_`var'==0 & above_cutoff==0 
        local mean_`var': dis %05.3f r(mean)
		
		test above_cutoff = 0
		local p_val_`var': dis %05.3f r(p)
		local vlabel : variable label `var'
		
		tex `vlabel' & `mean_`var'' & `coef_`var'' & `p_val_`var'' \\
		tex &  & (`se_`var'') & \\
		
	}

    *** Joint test using RI ***
	/*
	ritest above_cutoff e(F), reps(1000)  : reghdfe above_cutoff female english_learner special_edu poverty college_grad eng_home esp_home L1numsuspensions L1math L1ela app_nearest_mag_dist  ms_L1math ms_L1ela ms_L1numsuspensions   if inData==1, absorb(lottery_ID ) cluster(lottery_ID)
	
	*saving p value. exported as a matrix
	scalar p_temp = r(p)[1,1]  
	local pvalue_ri : display %05.3f p_temp
    */
	*** Joint test using SUR ***

	* due to so many fe, we residualize first, then do the sur stuff
	reghdfe above_cutoff if inData==1, absorb(lottery_ID) vce(cluster lottery_ID) residuals(r_above_cutoff)
	foreach var of local chars {
			reghdfe `var' if ms_`var'==0, absorb(lottery_ID) vce(cluster lottery_ID) residuals(r_`var')
	}
 

	foreach var of local chars {
		reg r_`var' r_above_cutoff
		estimates store m_`var'
		}
    local obs = e(N)

	suest m_*, cluster(lottery_ID)
	test r_above_cutoff
	local pvalue_sur: dis %05.3f r(p) 
	
	tex \addlinespace \hline \addlinespace
	tex Joint test p-value & \multicolumn{3}{c}{`pvalue_sur'} \\ 
    tex Observations & \multicolumn{3}{c}{`obs'} \\ \addlinespace \hline\bottomrule
	tex  \end{tabular}
	texdoc close 
	
end


capture program drop distance_analysis 
program define distance_analysis
    * ssc install ritest
    *seed for randomization inference
    set seed 2025 
    
    use "$DTAFINAL/lotteries.dta", clear

    ********************************* First Stage Figure *****************************************************
    * App distanace 
	reghdfe enrolled_app above_dist*  i.app_dist_rel_quint ///
							if  endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)

   
	test above_dist1=above_dist5
	local pval : di %9.2f round(r(p), .001)
	test above_dist1=above_dist2= above_dist3=above_dist4=above_dist5=0
	local pval_joint : di %9.2f round(r(p), .001)
	preserve
		parmest, norestore level(90)

		keep if _n<=5
		gen distquint = substr(parm, -1,1)
		destring distquint, replace
		twoway (line estimate distquint, color(black)) (rcap min90 max90 distquint, color(gs10) ), ///
			xtitle("Distance Quintile") ytitle("First Stage Effect") ///
			legend(off)  note("H0 Q5=Q1 p-value: `pval'"         "H0 Jointly Zero p-value:   `pval_joint'")
		graph export "$figures/dist_hetero_fs_app_dist_rel_quint.pdf", replace
	restore

    ********************************* First Stage Figure *****************************************************
    * relative distance 
	reghdfe enrolled_app above_rel_dist*  i.dist_rel_quint ///
							if  endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)

    test above_rel_dist1=above_rel_dist5
	local pval : di %9.2f round(r(p), .001)
	test above_rel_dist1=above_rel_dist2= above_rel_dist3=above_rel_dist4=above_rel_dist5=0
	local pval_joint : di %9.2f round(r(p), .001)
	preserve
		parmest, norestore level(90)

		keep if _n<=5
		gen distquint = substr(parm, -1,1)
		destring distquint, replace
		twoway (line estimate distquint, color(black)) (rcap min90 max90 distquint, color(gs10) ), ///
			xtitle("Distance Quintile") ytitle("First Stage Effect") ///
			legend(off)  note("H0 Q5=Q1 p-value: `pval'"         "H0 Jointly Zero p-value:   `pval_joint'")
		graph export "$figures/dist_hetero_fs_dist_rel_quint.pdf", replace
	restore  

    ********************************* Reduced Form Figure *****************************************************
    * Reduced form for app distance
    local exam_vars "math ela"
    foreach exam of local exam_vars {
        ivreghdfe avg_F`exam' above_dist*  ///
                      i.app_dist_rel_quint   ///
                        if  endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)

        
        test above_dist1=above_dist2= above_dist3=above_dist4=above_dist5
        local pval : di %9.2f round(r(p), .001)
        test above_dist1=above_dist2= above_dist3=above_dist4=above_dist5=0
        local pval_joint : di %9.2f round(r(p), .001)
        preserve 
        parmest, norestore level(90)

        keep if _n<=5
        gen distquart = substr(parm, -1,1)
        destring distquart, replace 
        twoway (line estimate distquart, color(black)) (rcap min90 max90 distquart, color(gs10) ), ///
            xtitle("Distance Quintile") ytitle("Reduced Form Effect") ///
            legend(off)  note("H0 Equal Effects p-value: `pval'"         "H0 Jointly Zero p-value:   `pval_joint'")
        graph export "$figures/dist_hetero_rf_`exam'_app_distance.pdf", replace 
        restore 

    }

    ********************************* Reduced Form Figure *****************************************************
    * Reduced form for relative distance
    local exam_vars "math ela"
    foreach exam of local exam_vars {
        ivreghdfe avg_F`exam' above_rel_dist*  ///
                     i.dist_rel_quint   ///
                        if  endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)

        test above_rel_dist1=above_rel_dist2= above_rel_dist3=above_rel_dist4=above_rel_dist5
        local pval : di %9.2f round(r(p), .001)
        test above_rel_dist1=above_rel_dist2= above_rel_dist3=above_rel_dist4=above_rel_dist5=0
        local pval_joint : di %9.2f round(r(p), .001)
        preserve 
        parmest, norestore level(90)

        keep if _n<=5
        gen distquart = substr(parm, -1,1)
        destring distquart, replace 
        twoway (line estimate distquart, color(black)) (rcap min90 max90 distquart, color(gs10) ), ///
            xtitle("Distance Quintile") ytitle("Reduced Form Effect") ///
            legend(off)  note("H0 Equal Effects p-value: `pval'"         "H0 Jointly Zero p-value:   `pval_joint'")
        graph export "$figures/dist_hetero_rf_`exam'_rel_distance.pdf", replace 
        restore 

    }


end 



 capture program drop dist_hetero_robust
 program define dist_hetero_robust

	syntax, [gradetype(string) dist_var(string) college(string) group(string) ]


    use "$DTAFINAL/lotteries.dta", clear

    ********************************** Make Figure *****************************************************
    * Restrict to cohorts on or before 2017 (enrollment year) because of F1-F3 average 
    local vars "lag_math lag_suspensions female poverty special_edu gifted english_learner eng_home esp_home born_usa"
    foreach var of local vars{
        tempfile temp_`var'
    }
    tempfile all

    ivreghdfe avg_Fmath above_dist*  i.app_dist_rel_quint ///
                    ///
                    if app_dist_rel_quint!=. & endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)
    nlcom (diff_all: _b[above_dist1] - _b[above_dist5]), post 
    parmest, saving(`all')
    foreach var of local vars{
        ivreghdfe avg_Fmath above_dist* c.above_cutoff#c.`var' `var' i.app_dist_rel_quint ///
                    if app_dist_rel_quint!=. & endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)
         nlcom (diff_`var': _b[above_dist1] - _b[above_dist5]), post
        parmest, saving(`temp_`var'')
    }
    preserve 
    use `all', clear 
    foreach var of local vars{
        append using `temp_`var'', force
    }
    split parm, parse("diff_")
    gen id = _n 
    twoway (scatter estimate id, color(black)) ///
            (rspike max95 min95 id, color(gs10)), ///
            legend(off) ylabel(0(0.05)0.2) ///
            xlabel( 0 " " 1 "Baseline" 2 "Achievement" ///
                    3 "Suspensions"  4 "Female" ///
                    5 "Poverty" 6 "Special Education" ///
                    7 "Gifted" 8 "English Learner" ///
                    9 "English at Home" 10 "Spanish at Home" ///
                    11 "Born USA" 12 " " , angle(45)) ///
            ytitle("Q1 - Q5 Difference") ///
            xtitle("")
    graph export $figures/dist_hetero_robust_math.pdf, replace
    restore 
  

    ivreghdfe avg_Fela above_dist*  i.app_dist_rel_quint ///
                    ///
                    if app_dist_rel_quint!=. & endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)
    nlcom (diff_all: _b[above_dist1] - _b[above_dist5]), post 
    parmest, saving(`all', replace)
    foreach var of local vars{
        ivreghdfe avg_Fela above_dist* c.above_cutoff#c.`var' `var' i.app_dist_rel_quint ///
                    if app_dist_rel_quint!=. & endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)
         nlcom (diff_`var': _b[above_dist1] - _b[above_dist5]), post
        parmest, saving(`temp_`var'', replace)
    }
    use `all', clear 
    foreach var of local vars{
        append using `temp_`var'', force
    }
    split parm, parse("diff_")
    gen id = _n 
    twoway (scatter estimate id, color(black)) ///
            (rspike max95 min95 id, color(gs10)), ///
            legend(off) ylabel(0(0.05)0.2) ///
            xlabel( 0 " " 1 "Baseline" 2 "Achievement" ///
                    3 "Suspensions"  4 "Female" ///
                    5 "Poverty" 6 "Special Education" ///
                    7 "Gifted" 8 "English Learner" ///
                    9 "English at Home" 10 "Spanish at Home" ///
                    11 "Born USA" 12 " " , angle(45)) ///
            ytitle("Q1 - Q5 Difference") ///
            xtitle("")
    graph export $figures/dist_hetero_robust_ela.pdf, replace

end 




 capture program drop pref_hetero
 program define pref_hetero
	syntax, [gradetype(string) dist_var(string) college(string) group(string) ]


	use "$DTAFINAL/lotteries.dta", clear
	local vars "lag_math lag_suspensions female poverty special_edu gifted english_learner eng_home esp_home born_usa"

	foreach var of local vars{
		egen app_`var' = mean(`var')
		bys lottery_ID: egen app_l_`var' = mean(`var')
		bys mag_cd school_yr: egen app_j_`var' = mean(`var')
		* create differences 
		gen stu_app_`var' = (`var' - app_`var')^2 
		gen stu_j_`var' = (`var' - app_j_`var')^2		
		gen stu_l_`var' = (`var' - app_l_`var')^2
	}
/*    local vars dist_rel
	foreach var of local vars{
		bys lottery_ID: egen app_l_`var' = mean(`var')
		bys mag_cd school_yr: egen app_j_`var' = mean(`var')
		* create differences 
		gen stu_app_`var' = (`var' - app_`var')^2 
		gen stu_j_`var' = (`var' - app_j_`var')^2		
		gen stu_l_`var' = (`var' - app_l_`var')^2
	}    
*/
	egen app_norm = rowtotal(stu_app_*)

	egen app_j_norm = rowtotal(stu_j_*)
	egen app_l_norm = rowtotal(stu_l_*)
	replace app_norm = sqrt(app_norm)
	replace app_j_norm = sqrt(app_j_norm)
	replace app_l_norm = sqrt(app_l_norm)
	egen norm_quart=xtile(app_j_norm ), n(5) by(endyear)
	gen above_norm1 = above_cutoff*(norm_quart==1)
	gen above_norm2 = above_cutoff*(norm_quart==2)
	gen above_norm3 = above_cutoff*(norm_quart==3)
	gen above_norm4 = above_cutoff*(norm_quart==4)
	gen above_norm5 = above_cutoff*(norm_quart==5)

	ivreghdfe avg_Fmath above_norm* i.norm_quart           if dist_rel_quint!=. & endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)

	test above_norm1=above_norm2= above_norm3=above_norm4=above_norm5
	local pval : di %9.2f round(r(p), .001)
	test above_norm1=above_norm2= above_norm3=above_norm4=above_norm5=0
	local pval_joint : di %9.2f round(r(p), .001)
	preserve 
	parmest, norestore level(90)

	keep if _n<=5
	gen normquart = substr(parm, -1,1)
	destring normquart, replace 
	twoway (line estimate normquart, color(black)) (rcap min90 max90 normquart, color(gs10) ), ///
		xtitle("Preference Index Quintile") ytitle("Reduced Form Effect") ///
		legend(off)  note("H0 Equal Effects p-value: `pval'"         "H0 Jointly Zero p-value:   `pval_joint'")
stop;
	graph export "$figures/pref_hetero_rf_math.pdf", replace
    restore

    ivreghdfe avg_Fela above_norm* i.norm_quart      if dist_rel_quint!=. & endyear<=2017, absorb(lottery_ID ) cluster(lottery_ID)
    	test above_norm1=above_norm2= above_norm3=above_norm4=above_norm5
	local pval : di %9.2f round(r(p), .001)
	test above_norm1=above_norm2= above_norm3=above_norm4=above_norm5=0
	local pval_joint : di %9.2f round(r(p), .001)
	preserve 
	parmest, norestore level(90)

	keep if _n<=5
	gen normquart = substr(parm, -1,1)
	destring normquart, replace 
	twoway (line estimate normquart, color(black)) (rcap min90 max90 normquart, color(gs10) ), ///
		xtitle("Preference Index Quintile") ytitle("Reduced Form Effect") ///
		legend(off)  note("H0 Equal Effects p-value: `pval'"         "H0 Jointly Zero p-value:   `pval_joint'")
	graph export "$figures/pref_hetero_rf_ela.pdf", replace
    restore
end


capture program drop peer_effects
program define peer_effects 


	use "$DTAFINAL/lotteries.dta", clear
	bys mag_cd: gen cellsize = _N

	* Drop tiny cells 
	keep if cellsize >=20

	* make lottery offers real quick
	levelsof mag_cd , local(lotteries) clean
	foreach l of local lotteries{
		gen offer_l`l' = (mag_cd==`l')*(above_cutoff==1)
	}

	* first ensure we use same sample across both regs
	gen math_sample = (!missing(sch_z_math_all_enr)&!missing(avg_Fmath)&endyear<=2017)
	gen ela_sample = (!missing(sch_z_ela_all_enr)&!missing(avg_Fela)&endyear<=2017)
	
	***** math peer effect *****
	reghdfe sch_z_math_all_enr offer_l*  lag_ela lag_math /// 
		ms_lag_ela ms_lag_math lag_suspensions special_edu gifted female if math_sample , absorb(lottery_ID) vce(robust)
	preserve
		parmest, norestore
		drop if strpos(parm, "o.offer_")
		keep if strpos(parm, "offer_")
		keep if stderr > 0
		rename estimate math_pq_coef
		rename stderr math_pq_se
		rename parm lottery
		replace lottery = substr(lottery, 8, .)
		keep lottery math_pq_coef math_pq_se
		tempfile math_pq
		save `math_pq'
	restore

	***** math achievement effect *****
	reghdfe avg_Fmath offer_l*  lag_ela lag_math  ///
		 ms_lag_ela ms_lag_math lag_suspensions special_edu gifted female  if math_sample , absorb(lottery_ID) vce(robust)
	preserve
		parmest, norestore 
		drop if strpos(parm, "o.offer_")
		keep if strpos(parm, "offer_")
		keep if stderr > 0
		rename estimate math_ae_coef
		rename stderr math_ae_se
		rename parm lottery
		replace lottery = substr(lottery, 8, .)
		keep lottery math_ae_coef math_ae_se
		tempfile math_ae
		save `math_ae'
	restore
	
	***** ela peer effect *****
	reghdfe sch_z_ela_all_enr offer_l*   lag_ela lag_math  ///
		ms_lag_ela ms_lag_math lag_suspensions special_edu gifted female if ela_sample, absorb(lottery_ID) vce(robust)
	preserve
		parmest, norestore
		drop if strpos(parm, "o.offer_")
		keep if strpos(parm, "offer_")
		rename estimate ela_pq_coef
		rename stderr ela_pq_se
		rename parm lottery
		replace lottery = substr(lottery, 8, .)
		keep lottery ela_pq_coef ela_pq_se
		tempfile ela_pq
		save `ela_pq'
	restore

	***** ela achievement effect *****
	reghdfe avg_Fmath offer_l*  lag_ela lag_math  ///
		ms_lag_ela ms_lag_math lag_suspensions special_edu gifted female if ela_sample, absorb(lottery_ID) vce(robust)
	parmest, norestore 
	drop if strpos(parm, "o.offer_")
	keep if strpos(parm, "offer_")
	keep if stderr > 0
	rename estimate ela_ae_coef
	rename stderr ela_ae_se
	rename parm lottery
	replace lottery = substr(lottery, 8, .)
	keep lottery ela_ae_coef ela_ae_se
	tempfile ela_ae
	save `ela_ae'
	
	* merge all together to plot
	merge 1:1 lottery using `math_pq', keep(1 3) nogen
	merge 1:1 lottery using `math_ae', keep(1 3) nogen
	merge 1:1 lottery using `ela_pq', keep(1 3) nogen
	
	* Keep only schools with complete data for at least one outcome
	keep if (math_pq_coef != . & math_ae_coef != .) | (ela_pq_coef != . & ela_ae_coef != .)
	
	reg math_ae_coef math_pq_coef [aweight=1/math_ae_se^2], r
	local slope_math: dis %05.3f _b[math_pq_coef]
	local se_math: dis %05.3f _se[math_pq_coef]
	
	reg ela_ae_coef ela_pq_coef [aweight=1/ela_ae_se^2], r
	local slope_ela: dis %05.3f _b[ela_pq_coef]
	local se_ela: dis %05.3f _se[ela_pq_coef]

	* Math VIV Plot
	twoway (scatter math_ae_coef math_pq_coef  [aweight=1/math_ae_se], msymbol(oh) mcolor(black)) ///
	    (lfit math_ae_coef math_pq_coef [aweight=1/math_ae_se], lcolor(maroon) lpattern(dash)), ///
	    xtitle("Peer Quality") ///
	    ytitle("Math Achievement Effect") ///
	    name(viv_math, replace) legend(off) note("Slope: `slope_math' (`se_math')")
	graph export "$figures/viv_math.pdf", replace
	
	* ela VIV Plot
	twoway (scatter ela_ae_coef ela_pq_coef [aweight=1/ela_ae_se], msymbol(oh) mcolor(black)) ///
	    (lfit ela_ae_coef ela_pq_coef [aweight=1/ela_ae_se], lcolor(maroon) lpattern(dash)), ///
	    xtitle("Peer Quality") ///
	    ytitle("ELA Achievement Effect") ///
	    name(viv_ela, replace) legend(off) note("Slope: `slope_ela' (`se_ela')")
	graph export "$figures/viv_ela.pdf", replace


end 


capture program drop model_validate 
program define model_validate 
syntax, [eta(string) homog(string) types(int 1) subject(string) blockfe(string) app_2013(int 0)]

    * Read in demand model estimates to define preference index below 
    local last_yr = cond(`app_2013', 2013, 2008)
    local app_suffix = cond(`app_2013', "_2013", "")
	use $ESTIMATES/y1y0_`eta'_`homog'_K`types'_math_blockfe`blockfe'`app_suffix'.dta, clear 
	destring studentpseudoid, replace 
	rename magnet magnet_enroll_date 
	tempfile temp 
	save `temp'

	use $ESTIMATES/y1y0_`eta'_`homog'_K`types'_ela_blockfe`blockfe'`app_suffix'.dta, clear 
	destring studentpseudoid, replace 
	drop mu_theta magnet theta_percentile overall_percentile 
	rename (y0 y1 y delta matcheffect) (y0_ela y1_ela y_ela delta_ela matcheffect_ela)
	rename y1_novam y1_novam_ela
	rename yhat yhat_ela
	rename yhat_vam yhat_vam_ela
	tempfile tempela 
	save `tempela'
	use "$DTAFINAL/lotteries.dta", clear
	bys mag_cd: gen cellsize = _N

	* Drop tiny cells 
	*keep if cellsize >=20
	drop cellsize 
	bys lottery_ID: gen cellsize = _N
	keep if cellsize >=100

	* make lottery offers real quick
	levelsof lottery_ID , local(lotteries) clean
	foreach l of local lotteries{
		gen offer_l`l' = (lottery_ID==`l')*(above_cutoff==1)
		sum offer_l`l' if avg_Fmath!=. 
		if r(mean)==0 | r(mean)==.{
			drop offer_l`l'
		}
	}
	merge 1:1 studentpseudoid endyear using `temp', keep(1 3) nogen
	merge 1:1 studentpseudoid endyear using `tempela', keep(1 3) nogen
	rename offer_magnet out_offer 
	ivreg2 avg_Fmath (y = offer_*) i.lottery_ID , partial(i.lottery_ID) robust
	stop;
	gen sample_math = e(sample)
	tempfile math_rf math_fs
	reghdfe avg_Fmath offer_* if sample_math==1, absorb(lottery_ID) cluster(lottery_ID)
	parmest, saving(`math_rf', replace)
	reghdfe y offer_* if sample_math==1, absorb(lottery_ID) cluster(lottery_ID)
	parmest, saving(`math_fs', replace)
	tempfile ela_rf ela_fs
	ivreg2 avg_Fela (y_ela = offer_*) i.lottery_ID , partial(i.lottery_ID) robust
	gen sample_ela = e(sample)
	reghdfe avg_Fela offer_* if sample_ela==1, absorb(lottery_ID) cluster(lottery_ID)
	parmest, saving(`ela_rf', replace)
	reghdfe y_ela offer_* if sample_ela==1, absorb(lottery_ID) cluster(lottery_ID)
	parmest, saving(`ela_fs', replace)

	use `math_rf', clear
	drop if regexm(parm, "o.offer_") | regexm(parm, "cons")
	keep estimate parm stderr 
	gen subject = "math"
	gen type = "rf"
	save `math_rf', replace
	use `math_fs', clear
	drop if regexm(parm, "o.offer_") | regexm(parm, "cons")
	keep estimate parm stderr 
	gen subject = "math"
	gen type = "fs"
	save `math_fs', replace
	use `ela_rf', clear
	drop if regexm(parm, "o.offer_") | regexm(parm, "cons")
	keep estimate parm stderr 
	gen subject = "ela"
	gen type = "rf"
	save `ela_rf', replace
	use `ela_fs', clear
	drop if regexm(parm, "o.offer_") | regexm(parm, "cons")
	keep estimate parm stderr 
	gen subject = "ela"
	gen type = "fs"
	save `ela_fs', replace

	use `math_rf', clear
	append using `math_fs'
	append using `ela_rf'
	append using `ela_fs'
	stop;
end 

capture program drop model_lates 
program define model_lates 
syntax, [eta(string) homog(string) types(int 1) subject(string) blockfe(string) app_2013(int 0)]

    * Read in demand model estimates to define preference index below 
    local last_yr = cond(`app_2013', 2013, 2008)
    local app_suffix = cond(`app_2013', "_2013", "")
	use $ESTIMATES/y1y0_`eta'_`homog'_K`types'_math_blockfe`blockfe'`app_suffix'.dta, clear 
	destring studentpseudoid, replace 
	rename magnet magnet_enroll_date 
	tempfile temp 
	save `temp'

	use $ESTIMATES/y1y0_`eta'_`homog'_K`types'_ela_blockfe`blockfe'`app_suffix'.dta, clear 
	destring studentpseudoid, replace 
	drop mu_theta magnet theta_percentile overall_percentile 
	rename (y0 y1 y delta matcheffect) (y0_ela y1_ela y_ela delta_ela matcheffect_ela)
	rename y1_novam y1_novam_ela
	rename yhat yhat_ela
	tempfile tempela 
	save `tempela'
	use "$DTAFINAL/lotteries.dta", clear
	bys mag_cd: gen cellsize = _N

	* Drop tiny cells 
	*keep if cellsize >=20
	*drop cellsize 
	*bys lottery_ID: gen cellsize = _N
	*keep if cellsize >=20

	merge 1:1 studentpseudoid endyear using `temp', keep(1 3) nogen
	merge 1:1 studentpseudoid endyear using `tempela', keep(1 3) nogen
	gen low_score = lag_math <0 
	gen high_score = lag_math >=0

	gen nonpoverty = poverty==0
	local demos "app_stu_asian app_stu_black app_stu_latino app_stu_white female english_learner eng_home born_usa college_grad esp_home poverty nonpoverty low_score high_score"

	* Create a matrix to store IV and RF estimates for each covariate group above
	matrix results = J(14, 16, .)
	local i = 1
	foreach d of local demos {
		* Math RF
		ivreghdfe avg_Fmath above_cutoff   if `d'==1 & endyear<=2013, cluster( lottery_ID) absorb(lottery_ID)
		local rf_math: dis %05.3f _b[above_cutoff]
		local rf_math_se: dis %05.3f _se[above_cutoff]
		* Math IV
		ivreghdfe avg_Fmath (magnet_enroll_date = above_cutoff) if `d'==1 & endyear<=2013, absorb(lottery_ID) cluster(lottery_ID)
		local iv_math: dis %05.3f _b[magnet_enroll_date]
		local iv_math_se: dis %05.3f _se[magnet_enroll_date]
		* ELA RF
		ivreghdfe avg_Fela above_cutoff   if `d'==1 & endyear<=2013, cluster(lottery_ID) absorb(lottery_ID)
		local rf_ela: dis %05.3f _b[above_cutoff]
		local rf_ela_se: dis %05.3f _se[above_cutoff]
		* ELA IV
		ivreghdfe avg_Fela (magnet_enroll_date = above_cutoff)  if `d'==1 & endyear<=2013, absorb(lottery_ID) cluster(lottery_ID)
		local iv_ela: dis %05.3f _b[magnet_enroll_date]
		local iv_ela_se: dis %05.3f _se[magnet_enroll_date]
		* Store in matrix
		matrix results[`i', 1] = `rf_math'
		matrix results[`i', 2] = `rf_math_se'
		matrix results[`i', 3] = `iv_math'
		matrix results[`i', 4] = `iv_math_se'
		matrix results[`i', 5] = `rf_ela'
		matrix results[`i', 6] = `rf_ela_se'
		matrix results[`i', 7] = `iv_ela'
		matrix results[`i', 8] = `iv_ela_se'

		* now estimate again on the model-based outcomes 
		* Math RF
		ivreghdfe yhat above_cutoff   if `d'==1 & endyear<=2013, cluster( lottery_ID) absorb(lottery_ID)
		local rf_math: dis %05.3f _b[above_cutoff]
		local rf_math_se: dis %05.3f _se[above_cutoff]
		* Math IV
		ivreghdfe yhat (magnet_enroll_date = above_cutoff) if `d'==1 & endyear<=2013, absorb(lottery_ID) cluster(lottery_ID)
		local iv_math: dis %05.3f _b[magnet_enroll_date]
		local iv_math_se: dis %05.3f _se[magnet_enroll_date]
		* ELA RF
		ivreghdfe yhat_ela above_cutoff   if `d'==1 & endyear<=2013, cluster(lottery_ID) absorb(lottery_ID)
		local rf_ela: dis %05.3f _b[above_cutoff]
		local rf_ela_se: dis %05.3f _se[above_cutoff]
		* ELA IV
		ivreghdfe yhat_ela (magnet_enroll_date = above_cutoff)  if `d'==1 & endyear<=2013, absorb(lottery_ID) cluster(lottery_ID)
		local iv_ela: dis %05.3f _b[magnet_enroll_date]
		local iv_ela_se: dis %05.3f _se[magnet_enroll_date]

		* Store in matrix
		matrix results[`i', 9] = `rf_math'
		matrix results[`i', 10] = `rf_math_se'
		matrix results[`i', 11] = `iv_math'
		matrix results[`i', 12] = `iv_math_se'
		matrix results[`i', 13] = `rf_ela'
		matrix results[`i', 14] = `rf_ela_se'
		matrix results[`i', 15] = `iv_ela'
		matrix results[`i', 16] = `iv_ela_se'

		local i = `i' + 1
	}
	matrix list results 

end 
main 