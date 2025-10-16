/********************************************************************************
3_outcome_analysis.do
- Author: Chris Campos 
- Description: Conducts lottery and outcome analysis for the decentralized choice project
- Date created: 08/23/2025 


* Installs: unique, ritest
*******************************************************************************/
clear all
set trace on 
set tracedepth 2
set maxvar 32000
*Packages needed
*ssc inst egenmore
*ssc install coefplot, replace
*ssc install seq
*ssc install unique
*net install grc1leg, from( http://www.stata.com/users/vwiggins/)

capture program drop main 
program define main 

	set_paths 

    * Event study 
*    event_study, type("tract") relative("no")
*    event_study, type("tract") relative("yes")
    * stop;
    *event_study, type("block")



    * Table 4 
    forval k = 3/3{
    	foreach c in  "eta_out" {
		foreach sbjt in "math" "ela"{
			*demand_for_effectiveness, eta("`c'") homog("hetero") types(`k') subject("`sbjt'") app_2013(1)
		
			
		 	cf_regressions, eta("`c'") homog("hetero") types(`k') subject("`sbjt'") blockfe("TRUE") app_2013(1)
            stop;
		}
		*cf_regressions_noncog, eta("`c'") homog("hetero") types(`k')  blockfe("TRUE") app_2013(1)

		*demand_for_effectiveness_noncog, eta("`c'") homog("hetero") types(`k') app_2013(1)
        
		cf_table2, eta("`c'") homog("hetero") types(`k') blockfe("TRUE") app_2013(1)
	}
    }
    *counterfactuals
stop;

    *lottery_model_hetero
	
	*maxAllocation
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
    global ESTIMATES "$ROOT/estimates"

	* Where we save output 
	global tables "$ROOT/output/tables"
	global figures "$ROOT/output/figures"
	global slidefigures "$ROOT/output/slidefigures"

	* Old directories, to be updated 
	global LOGS "$ROOT/logs"

	* Raw source data
	global UE "$BUILDROOT/rawdata/Campos_Master_DUA/Unified_Enrollment"
end

capture program drop event_study 
program define event_study 
syntax, type(string) [relative(string)]

    use $DTAFINAL/event_study_data_`type'_relative`relative'.dta, clear 
    if "`type'"=="tract" local geoid censustractid 
    if "`type'"=="block" local geoid censusblockid

    * Compositional changes 
    reghdfe app_yhat_wt         b3_d2 b2_d2  a*_d2 if stack<=2013 , absorb(stack_group stack_yr ) vce(cluster `geoid' )  

      preserve
    parmest , norestore
    expand 2 in 1
    split parm, parse("_")
    drop if parm=="_cons"
    replace parm1 = "b1" if _n==_N
    replace estimate = 0 if parm1=="b1"
    replace min95 = 0 if parm1=="b1"
    replace max95 = 0 if parm1=="b1"

    drop t
    gen t = substr(parm, 2,1)
    destring t , replace
    replace t = -t if substr(parm1, 1,1)=="b"
    replace t = -1 if parm1=="b1"
    sort t
    twoway (connected estimate t, color(black) ) (rcap min95 max95 t, color(gs10)), legend(off) ytitle("Change in Composition Index") xtitle("Time Relative to Entry") ylabel(-0.025(0.025)0.025)
    graph export $figures/event_study_app_yhat_wt_`type'_relative`relative'.pdf, replace
    restore  

    reghdfe appliedMag         b3_d2 b2_d2  a*_d2 if stack<=2013 , absorb(stack_group stack_yr ) vce(cluster `geoid' )  

    preserve
    parmest , norestore
    expand 2 in 1
    split parm, parse("_")
    drop if parm=="_cons"
    replace parm1 = "b1" if _n==_N
    replace estimate = 0 if parm1=="b1"
    replace min95 = 0 if parm1=="b1"
    replace max95 = 0 if parm1=="b1"

    drop t
    gen t = substr(parm, 2,1)
    destring t , replace
    replace t = -t if substr(parm1, 1,1)=="b"
    replace t = -1 if parm1=="b1"
    sort t
    twoway (connected estimate t, color(black) ) (rcap min95 max95 t, color(gs10)), legend(off) ytitle("Change in Magnet Applications") xtitle("Time Relative to Entry") 
    graph export $figures/event_study_apps_`type'_relative`relative'.pdf, replace
    restore

    if("`type'"=="tract"){

        reghdfe lnapplied         b3_d2 b2_d2  a*_d2 if stack<=2013 , absorb(stack_group stack_yr ) vce(cluster `geoid' )  

        preserve
        parmest , norestore
        expand 2 in 1
        split parm, parse("_")
        drop if parm=="_cons"
        replace parm1 = "b1" if _n==_N
        replace estimate = 0 if parm1=="b1"
        replace min95 = 0 if parm1=="b1"
        replace max95 = 0 if parm1=="b1"

        drop t
        gen t = substr(parm, 2,1)
        destring t , replace
        replace t = -t if substr(parm1, 1,1)=="b"
        replace t = -1 if parm1=="b1"
        sort t
        twoway (connected estimate t, color(black) ) (rcap min95 max95 t, color(gs10)), legend(off) ytitle("Change in Log Magnet Applications") xtitle("Time Relative to Entry")
        graph export $figures/event_study_lnapps_base_`type'_relative`relative'.pdf, replace
        restore
        reghdfe lnapplied         b3_d b2_d  a*_d if stack<=2013 , absorb(stack_group stack_yr ) vce(cluster `geoid' )  

        preserve
        parmest , norestore
        expand 2 in 1
        split parm, parse("_")
        drop if parm=="_cons"
        replace parm1 = "b1" if _n==_N
        replace estimate = 0 if parm1=="b1"
        replace min95 = 0 if parm1=="b1"
        replace max95 = 0 if parm1=="b1"

        drop t
        gen t = substr(parm, 2,1)
        destring t , replace
        replace t = -t if substr(parm1, 1,1)=="b"
        replace t = -1 if parm1=="b1"
        sort t
        twoway (connected estimate t, color(black) ) (rcap min95 max95 t, color(gs10)), legend(off) ytitle("Change in Log Magnet Applications") xtitle("Time Relative to Entry")
        graph export $figures/event_study_lnapps_`type'_relative`relative'.pdf, replace
        restore
    }
    reghdfe appliedMag         b3_d b2_d  a*_d if stack<=2013 , absorb(stack_group stack_yr ) vce(cluster `geoid' )  

    preserve
    parmest , norestore
    expand 2 in 1
    split parm, parse("_")
    drop if parm=="_cons"
    replace parm1 = "b1" if _n==_N
    replace estimate = 0 if parm1=="b1"
    replace min95 = 0 if parm1=="b1"
    replace max95 = 0 if parm1=="b1"

    drop t
    gen t = substr(parm, 2,1)
    destring t , replace
    replace t = -t if substr(parm1, 1,1)=="b"
    replace t = -1 if parm1=="b1"
    sort t
    twoway (connected estimate t, color(black) ) (rcap min95 max95 t, color(gs10)), legend(off) ytitle("Change in Magnet Applications") xtitle("Time Relative to Entry")
    graph export $figures/event_study_apps_base_`type'_relative`relative'.pdf, replace
    restore


    * Make a table showing the estimates of the subindex 
    gen pre_treat = (t <-1)*(treat==1)
    gen post_treat = (t>-1)* (treat==1)
    local vars "app_math_wt app_ela_wt app_black_wt app_hispanic_wt app_white_wt app_eng_home_wt app_suspensions_wt  app_poverty_wt app_yhat_wt"
	texdoc init "$tables/event_study_composition_table_relative`relative'.tex", replace force
	tex \begin{tabular}{lccc} \toprule \hline
    tex & Pre & Post \\ 
    tex & (1) & (2)  \\ \hline \hline 
    tex & &  \\
    foreach var of local vars{ 
        if "`var'"=="app_yhat_wt" local label = "Composition Index"
        if "`var'"=="app_math_wt" local label = "Math Score"
        if "`var'"=="app_ela_wt" local label = "ELA Score"
        if "`var'"=="app_black_wt" local label = "Black Share"
        if "`var'"=="app_hispanic_wt" local label = "Hispanic Share"
        if "`var'"=="app_white_wt" local label = "White Share"
        if "`var'"=="app_eng_home_wt" local label = "English Home Share"
        if "`var'"=="app_suspensions_wt" local label = "Suspension Rate"
        if "`var'"=="app_poverty_wt" local label = "Poverty Rate"

        reghdfe `var' pre_treat post_treat if stack<=2013 , absorb(stack_group stack_yr ) vce(cluster `geoid' )
        local b : di  %9.3f _b[pre_treat]  
        local se : di  %9.3f _se[pre_treat]
        local b2 : di  %9.3f _b[post_treat]
        local se2 : di  %9.3f _se[post_treat]
        texdoc write  `label' &  `b'  &  `b2' \\ 
        texdoc write  & (`se') & (`se2') \\ 

    }   
    count if e(sample)==1
    local n = r(N)
    texdoc write & & \\ \hline 
    texdoc write N & \multicolumn{2}{c}{`n'} \\ 
    texdoc write \hline \bottomrule
    texdoc write \end{tabular}
    texdoc close

end


capture program drop demand_for_effectiveness 
program define demand_for_effectiveness
    syntax, [eta(string) homog(string) types(int 1) subject(string) app_2013(int 0)]

    * Read in demand model estimates to define preference index below 
    local last_yr = cond(`app_2013', 2013, 2008)
    local app_suffix = cond(`app_2013', "_2013", "")
    import delimited "$ESTIMATES/in_mag_ind_`eta'_estimates_choice_`homog'_K`types'`app_suffix'.csv", clear
    sum omega if names == "Distance Cost"
    local distance_cost = r(mean)
    
    local vars "female black white hispanic poverty el eng_home lag_math lag_ela median_income curr"
    sum omega if names =="U x Female"
    local p_female = r(mean)
    sum omega if names =="U x Black"
    local p_black = r(mean)
    sum omega if names =="U x White"
    local p_white = r(mean)
    sum omega if names =="U x Hispanic"
    local p_hispanic = r(mean)
    sum omega if names =="U x Poverty"
    local p_poverty = r(mean)
    sum omega if names =="U x LEP"
    local p_el = r(mean)
    sum omega if names =="U x English Home"
    local p_eng_home = r(mean)
    sum omega if names =="U x Lag Math"
    local p_lag_math = r(mean)
    sum omega if names =="U x Lag ELA"
    local p_lag_ela = r(mean)
    sum omega if names =="U x Income"
    local p_median_income = r(mean)
    sum omega if names =="U x In Magnet"
    local p_in_magnet = r(mean)
    sum omega if names =="U x Peer Q"
    local p_peer_q = r(mean)
    
    keep if regexm(names, "Mean Utility \d+")
    
    rename v1 sid 
    gen wtt = (-1)*(omega/`distance_cost')
    tempfile meanUtilities 
    save `meanUtilities', replace

 * Read in posterior mean estimates 
    import delimited using "$ESTIMATES/in_mag_ind_`eta'_posteriors_`homog'_K`types'`app_suffix'.csv", clear stringcols(3)

    drop v1 
    rename (v4 v5) (mu_theta var_theta)
    isid studentpseudoid endyear 
  
    tempfile posteriors 
    save `posteriors'


    use "$DTAFINAL/structural_data_2004_`last_yr'.dta", clear 
        rename current_in_mag in_magnet 
    rename current_peer_quality peer_q 

    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela peer_q"
    merge m:1 studentpseudoid endyear using `posteriors', keepusing(mu_theta var_theta) gen(mergePosteriors) keep(1 3)
    gen magnet_enroll = D_0==0
    foreach var of local vars{
        gen magnet_`var' = magnet_enroll*`var'
    }
    gen theta_magnet = mu_theta * magnet_enroll
    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela in_magnet peer_q"

    foreach var of local vars{
        egen meantemp = mean(`var')
        gen orig_`var' = `var'
        replace `var' = `var' - meantemp
        drop meantemp
    }

    rename magnet_offer offer_var
    tab endyear, gen(yrdum)
    tab localdistrictcode, gen(districtdum)
    rename yrdum1 outyr 
    rename districtdum1 outdist

    rename D_0 neighborhood_school 


    local vars "female black white hispanic asian poverty el eng_home lag_math lag_ela peer_q"
    reghdfe  avg_F`subject'   D_* `vars'   magnet_*  ///  
        mu_theta  theta_magnet  yrdum*  districtdum*, absorb(censusblockid, savefe)  vce(cluster schoolcostcentercode)


    preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
    restore 

    use `mag_data', clear
    gen sid = _n 
    merge 1:1 sid using `meanUtilities', gen(mergeMeanUtilities) keep(1 3)

    twoway (scatter   estimate wtt, mcolor(maroon) lcolor(black) )  ///
        (lfit estimate wtt, lcolor(gs10) lwidth(medium)) , ///
        ytitle("Average Treatment Effect") xtitle("Willingness to Travel") ///
        legend(off) 

    graph export "$figures/demand_for_effectiveness_`eta'_`homog'_K`types'_`subject'`app_suffix'.pdf", replace
end 



capture program drop demand_for_effectiveness_noncog 
program define demand_for_effectiveness_noncog
    syntax, [eta(string) homog(string) types(int 1) subject(string) app_2013(int 0)]

    * Read in demand model estimates to define preference index below 
    local last_yr = cond(`app_2013', 2013, 2008)
    local app_suffix = cond(`app_2013', "_2013", "")
    import delimited "$ESTIMATES/in_mag_ind_`eta'_estimates_choice_`homog'_K`types'`app_suffix'.csv", clear
    sum omega if names == "Distance Cost"
    local distance_cost = r(mean)
    
    local vars "female black white hispanic poverty el eng_home lag_math lag_ela median_income curr"
    sum omega if names =="U x Female"
    local p_female = r(mean)
    sum omega if names =="U x Black"
    local p_black = r(mean)
    sum omega if names =="U x White"
    local p_white = r(mean)
    sum omega if names =="U x Hispanic"
    local p_hispanic = r(mean)
    sum omega if names =="U x Poverty"
    local p_poverty = r(mean)
    sum omega if names =="U x LEP"
    local p_el = r(mean)
    sum omega if names =="U x English Home"
    local p_eng_home = r(mean)
    sum omega if names =="U x Lag Math"
    local p_lag_math = r(mean)
    sum omega if names =="U x Lag ELA"
    local p_lag_ela = r(mean)
    sum omega if names =="U x Income"
    local p_median_income = r(mean)
    sum omega if names =="U x In Magnet"
    local p_in_magnet = r(mean)
    sum omega if names =="U x Peer Q"
    local p_peer_q = r(mean)
    
    keep if regexm(names, "Mean Utility \d+")
    
    rename v1 sid 
    gen wtt = (-1)*(omega/`distance_cost')
    tempfile meanUtilities 
    save `meanUtilities', replace

 * Read in posterior mean estimates 
    import delimited using "$ESTIMATES/in_mag_ind_`eta'_posteriors_`homog'_K`types'`app_suffix'.csv", clear stringcols(3)

    drop v1 
    rename (v4 v5) (mu_theta var_theta)
    isid studentpseudoid endyear 
  
    tempfile posteriors 
    save `posteriors'


    use "$DTAFINAL/structural_data_2004_`last_yr'.dta", clear 
        rename current_in_mag in_magnet 
    rename current_peer_quality peer_q 

    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela peer_q"
    merge m:1 studentpseudoid endyear using `posteriors', keepusing(mu_theta var_theta) gen(mergePosteriors) keep(1 3)
    gen magnet_enroll = D_0==0
    foreach var of local vars{
        gen magnet_`var' = magnet_enroll*`var'
    }
    gen theta_magnet = mu_theta * magnet_enroll
    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela in_magnet peer_q"

    foreach var of local vars{
        egen meantemp = mean(`var')
        gen orig_`var' = `var'
        replace `var' = `var' - meantemp
        drop meantemp
    }

    rename magnet_offer offer_var
    tab endyear, gen(yrdum)
    tab localdistrictcode, gen(districtdum)
    rename yrdum1 outyr 
    rename districtdum1 outdist

    rename D_0 neighborhood_school 


    local vars "female black white hispanic asian poverty el eng_home lag_math lag_ela peer_q"

    reghdfe  z_socio   D_* `vars'   magnet_*  ///  
        mu_theta  theta_magnet  yrdum*  districtdum*, absorb(censusblockid, savefe)  vce(cluster schoolcostcentercode)


    preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
    restore 

    use `mag_data', clear
    gen sid = _n 
    merge 1:1 sid using `meanUtilities', gen(mergeMeanUtilities) keep(1 3)

    twoway (scatter   estimate wtt, mcolor(maroon) lcolor(black) )  ///
        (lfit estimate wtt, lcolor(gs10) lwidth(medium)) , ///
        ytitle("Average Treatment Effect") xtitle("Willingness to Travel") ///
        legend(off) 

    graph export "$figures/demand_for_effectiveness_noncog_`eta'_`homog'_K`types'_`app_suffix'.pdf", replace
end 


********************************************************************************
* cf_regressions 
********************************************************************************/
capture program drop cf_regressions
program define cf_regressions
syntax, [eta(string) homog(string) types(int 1) subject(string) blockfe(string) app_2013(int 0)]

  
    
    ********************************* Prepare some inputs for the analysis *********************************

    local last_yr = cond(`app_2013', 2013, 2008)
    local app_suffix = cond(`app_2013', "_2013", "")
    
    * Read in demand model estimates to define preference index below 
    import delimited "$ESTIMATES/in_mag_ind_`eta'_estimates_choice_`homog'_K`types'`app_suffix'.csv", clear
    local vars "female black white hispanic poverty el eng_home lag_math lag_ela median_income curr"
    sum omega if names =="U x Female"
    local p_female = r(mean)
    sum omega if names =="U x Black"
    local p_black = r(mean)
    sum omega if names =="U x White"
    local p_white = r(mean)
    sum omega if names =="U x Hispanic"
    local p_hispanic = r(mean)
    sum omega if names =="U x Poverty"
    local p_poverty = r(mean)
    sum omega if names =="U x LEP"
    local p_el = r(mean)
    sum omega if names =="U x English Home"
    local p_eng_home = r(mean)
    sum omega if names =="U x Lag Math"
    local p_lag_math = r(mean)
    sum omega if names =="U x Lag ELA"
    local p_lag_ela = r(mean)
    sum omega if names =="U x Income"
    local p_median_income = r(mean)
    sum omega if names =="U x In Magnet"
    local p_in_magnet = r(mean)
    sum omega if names =="U x Peer Q"
    local p_peer_q = r(mean)
    
    keep if regexm(names, "Mean Utility \d+")
    
    rename v1 sid 
    tempfile meanUtilities 
    save `meanUtilities', replace

    * Flag students in the lottery sample
    use "$DTAFINAL/lotteries.dta", clear
    tostring studentpseudoid , gen(sid) force usedisplay
    drop studentpseudoid 
    rename sid studentpseudoid
    keep if stu_grade ==6
    replace endyear = endyear -1 
    keep if inrange(endyear, 2004, 2008)
    keep studentpseudoid endyear
    isid studentpseudoid endyear
    tempfile stuids 
    save `stuids'

    * Read in posterior mean estimates 
    import delimited using "$ESTIMATES/in_mag_ind_`eta'_posteriors_`homog'_K`types'`app_suffix'.csv", clear stringcols(3)
    drop v1 
    rename (v4 v5) (mu_theta var_theta)
    isid studentpseudoid endyear 
    tempfile posteriors 
    save `posteriors'


    ********************************* Read in data and merge inputs  *********************************

    use "$DTAFINAL/structural_data_2004_`last_yr'.dta", clear 
    rename current_in_mag in_magnet 
    rename current_peer_quality peer_q 

    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela peer_q"
    merge m:1 studentpseudoid endyear using `posteriors', keepusing(mu_theta var_theta) gen(mergePosteriors) keep(1 3)
    merge m:1 studentpseudoid using `stuids', keepusing(endyear) gen(mergeStuids) keep(1 3)
    gen inLottery = mergeStuids==3
    gen magnet_enroll = D_0==0

    rename magnet_offer offer_var
    tab endyear, gen(yrdum)
    tab localdistrictcode, gen(districtdum)
    rename yrdum1 outyr 
    rename districtdum1 outdist

    local vars "female black white hispanic asian poverty el eng_home median_income  born_usa gifted missing_ela missing_math lag_math lag_ela  districtdum2 districtdum3 districtdum4 districtdum5 districtdum6 districtdum7"

    foreach var of local vars{
        egen meantemp = mean(`var')
        gen orig_`var' = `var'
        replace `var' = `var' - meantemp
        drop meantemp
    }
    
    gen orig_in_magnet = in_magnet
    gen orig_peer_q = peer_q
    gen theta_magnet = mu_theta * magnet_enroll
    *gen orig_median_income = median_income
    *drop magnet_mu_theta 

    * Don't have TE heterogeneity wrt to district dummies 
    local vars "female black white hispanic asian poverty el eng_home median_income  born_usa gifted missing_ela missing_math lag_math lag_ela"
    foreach var of local vars{
        gen magnet_`var' = magnet_enroll*`var'
    }

    ********************************* Create preference indices and percentiles *********************************

    gen theta_index = (mu_theta) 
    egen theta_percentile = rank(theta_index), track
    replace theta_percentile = theta_percentile / _N
    rename magnet_enroll magnet 


    local pref_vars "black white hispanic poverty el eng_home lag_math lag_ela median_income in_magnet peer_q"
    local pref_index "(`p_female') * orig_female"
    foreach var of local pref_vars{
        gen u_`var' = `p_`var'' 
        local pref_index "`pref_index' + (u_`var') * orig_`var'"
    }
    gen pref_index = (1)*(`pref_index')
    egen pref_percentile = rank(pref_index), track
    replace pref_percentile = pref_percentile / _N

    gen overall_index = mu_theta + pref_index 
    egen overall_percentile = rank(overall_index), track
    replace overall_percentile = overall_percentile / _N


     ********************************* Estimate the outcome model *********************************

    if("`blockfe'"=="TRUE"){
        local blockstub "absorb(censusblockid, savefe)"
    }
    else{
        local blockstub "noabsorb"
    }
    rename D_0 neighborhood_school 


    *  
    reghdfe  avg_F`subject'   D_* `vars'   magnet_*  ///  
        mu_theta  theta_magnet  yrdum*  districtdum*, absorb(censusblockid)  vce(cluster schoolcostcentercode) residuals(yres)
    gen yhat = avg_F`subject' - yres

    reghdfe avg_F`subject'   D_* `vars'  yrdum*  districtdum*, absorb(censusblockid)  vce(cluster schoolcostcentercode) residuals(yresvam)
    gen yhat_vam = avg_F`subject' - yresvam
    *reghdfe  avg_F`subject'   magnet `vars'  magnet_*  ///
        mu_theta  theta_magnet  yrdum*  districtdum*, absorb(censusblockid, savefe)  vce(cluster schoolcostcentercode)
*    reghdfe  z_socio   magnet `vars'  magnet_*  ///
        mu_theta  theta_magnet  yrdum*  districtdum*, absorb(censusblockid, savefe)  vce(cluster schoolcostcentercode)
    *local mag_effects_avg = _b[magnet]
    reghdfe  avg_F`subject'  magnet `vars'   magnet_*  ///  
        mu_theta  theta_magnet  yrdum* , `blockstub'  vce(cluster schoolcostcentercode)
    local mag_effects_avg = _b[magnet]


    ********************************* Calculate model implied Y1 and Y0 *********************************
    if("`blockfe'"=="TRUE"){
        rename __hdfe1__ blockfe
        local blockstub "+ blockfe"
    }
    else{
        local blockstub ""
    }
    /*
    preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
    restore 
   
    preserve 
        *rename neighborhood_school D_0
        collapse (mean) D_* 
        
        egen rowsum = rowtotal(D_*) 
        gen i =1 
        reshape long D_, i(i)
        replace D_ = D_ / rowsum
        rename _j D 
        keep D D_ 
        rename D_ enrollment_share 
        tempfile shares 
        save `shares', replace
    restore 
    preserve
        use `mag_data', clear
        merge m:1 D using `shares', keepusing(enrollment_share) gen(mergeShares) keep(1 3)
        gen magnet_effect = estimate * enrollment_share
        egen avg_effect = total(magnet_effect)
        sum avg_effect
        local mag_effects = r(mean)
        egen avg_effect2 = mean(estimate)
    restore 
    */


    local mag_effects = `mag_effects_avg'
    ds D_*   
    local magnet_schools `r(varlist)'
    local schooleffectstub ""
    foreach school of local magnet_schools{
        local schooleffectstub = "`schooleffectstub' + _b[`school']*`school'"
    }
    local maineffectstub ""
    local matcheffectstub ""
    foreach var of local vars{
        local maineffectstub = "`maineffectstub' + _b[`var']*`var'"
        local matcheffectstub = "`matcheffectstub' + _b[magnet_`var']*`var'"
    }
    local yearstub "+ _b[yrdum2]*yrdum2 + _b[yrdum3]*yrdum3 + _b[yrdum4]*yrdum4 + _b[yrdum5]*yrdum5 + _b[yrdum6]*yrdum6 + _b[yrdum7]*yrdum7 + _b[yrdum8]*yrdum8 + _b[yrdum9]*yrdum9 + _b[yrdum10]*yrdum10"
    *local districtstub "+ _b[districtdum2]*districtdum2 + _b[districtdum3]*districtdum3 + _b[districtdum4]*districtdum4 + _b[districtdum5]*districtdum5 + _b[districtdum6]*districtdum6 + _b[districtdum7]*districtdum7"
	
	
    gen y0 = _b[_cons] + _b[mu_theta]*mu_theta  `maineffectstub' `yearstub' `districtstub' `blockstub'
    gen y1 = y0 + _b[magnet] + _b[theta_magnet]*mu_theta `matcheffectstub'

    gen delta = y1 - y0 
	gen matcheffect = 0 `matcheffectstub'
	gen y1_novam = y0 + _b[theta_magnet]*mu_theta `matcheffectstub'
    gen y = y0 + magnet*(y1-y0)
	
	tempfile y1y0


	preserve 

	keep studentpseudoid endyear y0 y1 delta matcheffect y1_novam y overall_percentile theta_percentile mu_theta magnet  yhat yhat_vam
	save $ESTIMATES/y1y0_`eta'_`homog'_K`types'_`subject'_blockfe`blockfe'`app_suffix'.dta, replace 
	restore 
  
    stop;
   /*
    ********************************* Visual MTEs and Y1 and Y0 curves *********************************
    * WRT preference index 
    twoway lpoly delta overall_percentile if overall_percentile>0.05 & overall_percentile<0.95, ytitle("{&Epsilon}[Y{sub:1} - Y{sub:0} ]") xtitle("Inclination to Treatment (Percentile)")  lc(black)  
    graph export "$figures/cf_regressions_delta_`eta'_`homog'_K`types'_`subject'_blockfe`blockfe'`app_suffix'.pdf", replace

    twoway lpoly delta overall_percentile if overall_percentile>0.05 & overall_percentile<0.95, ytitle("{&Epsilon}[Y{sub:1} - Y{sub:0} ]") xtitle("Inclination to Treatment (Percentile)")  lc(white)
    graph export "$figures/cf_regressions0_delta_`eta'_`homog'_K`types'_`subject'_blockfe`blockfe'`app_suffix'.pdf", replace
    
    twoway lpoly delta overall_percentile if overall_percentile>0.05 & overall_percentile<0.95, ytitle("{&Epsilon}[Y{sub:1} - Y{sub:0} ]") xtitle("Inclination to Treatment (Percentile)")  lc(black)   
    graph export "$figures/cf_regressions1_delta_`eta'_`homog'_K`types'_`subject'_blockfe`blockfe'`app_suffix'.pdf", replace

    twoway (lpoly y1 overall_percentile if overall_percentile>0.05 & overall_percentile<0.95, lcolor(black) ) (lpoly y0 overall_percentile if overall_percentile>0.05 & overall_percentile<0.95, lcolor(maroon) ) , xtitle("Preference Index (Percentile)")  legend(order(1 "Y{sub:0}" 2 "Y{sub:1}") pos(6) row(1)) 
    graph export "$figures/cf_regressions_y1_y0_`eta'_`homog'_K`types'_`subject'_blockfe`blockfe'`app_suffix'.pdf", replace

    * WRT theta index
    twoway lpoly delta theta_percentile if theta_percentile>0.05 & theta_percentile<0.95, ytitle("{&Epsilon}[Y{sub:1} - Y{sub:0} ]") xtitle("Unobservable Inclination to Treatment (Percentile)")  lc(black)  
    graph export "$figures/cf_regressions_delta_theta_`eta'_`homog'_K`types'_`subject'_blockfe`blockfe'`app_suffix'.pdf", replace
	*/
end


********************************************************************************
* cf_regressions 
********************************************************************************/
capture program drop cf_regressions_noncog
program define cf_regressions_noncog
syntax, [eta(string) homog(string) types(int 1)  blockfe(string) app_2013(int 0)]

  
    
    ********************************* Prepare some inputs for the analysis *********************************

    local last_yr = cond(`app_2013', 2013, 2008)
    local app_suffix = cond(`app_2013', "_2013", "")
    
    * Read in demand model estimates to define preference index below 
    import delimited "$ESTIMATES/in_mag_ind_`eta'_estimates_choice_`homog'_K`types'`app_suffix'.csv", clear
    local vars "female black white hispanic poverty el eng_home lag_math lag_ela median_income curr"
    sum omega if names =="U x Female"
    local p_female = r(mean)
    sum omega if names =="U x Black"
    local p_black = r(mean)
    sum omega if names =="U x White"
    local p_white = r(mean)
    sum omega if names =="U x Hispanic"
    local p_hispanic = r(mean)
    sum omega if names =="U x Poverty"
    local p_poverty = r(mean)
    sum omega if names =="U x LEP"
    local p_el = r(mean)
    sum omega if names =="U x English Home"
    local p_eng_home = r(mean)
    sum omega if names =="U x Lag Math"
    local p_lag_math = r(mean)
    sum omega if names =="U x Lag ELA"
    local p_lag_ela = r(mean)
    sum omega if names =="U x Income"
    local p_median_income = r(mean)
    sum omega if names =="U x In Magnet"
    local p_in_magnet = r(mean)
    sum omega if names =="U x Peer Q"
    local p_peer_q = r(mean)
    
    keep if regexm(names, "Mean Utility \d+")
    
    rename v1 sid 
    tempfile meanUtilities 
    save `meanUtilities', replace

    * Flag students in the lottery sample
    use "$DTAFINAL/lotteries.dta", clear
    tostring studentpseudoid , gen(sid) force usedisplay
    drop studentpseudoid 
    rename sid studentpseudoid
    keep if stu_grade ==6
    replace endyear = endyear -1 
    keep if inrange(endyear, 2004, 2008)
    keep studentpseudoid endyear
    isid studentpseudoid endyear
    tempfile stuids 
    save `stuids'

    * Read in posterior mean estimates 
    import delimited using "$ESTIMATES/in_mag_ind_`eta'_posteriors_`homog'_K`types'`app_suffix'.csv", clear stringcols(3)
    drop v1 
    rename (v4 v5) (mu_theta var_theta)
    isid studentpseudoid endyear 
    tempfile posteriors 
    save `posteriors'


    ********************************* Read in data and merge inputs  *********************************

    use "$DTAFINAL/structural_data_2004_`last_yr'.dta", clear 
    rename current_in_mag in_magnet 
    rename current_peer_quality peer_q 

    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela peer_q"
    merge m:1 studentpseudoid endyear using `posteriors', keepusing(mu_theta var_theta) gen(mergePosteriors) keep(1 3)
    merge m:1 studentpseudoid using `stuids', keepusing(endyear) gen(mergeStuids) keep(1 3)
    gen inLottery = mergeStuids==3
    gen magnet_enroll = D_0==0

    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela peer_q mu_theta"

    foreach var of local vars{
        egen meantemp = mean(`var')
        gen orig_`var' = `var'
        replace `var' = `var' - meantemp
        drop meantemp
    }
        foreach var of local vars{
        gen magnet_`var' = magnet_enroll*`var'
    }
    gen orig_in_magnet = in_magnet
    gen theta_magnet = mu_theta * magnet_enroll
    drop magnet_mu_theta 

    rename magnet_offer offer_var
    tab endyear, gen(yrdum)
    tab localdistrictcode, gen(districtdum)
    rename yrdum1 outyr 
    rename districtdum1 outdist

    ********************************* Create preference indices and percentiles *********************************

    gen theta_index = (mu_theta) 
    egen theta_percentile = rank(theta_index), track
    replace theta_percentile = theta_percentile / _N
    rename magnet_enroll magnet 


    local pref_vars "black white hispanic poverty el eng_home lag_math lag_ela median_income in_magnet peer_q"
    local pref_index "(`p_female') * orig_female"
    foreach var of local pref_vars{
        gen u_`var' = `p_`var'' 
        local pref_index "`pref_index' + (u_`var') * orig_`var'"
    }
    gen pref_index = (1)*(`pref_index')
    egen pref_percentile = rank(pref_index), track
    replace pref_percentile = pref_percentile / _N

    gen overall_index = mu_theta + pref_index 
    egen overall_percentile = rank(overall_index), track
    replace overall_percentile = overall_percentile / _N


     ********************************* Estimate the outcome model *********************************

    if("`blockfe'"=="TRUE"){
        local blockstub "absorb(censusblockid, savefe)"
    }
    else{
        local blockstub "noabsorb"
    }
    rename D_0 neighborhood_school 
    local vars "female black white hispanic asian poverty el eng_home lag_math lag_ela peer_q"

    reghdfe  z_socio   magnet `vars'  magnet_*  ///
        mu_theta  theta_magnet  yrdum*  districtdum*, absorb(censusblockid, savefe)  vce(cluster schoolcostcentercode)
    local mag_effects_avg = _b[magnet]


    ********************************* Calculate model implied Y1 and Y0 *********************************
    if("`blockfe'"=="TRUE"){
        rename __hdfe1__ blockfe
        local blockstub "+ blockfe"
    }
    else{
        local blockstub ""
    }
    /*
    preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
    restore 
   
    preserve 
        *rename neighborhood_school D_0
        collapse (mean) D_* 
        
        egen rowsum = rowtotal(D_*) 
        gen i =1 
        reshape long D_, i(i)
        replace D_ = D_ / rowsum
        rename _j D 
        keep D D_ 
        rename D_ enrollment_share 
        tempfile shares 
        save `shares', replace
    restore 
    preserve
        use `mag_data', clear
        merge m:1 D using `shares', keepusing(enrollment_share) gen(mergeShares) keep(1 3)
        gen magnet_effect = estimate * enrollment_share
        egen avg_effect = total(magnet_effect)
        sum avg_effect
        local mag_effects = r(mean)
        egen avg_effect2 = mean(estimate)
    restore 
    */

    local mag_effects = `mag_effects_avg'
    ds D_*   
    local magnet_schools `r(varlist)'
    local schooleffectstub ""
    foreach school of local magnet_schools{
        local schooleffectstub = "`schooleffectstub' + _b[`school']*`school'"
    }
    local maineffectstub ""
    local matcheffectstub ""
    foreach var of local vars{
        local maineffectstub = "`maineffectstub' + _b[`var']*`var'"
        local matcheffectstub = "`matcheffectstub' + _b[magnet_`var']*`var'"
    }
    local yearstub "+ _b[yrdum2]*yrdum2 + _b[yrdum3]*yrdum3 + _b[yrdum4]*yrdum4 + _b[yrdum5]*yrdum5 + _b[yrdum6]*yrdum6 + _b[yrdum7]*yrdum7 + _b[yrdum8]*yrdum8 + _b[yrdum9]*yrdum9 + _b[yrdum10]*yrdum10"
    local districtstub "+ _b[districtdum2]*districtdum2 + _b[districtdum3]*districtdum3 + _b[districtdum4]*districtdum4 + _b[districtdum5]*districtdum5 + _b[districtdum6]*districtdum6 + _b[districtdum7]*districtdum7"
	
	
    gen y0 = _b[_cons] + _b[mu_theta]*mu_theta  `maineffectstub' `yearstub' `districtstub' `blockstub'
    gen y1 = y0 + _b[magnet] + _b[theta_magnet]*mu_theta `matcheffectstub'

    gen delta = y1 - y0 
	gen matcheffect = 0 `matcheffectstub'
	gen y1_novam = y0 + _b[theta_magnet]*mu_theta `matcheffectstub'
	
	tempfile y1y0

    /*
	preserve 

	keep studentpseudoid endyear y0 y1 delta matcheffect y1_novam
	save $ESTIMATES/y1y0_`eta'_`homog'_K`types'_`subject'_blockfe`blockfe'`app_suffix'.dta, replace 
	restore 
    */


    ********************************* Visual MTEs and Y1 and Y0 curves *********************************
    * WRT preference index 
    twoway lpoly delta overall_percentile if overall_percentile>0.05 & overall_percentile<0.95, ytitle("{&Epsilon}[Y{sub:1} - Y{sub:0} ]") xtitle("Inclination to Treatment (Percentile)")  lc(black)  
    graph export "$figures/cf_regressions_noncog_delta_`eta'_`homog'_K`types'_blockfe`blockfe'`app_suffix'.pdf", replace

end



capture program drop cf_table 
program define cf_table 
    syntax, [eta(string) homog(string) types(int 1)  blockfe(string) app_2013(int 0)]

 
    local last_yr = cond(`app_2013', 2013, 2008)
    local app_suffix = cond(`app_2013', "_2013", "")
    local n_schools = cond(`app_2013', 53, 40)
    
    ********************************* Prepare some inputs for the analysis *********************************

    * Read in posterior mean estimates 
    import delimited using "$ESTIMATES/in_mag_ind_`eta'_posteriors_`homog'_K`types'`app_suffix'.csv", clear stringcols(3)
    drop v1 
    rename (v4 v5) (mu_theta var_theta)
    isid studentpseudoid endyear 
    tempfile posteriors 
    save `posteriors'


    ********************************* Read in data and merge inputs  *********************************

    use $DTAFINAL/structural_data_2004_`last_yr'.dta, clear 
        rename current_in_mag in_magnet 
    rename current_peer_quality peer_q 

    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela peer_q"
    merge m:1 studentpseudoid endyear using `posteriors', keepusing(mu_theta var_theta) gen(mergePosteriors) keep(1 3)
    gen magnet_enroll = D_0==0

    gen theta_magnet = mu_theta * magnet_enroll
    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela peer_q"

    foreach var of local vars{
        egen meantemp = mean(`var')
        gen orig_`var' = `var'
        replace `var' = `var' - meantemp
        drop meantemp
    }
        foreach var of local vars{
        gen magnet_`var' = magnet_enroll*`var'
    }

    rename magnet_offer offer_var
    tab endyear, gen(yrdum)
    tab localdistrictcode, gen(districtdum)
    rename yrdum1 outyr 
    rename districtdum1 outdist
    rename magnet_enroll magnet
    ********************************* Estimate the outcome model *********************************
    if("`blockfe'"=="TRUE"){
        local blockstub "absorb(censusblockid, savefe)"
    }
    else{
        local blockstub "noabsorb"
    }
    rename D_0 neighborhood_school 
    local vars "female black white hispanic asian poverty el eng_home  lag_math lag_ela peer_q"
    reghdfe  avg_Fmath   D_* `vars'   magnet_*  ///  
        mu_theta  theta_magnet  yrdum*  districtdum*, `blockstub'  vce(cluster schoolcostcentercode)


    test mu_theta = theta_magnet=0
    local p_math : di %9.3f round(r(p), 0.001)
    di `p_math'
    local obs_math : di %12.0fc e(N)
    test magnet_female = magnet_black = magnet_white = magnet_hispanic = magnet_asian = magnet_poverty = magnet_el = magnet_eng_home = magnet_lag_math = magnet_lag_ela = theta_magnet=0 
    local p_hetero_math : di %9.3f round(r(p), 0.001)
  
    preserve 
        parmest, norestore

        drop if _n<=`n_schools'
        drop if _n>=26
	

        keep parm estimate stderr 
        rename (estimate stderr) (est_math se_math)
        tempfile math_estimates
        save `math_estimates', replace
        export delimited "$ESTIMATES/outcome_model_betas_`eta'_`homog'_K`types'`app_suffix'_math.csv", replace
    restore    


    * Estimate enrollment-weighted average of magnet effects 
      preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
        export delimited "$ESTIMATES/outcome_model_ates_`eta'_`homog'_K`types'`app_suffix'_math.csv", replace
    restore 
    preserve 
        collapse (mean) D_*  
        egen rowsum = rowtotal(D_*) 
        gen i =1 
        reshape long D_, i(i)
        replace D_ = D_ / rowsum
        rename _j D 
        keep D D_ 
        rename D_ enrollment_share 
        tempfile shares 
        save `shares', replace
    restore 
    preserve
        use `mag_data', clear
        merge m:1 D using `shares', keepusing(enrollment_share) gen(mergeShares) keep(1 3)
        gen magnet_effect = estimate * enrollment_share
        egen avg_effect = total(magnet_effect)
        sum avg_effect
        local mag_effects_math : di %9.3f round(r(mean), 0.001)
        egen mean_effect = mean(estimate)
        gen var_component = (estimate - mean_effect)^2 - stderr^2
        replace var_component = var_component * enrollment_share
        egen var = total(var_component)
        gen sd = sqrt(var)
        sum sd
        local sd_math : di %9.3f round(r(mean), 0.001)
    restore 

    * Save the parameter estimates and fixed effects; they are inputs for the structural estimation 
    reghdfe  avg_Fmath   D_* `vars'   magnet_*  ///  
        mu_theta  theta_magnet  , absorb(endyear censusblockid localdistrictcode, savefe)  vce(cluster schoolcostcentercode)
    rename (__hdfe1__ __hdfe2__ __hdfe3__) (yearfe blockfe districtfe)
	replace yearfe = _b[_cons] + yearfe + districtfe 
    preserve
    	keep studentpseudoid endyear yearfe blockfe 
        export delimited "$ESTIMATES/outcome_model_fes_`eta'_`homog'_K`types'`app_suffix'_math.csv", replace
    restore 
    drop yearfe blockfe 
    

    reghdfe  avg_Fela   D_* `vars'   magnet_*  ///  
        mu_theta  theta_magnet  yrdum*  districtdum*, `blockstub'  vce(cluster schoolcostcentercode)
    test mu_theta = theta_magnet=0
    local p_ela : di %9.3f round(r(p), 0.001)
    local obs_ela : di %12.0fc e(N)
    test magnet_female = magnet_black = magnet_white = magnet_hispanic = magnet_asian = magnet_poverty = magnet_el = magnet_eng_home = magnet_lag_math = magnet_lag_ela = theta_magnet=0 
    local p_hetero_ela : di %9.3f round(r(p), 0.001)

    preserve 
        parmest, norestore
        drop if _n<=`n_schools'
        drop if _n>=26
        keep parm estimate stderr 
        rename (estimate stderr) (est_ela se_ela)
        tempfile ela_estimates
        save `ela_estimates', replace
        export delimited "$ESTIMATES/outcome_model_betas_`eta'_`homog'_K`types'`app_suffix'_ela.csv", replace
    restore

      preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
        export delimited "$ESTIMATES/outcome_model_ates_`eta'_`homog'_K`types'`app_suffix'_ela.csv", replace
    restore 
    preserve
        use `mag_data', clear
        merge m:1 D using `shares', keepusing(enrollment_share) gen(mergeShares) keep(1 3)
        gen magnet_effect = estimate * enrollment_share
        egen avg_effect = total(magnet_effect)
        sum avg_effect
        local mag_effects_ela : di %9.3f round(r(mean), 0.001)
        egen mean_effect = mean(estimate)
        gen var_component = (estimate - mean_effect)^2 - stderr^2
        replace var_component = var_component * enrollment_share
        egen var = total(var_component)
        gen sd = sqrt(var)
        sum sd
        local sd_ela : di %9.3f round(r(mean), 0.001)
    restore 

    * Save the parameter estimates and fixed effects; they are inputs for the structural estimation 
    reghdfe  avg_Fela   D_* `vars'   magnet_*  ///  
        mu_theta  theta_magnet  districtdum*, absorb(endyear censusblockid, savefe)  vce(cluster schoolcostcentercode)
    rename (__hdfe1__ __hdfe2__) (yearfe blockfe)
	replace yearfe = _b[_cons] + yearfe 
    preserve
    	keep studentpseudoid endyear yearfe blockfe 
        export delimited "$ESTIMATES/outcome_model_fes_`eta'_`homog'_K`types'`app_suffix'_ela.csv", replace
    restore 

    use `math_estimates', clear
    merge 1:1 parm using `ela_estimates', gen(mergeEstimates) keep(1 3)

	texdoc init "$tables/`eta'_`homog'_K`types'`app_suffix'_cf_table.tex" , replace force 
	tex \begin{tabular}{lcccc} \toprule \hline
    tex & \multicolumn{2}{c}{Math} & \multicolumn{2}{c}{ELA} \\ \cline{2-5}
	tex  & Neighborhood School & Choice School & Neighborhood School & Choice School \\ 
    tex & (1) & (2) & (3) & (4) \\\hline \hline 
    tex Main Effects & & `mag_effects_math' & & `mag_effects_ela' \\
    tex &  & [`sd_math']  &  & [`sd_ela']  \\
    foreach var of local vars{
        if "`var'"=="female" local rowname "Female"
        else if "`var'"=="black" local rowname "Black"
        else if "`var'"=="white" local rowname "White"
        else if "`var'"=="hispanic" local rowname "Hispanic"
        else if "`var'"=="asian" local rowname "Asian"
        else if "`var'"=="poverty" local rowname "Poverty"
        else if "`var'"=="el" local rowname "LEP"
        else if "`var'"=="eng_home" local rowname "English Home"
        else if "`var'"=="lag_math" local rowname "Baseline Math"
        else if "`var'"=="lag_ela" local rowname "Baseline ELA"
        else if "`var'"=="peer_q" local rowname "Baseline Peer Quality"

        sum est_math if parm == "`var'"
        local est_math : di %9.3f round(r(mean), 0.001)
        sum se_math if parm == "`var'"
        local se_math : di %9.3f round(r(mean), 0.001)
        sum est_math if parm == "magnet_`var'"
        local est_math2 : di %9.3f round(r(mean), 0.001)
        sum se_math if parm == "magnet_`var'"
        local se_math2 : di %9.3f round(r(mean), 0.001)
        sum est_ela if parm == "`var'"
        local est_ela : di %9.3f round(r(mean), 0.001)
        sum se_ela if parm == "`var'"
        local se_ela : di %9.3f round(r(mean), 0.001)
        sum est_ela if parm == "magnet_`var'"
        local est_ela2 : di %9.3f round(r(mean), 0.001)
        sum se_ela if parm == "magnet_`var'"
        local se_ela2 : di %9.3f round(r(mean), 0.001)
        tex `rowname' & `est_math' & `est_math2' & `est_ela' & `est_ela2' \\ 
        tex & (`se_math') & (`se_math2') & (`se_ela') & (`se_ela2') \\ 
    }
    sum est_math if parm == "mu_theta" 
    local est_math : di %9.3f round(r(mean), 0.001)
    sum se_math if parm == "mu_theta"
    local se_math : di %9.3f round(r(mean), 0.001)
    sum est_math if parm == "theta_magnet"
    local est_math2 : di %9.3f round(r(mean), 0.001)
    sum se_math if parm == "theta_magnet"
    local se_math2 : di %9.3f round(r(mean), 0.001)
    sum est_ela if parm == "mu_theta"
    local est_ela : di %9.3f round(r(mean), 0.001)
    sum se_ela if parm == "mu_theta"
    local se_ela : di %9.3f round(r(mean), 0.001)
    sum est_ela if parm == "theta_magnet"
    local est_ela2 : di %9.3f round(r(mean), 0.001)
    sum se_ela if parm == "theta_magnet"
    local se_ela2 : di %9.3f round(r(mean), 0.001)
    tex Choice School Preference $\theta_i$ & `est_math' & `est_math2' & `est_ela' & `est_ela2' \\
    tex & (`se_math') & (`se_math2') & (`se_ela') & (`se_ela2') \\ 
	tex \addlinespace \hline \addlinespace
    tex Neighborhood Effects & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} \\
    tex Year Effects & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} \\
    tex District Effects & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} \\
	tex $H_0:$ No selection on unobservables  (p-values)& \multicolumn{2}{c}{`p_math'} & \multicolumn{2}{c}{`p_ela'} \\ 
    tex $H_0:$ No treatment effect heterogeneity (p-values) & \multicolumn{2}{c}{`p_hetero_math'} & \multicolumn{2}{c}{`p_hetero_ela'} \\
    tex Observations & \multicolumn{2}{c}{`obs_math'} & \multicolumn{2}{c}{`obs_ela'} \\ \addlinespace \hline\bottomrule
	tex  \end{tabular}
	texdoc close 
end

capture program drop cf_table2
program define cf_table2
    syntax, [eta(string) homog(string) types(int 1)  blockfe(string) app_2013(int 0)]

    local last_yr = cond(`app_2013', 2013, 2008)
    local app_suffix = cond(`app_2013', "_2013", "")
    local n_schools = cond(`app_2013', 53, 40)
    
    ********************************* Prepare some inputs for the analysis *********************************
    import delimited using "$ESTIMATES/in_mag_ind_`eta'_posteriors_`homog'_K`types'`app_suffix'.csv", clear stringcols(3)
    drop v1 
    rename (v4 v5) (mu_theta var_theta)
    isid studentpseudoid endyear 
    tempfile posteriors 
    save `posteriors'

    ********************************* Read in data and merge inputs  *********************************
    use $DTAFINAL/structural_data_2004_`last_yr'.dta, clear 
    rename current_in_mag in_magnet 
    rename current_peer_quality peer_q 

    local vars "female black white hispanic asian poverty el eng_home median_income lag_math lag_ela peer_q"
    merge m:1 studentpseudoid endyear using `posteriors', keepusing(mu_theta var_theta) gen(mergePosteriors) keep(1 3)
    gen magnet_enroll = D_0==0
    gen theta_magnet = mu_theta * magnet_enroll

    foreach var of local vars{
        egen meantemp = mean(`var')
        gen orig_`var' = `var'
        replace `var' = `var' - meantemp
        drop meantemp
    }
    foreach var of local vars{
        gen magnet_`var' = magnet_enroll*`var'
    }

    rename magnet_offer offer_var
    tab endyear, gen(yrdum)
    tab localdistrictcode, gen(districtdum)
    rename yrdum1 outyr 
    rename districtdum1 outdist
    rename magnet_enroll magnet

    ********************************* Estimate the outcome model *********************************
    if("`blockfe'"=="TRUE"){
        local blockstub "absorb(censusblockid, savefe)"
    }
    else{
        local blockstub "noabsorb"
    }
    rename D_0 neighborhood_school 
    local vars "female black white hispanic asian poverty el eng_home  lag_math lag_ela peer_q"

    * -------- MATH --------
    reghdfe  avg_Fmath   D_* `vars'   magnet_*  ///
        mu_theta  theta_magnet  yrdum*  districtdum*, `blockstub'  vce(cluster schoolcostcentercode)
    test mu_theta = theta_magnet=0
    local p_math : di %9.3f round(r(p), 0.001)
    local obs_math : di %12.0fc e(N)
    test magnet_female = magnet_black = magnet_white = magnet_hispanic = magnet_asian = magnet_poverty = magnet_el = magnet_eng_home = magnet_lag_math = magnet_lag_ela = theta_magnet=0 
    local p_hetero_math : di %9.3f round(r(p), 0.001)
  
    preserve 
        parmest, norestore
        drop if _n<=`n_schools'
        drop if _n>=26
        keep parm estimate stderr 
        rename (estimate stderr) (est_math se_math)
        tempfile math_estimates
        save `math_estimates', replace
        export delimited "$ESTIMATES/outcome_model_betas_`eta'_`homog'_K`types'`app_suffix'_math.csv", replace
    restore    

    * MATH ATEs
    preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
        export delimited "$ESTIMATES/outcome_model_ates_`eta'_`homog'_K`types'`app_suffix'_math.csv", replace
    restore 
    preserve 
        collapse (mean) D_*  
        egen rowsum = rowtotal(D_*) 
        gen i =1 
        reshape long D_, i(i)
        replace D_ = D_ / rowsum
        rename _j D 
        keep D D_ 
        rename D_ enrollment_share 
        tempfile shares 
        save `shares', replace
    restore 
    preserve
        use `mag_data', clear
        merge m:1 D using `shares', keepusing(enrollment_share) gen(mergeShares) keep(1 3)
        gen magnet_effect = estimate * enrollment_share
        egen avg_effect = total(magnet_effect)
        sum avg_effect
        local mag_effects_math : di %9.3f round(r(mean), 0.001)
        egen mean_effect = mean(estimate)
        gen var_component = (estimate - mean_effect)^2 - stderr^2
        replace var_component = var_component * enrollment_share
        egen var = total(var_component)
        gen sd = sqrt(var)
        sum sd
        local sd_math : di %9.3f round(r(mean), 0.001)
    restore 

    * Save FEs for structural inputs (math)
    reghdfe  avg_Fmath   D_* `vars'   magnet_*  ///
        mu_theta  theta_magnet  , absorb(endyear censusblockid localdistrictcode, savefe)  vce(cluster schoolcostcentercode)
    rename (__hdfe1__ __hdfe2__ __hdfe3__) (yearfe blockfe districtfe)
    replace yearfe = _b[_cons] + yearfe + districtfe 
    preserve
        keep studentpseudoid endyear yearfe blockfe 
        export delimited "$ESTIMATES/outcome_model_fes_`eta'_`homog'_K`types'`app_suffix'_math.csv", replace
    restore 
    drop yearfe blockfe 

    * -------- ELA --------
    reghdfe  avg_Fela   D_* `vars'   magnet_*  ///
        mu_theta  theta_magnet  yrdum*  districtdum*, `blockstub'  vce(cluster schoolcostcentercode)
    test mu_theta = theta_magnet=0
    local p_ela : di %9.3f round(r(p), 0.001)
    local obs_ela : di %12.0fc e(N)
    test magnet_female = magnet_black = magnet_white = magnet_hispanic = magnet_asian = magnet_poverty = magnet_el = magnet_eng_home = magnet_lag_math = magnet_lag_ela = theta_magnet=0 
    local p_hetero_ela : di %9.3f round(r(p), 0.001)

    preserve 
        parmest, norestore
        drop if _n<=`n_schools'
        drop if _n>=26
        keep parm estimate stderr 
        rename (estimate stderr) (est_ela se_ela)
        tempfile ela_estimates
        save `ela_estimates', replace
        export delimited "$ESTIMATES/outcome_model_betas_`eta'_`homog'_K`types'`app_suffix'_ela.csv", replace
    restore

    * ELA ATEs
    preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
        export delimited "$ESTIMATES/outcome_model_ates_`eta'_`homog'_K`types'`app_suffix'_ela.csv", replace
    restore 
    preserve
        use `mag_data', clear
        merge m:1 D using `shares', keepusing(enrollment_share) gen(mergeShares) keep(1 3)
        gen magnet_effect = estimate * enrollment_share
        egen avg_effect = total(magnet_effect)
        sum avg_effect
        local mag_effects_ela : di %9.3f round(r(mean), 0.001)
        egen mean_effect = mean(estimate)
        gen var_component = (estimate - mean_effect)^2 - stderr^2
        replace var_component = var_component * enrollment_share
        egen var = total(var_component)
        gen sd = sqrt(var)
        sum sd
        local sd_ela : di %9.3f round(r(mean), 0.001)
    restore 

    * Save FEs for structural inputs (ela)
    reghdfe  avg_Fela   D_* `vars'   magnet_*  ///
        mu_theta  theta_magnet  districtdum*, absorb(endyear censusblockid, savefe)  vce(cluster schoolcostcentercode)
    rename (__hdfe1__ __hdfe2__) (yearfe blockfe)
    replace yearfe = _b[_cons] + yearfe 
    preserve
        keep studentpseudoid endyear yearfe blockfe 
        export delimited "$ESTIMATES/outcome_model_fes_`eta'_`homog'_K`types'`app_suffix'_ela.csv", replace
    restore 

    * -------- NON-COGNITIVE (z_socio) --------
    reghdfe  z_socio   D_* `vars'  magnet_*  ///
        mu_theta  theta_magnet  yrdum*  districtdum*, `blockstub'  vce(cluster schoolcostcentercode)
    test mu_theta = theta_magnet=0
    local p_socio : di %9.3f round(r(p), 0.001)
    local obs_socio : di %12.0fc e(N)
    test magnet_female = magnet_black = magnet_white = magnet_hispanic = magnet_asian = magnet_poverty = magnet_el = magnet_eng_home = magnet_lag_math = magnet_lag_ela = theta_magnet=0 
    local p_hetero_socio : di %9.3f round(r(p), 0.001)

    preserve 
        parmest, norestore
        drop if _n<=`n_schools'
        drop if _n>=26
        keep parm estimate stderr 
        rename (estimate stderr) (est_socio se_socio)
        tempfile socio_estimates
        save `socio_estimates', replace
        export delimited "$ESTIMATES/outcome_model_betas_`eta'_`homog'_K`types'`app_suffix'_socio.csv", replace
    restore

    * SOCIO ATEs
    preserve 
        parmest, norestore 
        keep if regexm(parm, "D_")
        sum estimate 
        local mag_effects = r(mean)
        tempfile mag_data 
        split parm, parse("_")
        rename parm2 D 
        destring D, replace 
        keep D estimate stderr 
        save `mag_data', replace
        export delimited "$ESTIMATES/outcome_model_ates_`eta'_`homog'_K`types'`app_suffix'_socio.csv", replace
    restore 
    preserve
        use `mag_data', clear
        merge m:1 D using `shares', keepusing(enrollment_share) gen(mergeShares) keep(1 3)
        gen magnet_effect = estimate * enrollment_share
        egen avg_effect = total(magnet_effect)
        sum avg_effect
        local mag_effects_socio : di %9.3f round(r(mean), 0.001)
        egen mean_effect = mean(estimate)
        gen var_component = (estimate - mean_effect)^2 - stderr^2
        replace var_component = var_component * enrollment_share
        egen var = total(var_component)
        gen sd = sqrt(var)
        sum sd
        local sd_socio : di %9.3f round(r(mean), 0.001)
    restore 

    * -------- Combine estimates for table --------
    use `math_estimates', clear
    merge 1:1 parm using `ela_estimates', gen(mergeEstimates1) keep(1 3)
    merge 1:1 parm using `socio_estimates', gen(mergeEstimates2) keep(1 3)

    * -------- LaTeX table (now 6 columns) --------
    texdoc init "$tables/`eta'_`homog'_K`types'`app_suffix'_cf_table.tex" , replace force 
    tex \begin{tabular}{lcccccc} \toprule \hline
    tex & \multicolumn{2}{c}{Math} & \multicolumn{2}{c}{ELA} & \multicolumn{2}{c}{Non-Cognitive Index} \\ \cline{2-7}
    tex  & Neighborhood School & Choice School & Neighborhood School & Choice School & Neighborhood School & Choice School \\ 
    tex & (1) & (2) & (3) & (4) & (5) & (6) \\\hline \hline 
    tex Main Effects & & `mag_effects_math' & & `mag_effects_ela' & & `mag_effects_socio' \\
    tex &  & [`sd_math']  &  & [`sd_ela']  &  & [`sd_socio'] \\
    foreach var of local vars{
        if "`var'"=="female" local rowname "Female"
        else if "`var'"=="black" local rowname "Black"
        else if "`var'"=="white" local rowname "White"
        else if "`var'"=="hispanic" local rowname "Hispanic"
        else if "`var'"=="asian" local rowname "Asian"
        else if "`var'"=="poverty" local rowname "Poverty"
        else if "`var'"=="el" local rowname "LEP"
        else if "`var'"=="eng_home" local rowname "English Home"
        else if "`var'"=="lag_math" local rowname "Baseline Math"
        else if "`var'"=="lag_ela" local rowname "Baseline ELA"
        else if "`var'"=="peer_q" local rowname "Baseline Peer Quality"

        * Math
        sum est_math if parm == "`var'"
        local est_math : di %9.3f round(r(mean), 0.001)
        sum se_math if parm == "`var'"
        local se_math : di %9.3f round(r(mean), 0.001)
        sum est_math if parm == "magnet_`var'"
        local est_math2 : di %9.3f round(r(mean), 0.001)
        sum se_math if parm == "magnet_`var'"
        local se_math2 : di %9.3f round(r(mean), 0.001)

        * ELA
        sum est_ela if parm == "`var'"
        local est_ela : di %9.3f round(r(mean), 0.001)
        sum se_ela if parm == "`var'"
        local se_ela : di %9.3f round(r(mean), 0.001)
        sum est_ela if parm == "magnet_`var'"
        local est_ela2 : di %9.3f round(r(mean), 0.001)
        sum se_ela if parm == "magnet_`var'"
        local se_ela2 : di %9.3f round(r(mean), 0.001)

        * Socio
        sum est_socio if parm == "`var'"
        local est_socio : di %9.3f round(r(mean), 0.001)
        sum se_socio if parm == "`var'"
        local se_socio : di %9.3f round(r(mean), 0.001)
        sum est_socio if parm == "magnet_`var'"
        local est_socio2 : di %9.3f round(r(mean), 0.001)
        sum se_socio if parm == "magnet_`var'"
        local se_socio2 : di %9.3f round(r(mean), 0.001)

        tex `rowname' & `est_math' & `est_math2' & `est_ela' & `est_ela2' & `est_socio' & `est_socio2' \\ 
        tex & (`se_math') & (`se_math2') & (`se_ela') & (`se_ela2') & (`se_socio') & (`se_socio2') \\ 
    }

    * Theta rows
    sum est_math if parm == "mu_theta" 
    local est_math : di %9.3f round(r(mean), 0.001)
    sum se_math if parm == "mu_theta"
    local se_math : di %9.3f round(r(mean), 0.001)
    sum est_math if parm == "theta_magnet"
    local est_math2 : di %9.3f round(r(mean), 0.001)
    sum se_math if parm == "theta_magnet"
    local se_math2 : di %9.3f round(r(mean), 0.001)

    sum est_ela if parm == "mu_theta"
    local est_ela : di %9.3f round(r(mean), 0.001)
    sum se_ela if parm == "mu_theta"
    local se_ela : di %9.3f round(r(mean), 0.001)
    sum est_ela if parm == "theta_magnet"
    local est_ela2 : di %9.3f round(r(mean), 0.001)
    sum se_ela if parm == "theta_magnet"
    local se_ela2 : di %9.3f round(r(mean), 0.001)

    sum est_socio if parm == "mu_theta"
    local est_socio : di %9.3f round(r(mean), 0.001)
    sum se_socio if parm == "mu_theta"
    local se_socio : di %9.3f round(r(mean), 0.001)
    sum est_socio if parm == "theta_magnet"
    local est_socio2 : di %9.3f round(r(mean), 0.001)
    sum se_socio if parm == "theta_magnet"
    local se_socio2 : di %9.3f round(r(mean), 0.001)

    tex Choice School Preference $\theta_i$ & `est_math' & `est_math2' & `est_ela' & `est_ela2' & `est_socio' & `est_socio2' \\
    tex & (`se_math') & (`se_math2') & (`se_ela') & (`se_ela2') & (`se_socio') & (`se_socio2') \\ 
    tex \addlinespace \hline \addlinespace
    tex Neighborhood Effects & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} \\
    tex Year Effects & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} \\
    tex District Effects & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} & \multicolumn{2}{c}{\checkmark} \\
    tex $H_0:$ No selection on unobservables  (p-values)& \multicolumn{2}{c}{`p_math'} & \multicolumn{2}{c}{`p_ela'} & \multicolumn{2}{c}{`p_socio'} \\ 
    tex $H_0:$ No treatment effect heterogeneity (p-values) & \multicolumn{2}{c}{`p_hetero_math'} & \multicolumn{2}{c}{`p_hetero_ela'} & \multicolumn{2}{c}{`p_hetero_socio'} \\
    tex Observations & \multicolumn{2}{c}{`obs_math'} & \multicolumn{2}{c}{`obs_ela'} & \multicolumn{2}{c}{`obs_socio'} \\ \addlinespace \hline\bottomrule
    tex  \end{tabular}
    texdoc close 
end


capture program drop maxAllocation
program define maxAllocation 
	tempfile seats
	tempfile va 

	

	* Read in vam estimates 
	import delimited "/Volumes/lausd/decentralized_choice/estimates/outcome_model_ates_eta_out_hetero_K3_2013_math.csv", clear 
	tempfile va 
	save `va'

	use "/Volumes/lausd/decentralized_choice/data/prog_seats_2004_2013.dta"
	reshape long seats_, i(endyear stu_grade phbao ) 
	rename _j d
	save `seats', replace
	
	merge m:1 d using `va', gen(mergeVAM) keep(1 3) 
	
	reshape wide seats_ , i(endyear d stu_grade estimate stderr ) j(phbao )
	
	gsort endyear -estimate  
	
	bys endyear: gen vamrank = _n
	
	tempfile maxalloc 
	save `maxalloc', replace 
	
	use $DTAFINAL/structural_data_2004_2013.dta, clear 
	* merge in match effects 
	merge 1:1 studentpseudoid endyear using $ESTIMATES/y1y0_eta_out_hetero_K3_math_blockfeTRUE_2013.dta, gen(mergeMatchEffects) keep(1 3)
	gen phbao = 1- white
	gsort endyear phbao -matcheffect 

	bys endyear phbao: gen stu_rank = _n
	
	gsort endyear -matcheffect
	bys endyear: gen stu_rank2 = _n 
	
	gen y1_max = .
	
	gen y1_max2 = .
	

	* Loop through each one of the 53 schools in each year 
	forvalues yr = 2004/2013{
		local phbao0_ul = 0 
		local phbao0_ll = 0
		local phbao1_ul = 0
		local phbao1_ll = 0 
		
		local all_ll = 0
		forvalues j=1/53{
			preserve 
				use `maxalloc', clear 
				* vam estimate to assign 
				sum estimate if vamrank==`j' & endyear == `yr'
				local vam = r(mean)
				* seats to allocate to phbao 
				sum seats_1 if vamrank==`j' & endyear == `yr'
				local s1 = r(mean)
				* seats to allocate to non-phbao 
				sum seats_0 if vamrank==`j' & endyear == `yr'
				local s0 = r(mean) 
			restore 
			if `s1' >0{
				replace y1_max = y1_novam + `vam' if endyear == `yr' & phbao==1 & stu_rank > `phbao1_ll' & stu_rank <= `phbao1_ll' + `s1'
				local phbao1_ll = `phbao1_ll' + `s1'  
				}
			if `s0' >0 {
				replace y1_max = y1_novam + `vam' if  endyear == `yr' & phbao==0 & stu_rank > `phbao0_ll' & stu_rank <= `phbao0_ll' + `s0'
				local phbao0_ll = `phbao0_ll' + `s0'  	
				}
				* Create a third that ignores the race-specific seats 
				if (`s1' > 0 | `s0' > 0)  {
					replace y1_max2 = y1_novam + `vam' if  endyear == `yr' & stu_rank2 > `all_ll' & stu_rank2 <= `all_ll' + `s0' + `s1'
				local all_ll = `all_ll' + `s0' + `s1'
				}
			
				
				
		}
	}
	
	gen delta_max2 = y1_max2 - y0 
	replace delta_max2 = 0 if missing(delta_max2)
	gen delta_max = y1_max - y0 
	replace delta_max = 0 if missing(delta_max )
	rename poverty pov 
	gen nonpov = 1- pov 
		gen treat1 = delta_max !=0
	gen treat2 = delta_max2 !=0 
	local groups "black white hisp asian pov nonpov"
	local subject "math"

	foreach g of local groups{
        gen y_`g'_`subject'1 = y1_max if treat1==1 & `g'==1
		gen y_`g'_`subject'2 = y1_max2 if treat2==1 & `g'==1
		replace y_`g'_`subject'1 = y0 if treat1==0 & `g'==1
		replace y_`g'_`subject'2 = y0 if treat2==0 & `g'==1
	}
	gen yobs_math_1 = y1_max if treat1==1
	replace yobs_math_1 = y0 if treat1==0 
	gen yobs_math_2 = y1_max2 if treat2==1
	replace yobs_math_2 = y0 if treat2==0
	
	gen tot_math_1 = delta_max if treat1==1
	gen tot_math_2 = delta_max2 if treat2==1 
	
	collapse (mean) tot_math_* yobs_math_* y_*_math*  delta_max delta_max2 ///
		(sd) sd_yobs_math_1 = yobs_math_1 sd_yobs_math2 = yobs_math_2 
			
	rename delta_max delta_max1
	rename sd_yobs_math2 sd_yobs_math_2
	gen i = 1	
	reshape long tot_math_ yobs_math_ y_black_math y_white_math y_hisp_math y_asian_math y_pov_math y_nonpov_math delta_max sd_yobs_math_, i(i) 
	tempfile math 
	save `math'
	
	

		

	* Read in vam estimates 
	import delimited "/Volumes/lausd/decentralized_choice/estimates/outcome_model_ates_eta_out_hetero_K3_2013_ela.csv", clear 
	tempfile va 
	save `va'

	use "/Volumes/lausd/decentralized_choice/data/prog_seats_2004_2013.dta"
	reshape long seats_, i(endyear stu_grade phbao ) 
	rename _j d
	save `seats', replace
	
	merge m:1 d using `va', gen(mergeVAM) keep(1 3) 
	
	reshape wide seats_ , i(endyear d stu_grade estimate stderr ) j(phbao )
	
	gsort endyear -estimate  
	
	bys endyear: gen vamrank = _n
	
	tempfile maxalloc 
	save `maxalloc', replace 
	
	use $DTAFINAL/structural_data_2004_2013.dta, clear 
	* merge in match effects 
	merge 1:1 studentpseudoid endyear using $ESTIMATES/y1y0_eta_out_hetero_K3_ela_blockfeTRUE_2013.dta, gen(mergeMatchEffects) keep(1 3)
	gen phbao = 1- white
	gsort endyear phbao -matcheffect 

	bys endyear phbao: gen stu_rank = _n
	
	gsort endyear -matcheffect
	bys endyear: gen stu_rank2 = _n 
	
	gen y1_max = .
	
	gen y1_max2 = .
	

	* Loop through each one of the 53 schools in each year 
	forvalues yr = 2004/2013{
		local phbao0_ul = 0 
		local phbao0_ll = 0
		local phbao1_ul = 0
		local phbao1_ll = 0 
		
		local all_ll = 0
		forvalues j=1/53{
			preserve 
				use `maxalloc', clear 
				* vam estimate to assign 
				sum estimate if vamrank==`j' & endyear == `yr'
				local vam = r(mean)
				* seats to allocate to phbao 
				sum seats_1 if vamrank==`j' & endyear == `yr'
				local s1 = r(mean)
				* seats to allocate to non-phbao 
				sum seats_0 if vamrank==`j' & endyear == `yr'
				local s0 = r(mean) 
			restore 
			if `s1' >0{
				replace y1_max = y1_novam + `vam' if endyear == `yr' & phbao==1 & stu_rank > `phbao1_ll' & stu_rank <= `phbao1_ll' + `s1'
				local phbao1_ll = `phbao1_ll' + `s1'  
				}
			if `s0' >0 {
				replace y1_max = y1_novam + `vam' if  endyear == `yr' & phbao==0 & stu_rank > `phbao0_ll' & stu_rank <= `phbao0_ll' + `s0'
				local phbao0_ll = `phbao0_ll' + `s0'  	
				}
				* Create a third that ignores the race-specific seats 
				if (`s1' > 0 | `s0' > 0)  {
					replace y1_max2 = y1_novam + `vam' if  endyear == `yr' & stu_rank2 > `all_ll' & stu_rank2 <= `all_ll' + `s0' + `s1'
				local all_ll = `all_ll' + `s0' + `s1'
				}
			
				
				
		}
	}
	
	gen delta_max2 = y1_max2 - y0 
	replace delta_max2 = 0 if missing(delta_max2)
	gen delta_max = y1_max - y0 
	replace delta_max = 0 if missing(delta_max )
	rename poverty pov 
	gen nonpov = 1- pov 
		gen treat1 = delta_max !=0
	gen treat2 = delta_max2 !=0 
	local groups "black white hisp asian pov nonpov"
	local subject "ela"

	foreach g of local groups{
        gen y_`g'_`subject'1 = y1_max if treat1==1 & `g'==1
		gen y_`g'_`subject'2 = y1_max2 if treat2==1 & `g'==1
		replace y_`g'_`subject'1 = y0 if treat1==0 & `g'==1
		replace y_`g'_`subject'2 = y0 if treat2==0 & `g'==1
	}
	gen yobs_ela_1 = y1_max if treat1==1
	replace yobs_ela_1 = y0 if treat1==0 
	gen yobs_ela_2 = y1_max2 if treat2==1
	replace yobs_ela_2 = y0 if treat2==0
	
	gen tot_ela_1 = delta_max if treat1==1
	gen tot_ela_2 = delta_max2 if treat2==1 
	
	gen attend_rate_black1 = delta_max !=0 if black==1
	gen attend_rate_hisp1 = delta_max !=0 if hispanic==1
	gen attend_rate_asian1 = delta_max !=0 if asian ==1 
	gen attend_rate_white1 = delta_max !=0 if white==1
	
	gen attend_rate_black2 = delta_max2 !=0 if black==1
	gen attend_rate_hisp2 = delta_max2 !=0 if hispanic==1
	gen attend_rate_asian2 = delta_max2 !=0 if asian ==1 
	gen attend_rate_white2 = delta_max2 !=0 if white==1
	
	collapse (mean) tot_ela_* yobs_ela_* y_*_ela*  delta_max delta_max2 ///
			attend_rate_* ///
		(sd) sd_yobs_ela_1 = yobs_ela_1 sd_yobs_ela2 = yobs_ela_2 


	rename delta_max delta_max1
	rename sd_yobs_ela2 sd_yobs_ela_2
	gen i = 1	
	reshape long tot_ela_ yobs_ela_ y_black_ela y_white_ela y_hisp_ela y_asian_ela y_pov_ela y_nonpov_ela delta_max sd_yobs_ela_ attend_rate_asian attend_rate_black attend_rate_hisp attend_rate_white, i(i) 
	rename delta_max delta_max_ela
	

	merge 1:1 i _j using `math'
	
	rename (tot_ela_ yobs_ela_ sd_yobs_ela_ tot_math_ yobs_math_ sd_yobs_math_) ///
		(tot_ela yobs_ela sd_yobs_ela tot_math yobs_math sd_yobs_math)
	rename delta_max delta_max_math
	stop;
	
	export delimited $ESTIMATES/counterfactual_max_policy.csv, replace 
end 
main 
