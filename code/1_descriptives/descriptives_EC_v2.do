/********************************************************************************
1_descriptives.do
- Author: Chris Campos and Antonia Vazquez
- Description: Conducts descriptive and motivating analysis  

- Date started: 07/10/2025

* Installs: unique
*******************************************************************************/
clear all
set trace on 
set tracedepth 2
*set maxvar 120000
*Packages needed
*ssc inst egenmore
*ssc install coefplot, replace
*ssc install seq
*ssc install unique
*net install grc1leg, from( http://www.stata.com/users/vwiggins/)

capture program drop main 
program define main 

	set_paths 

	* Appendix Table
	*sampleComparisons
	
	* Figure 1 
	*motivation
	* Figure 3 
	*application_distance_descriptive
	* Table 1 
	*whoApplies
	* Figure 4 
	distance_correlates
	* Figure 2 
	*distance_maps
end 



********************************************************************************
* set_paths 
* 
* Description: Establish paths 
********************************************************************************
capture program drop set_paths
program define set_paths

	/*
	if "`c(os)'"=="MacOSX" & "`c(username)'"=="cqcampos"{
		global ROOT "/Volumes/lausd/decentralized_choice"
		global BUILD "/Volumes/lausd/build/output/data"
    	global MAGNET "/Volumes/lausd/magnet"
		global RAWDATA "$BUILD/rawdata/Campos_Master_DUA" 
		global BUILDROOT "/Volumes/lausd/build/" 
	}
	else{
	*/
			global ROOT "/Volumes/lausd/decentralized_choice"
			global BUILD "/Volumes/lausd/build/output/data"
			global MAGNET "/Volumes/lausd/magnet"
			global RAWDATA "$BUILD/rawdata/Campos_Master_DUA" 
			global BUILDROOT "/Volumes/lausd/build/" 

	*}

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
	global figures "$ROOT/output/figures"
	global slidefigures "$ROOT/output/slidefigures"

	* Old directories, to be updated 
	global LOGS "$ROOT/logs"

	* Raw source data
	global UE "$BUILDROOT/rawdata/Campos_Master_DUA/Unified_Enrollment"
end

********************************************************************************
* motivation
*
* Description: Creates a time series plot of number of choice programs and 
* declining enrollment 
*
********************************************************************************
capture program drop motivation
program define motivation

	/* Code to create aggregates using NSC data
	use "/Users/cqcampos/Library/CloudStorage/Dropbox/chyn_campos_bruhn/Magnet/data/national level/ccd_directory_1986_2020.dta", clear
	replace state_leaid = "1964733" if regexm(state_leaid, "1964733")
	replace state_leaid = "1964733" if regexm(state_leaid, "19647")
	replace lea_name = proper(lea_name )
	tab charter 
	sum charter 
	tab charter if charter==1
	drop if charter==1
	collapse (sum) enrollment , by(leaid lea_name state_leaid year) 
	keep if leaid =="0622710"
	*/
    
	*create an empty dataset to store our results
	clear
	set obs 23
	gen endyear = 2000 + _n
	tempfile results
	save `results'
	
    * Open LAUSD dataset and keep the n of students per year 
    use  endyear using "$BUILD/lausd2001_2023_encoded.dta", clear
    gen stu = 1
    collapse (count) stu, by(endyear)
    rename stu n_students
    sort endyear 
    tempfile temp 
    save `temp'

    * Describe the lotteries dataset
    use "$DTAFINAL/magnet_lotteries_1perc_b.dta", clear
    keep if choice_ranking ==1
    contract mag_cd_ endyear
    collapse (count) mag_cd_ , by(endyear)
    tempfile temp2 
    save `temp2'
	
	*Magnet *********
	use "$ROOT/rawdata/magnet_rank.dta", clear 
	keep if choice_ranking ==1
    contract mag_cd_ endyear
	*AV: drop if missing or code is 0
	drop if mag_cd_==. |  mag_cd_==0
	*AV: kepp only the programs that appear more than 10 times
    drop if _freq<10
	collapse (count) mag_cd_ , by(endyear)  
    rename mag_cd_ n_mag
	tempfile mag 
    save `mag'
	
	*sas ********* Check how this look like in the build before adding them 
	use "$BUILD/sas_rank.dta", clear
	keep if !missing(studentpseudoid)
	gen endyear= school_yr+1
	drop if sas_ch_=="NA"
	contract sas_ch_ endyear
	destring sas_ch_, replace
	collapse (count) sas_ch_ , by(endyear)
	rename sas_ch_ n_sas
	
	tempfile sas 
    save `sas'
	
	*acs  ********* these are 10 and stable acrross time 
	use "$BUILD/acs_rank.dta", clear 
	keep if !missing(studentpseudoid)
	gen endyear= school_yr+1
	contract acs_ch_1 endyear
	destring acs_ch_1, replace
	collapse (count) acs_ch_1 , by(endyear)  
	rename acs_ch_1 n_acs

	tempfile acs 
    save `acs'
	
	*afc ********* 
	use "$BUILD/afc_app.dta", clear 
	keep if !missing(studentpseudoid)
	gen endyear= school_yr+1
	contract cost_center_code endyear
	*There is a DLE program here, I remove the letters 
	replace cost_center_code ="1511102" if cost_center_code=="1511102TWS"
	destring cost_center_code, replace
	collapse (count) cost_center_code , by(endyear)  
	rename cost_center_code n_afc

	tempfile afc 
    save `afc'
   
   *merge all results 
	use `results', clear
    merge 1:1 endyear using `temp', keep(1 3) nogen
    merge 1:1 endyear using `temp2', keep(1 3) nogen
	replace n_students = n_students/100000
	
	merge 1:1 endyear using `mag', keep(1 3) nogen
	merge 1:1 endyear using `sas', keep(1 3) nogen
	merge 1:1 endyear using `acs', keep(1 3) nogen
	merge 1:1 endyear using `afc', keep(1 3) nogen
	
	tempfile results_filled
	save `results_filled'
	

	
	*Now we try different ways to identify DLE programs
	*dle ********* From the applications data 
	use "$BUILD/dle_rank.dta" , clear
	keep if choice_ranking ==1
	contract choice endyear
	encode choice, gen(choice_code)
	collapse (count) choice_code , by(endyear)  
	rename choice_code n_dle
	tempfile dle 
    save `dle'
   
	*dle ******* From the LAUSD 2001 2023 data using the name of the program

	**Get magnet codes
	use "$ROOT/rawdata/magnet_rank.dta", clear
    keep if choice_ranking ==1
	 contract mag_cd_ endyear
	drop if mag_cd_==. | mag_cd_==0
	*There is an issue with the magnet codes, so we collapse at the costcentercode (mag_cd_) and drop programs that are applied to less than 10 times across all years
	drop if _freq <10
	drop _freq
	contract mag_cd_
	rename mag_cd_ magnet_code 
	tempfile tempcodes 
	save `tempcodes'
    
	use  endyear preferredlocationcode schoollocationcode schoolcostcentercode studentpseudoid gradecode schoollocationname using "$BUILD/lausd2001_2023_encoded.dta", clear
		
	*unique endyear studentpseudoid gradecode 
	*unique preferredlocationcode
	*unique schoollocationcode
	*unique schoolcostcentercode
	*unique preferredlocationcode schoollocationcode schoolcostcentercode schoollocationname
	*1676 unique combinations
	gen magnet_code = schoolcostcentercode
	merge m:1 magnet_code using `tempcodes'
	 drop if _merge==2
	 *I keep the schools that are not magnet 
	drop if  _merge==3
	collapse gradecode, by(preferredlocationcode schoollocationcode schoolcostcentercode schoollocationname )
	drop gradecode
   
   sort preferredlocationcode schoollocationcode schoolcostcentercode schoollocationname
   
   drop if missing(schoollocationname)
   gen last_two_digits_costcode = mod(schoolcostcentercode, 100)
   
   *Drop the codes that are likely regular schools identify by the last 2 digits of schoolcostcentercode
   drop if last_two_digits_costcode==1
 
   decode schoollocationname, gen(schoollocationname_string)
   
   *Identify potential DLE 
   gen potential_dle = regexm(schoollocationname_string, " EL ")
   replace potential_dle = 1 if regexm(schoollocationname_string, " DT")
   replace potential_dle = 1 if regexm(schoollocationname_string, " DOS")
   replace potential_dle = 1 if regexm(schoollocationname_string, " DLC ")
   replace potential_dle = 1 if regexm(schoollocationname_string, " DL")

   gen has_slc = regexm(schoollocationname_string, " SLC ")
   drop if has_slc==1 
   drop if potential_dle==0
   contract schoolcostcentercode
   keep schoolcostcentercode
   
   tempfile dle_codes
   save `dle_codes', replace
   
   **Now i use this to count the programs in LAUSD 2001 to 2023
   use  endyear preferredlocationcode schoollocationcode schoolcostcentercode studentpseudoid gradecode schoollocationname using "$BUILD/lausd2001_2023_encoded.dta", clear
   merge m:1 schoolcostcentercode using `dle_codes', gen(DLE_prog)
	keep if DLE_prog==3
   contract schoolcostcentercode endyear
	
	collapse (count) schoolcostcentercode , by(endyear)  
	rename schoolcostcentercode n_dle2
	tempfile dle2
    save `dle2', replace 
	
	
	*dle ******* From the data on DLE enrollment that Pablo shared 
	*This is the most conservative and the one we are using 
	 import excel using  "$BUILDROOT/rawdata/ZOC Obfuscation 20250214/DUAL_LANGUAGE_STUDENTS_DEIDENTIFIED.xlsx", first clear
	 drop FIRSTNAME LASTNAME MIDDLE_INIT HOME_ADDRESS HOME_CITY HOME_ZIP
	 rename *, lower
	 destring school_yr, replace
	 gen endyear=school_yr+1
	 *Delete magnets
	    gen has_mag = regexm(current_school_name, "MAG")
		drop if has_mag==1
		 gen has_mg = regexm(current_school_name, "MG")
		drop if has_mg==1
		
	contract current_school_code endyear
	encode current_school_code, gen(n_dle3)
	collapse (count) n_dle3 , by(endyear)  
	
	
	tempfile dle3
    save `dle3', replace 
		
	use `results_filled', clear 
	merge 1:1 endyear using `dle3', keep(1 3) nogen	
	
	* For 2017, inferred from this page: https://www.laschoolreport.com/dual-language-programs-are-so-popular-that-lausd-plans-to-double-the-number-of-schools-offering-them-by-next-year/
	replace n_dle3 = 13 if endyear==2005
	replace n_dle3 = 44 if endyear==2013
	replace n_dle3 = 57 if endyear==2015
	replace n_dle3 = 56 if endyear==2016
	replace n_dle3 = 101 if endyear==2017
	
	ipolate n_dle3 endyear, generate(n_dle_final)
	replace n_dle_final =round(n_dle_final)
	
	* Numbers taken from output produced by independent research
	gen n_charter = 10 if endyear==2001 
	replace n_charter = 10 if endyear==2002
	replace n_charter = 12 if endyear==2003
	replace n_charter = 14 if endyear==2004
	replace n_charter = 14 if endyear==2005
	replace n_charter = 14 if endyear==2006
	replace n_charter = 14 if endyear==2007
	replace n_charter = 14 if endyear==2008
	replace n_charter = 16 if endyear==2009
	replace n_charter = 17 if endyear==2010
	replace n_charter = 18 if endyear==2011
	replace n_charter = 24 if endyear==2012
	replace n_charter = 43 if endyear==2013
	replace n_charter = 47 if endyear==2014
	replace n_charter = 48 if endyear==2015
	replace n_charter = 48 if endyear==2016
	replace n_charter = 49 if endyear==2017
	replace n_charter = 49 if endyear==2018
	replace n_charter = 51 if endyear==2019
	replace n_charter = 51 if endyear==2020
	replace n_charter = 52 if endyear==2021
	replace n_charter = 52 if endyear==2022
	replace n_charter = 52 if endyear==2023
	egen n_programs = rowtotal(n_dle_final n_mag n_charter)


	*Plot enrollment in the district against n of choice programs 
   twoway (line n_students endyear , yaxis(1) color(black) ylab(3(2)8)) ///
        (line n_programs endyear, yaxis(2) color(black) lpattern(dash)), ///
        ytitle("Number of Students (100,000s)", axis(1)) ///
        ytitle("Number of Choice Programs", axis(2)) ///
        xtitle("Year") ///
        legend(order(1 "Number of students" 2 "Number of programs" ) row(1) pos(6)) 

    graph export $figures/motivation.pdf, replace

   * Motivating figure 
   twoway (line n_students endyear , yaxis(1) color(black) ylab(3(2)8)) ///
        (line n_programs endyear, yaxis(2) color(white) lpattern(dash)), ///
        ytitle("Number of Students (100,000s)", axis(1)) ///
        ytitle("Number of Choice Programs", axis(2) color(white)) ///
        xtitle("Year") ///
        legend(order(1 "Number of students" 2 "Number of programs" ) row(1) pos(6)) 

	* Dynamic version for slides 
   twoway (line n_students endyear , yaxis(1) color(maroon) ylab(3(2)8))  (line n_programs endyear, yaxis(2) color(white) lpattern(dash)), ytitle("Number of Students (100,000s)", axis(1))  ytitle("Number of Choice Programs", axis(2)  color(white) ) xtitle("Year") legend(order(1 "Number of students" 2 "Number of programs" ) row(1) pos(6))
   graph export $slidefigures/motivation_1.pdf, replace

   twoway (line n_students endyear , yaxis(1) color(maroon) ylab(3(2)8))  (line n_programs endyear, yaxis(2) color(maroon) lpattern(dash)), ytitle("Number of Students (100,000s)", axis(1))  ytitle("Number of Choice Programs", axis(2)  ) xtitle("Year") legend(order(1 "Number of students" 2 "Number of programs" ) row(1) pos(6))
   graph export $slidefigures/motivation_2.pdf, replace
	
end


********************************************************************************
* application_distance_descriptive 
*
* Description:  Gradient of applying wrt to relative distance 
*
********************************************************************************
capture program drop application_distance_descriptive
program define application_distance_descriptive
	use $DTAFINAL/5th_grade_cohorts.dta, clear

	encode localdistrictcode , gen(districtcode)
	egen district_year = group(districtcode endyear)
	binscatter appliedAny rel_nearest_choice_dist, ///
		absorb(district_year ) linetype(qfit) ///
		lcolor(black) mcolor(maroon) ///
		xlabel(-0.5(1)2.5) ///
		ytitle("Share Applying to Any Choice School") ///
		xtitle("Relative Distance to Nearest Choice School")
	graph export $figures/rel_dist_app_any.pdf, replace
end 

capture program drop sampleComparisons 
program define sampleComparisons 

	use "$DTAFINAL/lotteries.dta", clear
	keep if endyear <=2017
	tempfile lotteries 
	save `lotteries', replace 

	use "$DTAFINAL/structural_data_2004_2013.dta", clear 
	destring studentpseudoid, replace 
	tempfile structural
	save `structural', replace

	use $DTAFINAL/5th_grade_cohorts.dta, clear
	keep if endyear <=2017
	replace endyear = endyear + 1 
	merge m:1 studentpseudoid endyear using `lotteries', keep(1 3) gen(mergeLotteries) keepusing(studentpseudoid) 
	replace endyear = endyear - 1 

	merge m:1 studentpseudoid endyear using `structural', keep(1 3) gen(mergeStructural) keepusing(studentpseudoid)


	********************************** Make a table *****************************************************	
	texdoc init $tables/sample_comparison.tex, replace force
	texdoc write \begin{tabular}{lccc} \hline \hline
	texdoc write & Fifth-Grade Students & Lottery Sample & Structural Sample \\
	texdoc write  & (1) & (2) & (3) \\ \hline
	texdoc write & & & \\
	texdoc write \multicolumn{3}{l}{\textbf{Panel A: Student Demographics}} \\
	texdoc write & & & \\
	local demo_vars "hispanic black white asian female poverty el "
	foreach var of local demo_vars {
		
		if "`var'"=="hispanic" local rowname "Hispanic"
		if "`var'"=="black" local rowname "Black"
		if "`var'"=="white" local rowname "White"
		if "`var'"=="asian" local rowname "Asian"

		if "`var'"=="female" local rowname "Female"
		if "`var'"=="poverty" local rowname "Poverty"
		if "`var'"=="el" local rowname "English Learner"
		
		sum `var' 
		local mu_all : di %9.2f round(r(mean), .01)
		sum `var' if mergeLotteries==3
		local mu_m : di %9.2f round(r(mean), .01)
		sum `var' if mergeStructural==3
		local mu_s : di %9.2f round(r(mean), .01)
		texdoc write `rowname' & `mu_all' & `mu_m' & `mu_s' \\

	}

	texdoc write \multicolumn{3}{l}{\textbf{Panel B: Standardized Test scores}} \\
	texdoc write & & & \\
	local tests "z_ela_all z_math_all"
	foreach var of local tests{
		if "`var'"=="z_ela_all" local rowname "Baseline ELA"
		if "`var'"=="z_math_all" local rowname "Baseline Math"

		sum `var' 
		local mu_all : di %9.2f round(r(mean), .01)
		sum `var' if mergeLotteries==3
		local mu_m : di %9.2f round(r(mean), .01)
		sum `var' if mergeStructural==3
		local mu_s : di %9.2f round(r(mean), .01)
		texdoc write `rowname' & `mu_all' & `mu_m' & `mu_s' \\

	}
	
	texdoc write \multicolumn{3}{l}{\textbf{Panel C: Distance to Schools}} \\
	texdoc write & & & \\
	local dist_vars "nearest_schl_dist rel_nearest_choice_dist "
	foreach var of local dist_vars{
		if "`var'"=="nearest_schl_dist" local rowname "Nearest District School"
		if "`var'"=="rel_nearest_choice_dist" local rowname "Relative Dist: Nearest Choice School"
		
		sum `var' 
		local mu_all : di %9.2f round(r(mean), .01)
		sum `var' if mergeLotteries==3
		local mu_m : di %9.2f round(r(mean), .01)
		sum `var' if mergeStructural==3
		local mu_s : di %9.2f round(r(mean), .01)
		texdoc write `rowname' & `mu_all' & `mu_m' & `mu_s' \\

	}
	texdoc write & & & \\
	texdoc write \hline \hline 
	texdoc write \end{tabular}
	texdoc close

end
capture program drop whoApplies
program define whoApplies

	use $DTAFINAL/5th_grade_cohorts.dta, clear
	keep if endyear<=2017
	********************************** Make a table *****************************************************	
	texdoc init $tables/descriptives_v2.tex, replace force
	texdoc write \begin{tabular}{lcccc} \hline \hline
	texdoc write & & & Estimated \\
	texdoc write & Non-Applicants & Applicants & Difference \\
	texdoc write  & (1) & (2) & (3) \\ \hline
	texdoc write & & & \\
	*texdoc write \multicolumn{4}{l}{\textbf{Panel A: Choice School Application and Attendance}} \\
	*texdoc write & & & \\
	/*local panela_vars "appliedAny anyOffer next_schoice"
	*local panela_vars "appliedMagnet magnet_offer"
	foreach var of local panela_vars {
		if "`var'"=="appliedAny" local rowname "Applied to a choice school"
		if "`var'"=="anyOffer" local rowname "Received a choice offer"
		if "`var'"=="next_schoice" local rowname "Enrolled in choice school"

		sum `var' if appliedAny==0
		local mu_all : di %9.2f round(r(mean), .01)
		sum `var' if appliedAny==1
		local mu_m : di %9.2f round(r(mean), .01)
		sum next_schoice if appliedAny==1
		texdoc write `rowname' & `mu_all' & `mu_m' & \\
	}
	*/
	
	count 
	local allStu = r(N)
	count if appliedAny ==0
	local nonMagStu = r(N)
	count if appliedAny ==1 
	local magStu = r(N)
	
	local sharenon : di %4.2f `nonMagStu'/(`nonMagStu' + `magStu')
	local shareapp : di %4.2f `magStu'/(`nonMagStu' + `magStu')
	
	local nonMagStu_comma : di %12.0fc `nonMagStu'
	local magStu_comma : di %12.0fc `magStu'
	local allStu_comma : di %12.0fc `allStu'

	texdoc write Observations & `nonMagStu_comma' & `magStu_comma' & `allStu_comma' \\
	texdoc write Student Share & `sharenon' & `shareapp' & -- \\

	texdoc write \multicolumn{4}{l}{\textbf{Panel A: Student Demographics}} \\
	texdoc write & & & &\\
	local demo_vars "hispanic black white asian female poverty el "
	foreach var of local demo_vars {
		
		if "`var'"=="hispanic" local rowname "Hispanic"
		if "`var'"=="black" local rowname "Black"
		if "`var'"=="white" local rowname "White"
		if "`var'"=="asian" local rowname "Asian"

		if "`var'"=="female" local rowname "Female"
		if "`var'"=="poverty" local rowname "Poverty"
		if "`var'"=="el" local rowname "English Learner"
		
		reghdfe `var' appliedAny, absorb(endyear) vce(cluster preferredlocationcode)
		local df = e(df_r)
		local d : di %9.2f round(_b[appliedAny], .01)
		local d_se : di trim(string(round(_se[appliedAny], .01), "%4.2f"))
		local d_p = (2*ttail(`df', abs(_b[appliedAny]/_se[appliedAny]) ) )
		local d_star ""
		if `d_p' <=.01 local d_star "***"
		else if `d_p'<=.05 & `d_p ' >.01 local d_star "**"
		else if `d_p'<=.10 & `d_p ' >.05 local d_star "*"

		*sum `var' if e(sample)==1 & appliedAny==0
		local mu_all : di %9.2f round(_b[_cons] , .01)
		*sum `var' if appliedAny==1 & e(sample)==1
		local mu_m : di %9.2f round(_b[_cons] + _b[appliedAny], .01)
		
		*reghdfe `var' next_schoice, absorb(endyear) vce(cluster preferredlocationcode)
		*local mu_e : di %9.2f round(_b[_cons] + _b[next_schoice], .01) 

		texdoc write `rowname' & `mu_all' & `mu_m' &  `d'`d_star' \\[-0.5ex]
		texdoc write  &  &  & (`d_se') \\
	}	
	 
	texdoc write \multicolumn{4}{l}{\textbf{Panel B: Standardized Test scores}} \\
	texdoc write & & & \\
	local tests "z_ela_all z_math_all"
	foreach var of local tests{
		if "`var'"=="z_ela_all" local rowname "ELA"
		if "`var'"=="z_math_all" local rowname "Math"

		reghdfe `var' appliedAny, absorb(endyear) vce(cluster preferredlocationcode)
		local df = e(df_r)
		local d : di %9.2f round(_b[appliedAny], .01)
		local d_se : di trim(string(round(_se[appliedAny], .01), "%4.2f"))
		local d_p = (2*ttail(`df', abs(_b[appliedAny]/_se[appliedAny]) ) )
		local d_star ""
		if `d_p' <=.01 local d_star "***"
		else if `d_p'<=.05 & `d_p ' >.01 local d_star "**"
		else if `d_p'<=.10 & `d_p ' >.05 local d_star "*"
		

		*sum `var' if e(sample)==1 & appliedAny==0
		local mu_all : di %9.2f round(_b[_cons] , .01)
		*sum `var' if appliedAny==1 & e(sample)==1
		local mu_m : di %9.2f round(_b[_cons] + _b[appliedAny], .01)

		*reghdfe `var' next_schoice, absorb(endyear) vce(cluster preferredlocationcode)
		*local mu_e : di %9.2f round(_b[_cons] + _b[next_schoice], .01) 

		texdoc write `rowname' & `mu_all' & `mu_m' &  `d'`d_star' \\[-0.5ex]
		texdoc write  &  &  & (`d_se') \\
	}
	
	texdoc write \multicolumn{4}{l}{\textbf{Panel C: Distance to Schools}} \\
	texdoc write & & & \\
	local dist_vars "nearest_schl_dist rel_nearest_choice_dist "
	foreach var of local dist_vars{
		if "`var'"=="nearest_schl_dist" local rowname "Nearest District School"
		if "`var'"=="rel_nearest_choice_dist" local rowname "Relative Dist: Nearest Choice School"
		
		
		reghdfe `var' appliedAny, absorb(gradecode endyear) vce(cluster preferredlocationcode)
		local df = e(df_r)
		local d : di %9.2f round(_b[appliedAny], .01)
		local d_se : di %9.2f round(_se[appliedAny], .01)
		local d_p = (2*ttail(`df', abs(_b[appliedAny]/_se[appliedAny]) ) )
		local d_star ""
		if `d_p' <=.01 local d_star "***"
		else if `d_p'<=.05 & `d_p ' >.01 local d_star "**"
		else if `d_p'<=.10 & `d_p ' >.05 local d_star "*"
		

		*sum `var' if e(sample)==1 & appliedAny==0
		local mu_all : di %9.2f round(_b[_cons] , .01)
		*sum `var' if appliedAny==1 & e(sample)==1
		local mu_m : di %9.2f round(_b[_cons] + _b[appliedAny], .01)

		*reghdfe `var' next_schoice, absorb(endyear) vce(cluster preferredlocationcode)
		*local mu_e : di %9.2f round(_b[_cons] + _b[next_schoice], .01) 

		texdoc write `rowname' & `mu_all' & `mu_m' &  `d'`d_star' \\[0.5ex]
		texdoc write  &  &  & (`d_se') \\
	}

	texdoc write \hline \hline 
	texdoc write \end{tabular}
	texdoc close

end 


********************************************************************************
* distance_correlates 
*
* Description:  Correlates of distance 
*
********************************************************************************
capture program drop distance_correlates
program define distance_correlates

	use $DTAFINAL/5th_grade_cohorts.dta, clear

	********************************** Distance Correlates *****************************************************
	* Distance to nearest choice school 
	local vars "z_math_all z_ela_all female black hispanic white asian poverty el eng_home esp_home"
	*quietly eststo Multivariate: reghdfe rel_nearest_mag_dist `vars' if endyear <=2013 & endyear>=2004, absorb(endyear censusblockid) vce(r) 
  
	local vars_nocontrol "z_math_all_nocontrol z_ela_all_nocontrol female_nocontrol black_nocontrol hispanic_nocontrol white_nocontrol asian_nocontrol poverty_nocontrol el_nocontrol eng_home_nocontrol esp_home_nocontrol"
    stop;
	foreach var in `vars' {
		quietly  eststo `var' : reghdfe nearest_choice_dist    `var',  absorb(endyear censusblockid) vce(cluster censusblockid)
		quietly eststo `var'_nocontrol: reghdfe nearest_choice_dist    `var', noabsorb  vce(robust)
		
	}
	coefplot ( `vars' , label(Changes in Distance) color (black) ciopts(lcolor(black) ) ) ///
		(`vars_nocontrol', label(Distance) color(maroon) ciopts(lcolor(maroon) ) ), drop(_cons) xline(0) legend(pos(6) row(1)) ///
		xscale(range(-0.25 0.45))  xtick(-0.2 0 0.2 0.4) ///
		rename( z_ela_all = "ELA Scores" ///
				z_math_all = "Math Scores" ///
				female = "Female" ///
				black = "Black" ///
				hispanic = "Hispanic" ///
				white = "White" ///
				asian = "Asian" ///
				poverty = "Poverty" ///
				el = "English Learner" ///
				eng_home = "English Speaking Home" ///
				esp_home = "Spanish Speaking Home" ///
				born_usa = "Born in USA" )
	graph export $figures/distance_correlates_nearest.pdf, replace

	********************************** Distance Correlates *****************************************************
	* Relative distance to nearest choice school 
	local vars "z_math_all z_ela_all female black hispanic white asian poverty el eng_home esp_home"
	*quietly eststo Multivariate: reghdfe rel_nearest_mag_dist `vars' if endyear <=2013 & endyear>=2004, absorb(endyear censusblockid) vce(r) 
  
	local vars_nocontrol "z_math_all_nocontrol z_ela_all_nocontrol female_nocontrol black_nocontrol hispanic_nocontrol white_nocontrol asian_nocontrol poverty_nocontrol el_nocontrol eng_home_nocontrol esp_home_nocontrol"
    
	foreach var in `vars' {
		quietly  eststo `var' : reghdfe rel_nearest_choice_dist     `var',  absorb(endyear censusblockid) vce(cluster censusblockid)
		quietly eststo `var'_nocontrol: reghdfe rel_nearest_choice_dist    `var', noabsorb  vce(robust)

	}
	coefplot ( `vars' , label(Changes in Distance) color (black) ciopts(lcolor(black) ) ) ///
		(`vars_nocontrol', label(Distance) color(maroon) ciopts(lcolor(maroon) ) ), drop(_cons) xline(0) legend(pos(6) row(1)) ///
		xscale(range(-0.25 0.45)) xlabel(-0.2 0 0.2 0.4) ///
		rename( z_ela_all = "ELA Scores" ///
				z_math_all = "Math Scores" ///
				female = "Female" ///
				black = "Black" ///
				hispanic = "Hispanic" ///
				white = "White" ///
				asian = "Asian" ///
				poverty = "Poverty" ///
				el = "English Learner" ///
				eng_home = "English Speaking Home" ///
				esp_home = "Spanish Speaking Home" ///
				born_usa = "Born in USA" )
	graph export $figures/distance_correlates_nearest_relative.pdf, replace
end 



********************************************************************************
* distance_maps  
*
* Description:  Calls R script to make distance maps 
*
********************************************************************************
capture program drop distance_maps
program define distance_maps	
	local R "/usr/local/bin/Rscript"
	local Rfile "$ROOT/code/distance_to_choice_maps.R"

	shell `R' --vanilla `Rfile'

end 




capture program drop nextProgram 
program define nextProgram
	********************************** Make a figure based on the table **************************************	
	** Comparison between all students versus applicants 
	
	local demo_vars "hispanic black white asian female poverty el z_ela_all z_math_all"
	foreach var of local demo_vars{
		
		reghdfe `var' appliedMagnet, absorb(gradecode endyear) vce(cluster preferredlocationcode)
		local d_`var' : di %9.2f round(_b[appliedMagnet], .01)
		local d_se_`var' : di %9.2f round(_se[appliedMagnet], .01)
		local mu_all_`var' : di %9.2f round(_b[_cons] , .01)
		local mu_m_`var' : di %9.2f round(_b[_cons] + _b[appliedMagnet], .01)

	}
	
	
	preserve 
	
	* store results as a separate file
	clear 
	set obs 18
	
	* sequence variables
	seq group, f(1) t(2)
	 * all student and applicants place holders
	seq outcome, b(2) f(1) t(9)

	 *generates repeated sequence... 1 1, 2 2,...

	* variables for outcomes and standard errors
	gen y = . 
	gen se = .

	local l = 1

	* bars
	*order of vars matter
	foreach var in hispanic black white asian female poverty el z_ela_all z_math_all   { 
	*all students bar
	  replace y = `mu_all_`var'' if group==1 & outcome==`l'
	 *applicants bar
	  replace y = `mu_m_`var'' if group==2 & outcome==`l' 
	  *SE
	  replace se = `d_se_`var'' if group==2 & outcome==`l' 
	  local l = `l'+1
	  }

	gen y_upper = y+1.96*se
	gen y_lower = y-1.96*se
	
	gen outcomegroup=group if outcome==1
	replace outcomegroup=group+3 if outcome==2
	replace outcomegroup=group+6 if outcome==3
	replace outcomegroup=group+9 if outcome==4
	replace outcomegroup=group+12 if outcome==5
	replace outcomegroup=group+15 if outcome==6
	replace outcomegroup=group+18 if outcome==7
	replace outcomegroup=group+21 if outcome==8
	replace outcomegroup=group+24 if outcome==9

	list
	
	* create components for the figure	
	 graph twoway (bar y outcomegroup if group==1 & outcome<=7, yaxis(1) color(gs10)) (bar y outcomegroup if group==2 & outcome<=7, yaxis(1) color(maroon)), legend(order (1 "All students" 2 "Applicants") row(1) position(6)) fxsize(100)   ylabel(0(0.2)1, angle(0) format(%03.1f))  xlabel(1.5 `"Hispanic"' 4.5 `" Black"' 7.5 `" White"' 10.5 `" Asian"' 13.5 `" Female"' 16.5 `"Poor"' 19.5 `" "English" "Learner" "' , labgap(*1) noticks)  ytitle("Share") xtitle("") 
	 graph save "$figures/figure_stats_part1.gph", replace
	 
		
	graph twoway (bar y outcomegroup if group==1 & outcome>=8, yaxis(1) color(gs10))  (bar y outcomegroup if group==2 & outcome>=8, yaxis(1) color(maroon)), legend(order (1  "All students" 2 "Applicants") row(1) position(6)) fxsize(60)  ylabel(, angle(0) format(%15.0gc))xlabel(22.5 `" "ELA" "' 25.5 `" "Math" "') ytitle("Z score") xtitle("") 
	
	graph save "$figures/figure_stats_part2.gph", replace
	

	*Combine graphs 
	grc1leg "$figures/figure_stats_part1.gph" "$figures/figure_stats_part2.gph", legendfrom("$figures/figure_stats_part1.gph") rows(1)
	graph export "$figures/figure_stats_all_vsapp.pdf", replace as(pdf)
	
	restore
	
	
	********************************** Make a figure on the schools characteristics **************************************	
	preserve
	
	merge 1:1 studentpseudoid endyear using "$BUILD/ses_2016_2023.dta", keepus(  z_happy_index ) gen(mergeSES) 
	drop if mergeSES==2

	keep endyear preferredlocationcode hispanic black white  asian other female  poverty el gifted z_ela_all z_math_all next_smagnet next_scode z_happy_index  
	   
	collapse (mean)  hispanic black white  asian other female  poverty el gifted z_ela_all z_math_all next_smagnet (first) preferredlocationcode z_happy_index, by(endyear next_scode)
	
	collapse (mean)  hispanic black white  asian other female  poverty el gifted z_ela_all z_math_all next_smagnet (first) preferredlocationcode z_happy_index ,by(next_scode)
	
	** Comparison between neighborhood schools and magnet schools
	
	local demo_vars "hispanic black white asian female poverty el z_ela_all z_math_all z_happy_index "
	
	foreach var of local demo_vars{
		
		reghdfe `var' next_smagnet, vce(cluster preferredlocationcode)
		local d_`var' : di %9.2f round(_b[next_smagnet], .01)
		local d_se_`var' : di %9.2f round(_se[next_smagnet], .01)
		local mu_all_`var' : di %9.2f round(_b[_cons] , .01)
		local mu_m_`var' : di %9.2f round(_b[_cons] + _b[next_smagnet], .01)

	}
	
	* store results as a separate file
	clear 
	set obs 20
	
	* sequence variables
	seq group, f(1) t(2)
	 * all student and applicants place holders
	seq outcome, b(2) f(1) t(10)

	 *generates repeated sequence... 1 1, 2 2,...

	* variables for outcomes and standard errors
	gen y = . 
	gen se = .

	local l = 1

	* bars
	*order of vars matter
	foreach var in hispanic black white asian female poverty el z_ela_all z_math_all  z_happy_index  { 
	*all students bar
	  replace y = `mu_all_`var'' if group==1 & outcome==`l'
	 *applicants bar
	  replace y = `mu_m_`var'' if group==2 & outcome==`l' 
	  *SE
	  replace se = `d_se_`var'' if group==2 & outcome==`l' 
	  local l = `l'+1
	  }

	gen y_upper = y+1.96*se
	gen y_lower = y-1.96*se
	
	gen outcomegroup=group if outcome==1
	replace outcomegroup=group+3 if outcome==2
	replace outcomegroup=group+6 if outcome==3
	replace outcomegroup=group+9 if outcome==4
	replace outcomegroup=group+12 if outcome==5
	replace outcomegroup=group+15 if outcome==6
	replace outcomegroup=group+18 if outcome==7
	replace outcomegroup=group+21 if outcome==8
	replace outcomegroup=group+24 if outcome==9
	replace outcomegroup=group+27 if outcome==10
	

	list
	
	
	* create components for the figure	
	 graph twoway (bar y outcomegroup if group==1 & outcome<=7, yaxis(1) color(gs10))  (bar y outcomegroup if group==2 & outcome<=7, yaxis(1) color(maroon)) , legend(order (1 "Neighborhood schools" 2 "Choice schools") row(1) position(6)) fxsize(100)   ylabel(0(0.2)1, angle(0) format(%03.1f))  xlabel(1.5 `"Hispanic"' 4.5 `" Black"' 7.5 `" White"' 10.5 `" Asian"' 13.5 `" Female"' 16.5 `"Poor"' 19.5 `" "English" "Learner" "' , labgap(*1) noticks)  ytitle("Share") xtitle("")
	 
	 graph save "$figures/figure_stats_schools_part1.gph", replace
		
	graph twoway (bar y outcomegroup if group==1 & outcome>=8 & outcome<=9, yaxis(1) color(gs10))  (bar y outcomegroup if group==2 & outcome>=8 & outcome<=9 , yaxis(1) color(maroon)) , legend(order (1 "Neighborhood schools" 2 "Choice schools") row(1) position(6)) fxsize(60)  ylabel(, angle(0) format(%15.0gc)) xlabel(22.5 `" "ELA" "' 25.5 `" "Math" "') ytitle("Z score") xtitle("")
	 
	graph save "$figures/figure_stats_schools_part2.gph", replace
	
	graph twoway (bar y outcomegroup if group==1 & outcome>=10, yaxis(1) color(gs10))  (bar y outcomegroup if group==2 & outcome>=10, yaxis(1) color(maroon)), legend(order (1 "Neighborhood schools" 2 "Choice schools") row(1) position(6)) fxsize(60)  ylabel(0(0.2)0.6, angle(0) format(%15.0gc)) xlabel( 28.5 `" "Happy" "') ytitle("Z score") xtitle("")
	 
	graph save "$figures/figure_stats_schools_part3.gph", replace

	*Combine graphs 
	grc1leg "$figures/figure_stats_schools_part1.gph" "$figures/figure_stats_schools_part2.gph" "$figures/figure_stats_schools_part3.gph", legendfrom("$figures/figure_stats_schools_part1.gph") rows(1)
	graph export "$figures/figure_stats_choices.pdf", replace as(pdf)
	
	restore
	

end 



	
main 
