/********************************************************************************
0_build.do 

- Authors: Antonia, Vazquez and Chris Campos
- Description: Creates analyses datasets for decentralized choice project. 
    (i) Lottery data - This includes information at the student level on 
    demographics, attendance, GPA, test scores, suspensions, and SES. 
    Also includes information at the school level on HHI, VAM, and peer quality 
    both for the school the student applied to and the one the student enrolled in 
    (ii) Post-model estimation analyses data - after estimating the structural 
    model, this dataset includes information on all students, not just lottery
    applicants along with relevant information from the model 
	(iii) Block- and tract-by-year data for quasi-experimental event study analysis 
- Dependencies: This project uses the following raw data 
    (i) Lottery records from magnet project
    (ii) Distance calculations from aux_data_make in the build folder 
    (iii) TBD 
- Date started: 08/13/2024
- Last update: 8/23/2025
* Packages to install: getcensus, geonear
*******************************************************************************/
clear all 
set trace on 
set tracedepth 2
set maxvar 120000

capture program drop main 
program define main 

	set_paths 
	
	/* Prepare Lottery Data */ 
	*prep_relative_distance
	*prep_relative_distance_choices
	
	*prep_st_info

	*school_demo
	*hhi
	*combineVA
    
	*PeerQ_school
	*ses_school
	
    *merge_st_info
	*merge_other_choices
	
	/* Prepare data for descriptives */ 
	*prep_descriptives

	prep_lottery

	*prep_event_study_tract, relative(yes)
	*prep_event_study_tract, relative(no)
	*prep_event_study_block
end

********************************************************************************
* set_paths 
* 
* Description: Establish paths 
********************************************************************************
capture program drop set_paths
program define set_paths

	* Establish directories 
	if "`c(os)'"=="MacOSX" & "`c(username)'"=="ccampos"{
		global ROOT "/Volumes/lausd/decentralized_choice"
		global BUILD "/Volumes/lausd/build"

	}
	else{
		global ROOT "/project/lausd/decentralized_choice"
		global BUILD "/project/lausd/build"
	}
	
	global RAWDATA "$BUILD/rawdata/Campos_Master_DUA" 
	
    * Data directories for decentralized choice folder 
	global DTARAW "$ROOT/data/raw" 
	global DTAVAM "$ROOT/rawdata/vam/data/intermediate"
	global DTAINT "$ROOT/data/intermediate" 
	global DTAFINAL "$ROOT/data"
	global LOGS		"$ROOT/logs"
  
	
end


********************************************************************************
* prep_relative_distance
* 
* Description: Prepares a relative distance dataset that is merged in later 
* 
* Dependencies: This project depends on relative distances calculated for 
* every census block in the county. You can find the program in 
* aux_data_make in the lausd build folder. 
********************************************************************************/
capture program drop prep_relative_distance 
program define prep_relative_distance

	use $BUILD/rawdata/aux_data/relative_magnet_distances_grade2.dta, clear
	append using $BUILD/rawdata/aux_data/relative_magnet_distances_grade5
	append using $BUILD/rawdata/aux_data/relative_magnet_distances_grade8
	
	keep censusblockid grade_level  rank_closest_mag min_distance min_mag_dist mean_mag_dist20 rel_mean_mag_dist20 rel_mag_distance rel_mag2_distance min2_mag_dist endyear os_total_150 os2_total_150 os3_total_150 min_os_mag_dist min_os2_mag_dist min_os3_mag_dist rel_os_mag_dist rel_os2_mag_dist rel_os3_mag_dist os_total_20 os2_total_20 os3_total_20 

	rename (rank_closest_mag min_distance min_mag_dist  ) (nearest_mag_rank nearest_schl_dist nearest_mag_dist) 
	rename (mean_mag_dist20 rel_mean_mag_dist20 rel_mag_distance) (avg_mag_dist20 rel_avg_mag_dist20 rel_nearest_mag_dist)
	rename (rel_mag2_distance min2_mag_dist) (rel_second_mag_dist second_mag_dist )

    label var nearest_mag_rank "Distance rank of nearest school"	
    label var nearest_schl_dist "Distance between block and nearest school (miles)"
    label var nearest_mag_dist "Distance between block and nearest magnet (miles)"
    label var avg_mag_dist20 "Average distance to all magnets among the 20 closest schools (miles)"
    label var rel_avg_mag_dist20 "avg_mag_dist20 relative to nearest school (miles)"
    label var rel_nearest_mag_dist "Nearest magnet relative distance (miles)"
    label var rel_second_mag_dist "Relative distance to second-nearest magnet (miles)"
    label var second_mag_dist "Distance between block and second-nearest magnet (miles)"
    label var os_total_150 "Number of oversubscribed magnets among the 50 closest schools to census block"
    label var os2_total_150 "Number of oversubscribed w 2x apps to seats among the 150 closest schools to census block"
    label var os3_total_150 "Number of oversubscribed w 3x apps to seats among the 150 closest schools to census block"
    label var min_os_mag_dist "Distance to nearest oversubscribed magnet (miles)"
    label var min_os2_mag_dist "Distance to nearest oversubscribed w 2x apps to seats magnet (miles)"
    label var min_os3_mag_dist "Distance to nearest oversubscribed w 3x apps to seats magnet (miles)"
    label var rel_os_mag_dist "Relative distance to nearest oversubscribed magnet (miles)"
    label var rel_os2_mag_dist "Relative distance to nearest oversubscribed w 2x apps to seats magnet (miles)"
    label var rel_os3_mag_dist "Relative distance to nearest oversubscribed w 3x apps to seats magnet (miles)"
    label var os_total_20 "Number of oversubscribed magnets among the 20 closest schools to census block"
    label var os2_total_20 "Number of oversubscribed w 2x apps to seats among the 20 closest schools to census block"
    label var os3_total_20 "Number of oversubscribed w 3x apps to seats among the 20 closest schools to census block"

	gen schooltype = "ele" if grade_level==2
	replace schooltype = "ms" if grade_level==5
	replace schooltype= "hs" if grade_level==8
	drop grade_level
	save $DTAINT/magnet_distances.dta , replace 

end 


*****************************************************************************
* prep_relative_distance_choices
* 
* Description: Prepares a relative distance dataset that is merged in later 
********************************************************************************/
capture program drop prep_relative_distance_choices 
program define prep_relative_distance_choices

	use $BUILD/rawdata/aux_data/relative_choice_distances_grade2.dta, clear
	append using $BUILD/rawdata/aux_data/relative_choice_distances_grade5
	append using $BUILD/rawdata/aux_data/relative_choice_distances_grade8
	
	keep censusblockid grade_level  rank_closest_choice min_distance min_choice_dist mean_choice_dist20 rel_mean_choice_dist20 rel_choice_distance rel_choice2distance min2_choice_dist endyear             
	
	rename (rank_closest_choice min_distance min_choice_dist  ) (nearest_choice_rank nearest_schl_dist nearest_choice_dist) 
	rename (mean_choice_dist20 rel_mean_choice_dist20 rel_choice_distance) (avg_choice_dist20 rel_avg_choice_dist20 rel_nearest_choice_dist)
	rename (rel_choice2distance min2_choice_dist) (rel_second_choice_dist second_choice_dist)
                                                   
		label var nearest_choice_rank "Distance rank of nearest school"	
		label var nearest_schl_dist "Distance between block and nearest school (miles)"
		label var nearest_choice_dist "Distance between block and nearest choice (miles)"
		label var avg_choice_dist20 "Average distance to all choice among the 20 closest schools (miles)"
		label var rel_avg_choice_dist20 "avg_choice_dist20 relative to nearest school (miles)"
		label var rel_nearest_choice_dist "Nearest choice relative distance (miles)"
		label var rel_second_choice_dist "Relative distance to second-nearest choice (miles)"
		label var second_choice_dist "Distance between block and second-nearest choice (miles)"
		

	gen schooltype = "ele" if grade_level==2
	replace schooltype = "ms" if grade_level==5
	replace schooltype= "hs" if grade_level==8
	drop grade_level
	
	save $DTAINT/choice_distances.dta , replace 


end 


********************************************************************************
* prep_st_info
* 
* Description: prepares student data for merge to lotteries with cut-off data 
********************************************************************************/
capture program drop prep_st_info 
program define prep_st_info 
	cap log close
	log using "$LOGS/prep_st_info.log", text replace
	
	use endyear studentpseudoid preferredlocationcode preferredlocationname schoollocationcode schoollocationname schoolcostcentercode  gradecode gendercode studentclassofname studentbirthcountry parentedulevelname studentusschoolfirstattenddate ethnicitydescription languageclasscode homelanguagedescription studentspedflag studentgiftedprogramdescription studentpovertyindicator eng_home esp_home born_usa using $BUILD/output/data/lausd2001_2023, clear 

	* Want to keep track of whether or not we observe students in high school
	gen inHS = gradecode==9 | gradecode==10 | gradecode == 11 | gradecode == 12
	bys studentpseudoid: egen observe_in_HS = max(inHS)
	label var observe_in_HS "Observed in high school"
	drop inHS 
	
	gen in9 = gradecode==9
	bys studentpseudoid: egen observe_in_9 = max(in9)
	label var observe_in_9 "Observed in 9th grade"
	drop in9
	
	label var endyear "School year"
	
	destring  studentpseudoid, replace
	bysort endyear studentpseudoid: gen dup = cond(_N==1,0,_n)
	drop if dup>0
	drop dup
	
	
	save "$DTAINT/demographics.dta", replace
	
	*Attendance
	use studentpseudoid endyear ytdofattendeddays ytdofenrolleddays ytdofabsentdays otherabsencecount ytdofattendance appabsencecount using $BUILD/output/data/lausd2001_2023, clear 	
	label var endyear "School year"
	
	destring  studentpseudoid, replace
	bysort endyear studentpseudoid: gen dup = cond(_N==1,0,_n)
	drop if dup>0
	drop dup
	save "$DTAINT/att.dta", replace
	

	*GPA 
	use studentpseudoid gradecode endyear gpa_fall gpa_spring ///
		using $BUILD/output/data/lausd2001_2023, clear 
	
	label var endyear "School year"
	
	destring  studentpseudoid, replace
	bysort endyear studentpseudoid gradecode: gen dup = cond(_N==1,0,_n)
	drop if dup>0
	drop dup
	destring gpa_fall gpa_spring, replace
	format studentpseudoid %20.0g
	*egen gpa = rowmean(gpa_fall gpa_spring)
	
	/*
	count if missing( gpa_fall )
	count if missing( gpa_spring )
	reshape long gpa, i(studentpseudoid gradecode endyear) j(term) string 
	format %ty endyear
	gen hy_gpa = .
	replace hy_gpa = 1 if term=="_fall"
	replace hy_gpa = 2 if term=="_spring"
	gen endyear_halfyear = yh(endyear,hy_gpa)
	format endyear_halfyear %th
	drop if missing(gpa)
	*/
	
	* The year of reference here is the academic year the student is on  (endyear)
	* The match after is going to be using the variable endyear
	* within an endyear,  fall GPA is measured in december and spring GPA is measured in june
	*We won't consider the fall GPA or spring GPA of the application year because the GPA is measured after applying
	* That is why the lag is two periods before 
	
	xtset studentpseudoid endyear 
	
	gen F1GPA_s = gpa_spring // leads
	gen F2GPA_s = F1.gpa_spring
	gen F3GPA_s = F2.gpa_spring 
	gen F4GPA_s = F3.gpa_spring
	gen L1GPA_s = L2.gpa_spring // lag
	gen L2GPA_s = L3.gpa_spring
	
	gen F1GPA_f = gpa_fall // lead
	gen F2GPA_f = F1.gpa_fall
	gen F3GPA_f = F2.gpa_fall 
	gen F4GPA_f = F3.gpa_fall
	gen L1GPA_f = L2.gpa_fall // lag
	gen L2GPA_f = L3.gpa_fall

	drop gpa_fall gpa_spring 

	label var F1GPA_s "Spring GPA (year after application)"
	label var F2GPA_s "Spring GPA (2 years after application)"
	label var F3GPA_s "Spring GPA (3 years after application)"
	label var F4GPA_s "Spring GPA (4 years after application)"
	label var L1GPA_s "Spring GPA (year before application)"
	label var L2GPA_s "Spring GPA (2 years before application)"

	save "$DTAINT/gpa.dta", replace
	
	
	* Student outcomes (test scores and suspensions)
	use studentpseudoid gradecode endyear z_ela_all z_math_all ///
		total_num_suspended_days total_num_suspensions ///
		using $BUILD/output/data/lausd2001_2023, clear
	
	destring  studentpseudoid, replace
	format studentpseudoid %14.0g
	replace total_num_suspended_days =0 if missing(total_num_suspended_days)
	replace total_num_suspensions	=0 if missing(total_num_suspensions)
	
	* The year of reference here is the academic year the student is on  (endyear)
	* The match after is going to be using the variable endyear
	*We won't consider the score of the application year because students are tested right after applying
	* That is why the lag is two periods before 
	xtset studentpseudoid endyear 
	gen F1math = z_math_all // lead
	gen F2math = F1.z_math_all
	gen F3math = F2.z_math_all 
	gen L1math = L2.z_math_all // lag
	gen L2math = L3.z_math_all 
	gen L3math = L4.z_math_all 
	
	gen F1ela = z_ela_all 
	gen F2ela = F1.z_ela_all 
	gen F3ela = F2.z_ela_all
	gen L1ela = L2.z_ela_all
	gen L2ela = L3.z_ela_all 
	gen L3ela = L4.z_ela_all 
	
	gen F1numsuspensions = total_num_suspensions
	gen F2numsuspensions = F1.total_num_suspensions
	gen F3numsuspensions = F2.total_num_suspensions 
	gen L1numsuspensions = L2.total_num_suspensions
	gen L2numsuspensions = L3.total_num_suspensions 
	gen L3numsuspensions = L4.total_num_suspensions
	
	drop total_num_suspended_days total_num_suspensions z_ela_all z_math_all

	label var F1math "Math score (year after application)"
	label var F2math "Math score (2 years after application)"
	label var F3math "Math score (3 years after application)"
	label var L1math "Math score (year before application)"
	label var L2math "Math score (2 years before application)"
	label var L3math "Math score (3 years before application)"

	label var F1ela "ELA score (year after application)"
	label var F2ela "ELA score (2 years after application)"
	label var F3ela "ELA score (3 years after application)"
	label var L1ela "ELA score (year before application)"
	label var L2ela "ELA score (2 years before application)"
	label var L3ela "ELA score (3 years before application)"

	label var F1numsuspensions "Number of suspensions (year after application)"
	label var F2numsuspensions "Number of suspensions (2 years after application)"
	label var F3numsuspensions "Number of suspensions (3 years after application)"
	label var L1numsuspensions "Number of suspensions (year before application)"
	label var L2numsuspensions "Number of suspensions (2 years before application)"
	label var L3numsuspensions "Number of suspensions (3 years before application)"

	save "$DTAINT/test_and_suspensions.dta", replace
	
	***SES data student 
	use "$BUILD/output/data/ses_2016_2023.dta", clear 
	
	keep endyear studentpseudoid z_happy_index z_happy_index_qs z_bully_index z_bully_index_qs z_no_bully_index z_no_bully_index_qs z_climate_index z_climate_index_qs z_extra_index z_extra_index_qs z_expectation_index z_expectation_index_qs z_perceptions_index z_perceptions_index_qs z_interpersonal_index z_interpersonal_index_qs z_connectedness_index z_connectedness_index_qs z_effort_index z_effort_index_qs z_grit_index z_grit_index_qs
	
	* The year of reference here is the academic year the student is on  (endyear)
	* The match after is going to be using the variable endyear
	*We won't consider the survey of the application year 
	* That is why the lag is two periods before 
	xtset studentpseudoid endyear 
	
	
	local vars "z_happy_index z_bully_index  z_no_bully_index  z_climate_index  z_extra_index  z_expectation_index  z_perceptions_index  z_interpersonal_index  z_connectedness_index  z_effort_index  z_grit_index"
	
	foreach var of local vars {
	gen F1`var' = `var' // lead
	gen F2`var' = F1.`var'
	gen F3`var' = F2.`var' 
	gen L1`var' = L2.`var' // lag
	gen L2`var' = L3.`var'
	gen L3`var' = L4.`var'

	}

	drop z_happy_index z_happy_index_qs z_bully_index z_bully_index_qs z_no_bully_index z_no_bully_index_qs z_climate_index
	
	label var F1z_happy_index "Happy index (year after application)"
	label var F2z_happy_index "Happy index (2 years after application)"
	label var F3z_happy_index "Happy index (3 years after application)"
	label var L1z_happy_index "Happy index (year before application)"
	label var L2z_happy_index "Happy index (2 years before application)"
	label var L3z_happy_index "Happy index (3 years before application)"
	
	label var F1z_bully_index "Bully index (year after application)"
	label var F2z_bully_index "Bully index (2 years after application)"
	label var F3z_bully_index "Bully index (3 years after application)"
	label var L1z_bully_index "Bully index (year before application)"
	label var L2z_bully_index "Bully index (2 years before application)"
	label var L3z_bully_index "Bully index (3 years before application)"

	label var F1z_no_bully_index "No bully index (year after application)"
	label var F2z_no_bully_index "No bully index (2 years after application)"
	label var F3z_no_bully_index "No bully index (3 years after application)"
	label var L1z_no_bully_index "No bully index (year before application)"
	label var L2z_no_bully_index "No bully index (2 years before application)"
	label var L3z_no_bully_index "No bully index (3 years before application)"

	label var F1z_climate_index "School Climate index (year after application)"
	label var F2z_climate_index "School Climate index (2 years after application)"
	label var F3z_climate_index "School Climate index (3 years after application)"
	label var L1z_climate_index "School Climate index (year before application)"
	label var L2z_climate_index "School Climate index (2 years before application)"
	label var L3z_climate_index "School Climate index (3 years before application)"

	label var F1z_extra_index "Extracurricular Activities index (year after application)"
	label var F2z_extra_index "Extracurricular Activities index (2 years after application)"
	label var F3z_extra_index "Extracurricular Activities index (3 years after application)"
	label var L1z_extra_index "Extracurricular Activities index (year before application)"
	label var L2z_extra_index "Extracurricular Activities index (2 years before application)"
	label var L3z_extra_index "Extracurricular Activities index (3 years before application)"

	label var F1z_expectation_index "Expectation index (year after application)"
	label var F2z_expectation_index "Expectation index (2 years after application)"
	label var F3z_expectation_index "Expectation index (3 years after application)"
	label var L1z_expectation_index "Expectation index (year before application)"
	label var L2z_expectation_index "Expectation index (2 years before application)"
	label var L3z_expectation_index "Expectation index (3 years before application)"

	label var F1z_perceptions_index "Perceptions index (year after application)"
	label var F2z_perceptions_index "Perceptions index (2 years after application)"
	label var F3z_perceptions_index "Perceptions index (3 years after application)"
	label var L1z_perceptions_index "Perceptions index (year before application)"
	label var L2z_perceptions_index "Perceptions index (2 years before application)"
	label var L3z_perceptions_index "Perceptions index (3 years before application)"

	label var F1z_interpersonal_index "Interpersonal skills index (year after application)"
	label var F2z_interpersonal_index "Interpersonal skills index (2 years after application)"
	label var F3z_interpersonal_index "Interpersonal skills index (3 years after application)"
	label var L1z_interpersonal_index "Interpersonal skills index (year before application)"
	label var L2z_interpersonal_index "Interpersonal skills index (2 years before application)"
	label var L3z_interpersonal_index "Interpersonal skills index (3 years before application)"

	label var F1z_connectedness_index "School Connectedness index (year after application)"
	label var F2z_connectedness_index "School Connectedness index (2 years after application)"
	label var F3z_connectedness_index "School Connectedness index (3 years after application)"
	label var L1z_connectedness_index "School Connectedness index (year before application)"
	label var L2z_connectedness_index "School Connectedness index (2 years before application)"
	label var L3z_connectedness_index "School Connectedness index (3 years before application)"


	save "$DTAINT/ses_2016_2023.dta", replace
	
	log close
	
	
end


********************************************************************************
* school_demo
* 
* Description: Prepares school level demographics
*
* Inspected by Chris on 02/03/2025 
*
********************************************************************************/
capture program drop school_demo 
program define school_demo 
	
	use  schoollocationname schoolcostcentercode endyear parentedulevelname studentspedflag studentgiftedprogramdescription studentpovertyindicator born_usa gendercode languageclasscode using "$BUILD/output/data/lausd2001_2023", clear 
		
	*labels 1 "Not HS Grad" 2 "HS Grad" 3 "Some College" 4 "College Grad" 5 "Grad Sch/ Post Grad" 6 "Decline to Answer" 7 " "
	gen college_grad = 0
	replace college_grad=1 if parentedulevelname==4 | parentedulevelname==5
	
	gen special_edu = studentspedflag=="Y"
	
	gen gifted = !missing(studentgiftedprogramdescription) 
	
	gen poverty =  studentpovertyindicator=="Y"
	
	gen immigrant = born_usa==0
	
	*label: 1 "F" 2 "M" 3 "UNKNOWN" 4 "N" 5 " "
	gen female = gendercode==1
	
	gen english_learner=0
	replace english_learner = 1 if languageclasscode == "LEP"


	collapse (first) schoollocationname (mean) college_grad special_edu gifted poverty immigrant female english_learner, by(endyear schoolcostcentercode)
	drop if missing(schoolcostcentercode)
	
	unique schoolcostcentercode
	* There are 1,566 programs 
	
	*Calculate the first year, age of a program, and if it is an entrant program
	sort schoolcostcentercode endyear, stable 
	bysort schoolcostcentercode: gen first_year = endyear if _n == 1
	bysort schoolcostcentercode: egen fst_yr_LAUSD =max(first_year)
	drop first_year
	
	
	gen mag_cd_ = schoolcostcentercode
	gen app_endyear= endyear-1
	
	*Bringing in the information on openings of magnet 
	merge m:1 schoolcostcentercode using "$ROOT/rawdata/magnet_chronology.dta", gen(mergeAge) keepusing(year_opened loc_name)
	drop if mergeAge==2
	
	rename year_opened fst_yr_mag_info 

	*rename variables for enrolled adn applied school
	local vars "college_grad special_edu gifted poverty immigrant female english_learner fst_yr_LAUSD fst_yr_mag_info"
	foreach var of local vars{
		gen sch_`var'_app = `var'
		label var sch_`var'_app "Applied school `var'" 
		rename `var' sch_`var'_enr
		label var sch_`var'_enr "Enrolled school `var'" 
	}
	
	
	save "$DTAINT/school_demo.dta", replace
	
end 

********************************************************************************
* hhi
* 
* Description: Prepares hhi measure
********************************************************************************/
capture program drop hhi 
program define hhi 
	
	use   schoolcostcentercode endyear ethnicitydescription using "$BUILD/output/data/lausd2001_2023", clear 
	gen hispanic = ethnicitydescription ==1 
	gen black = ethnicitydescription ==2
	gen white = ethnicitydescription ==3 
	gen asian = ethnicitydescription==4 
	gen filipino = ethnicitydescription==5
	gen other = ethnicitydescription>5

	collapse (mean) hispanic black white asian filipino other, by(endyear schoolcostcentercode)
	drop if missing(schoolcostcentercode)
	local vars "hispanic black white asian filipino other"
	foreach var of local vars {
		gen `var'2 = `var'^2 
	}
	egen hhi = rowtotal(*2)
	duplicates drop schoolcostcentercode endyear, force
	keep schoolcostcentercode hhi hispanic black white asian filipino other endyear 
	
	*rename variables for enrolled school
	local vars "hispanic black white asian filipino other hhi"
	foreach var of local vars{
		gen sch_`var'_app = `var'
		label var sch_`var'_app "Applied school `var'" 
		rename `var' sch_`var'_enr
		label var sch_`var'_enr "Enrolled school `var'" 
	}
	
	gen mag_cd_ = schoolcostcentercode
	gen app_endyear= endyear-1
	
	save "$DTAINT/hhi.dta", replace
	
end 
********************************************************************************
* combineVA
* 
* Description: Prepares VA and peer quality data created by Connor Fogal
********************************************************************************/
capture program drop combineVA 
program define combineVA 
	
	*Vam elementary using data Connor prepared
	use "$DTAVAM/ele_vam_raw.dta", clear 
	sort year schoolcostcentercode
	*586 unique schools 
	rename year endyear 
	label var endyear "School year"
	
	gen schooltype = "ele" 
	*reshape wide vam stderr, i(schoolcostcentercode) j(subject) string	
	save "$DTAINT/vam.dta", replace 
	
	*Vam middle school using data Connor prepared
	use "$DTAVAM/ms_vam_raw.dta", clear
	sort year schoolcostcentercode
	*301 unique schools 
	rename year endyear 
	label var endyear "School year"
	
	gen schooltype = "ms" 
		
	append using "$DTAINT/vam.dta"
	
	save "$DTAINT/vam.dta", replace
	
	*Vam high school using data Connor prepared
	use "$DTAVAM/hs_vam_raw.dta", clear
	sort year schoolcostcentercode
	*358 unique schools 
	
	rename year endyear 
	label var endyear "School year"
	
	gen schooltype = "hs"
	*drop stderr  
	*reshape wide vam ahat, i(schoolcostcentercode endyear) j(subject) string 
	*rename (vamhatmath vamhatela) (vamhat_math_hs vamhat_ela_hs) 
	*rename (ahatmath ahatela) (ahat_math_hs ahat_ela_hs)
	append using "$DTAINT/vam.dta"
	
	save "$DTAINT/vam.dta", replace
	
	*Special VAM years
	* We don't have the vam and PeerQ estimates for 2014 (change of exams), 2020, and 2021 (covid) due to missing exams those years. Solution: for 2014, let' add the measures from the same school in 2013. For 2020 and 2021, assign 2019.  
	use  "$DTAINT/vam.dta", clear
	
	keep if inlist(endyear, 2013, 2019)
	replace endyear = 2014 if endyear==2013
	
	replace endyear = 2020 if endyear==2019
	
	save "$DTAINT/vam_spec_years.dta", replace
	
	use  "$DTAINT/vam.dta", clear
	keep if inlist(endyear, 2019)
	replace endyear = 2021 if endyear==2019
	append using "$DTAINT/vam_spec_years.dta"
	save "$DTAINT/vam_spec_years.dta", replace
	
	*All years together including special years 
	use  "$DTAINT/vam.dta", clear
	append using "$DTAINT/vam_spec_years.dta"
	save "$DTAINT/vam_WITHspec_years.dta", replace
	
	*bring the data on the school center code
	* Don't need this anymore because VAM at schoolcostecentercode level
	*use  schoollocationcode schoolcostcentercode endyear using "$BUILD/output/data/lausd2001_2023.dta", clear 
	*duplicates drop endyear schoollocationcode schoolcostcentercode, force
	*drop if missing(schoolcostcentercode)
	*duplicates drop schoollocationcode schoolcostcentercode, force
	*drop endyear
	*tempfile sc_codes 
	*save `sc_codes'
	
	*use "$DTAINT/vam_WITHspec_years.dta", clear 
	*merge m:1 schoollocationcode using `sc_codes'
	*drop if _merge==2 |  _merge==1
	*drop _merge
	
	*gen variables to do the merge for the magnet application school
	local vars "vamhat_ela_hs ahat_ela_hs vamhat_math_hs ahat_math_hs vamhat_ela_ms ahat_ela_ms vamhat_math_ms ahat_math_ms vamhat_ela_ele ahat_ela_ele vamhat_math_ele ahat_math_ele"
	foreach var of local vars{
		gen sch_`var'_app = `var'
		label var sch_`var'_app "Applied school `var'" 
		rename `var' sch_`var'_enr
		label var sch_`var'_enr "Enrolled school `var'" 
	}
	
	gen mag_cd_ = schoolcostcentercode
	gen app_endyear= endyear-1
	save "$DTAINT/vam_WITHspec_years.dta", replace
	
end 

	
********************************************************************************
* PeerQ_school
* 
* Description: Prepares peer quality at the school level 
********************************************************************************/
capture program drop PeerQ_school 
program define PeerQ_school 

* Calculation of Peer quality
	use schoollocationcode schoollocationname schoolcostcentercode endyear z_ela_all z_math_all using "$BUILD/output/data/lausd2001_2023", clear 
	collapse (mean) z_ela_all z_math_all, by(endyear schoolcostcentercode)
	drop if missing(schoolcostcentercode)
	label var z_ela_all "PeerQ ELA - mean school's test score"
	label var z_math_all "PeerQ Math - mean school's test score"
	
	*gen a variable to do the merge for the magnet application school
	local vars "z_ela_all z_math_all"
	foreach var of local vars{
		gen sch_`var'_app = `var'
		label var sch_`var'_app "Applied school `var'" 
		rename `var' sch_`var'_enr
		label var sch_`var'_enr "Enrolled school `var'" 
	}
	gen mag_cd_ = schoolcostcentercode
	gen app_endyear=endyear-1
	
	save "$DTAINT/peerQ.dta", replace

end 

********************************************************************************
* ses_school
* 
* Description: Prepares SES data at the school level 
********************************************************************************/
capture program drop ses_school 
program define ses_school 
	use "$BUILD/output/data/ses_2016_2023.dta", clear 
	
	*Keep relevant variables 
	keep endyear schoolcostcentercode z_happy_index z_bully_index z_no_bully_index z_climate_index z_extra_index z_expectation_index z_perceptions_index z_interpersonal_index z_connectedness_index z_effort_index z_grit_index
	
	*Collapse infotmation at the year and school level 
	collapse  z_happy_index z_bully_index z_no_bully_index z_climate_index z_extra_index z_expectation_index z_perceptions_index z_interpersonal_index z_connectedness_index z_effort_index z_grit_index, by(endyear schoolcostcentercode)
	
	drop if missing(schoolcostcentercode)
	
	local vars "z_happy_index z_bully_index z_no_bully_index z_climate_index z_extra_index z_expectation_index z_perceptions_index z_interpersonal_index z_connectedness_index z_effort_index z_grit_index"
	foreach var of local vars {
		gen sch_`var'_app = `var'
		label var sch_`var'_app "Applied school `var'" 
		rename `var'  sch_`var'_enr
		label var sch_`var'_enr "Enrolled school `var'"
	}
	
	
	gen mag_cd_ = schoolcostcentercode
	gen app_endyear=endyear-1
	save "$DTAINT/ses_school.dta", replace


end
		
	
********************************************************************************
* merge_st_info
* 
* Description: merges clean student data to lotteries with cut-off data 
* The year of merging is the enrollment year, except for the information regarding
* the application school 
********************************************************************************
capture program drop merge_st_info 
program define merge_st_info 
	cap log close
	log using "$LOGS/merge_st_info.log", text replace
	
	use "$BUILD/output/data/magnet_cutoffs/potential_lotteries_WITHcutoffs.dta", clear 
	
	*Keep only lotteries with an estimated cutoff
	keep if lottery_nocutoff==0
	
	gen schooltype = "ele" if stu_grade<=5
	replace schooltype = "ms" if inrange(stu_grade, 6, 8)
	replace schooltype = "hs" if stu_grade>=9

	
	
	destring studentpseudoid, replace 
	
	gen endyear=app_endyear+1
	label var endyear "End year of enrrollment year"
	
	****** STUDENT INFO *******
	
	****DEMOGRAPHICS
	*The merge is using endyear, which is the year when the studen would start in LAUSD
	merge 1:1 endyear studentpseudoid using "$DTAINT/demographics.dta", gen(mergeDem)
	drop if mergeDem==2 // students that were not applicants
	/*Check mergeDem==1 these are students who applied but are not is lausd data. Most of them correspond 
	to application year 2023 and are coming from outside district or have missing current school info
	*/
	destring studentpseudoid preferredlocationcode schoollocationcode schoolcostcentercode gradecode, replace 


	***ATTENDANCE
	merge m:1 studentpseudoid endyear  using "$DTAINT/att.dta", gen(mergeATT)
	drop if mergeATT==2 // students that were not applicants
	
	***GPA
	*Because some students (102,746) change grades in the same year, we also use the grade to keep the one that matches the demographics data
	destring gradecode, replace
	merge m:1 studentpseudoid endyear gradecode using "$DTAINT/gpa.dta", gen(mergeGPA)
	
	drop if mergeGPA==2 // students that were not applicants
	
	***Test scores and suspensions data 
	merge m:1 studentpseudoid endyear  using "$DTAINT/test_and_suspensions.dta", gen(mergeTests)
	drop if mergeTests==2
	
	***SES data student 
	merge 1:1 studentpseudoid endyear using "$DTAINT/ses_2016_2023.dta", gen(mergeSES) keepus( F1* F2* F3* L1* L2* L3*)
	drop if mergeSES==2
	
	
	****** SCHOOL INFO *******
	
	*DEMOGRAPHICS
	merge m:1 schoolcostcentercode endyear using "$DTAINT/school_demo.dta", gen(mergeDemEnr) keepusing(sch_female_enr sch_poverty_enr sch_college_grad_enr sch_gifted_enr sch_english_learner_enr sch_immigrant_enr sch_special_edu_enr sch_fst_yr_LAUSD_enr  sch_fst_yr_mag_info_enr  ) 
	drop if mergeDemEnr==2
	
	merge m:1 mag_cd_ app_endyear using "$DTAINT/school_demo.dta", gen(mergeDemApp) keepusing(sch_english_learner_app sch_immigrant_app sch_college_grad_app  sch_gifted_app sch_special_edu_app sch_female_app  sch_poverty_app sch_fst_yr_LAUSD_app sch_fst_yr_mag_info_app   )
	drop if mergeDemApp==2
	
	****VAM
	* Merge school vam for the school the st enrolled to and for school the st applied to
	* We don't have the vam and PeerQ estimates for 2014 (change of exams), 2020, and 2021 (covid) due to missing exams those years. Solution: for 2014, let' add the measures from the same school in 2013. For 2020 and 2021, assign 2019.  
	merge m:1 schoolcostcentercode endyear schooltype using "$DTAINT/vam_WITHspec_years.dta", gen(mergeVAM) keepusing(sch_vamhat_ela_hs_enr sch_ahat_ela_hs_enr sch_vamhat_math_hs_enr sch_ahat_math_hs_enr sch_vamhat_ela_ms_enr sch_ahat_ela_ms_enr sch_vamhat_math_ms_enr sch_ahat_math_ms_enr sch_vamhat_ela_ele_enr sch_ahat_ela_ele_enr sch_vamhat_math_ele_enr sch_ahat_math_ele_enr)
	drop if mergeVAM==2
	
	merge m:1 mag_cd_ app_endyear schooltype using "$DTAINT/vam_WITHspec_years.dta", gen(mergeVAM_app) keepusing(sch_vamhat_ela_hs_app sch_ahat_ela_hs_app sch_vamhat_math_hs_app sch_ahat_math_hs_app sch_vamhat_ela_ms_app sch_ahat_ela_ms_app sch_vamhat_math_ms_app sch_ahat_math_ms_app sch_vamhat_ela_ele_app sch_ahat_ela_ele_app sch_vamhat_math_ele_app sch_ahat_math_ele_app)
	drop if mergeVAM_app==2
	
	***PEERQ
	* Merge in school information, peer quality using mean of school's test score, for the enrolled and application school 
	merge m:1 endyear schoolcostcentercode using "$DTAINT/peerQ.dta", gen(mergePeerQ) keepusing(sch_z_ela_all_enr sch_z_math_all_enr)
	drop if mergePeerQ==2
	
	merge m:1 mag_cd_ app_endyear using "$DTAINT/peerQ.dta", gen(mergePeerQ_app) keepusing(sch_z_ela_all_app sch_z_math_all_app)
	drop if mergePeerQ_app==2
	
	***HHI
	*Merge school information, diversity measure (HHI) for the enrolled and application school 
	merge m:1 schoolcostcentercode endyear  using "$DTAINT/hhi",  gen(mergeHHI) keepusing(sch_hhi_enr sch_hispanic_enr sch_black_enr sch_white_enr sch_asian_enr sch_filipino_enr sch_other_enr)
	drop if mergeHHI==2
	
	merge m:1 mag_cd_ app_endyear using "$DTAINT/hhi",  gen(mergeHHI_app) keepusing(sch_hhi_app sch_hispanic_app sch_black_app sch_white_app sch_asian_app sch_filipino_app sch_other_app)
	drop if mergeHHI_app==2
	
	* Normalize so that higher HHI means more diversity
	replace sch_hhi_enr = 1-sch_hhi_enr
	replace sch_hhi_app= 1-sch_hhi_app
	
	***SES at school level for the enrolled and application school 
	merge m:1 schoolcostcentercode endyear  using "$DTAINT/ses_school.dta",  gen(mergeSES_enr) keepusing(sch_z_happy_index_enr sch_z_bully_index_enr sch_z_no_bully_index_enr sch_z_climate_index_enr sch_z_extra_index_enr sch_z_expectation_index_enr sch_z_perceptions_index_enr sch_z_interpersonal_index_enr sch_z_connectedness_index_enr sch_z_effort_index_enr sch_z_grit_index_enr)
	drop if mergeSES_enr==2
	
	merge m:1 mag_cd_ app_endyear using "$DTAINT/ses_school.dta",  gen(mergeSES_app) keepusing(sch_z_happy_index_app sch_z_bully_index_app sch_z_no_bully_index_app sch_z_climate_index_app sch_z_extra_index_app sch_z_expectation_index_app sch_z_perceptions_index_app sch_z_interpersonal_index_app sch_z_connectedness_index_app sch_z_effort_index_app sch_z_grit_index_app)
	drop if mergeSES_app==2
	

	***
	**Some variables cleaning and creation
	***
	
	format studentpseudoid %14.0g
	sort   studentpseudoid app_endyear, stable 
	
	*Identify the students who enroll in magnet.
	*Gen enrolled in applied program variable 
	gen enrolled_app = 0
	replace enrolled_app = 1 if mag_cd_== schoolcostcentercode
	replace enrolled_app = . if schoolcostcentercode==.
	
	*Gen standardize cutoff
	gen std_cutoff= sub_prio_pts - lottery_cutoff
	
	*Gen dummy variable for being above or below the cutoff
	gen above_cutoff = 0
	replace above_cutoff = 1 if sub_prio_pts>= lottery_cutoff
	replace above_cutoff = . if missing(sub_prio_pts)
	replace above_cutoff = . if missing(lottery_cutoff)
	
	* treatment indicator: above_cutoff
	label define treat 0 "Below" 1 "Equal or above" 
	label values  above_cutoff treat
 
	
	****Baseline covariates and outcomes 
	
	*Indicator of st who has info in LAUSD 
	gen inData = mergeDem==3 

	*Lags with three years previous to the application year (we don't include the application year because it is measured after the app)
	egen lag_ela = rowmean(L1ela L2ela L3ela)
	egen lag_math = rowmean(L1math L2math L3math)
	egen lag_suspensions = rowmean(L1numsuspensions L2numsuspensions L3numsuspensions)
	
	replace lag_suspensions = 0 if lag_suspensions==. 
	
	* create squared and cubed lag scores
	gen lag_math_2 = lag_math^2
	gen lag_math_3 = lag_math^3
	gen lag_ela_2 = lag_ela^2
	gen lag_ela_3 = lag_ela^3
	
	gen college_grad = 0
	replace college_grad=1 if parentedulevelname==4 | parentedulevelname==6
	
	gen special_edu = studentspedflag=="Y"
	
	gen gifted = !missing(studentgiftedprogramdescription) 
	
	gen poverty =  studentpovertyindicator=="Y"
	
	gen immigrant = born_usa==0
	
	gen female = gendercode==1
	
	gen english_learner=0
	replace english_learner = 1 if languageclasscode == "LEP"
	
	*For test scores. Imputing zeros and generating the mis variable for all  
	foreach var of varlist L1math L1ela lag_ela lag_math lag_math_2 lag_math_3 lag_ela_2 lag_ela_3   	  	{ 
		gen ms_`var' = (`var'==.)
		replace `var' = 0 if ms_`var'==1
		replace `var'=. if mergeDem!=3
		label var ms_`var' "Missing `var'"
	}
	
	gen anyScore_ms = (ms_lag_ela ==1) | (ms_lag_math ==1)
	label var anyScore_ms "Missing lag_ela or lag_math"
	
	
	*Putting missing when the student didn't merge to LAUSD and generating missing values dummy
	foreach var of varlist lag_suspensions college_grad special_edu gifted poverty immigrant female english_learner eng_home esp_home anyScore_ms { 
		replace `var'=. if mergeDem!=3
		gen ms_`var' = (`var'==.)
		label var ms_`var' "Missing `var'"
		
	}
	

	*These are not relevant variables, it is just to make the balance tables easier to run
	gen ms_inData = 0
	
	gen sch_pooled_vam_math_enr = 0
	replace sch_pooled_vam_math_enr = sch_vamhat_math_ele_enr if schooltype=="ele"
	replace sch_pooled_vam_math_enr = sch_vamhat_math_ms_enr if schooltype=="ms"
	replace sch_pooled_vam_math_enr = sch_vamhat_math_hs_enr if schooltype=="hs"
	replace sch_pooled_vam_math_enr =. if sch_vamhat_math_ele_enr ==. & schooltype=="ele"
	replace sch_pooled_vam_math_enr =. if sch_vamhat_math_ms_enr ==. & schooltype=="ms"
	replace sch_pooled_vam_math_enr =. if sch_vamhat_math_hs_enr ==. & schooltype=="hs"
	
	gen sch_pooled_vam_ela_enr=0
	replace sch_pooled_vam_ela_enr = sch_vamhat_ela_ele_enr if schooltype=="ele"
	replace sch_pooled_vam_ela_enr = sch_vamhat_ela_ms_enr if schooltype=="ms"
	replace sch_pooled_vam_ela_enr = sch_vamhat_ela_hs_enr if schooltype=="hs"
	replace sch_pooled_vam_ela_enr =. if sch_vamhat_ela_ele_enr==. & schooltype=="ele"
	replace sch_pooled_vam_ela_enr =. if sch_vamhat_ela_ms_enr==. & schooltype=="ms"
	replace sch_pooled_vam_ela_enr =. if sch_vamhat_ela_hs_enr==. & schooltype=="hs"
	
	*Measures for SES taking the average of three periods forward
	
	foreach var in "z_happy_index" "z_bully_index"  "z_no_bully_index" "z_climate_index"  "z_extra_index"  "z_expectation_index"  "z_perceptions_index"  "z_interpersonal_index"  "z_connectedness_index"  "z_effort_index"  "z_grit_index" {
		egen avg_`var' = rowmean(F1`var' F2`var' F3`var')
		label var avg_`var' "Average of `var' for 3 periods foward"
	}
	
	
	*Measures for high school testscore
	
	egen avg_Fela = rowmean(F1ela F2ela F3ela)
	egen avg_Fmath = rowmean(F1math F2math F3math)
	
	label var english_learner "Student is english learner"
	label var college_grad	"Parent reports going to college"
	label var special_edu	"Special education student"
	label var gifted	"Gifted student"
	label var poverty 	"Under poverty"
	label var immigrant	"Immigrant student"
	label var female	"Female"
	label var L1math  	"Baseline Math Score"
	label var L1ela 	"Baseline ELA Score"
	label var ms_anyScore	"Missing any ELA or MATH - avg 3 ys"
	label var esp_home	"Spanish at home"
	label var eng_home	"English at home"
	label var lag_ela 	"Baseline ELA Score - avg 3 ys"
	label var lag_math 	"Baseline Math Score - avg 3 ys"
	label var inData 	"Stays in District"
	label var lag_suspensions "Baseline Suspension Days"
	label var avg_Fela 	"Average of ELA for 3 periods foward"
	label var avg_Fmath 	"Average of Math for 3 periods foward"
	
	*Gen an indicator if the program is in a magnet campus
	gen magnet_campus = inlist(schoolcostcentercode, 1873801, 1855801, 1480801, 1519801, 1874101, 1891701, 1884201, 1789501, 1226901, 1760401, 1432201, 1274101, 1739001, 1588901, 1712301, 1875401, 1328801, 1627402, 1872701, 1230701, 1823002, 1885301, 1777701, 1756701, 1894301, 1225001, 1761501, 1250701, 1806601, 1284901, 1302701, 1311001, 1859601, 1330201, 1408201, 1412301, 1816801, 1445201, 1818901, 1493201, 1498601, 1521901, 1350001, 1605201, 1613701, 1775101, 1638401, 1649301, 1857701, 1694501, 1839601, 1839602, 1839603, 1765801, 1756201, 1811703, 1894310, 1767101, 1769901, 1331101, 1849001, 1782201, 1849301, 1661601, 1713701, 1253002, 1771501, 1357501, 1820801, 1501401, 1632901, 1687001, 1860601 )
	*Gen an indicator if the program is a highly gifted program
	gen highly_gifted = inlist(schoolcostcentercode,1350702, 1878602, 1810702, 1647902)
	
	
	* Merge in school-level traact info
	****
	/*preserve 
		import delimited "$DTARAW/censustract_2010_data.csv", clear stringcols(1)
		keep geoid total_population median_income  poverty_share college_educated_share black_share hispanic_share
		rename geoid censustract
		rename total_population schl_ctract_pop
		rename median_income schl_ctract_median_income
		rename college_educated_share schl_ctract_college_share
		rename black_share schl_ctract_black_share
		rename hispanic_share schl_ctract_hispanic_share
		rename poverty_share schl_ctract_poverty
		tempfile ctract_temp 
		save `ctract_temp'
	restore 
	*/
	* For the school the student is enrolled in
	preserve 
		import delimited "$ROOT/rawdata/costcenters_with_tractinfo.csv", clear
		rename geoid censustract
		keep schoolcostcentercode censustract magnetyesno latitude longitude
		rename (latitude longitude ) (schl_lat schl_lon)
		drop if schoolcostcentercode=="NA"
		destring schoolcostcentercode, replace
		tempfile schl_temp 
		save `schl_temp'
	restore 
	* For the school the student applied to 
		preserve 
		import delimited "$ROOT/rawdata/costcenters_with_tractinfo.csv", clear
		rename geoid censustract
		keep schoolcostcentercode censustract  latitude longitude
		rename (latitude longitude ) (app_lat app_lon)
		rename censustract app_censustrat 
		drop if schoolcostcentercode=="NA"
		destring schoolcostcentercode, replace
		rename schoolcostcentercode mag_cd_
		tempfile app_temp 
		save `app_temp'
	restore 
	preserve 
		use "$ROOT/rawdata/magnet_costcentercodes.dta", clear 
		rename mag_cd_ schoolcostcentercode
		drop if missing(schoolcostcentercode)
		tempfile mag_tmp
		save `mag_tmp'
	restore 
	preserve 
		import delimited "$ROOT/rawdata/all_ca_censusblock_grids.csv", delimiter("", collapse) varnames(1) emptylines(include) colrange(:3) stringcols(1) numericcols(2) clear
		drop if missing(geoid)
		rename geoid censusblockid 
		tempfile allblocks
		save `allblocks'
	restore
	* Merge in addresses at the time the student applies 
	preserve 
		use studentpseudoid  endyear censusblockid1 using "$BUILD/output/data/lausd2001_2023", clear
		rename censusblockid1 censusblockid 
		merge m:1 censusblockid using `allblocks', gen(mergeStuBlockCoords) keep(1 3)
		drop mergeStuBlockCoords 
		* make year endyear-1 so that we merge in addresses at the time students apply 
		replace endyear = endyear - 1
		rename (lat lon) (stu_lat stu_lon)
		destring studentpseudoid, replace
		tempfile allstublocks 
		save `allstublocks'
	restore

	* Merge in address the student is at the time of enrollment 
	preserve 
		use studentpseudoid  endyear censusblockid1 using "$BUILD/output/data/lausd2001_2023", clear
		rename censusblockid1 censusblockid 
		merge m:1 censusblockid using `allblocks', gen(mergeStuBlockCoords) keep(1 3)
		drop mergeStuBlockCoords 
		* make year endyear-1 so that we merge in addresses at the time students apply 
		replace endyear = endyear
		rename (lat lon) (stu_enroll_yr_lat stu_enroll_yr_lon)
		rename censusblockid censusblockid_enroll_yr 
		destring studentpseudoid, replace
		tempfile allstublocks2
		save `allstublocks2'
	restore	
	* Read in and prepare neighborhood school assignments 
	preserve 
		import delimited "$ROOT/rawdata/censusblock_ms_with_zones_and_schools.csv", stringcols(1) clear
		keep geoid school_id1 name1 lo_grd1 hi_grd1 cds1 school_id2 name2 lo_grd2 hi_grd2 cds2 school_id3 name3 lo_grd3 hi_grd3 cds3 school_id4 name4 lo_grd4 hi_grd4 cds4
		isid geoid 
		rename geoid censusblockid 
		tempfile schools
		save `schools'
	restore


	
	
	* Most of the not matching are those that leave the district (253 with costcentercodes without geo info)
	merge m:1  schoolcostcentercode using `schl_temp', gen(mergeSchoolGeo) keep(1 3)
	*merge m:1 censustract using `ctract_temp', gen(mergeSchoolCtractInfo) keep(1 3)
	merge m:1  schoolcostcentercode using `mag_tmp', gen(mergeMagFlag) keep(1 3)
	merge m:1 studentpseudoid endyear using `allstublocks', gen(mergeStuBlockCoords) keep(1 3)
	merge m:1 studentpseudoid endyear using `allstublocks2', gen(mergeStuBlockCoords2) keep(1 3)
	rename censustract schl_censustract
	rename magnetyesno magnetcampusflag 
	gen magnetprogramflag = mergeMagFlag==3 

	* Merge geo info for schools you applied to 
	merge m:1 mag_cd_ using `app_temp', gen(mergeAppGeo) keep(1 3)
	
	geodist stu_lat stu_lon schl_lat schl_lon, gen(schl_dist)
	geodist stu_lat stu_lon app_lat app_lon, gen(app_dist)
	geodist stu_enroll_yr_lat stu_enroll_yr_lon app_lat app_lon, gen(app_enroll_dist)

	label var schl_dist "Distance from residence to enrolled school"
	label var app_dist "Distance from residence to applied school"
	label var app_enroll_dist "Distance from residence (at the time of enrollment) to applied school (useful for this with missing addresses at the time of app)"
	label var stu_lat "Student residence latitude (time of applying)"
	label var stu_lon "Student residence longitude (time of applying)"
	label var schl_lat "School latitude (enrolled school)"
	label var schl_lon "School longitude (enrolled school)"
	label var app_lat "School latitude (applied school)"
	label var app_lon "School longitude (applied school)"
	label var stu_enroll_yr_lat "Student residence latitude (year of enrollment)"
	label var stu_enroll_yr_lon "Student residence longitude (year of enrollment)"

	label var magnetprogramflag "Flag if enrolled school is a magnet program"
	
	* merge neighborhood school info 
	merge m:1 censusblockid using `schools', gen(mergeNeighborhoodSchoolInfo) keep( 1 3)
	local schoolids "school_id1 school_id2 school_id3 school_id4"
	gen in_neighborhood_school = 0
	gen in_neighborhood_campus = 0
	foreach id of local schoolids{
		if "`id'" != "school_id1"{
			replace `id' = "" if `id'=="NA"
			destring `id' , replace
		}
		replace in_neighborhood_campus = 1 if `id' == preferredlocationcode & !missing(preferredlocationcode)
		replace in_neighborhood_school = 1 if `id' == schoollocationcode  & !missing(schoollocationcode)
	}
	* Note that if student lives outside LAUSD boundaries, we set them to not being enrolled in their neighborhood school/campus
	replace res_code = "" if res_code=="NA"
	destring res_code, replace
	gen in_neighborhood_campus_res_code = 0 
	replace in_neighborhood_campus_res_code = 1 if preferredlocationcode==res_code & !missing(preferredlocationcode)
	replace in_neighborhood_campus_res_code = 1 if schoollocationcode==res_code & !missing(schoollocationcode)
	label var in_neighborhood_campus "Indicator for being enrolled in neighborhood campus based on address"
	label var in_neighborhood_school "Indicator for being enrolled in neighborhood school based on address"
	label var in_neighborhood_campus_res_code "Indictor for being enrolled in neighborhood campus based on application data"

	drop school_id1 name1 lo_grd1 hi_grd1 cds1 school_id2 name2 lo_grd2 hi_grd2 cds2 school_id3 name3 lo_grd3 hi_grd3 cds3 school_id4 name4 lo_grd4 hi_grd4 cds4 mergeNeighborhoodSchoolInfo
	* Merge relative distances to neighborhood and magnet schools
	* Make sure the distance is at time of application 
	replace endyear = endyear - 1
	merge m:1 censusblockid schooltype endyear  using $DTAINT/magnet_distances.dta, keep(1 3) gen(mergeRelativeDistances)
	
	* temporarily rename censusblockid at time of app to at time of enrollment
	rename (censusblockid censusblockid_enroll_yr) (censusblockid_enroll_yr censusblockid) 
	rename 	(nearest_mag_rank nearest_schl_dist nearest_mag_dist avg_mag_dist20 rel_nearest_mag_dist) ///
		 	(app_nearest_mag_rank app_nearest_schl_dist app_nearest_mag_dist app_avg_mag_dist20 app_rel_nearest_mag_dist) 
	rename	(os_total_20 os2_total_20 os3_total_20 ) (app_os_total_20 app_os2_total_20 app_os3_total_20) 
	rename	(os_total_150 os2_total_150 os3_total_150 ) (app_os_total_150 app_os2_total_150 app_os3_total_150) 
	rename	(min_os_mag_dist min_os2_mag_dist min_os3_mag_dist) (app_min_os_mag_dist app_min_os2_mag_dist app_min_os3_mag_dist) 
	rename	(rel_os_mag_dist rel_os2_mag_dist rel_os3_mag_dist) (app_rel_os_mag_dist app_rel_os2_mag_dist app_rel_os3_mag_dist) 
	rename	(second_mag_dist rel_second_mag_dist ) (app_second_mag_dist app_rel_second_mag_dist) 
	merge m:1 censusblockid schooltype endyear using $DTAINT/magnet_distances.dta, keep(1 3) gen(mergeRelativeDistancesEnroll)
	* Change back to application year 
	replace endyear = endyear + 1 
	rename (nearest_mag_rank nearest_schl_dist nearest_mag_dist avg_mag_dist20 rel_nearest_mag_dist) ///
	       (enr_nearest_mag_rank enr_nearest_schl_dist enr_nearest_mag_dist enr_avg_mag_dist20 enr_rel_nearest_mag_dist) 
	rename	(os_total_20 os2_total_20 os3_total_20 ) (enr_os_total_20 enr_os2_total_20 enr_os3_total_20) 
	rename	(os_total_150 os2_total_150 os3_total_150 ) (enr_os_total_150 enr_os2_total_150 enr_os3_total_150) 
	rename	(min_os_mag_dist min_os2_mag_dist min_os3_mag_dist) (enr_min_os_mag_dist enr_min_os2_mag_dist enr_min_os3_mag_dist) 
	rename	(rel_os_mag_dist rel_os2_mag_dist rel_os3_mag_dist) (enr_rel_os_mag_dist enr_rel_os2_mag_dist enr_rel_os3_mag_dist) 
	rename	(second_mag_dist rel_second_mag_dist ) (enr_second_mag_dist enr_rel_second_mag_dist)
	* Rename back 
	rename (censusblockid censusblockid_enroll_yr) (censusblockid_enroll_yr censusblockid)


	* Add in college enrollment/grad outcomes and HS graduation
	preserve 
		use "$BUILD/output/data/intermediate/nsc_student_college.dta", clear 
		gen enroll_college = college_code_branch1!="NA" 
		gen enroll_two_year = two_year_four_year1=="2-year"
		gen enroll_four_year = two_year_four_year1=="4-year" 
		rename student_pseudo_id studentpseudoid
		keep studentpseudoid enroll_college enroll_two_year enroll_four_year 
		destring studentpseudoid, replace
		tempfile tempcoll 
		save `tempcoll'
	restore 
	
	merge m:1 studentpseudoid using `tempcoll', gen(mergeColege) keep( 1 3)
	label var enroll_college "Student enrolled in any college"
	label var enroll_two_year "Student enrolled in any two-year college"
	label var enroll_four_year "Student enrolled in any four-year college"



	* drop unnecessary variables 
	drop correc_imge_id validation_number formatted_timestamp ethnic_cd lottery_potential lottery_nocutoff

	* Some additional cleaning of variables and labeling 

	label var choice_ranking "Position of program on ROL (should be 1 for all)"
	label var school_yr "End year of application year"
	label var studentpseudoid "Student identifier"
	label var application_id "Application identifier"
	label var stu_grade "Grade the student is applying to"
	label var mag_cd_ "Magnet program code student is applying to"
	label var mag_name_ "Magnet program name student is applying to"
	label var mag_action_ "Assignment status"
	
	label var pwt_mgt_slctd_flag "Application: Magnet or pwt flag"
	label var vrfy_gft_flag "Application: Gifted flag"
	label var loc_cd_prsnt "Application: Current school location code"
	label var loc_name_prsnt "Application: Current school location name"
	label var sib_purls_id "Application: Sibling priority"
	label var stu_wait1_flag "Application: st applied in previous year priority"
	label var stu_wait2_flag "Application: st applied 2 years prior priority"
	label var stu_wait3_flag "Application: st applied 3 years prior priority"
	label var matriculation_flag "Application: st is in transition year between magnet programs"
	label var phbao_flag "Application: st assigned to a primarily Hisp, Black, or other non-white school"
	label var over_flag "Application: st assigned overcrowded school"
	label var priority_pts "Application: Calculated points by the district"
	label var sub_prio_pts "Application: Lottery number"
	label var res_code "Application: School code student is assigned if no action is taken"
	label var res_name "Application: School name student is assigned if no action is taken"
	label var application_date "Application: Date of application" 
	label var app_yr "Application: Application calendar year"
	label var app_mnth "Application: Application calendar month"
	label var app_endyear "Application: Endyear of application"
	
	rename (black white latino asian phbao) (app_stu_black app_stu_white app_stu_latino app_stu_asian app_stu_phbao)
	label var app_stu_black "Application: St is black"
	label var app_stu_latino "Application: St is latino"
	label var app_stu_asian "Application: St is asian"
	label var app_stu_white "Application: St is white"
	label var app_stu_phbao "Application: St is phbao"

	label var lottery_cutoff_rank "Lottery: Cutoff Rank"
	label var lottery_cutoff "Lottery: Estimated cutoff"
	label var lottery_coeff "Lottery: Estimated jump at cutoff"
	label var lottery_p_val "Lottery: P-value associated with estimated jump at cutoff"
	label var lottery_N_below "Lottery: Number of apps below estimated cutoff"
	label var lottery_N_above "Lottery: Number of apps above estimated cutoff"

	label var schooltype "School level (ele, ms, hs) st is applying to"

	drop studentclassofname  
	label var endyear "Demographics: Endyear of enrollment year"
	label var preferredlocationcode "Demographics: Preferred location code"
	label var preferredlocationname "Demographics: Preferred location name"
	label var schoollocationcode "Demographics: School location code"
	label var schoollocationname "Demographics: School location name"
	label var schoolcostcentercode "Demographics: School cost center code"
	label var gradecode "Demographics: Grade code (enrolled year)"
	label var gendercode "Demographics: Gender"
	label var studentbirthcountry "Demographics: St birth country"
	label var parentedulevelname "Demographics: Parent education level"
	label var studentusschoolfirstattenddate "Demographics: Date st first attended US School"
	label var languageclasscode "Demographics: English Learner class code"
	label var studentspedflag "Demographics: St sped flag"
	label var studentgiftedprogramdescription "Demographics: St gifted program flag"
	label var studentpovertyindicator "Demographics: St poverty indicator"
	label var eng_home "Demographics: St speaks english at home"
	label var esp_home "Demographics: St speaks spanish at home"
	label var born_usa "Demographics: St born in USA"

	label var ytdofattendeddays "Attendance: Attended Days"
	label var ytdofenrolleddays "Attendance: Enrolled days"
	label var ytdofabsentdays "Attendance: Absent days"
	drop otherabsencecount appabsencecount  
	label var ytdofattendance "Attendance: Attendance rate"

	label var enrolled_app "Lottery: Enrolled in top-ranked program"
	label var std_cutoff "Lottery: Standardized cutoff"
	label var above_cutoff "Lottery: Random number above cutoff"
	order enrolled_app std_cutoff above_cutoff, after(lottery_ID)

	save "$DTAFINAL/magnet_lotteries_1perc.dta", replace
	
	log close
	
	
end



********************************************************************************
* merge_other_choices
* 
* Description: merges 2nd and 3rd choices options to generate another datasets 
*that includes the rest of the options and the result of the 2nd and 3rd lotteries

*The year of merging is the application year 
********************************************************************************/
capture program drop merge_other_choices 
program define merge_other_choices 
	cap log close
	log using "$LOGS/merge_other_choices.log", text replace

**Create second database with the information from the second and third choice
	
	use "$DTAFINAL/magnet_lotteries_1perc.dta", clear
	
	gen has_1st=1
	*Bring in information if the applicant selected a second and third choice
	
	merge 1:1 studentpseudoid app_endyear using "$BUILD/output/data/magnet_cutoffs/raw_all_app_choice2Info.dta", gen(mergeChoice2) keepusing(has_2nd offer_2nd) 
	drop if mergeChoice2==2
	
	merge 1:1 studentpseudoid app_endyear using "$BUILD/output/data/magnet_cutoffs/raw_all_app_choice3Info.dta", gen(mergeChoice3) keepusing(has_3rd offer_3rd) 
	drop if mergeChoice3==2
	
	*Bring in information if the applicant participated in a second and third choice lottery
	
	merge 1:1 studentpseudoid app_endyear using "$BUILD/output/data/magnet_cutoffs/potential_lotteries_WITHcutoffs_choice2.dta", gen(mergeL2) keepusing(mag_cd_2nd mag_name_2nd offer_magnet_2nd lottery_ID2nd lottery_cutoff_rank2nd lottery_cutoff2nd above_cutoff_2nd lottery_nocutoff2nd)       
	drop if mergeL2==2
	drop mergeL2
	
	*Identify the students who enroll in magnet option 2.
	*Gen enrolled in applied program variable 
	gen enrolled_app_2nd = 0
	replace enrolled_app_2nd = 1 if mag_cd_2nd== schoolcostcentercode
	replace enrolled_app_2nd = . if schoolcostcentercode==.	
	
	merge 1:1 studentpseudoid app_endyear using "$BUILD/output/data/magnet_cutoffs/potential_lotteries_WITHcutoffs_choice3.dta", gen(mergeL3) keepusing(mag_cd_3rd mag_name_3rd offer_magnet_3rd lottery_ID3rd lottery_cutoff_rank3rd lottery_cutoff3rd above_cutoff_3rd lottery_nocutoff3rd)
	drop if mergeL3==2
	drop mergeL3

	save "$DTAFINAL/magnet_lotteries_1perc_allchoices.dta", replace
	

	
end

capture program drop prep_descriptives
program define prep_descriptives 

	* Determine magnet, charter, and dle codes 
	use "$BUILD/output/data/magnet_rank.dta", clear
	keep if stu_grade==6 
    keep if choice_ranking ==1
    contract mag_cd_ mag_name
	drop if mag_cd_==. | mag_cd_==0
	*There is an issue with the magnet codes, so we collapse at the costcentercode (mag_cd_) and drop programs that are applied to less than 15 times across all years
	drop if _freq <15
	drop _freq
	contract mag_cd_
	rename mag_cd_ magnet_code 
	rename magnet_code code 
	gen magnetenroll=1 
	tempfile tempcodes 
	save `tempcodes'
	
	use $BUILD/output/data/afc_app.dta, clear 
	keep if grade_applied=="06"
	contract cost_center_code 
	replace cost_center_code = subinstr(cost_center_code, "TWS", "", .) 
	destring cost_center_code, replace
	drop if cost_center_code==. | cost_center_code==0
	drop _freq
	rename cost_center_code afc_code 
	rename afc_code code
	gen afcenroll=1
	* four duplictes 

	append using `tempcodes'
	duplicates drop code , force 
	save `tempcodes', replace 


	*******
	* Create datasets that tracks who gets offers 
	*******
	/*
	import delimited "$UE/Deidentified_NEW_MAG_ALL.csv", clear stringcols(_all)
	destring application_id school_yr, replace
	merge 1:1 application_id using "$BUILDROOT/output/data/apps_time_stamp.dta"
	drop _merge
	replace app_endyear = school_yr if missing(app_endyear) // so we don;t lose the 2023, 2024 people
	sort student_pseudo_id
	rename (student_pseudo_id ) (studentpseudoid) 
	gen magnet_offer = mag_action_1=="B" | mag_action_2 =="B" | mag_action_3=="B"
	contract studentpseudoid app_endyear magnet_offer
	gen endyear = app_endyear //just to use for the match in the build code 
	drop _freq 
	duplicates drop studentpseudoid endyear, force 
	save "$DTAINT/magnet_offers.dta", replace 
	*/
	
	
	*******
	*MAGNET Application data available from years 2000 to 2024   
	*******
	
	import delimited "$BUILD/rawdata/Campos_Master_DUA/Unified_Enrollment/Deidentified_NEW_MAG_ALL.csv", clear stringcols(_all)

	destring application_id school_yr, replace
	merge 1:1 application_id using "$BUILD/output/data/apps_time_stamp.dta"
	drop if _merge==2
	drop _merge
	*To avoid losing the 2023,2024 people:
	replace app_endyear = school_yr if missing(app_endyear)
	sort student_pseudo_id
	rename (student_pseudo_id ) (studentpseudoid) 
	*We don;t have offer info for 2024
	gen magnet_offer = mag_action_1=="B" | mag_action_2 =="B" | mag_action_3=="B"

	/*Checks for duplicates at the st/yr level (people who applied multiple times in the same year)
	sort studentpseudoid app_endyear
    quietly by studentpseudoid app_endyear:  gen dup = cond(_N==1,0,_n)
	
	*We some cases of 2 or more applications ~7,000. Careful for the dropping of dups
	drop dup
	*/
	
	*Create an application variable dataset
	contract studentpseudoid app_endyear magnet_offer
	*To use for the merge in the build code:
	gen endyear = app_endyear 
	bysort studentpseudoid app_endyear: egen magnet_offer2 = max(magnet_offer)
	replace magnet_offer = magnet_offer2
	drop _freq magnet_offer2 
	duplicates drop studentpseudoid endyear, force 
	save "$DTARAW/magnet_offers.dta", replace 
	
	
	
	*******
	* DLE
	*******
	
	import delimited  "$BUILD/rawdata/Campos_Master_DUA/Unified_Enrollment/Deidentified_NEW_DLE_ALL.csv", clear delimiter(comma) bindquote(strict)  stripquote(no) stringcols(_all)
	
	destring application_id school_yr, replace
	merge 1:1 application_id using "$BUILD/output/data/apps_time_stamp.dta"
	drop if _merge==2
	drop _merge
	*To avoid losing the 2023,2024 people:
	replace app_endyear = school_yr if missing(app_endyear) 
	sort student_pseudo_id
	rename (student_pseudo_id ) (studentpseudoid)
	keep if app_late_flag=="FALSE"
	*We don;t have offer info for 2024
	gen dle_offer = first_choice_action=="B" | second_choice_action =="B" | third_choice_action=="B"
	
	/*Checks for duplicates at the st/yr level (people who applied multiple times in the same year)
	sort studentpseudoid app_endyear
    quietly by studentpseudoid app_endyear:  gen dup = cond(_N==1,0,_n)
	
	*We some cases of 2 or more applications ~9,000. Careful for the dropping of dups
	drop dup
	*/
	
	*Create an application variable dataset
	contract studentpseudoid app_endyear dle_offer
	*To use for the merge in the build code:
	gen endyear = app_endyear 
	bysort studentpseudoid app_endyear: egen dle_offer2 = max(dle_offer)
	replace dle_offer = dle_offer2
	drop _freq dle_offer2 
	duplicates drop studentpseudoid endyear, force 
	save "$DTARAW/dle_offers.dta", replace 
	

	
	*******
	*SAS * Application data available from years 2019 to 2024   
	*******
	
	import delimited "$BUILD/rawdata/Campos_Master_DUA/Unified_Enrollment/Deidentified_NEW_SAS_ALL.csv", clear bindquote(strict)  stripquote(no) stringcols(_all)
	destring application_id school_yr, replace
	merge 1:1 application_id using "$BUILD/output/data/apps_time_stamp.dta"
	drop if _merge==2
	drop _merge
	*To avoid losing the 2023,2024 people:
	replace app_endyear = school_yr if missing(app_endyear) 
	sort student_pseudo_id
	rename (student_pseudo_id ) (studentpseudoid)
	*We don't have offer info for 2024
	gen sas_offer = sas_action_1=="B" | sas_action_2 =="B" | sas_action_3=="B" 
	
	/*Checks for duplicates at the st/yr level (people who applied multiple times in the same year)
	sort studentpseudoid app_endyear
    quietly by studentpseudoid app_endyear:  gen dup = cond(_N==1,0,_n)
	
	*We have only 18 cases of 2 or more applications. Careful for the dropping of dups
	drop dup
	*/
	
	*Create an application variable dataset
	contract studentpseudoid app_endyear sas_offer
	*To use for the merge in the build code:
	gen endyear = app_endyear
	bysort studentpseudoid app_endyear: egen sas_offer2 = max(sas_offer)
	replace sas_offer = sas_offer2
	drop _freq sas_offer2 
	duplicates drop studentpseudoid endyear, force 
	save "$DTARAW/sas_offers.dta", replace 
	
	*******
	*ACS  ********* these are 10 programs and stable acrross time 
	*******
	* Application data available from years 2019 to 2024   
	import delimited "$BUILD/rawdata/Campos_Master_DUA/Unified_Enrollment/Deidentified_NEW_ACS_ALL.csv", clear stringcols(_all)
	destring application_id school_yr, replace
	merge 1:1 application_id using "$BUILD/output/data/apps_time_stamp.dta" 
	drop if _merge==2
	drop _merge
	*To avoid losing the 2023,2024 people:
	replace app_endyear = school_yr if missing(app_endyear) 
	sort student_pseudo_id
	rename (student_pseudo_id ) (studentpseudoid) 
	*We don;t have offer info for 2024
	gen acs_offer = acs_action_1=="B" 
	/*Checks for duplicates at the st/yr level (people who applied multiple times in the same year)
	sort studentpseudoid app_endyear
    quietly by studentpseudoid app_endyear:  gen dup = cond(_N==1,0,_n)
	
	*We have only 16 cases of 2 or more applications. Careful for the dropping of dups
	drop dup
	*/
	
	*Create an application variable dataset
	contract studentpseudoid app_endyear acs_offer
	*To use for the merge in the build code:
	gen endyear = app_endyear 
	bysort studentpseudoid app_endyear: egen acs_offer2 = max(acs_offer)
	replace acs_offer = acs_offer2
	drop _freq acs_offer2 
	duplicates drop studentpseudoid endyear, force 
	save "$DTARAW/acs_offers.dta", replace 
	
	*******
	*AFC  
	*******
	* Application data available from years 2020 to 2024   
	* When they apply to charter they don't put a ranking. It is a lottery. The lottery is public for people who attend the lottery. Relevant variables: school_decision (B accepted, W waitlisted, P pending -probably a late application), parent_decision_description
	* They can apply to more than one school. That is why we have 70,556 obs but 33,469 unique at the st year level
	import delimited "$BUILD/rawdata/Campos_Master_DUA/Unified_Enrollment/Deidentified_NEW_AFC_ALL.csv", clear stringcols(_all)
	destring application_id school_yr, replace
	*application_id is not unique for this program
	merge m:1 application_id using "$BUILD/output/data/apps_time_stamp.dta"  
	drop if _merge==2
	drop _merge
	*To avoid losing the 2023,2024 people:
	replace app_endyear = school_yr if missing(app_endyear)
	sort student_pseudo_id
	rename (student_pseudo_id ) (studentpseudoid) 
	*We don;t have offer info for 2024
	gen afc_offer = school_decision=="B" 
	
	/*Checks for duplicates at the st/yr level (people who applied to multiple affiliated charters)
	sort studentpseudoid app_endyear 
    quietly by studentpseudoid app_endyear:  gen dup = cond(_N==1,0,_n)
	*We have only 20,495 unique applications
	
	
	*Check for dups at the st/yr charter level. There are only 3 
	sort studentpseudoid app_endyear cost_center_code
	 quietly by studentpseudoid app_endyear cost_center_code:  gen dup2 = cond(_N==1,0,_n)
	 drop dup dup2
	 
	*From the data it seems students can get more than one offer from diff schools. 
	*Some examples: students 154112971241, 829280413406, 819611989996. 
	*/
	
	*Create an application variable dataset. Keeping in mind the feature of multiple offers 
	contract studentpseudoid app_endyear afc_offer
	*To use for the merge in the build code:
	gen endyear = app_endyear
	bysort studentpseudoid app_endyear: egen afc_offer2 = max(afc_offer)
	replace afc_offer = afc_offer2
	drop _freq afc_offer2 
	duplicates drop studentpseudoid endyear, force 
	save "$DTARAW/afc_offers.dta", replace 
	
	

	* Starting with LA panel add info about who applies to the various choice programs 
	use if endyear >=2004 & endyear <=2020  using $BUILD/output/data/lausd2001_2023_encoded.dta, clear 
	merge m:1 studentpseudoid endyear using $DTARAW/magnet_offers, gen(mergeOffers) keep( 1 3)
	drop mergeOffers 
	gen appliedMagnet = mergeMagnet ==3
	merge m:1 studentpseudoid endyear using $DTARAW/dle_offers, gen(mergeDLEoffers) keep(1 3)
	drop if mergeDLEoffers==2 
	drop mergeDLEoffers 
	gen appliedDLE = mergeDLE ==3
	merge m:1 studentpseudoid endyear using $DTARAW/afc_offers, gen(mergeAFCoffers) keep(1 3)
	drop if mergeAFCoffers==2 
	drop mergeAFCoffers 
	gen appliedAFC = mergeAFC ==3
	replace magnet_offer = 0 if missing(magnet_offer)
	replace dle_offer = 0 if missing(dle_offer)
	replace afc_offer = 0 if missing(afc_offer)
	gen appliedAny = (appliedMagnet==1 | appliedDLE==1 | appliedAFC==1)
	gen anyOffer = (magnet_offer==1 | dle_offer==1 | afc_offer==1)
	label var appliedAny "Student applied to any choice program"
	label var anyOffer "Student received an offer to any choice program"
	egen sid = group(studentpseudoid)
	xtset sid endyear 
	decode schoollocationname, gen(sname) 
	gen next_schoolname = F1.schoollocationname
	gen next_scode = F1.schoolcostcentercode 
	label values next_schoolname schoollocationname_
	decode next_schoolname, gen(next_sname)
	label var next_schoolname "Name of school student enrolls in in sixth grade"
	label var next_scode "Cost center code of school student enrolls in in sixth grade"
	label var next_sname "Name of school student enrolls in in sixth grade"

	
	* Add information about current school been a magnet
	gen code = schoolcostcentercode
	merge m:1 code using `tempcodes', gen(merge_enroll_choice) keep (1 3)
	gen cur_enroll_choice = merge_enroll_choice==3
	replace cur_enroll_choice = 0 if afcenroll==1 & endyear <2020
	drop code merge_enroll_choice
	label var cur_enroll_choice "Indicator for being enrolled in a choice program at the time of application"

	* Add information about next year school been a magnet
	gen code = next_scode
	merge m:1 code using `tempcodes', gen(merge_next_sch_choice) keep(1 3)
	drop if merge_next_sch_choice==2
	gen next_schoice= merge_next_sch_choice==3 
	replace next_schoice = 0 if afcenroll==1 & endyear <2020
	* Address issues with conversion charters -- they keep the same code so being flagged as if students are enrolling in choice programs but not really 
	replace next_schoice = 0 if !(regexm(next_sname , "MAG") | regexm(next_sname, "MG") | regexm(next_sname, "CS") | regexm(next_sname, "CMS") | regexm(next_sname, "CM"))
	replace next_schoice = 0 if  ( (regexm(next_sname, "CS") | regexm(next_sname, "CMS") | regexm(next_sname, "CM") ) & endyear <2020)
	* Charters are only there for 2021 
	drop code merge_next_sch_choice
	label var next_schoice "Indicator for being enrolled in a choice program at the time of sixth grade enrollment"

	unique next_scode if next_schoice==1 // this identifies the number of unique choice programs for this sample

	* Create some demographic flags
	gen hispanic = ethnicitydescription==1
	gen black = ethnicitydescription==2
	gen white = ethnicitydescription==3
	gen asian= ethnicitydescription==4 
	gen other = hispanic==0 & black==0 & white==0 & asian==0 
	gen female = gendercode==1
	gen poverty = studentpovertyindicator==3
	gen el = languageclasscode==3
	gen gifted = !missing(studentgiftedprogramdescription) 

	* Weird issue where elementary school students are more likely to be missing addresses
	* Fix by assigning fifth grade students their sixth-grade address
	gen sixthGradeAddress = censusblockid1 if gradecode ==6 
	bys studentpseudoid : egen addressSix = mean(sixthGradeAddress )
	replace censusblockid1 = addressSix if missing(censusblockid1 ) & gradecode ==5
		
	* Keep cohorts of fifth graders 
	keep if gradecode==5
	destring studentpseudoid, replace 
	* Create school type variable (for future merges)
	gen schooltype = "ms" if gradecode ==5

	* Merge distance data 
	rename censusblockid1 censusblockid
	decode censusblockid, gen(blockid)
	rename (censusblockid blockid) (blockid censusblockid)

	merge m:1 censusblockid schooltype endyear using $DTAINT/choice_distances.dta   , gen(mergeDistances) keep(1 3) 
	keep if mergeDistances==3 
	replace magnet_offer = 0 if missing(magnet_offer)
	* Calculate quintiles of relative distance to nearest choice school 
	egen dist_quint = xtile(rel_nearest_choice_dist), n(5) by(schooltype endyear)

	save $DTAFINAL/5th_grade_cohorts.dta, replace 
end 

capture program drop prep_lottery 
program define prep_lottery 

	use $DTAFINAL/lausd_5thgrade_cohorts_2004_2017.dta, clear 
	destring studentpseudoid, replace
	tempfile temp 
	save `temp'

	use if schooltype=="ms" using "$DTAFINAL/magnet_lotteries_1perc.dta", clear
	/*
	levelsof mag_cd_, local(lotteryIDs)
	foreach id of local lotteryIDs {
			cap reghdfe enrolled_app above_cutoff if mag_cd_ == `id', noabsorb
			if _rc==0{
			if _b[above_cutoff] < 0  {
				drop if mag_cd_==`id'
			}
			}
	} 
	*/

	* Merge in socioemotional outcomes 
	replace endyear = endyear -1  
	merge m:1 studentpseudoid endyear using `temp', keepusing(z_suspensions z_absent z_gpa z_socio F1suspensions F2suspensions F3suspensions F1absent F2absent F3absent F1gpa F2gpa F3gpa avg_Fmath avg_Fela avg_Fsuspensions avg_Fabsent avg_Fgpa ytdofabsentdays numSuspensions) keep(1 3)

	
	* Before removing negative FS, we had 1908 lottery IDs; after we have 1574 
	unique lottery_ID 


	rename mag_cd_ mag_cd
	rename mag_name mag_nam

	* additional variables
	gen hispanic = (ethnicity==1)
	gen black = (ethnicity==2)
	gen white = (ethnicity==3)
	gen other = (app_stu_black==0 & app_stu_white==0 & hispanic==0)
	
	
	* balance table labels
    replace L1numsuspensions = 0 if missing(L1numsuspensions)
	local chars "female english_learner special_edu poverty college_grad eng_home esp_home L1numsuspensions L1math L1ela app_nearest_mag_dist anyScore_ms "
    gen ms_app_nearest_mag_dist = missing(app_nearest_mag_dist)
    gen ms_ytdofattendance = missing(ytdofattendance)
    gen ms_L1numsuspensions = missing(L1numsuspensions)


	* Neighborhood school and campuse indicators 
	rename in_neighborhood_school in_nghd_sch
	rename in_neighborhood_campus in_nghd_camp

	gen app_dist_rel = app_dist - app_nearest_schl_dist
	
	* Add some jitter so that quintiles are balanced (mass at zero messes things up)
	gen dist_rel = app_nearest_mag_dist - app_nearest_schl_dist+ runiform(0,0.25)

	egen app_dist_rel_quint=xtile(app_dist_rel), n(5) by(schooltype endyear)
	egen dist_rel_quint=xtile(dist_rel), n(5) by(schooltype endyear)
	

	* Interact offers with quintiles of distance to applied school
	gen above_dist1 = above_cutoff*(app_dist_rel_quint==1)
	gen above_dist2 = above_cutoff*(app_dist_rel_quint==2)
	gen above_dist3 = above_cutoff*(app_dist_rel_quint==3)
	gen above_dist4 = above_cutoff*(app_dist_rel_quint==4)
	gen above_dist5 = above_cutoff*(app_dist_rel_quint==5)
	gen enroll_dist1 = enrolled_app*(app_dist_rel_quint==1)
	gen enroll_dist2 = enrolled_app*(app_dist_rel_quint==2)
	gen enroll_dist3 = enrolled_app*(app_dist_rel_quint==3)
	gen enroll_dist4 = enrolled_app*(app_dist_rel_quint==4)
	gen enroll_dist5 = enrolled_app*(app_dist_rel_quint==5)

	* Interact offers with quintiles of relative distance to nearest school 
	gen above_rel_dist1 = above_cutoff*(dist_rel_quint==1)
	gen above_rel_dist2 = above_cutoff*(dist_rel_quint==2)
	gen above_rel_dist3 = above_cutoff*(dist_rel_quint==3)
	gen above_rel_dist4 = above_cutoff*(dist_rel_quint==4)
	gen above_rel_dist5 = above_cutoff*(dist_rel_quint==5)
	gen enroll_rel_dist1 = enrolled_app*(dist_rel_quint==1)
	gen enroll_rel_dist2 = enrolled_app*(dist_rel_quint==2)
	gen enroll_rel_dist3 = enrolled_app*(dist_rel_quint==3)
	gen enroll_rel_dist4 = enrolled_app*(dist_rel_quint==4)	
	gen enroll_rel_dist5 = enrolled_app*(dist_rel_quint==5)

	save $DTAFINAL/lotteries.dta, replace 



end 


capture program drop prep_event_study_tract 
program define prep_event_study_tract
syntax, [pscore(string) relative(string)]

	ssc install getcensus
    *************************************************************************
    * Download Census tract data for matching exercise 
    *************************************************************************
    local years 2022               
    local state "06"               
    local counties "037"     

    * 1) SUBJECT TABLES (ST): poverty rate, school-aged children
    *    S1701_C03_001 = % below poverty level (all people)
    *    S0101_C01_020 = pop. count 5–14; S0101_C01_021 = count 15–17
    getcensus S1701_C03_001  S0101_C01_020 S0101_C01_021 ///
        , years(`years') sample(5) ///
        geography(tract) statefips(`state') countyfips(`counties') ///
        clear

    rename s1701_c03_001e pov_rate
    rename s0101_c01_020e n_5_14
    rename s0101_c01_021e n_15_17
    gen n_school_age = n_5_14 + n_15_17
    tempfile st
    save `st', replace

    * 2) DETAILED TABLE (DT): median household income
    *    B19013_001 = Median household income (dollars)
    getcensus B19013_001 ///
        , years(`years') sample(5) ///
        geography(tract) statefips(`state') countyfips(`counties') ///
        clear

    rename b19013_001e medhhinc
    tempfile dt
    save `dt', replace

    * 3) DATA PROFILE (DP): race/ethnicity shares (percent)
    *    B03002_003 = White (alone or in combination)
    *    B03002_004 = Black or African American (alone or in combination)
    *    B03002_012 = Hispanic or Latino (of any race)
    getcensus B03002_003 B03002_004 B03002_012 ///
        , years(`years') sample(5) ///
        geography(tract) statefips(`state') countyfips(`counties') ///
        clear

    rename b03002_003e white_pct
    rename b03002_004e black_pct
    rename b03002_012e hispanic_pct
    tempfile dp
    save `dp', replace

    * Merge all three pulls by tract and keep just what we need
    use `st', clear
    merge 1:1 state county tract using `dt', nogen
    merge 1:1 state county tract using `dp', nogen
    keep state county tract name   ///
        pov_rate  n_school_age medhhinc ///
        white_pct black_pct hispanic_pct
    order state county tract name pov_rate medhhinc n_school_age ///
         white_pct black_pct hispanic_pct
    tempfile tracts 
    save `tracts', replace 
    gen censustract = state + county + tract 

    save $DTAINT/census_tracts.dta, replace 
  
    clear 
    tempfile stacks 
    save `stacks', emptyok

    *************************************************************************
    * Create a block-year level dataset with mag apps and other vars 
    *************************************************************************
    use "$BUILD/output/data/lausd2001_2023_encoded.dta", clear
    gen appliedMag = mergeMagnet ==3  // Applied to any magnet 
    keep if gradecode==4 | gradecode ==5 | gradecode ==6 | gradecode ==7 // students applying to middle schools 
	egen sid = group(studentpseudoid)
	xtset sid endyear 
	gen lag_math = L1.z_math_all 
	gen lag_ela = L1.z_ela_all
    gen stu=1
	drop if gradecode==4
    decode censusblockid1 , gen(blockid) 
    drop if missing(blockid)
	gen poverty = studentpovertyindicator==3
	reg z_ela_all i.gendercode i.ethnicitydescription eng_home esp_home i.gradecode i.endyear ///
		total_num_suspensions poverty z_math_all if endyear<=2013, vce(robust)  
	predict yhat, xb 
	predict yres, res 
    gen app_math = z_math_all if appliedMag ==1 // math scores of applicants
    gen app_ela = z_ela_all if appliedMag ==1 // ela scores of applicants 

	gen app_black = (ethnicitydescription == 2) if appliedMag ==1
	gen app_hispanic = (ethnicitydescription == 1) if appliedMag ==1
	gen app_white = (ethnicitydescription == 3) if appliedMag ==1
	gen app_eng_home = (eng_home==1) if appliedMag ==1
	gen app_female = (gendercode==1) if appliedMag ==1
	gen app_suspensions = (total_num_suspensions) if appliedMag ==1
	gen app_poverty = poverty if appliedMag==1
	gen app_yhat = (yhat) if appliedMag ==1
	gen app_yres = yres if appliedMag==1

    collapse (sum) appliedMag (count) num_stu = stu ///
		(mean) z_math_all z_ela_all app_math app_ela app_yhat app_yres ///
			app_black app_hispanic app_white app_eng_home app_female app_suspensions app_poverty ///
		, by(blockid endyear) 
    destring blockid, replace 
    format blockid %20.0g
    tempfile temp 
    save `temp'

    * Keep a list of relevant census block ids 
    contract blockid 
    drop _freq 
    tempfile temp2 
    save `temp2'

    *************************************************************************************
    * Create census block centroid data (for distance calculations later)
    *************************************************************************************
	cd $DTARAW
    spshape2dta "tl_2020_06_tabblock20.shp", replace saving(blocks_ca)
    * Attribute file with GEOID and centroids
	    use blocks_ca.dta, clear
	       keep if COUNTYFP20 =="037"
    keep GEOID20 _CX _CY
    rename (GEOID20 _CX _CY) (blockid lon lat)    // lon/lat in decimal degrees (NAD83/WGS84)
    tostring blockid, replace force  // if needed
    tempfile blk_centroids
    save `blk_centroids'
	cd $ROOT

    *************************************************************************************
    * Read in block-to-school distance data (created in the LAUSD master build directory)
    *************************************************************************************
    use "$BUILD/rawdata/aux_data/relative_choice_distances_grade5.dta", clear
    rename censusblockid blockid 
    merge m:1 blockid using `blk_centroids', keep(3) nogen
    assert !missing(lat, lon)
    rename blockid censusblockid
    destring censusblockid , gen(blockid)
    format blockid %20.0g
    drop if missing(blockid)
    xtset blockid endyear

    * Create a variable that tracks changes in relative distance to a block's nearest choice school (used to define events)
	if "`relative'" == "yes" {
	    gen change_rel = rel_choice_dist  - L1.rel_choice_dist 
	}
	else {
		gen change_rel = min_choice_dist  - L1.min_choice_dist
	}
    * Those that do not merge (==1) are blocks in LA county that LAUSD students do not live in
    merge m:1 blockid using `temp2', gen(mergeExistBlock) keep(3) // Keep relevant census block IDs
    merge 1:1 blockid endyear using `temp', gen(mergeApps) keep(1 3) // Merge block application data 
    replace appliedMag = 0 if missing(appliedMag)
    replace num_stu = 0 if missing(num_stu)

    *************************************************************************************
    * Collapse to tract level and create additional vars to find events 
    * Note that for distance variables, we are creating tract-level averages 
    *************************************************************************************
    gen censustract = substr(censusblockid, 1,11)
    bys censustract endyear : egen total_stu = total(num_stu)
    gen stu_share = num_stu /total_stu 

    * Weight test scores by block shares before collapsing 
    gen z_math_wt = z_math_all * stu_share 
    gen z_ela_wt = z_ela_all * stu_share
    gen app_math_wt = app_math * stu_share 
    gen app_ela_wt = app_ela * stu_share 
	gen app_black_wt = app_black * stu_share
	gen app_hispanic_wt = app_hispanic * stu_share
	gen app_white_wt = app_white * stu_share
	gen app_eng_home_wt = app_eng_home * stu_share	
	gen app_suspensions_wt = app_suspensions * stu_share
	gen app_yhat_wt = app_yhat * stu_share
	gen app_poverty_wt = app_poverty * stu_share
	gen app_yres_wt = app_yres * stu_share

    encode censustract, gen(censustractid)
    collapse (mean) dist min_distance min_choice_dist choice_dist rel_choice_distance change_rel lat lon ///
            (sum) z_math_wt z_ela_wt app_math_wt app_ela_wt ///
				app_black_wt app_hispanic_wt app_white_wt app_eng_home_wt app_suspensions_wt ///
				app_yhat_wt app_poverty_wt app_yres_wt ///
                num_stu appliedMag, ///
            by(censustractid censustract endyear )

    * Merge in census tract information
    merge m:1 censustract using $DTAINT/census_tracts.dta, gen(mergeTractInfo) keep(1 3) 
    xtset censustractid endyear

    * Create flags for blocks not experiencing a change in distance to nearest choice school 
    gen L1null = L2.min_choice_dist  == L1.min_choice_dist  
    gen L2null = L3.min_choice_dist  == L2.min_choice_dist  
    gen L3null = L4.min_choice_dist  == L3.min_choice_dist 
    
    * Create variables that keep track of tract-level enrolled LAUSD students 
    gen L1num_stu = L1.num_stu
    gen L2num_stu = L2.num_stu
    gen L3num_stu = L3.num_stu

    * Drop tracts that have low student enrollment 
    bys censustractid: egen maxStu = max(num_stu)
    drop if maxStu <10 

    * Identify potential events 
    sum change_rel if change_rel <0 , det 
    gen event_occur = change_rel < -0.5 & endyear !=2004 if !missing(change_rel)

    *************************************************************************************
    * Identify candidate event years
    * Restrictions: 
    * 1) Avg. tract-level reduction in relative distance is greater than 1 mile (event_occur==1)
    * 2) No changes in proximity to choice schools in t-1, t-2, and t-3 (L1null==1 & L2null==1)
    * 3) The tract has enrolled students at t, t-1, t-2, and t-3 (L1num_stu>0 & L2num_stu>0)
    * 4) Restrict to years that we use in the structural estimation (+ 2)
    *************************************************************************************
    levelsof endyear if event_occur==1 ///
        & L1null==1 & L2null==1 & num_stu>0 & L1num_stu>0 & L2num_stu>0 ///
        & inrange(endyear, 2007, 2015), local(event_years) clean
    tempfile basedata
    save `basedata', replace


    *************************************************************************************
    * Loop through all years with identified events and create a stack of data for each
    *************************************************************************************
    foreach yr of local event_years {

        use `basedata', clear
        * Event study will consider t-3, t-2, t-1, t, t+1, t+2 
        keep if inrange(endyear, `yr'-3, `yr'+2)

        * Indicator for those treated (same restrictions as above) in year == `yr'
        gen core_treat0 = (event_occur==1  ///
                            & L1null==1 & L2null==1  ///
                            & num_stu>0 & L1num_stu>0 & L2num_stu>0 ///
                            & endyear==`yr')

        * Extract their centroids (note this is an average across blocks since we collapsed to tract level)
        preserve
            keep if core_treat0==1 & endyear==`yr'
            keep censustractid lon lat
            rename (censustractid lon lat) (treatblockid tr_lon tr_lat)
            tempfile core
            save `core'
        restore

        * Identify all blocks that are within a 2.5-mile radius of the core treated blocks 
        * The complement of this list will be used as potential controls
        tempfile near
        preserve 
        keep if endyear==   `yr'
        geonear censustractid  lat lon using `core', neighbors(treatblockid tr_lat tr_lon)  miles long ignoreself
        keep if mi_to_treatblockid <= 2.5
        duplicates drop censustractid, force 
        save `near', replace
        restore
        * Merge in the nearby blocks
        merge m:1 censustractid using `near',  keep(1 3) gen(mergeTreat)

        * Control blocks are blocks farther away than the 2.5 mile radius
        tempfile control_blocks
        preserve 
            merge m:1 censustractid using `near', keep(1) nogen 
            drop treatblockid  mi_to_treatblockid
            rename censustractid treatblockid
            merge m:1 treatblockid using `core', keep(1) nogen 
            rename treatblockid censustractid 
            keep if endyear==`yr' 
            gen  control_obs = ( L1null==1 & L2null==1 ///
                                  & num_stu>0 & L1num_stu>0 & L2num_stu>0 ///
                                  & endyear==`yr')
            keep if control_obs==1 
            duplicates drop censustractid , force 
            keep censustractid  
            save `control_blocks', replace
        restore

        * Create stack-specific treatment and control indicators
        bys censustractid: egen treat = max(core_treat0)
        merge m:1 censustractid using `control_blocks', keep(1 3) gen(mergeControlTracts)
        gen control = mergeControl==3 
        drop if control == 0 & treat==0

        * Further restrict based on propensity scores if option is enabled 
        if("`pscore'"=="TRUE"){
            * Estimate logit model predicting treatment assignment 
            logit treat rel_choice_dist num_stu L1num_stu L2num_stu /// 
                        z_math_wt  ///
                        pov_rate medhhinc n_school_age white_pct black_pct  ///
                        if endyear==`yr'
            * Generate predicted probabilities
            predict phat if endyear==`yr' & (inlist(treat,1) | inlist(control,1)), pr
            sum phat if endyear==`yr' & treat==1, detail 
            * Restrict to control obs whose pscore is above the median of the treated pscore distribution 
            gen good_control =  phat > r(p50) & treat==0 & endyear==`yr'
            bys censustractid: egen keep_control = total(good_control)
            keep if keep_control==1 | treat==1 
            tab treat keep_control
        }        

        * Create event-time variable 
        gen t = endyear - `yr'
        gen stack = `yr'
        append using `stacks' 
        save `stacks', replace 
    }

    *************************************************************************************
    * Create stacked data 
    *************************************************************************************
    use `stacks', clear 

    * Event time indicators 
    gen b3_d = (t==-3)*(treat==1)
    gen b2_d = (t==-2)*(treat==1)
    gen a0_d = (t==0)*(treat==1)
    gen a1_d = (t==1)*(treat==1)
    gen a2_d = (t==2)*(treat==1)

    * Event time indicators interacted with a variable measuring intensity of treatment 
    gen focal_change = change_rel if t==0
    bys stack censustractid : egen focal_change_stack = mean(focal_change )
    replace focal_change_stack =  focal_change_stack*(-1)
    gen b3_d2 = (t==-3)*(treat==1)*focal_change_stack
    gen b2_d2 = (t==-2)*(treat==1)*focal_change_stack
    gen a0_d2 = (t==0)*(treat==1)*focal_change_stack
    gen a1_d2 = (t==1)*(treat==1)*focal_change_stack
    gen a2_d2 = (t==2)*(treat==1)*focal_change_stack
    
    * Stack-specific fixed effects 
    egen stack_yr = group(stack t) 
    egen stack_group = group(stack censustractid )

    * Create some variables to assess robustness 
    gen share_applied = appliedMag /num_stu 
    gen num_stu_mag_t = num_stu *t
    gen num_stu_mag_t2 = num_stu *t^2
    gen num_stu_mag_t3 = num_stu* t^3
    local vars "pov_rate medhhinc n_school_age white_pct black_pct hispanic_pct "
    foreach var of local vars{
        gen t_`var' = `var' * t
        gen t2_`var' = `var' * t^2
        gen t3_`var' = `var' * t^3
        gen t4_`var' = `var' * t^4
    }

    * Note a few zeros are dropped but very few, so just use logs 
    gen lnapplied = ln(applied)

    save $DTAFINAL/event_study_data_tract_relative`relative'.dta, replace 

end 



capture program drop prep_event_study_block 
program define prep_event_study_block

    clear 
    tempfile stacks 
    save `stacks', emptyok

    *************************************************************************
    * Create a block-year level dataset with mag apps and other vars 
    *************************************************************************
    use "$BUILD/output/data/lausd2001_2023_encoded.dta", clear    
	gen appliedMag = mergeMagnet ==3 
    keep if gradecode ==5 | gradecode ==6 | gradecode ==7
    gen stu=1
    decode censusblockid1 , gen(blockid) 
    drop if missing(blockid)
    gen app_math = z_math_all if appliedMag ==1
    gen app_ela = z_ela_all if appliedMag ==1
	gen app_black = (ethnicitydescription == 2) if appliedMag ==1
	gen app_hispanic = (ethnicitydescription == 1) if appliedMag ==1
	gen app_white = (ethnicitydescription == 3) if appliedMag ==1
	gen app_eng_home = (eng_home==1) if appliedMag ==1
	gen app_female = (gendercode==1) if appliedMag ==1
	gen app_suspensions = (total_num_suspensions) if appliedMag ==1

    collapse (sum) appliedMag (count) num_stu = stu ///
		(mean) z_math_all z_ela_all app_math app_ela ///
			app_black app_hispanic app_white app_eng_home app_female app_suspensions , by(blockid endyear) 
    destring blockid, replace 
    format blockid %20.0g
    tempfile temp 
    save `temp'

    * List of relevant block ids 
    contract blockid 
    drop _freq 
    tempfile temp2 
    save `temp2'

    *************************************************************************************
    * Create census block centroid data (for distance calculations later)
    *************************************************************************************
    *cd /Users/cqcampos/Library/CloudStorage/Dropbox/research/magnet_structural/rawdata/
    *spshape2dta "tl_2020_06_tabblock20.shp", replace saving(blocks_ca)
    * Attribute file with GEOID and centroids
    use $DTARAW/blocks_ca.dta, clear
        keep if COUNTYFP20 =="037"
    keep GEOID20 _CX _CY
    rename (GEOID20 _CX _CY) (blockid lon lat)    // lon/lat in decimal degrees (NAD83/WGS84)
    tostring blockid, replace force  // if needed
    tempfile blk_centroids
    save `blk_centroids'



    *************************************************************************************
    * Read in block-to-school distance data (created in the LAUSD master build directory)
    *************************************************************************************    
    use "$BUILD/rawdata/aux_data/relative_choice_distances_grade5.dta", clear    
	rename censusblockid blockid 
    merge m:1 blockid using `blk_centroids', keep(3) nogen
    assert !missing(lat, lon)
    rename blockid censusblockid

    destring censusblockid , gen(blockid)
    format blockid %20.0g
    drop if missing(blockid)
    xtset blockid endyear

    gen change_rel = rel_choice_dist  - L1.rel_choice_dist 

    * replace change_rel = 0 if change_rel==.
    * Those that do not merge (==1) are blocks in LA county that LAUSD students do not live in
    merge m:1 blockid using `temp2', gen(mergeExistBlock) keep(3)
    merge 1:1 blockid endyear using `temp', gen(mergeApps) keep(1 3)
    replace appliedMag = 0 if missing(appliedMag)
    replace num_stu = 0 if missing(num_stu)

    *************************************************************************************
    * Create additional vars to find events 
    *************************************************************************************
    xtset blockid endyear
    gen L1null = L2.min_choice_dist  == L1.min_choice_dist  
    gen L2null = L3.min_choice_dist  == L2.min_choice_dist  
    gen L3null = L4.min_choice_dist  == L3.min_choice_dist 
    
    gen L1num_stu = L1.num_stu
    gen L2num_stu = L2.num_stu
    gen L3num_stu = L3.num_stu

    * Identify potential events 
    gen event_occur = change_rel < -0.5 & endyear !=2004 if !missing(change_rel)


     *************************************************************************************
    * Identify candidate event years
    * Restrictions: 
    * 1) Avg. tract-level reduction in relative distance is greater than 1 mile (event_occur==1)
    * 2) No changes in proximity to choice schools in t-1, t-2, and t-3 (L1null==1 & L2null==1)
    * 3) The tract has enrolled students at t, t-1, t-2, and t-3 (L1num_stu>0 & L2num_stu>0)
    * 4) Restrict to years that we use in the structural estimation (+ 2)
    *************************************************************************************
    levelsof endyear if event_occur==1 ///
        & L1null==1 & L2null==1 & num_stu>0 & L1num_stu>0 & L2num_stu>0 ///
        & inrange(endyear, 2007, 2017), local(event_years) clean
    tempfile basedata
    save `basedata', replace

    *************************************************************************************
    * Loop through all years with identified events and create a stack of data for each
    *************************************************************************************
    foreach yr of local event_years {
     use `basedata', clear
        * Event study will consider t-3, t-2, t-1, t, t+1, t+2 
        keep if inrange(endyear, `yr'-3, `yr'+2)

        * Indicator for those treated (same restrictions as above) in year == `yr'
        gen core_treat0 = (event_occur==1  ///
                            & L1null==1 & L2null==1  ///
                            & num_stu>0 & L1num_stu>0 & L2num_stu>0 ///
                            & endyear==`yr')

        * Extract their centroids (note this is an average across blocks since we collapsed to tract level)
        preserve
            keep if core_treat0==1 & endyear==`yr'
            keep censusblockid lon lat
            rename (censusblockid lon lat) (treatblockid tr_lon tr_lat)
            tempfile core
            save `core'
        restore

        * Identify all blocks that are within a 2.5-mile radius of the core treated blocks 
        * The complement of this list will be used as potential controls
        tempfile near
        preserve 
        keep if endyear==   `yr'
        geonear censusblockid  lat lon using `core', neighbors(treatblockid tr_lat tr_lon)  miles long ignoreself
        keep if mi_to_treatblockid <= 5
        duplicates drop censusblockid, force 
        save `near', replace
        restore
        * Merge in the nearby blocks
        merge m:1 censusblockid using `near',  keep(1 3) gen(mergeTreat)

        * Control blocks are blocks farther away than the 2.5 mile radius
        tempfile control_blocks
        preserve 
            merge m:1 censusblockid using `near', keep(1) nogen 
            drop treatblockid  mi_to_treatblockid
            rename censusblockid treatblockid
            merge m:1 treatblockid using `core', keep(1) nogen 
            rename treatblockid censusblockid
            keep if endyear==`yr' 
            gen  control_obs = ( L1null==1 & L2null==1 ///
                                  & num_stu>0 & L1num_stu>0 & L2num_stu>0 ///
                                  & endyear==`yr')
            keep if control_obs==1 
            duplicates drop censusblockid , force 
            keep censusblockid  
            save `control_blocks', replace
        restore

        * Create stack-specific treatment and control indicators
        bys censusblockid: egen treat = max(core_treat0)
        merge m:1 censusblockid using `control_blocks', keep(1 3) gen(mergeControlTracts)
        gen control = mergeControl==3 
        drop if control == 0 & treat==0

        * Create event-time variable 
        gen t = endyear - `yr'
        gen stack = `yr'
        append using `stacks' 
        save `stacks', replace 
    }

    *************************************************************************************
    * Create stacked data 
    *************************************************************************************
    use `stacks', clear 

    * Event time indicators 
    gen b3_d = (t==-3)*(treat==1)
    gen b2_d = (t==-2)*(treat==1)
    gen a0_d = (t==0)*(treat==1)
    gen a1_d = (t==1)*(treat==1)
    gen a2_d = (t==2)*(treat==1)

    * Event time indicators interacted with a variable measuring intensity of treatment 
    gen focal_change = change_rel if t==0
    bys stack censusblockid : egen focal_change_stack = mean(focal_change )
    replace focal_change_stack =  focal_change_stack*(-1)
    gen b3_d2 = (t==-3)*(treat==1)*focal_change_stack
    gen b2_d2 = (t==-2)*(treat==1)*focal_change_stack
    gen a0_d2 = (t==0)*(treat==1)*focal_change_stack
    gen a1_d2 = (t==1)*(treat==1)*focal_change_stack
    gen a2_d2 = (t==2)*(treat==1)*focal_change_stack
    
    * Stack-specific fixed effects 
    egen stack_yr = group(stack t) 
    egen stack_group = group(stack censusblockid )

    save $DTAFINAL/event_study_data_block.dta, replace 

end 

main
















