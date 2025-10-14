/********************************************************************************
- Author: Chris Campos
- Antonia Edited the path files to follow the folders restructuration 
- Description: some preliminary analysis of magnet programs 
- Date started: 12/17/2024
- Last update: 07/24/2025
*******************************************************************************/
clear all
set trace on
set tracedepth 2
set maxvar 120000

capture program drop main
program define main

        set_paths
	
	local minY = 2004
	local maxY = 2013	
	
       * prepStudentData, minYr(`minY') maxYr(`maxY')
    local maxY = 2017
        prepStudentData, minYr(`minY') maxYr(`maxY')
        stop;
    /*
        * Extrat subset of applications for relevant years 
        prepApplications, minYr(`minY') maxYr(`maxY')

        /* Prep magnet schools */ 
        isolateSchools, minYr(`minY') maxYr(`maxY')
        isolateApplications, minYr(`minY') maxYr(`maxY')
	
        /* Calculate distances */
        calculateDistances, minYr(`minY') maxYr(`maxY') 
    */
        /* Prep structural data */
       prepStructuralData, minYr(`minY') maxYr(`maxY')
    /*
        /* Admission probabilities */ 
	   prepareSeats, minYr(`minY') maxYr(`maxY')
       student_p_i, minYr(`minY') maxYr(`maxY')
	*/
	   
	   /* Make x-walk between program IDs from 2004-2008 sample and 2004-2013 sample */
	   *cc_xwalk_08_13
	   
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
        global ROOT /Volumes/lausd/decentralized_choice 
        global BUILD /Volumes/lausd/build/output/data
        global ROOTSOURCE "/Volumes/lausd/build/" 
        global magnet "/Volumes/lausd/magnet/"
    }
	else if "`c(os)'" == "Windows"{
		global ROOT "Z:/decentralized_choice"
        global BUILD "Z:/build/output/data"
        global ROOTSOURCE "Z:/build/" 
        global magnet "Z:/magnet/"
	}
    else{
            global ROOT "/project/lausd/decentralized_choice"
            global BUILD "/project/lausd/build/output/data"
            global ROOTSOURCE "/project/lausd/build/" 
            global magnet "/project/lausd/magnet"

    }
	
    global DTARAW "$ROOT/rawdata" 
    global RAW $ROOT/data/raw
    global DTAFINAL "$ROOT/data" 
    global tables "$ROOT/tables"
    global figures "$ROOT/figures"

    * Raw source data
    global UE "$ROOTSOURCE/rawdata/Campos_Master_DUA/Unified_Enrollment"
    global RAWOUT "$ROOTSOURCE/output/data" 


end


********************************************************************************
* prepStudentData 
* 
* Description: Preps relevant fifth grade cohort student data 
********************************************************************************
capture program drop prepStudentData 
program define prepStudentData 
syntax, [minYr(int 2004) maxYr(int 2017)]
	/*
    use "$BUILD/magnet_rank.dta", clear
    keep if choice_ranking ==1
    contract mag_cd_ endyear
    drop if _freq<=10 // those that have less than 10 apps are usually data entry errors
    drop _freq 
    rename mag_cd_ next_scode 
    drop if next_scode==. 
    tempfile tempcodes 
    save `tempcodes'


    * Create offer dummies (not conditioning on lottery) THIS DATA IS NOW CREATED IN 0_build.do and store in "./decentralized_choice/data/raw/magnet_offers.dta"
    import delimited "$UE/Deidentified_NEW_MAG_ALL.csv", clear stringcols(_all)
    destring application_id school_yr, replace
    merge 1:1 application_id using "$RAWOUT/apps_time_stamp.dta"
    drop _merge
    replace app_endyear = school_yr if missing(app_endyear) // so we don;t lose the 2023, 2024 people
    sort student_pseudo_id
    rename (student_pseudo_id ) (studentpseudoid) 
    gen magnet_offer = mag_action_1=="B" | mag_action_2 =="B" | mag_action_3=="B"
    contract studentpseudoid app_endyear magnet_offer
    gen endyear = app_endyear //just to use for the match in the build code 
    drop _freq 
    duplicates drop studentpseudoid endyear, force 
    save "$magnetINT/magnet_offers.dta", replace 
	*/
	
	local minYr = 2004
	local maxYr = `maxYr'
	
    * Read in main sample, merge magnet offers, and create necessary variables 
    use if endyear >=`minYr' & endyear <=(`maxYr' +4)  using $RAWOUT/lausd2001_2023_encoded.dta, clear 
    merge m:1 studentpseudoid endyear using $RAW/magnet_offers, gen(mergeOffers) keep( 1 3)
    drop mergeOffers 
    gen appliedMagnet = mergeMagnet ==3
    egen sid = group(studentpseudoid)

    xtset sid endyear 
    decode schoollocationname, gen(sname)
    gen next_schoolname = F1.schoollocationname
    label values next_schoolname schoollocationname_
    gen next_scode = F1.schoolcostcentercode 
    decode next_schoolname, gen(next_sname)
    gen F1math = F1.z_math_all 
    gen F1ela = F1.z_ela_all
    gen F2math = F2.z_math_all
    gen F2ela = F2.z_ela_all
    gen F3math = F3.z_math_all
    gen F3ela = F3.z_ela_all
    egen numSuspensions = rowtotal(num_of_suspensions_Class num_of_suspensions_In_school num_of_suspensions_Out_of_school)
    replace numSuspensions =0 if missing(numSuspensions)
    gen F1suspensions = F1.numSuspensions
    gen F2suspensions = F2.numSuspensions
    gen F3suspensions = F3.numSuspensions
    replace ytdofabsentdays =0 if missing(ytdofabsentdays)
    gen F1absent = F1.ytdofabsentdays
    gen F2absent = F2.ytdofabsentdays
    gen F3absent = F3.ytdofabsentdays
    gen F1gpa = F1.gpa_fall 
    gen F2gpa = F2.gpa_fall
    gen F3gpa = F3.gpa_fall


    egen avg_Fmath = rowmean(F1math F2math F3math)
    egen avg_Fela = rowmean(F1ela F2ela F3ela)
	egen avg_Fsuspensions = rowmean(F1suspensions F2suspensions F3suspensions)
    egen avg_Fabsent = rowmean(F1absent F2absent F3absent)
    egen avg_Fgpa = rowmean(F1gpa F2gpa F3gpa)
    bys endyear gradecode: egen mean_suspensions = mean(avg_Fsuspensions)
    bys endyear gradecode: egen mean_absent = mean(avg_Fabsent)
    bys endyear gradecode: egen mean_gpa = mean(avg_Fgpa)
    bys endyear gradecode: egen sd_suspensions = sd(avg_Fsuspensions)
    bys endyear gradecode: egen sd_absent = sd(avg_Fabsent)
    bys endyear gradecode: egen sd_gpa = sd(avg_Fgpa)

    gen z_suspensions = (avg_Fsuspensions - mean_suspensions) / sd_suspensions

    gen z_absent = (avg_Fabsent - mean_absent) / sd_absent
    gen z_gpa = (avg_Fgpa - mean_gpa) / sd_gpa
    replace z_suspensions = - z_suspensions
    replace z_absent = - z_absent

    
    egen z_socio = rowmean(z_suspensions z_absent z_gpa)



	keep if endyear >=`minYr' & endyear <=`maxYr'

    gen hispanic = ethnicitydescription==1
    gen black = ethnicitydescription==2
    gen white = ethnicitydescription==3
    gen asian= ethnicitydescription==4 
    gen other = hispanic==0 & black==0 & white==0 & asian==0 
    gen female = gendercode==1
    gen poverty = studentpovertyindicator==3
    gen el = languageclasscode==3
    gen gifted = !missing(studentgiftedprogramdescription) 

    gen sixthGradeAddress = censusblockid1 if gradecode ==6
    bys studentpseudoid : egen addressSix = mean(sixthGradeAddress )
    replace censusblockid1 = addressSix if missing(censusblockid1 ) & gradecode ==5
    keep if gradecode==5

    gen schooltype = "ms" if gradecode ==5
    
    rename censusblockid1 censusblockid
    decode censusblockid, gen(blockid)  
    rename (censusblockid blockid) (blockid censusblockid)
    merge m:1 censusblockid schooltype endyear using $DTAFINAL/intermediate/magnet_distances.dta  , gen(mergeDistances) keep(1 3)

    gen next_mag = regexm(next_sname , "MAG") | regexm(next_sname, "MG")    


    drop city2 api_city2 api_state2 api_zip_code2 censusblockid2 city3 api_city3 api_state3 api_zip_code3 censusblockid3 ///
        city4 api_city4 api_state4 api_zip_code4 censusblockid4 city5 api_city5 api_state5 api_zip_code5 censusblockid5 ///
        city6 api_city6 api_state6 api_zip_code6 censusblockid6 city7 api_city7 api_state7 api_zip_code7 censusblockid7 ///
        city8 api_city8 api_state8 api_zip_code8 censusblockid8 city9 api_city9 api_state9 api_zip_code9 censusblockid9 ///
        sbac_ela sbac_math mergeSBAC ss_cst_ela1 ss_cst_math1 ss_cst_ela2 ss_cst_math2 numTests testid mergeCST ///
        effective_date2 zip_code2 effective_date3 zip_code3 effective_date4 zip_code4 effective_date5 zip_code5 effective_date6 zip_code6 effective_date7 zip_code7 effective_date8 zip_code8 effective_date9 zip_code9 mergeAddresses ///
        ethnicitydescription gendercode parentedulevelname preferredlocationname schoollocationname studentbirthcountry studentusschoolfirstattenddate languageclasscode homelanguagedescription studentspedflag studentgiftedprogramdescription studentpovertyindicator overallperformancelevelnameELA overallperformancelevelnameMath testnameela1 testnamemath1 testnameela2 testnamemath2 perf_cst_ela1 perf_cst_math1 perf_cst_ela2 perf_cst_math2 ///
        sixthGradeAddress addressSix 

    save $DTAFINAL/lausd_5thgrade_cohorts_`minYr'_`maxYr'.dta, replace

end 


capture program drop prepApplications 
program define prepApplications
syntax, [minYr(int 2004) maxYr(int 2017)]


    use "$BUILD/magnet_rank.dta", clear
    **In this case these are the same because we want to know who apply to a magnet when they were in 5th year
    gen endyear = school_yr
    keep if choice_ranking ==1 
    keep if endyear >=`minYr' & endyear <= `maxYr'
    keep if stu_grade==6


    save $DTAFINAL/applications_5thgrade_cohorts_`minYr'_`maxYr'.dta, replace
end 



/********************************************************************************
* isolateSchools 
* Description: There are many programs middle school students apply to that are 
* elementary programs or other highly gifted programs with additional screening 
* criteria. Drop these programs and keep only those that are middle school magnets. 
* This will be the list of magnets that we use moving forward. 
********************************************************************************/
capture program drop isolateSchools
program define isolateSchools
syntax, [minYr(int 2004) maxYr(int 2017)]

    use $DTAFINAL/applications_5thgrade_cohorts_`minYr'_`maxYr'.dta, clear 
    drop if regexm(mag_name_ , "MS")
    drop if regexm(mag_name_ , " EL") | regexm(mag_name_ , " CM") | regexm(mag_name_ , "ACAD")
    drop if regexm(mag_name_ , " CA") | regexm(mag_name_ , " LC") | regexm(mag_name_ , "ACAD")
    drop if regexm(mag_name_, "SVCS-CO")
    drop if regexm(mag_name_, "CONTRACT ISSUED")
    drop if regexm(mag_name_, "CHATSWORTH LAAMP")
    gen obs = 1 
    bys mag_cd_: egen numApps = total(obs)
    tab numApps
    keep if numApps>15

    bys endyear: egen nschools = nvals(mag_cd_ )
    egen nschoolsall = nvals(mag_cd_ )

    tab nschools 
    tab nschoolsall 

    * These magnets do not exist in enrollment data -- all are elementary schools
    drop if mag_cd_== 1226901 | mag_cd_== 1493201 | mag_cd_== 1713702 | mag_cd_== 1861402 | mag_cd_== 1875401 | mag_cd == 1206802  
    
    * The following are highly gifted and have other screening criteria 
    drop if mag_cd == 1350702 | mag_cd == 1350703
    
    * The following are elementary schools (checked manually)
    drop if inlist(mag_cd, 1250701, 1258902, 1274101, 1367102, 1461602, 1531502, 1687502)
    
    * Ryan: Between 2009 and 2013, mag_name_has trailing spaces (e.g., "LA MAG " instead of "LA MAG")
    * Remove the trailing spaces
    replace mag_name_ =  trim(mag_name_)	
    contract mag_cd_ mag_name_ 
    rename mag_cd_ schoolcostcentercode
   
     drop if missing(schoolcostcentercode)
    save $DTAFINAL/magnet_codes_`minYr'_`maxYr'.dta, replace
end 


capture program drop isolateApplications 
program define isolateApplications
syntax, [minYr(int 2004) maxYr(int 2017)]

    set sortseed 123456 

    use $DTAFINAL/applications_5thgrade_cohorts_`minYr'_`maxYr'.dta, clear 
    drop if regexm(mag_name_ , "MS")
    drop if regexm(mag_name_ , " EL") | regexm(mag_name_ , " CM") | regexm(mag_name_ , "ACAD")
    drop if regexm(mag_name_ , " CA") | regexm(mag_name_ , " LC") | regexm(mag_name_ , "ACAD")

    duplicates drop studentpseudoid endyear, force 
    rename mag_cd_ schoolcostcentercode 
    merge m:1 schoolcostcentercode using $DTAFINAL/magnet_codes_`minYr'_`maxYr'.dta, gen(mergeMagnet) keep(3)
    rename schoolcostcentercode mag_cd_
    save $DTAFINAL/applications_`minYr'_`maxYr'.dta, replace 

    gen offer = mag_action_=="B"
    keep if offer ==1 
    rename mag_cd_ mag_offered_code
    keep studentpseudoid endyear mag_offered
    duplicates drop studentpseudoid endyear, force
    save $DTAFINAL/offers_`minYr'_`maxYr'.dta, replace
end 


capture program drop calculateDistances
program define calculateDistances
syntax, [cyear(string) minYr(int 2004) maxYr(int 2017)]

	cap  log close 

	import delimited "$ROOTSOURCE/rawdata/raw_base/pubschls.csv", clear 
	tempfile cde 
	save `cde' 

	import delimited "$ROOTSOURCE/rawdata/raw_base/blocks_la.csv", clear stringcols(_all) varnames(1)
	destring lon lat, replace  float 
	sort censusblockid 
	gen blockid = _n 
	count 
	local Nblocks = r(N)
	format lat %20.0g
	format lon %20.0g
	cap drop block_name
	rename (lat lon) (block_lat block_lon)
	tempfile blocks 
	save `blocks'
	
	*2

		use if endyear>=`minYr' & endyear<=(`maxYr'+1) using $RAWOUT/lausd2001_2023_encoded.dta, clear 
		keep if gradecode==6 
		 
		collapse (count) enrollment = endyear, by(preferredlocationcode schoolcostcentercode gradecode) 
		merge m:1 preferredlocationcode using "$ROOTSOURCE/rawdata/raw_base/plocn_cde_xwalk.dta", gen(mergeCDEcodes) keep(3) 		
		
		* Four obs are unmatched (check website to manually assign)
		replace cdecode = 0140046 if preferredlocationcode==7895
		replace cdecode = 6058028 if preferredlocationcode==7567
		replace cdecode = 6057962 if preferredlocationcode==7566
		tostring cdecode , replace
		replace cdecode = "0" + cdecode if length(cdecode )==6
		replace cdecode = "1964733" + cdecode
		gen cdscode = cdecode
		merge m:1 cdscode using `cde', keep(1 3) gen(mergeCDEInfo)
  
		* Keep only magnet programs that we have flagged 
		* If _merge == 2 exists, some schools that are applied by students are problematic 
		merge m:1 schoolcostcentercode using "$DTAFINAL/magnet_codes_`minYr'_`maxYr'.dta", gen(mergeMagnet) keep(2 3)
	

		* create an alternate identifier to use in expand 
		*drop preferredlocationcode
		gen grade_level = gradecode 
		rename schoolcostcentercode school_number
		sort grade_level school_number
		bys grade_level: gen auxid = _n

		* keep relevent variables to merge later 
		keep school_number  grade_level cdscode  latitude longitude  
		order  grade_level, after(school_number)
		foreach var of varlist cdscode latitude longitude {
			rename `var' receiving_`var'
		}
		gen feeder_grade = grade_level 
		replace grade_level = grade_level-1

		* drop additional Special Ed, Continuation, Opportunity and adult schools 
		*drop if inlist(receiving_edopsname, "Community Day School", "Continuation School", "Opportunity School", "Special Education School")
		tempfile receiving_schools 
		drop if receiving_latitude=="No Data" | receiving_longitude=="No Data"
		drop if receiving_latitude=="NA" | receiving_longitude=="NA"
		destring  receiving_latitude receiving_longitude , replace 
		save `receiving_schools' 
		

        use `receiving_schools', clear 
        keep if grade_level==5
        count 
        
        local Nexpand = r(N) 
        expand `Nblocks'
        sort school_number
        bys school_number: gen blockid = _n 
        merge m:1 blockid using `blocks', gen(mergeBlocks) 
        geodist block_lat block_lon receiving_latitude receiving_longitude, gen(dist) miles  
        sort censusblockid dist 

        * Identify magnet schools

        cap destring receiving_cdscode , replace
        cap format receiving_cdscode %20.0g 
        drop feeder_grade pop receiving_cdscode
        
        keep censusblockid block_lat block_lon dist school_number 
        rename dist d_ 
        reshape wide d_, i(censusblockid) j(school_number)
        save $DTAFINAL/distance_to_magnets_`minYr'_`maxYr'.dta, replace 




end 



capture program drop prepStructuralData
program define prepStructuralData
syntax, [minYr(int 2004) maxYr(int 2017)]	

	*AV 07/24/2025 when using the command getcensus i was getting a weird error about the 
	* directoy, this was my solution but may need to be adapted for different users
	*! mkdir -p /home/avazque0/.local/share/getcensus/
	* ! ls -la /home/avazque0/.local/share/getcensus/
	***
	getcensus B19013, geography(tract) sample(5) state("06") clear year(2020)
    gen censustract = state + county + tract 
    rename b19013_001e median_income 
    keep censustract median_income
    replace median_income = 0 if missing(median_income)
    replace median_income = median_income/10000
    tempfile census_tracts
    save `census_tracts', replace

    
    
    
	local minYr =2004 
	local maxYr = 2013
    use $DTAFINAL/magnet_codes_`minYr'_`maxYr'.dta, clear 
    levelsof schoolcostcentercode, local(mag_codes) clean 
    
    use $DTAFINAL/lausd_5thgrade_cohorts_`minYr'_`maxYr'.dta, clear 
    merge m:1 studentpseudoid endyear using $DTAFINAL/applications_`minYr'_`maxYr'.dta, keepusing(mag_cd_) gen(mergeApps) keep(1 3)
    rename mag_cd_ mag_applied_code 
    merge m:1 studentpseudoid endyear using $DTAFINAL/offers_`minYr'_`maxYr'.dta, keepusing(mag_offered_code) gen(mergeOffers) keep(1 3)
     * Organize the data in a particular order for the structural code 
    local covs "female black white hispanic poverty el eng_home"
    order endyear `covs', first 
    local scores "F1math F1ela F2math F2ela F3math F3ela avg_Fmath avg_Fela"
    order `scores', after(eng_home)
    rename (z_math_all z_ela_all) (lag_math lag_ela)
    bys endyear: egen mean_math = mean(lag_math)
    bys endyear: egen mean_ela = mean(lag_ela)
    replace lag_math = lag_math - mean_math
    replace lag_ela = lag_ela - mean_ela
    gen missing_math = missing(lag_math)
    gen missing_ela = missing(lag_ela)

    replace lag_math = 0 if missing_math == 1
    replace lag_ela = 0 if missing_ela == 1
    local lags "lag_math lag_ela missing_math missing_ela"
    order `lags', after(avg_Fela)

    gen Dchoice = next_scode 
    order Dchoice, after(missing_ela)

    * Create enrollment and application dummies 
    foreach code of local mag_codes{
        gen D_`code' = next_scode == `code'
    }
    foreach code of local mag_codes{
        gen A_`code' = mag_applied_code == `code'
    }
    foreach code of local mag_codes{
        gen Z_`code' = mag_offered_code == `code'
    }
    gen Achoice= mag_applied_code
    replace Achoice = 0 if missing(Achoice)


    egen enroll_one = rowtotal(D_*)
    gen D_0 = enroll_one ==0 
    drop enroll_one
    
    order Achoice, after(Dchoice)
    order D_* A_* Z_*, after(Dchoice)
    order D_0, after(Dchoice)

    * Keep students with addresses 
    merge m:1 censusblockid using $DTAFINAL/distance_to_magnets_`minYr'_`maxYr'.dta, gen(mergeDistance) keep(3)

    order d_*, after(Z_1884201)
    foreach code of local mag_codes{
        replace d_`code' = d_`code' - nearest_schl_dist
    }

    order block_lat block_lon, after(d_1884201)

    drop if missing(avg_Fela) | missing(avg_Fmath)

     * DROP STUDENTS WHO ENROLLED IN A SCHOOL THEY DIDN'T RECEIVE AN OFFER FROM
    gen inconsistent_enrollment = 0
    foreach code of local mag_codes {
        replace inconsistent_enrollment = 1 if D_`code' == 1 & Z_`code' == 0
    }
    count if inconsistent_enrollment == 1
    local dropped = r(N)
    drop if inconsistent_enrollment == 1
    drop inconsistent_enrollment
    display "Dropped `dropped' students who enrolled in schools they didn't receive offers from"


    drop studentclassofname esp_home ytdofattendeddays ytdofenrolleddays ytdofabsentdays otherabsencecount ytdofattendance appabsencecount mergeATT num_of_suspensions_Class num_suspended_days_Class num_of_suspensions_In_school num_suspended_days_In_school num_of_suspensions_Out_of_school num_suspended_days_Out_of_school total_num_suspended_days total_num_suspensions mergeSuspensions gpa_fall gpa_spring mergeGPA

    drop mergeAFC mergeMagnet mergeDLE mergeACS mergePWT mergeSAS cst_ela cst_math z_cst_ela z_cst_math z_sbac_ela z_sbac_math
    
    
    gen censustract = substr(censusblockid, 1, 11) 
    merge m:1 censustract using `census_tracts', gen(mergeCensus) keep(1 3)

    gen any_missing_score = missing_math==1 | missing_ela ==1
    gen current_in_mag = regexm(sname, "MAG") | regexm(sname, "MG") | regexm(sname, "SOCES") | regexm(sname, "LACES")
    bys sname endyear : egen current_peer_quality =mean(lag_math)
    
    * For now have a subset of applicants that apply to magnets that we have not flagged as magnets (3%)
    * Drop for now 
    gen applied = Achoice !=0
    drop if applied==0 & appliedMagnet==1

    sort endyear studentpseudoid
    save $DTAFINAL/structural_data_`minYr'_`maxYr'.dta, replace

    gen rand = runiform()
    sort rand 
    keep if rand <=0.25
    drop rand 
    sort endyear studentpseudoid
    save $DTAFINAL/structural_data_`minYr'_`maxYr'_sample.dta, replace

end 



capture program drop prepareSeats
program define prepareSeats 
syntax, [minYr(int 2004) maxYr(int 2017)]

    use $BUILD/magnet_cutoffs/prog_seats_app.dta, clear 
    rename mag_cd_ schoolcostcentercode 
    rename app_endyear endyear
    keep if endyear <=`maxYr' & endyear>=`minYr' 
    keep if stu_grade==6 
   

    merge m:1 schoolcostcentercode using $DTAFINAL/magnet_codes_`minYr'_`maxYr'.dta, gen(mergeSchools) keep(3)
    preserve 
        keep endyear schoolcostcentercode stu_grade phbao seats_available
        rename seats_available seats_ 
        reshape wide seats_, i(endyear stu_grade phbao) j(schoolcostcentercode)
        ds seats_*
        foreach var of varlist `r(varlist)'{
            replace `var' = 0 if missing(`var')
        }
        save $DTAFINAL/prog_seats_`minYr'_`maxYr'.dta, replace
    restore
    preserve 
        collapse (sum) applications seats_available, by(schoolcostcentercode)
        gen p = seats_available / applications
        replace p = 1 if seats_available> applications & !missing(seats_available)
        save $DTAFINAL/aggregate_prog_admission_prob_`minYr'_`maxYr'.dta, replace
    restore 
 
    gen p_i = seats_available /applications 
    replace p_i = 1 if seats_available> applications & !missing(seats_available)
    replace p_i = 0 if missing(p_i)
    

    keep endyear schoolcostcentercode stu_grade phbao p_i

    rename p_i p_i_
    reshape wide p_i_ , i(endyear stu_grade phbao ) j(schoolcostcentercode )
    ds p_i_*
    foreach var of varlist `r(varlist)'{
        replace `var' = 0 if missing(`var')
    }
    

    save $DTAFINAL/prog_seats_app_`minYr'_`maxYr'.dta, replace

end 


capture program drop student_p_i 
program define student_p_i 
syntax, [minYr(int 2004) maxYr(int 2017)]

    use $DTAFINAL/structural_data_`minYr'_`maxYr'.dta, clear
    gen phbao = (black==1 | asian==1 | hispanic==1 | other==1)
    gen stu_grade = 6 
    merge m:1 phbao stu_grade endyear using $DTAFINAL/prog_seats_app_`minYr'_`maxYr'.dta, gen(mergeP) keep( 1 3)
    sort endyear studentpseudoid
    keep studentpseudoid endyear p_i_*

    save $DTAFINAL/student_p_i_`minYr'_`maxYr'.dta, replace

    use $DTAFINAL/structural_data_`minYr'_`maxYr'_sample.dta, clear
    gen phbao = (black==1 | asian==1 | hispanic==1 | other==1)
    gen stu_grade = 6 
    merge m:1 phbao stu_grade endyear using $DTAFINAL/prog_seats_app_`minYr'_`maxYr'.dta, gen(mergeP) keep(1 3)
    sort endyear studentpseudoid
    keep studentpseudoid endyear p_i_*
    save $DTAFINAL/student_p_i_`minYr'_`maxYr'_sample.dta, replace
end 


capture program drop cc_xwalk_08_13
program define cc_xwalk_08_13

    use $DTAFINAL/magnet_codes_2004_2013.dta, clear 
	
	rename _freq n_apps
	rename mag_name_ mag_name 
	rename schoolcostcentercode costcentercode
	
	tempfile magnames 
	save `magnames'
	
	* Grab cost center ID from 2004-2008(2013) applications 
	forval year = 2008(5)2013{
		use "$DTAFINAL/student_p_i_2004_`year'_sample.dta", clear
		keep if _n == 1
		keep p_i_*
		
		gen last_yr = `year'
		
		reshape long p_i_, i(last_yr) j(costcentercode)
		drop p_i*
		gen j_`year' = _n
			
		tempfile j_`year'
		save `j_`year''
	}
	
	merge 1:1 costcentercode using `j_2008', nogen
	drop last_yr
	
	* Get the name for these costcentercode
	
	merge 1:1 costcentercode using `magnames', nogen
	export delimited "$DTAFINAL/cc_xwalk_08_13.csv", replace
	
end 

main

