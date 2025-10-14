/********************************************************************************
5_outcome_analysis.do
- Author: Chris Campos 
- Description: Summarizes results from counterfactual simulations
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

	*set_paths 
 
    counterfactuals

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

capture program drop counterfactuals
program define counterfactuals

    * Total number of seats in the market 
    use "/Volumes/lausd/decentralized_choice/data/prog_seats_2004_2013.dta", clear 

    reshape long seats_, i(endyear phbao stu_grade )
    collapse (sum) seats_ 
    sum seats_ 
    local total_seats = r(mean)

    * Total number of students 
    use "/Volumes/lausd/decentralized_choice/data/structural_data_2004_2013.dta", clear
    sum F1math  
    local ybar_math_base = r(mean)
    sum F1ela
    local ybar_ela_base = r(mean)
    count 
    local total_students = r(N)

    * Read in counterfactual simulation estimates 
    tempfile ests 
    * General info campaign
    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_2_used_eq_pi_2013_maxapp_1.csv", clear
    save `ests', replace
    * Targeted info campaign
    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_3_used_eq_pi_2013_maxapp_1.csv", clear
    append using `ests'
    save `ests', replace
    * Centralized UE 
    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_4_2013_maxapp_1.csv", clear
    replace n_over_sub_j = ""  
    destring n_over_sub_j, replace 
    append using `ests'
    save `ests', replace
    * Baseline allocation
    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_0_2013_maxapp_1.csv", clear
    append using `ests'
    save `ests', replace
    * Centralized UE + no travel costs
    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_6_used_eq_pi_2013_maxapp_1.csv", clear
    replace amte_math = ""
    replace amte_ela = ""
    destring amte_math amte_ela, replace 
    replace n_over_sub_j = ""  
    destring n_over_sub_j, replace 
    append using `ests' 
    save `ests', replace

    * No travel costs 
    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_5_used_eq_pi_2013_maxapp_1.csv", clear
    replace amte_math = ""
    replace amte_ela = ""
    destring amte_math amte_ela, replace 
    append using `ests' 

    save `ests', replace 

    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_7_used_eq_pi_2013_maxapp_1.csv", clear
    replace app_rate = ""
    replace offer_rate = ""
    replace n_over_sub_j = ""  
    destring n_over_sub_j app_rate offer_rate, replace 
 
    append using `ests'

    save `ests', replace 


    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_8_used_eq_pi_2013_maxapp_1.csv", clear
  
    replace n_over_sub_j = ""  
    destring n_over_sub_j , replace 
    append using `ests'

    save `ests', replace 

    import delimited "/Volumes/lausd/decentralized_choice/estimates/counterfactuals/eta_out_K3_simulate_counterfactual_9_used_eq_pi_2013_maxapp_1.csv", clear
    replace amte_math = ""
    replace amte_ela = ""
    destring amte_math amte_ela, replace
    replace n_over_sub_j = ""  
    destring n_over_sub_j , replace 
    append using `ests'
    rename attned_rate_black attend_rate_black


    collapse (mean) app_rate offer_rate attend_rate ///
                    tot_* amte_* bw_gap_* hw_gap_* pov_gap_* y_* sd_y_* ///
                    attend_rate_* n_*,  by(preference_model)
					
	gen preference_model2 = preference_model 
    replace preference_model2 = 1 if preference_model==0 // baseline
    replace preference_model2 = 2 if preference_model==2 // information
    replace preference_model2 = 3 if preference_model==3
    replace preference_model2 = 4 if preference_model==4 // centralized ue 
    replace preference_model2 = 5 if preference_model==5 // no travel costs
    replace preference_model2 = 6 if preference_model==6 // centralized ue + no travel costs 
    replace preference_model2 = 7 if preference_model==8 // centralized_ue + info  
    replace preference_model2 = 8 if preference_model==9 // centralized_ue + info  + no travel costs
    replace preference_model2 = 9 if preference_model==7 // maximizing achievement 

    * Implied seats filled 
    gen seats_filled = attend_rate * `total_students'
    gen share_seats_filled = seats_filled/`total_seats'


    * Calculate changes in group achievement relative to baseline
    local groups "black white hisp asian pov nonpov"
    foreach g of local groups{
        gen temp = y_`g'_math if preference_model==0
        egen baseline_`g'_math = mean(temp)
        drop temp
        gen change_`g'_math = y_`g'_math - baseline_`g'_math
        drop baseline_`g'

        gen temp = y_`g'_ela if preference_model==0
        egen baseline_`g'_ela = mean(temp)
        drop temp
        gen change_`g'_ela = y_`g'_ela - baseline_`g'_ela
        drop baseline_`g'

 
    }
    gen temp = share_seats_filled if preference_model==0
    egen baseline_share_seats_filled = mean(temp)
    drop temp
    gen change_share_seats_filled = share_seats_filled - baseline_share_seats_filled
    drop baseline_share_seats_filled

    gen temp = app_rate if preference_model==0
    egen baseline_app_rate = mean(temp)
    drop temp
    gen change_app_rate = app_rate - baseline_app_rate
    drop baseline_app_rate

    gen temp = tot_math if preference_model==0
    egen baseline_tot_math = mean(temp)
    drop temp
    gen change_tot_math = tot_math - baseline_tot_math
    drop baseline_tot_math

    gen temp = tot_ela if preference_model==0
    egen baseline_tot_ela = mean(temp)
    drop temp
    gen change_tot_ela = tot_ela - baseline_tot_ela
    drop baseline_tot_ela


    * Calculate changes in ybar (using poverty shares from the main sample) 
    gen delta_ybar_math = change_pov_math * 0.6960 + change_nonpov_math * 0.3040
    gen delta_ybar_ela = change_pov_ela * 0.6960 + change_nonpov_ela * 0.3040

    * Make a table summarizing results 
    texdoc init "/Volumes/lausd/decentralized_choice/output/tables/counterfactual_table.tex" , replace force 
	tex \begin{tabular}{lcccccccccc} \toprule \hline \hline 
    tex & \multicolumn{6}{c}{Demand and Enrollment} & \multicolumn{4}{c}{Achievement} \\ \cline{2-11}
    tex & & & & & & & & & & \\
	tex  & Asian & Black & Hispanic & White & Seats Filled & Apply & \multicolumn{2}{c}{$\bar{Y}$} & \multicolumn{2}{c}{TOT} \\ \cline{8-11}
    tex & & & & & & & Math & ELA & Math & ELA  \\ \hline 
    tex & & & & & & & & & & \\
    local models "1 2 3 4 5 6 7 8 9 "
    foreach model of local models{
        if "`model'"=="1"{
            local model_name "Baseline"
        }
        if "`model'"=="2"{
            local model_name "General Information"
        }
        if "`model'"=="3"{
            local model_name "Targeted Information"
        }
        if "`model'"=="4"{
            local model_name "Centralized UE"
        }
        if "`model'"=="5"{
            local model_name "No Travel Costs"
        }
        if "`model'"=="6"{
            local model_name "Centralized UE + No Travel Costs"
        }
        if "`model'"=="7"{
            local model_name "Centralized UE + Info"
        }
        if "`model'"=="8"{
            local model_name "Centralized UE + Info + No Travel Costs"
        }
        if "`model'"=="9"{
            local model_name "Achievement Optimal"
        }
        sum attend_rate_black if preference_model2==`model'
        local attend_rate_black = r(mean)
        sum attend_rate_white if preference_model2==`model'
        local attend_rate_white = r(mean)
        sum attend_rate_hisp if preference_model2==`model'
        local attend_rate_hisp = r(mean)
        sum attend_rate_asian if preference_model2==`model'
        local attend_rate_asian = r(mean)
        sum change_share_seats_filled if preference_model2==`model'
        local share_seats_filled = r(mean)
        sum change_app_rate if preference_model2==`model'
        local app_rate = r(mean)
        sum delta_ybar_math if preference_model2==`model'
        local delta_ybar_math = r(mean)
        sum delta_ybar_ela if preference_model2==`model'
        local delta_ybar_ela = r(mean)
        sum change_tot_math if preference_model2==`model'
        local change_tot_math = r(mean)
        sum change_tot_ela if preference_model2==`model'
        local change_tot_ela = r(mean)

        * For baseline model, do not uses changes 
        if `model'==1{
            local delta_ybar_math = .
            local delta_ybar_ela = .
            sum tot_math if preference_model2==`model'
            local change_tot_math = r(mean)
            sum tot_ela if preference_model2==`model'
            local change_tot_ela = r(mean)

            sum share_seats_filled if preference_model2==`model'
            local share_seats_filled = r(mean)
            sum app_rate if preference_model2==`model'
            local app_rate = r(mean)

        }
        * add the panel label 
        if `model'==2{
            texdoc write & & & & & & & & & & \\
            texdoc write & \multicolumn{10}{c}{\textbf{ Change Relative to Baseline} } \\
            texdoc write & & & & & & & & & & \\
        }

        texdoc write `model_name' & ///
            `=string(`attend_rate_asian'*100, "%9.2f")' & ///
            `=string(`attend_rate_black'*100, "%9.2f")' & ///
            `=string(`attend_rate_hisp'*100, "%9.2f")' & ///
            `=string(`attend_rate_white'*100, "%9.2f")' & ///
            `=string(`share_seats_filled'*100, "%9.2f")' & ///
            `=string(`app_rate'*100, "%9.2f")' & ///
            `=string(`delta_ybar_math', "%9.3f")' & ///
            `=string(`delta_ybar_ela', "%9.3f")' & ///
            `=string(`change_tot_math', "%9.3f")' & ///
            `=string(`change_tot_ela', "%9.3f")' \\

    }
    texdoc write & & & & & & & & & & \\ \hline \hline 
    tex  \end{tabular}


    graph bar  change_black_math change_white_math change_hisp_math change_asian_math  if preference_model!=0 , ///
        over(preference_model2, relabel(1 "General Info." 2 "Targeted Info." 3 "Mandatory" 4 "Busing" 5 "Mandatory + Busing" 6 "Mandatory + Info"  7 "Mandatory + Info + Busing" 8 "Achievement Optimal") label(angle(45))) ///
        bar(1, color(black)) bar(2, fcolor(maroon) lcolor(black)) bar(3, fcolor(gs10) lcolor(black)) bar(4, lcolor(black) fcolor(none)) ///
        ytitle("Change in Student Achievement (SD)") legend(order(1 "Black" 2 "White" 3 "Hispanic" 4 "Asian") row(1) pos(6))
        local groups "black white hisp asian pov nonpov"
    graph export $figures/counterfactuals_change_group_math.pdf, replace

    graph bar  change_black_ela change_white_ela change_hisp_ela change_asian_ela  if preference_model!=0 , ///
        over(preference_model2, relabel(1 "General Info." 2 "Targeted Info." 3 "Mandatory" 4 "Busing" 5 "Mandatory + Busing" 6 "Mandatory + Info"  7 "Mandatory + Info + Busing" 8 "Achievement Optimal") label(angle(45))) ///
        bar(1, color(black)) bar(2, fcolor(maroon) lcolor(black)) bar(3, fcolor(gs10) lcolor(black)) bar(4, lcolor(black) fcolor(none)) ///
        ytitle("Change in Student Achievement (SD)") legend(order(1 "Black" 2 "White" 3 "Hispanic" 4 "Asian") row(1) pos(6))
        local groups "black white hisp asian pov nonpov"
    graph export $figures/counterfactuals_change_group_ela.pdf, replace


/* 
    graph bar  tot_math tot_ela  , ///
        over(preference_model2, relabel(1 "Baseline" 2 "General Information" 3 "Targeted Information" 4 "Centralized UE") ) ///
        bar(1, color(black)) bar(2, fcolor(maroon) lcolor(black))  ///
        ytitle("Treatment on the Treated") legend(order(1 "Math" 2 "ELA") row(1) pos(6))
    graph export $figures/counterfactuals_tot_math_ela.pdf, replace */


    *gen ybar_math = change_asian_math*.0368641  + change_white_math * .0731182  + change_black_math * .0959533  + change_hisp_math *.7662734 
    *gen ybar_ela = change_asian_ela*.0368641  + change_white_ela * .0731182  + change_black_ela * .0959533  + change_hisp_ela *.7662734



    graph bar  delta_ybar_math delta_ybar_ela  if preference_model2!=1 , ///
        over(preference_model2, relabel(1 "General Info." 2 "Targeted Info." 3 "Mandatory" 4 "Busing" 5 "Mandatory + Busing" 6 "Mandatory + Info"  7 "Mandatory + Info + Busing" 8 "Achievement Optimal") label(angle(45))) ///
        bar(1, color(black)) bar(2, fcolor(maroon) lcolor(black))  ///
        ytitle("Change in Aggregate Achievement (SD)") legend(order(1 "Math" 2 "ELA") row(1) pos(6))

 

    graph export $figures/counterfactuals_aggregate_achievement.pdf, replace


end 


capture program drop modelAccuracy 
program define modelAccuracy
    import delimited "/Volumes/lausd/decentralized_choice/estimates/eta_out_k3_group_stats2.csv", clear 
    rename (*_pred) (pred_*)
    rename (*_actual) (act_*)
    rename (n_offers_high_score_school_actua) (act_n_offers_high_score_school)



    reshape long pred_ act_ , i(program) string
    twoway (scatter act_ pred_ if program !=0 & _j!="n_enr", color(gs10)  ) (lfit act_ pred_ if program !=0) (function y= x if program !=0) (scatter act_ pred_ if _j=="n_enr" & program!=0, color(maroon)), legend(off) ytitle("Enrollment") xtitle("Model Predicted Enrollment")
    graph export /Volumes/lausd/decentralized_choice/output/figures/model_accuracy_enrollment.pdf, replace

    import delimited "/Volumes/lausd/decentralized_choice/estimates/eta_out_k3_group_stats.csv", clear 
    rename (*_pred) (pred_*)
    rename (*_actual) (act_*)
    *rename (n_offers_high_score_school_actua) (act_n_offers_high_score_school)
    drop if program==0


    reshape long pred_ act_ , i(program) string
    twoway (scatter act_ pred_ if program !=0 & _j!="n_apps", color(gs10)  ) (lfit act_ pred_ if program !=0) (function y= x if program !=0, color(black)) (scatter act_ pred_ if _j=="n_apps" & program!=0, color(maroon)), legend(off) ytitle("Applications") xtitle("Model Predicted Applications")
    graph export /Volumes/lausd/decentralized_choice/output/figures/model_accuracy_applications.pdf, replace
end

main