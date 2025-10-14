# DOCUMENTATION -----------------------------------------------------------

#   Author:     Jesse Bruhn
#   Contact:    jesse_bruhn@brown.edu

# PREAMBLE ----------------------------------------------------------------

#Clear workspace
rm(list=ls())

#Load Packages
pckgs <- c("tidyverse")

lapply(pckgs, library, character.only=TRUE)

#Set Options
options(scipen=100)

#set seed
set.seed(2025)

#Clean Workspace
rm(pckgs)

#Session Info for Reproducibility
Sys.Date()
sessionInfo()

#root folder
root_dir <- "/Volumes/lausd/decentralized_choice"

# NOTES -------------------------------------------------------------------



# SET TOGGLES -------------------------------------------------------------


#end of sample period
sample_endyear <- 2017

# Variables to check for balance
balance <- tribble(
  
  ~vars, ~var_lables,
  
  #default list
  "female", "Female",
  "app_stu_black","Black",
  "app_stu_latino", "Hispanic",
  "english_learner", "English Learner",
  "special_edu", "Special Education",
  "poverty", "Poverty",
  "college_grad", "Parent reports going to college",
  "eng_home", "Speaks english at home",
  "esp_home", "Speaks spanish at home", 
  # "app_nearest_mag_dist", "Distance to nearest choice school (miles)",
  "L1numsuspensions", "Baseline Suspensions",
  "L1math", "Baseline Math Score",
  "L1ela","Baseline ELA Score",
  "anyScore_ms", "Missing Baseline Score" 
)

#variables that we require to be non-missing
require_non_missing <- c(
  "app_dist_rel_quint",
  "above_cutoff",
  "lottery_ID" 
)

#variables to check for heterogeneity
dist_het_vars <- tibble::tribble(
  ~variable,          ~label,
  "lag_math",         "Achievement",
  "lag_suspensions",  "Suspensions",
  "female",           "Female",
  "poverty",          "Poverty",
  "special_edu",      "Special Education",
  "gifted",           "Gifted",
  "english_learner",  "English Learner",
  "eng_home",         "English at Home",
  "esp_home",         "Spanish at Home",
  "born_usa",         "Born USA"
)

# #variables to use for preference heterogeneity exercise
# pref_het_vars <- c(
#   "lag_math",
#   "lag_suspensions",
#   "female",
#   "poverty",
#   "special_edu",
#   "gifted",
#   "english_learner",
#   "eng_home",
#   "esp_home",
#   "born_usa"
# )
# 

# #distance or relative distance ("above_dist" or "rel_dist")
# distance_var <- "above_dist"
# 
# #Do we want the "N" on our figures as meta-data? 
# figure_N <- F
# 
# #Do we want a new pull of census block characteristics? 
# new_census_pull <- F
# 

#output folder name
figures <- "/Volumes/lausd/decentralized_choice/output/figures/"
tables  <- "/Volumes/lausd/decentralized_choice/output/tables/"



# HELPER FUNCTIONS --------------------------------------------------------

# robust cluster vcov spec shortcut
vcov_cluster <- function(cluster_var) stats::as.formula(paste0("~", cluster_var))

# Pull control mean (for control = above_cutoff == 0) with non-missing mask
control_mean <- function(df, var, mask) {
  df %>%
    # filter(!!mask, .data$above_cutoff == 0) %>%
    summarise(m = mean(.data[[var]], na.rm = TRUE)) %>%
    pull(m)
}

#helps with building formula's programatically
backtick <- function(x) paste0("`", x, "`")

# LOAD DATA ---------------------------------------------------------------

#lottery data
lotto <- haven::read_dta("/Volumes/lausd/decentralized_choice/data/lotteries.dta") %>%
  mutate(across(everything(), ~ .))

# #legacy code to view variable labels
# var_labels <- map_chr(lotto, ~ attr(.x, "label") %||% NA_character_)
# var_mapping <- tibble(variable = names(lotto), label = var_labels)


#better distance data
distance <- read_csv("/Volumes/lausd/decentralized_choice/rawdata/Geocoded_Student_Addresses_2003_2025.csv")  





# MERGE IN DISTANCE DATA --------------------------------------------------


#clean up distance data
distance_mod <- distance %>%
  mutate(
    
    studentpseudoid = Student_Pseudo_Id,
    # school_yr = School_Year,
    new_lat = Lat, 
    new_lon = Long
    
  ) %>%
  mutate(
    school_yr = str_split_fixed(School_Year, "-", 2)[, 2] |> as.numeric()
  ) %>%
  select(studentpseudoid, 
         school_yr, 
         new_lat, 
         new_lon
  ) %>%
  filter(complete.cases(.))

#merge to lottery data
df <- lotto %>%
  left_join(distance_mod)

#now caluclate app_dist
df <- df %>%
  mutate(
    app_lon_num = vctrs::vec_data(app_lon),  # drop label attrs/classes
    app_lat_num = vctrs::vec_data(app_lat),
    new_lon_num = vctrs::vec_data(new_lon),
    new_lat_num = vctrs::vec_data(new_lat),
    app_dist_2  = geosphere::distHaversine(
      cbind(app_lon_num, app_lat_num),   # (lon, lat)
      cbind(new_lon_num, new_lat_num)
    ) / 1000  # kilometers
  ) %>%
  select(-ends_with("_num"))



df <- df %>%
  select(-app_dist_rel_quint, -above_dist1, -above_dist2, -above_dist3, 
         -above_dist4, -above_dist5) %>%
  mutate(
    app_dist_rel_quint = ntile(app_dist_2, 5), 
    above_dist1 = ifelse(app_dist_rel_quint == 1 & above_cutoff==1, 1, 0),
    above_dist2 = ifelse(app_dist_rel_quint == 2 & above_cutoff==1, 1, 0),
    above_dist3 = ifelse(app_dist_rel_quint == 3 & above_cutoff==1, 1, 0),
    above_dist4 = ifelse(app_dist_rel_quint == 4 & above_cutoff==1, 1, 0),
    above_dist5 = ifelse(app_dist_rel_quint == 5 & above_cutoff==1, 1, 0)
  )

  

# APPLY SAMPLE RESTRICTIONS -----------------------------------------------

df <- df %>%
  filter(
    if_all(all_of(require_non_missing), ~!is.na(.)),
    endyear <= sample_endyear)

# LOTTERY BALANCE ---------------------------------------------------------


#function to build table
buildTable <- function(data, chars, controls = character()) {
  # helper to quote *simple* variable names safely in formulas
  backtick <- function(x) paste0("`", x, "`")
  
  ctrl_part <- if (length(controls) > 0) {
    paste0(" + ", paste(backtick(controls), collapse = " + "))
  } else {
    ""
  }
  
  rows <- purrr::map_dfr(chars, function(v) {
    
    fml <- as.formula(
      paste0(backtick(v), " ~ above_cutoff", ctrl_part, " | lottery_ID")
    )
    
    m <- fixest::feols(
      fml,
      data = df,
      vcov = fixest::vcov_cluster("lottery_ID")
    )
    
    tt <- broom::tidy(m)
    co <- dplyr::filter(tt, term == "above_cutoff") |> dplyr::slice(1)
    
    ctrl_mean <- control_mean(df, v)  # your existing helper
    
    tibble::tibble(
      variable     = v,
      control_mean = ctrl_mean,
      diff         = co$estimate,
      se           = co$std.error,
      p_value      = co$p.value,
      N            = stats::nobs(m)  # number of rows actually used by the model
    )
  })
  
  rows
}

#build balance table
balance_table <- buildTable(
  data = df,
  chars = balance$vars
)

balance_table <- balance_table %>%
  mutate(variable_label = balance$var_lables)

#joint test
#will operationalize by creating a stacked data frame
df_stack <- df %>%
  select(all_of(require_non_missing),
         all_of(balance$vars)
  ) %>%
  pivot_longer(
    cols = all_of(balance$vars),
    names_to = "stack",
    values_to = "outcome"
  ) %>%
  mutate(
    lottery_ID_by_stack = str_glue("{lottery_ID}:{stack}")
  )


fml <- as.formula("outcome ~ above_cutoff:stack | lottery_ID_by_stack")


joint_mod <- fixest::feols(fml,
                           data = df_stack,
                           vcov = fixest::vcov_cluster("lottery_ID")
)

all_coef_names <- names(coef(joint_mod))
test_coef_names <- all_coef_names[str_detect(all_coef_names, "above_cutoff:")]

pvalue_sur <- fixest::wald(joint_mod, keep = test_coef_names)$p


# Write LaTeX table
obs <- nrow(df)

# helper: map p-values to stars
p_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) "$^{***}$"
  else if (p < 0.05) "$^{**}$"
  else if (p < 0.1) "$^{*}$"
  else ""
}

# keep only the columns we need, in the order we print them
balance_table <- balance_table %>%
  dplyr::select(variable_label, control_mean, diff, se, p_value)

lines <- c(
  "\\begin{tabular}{lccc} \\toprule \\hline",
  " & Control Mean & Lottery Offer Difference & p-value \\\\ \\hline \\addlinespace"
)

# helper to format p-values nicely
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf('%.3f', p)
}

add_row <- function(vlab, mean, coef, p){
  str_glue(
    "{vlab} & {sprintf('%.3f', mean)} & {sprintf('%.3f', coef)}{p_stars(p)} & {fmt_p(p)} \\\\"
  )
}

# add SE line (blank p-value column)
add_se <- function(se) str_glue("&  & ({sprintf('%.3f', se)}) &  \\\\")

for (i in seq_len(nrow(balance_table))) {
  vlabel <- balance_table$variable_label[i]
  lines <- c(
    lines,
    add_row(
      vlab = vlabel,
      mean = balance_table$control_mean[i],
      coef = balance_table$diff[i],
      p    = balance_table$p_value[i]
    ),
    add_se(balance_table$se[i])
  )
}

lines <- c(
  lines,
  "\\addlinespace \\hline \\addlinespace",
  str_glue("Joint test & \\multicolumn{{3}}{{c}}{{{round(pvalue_sur, digits = 3)}}} \\\\"),
  str_glue("Observations & \\multicolumn{{3}}{{c}}{{{formatC(obs, format='d', big.mark=',')}}} \\\\ \\addlinespace \\hline\\bottomrule"),
  "\\end{tabular}"
)

# Now write the table
out_path <- str_glue("{tables}/lottery_balance.tex")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(paste(lines, collapse = "\n"), con = out_path)




# DISTANCE GRADIENT -------------------------------------------------------


# Plot quintile effects 
plot_quintile <- function(model,
                          label_y,
                          note1,
                          note2,
                          outpath,
                          include_N = FALSE,
                          table_outpath = NULL,
                          table_group_header = "Quintile of Distance to Choice School") {
  # ---- N for caption (plot only) ----
  N <- tryCatch(stats::nobs(model), error = function(e) NA_integer_)
  note2_aug <- if (isTRUE(include_N)) {
    paste0(note2, "\n", "(N = ", ifelse(is.na(N), "NA", format(N, big.mark=",")), ")")
  } else {
    note2
  }
  
  # ---- Tidy model with CI ----
  tt <- broom::tidy(model, conf.int = TRUE, conf.level = 0.95)
  qq <- tt %>%
    dplyr::filter(stringr::str_detect(term, "^above_.*[1-5]$")) %>%
    dplyr::mutate(quint = as.integer(stringr::str_extract(term, "[1-5]"))) %>%
    dplyr::arrange(quint)
  
  # Ensure we have exactly quintiles 1..5 and keep CI columns too
  df5 <- dplyr::tibble(quint = 1:5) %>%
    dplyr::left_join(
      qq[, c("quint", "estimate", "std.error", "conf.low", "conf.high")],
      by = "quint"
    )
  
  # ---- Plot ----
  whisker_half_width <- 0.02
  g <- ggplot2::ggplot(df5, ggplot2::aes(x = quint, y = estimate)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      width = 0
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = quint - whisker_half_width, xend = quint + whisker_half_width,
                   y = conf.high, yend = conf.high),
      linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = quint - whisker_half_width, xend = quint + whisker_half_width,
                   y = conf.low, yend = conf.low),
      linewidth = 0.4
    ) +
    ggplot2::scale_x_continuous(breaks = 1:5) +
    ggplot2::labs(
      x = "Distance Quintile",
      y = label_y,
      caption = paste0(note1, "\n", note2_aug)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor      = ggplot2::element_blank(),
      panel.grid.major.x    = ggplot2::element_blank(),
      axis.line             = ggplot2::element_line(linewidth = 0.5, colour = "#9aa0a6"),
      axis.ticks            = ggplot2::element_line(linewidth = 0.5, colour = "#b9bfc6"),
      axis.ticks.length     = grid::unit(3, "pt"),
      plot.caption.position = "plot",
      plot.caption          = ggplot2::element_text(hjust = 0, size = 9),
      legend.position       = "none"
    )
  
  if (!missing(outpath) && is.character(outpath) && nzchar(outpath)) {
    dir.create(dirname(outpath), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = outpath, plot = g, width = 12, height = 8, dpi = 300)
  }
  
  # ---- Write pure LaTeX tabular (no notes; ready to \input{}) ----
  f3 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))
  
  cell <- function(b, se) {
    b_  <- f3(b); se_ <- f3(se)
    if (!nzchar(b_) && !nzchar(se_)) return("")
    glue::glue("\\shortstack{{{b_}\\\\({se_})}}")
  }
  
  cells <- mapply(cell, df5$estimate, df5$std.error, USE.NAMES = FALSE)
  
  if (is.null(table_outpath) || !nzchar(table_outpath)) {
    if (!missing(outpath) && is.character(outpath) && nzchar(outpath)) {
      base <- sub("\\.[A-Za-z]+$", "", outpath)
      table_outpath <- paste0(base, "_quintile_tabular.tex")
    } else {
      table_outpath <- file.path(getwd(), "quintile_tabular.tex")
    }
  }
  dir.create(dirname(table_outpath), recursive = TRUE, showWarnings = FALSE)
  
  latex_lines <- c(
    "\\begin{tabular}{@{} l *{5}{c} @{}}",
    "\\toprule",
    glue::glue(" & \\multicolumn{{5}}{{c}}{{\\textit{{{table_group_header}}}}} \\\\"),
    "\\cmidrule(lr){2-6}",
    " & (1) & (2) & (3) & (4) & (5) \\\\",
    "\\midrule",
    glue::glue("Lottery Offer & {cells[1]} & {cells[2]} & {cells[3]} & {cells[4]} & {cells[5]} \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )
  
  writeLines(latex_lines, con = table_outpath)
  
  structure(g, N = N, table_outpath = table_outpath)
}



fml_rhs <- paste0( " ~ 0 + above_dist1 + above_dist2 + above_dist3 + above_dist4 + above_dist5 +
                            i(app_dist_rel_quint) ",  " | lottery_ID")  

joint_q5eqq1 <- c(
  "above_dist1 = above_dist2",
  "above_dist1 = above_dist3",
  "above_dist1 = above_dist4",
  "above_dist1 = above_dist5"
)

joint_all0 <- c("above_dist1", "above_dist2", "above_dist3", "above_dist4", "above_dist5")




###
# FS Distance Gradient
###

#create formula
fml <- as.formula(
  paste0("enrolled_app", fml_rhs)
)

#Estimate distance gradient
fs_gradient <- fixest::feols(fml,
                             data = df,
                             vcov = vcov_cluster("lottery_ID"))



p_q5eqq1 <- car::linearHypothesis(fs_gradient,
                                  joint_q5eqq1,
                                  vcov. = vcov(fs_gradient),
                                  test = "F"
                                  
)$`Pr(>F)`[2]

p_all_0 <- fixest::wald(fs_gradient, joint_all0)$p

plot_quintile(
  model = fs_gradient, 
  label_y = "First Stage Effect",
  note1 = str_glue("H0 All Equal p-value: {sprintf('%.3f', p_q5eqq1)}"),
  note2 = str_glue("H0 Jointly Zero p-value: {sprintf('%.3f', p_all_0)}"),
  outpath = str_glue("{figures}/fs_dist_gradient.pdf"), 
)


###
# RF Distance Gradient -- Math
###

#create formula
fml <- as.formula(
  paste0("avg_Fmath", fml_rhs)
)

#Estimate distance gradient
rf_math_gradient <- fixest::feols(fml,
                                  data = df,
                                  vcov = vcov_cluster("lottery_ID"))



p_q5eqq1_m <- car::linearHypothesis(rf_math_gradient,
                                  joint_q5eqq1,
                                  vcov. = vcov(rf_math_gradient),
                                  test = "F"
                                  
)$`Pr(>F)`[2]

p_all_0_m <- fixest::wald(rf_math_gradient, joint_all0)$p

plot_quintile(
  model = rf_math_gradient, 
  label_y = "Math (Reduced Form)",
  note1 = str_glue("H0 All Equal p-value: {sprintf('%.3f', p_q5eqq1_m)}"),
  note2 = str_glue("H0 Jointly Zero p-value: {sprintf('%.3f', p_all_0_m)}"),
  outpath = str_glue("{figures}/math_dist_gradient.pdf")
)


###
# RF Distance Gradient -- ELA
###

#create formula
fml <- as.formula(
  paste0("avg_Fela", fml_rhs)
)

#Estimate distance gradient
rf_ela_gradient <- fixest::feols(fml,
                                 data = df,
                                 vcov = vcov_cluster("lottery_ID"))



p_q5eqq1_e <- car::linearHypothesis(rf_ela_gradient,
                                  joint_q5eqq1,
                                  vcov. = vcov(rf_ela_gradient),
                                  test = "F"
                                  
)$`Pr(>F)`[2]

p_all_0_e <- fixest::wald(rf_ela_gradient, joint_all0)$p

plot_quintile(
  model = rf_ela_gradient, 
  label_y = "ELA (Reduced Form)",
  note1 = str_glue("H0 All Equal p-value: {sprintf('%.3f', p_q5eqq1_e)}"),
  note2 = str_glue("H0 Jointly Zero p-value: {sprintf('%.3f', p_all_0_e)}"),
  outpath = str_glue("{figures}/ela_dist_gradient.pdf")
)



# DISTANCE HETEROGENEITY --------------------------------------------------

#function to estimate CI on difference
coef_diff_ci <- function(m,
                         term_hi = "above_dist1",
                         term_lo = "above_dist5",
                         level = 0.95,
                         df_method = c("auto", "resid", "normal")) {
  
  b  <- stats::coef(m)
  V  <- stats::vcov(m)
  
  if (!(term_hi %in% names(b) && term_lo %in% names(b))) {
    stop("One or both terms not found in model coefficients.")
  }
  
  est <- b[[term_hi]] - b[[term_lo]]
  se  <- sqrt(V[term_hi, term_hi] + V[term_lo, term_lo] - 2 * V[term_hi, term_lo])
  
  alpha <- 1 - level
  
  get_crit <- function() {
    if (df_method == "normal") {
      return(list(crit = stats::qnorm(1 - alpha/2), df = NA_real_, dist = "normal"))
    }
    if (df_method == "resid") {
      df_res <- suppressWarnings(tryCatch(stats::df.residual(m), error = function(e) NA_real_))
      if (is.na(df_res)) {
        # fallback to normal if we can't get residual df
        return(list(crit = stats::qnorm(1 - alpha/2), df = NA_real_, dist = "normal"))
      } else {
        return(list(crit = stats::qt(1 - alpha/2, df = df_res), df = df_res, dist = "t(resid)"))
      }
    }
    # df_method == "auto": try to borrow df from fixest::wald() for the equality test
    w <- suppressWarnings(tryCatch(
      fixest::wald(m, c(term_hi)),
      error = function(e) NULL
    ))
    # If w exists and contains a finite df, use t; else normal
    df_auto <- tryCatch(w$df2, error = function(e) NA_real_)
    if (!is.null(w) && is.finite(df_auto)) {
      list(crit = stats::qt(1 - alpha/2, df = df_auto), df = df_auto, dist = "t(fixest)")
    } else {
      list(crit = stats::qnorm(1 - alpha/2), df = NA_real_, dist = "normal")
    }
  }
  
  crit_info <- get_crit()
  ci <- c(est - crit_info$crit * se, est + crit_info$crit * se)
  
  # two-sided p-value using the same distribution choice
  tstat <- est / se
  p_val <- if (identical(crit_info$dist, "normal")) {
    2 * stats::pnorm(abs(tstat), lower.tail = FALSE)
  } else {
    2 * stats::pt(abs(tstat), df = crit_info$df, lower.tail = FALSE)
  }
  
  tibble::tibble(
    contrast  = paste0(term_hi, " - ", term_lo),
    estimate  = est,
    se        = se,
    conf.low  = ci[1],
    conf.high = ci[2],
    p.value   = p_val,
    df_used   = crit_info$df,
    dist_used = crit_info$dist
    # also return fixest's equality-test p-value for reference
    
  )
}

# function to help build interaction formulas
augment_rhs <- function(fml_rhs, var, mode = c("all", "cutoff")) {
  mode <- match.arg(mode)
  
  # Split RHS into the part before and after the "|" (fixest FE separator)
  parts <- strsplit(fml_rhs, "\\|")[[1]]
  pre   <- trimws(parts[1])                 # e.g. "~ 0 + above_rel_dist1 + ... + i(...) + `L1math`"
  post  <- if (length(parts) > 1) paste0(" |", parts[2]) else ""
  
  # Helper: does token already appear as a main effect?
  has_term <- function(txt, term) {
    pattern <- paste0("(^|[+~\\s\\(])`?", term, "`?(?=\\b|\\s)")
    str_detect(txt, pattern)
  }
  
  additions <- character(0)
  
  # 1) main effect of var (if not already present)
  if (!has_term(pre, var)) additions <- c(additions, var)
  
  # 2) interactions
  if (mode == "all") {
    # Find unique "above_*" variable names in the regressors part
    above_terms <- unique(str_extract_all(pre, "\\babove_[A-Za-z0-9_\\.]+\\b")[[1]])
    if (length(above_terms)) {
      inter_terms <- paste(var, above_terms, sep = ":")
      # drop any that already exist verbatim
      keep <- !str_detect(pre, paste0("\\b", str_replace_all(inter_terms, "([:])", "\\\\$1"), "\\b"))
      additions <- c(additions, inter_terms[keep])
    }
  } else if (mode == "cutoff") {
    cutoff <- "above_cutoff"
    inter  <- paste(var, cutoff, sep = ":")
    if (!str_detect(pre, paste0("\\b", var, "\\s*:\\s*", cutoff, "\\b")))
      additions <- c(additions, inter)
  }
  
  # Glue everything back together (before the "|")
  new_pre <- if (length(additions)) paste(pre, paste(additions, collapse = " + "), sep = " + ") else pre
  paste0(new_pre, post)
}



###
# Heterogeneity for math
###

math_het_results <- coef_diff_ci(m = rf_math_gradient,
                                 level = 0.95,
                                 df_method = "auto"
) %>%
  mutate(term = "baseline", 
         label = "Baseline")

for(var in dist_het_vars$variable){
  print(var)
  lab <- dist_het_vars %>%
    filter(variable == var) %>% 
    pull(label)
  
  het_rhs <- augment_rhs(fml_rhs, var, mode = "cutoff")
  fml <- as.formula(paste0("avg_Fmath ", het_rhs))
  het_reg <- fixest::feols(fml,
                           data = df,
                           vcov = vcov_cluster("lottery_ID"))
  
  new_het_results <- coef_diff_ci(m = het_reg,
                                  level = 0.95,
                                  df_method = "auto"
  ) %>%
    mutate(term = var, 
           label = lab
    )
  
  math_het_results <- bind_rows(math_het_results, new_het_results)
  
}

#make plot
math_het_results <- math_het_results %>%
  mutate(
    label_fctr = factor(label, levels = label, labels = label)
  )

g_math <- ggplot(math_het_results, aes(x = label_fctr, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(y = "Q1 - Q5 Difference", x = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())


ggsave(str_glue("{figures}/math_dist_het.pdf"), g_math, width = 10, height = 5)


###
# Heterogeneity for ela
###

ela_het_results <- coef_diff_ci(m = rf_ela_gradient,
                                level = 0.95,
                                df_method = "auto"
) %>%
  mutate(term = "baseline", 
         label = "Baseline")

for(var in dist_het_vars$variable){
  print(var)
  lab <- dist_het_vars %>%
    filter(variable == var) %>% 
    pull(label)
  
  het_rhs <- augment_rhs(fml_rhs, var, mode = "cutoff")
  fml <- as.formula(paste0("avg_Fela ", het_rhs))
  het_reg <- fixest::feols(fml,
                           data = df,
                           vcov = vcov_cluster("lottery_ID"))
  
  new_het_results <- coef_diff_ci(m = het_reg,
                                  level = 0.95,
                                  df_method = "auto"
  ) %>%
    mutate(term = var, 
           label = lab
    )
  
  ela_het_results <- bind_rows(ela_het_results, new_het_results)
  
}

#make plot
ela_het_results <- ela_het_results %>%
  mutate(
    label_fctr = factor(label, levels = label, labels = label)
  )

g_ela <- ggplot(ela_het_results, aes(x = label_fctr, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(y = "Q1 - Q5 Difference", x = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())


ggsave(str_glue("{figures}/ela_dist_het.pdf"), g_ela, width = 10, height = 5)





# PREFERENCE HETEROGENEITY ------------------------------------------------



# Compute app-level, lottery-level, school-year (magnet) means and squared diffs
augment_prefs <- function(df, v) {
  df %>%
    group_by(lottery_ID) %>% mutate(app_l = mean(.data[[v]], na.rm = TRUE)) %>% ungroup() %>%
    group_by(mag_cd, school_yr) %>% mutate(app_j = mean(.data[[v]], na.rm = TRUE)) %>% ungroup() %>%
    mutate(app_ = mean(.data[[v]], na.rm = TRUE)) %>%
    mutate(
      !!str_glue("stu_app_{v}") := (.data[[v]] - .data$app_)^2,
      !!str_glue("stu_j_{v}")   := (.data[[v]] - .data$app_j)^2,
      !!str_glue("stu_l_{v}")   := (.data[[v]] - .data$app_l)^2
    ) %>%
    select(-app_, -app_j, -app_l)
}

df_aug <- reduce(dist_het_vars$variable, ~ augment_prefs(.x, .y), .init = df)

# Norms (row totals of squared diffs)
df_aug <- df_aug %>%
  rowwise() %>%
  mutate(
    app_norm   = sqrt(sum(c_across(starts_with("stu_app_")), na.rm = TRUE)),
    app_j_norm = sqrt(sum(c_across(starts_with("stu_j_")),   na.rm = TRUE)),
    app_l_norm = sqrt(sum(c_across(starts_with("stu_l_")),   na.rm = TRUE))
  ) %>%
  ungroup()

# Quintiles by endyear (Stata xtile n(5) by(endyear))
df_aug <- df_aug %>%
  group_by(endyear) %>%
  mutate(norm_quart = ntile(app_j_norm, 5L)) %>%
  ungroup() %>%
  mutate(
    above_norm1 = as.integer(above_cutoff * as.numeric(norm_quart == 1)),
    above_norm2 = as.integer(above_cutoff * as.numeric(norm_quart == 2)),
    above_norm3 = as.integer(above_cutoff * as.numeric(norm_quart == 3)),
    above_norm4 = as.integer(above_cutoff * as.numeric(norm_quart == 4)),
    above_norm5 = as.integer(above_cutoff * as.numeric(norm_quart == 5))
  )

#create regression RHS and strings for test
fml_rhs <- paste0( " ~ 0 + above_norm1 + above_norm2 + above_norm3 + above_norm4 + above_norm5 +
                            i(norm_quart) ",  " | lottery_ID")  

joint_q5eqq1 <- c(
  "above_norm1 = above_norm2",
  "above_norm1 = above_norm3",
  "above_norm1 = above_norm4",
  "above_norm1 = above_norm5"
)

joint_all0 <- c("above_norm1", "above_norm2", "above_norm3", "above_norm4", "above_norm5")




###
# FS Preference Gradient
###

#create formula
fml <- as.formula(
  paste0("enrolled_app", fml_rhs)
)

#Estimate distance gradient
fs_pref_gradient <- fixest::feols(fml,
                             data = df_aug,
                             vcov = vcov_cluster("lottery_ID"))



p_q5eqq1 <- car::linearHypothesis(fs_pref_gradient,
                                  joint_q5eqq1,
                                  vcov. = vcov(fs_pref_gradient),
                                  test = "F"
                                  
)$`Pr(>F)`[2]

p_all_0 <- fixest::wald(fs_pref_gradient, joint_all0)$p

plot_quintile(
  model = fs_pref_gradient, 
  label_y = "First Stage Effect",
  note1 = str_glue("H0 All Equal p-value: {sprintf('%.3f', p_q5eqq1)}"),
  note2 = str_glue("H0 Jointly Zero p-value: {sprintf('%.3f', p_all_0)}"),
  outpath = str_glue("{figures}/fs_pref_gradient.pdf"), 
  table_group_header = "Quintile of Preference Index"
)


###
# RF Pref Gradient -- Math
###

#create formula
fml <- as.formula(
  paste0("avg_Fmath", fml_rhs)
)

#Estimate distance gradient
rf_math_pref_gradient <- fixest::feols(fml,
                                  data = df_aug,
                                  vcov = vcov_cluster("lottery_ID"))



p_q5eqq1 <- car::linearHypothesis(rf_math_pref_gradient,
                                  joint_q5eqq1,
                                  vcov. = vcov(rf_math_pref_gradient),
                                  test = "F"
                                  
)$`Pr(>F)`[2]

p_all_0 <- fixest::wald(rf_math_pref_gradient, joint_all0)$p

plot_quintile(
  model = rf_math_pref_gradient, 
  label_y = "Math (Reduced Form)",
  note1 = str_glue("H0 All Equal p-value: {sprintf('%.3f', p_q5eqq1)}"),
  note2 = str_glue("H0 Jointly Zero p-value: {sprintf('%.3f', p_all_0)}"),
  outpath = str_glue("{figures}/math_pref_gradient.pdf"), 
  table_group_header = "Quintile of Preference Index"
)


###
# RF Pref Gradient -- ELA
###

#create formula
fml <- as.formula(
  paste0("avg_Fela", fml_rhs)
)

#Estimate distance gradient
rf_ela_pref_gradient <- fixest::feols(fml,
                                 data = df_aug,
                                 vcov = vcov_cluster("lottery_ID"))



p_q5eqq1 <- car::linearHypothesis(rf_ela_pref_gradient,
                                  joint_q5eqq1,
                                  vcov. = vcov(rf_ela_pref_gradient),
                                  test = "F"
                                  
)$`Pr(>F)`[2]

p_all_0 <- fixest::wald(rf_ela_pref_gradient, joint_all0)$p

plot_quintile(
  model = rf_ela_pref_gradient, 
  label_y = "ELA (Reduced Form)",
  note1 = str_glue("H0 All Equal p-value: {sprintf('%.3f', p_q5eqq1)}"),
  note2 = str_glue("H0 Jointly Zero p-value: {sprintf('%.3f', p_all_0)}"),
  outpath = str_glue("{figures}/ela_pref_gradient.pdf"), 
  table_group_header = "Quintile of Preference Index"
)



