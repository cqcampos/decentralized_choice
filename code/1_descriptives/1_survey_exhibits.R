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

#Clean Workspace
rm(pckgs)

#Session Info for Reproducibility
Sys.Date()
sessionInfo()

#root folder
root_dir <- "/Volumes/lausd/decentralized_choice"

# NOTES -------------------------------------------------------------------



# READ IN DATA ------------------------------------------------------------


# original survey data -- 250 school sample
original_survey <- read_csv(str_glue("{root_dir}/rawdata/School Choice Data - 250 school sample.csv"))

#Chris's CCD directory data he uses for magnet definition
ccd_dir_86_20 <- haven::read_dta(str_glue("/Volumes/lausd/magnet/data/raw/ccd_directory_1986_2020.dta"))

# NOW CLEAN UP SOME STUFF IN SURVEY DATA ----------------------------------

#restrict to top 250 
df_survey <- original_survey %>% 
  filter(rank <= 152)

# CLEAN CCD DATA ----------------------------------------------------------


#drop school-year combos that are closed or that happen before 1999 or that 
#are not regular schools, and create magnet defintion which is: 
#(1) If listed as a magnet, you are a magnet
#(2) If you have magnet in your name,  you are a magnet. 
#(3) If you are listed as a magnet in more than 20% of observed years, we call you 
#    a magnet in all years after the "first" year you are listed as a magnet. 
df_ccd <- ccd_dir_86_20 %>%
  filter(year>=1998, 
         school_status != 2, 
         school_status != 6, 
         school_status != 7, 
         school_type == 1
         ) %>%
  group_by(ncessch) %>%
  mutate(
    magnet_dum = ifelse(!is.na(magnet) & magnet == 1, 1, 0), 
    mag_year = ifelse(magnet_dum ==1 , 
                      year, 
                      NA
                      ), 
    first_mag_year = min(mag_year, na.rm = T), 
    magnet_dum_improved = ifelse(!is.na(first_mag_year) & year >= first_mag_year, 1, 0), 
    # ever_mag = ifelse(any(magnet_dum == 1), 1, 0),
    mag_share = mean(magnet_dum == 1, na.rm = T), 
    mag_name = str_to_lower(school_name), 
    mag_in_name = str_detect(school_name, "magnet"),
    magnet_schl = ifelse((magnet_dum_improved==1 & mag_share >= 0.2 ) | isTRUE(mag_in_name), 1, 0)
    
  )


#NYC has tons of LEAIDs so need to recode them in NCES data to match the survey LEAID
nyc_codes <- c(
  "3600075", # NYC Alternative HS District
  "3600076", "3600077", "3600078", "3600079",
  "3600081", "3600083", "3600084", "3600085", "3600086",
  "3600087", "3600088", "3600090", "3600091", "3600092",
  "3600094", "3600095", "3600096", "3600097", "3600098",
  "3600099", "3600100", "3600101", "3600102", "3600103",
  "3600119", "3600120", "3600121", "3600122", "3600123",
  "3600135", # District 75 (city-wide special schools)
  "3600151", "3600152", "3600153"
)

#Also, EPIC in oklahoma has two LEAID's that were merged into 1, so need to recode
#that as well
epic_codes <- c("4000777", "4000800")

#clean some other useful variables
df_ccd <- df_ccd %>%
  ungroup() %>% 
  mutate(charter_schl = ifelse(!is.na(charter) & charter == 1, 1, 0), 
         free_lunch = ifelse(free_lunch < 0, NA, free_lunch),
         reduced_price_lunch = ifelse(reduced_price_lunch < 0, NA, reduced_price_lunch),
         free_or_reduced_price_lunch = ifelse(free_or_reduced_price_lunch < 0, NA, free_or_reduced_price_lunch),
         enrollment = ifelse(enrollment < 0, NA, enrollment),
         
         leaid = case_when(
                      leaid %in% nyc_codes ~ "3600076",
                      leaid %in% epic_codes ~ "4000777",
                      T ~ leaid
                    ),

         )

#Now fix charter schools
#strategy will be to grab every zip code contained in an LEAID, then find every
#charter school within those zip codes, and call that a charter that serves that LEAID.
non_charters <- df_ccd %>%
  filter(charter_schl == 0) %>%
  select(
    leaid, year, zip_location
  ) %>%
  filter(str_detect(zip_location, "^\\d{5}$"), #drop odd zip codes
         complete.cases(.)
         ) %>%
  unique(
  ) 
  

charters <- df_ccd %>%
  filter(charter_schl == 1, 
         !(leaid %in% unique(df_survey$`NCES ID`)), #removes charters where LEAID already puts them in one of the 120 biggest districts. 
         !(leaid %in% unique(non_charters$leaid))
         ) %>%
  select(year, zip_location, ncessch) %>%
  filter(str_detect(zip_location, "^\\d{5}$"), #drop odd zip codes
         complete.cases(.))

improved_charters <- charters %>%
  left_join(non_charters) %>%
  filter(!is.na(leaid)) %>%
  unique() %>%
  select(year, ncessch, leaid)

#after all this, there are 127 charter schools that serve more than one of the 
#LEAIDs in our survey in at least 1 year. For now, I'm going to decide which 
#district they serve based on the one they serve most frequently and after that
#by tiebreaking randomly
unique_sch_district_pairs <- improved_charters %>%
  select(leaid, ncessch) %>%
  group_by(ncessch) %>%
  count(leaid, name = "freq") %>%
  slice_max(freq, n = 1, with_ties = F) %>%
  select(-freq)

improved_charters <- improved_charters %>%
  semi_join(unique_sch_district_pairs) %>%
  mutate(new_charter_leaid = leaid) %>%
  select(-leaid)

df_ccd <- df_ccd %>%
  left_join(improved_charters) %>%
  mutate(
    dist_charter_schl = ifelse(charter_schl==1 & is.na(new_charter_leaid), 1, 0), 
    non_dist_charter_schl = ifelse(charter_schl == 1 & dist_charter_schl == 0, 1, 0),
    leaid = ifelse(is.na(new_charter_leaid), leaid, new_charter_leaid)
  ) 

#now summarise relevant variables to the district level

mode_chr <- function(x) names(which.max(table(x)))

df_ccd_districts <- df_ccd %>%
  arrange(leaid, year) %>%
  group_by(leaid, year) %>%
  summarise(
    lea_name = mode_chr(lea_name), 
    N = n(), 
    N_schl = n_distinct(ncessch), 
    N_magnet_schl = sum(magnet_schl, na.rm = T), 
    N_charter_schl = sum(charter_schl, na.rm = T),
    
    #NCES Enrollment
    NCES_free_lunch = sum(free_lunch*(1-non_dist_charter_schl), na.rm = T), 
    NCES_reduced_lunch = sum(reduced_price_lunch*(1-non_dist_charter_schl), na.rm = T), 
    NCES_frpl = sum(free_or_reduced_price_lunch*(1-non_dist_charter_schl), na.rm = T), 
    NCES_enrollment = sum(enrollment*(1-non_dist_charter_schl), na.rm = T), 
    
    #total charter enrollment
    charter_free_lunch = sum(charter_schl*free_lunch, na.rm = T), 
    charter_reduced_lunch = sum(charter_schl*reduced_price_lunch, na.rm = T), 
    charter_frpl = sum(charter_schl*free_or_reduced_price_lunch, na.rm = T), 
    charter_enrollment = sum(charter_schl*enrollment, na.rm = T), 
    
    #district charter enrollment
    dist_charter_free_lunch = sum(dist_charter_schl*free_lunch, na.rm = T), 
    dist_charter_reduced_lunch = sum(dist_charter_schl*reduced_price_lunch, na.rm = T), 
    dist_charter_frpl = sum(dist_charter_schl*free_or_reduced_price_lunch, na.rm = T), 
    dist_charter_enrollment = sum(dist_charter_schl*enrollment, na.rm = T),
    
    #nondistrict charter enrollment
    non_dist_charter_free_lunch = sum(non_dist_charter_schl*free_lunch, na.rm = T), 
    non_dist_charter_reduced_lunch = sum(non_dist_charter_schl*reduced_price_lunch, na.rm = T), 
    non_dist_charter_frpl = sum(non_dist_charter_schl*free_or_reduced_price_lunch, na.rm = T), 
    non_dist_charter_enrollment = sum(non_dist_charter_schl*enrollment, na.rm = T),
    
    #magnet magnet enrollment
    magnet_free_lunch = sum(magnet_schl*free_lunch, na.rm = T), 
    magnet_reduced_lunch = sum(magnet_schl*reduced_price_lunch, na.rm = T), 
    magnet_frpl = sum(magnet_schl*free_or_reduced_price_lunch, na.rm = T), 
    magnet_enrollment = sum(magnet_schl*enrollment, na.rm = T), 
    
    #district wide enrollment
    total_free_lunch = sum(free_lunch, na.rm = T), 
    total_reduced_lunch = sum(reduced_price_lunch, na.rm = T), 
    total_frpl = sum(free_or_reduced_price_lunch, na.rm = T), 
    total_enrollment = sum(enrollment, na.rm = T), 
    
    #TPS enrollment
    tps_free_lunch = total_free_lunch - charter_free_lunch - magnet_free_lunch, 
    tps_reduced_lunch = total_reduced_lunch - charter_reduced_lunch - magnet_reduced_lunch, 
    tps_frpl = total_frpl - charter_frpl - magnet_frpl, 
    tps_enrollment = total_enrollment - charter_enrollment - magnet_enrollment 
    
  )


# NOW MERGE EVERYTHING ----------------------------------------------------

df <- df_survey %>%
  rename(leaid = `NCES ID`) %>%
  inner_join(df_ccd_districts) 

#I do a bunch of spot checking on ones with a big discrepency between 
#my calculated enrollment in 2019 from CCD and the ones chris downloaded from an official CCD pub
#and that are included in the survey spreadsheet
# tmp <- df %>% 
#   filter(year == 2019) %>%
#   select(rank, leaid, state, district, enrollment, NCES_enrollment) %>%
#   mutate(
#     delta = enrollment - NCES_enrollment, 
#     delta_pct = abs(delta/enrollment)
#     
#   )
#looking at basically every district that has more than a 15p.p. change. 
#and also checking the top 10-20 districts in terms of overall size. 

#(1) Prosper, ISD has a huge delta, but this seems like a real decline
#(2) Epic charters had a huge delta, but this was due to having multiple LEAIDs in the 
#common core data which I've fixed.
#(3) Idea public schools looks like a real change. 
#(4) State Sponsored Charter schools looks like a real change. 
#(5) NEBO district looks like a real change. 
#(6) Lamar CISD looks like a real change. 
#(7) Northwest ISD looks like a real change.
#(8) Santa Ana looks like a real change. 
#(9) Brownsville looks like a real change. 
#(10) Kipp texas looks like a real change. 
#(11) Newark looks like a real change. 
#(12) LA looks like a real change
#(13) Houston ISD looks like a real change. 
#(14) Dallas ISD looks like a real change. 

#Find districts that only consist of charter schools
all_charter_districts <- df %>%
  filter(year == 2019) %>%
  mutate(charter_school_share = N_charter_schl/N_schl, 
         charter_student_share = dist_charter_enrollment/NCES_enrollment
         ) %>%
  select(leaid, district, charter_school_share, charter_student_share) %>%
  arrange(desc(charter_school_share)) %>%
  filter(charter_school_share == 1) %>%
  pull(leaid)

#apply some sample restrictions
df <- df %>%
  filter((!leaid %in% all_charter_districts)) #drops districts that are 100% charter. 

#now find districts that don't report FRPL
no_frpl_districts <- df %>%
  filter(year == 2019) %>% 
  mutate(frpl_pct = NCES_frpl/NCES_enrollment) %>%
  select(rank, leaid, district, frpl_pct) %>%
  arrange(frpl_pct) %>%
  filter(frpl_pct <0.1) %>%
  pull(leaid)

#now find places that report 0 enrollment
no_enrollment_districts <- df %>%
  filter(NCES_enrollment == 0) %>%
  pull(leaid)

# ENROLLMENT IN CHARTER V MAGNETS -----------------------------------------


#using national data (share all students)
time_series <- df_ccd_districts %>%
  group_by(year) %>%
  summarise(
    across(N_schl:tps_enrollment, sum)
    
  ) %>%
  filter(year %in% 1999:2019)


df_long <- time_series %>%
  transmute(
    year,
    `Charter` = charter_enrollment / total_enrollment,
     `District Option` = magnet_enrollment  / total_enrollment
  ) %>%
  pivot_longer(-year, names_to = "school_type", values_to = "share")

x_pad <- diff(range(df_long$year, na.rm = TRUE)) * 0.03

first_pts <- df_long %>%
  group_by(school_type) %>%
  filter(!is.na(share)) %>%
  slice_min(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    first_label = paste0(school_type, "  ",
                         scales::percent(share, accuracy = 0.1)),
    x_first = year - x_pad
  )

last_pts <- df_long %>%
  group_by(school_type) %>%
  filter(!is.na(share)) %>%
  slice_max(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    last_label = scales::percent(share, accuracy = 0.1),
    x_last = year + x_pad
  )


share_plots <-
  ggplot(df_long, aes(x = year, y = share, color = school_type)) +
  geom_line(linewidth = 1.1) +
  geom_point(data = first_pts, size = 1.6, show.legend = FALSE) +
  geom_point(data = last_pts,  size = 1.6, show.legend = FALSE) +
  geom_text(
    data = first_pts ,
    aes(x = x_first, label = first_label),
    hjust = 0.7, vjust = -1.3, show.legend = FALSE
  ) +
  geom_text(
    data = last_pts,
    aes(x = x_last, label = last_label),
    hjust = 0, vjust = 0.5, show.legend = FALSE
  ) +
  # === Axes back in: percent y, pretty breaks, small cushion ===
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.08)), 
    limits = c(0, NA)
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  scale_color_manual(values = c("Charter" = "#800000", "District Option" = "#000000")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "none",
    # Keep the panel clean; add only faint horizontal reference lines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.5),
    
    # AXES: thin lines + subtle ticks that fit the minimal look
    axis.line        = element_line(color = "#9aa0a6", linewidth = 0.5),
    axis.ticks       = element_line(color = "#b9bfc6", linewidth = 0.5),
    axis.ticks.length = unit(3, "pt"),
    
    # Axis text styling
    axis.text.x = element_text(size = 14, color = "#444444", margin = margin(t = 4)),
    axis.text.y = element_text(size = 12, color = "#444444", margin = margin(r = 6)),
    
    # Room on both sides for direct labels
    plot.margin       = margin(5.5, 40, 5.5, 40)
  )


ggsave(filename = str_glue("{root_dir}/output/figures/national_magnet_share_ts.pdf"), width = 14, height = 8)

#using 150 largest districts (share of total enrollment)
time_series <- df %>%
  group_by(year) %>%
  summarise(
    across(N_schl:tps_enrollment, sum)
    
  ) %>%
  filter(year %in% 1999:2019)

df_long <- time_series %>%
  transmute(
    year,
    Charter = charter_enrollment / total_enrollment,
    `District Option`  = magnet_enrollment  / total_enrollment
  ) %>%
  pivot_longer(-year, names_to = "school_type", values_to = "share")

x_pad <- diff(range(df_long$year, na.rm = TRUE)) * 0.03

first_pts <- df_long %>%
  group_by(school_type) %>%
  filter(!is.na(share)) %>%
  slice_min(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    first_label = paste0(school_type, "  ",
                         scales::percent(share, accuracy = 0.1)),
    x_first = year - x_pad
  )

last_pts <- df_long %>%
  group_by(school_type) %>%
  filter(!is.na(share)) %>%
  slice_max(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    last_label = scales::percent(share, accuracy = 0.1),
    x_last = year + x_pad
  )

share_plots <-
  ggplot(df_long, aes(x = year, y = share, color = school_type)) +
  geom_line(linewidth = 1.1) +
  geom_point(data = first_pts, size = 1.6, show.legend = FALSE) +
  geom_point(data = last_pts,  size = 1.6, show.legend = FALSE) +
  geom_text(
    data = first_pts,
    aes(x = x_first, label = first_label),
    hjust = 0.7, vjust = -1.3, show.legend = FALSE
  ) +
  geom_text(
    data = last_pts,
    aes(x = x_last, label = last_label),
    hjust = 0, vjust = 0.5, show.legend = FALSE
  ) +
  # Axes: percent y, start at 0, pretty breaks, small cushion
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.08)),
    limits = c(0, NA)
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  scale_color_manual(values = c("Charter" = "#800000", "District Option" = "#000000")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    # Grid: faint major y only
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.5),
    # Axes: thin lines + subtle ticks
    axis.line          = element_line(color = "#9aa0a6", linewidth = 0.5),
    axis.ticks         = element_line(color = "#b9bfc6", linewidth = 0.5),
    axis.ticks.length  = unit(3, "pt"),
    # Axis text styling
    axis.text.x        = element_text(size = 14, color = "#444444", margin = margin(t = 4)),
    axis.text.y        = element_text(size = 12, color = "#444444", margin = margin(r = 6)),
    # Room on both sides for direct labels
    plot.margin        = margin(5.5, 40, 5.5, 40)
  )

ggsave(filename = str_glue("{root_dir}/output/figures/top150_magnet_share_ts.pdf"), width = 14, height = 8)


# WHAT SHARE OF STUDENTS ARE IN OUR SAMPLE? -------------------------------

all_student_enrollment <- df_ccd %>%
  filter(year == 2019) %>%
  pull(enrollment) %>%
  sum(na.rm = T)

survey_sample_enrollment <- df %>%
  filter(year == 2019) %>%
  pull(total_enrollment) %>%
  sum(na.rm = T)

share_of_students_in_survey_sample <- survey_sample_enrollment/all_student_enrollment

# MANDATORY VERSUS OPT-IN -------------------------------------------------

df_cross_section <- df %>%
  filter(year == 2019) %>%
  mutate(
         choice_type = case_when(
           `centralized/algorithm` == 1 & mandatory == 0 ~ "Opt-in algorithm", 
           `centralized/algorithm` == 1 & mandatory == 1 ~ "Mandatory algorithm", 
           `centralized/algorithm` == 0 & school_choice == 1 ~ "Decentralized", 
           `centralized/algorithm` == 0 & school_choice == 0 ~ "No district choice", 
         ), 
         
         choice_type_fctr = factor(choice_type, levels = c("No district choice", "Decentralized", "Opt-in algorithm", "Mandatory algorithm"), labels = c("No district choice", "Decentralized", "Opt-in algorithm", "Mandatory algorithm"))
         )


#types of choice offerings
all_enrollment <- sum(df_cross_section$enrollment)
offerings <- df_cross_section %>%
  filter(!is.na(choice_type_fctr)) %>%
  group_by(
    choice_type_fctr
  ) %>%
  mutate(
    enrollment_share = total_enrollment/sum(total_enrollment)) %>%
  summarise(
    N = n(),
    share = sum(total_enrollment)/all_enrollment,
    weight_check = sum(enrollment_share),
    zoned_school = sum(enrollment_share*zoned_school, na.rm = T),
    across(easy_to_learn:virtual_schools, ~sum(enrollment_share*., na.rm = T)) 
    
  )

all_districts <- sum(!is.na(df_cross_section$choice_type))
offerings_v2 <- df_cross_section %>%
  filter(!is.na(choice_type_fctr)) %>%
  group_by(
    choice_type_fctr
  ) %>%
  summarise(
    N = n(),
    share = N/all_districts,
    zoned_school = mean(zoned_school, na.rm = T),
    across(easy_to_learn:virtual_schools, ~mean(., na.rm = T)) 
    
  )

#demographics
enrollment <- df_cross_section %>%
  filter(!(leaid %in% no_frpl_districts), 
         !is.na(choice_type_fctr)) %>% 
  group_by(choice_type_fctr) %>%
  summarise(
    N = n(), 
    dist_charter_frpl_share = sum(dist_charter_frpl)/sum(total_frpl), 
    non_dist_charter_frpl_share = sum(non_dist_charter_frpl)/sum(total_frpl), 
    magnet_frpl_share = sum(magnet_frpl)/sum(total_frpl), 
    all_choice_share = sum(dist_charter_frpl + non_dist_charter_frpl + magnet_frpl)/sum(total_frpl),
    tps_frpl_share = sum(tps_frpl)/sum(total_frpl)
  )

enrollment_v2 <- df_cross_section %>%
  filter(!(leaid %in% no_frpl_districts), 
         !is.na(choice_type_fctr)) %>% #drop columbus, boston, and DC
  group_by(choice_type_fctr) %>%
  mutate(
    dist_charter_frpl_share = dist_charter_frpl/total_frpl, 
    non_dist_charter_frpl_share = non_dist_charter_frpl/total_frpl, 
    magnet_frpl_share = magnet_frpl/total_frpl,
    all_choice_share = (dist_charter_frpl + non_dist_charter_frpl + magnet_frpl)/total_frpl,
    tps_frpl_share = tps_frpl/total_frpl
  ) %>%
  summarise(
    N = n(), 
    dist_charter_frpl_share = mean(dist_charter_frpl_share),
    non_dist_charter_frpl_share = mean(non_dist_charter_frpl_share), 
    magnet_frpl_share = mean(magnet_frpl_share), 
    all_choice_share = mean(all_choice_share),
    tps_frpl_share = mean(tps_frpl_share)
  )
  

####
#`district weighted` figures
####

#distribution of choice types
y_max <- max(offerings_v2$share, na.rm = TRUE)
y_pad <- y_max * 0.08

choice_share_plot_v2 <-
  ggplot(offerings_v2, aes(x = choice_type_fctr, y = share)) +
  geom_col(width = 0.65, fill = "#000000", alpha = 0.8) +
  geom_text(
    aes(y = share + y_pad * 0.25, label = scales::percent(share, accuracy = 0.1)),
    vjust = 0, size = 4.6
  ) +
  scale_y_continuous(
    limits = c(0, y_max + y_pad * 1.6),  # headroom for labels
    labels = NULL, breaks = NULL, expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  labs(title = "Choice system type (% of districts)") +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 12, 5.5),  # a touch more bottom space
    axis.text.x = element_text(size = 14), 
    plot.title.position = "plot",      # use the plot edge for alignment
    plot.title = element_text(size = 16, hjust = .045) 
  )

ggsave(choice_share_plot_v2, filename = str_glue("{root_dir}/output/figures/choice_type_shares_district_weighted.pdf"), width = 14, height = 8)


#decentralized systems tend to be harder to navigate
y_max <- max(1-offerings_v2$easy_to_learn, na.rm = TRUE)
y_pad <- y_max * 0.08

easy_to_learn_plot_v2 <-
  ggplot(offerings_v2, aes(x = choice_type_fctr, y = 1 - easy_to_learn)) +
  geom_col(width = 0.65, fill = "#000000", alpha = 0.8) +
  geom_text(
    aes(y = 1-easy_to_learn + y_pad * 0.25, label = scales::percent(1 - easy_to_learn, accuracy = 0.1)),
    vjust = 0, size = 4.6
  ) +
  scale_y_continuous(
    limits = c(0, y_max + y_pad * 1.6),  # headroom for labels
    labels = NULL, breaks = NULL, expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  labs(title = "Difficult to navigate system (% of districts)") +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 12, 5.5),  # a touch more bottom space
    axis.text.x = element_text(size = 14), 
    plot.title.position = "plot",      # use the plot edge for alignment
    plot.title = element_text(size = 16, hjust = .045) 
  )

ggsave(easy_to_learn_plot_v2, filename = str_glue("{root_dir}/output/figures/difficult_to_navigate_shares_district_weighted.pdf"), width = 14, height = 8)


#decentralized systems tend to enroll higher SES students in choice programs
y_max <- max(enrollment_v2$all_choice_share, na.rm = TRUE)
y_pad <- y_max * 0.08

all_choice_share_plot_v2 <-
  ggplot(enrollment_v2, aes(x = choice_type_fctr, y = all_choice_share)) +
  geom_col(width = 0.65, fill = "#000000", alpha = 0.8) +
  geom_text(
    aes(y = all_choice_share + y_pad * 0.25, label = scales::percent(all_choice_share, accuracy = 0.1)),
    vjust = 0, size = 4.6
  ) +
  scale_y_continuous(
    limits = c(0, y_max + y_pad * 1.6),  # headroom for labels
    labels = NULL, breaks = NULL, expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  labs(title = "Average share of low-SES students enrolled in a choice school") +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 12, 5.5),  # a touch more bottom space
    axis.text.x = element_text(size = 14), 
    plot.title.position = "plot",      # use the plot edge for alignment
    plot.title = element_text(size = 16, hjust = .045) 
  )

ggsave(all_choice_share_plot_v2, filename = str_glue("{root_dir}/output/figures/frpl_shares_district_weighted_v4.pdf"), width = 14, height = 8)



