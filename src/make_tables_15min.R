library(lqmm)
library(gtsummary)
library(lme4)
library(gt)
library(webshot)
library(car)
library(ggtext)
library(table1)
library(performance)
library(arrow)
#install.packages("table1")
#install.packages("metafor")
library(sjstats) #use for r2 functions
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(flextable)
library(ggplot2)

setwd('/Users/david/Dropbox/PhD/GitHub/15min_city/')
df <- read.csv("./data/gdf_final_1200.csv")
df_urban <- read.csv("./data/gdf_final_1200_urban.csv")
df$MVPA
result_folder <- './results'
df_wsleep <- df[!is.na(df$sleep_time), ]

df$overall_15min_city_proximity_time_cat_4 <- factor(df$overall_15min_city_proximity_time_cat_4, levels = c('<10 min','10-15 min', '15-30 min', '>30 min'), ordered=TRUE)
df_urban$overall_15min_city_proximity_time_cat_4 <- factor(df_urban$overall_15min_city_proximity_time_cat_4, levels = c('<5 min','5-10 min','10-15 min', '15-30 min', '>30 min'), ordered=TRUE)

t1 <- df %>% 
  select('sexe', 'age') %>% 
  tbl_summary(by = sexe, missing ='ifany',
              statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75})", #median and IQR
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 1,
              label = list(sexe = 'Sex',
                           age = 'Age'
              )) %>%
  add_p() %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()


t1 <- df %>% 
  select('sexe', 'age','bmi_category','employment_coded','education_coded','smoking_status','Mobility..commute...personal.time...standardized..min.day.','Leisure.time.MVPA..standardized..min.day.','Home..sedentary..standardized..min.day.','Total.adjusted.energy.kcal.day','overall_15min_city_proximity_time') %>% 
  tbl_summary(by = sexe, missing ='ifany',
              statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75})", #median and IQR
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 1,
              label = list(sexe = 'Sex',
                           age = 'Age (years)',
                           bmi_category = 'Body Mass Index (BMI)',
                           smoking_status = 'Smoking status',
                           education_coded = 'Education level',
                           employment_coded = 'Working status',
                           Leisure.time.MVPA..standardized..min.day. = 'Leisure-time moderate and vigorous PA time (min/day)',
                           Home..sedentary..standardized..min.day. = "Home sedentary time, non-occupational (min/day)",
                           Mobility..commute...personal.time...standardized..min.day. = "Active mobility time, non-occupational (min/day)",
                           Total.adjusted.energy.kcal.day = "Total adjusted energy expenditure (kcal/day)",
                           overall_15min_city_proximity_time = "Proximity Time - All POIs [min]"
              )) %>%
  add_p() %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()
t1

t1_light <-t1 %>%
  as_gt() %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(2),  # Adjust this value to change padding for data rows
    column_labels.padding = px(1),  # Adjust padding for column labels
    table.width = pct(80)  # Make table full width
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "white",
      weight = px(0)  # Adjust this value to change apparent cell padding
    ),
    locations = cells_body()
  ) %>% cols_width(
    everything() ~ px(40)  # Adjust this value as needed
  )

t1_light
gt::gtsave(as_gt(t1), file = file.path(result_folder,'Table 1.png'), vwidth=600, vheight=800)
save_as_docx(
  "Table: Descriptive statistics" = as_gt(t1),
  path = file.path(result_folder,'Table 1.docx'))
t1
save_as_docx(
  "Table: Descriptive statistics" = as_gt(t1),
  path = file.path(result_folder,'Table 1 - par sexe.docx'))

gt::gtsave(t1_light, file = file.path(result_folder,'Table 1.docx'), vwidth=600, vheight=800)



t1 <- df %>% 
  select('sexe', 'age','bmi_category','employment_coded','education_coded','smoking_status','Mobility..commute...personal.time...standardized..min.day.','Leisure.time.MVPA..standardized..min.day.','Home..sedentary..standardized..min.day.','Total.adjusted.energy.kcal.day','overall_15min_city_proximity_time') %>% 
  tbl_summary(by = sexe, missing ='ifany',
              statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75})", #median and IQR
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 1,
              label = list(sexe = 'Sex',
                           age = 'Age (years)',
                           bmi_category = 'Body Mass Index (BMI)',
                           smoking_status = 'Smoking status',
                           education_coded = 'Education level',
                           employment_coded = 'Working status',
                           Leisure.time.MVPA..standardized..min.day. = 'Leisure-time moderate and vigorous PA time (min/day)',
                           Home..sedentary..standardized..min.day. = "Home sedentary time, non-occupational (min/day)",
                           Mobility..commute...personal.time...standardized..min.day. = "Active mobility time, non-occupational (min/day)",
                           Total.adjusted.energy.kcal.day = "Total adjusted energy expenditure (kcal/day)",
                           overall_15min_city_proximity_time = "Proximity Time - All POIs [min]"
              )) %>%
  add_p() %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()
t1

t1_light <-t1 %>%
  as_gt() %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(2),  # Adjust this value to change padding for data rows
    column_labels.padding = px(1),  # Adjust padding for column labels
    table.width = pct(80)  # Make table full width
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "white",
      weight = px(0)  # Adjust this value to change apparent cell padding
    ),
    locations = cells_body()
  ) %>% cols_width(
    everything() ~ px(40)  # Adjust this value as needed
  )

t1_light
gt::gtsave(as_gt(t1), file = file.path(result_folder,'Table 1.png'), vwidth=600, vheight=800)
save_as_docx(
  "Table: Descriptive statistics" = as_gt(t1),
  path = file.path(result_folder,'Table 1.docx'))

gt::gtsave(t1_light, file = file.path(result_folder,'Table 1.docx'), vwidth=600, vheight=800)



#################

t1 <- df %>% 
  select('sexe','overall_15min_city_proximity_time_cat_4', 'age','bmi_category','employment_coded','education_coded','smoking_status','Mobility..commute...personal.time...standardized..min.day.','Leisure.time.MVPA..standardized..min.day.') %>% 
  tbl_summary(by = overall_15min_city_proximity_time_cat_4, missing ='ifany',
              statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75})", #median and IQR
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 1,
              label = list(sexe = 'Sex',
                           age = 'Age (years)',
                           bmi_category = 'Body Mass Index (BMI)',
                           smoking_status = 'Smoking status',
                           education_coded = 'Education level',
                           employment_coded = 'Working status',
                           Leisure.time.MVPA..standardized..min.day. = 'Leisure-time moderate and vigorous PA time (min/day)',
                           Mobility..commute...personal.time...standardized..min.day. = "Active mobility time, non-occupational (min/day)"
              )) %>%
  add_p() %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()
t1

t1_light <-t1 %>%
  as_gt() %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(2),  # Adjust this value to change padding for data rows
    column_labels.padding = px(1),  # Adjust padding for column labels
    table.width = pct(80)  # Make table full width
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "white",
      weight = px(0)  # Adjust this value to change apparent cell padding
    ),
    locations = cells_body()
  ) %>% cols_width(
    everything() ~ px(40)  # Adjust this value as needed
  )

t1_light
gt::gtsave(as_gt(t1), file = file.path(result_folder,'Table 1 - Proximity.png'), vwidth=1000, vheight=1000)
 save_as_docx(
  "Table: Descriptive statistics" = as_gt(t1),
  path = file.path(result_folder,'Table 1.docx'))

gt::gtsave(t1_light, file = file.path(result_folder,'Table 1 - Proximity.docx'), vwidth=600, vheight=800)



t1 <- df_urban %>% 
  select('sexe','overall_15min_city_proximity_time_cat_4', 'age','bmi_category','employment_coded','education_coded','smoking_status','Mobility..commute...personal.time...standardized..min.day.','Leisure.time.MVPA..standardized..min.day.') %>% 
  tbl_summary(by = overall_15min_city_proximity_time_cat_4, missing ='ifany',
              statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75})", #median and IQR
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 1,
              label = list(sexe = 'Sex',
                           age = 'Age (years)',
                           bmi_category = 'Body Mass Index (BMI)',
                           smoking_status = 'Smoking status',
                           education_coded = 'Education level',
                           employment_coded = 'Working status',
                           Leisure.time.MVPA..standardized..min.day. = 'Leisure-time moderate and vigorous PA time (min/day)',
                           Mobility..commute...personal.time...standardized..min.day. = "Active mobility time, non-occupational (min/day)"
              )) %>%
  add_p() %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()
t1

t1_light <-t1 %>%
  as_gt() %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(2),  # Adjust this value to change padding for data rows
    column_labels.padding = px(1),  # Adjust padding for column labels
    table.width = pct(80)  # Make table full width
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "white",
      weight = px(0)  # Adjust this value to change apparent cell padding
    ),
    locations = cells_body()
  ) %>% cols_width(
    everything() ~ px(40)  # Adjust this value as needed
  )

t1_light
gt::gtsave(as_gt(t1), file = file.path(result_folder,'Table 1 - Proximity.png'), vwidth=1000, vheight=1000)
save_as_docx(
  "Table: Descriptive statistics" = as_gt(t1),
  path = file.path(result_folder,'Table 1.docx'))

gt::gtsave(t1_light, file = file.path(result_folder,'Table 1 - Proximity - Urban.docx'), vwidth=600, vheight=800)




