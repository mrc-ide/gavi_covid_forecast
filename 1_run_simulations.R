library(nimue)
library(ggplot2)
library(dplyr)
library(purrr)
library(furrr)
library(tidyverse)

source("R/functions.R")

# read in the gavi forecast data
gavi_data <- read.csv("scenarios_2021-11.csv", header=TRUE) %>%
 # select(iso3c, country) %>% #scenario, 
  #filter(!iso3c %in% c("DMA", "FSM", "MHL", "TUV", "XKX", "KIR", "PRK", "SLB", "TON", "WSM", "VUT")) %>% #dma, mhl, tuv, xkx no input_params.json in global-lmic-reports
  filter(iso3c %in% c("IND", "NGA","COD","IDN","BOL","PAK")) %>% #dma, mhl, tuv, xkx no input_params.json in global-lmic-reports
  mutate(june_2022_cov = 0.8,
         dec_2022_cov_global_recov = 0.8,
         vacc_scenario = "Vaccine",
         scenario="scenario_18plus")

# add in scenario_18plus and scenario_50plus options
# gavi_data_scenario_18plus <- gavi_data %>%
#   mutate(scenario="scenario_18plus",
#          june_2022_cov = 0.7,
#          dec_2022_cov_global_recov = 0.7)
# gavi_data <- rbind(gavi_data, gavi_data_scenario_18plus)
# 
# gavi_data_scenario_50plus <- gavi_data %>%
#   filter(scenario %in% c("scenario_12plus"))%>%
#   mutate(scenario="scenario_50plus",
#          june_2022_cov = 0.85,
#          dec_2022_cov_global_recov = 0.85)
# gavi_data <- rbind(gavi_data, gavi_data_scenario_50plus)

# create 0 coverage option
gavi_data_0 <- gavi_data %>%
  mutate(june_2022_cov = 0,
         dec_2022_cov_global_recov = 0,
         vacc_scenario = "Counterfactual")

gavi_data <- rbind(gavi_data, gavi_data_0)

# include other parameters
gavi_data$future_Rt <- 2
end_date <- as.Date("2022-12-31") #set as end of June 2022-06-30 or end of Dec 2022-12-31
gavi_data$forecast <- as.integer(end_date - as.Date("2021-11-15"))

gavi_data_2 <- gavi_data %>%
  mutate(future_Rt = 3.5)


scenarios <- rbind(gavi_data, gavi_data_2)

# save the parameters
write_csv(scenarios, "scenarios_2021-11_JT.csv")

#### Run the model #############################################################
plan(multiprocess, workers = 2)
system.time({out <- future_pmap(select(scenarios, -country), create_vacc_fit, .progress = TRUE)}) #, -june_2022_cov

#### Format output #############################################################
out <- do.call(rbind,out)

### Save output ################################################################
saveRDS(out, "output/output.rds")

