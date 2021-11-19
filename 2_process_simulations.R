
library(patchwork)
library(scales)
library(tidyverse)
library(countrycode)

##############################################################################
# plotting stuff

col1 <- "#5da0b5"
col2 <- "#c59e96"
col3 <- "#747473"
col4 <- "#5c8e72"
col5 <- "#2a73bb" # Reff color

lw <- 1

##############################################################################

d <- readRDS("output/output.rds")

d_counter <- d %>%
  filter(vacc_scenario == "Counterfactual") %>%
  select(scenario, iso3c, future_Rt, total_deaths, total_infections, total_hospitalisations)
colnames(d_counter)[4:6] <- paste0(colnames(d_counter)[4:6], "_counter")

d2 <- d %>%
  left_join(d_counter, by = c("scenario", "iso3c", "future_Rt")) %>%
  mutate(total_infections_averted = total_infections_counter - total_infections,
         total_deaths_averted = total_deaths_counter - total_deaths)

d2_summary <- d2 %>%
  filter(vacc_scenario == "Vaccine") %>%
  group_by(scenario, vacc_scenario, future_Rt) %>%
  summarise(global_deaths_averted = sum(total_deaths_averted), 
            vaccines = sum(total_vaccines),
            .groups = "drop")

d2_summary

d2 <- d2 %>%
  unnest(timeseries)

###########################################################################
# Dataframe: trajectories

df1 <- d2 %>%  
  mutate(vacc_scenario = if_else(date < "2021-11-16", "Fitted", vacc_scenario)) %>% #last date with "deaths" in json
  unique() %>%
  mutate(scenario = if_else(vacc_scenario == "Fitted", "Fitted", scenario),
         scenario = if_else(vacc_scenario == "Counterfactual", "Counterfactual", scenario)) %>%
    unique()

###########################################################################
# Dataframe: vaccines delivered

df2 <- d2 %>%
  select(scenario, vacc_scenario, iso3c, population, future_Rt, date, vaccines, real) %>%
  filter(future_Rt == 2)

###########################################################################
# Dataframe: Rt

df3 <- d2 %>%
  select(scenario, vacc_scenario, iso3c, future_Rt, date, Rt) %>%
  filter(vacc_scenario == "Vaccine",
         scenario == "scenario_12plus")

###########################################################################
# Create plots

iso3c_list <- unique(select(d, c(iso3c)))
iso3c_list$country <- countrycode(iso3c_list$iso3c, "iso3c", "country.name")
plots <- list()

for(i in seq_along(iso3c_list$iso3c)){
  
  # plot trajectories
  p1 <- ggplot(filter(df1, iso3c == iso3c_list$iso3c[i])) +
    geom_point(aes(x = as.Date(date), y = real),
               col = "darkgrey") +
    geom_line(aes(x = as.Date(date), y = deaths, linetype = vacc_scenario, col = scenario), size = lw) +
    geom_vline(xintercept = as.Date("2021-11-16"), linetype = "dashed") +
    facet_wrap( ~ `future_Rt`, labeller = label_both, nrow = 2) +
    labs(x = "Time", y = "Deaths per day", col = "Allocation") +
    scale_color_manual(values = c("black", "black",col1, col2, col4)) +  #counterfactual, fitted, scenario_12plus, scenario_18plus, scenario_50plus
    scale_linetype_manual(values = c("dashed", "solid", "solid")) + #counterfactual=dashed line
    scale_y_continuous(labels = comma)+
    theme_bw() +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0,
          legend.position = "none") +
    ggtitle("Epidemic trajectories")
  
  # p1_inset <- ggplot(filter(deaths_dat, iso3c == iso3c_list$iso3c[i])) +
  #   geom_point(aes(x = as.Date(date), y = deaths), col = "darkgreen") +
  #   geom_line(data = filter(df1, iso3c == iso3c_list$iso3c[i], Scenario == "Counterfactual", `Rt Scenario` == "August 2021 lift", Allocation == "1.3 billion", date <= date_current), aes(x = date, y = value), col = "darkred", size = lw) +
  #   labs(x = "Time", y = "Deaths per day") +
  #   scale_y_continuous(labels = comma) +
  #   theme_bw() +
  #   theme(strip.background = element_rect(fill = NA),
  #         panel.border = element_blank(),
  #         axis.line = element_line(),
  #         legend.text.align = 0) +
  #   ggtitle("Fitted epidemic")
  
  # plot vaccines delivered over time
  p2 <- ggplot(filter(df2, iso3c == iso3c_list$iso3c[i], vacc_scenario == "Vaccine")) +
    geom_line(aes(x = as.Date(date), y = vaccines/population, col = scenario), size = lw) +
    geom_vline(xintercept = as.Date("2021-11-16"), linetype = "dashed") + 
    labs(x = "Time", y = "Proportion vaccinated", col = "Allocation") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_color_manual(values = c(col1, col2, col4)) +
    lims(y = c(0,1)) +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0#,
          #axis.text.x = element_text(angle = 315, vjust = 0.5, hjust = -0.1)
          )+
    ggtitle("Vaccine delivery")
  
  # plot R0
  p3 <- ggplot(filter(df3, iso3c == iso3c_list$iso3c[i])) +
    geom_line(aes(x = as.Date(date), y = Rt, col = factor(future_Rt)), size = lw) +
    geom_vline(xintercept = as.Date("2021-11-16"), linetype = "dashed") +
    labs(x = "Time", y = "Rt", col = "Rt scenario") +
    scale_color_manual(values = c("darkgrey", col5)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0#,
          #axis.text.x = element_text(angle = 315, vjust = 0.5, hjust = -0.1)
          )+
    ggtitle("Reproduction number over time")
  
  plots[[i]] <- p3 / p2 / p1 + plot_annotation(title = iso3c_list$country[i], theme = theme(plot.title = element_text(size = 18))) +  plot_layout(guides = 'collect', heights = unit(c(4, 4, 8), c('cm')))
}

pdf("GAVI_plots_covid_4.pdf", width = 11, height = 10)
for (i in seq_along(iso3c_list$iso3c)){
  print(plots[[i]])
}
dev.off()
