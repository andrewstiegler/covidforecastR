# Master script

library(EpiNow2) # for estimating Rt and forecasting infections
library(covid19us) # for COVID19 Tracking Project API
library(tidyverse) # dplyr data manipulation, ggplot
library(data.table) # for rbindlist
library(tsibble) 
library(ggsci) # NEJM color palett
library(leaflet)
library(RColorBrewer)
library(viridis)
library(albersusa)
library(zoo)
library(leaflet.extras)

# Helper functions------

daily_state_update <- function(state){
    daily_tibble <- get_states_daily(state) %>% filter(date > "2020-03-18")
    daily_tsibble <- as_tsibble(daily_tibble, index = "date")
    daily_tsibble$pct_pos <- NA
    daily_tsibble$pct_pos <- (daily_tsibble$positive_increase / 
                                  daily_tsibble$total_test_results_increase) * 100
    return(daily_tsibble)
}

select_pos_only <- function(data){
    data <- data %>% dplyr::select(date, state, positive_increase)
    return(data)
}

state_prep_for_epinow <- function(data){
    data <- data %>% dplyr::select(date, positive_increase)
    colnames(data)[2] <- "confirm"
    return(data)
}

state_forecast_extract <- function(data){
    data_extract <- data$estimates$summarised %>% tibble() %>%
        filter(date > Sys.Date() & variable == "reported_cases")
    return(data_extract)
}

state_rt_forecast_extract <- function(data){
    data_extract <- data$estimates$summarised %>% tibble() %>%
        filter(date > Sys.Date() & variable == "R")
    return(data_extract)
}

state_rt_extract <- function(data){
    data_extract <- data$estimates$summarised %>% tibble() %>%
        filter(date <= Sys.Date() & variable == "R")
    return(data_extract)
}

state_plots <- function(state_infections, state_forecast, state_rt, state_rt_forecast){
    
    state_infection_plot <- ggplot()+
        geom_col(data = state_infections, aes(x = date, y = positive_increase))+
        geom_line(data = state_infections, aes(x = date, y = avg_7_day),
                  color = plot_colors[2], size = 2)+
        ylab("New positive cases")+
        geom_ribbon(data = state_forecast, aes(x = date, ymin = lower_90, ymax = upper_90),
                    fill = plot_colors[1], alpha = 0.1)+
        geom_ribbon(data = state_forecast, aes(x = date, ymin = lower_50, ymax = upper_50),
                    fill = plot_colors[1], alpha = 0.25)+
        geom_ribbon(data = state_forecast, aes(x = date, ymin = lower_20, ymax = upper_20),
                    fill = plot_colors[1], alpha = 0.45)+
        geom_line(data = state_forecast, aes(x = date, y = median),
                  color = plot_colors[1], size = 1.1)+
        scale_x_date("", limits = c(as.Date("2020-03-25"), Sys.Date() + 14))+
        coord_cartesian(ylim = c(0, 1.05*max(state_forecast$upper_50, state_infections$positive_increase)))+
        annotate("segment", x = Sys.Date(), xend = Sys.Date(), y = 0, yend = 1000000,
                 linetype = "dashed")+
        labs(caption = "source: The COVID Tracking Project (https://covidtracking.com) / forecast with R/EpiNow2")+
        theme_bw()
    
    state_rt_plot <- ggplot()+
        geom_ribbon(data = state_rt, aes(x = date, ymin = lower_90, ymax = upper_90),
                    fill = plot_colors[2], alpha = 0.1)+
        geom_ribbon(data = state_rt, aes(x = date, ymin = lower_50, ymax = upper_50),
                    fill = plot_colors[2], alpha = 0.25)+
        geom_ribbon(data = state_rt, aes(x = date, ymin = lower_20, ymax = upper_20),
                    fill = plot_colors[2], alpha = 0.45)+
        geom_ribbon(data = state_rt_forecast, aes(x = date, ymin = lower_90, ymax = upper_90),
                    fill = plot_colors[1], alpha = 0.1)+
        geom_ribbon(data = state_rt_forecast, aes(x = date, ymin = lower_50, ymax = upper_50),
                    fill = plot_colors[1], alpha = 0.25)+
        geom_ribbon(data = state_rt_forecast, aes(x = date, ymin = lower_20, ymax = upper_20),
                    fill = plot_colors[1], alpha = 0.45)+
        geom_line(data = state_rt, aes(x = date, y = median), color = plot_colors[2], size = 1.1)+
        geom_line(data = state_rt_forecast, aes(x = date, y = median), color = plot_colors[1], size = 1.1)+
        geom_hline(yintercept = 1, linetype = "dashed")+
        ylab("Effective reproduction (R)")+
        xlab("")+
        coord_cartesian(ylim = c(.95*min(state_rt$lower_50, state_rt_forecast$lower_50), 
                                 1.05*max(state_rt$upper_50, state_rt_forecast$upper_50)))+
        annotate("segment", x = Sys.Date(), xend = Sys.Date(), y = 0, yend = 1000000,
                 linetype = "dashed")+
        # annotate("text", x = Sys.Date() + 0.5, y = 1.05*max(state_rt$upper_50, state_rt_forecast$upper_50),
        #          label = "Forecast", hjust = 0)+
        labs(caption = "source: The COVID Tracking Project (https://covidtracking.com) / Rt estimate and forecast with R/EpiNow2")+
        theme_bw()
    
    state_plots <- list(infection_plot = state_infection_plot, rt_plot = state_rt_plot)
    return(state_plots)
}

state_delta <- function(state){
    data <- state_estimates_data[[state]]
    
    data_initial <- filter(data[[1]], date >= (Sys.Date() - 7)) %>% tibble() %>%
        select(positive_increase) %>% colMeans(na.rm = TRUE)
    data_forecast <- filter(data[[2]], date > (Sys.Date() + 7)) %>% tibble() %>%
        select(median) %>% colMeans(na.rm = TRUE)
    inf_lower_50 <- filter(data[[2]], date >= (Sys.Date() - 7)) %>% tibble() %>%
        select(lower_50) %>% colMeans(na.rm = TRUE)
    inf_upper_50 <- filter(data[[2]], date >= (Sys.Date() - 7)) %>% tibble() %>%
        select(upper_50) %>% colMeans(na.rm = TRUE)
    
    rt_initial <- filter(data[[3]], date >= (Sys.Date() - 7)) %>% tibble() %>%
        select(median) %>% colMeans(na.rm = TRUE)
    rt_lower_50 <- filter(data[[3]], date >= (Sys.Date() - 7)) %>% tibble() %>%
        select(lower_50) %>% colMeans(na.rm = TRUE)
    rt_upper_50 <- filter(data[[3]], date >= (Sys.Date() - 7)) %>% tibble() %>%
        select(upper_50) %>% colMeans(na.rm = TRUE)
    
    output <- tibble(abb = state, 
                     forecast_increase = (data_forecast / data_initial - 1) * 100,
                     current_inf = data_initial,
                     inf_lower_50 = inf_lower_50,
                     inf_upper_50 = inf_upper_50,
                     current_rt = rt_initial,
                     rt_lower_50 = rt_lower_50,
                     rt_upper_50 = rt_upper_50 )
    return(output)
}
# Misc stuff -----
setDTthreads(0)

state_list <- list(
    "AL",
    "AK",
    "AZ",
    "AR",
    "CA",
    "CO",
    "CT",
    "DE",
    "DC",
    "FL",
    "GA",
    "HI",
    "ID",
    "IL",
    "IN",
    "IA",
    "KS",
    "KY",
    "LA",
    "ME",
    "MD",
    "MA",
    "MI",
    "MN",
    "MS",
    "MO",
    "MT",
    "NE",
    "NV",
    "NH",
    "NJ",
    "NM",
    "NY",
    "NC",
    "ND",
    "OH",
    "OK",
    "OR",
    "PA",
    "RI",
    "SC",
    "SD",
    "TN",
    "TX",
    "UT",
    "VT",
    "VA",
    "WA",
    "WV",
    "WI",
    "WY"
)

plot_colors <- pal_nejm("default")(8)

# Main update ------
# get daily data from COVID19 Tracking Project
state_data <- map(state_list, daily_state_update)
us_data <- get_us_daily()
us_data$avg_7_pos_inc <- floor(rollmean(us_data$positive_increase, k = 7,
                                  align = "left", fill = 0))

# transform data for epinow prediction
state_daily_infected <- lapply(state_data, select_pos_only)
state_daily_infected_tibble <- rbindlist(state_daily_infected) %>% tibble() %>%
    filter(date > (Sys.Date() - 60))
state_daily_7day <- map(state_daily_infected, 
                        ~rollmean(.x$positive_increase, k = 7, align = "left"))

for (i in 1:length(state_list)) {
    state_daily_infected[[i]] <- filter(state_daily_infected[[i]], 
                                        date > as.Date("2020-03-24"))
    state_daily_infected[[i]]$avg_7_day <- state_daily_7day[[i]]
    i <- i + 1
}

colnames(state_daily_infected_tibble)[3] <- "confirm"
us_daily_infected <- us_data %>% select(date, positive_increase) %>% 
    filter(date > (Sys.Date() - 60))
colnames(us_daily_infected)[2] <- "confirm"

# set epinow parameters
reporting_delay <- bootstrapped_dist_fit(rlnorm(1000, log(3), 1), max_value = 15,
                                         bootstraps = 1)

generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")

# epinow predict entire us
setup_logging("WARN")
us_estimates <- epinow(reported_cases = us_daily_infected, 
                       generation_time = generation_time,
                       delays = list(incubation_period, reporting_delay),
                       horizon = 14,
                       samples = 750,
                       max_execution_time = 1800)

# epinow predict states
state_epinow <- lapply(state_data, state_prep_for_epinow)
state_epinow <- map(state_epinow, ~filter(.x, date > (Sys.Date() - 30)))
names(state_epinow) <- state_list

state_estimates_epinow <- list()
state_estimates_data <- list()
state_estimates_plot <- list()

# repeat this because it crashes if run too many times
for(i in 1:25) {
    state_estimate <- try(epinow(reported_cases = state_epinow[[i]], 
                                 generation_time = generation_time,
                                 delays = list(incubation_period, reporting_delay),
                                 horizon = 14,
                                 samples = 700,
                                 max_execution_time = 1800))
    
    state_estimates_epinow[[i]] <- state_estimate
    
    i <- i + 1
}
# same repeat
for(i in 26:length(state_list)) {
    state_estimate <- try(epinow(reported_cases = state_epinow[[i]], 
                                 generation_time = generation_time,
                                 delays = list(incubation_period, reporting_delay),
                                 horizon = 14,
                                 samples = 700,
                                 max_execution_time = 1800))
    
    state_estimates_epinow[[i]] <- state_estimate
    
    i <- i + 1
}

for(i in 1:length(state_list)){
    
    state_infections <- state_daily_infected[[i]]
    
    if(is.list(state_estimates_epinow[[i]])) {
        state_forecast <- state_estimates_epinow[[i]] %>% state_forecast_extract()
        state_rt <- state_estimates_epinow[[i]] %>% state_rt_extract()
        state_rt_forecast <- state_estimates_epinow[[i]] %>% state_rt_forecast_extract()    
        
        state_plot_data <- list(state_infections, state_forecast, 
                                state_rt, state_rt_forecast)
        plots <- state_plots(state_infections, state_forecast, 
                             state_rt, state_rt_forecast)
        
        plots[[1]] <- plots[[1]] + ggtitle(paste(bquote(.(Sys.Date())),
                                                 bquote(.(state_list[[i]])), 
                                                 "new COVID19 cases + forecast",
                                                 sep = " - "))
        plots[[2]] <- plots[[2]] + ggtitle(paste(bquote(.(Sys.Date())), 
                                                 bquote(.(state_list[[i]])), 
                                                 "reproduction number + forecast",
                                                 sep = " - "))
    } else state_plot_data <- state_estimates_plot <- "Modeling Failed"
    
    
    state_estimates_data[[i]] <- state_plot_data
    state_estimates_plot[[i]] <- plots
    
    i <- i + 1
}

names(state_estimates_epinow) <- names(state_estimates_data) <- names(state_estimates_plot) <- state_list

# separate forecats from epinow output
us_forecast <- us_estimates$estimates$summarised %>% tibble() %>% 
    filter(date > Sys.Date() & variable == "reported_cases") 

# Label functions for plotly
pos_text <- function(data){
    text <- paste("Novel positive cases", data, sep = ": ")
    return(text)
}

avg_7_text <- function(data){
    text <- paste("7-day average", data, sep = ": ")
    return(text)
}

forecast_text <- function(data){
    text <- paste("Forecast", floor(data), sep = ": ")
    return(text)
}

# US cases plot
us_infection_plot <- ggplot()+
    geom_col(data = us_data, aes(x = date, y = positive_increase, text = pos_text(positive_increase)))+
    geom_line(data = us_data, aes(x = date, y = avg_7_pos_inc, group = 1, text = avg_7_text(avg_7_pos_inc)),
              color = plot_colors[2], size = 2)+
    ylab("New positive cases")+
    geom_ribbon(data = us_forecast, aes(x = date, ymin = lower_90, ymax = upper_90),
                fill = plot_colors[1], alpha = 0.1)+
    geom_ribbon(data = us_forecast, aes(x = date, ymin = lower_50, ymax = upper_50),
                fill = plot_colors[1], alpha = 0.25)+
    geom_ribbon(data = us_forecast, aes(x = date, ymin = lower_20, ymax = upper_20),
                fill = plot_colors[1], alpha = 0.45)+
    geom_line(data = us_forecast, aes(x = date, y = median, group = 1, text = forecast_text(median)),
              color = plot_colors[1], size = 1.1)+
    #ggtitle(paste(bquote(.(Sys.Date())), "new US COVID19 cases + forecast", sep = " - "))+
    scale_x_date("", limits = c(as.Date("2020-03-20"), Sys.Date() + 14))+
    annotate("segment", x = Sys.Date(), xend = Sys.Date(), y = 0, yend = 1000000,
             linetype = "dashed")+
    # annotate("text", x = Sys.Date() + 0.5, y = 1.05*max(us_forecast$upper_50),
    #          label = "Forecast",
    #          hjust = 0)+
    coord_cartesian(ylim = c(0, 1.05*max(us_forecast$upper_50)))+
    labs(caption = "source: The COVID Tracking Project (https://covidtracking.com) / forecast with R/EpiNow2")+
    theme_bw()

# US Rt plot
us_rt_forecast <- us_estimates$estimates$summarised %>% tibble() %>% 
    filter(date > Sys.Date() & variable == "R") 

us_rt <- us_estimates$estimates$summarised %>% tibble() %>% 
    filter(date <= Sys.Date() & variable == "R") 

est_r_label <- function(data){
    data <- round(data, 2)
    label <- paste("Estimated R: ", data, sep = "")
    return(label)
}

forecast_r_label <- function(data){
    data <- round(data, 2)
    label <- paste("Forecast R: ", data, sep = "")
    return(label)
}

us_rt_plot <- ggplot()+
    geom_ribbon(data = us_rt, aes(x = date, ymin = lower_90, ymax = upper_90),
                fill = plot_colors[2], alpha = 0.1)+
    geom_ribbon(data = us_rt, aes(x = date, ymin = lower_50, ymax = upper_50),
                fill = plot_colors[2], alpha = 0.25)+
    geom_ribbon(data = us_rt, aes(x = date, ymin = lower_20, ymax = upper_20),
                fill = plot_colors[2], alpha = 0.45)+
    geom_ribbon(data = us_rt_forecast, aes(x = date, ymin = lower_90, ymax = upper_90),
                fill = plot_colors[1], alpha = 0.1)+
    geom_ribbon(data = us_rt_forecast, aes(x = date, ymin = lower_50, ymax = upper_50),
                fill = plot_colors[1], alpha = 0.25)+
    geom_ribbon(data = us_rt_forecast, aes(x = date, ymin = lower_20, ymax = upper_20),
                fill = plot_colors[1], alpha = 0.45)+
    geom_line(data = us_rt, aes(x = date, y = median, group = 1, text = est_r_label(median)), color = plot_colors[2], size = 1.1)+
    geom_line(data = us_rt_forecast, aes(x = date, y = median, group = 1, text = forecast_r_label(median)), color = plot_colors[1], size = 1.1)+
    geom_hline(yintercept = 1, linetype = "dashed")+
    ylab("Effective reproduction no. (R)")+
    #ggtitle(paste(bquote(.(Sys.Date())), "US reproduction number + forecast", sep = " - "))+
    xlab("")+
    annotate("segment", x = Sys.Date(), xend = Sys.Date(), y = -10, yend = 10,
             linetype = "dashed")+
    # annotate("text", x = Sys.Date() + 0.5, y = 1.05*max(us_rt$upper_50, us_rt_forecast$upper_50),
    #          label = "Forecast", hjust = 0)+
    labs(caption = "source: The COVID Tracking Project (https://covidtracking.com) / Rt estimate and forecast with R/EpiNow2")+
    coord_cartesian(ylim = c(.95*min(us_rt$lower_50, us_rt_forecast$lower_50), 
                             1.05*max(us_rt$upper_50, us_rt_forecast$upper_50)))+
    theme_bw()

# Return data tables for future analysis
return_list <- list("US estimates" = us_estimates, "US infection forecast" = us_forecast,
                    "State estimate data" = state_estimates_data, 
                    "US daily Rt" = us_rt, "US Rt forecast" = us_rt_forecast,
                    "US infection plot" = us_infection_plot, "US Rt plot" = us_rt_plot,
                    "State plots" = state_estimates_plot)    

state_estimates_data <- return_list[[3]]

all_states_delta <- purrr::map(state_list, state_delta) %>% rbindlist()

all_states_delta$name <- c((state.name %>% tolower())[1:8], 
                           "district of columbia", 
                           (state.name %>% tolower())[9:50])


states_js <- usa_sf()

state_delta_tibble <- left_join(tibble(name = states_js$name %>% tolower()), 
                                all_states_delta, by = "name")

states_js$current_inf <- state_delta_tibble$current_inf
states_js$current_rt <- state_delta_tibble$current_rt
states_js$forecast_inf <- state_delta_tibble$forecast_increase
states_js$inf_100k <- state_delta_tibble$current_inf / (states_js$pop_2014 / 100000)
states_js$forecast_100k <- states_js$inf_100k * (1 + states_js$forecast_inf / 100)

basemap <- leaflet(states_js, options = leafletOptions(minZoom = 4, maxZoom = 7)) %>%
    setView(-96, 37.8, 4) %>%
    addProviderTiles("MapBox", options = providerTileOptions(
        id = "mapbox.light",
        accessToken = Sys.getenv('MAPBOX_ACCESS_TOKEN')))

pos_bins <- c(0, 5, 10, 20, 50, 100, 150, Inf)
pos_pal <- colorBin("inferno", domain = states_js$inf_100k, bins = pos_bins)

pos_labels <- sprintf(
    "<strong>%s</strong><br/>%g cases per 100k population",
    states_js$name, round(states_js$inf_100k,1)
) %>% lapply(htmltools::HTML)

cases_per_100k_map <- basemap %>% addPolygons(
    fillColor = ~pos_pal(inf_100k),
    weight = 2,
    opacity = 1,
    color = "white",
    dashArray = "3",
    fillOpacity = 0.7,
    highlight = highlightOptions(
        weight = 5,
        color = "#666",
        dashArray = "",
        fillOpacity = 0.7,
        bringToFront = TRUE),
    label = pos_labels,
    labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "15px",
        direction = "auto")) %>% 
    addLegend(pal = pos_pal, values = ~inf_100k, opacity = 0.7, title = "Cases per 100k",
              position = "bottomright", na.label = element_blank()) %>% suspendScroll(sleepOpacity = 0.9)

r_bins <- c(Inf, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7,-Inf)
r_pal <- colorBin(c("#e41a1c", "#377eb8", "#4daf4a"), domain = states_js$current_rt, bins = r_bins, reverse = TRUE)

r_labels <- sprintf(
    "<strong>%s</strong><br/>%g current R",
    states_js$name, round(states_js$current_rt,2)
) %>% lapply(htmltools::HTML)

current_r_map <- basemap %>% addPolygons(
    fillColor = ~r_pal(current_rt),
    weight = 2,
    opacity = 1,
    color = "white",
    dashArray = "3",
    fillOpacity = 0.7,
    highlight = highlightOptions(
        weight = 5,
        color = "#666",
        dashArray = "",
        fillOpacity = 0.7,
        bringToFront = TRUE),
    label = r_labels,
    labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "15px",
        direction = "auto")) %>% 
    addLegend(colors = colorBin(c("#e41a1c", "#377eb8", "#4daf4a"), 
                                domain = states_js$current_rt, bins = r_bins, reverse = TRUE)(r_bins)[-1], values = ~current_rt, opacity = 0.7, title = "Estimated R",
              position = "bottomright",
              labels = c("> 1.3", "1.2 - 1.3", "1.1 - 1.2", "1.0 - 1.1",
              "0.9 - 1.0", "0.8 - 0.9", "0.7 - 0.8", "< 0.7" )) %>% suspendScroll(sleepOpacity = 0.9)

forecast_labels <- sprintf(
    "<strong>%s</strong><br/>%g Forecast cases per 100k",
    states_js$name, round(states_js$forecast_100k,1)
) %>% lapply(htmltools::HTML)


forecast_map <- basemap %>% addPolygons(
    fillColor = ~pos_pal(forecast_100k),
    weight = 2,
    opacity = 1,
    color = "white",
    dashArray = "3",
    fillOpacity = 0.7,
    highlight = highlightOptions(
        weight = 5,
        color = "#666",
        dashArray = "",
        fillOpacity = 0.7,
        bringToFront = TRUE),
    label = forecast_labels,
    labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "15px",
        direction = "auto")) %>% 
    addLegend(pal = pos_pal, values = ~forecast_100k, opacity = 0.7, title = "Cases per 100k",
              position = "bottomright") %>% suspendScroll(sleepOpacity = 0.9)


good_map_list <- list(cases_per_100k_map, current_r_map, forecast_map)

us_estimates <- return_list[[1]]
state_estimates_plot <- return_list[[8]]
us_infection_pot <- return_list[[6]]
us_rt_plot <- return_list[[7]]

state_pop <- read.csv("https://raw.githubusercontent.com/ds4stats/r-tutorials/master/intro-maps/data/StatePopulation.csv", as.is = TRUE)

all_states_for_bar <- all_states_delta
all_states_for_bar$region <- all_states_for_bar$name
all_states_for_bar <- inner_join(all_states_for_bar, state_pop, by = "region")
all_states_for_bar$region <- all_states_for_bar$region %>% toupper()

all_states_for_bar$region <- factor(all_states_for_bar$region) %>% fct_reorder(all_states_for_bar$current_rt)

r_bar_label <- function(data){
    data <- round(data, 2)
    label <- paste("Current R: ", data, sep = "")
    return(label)
}

current_inf_label <- function(data){
    data <- floor(data)
    label <- paste("Current infections per 100k: ", data, sep = "")
    return(label)
}

forecast_inf_label <- function(data){
    data <- floor(data)
    label <- paste("Forecast infections per 100k: ", data, sep = "")
    return(label)
}

current_r_bar <- ggplot(all_states_for_bar)+
    geom_crossbar(aes(x = current_rt, xmin = rt_lower_50, xmax = rt_upper_50, 
                        y = region, color = current_rt, text = r_bar_label(current_rt)), 
                    size = 0.1, show.legend = FALSE, fatten = 10)+
#    geom_point(aes(x = current_rt, 
 #                  y = region, color = current_rt, text = r_bar_label(current_rt)), 
             #  size = 2.5, show.legend = FALSE)+
    annotate("segment", x = 1, xend = 1, y = 0, yend = 50, linetype = "dashed")+
    scale_color_gradient2("Rt", low = "#4daf4a", high = "#e41a1c", mid = "#377eb8",
                         midpoint = 1)+
    theme_bw(base_size = 12)+
    ylab("")+
    scale_x_continuous("", position = "top")

all_states_for_bar$region <- factor(all_states_for_bar$region) %>% fct_reorder(all_states_for_bar$current_inf / all_states_for_bar$population)

current_inf_bar <- ggplot(all_states_for_bar)+
    geom_col(aes(x = current_inf / population * 100000, y = region, 
                 fill = current_inf / population * 100000, 
                 text = current_inf_label(current_inf / population * 100000)), show.legend = FALSE)+
    # annotate("segment", x = 1, xend = 1, y = 0, yend = 60, linetype = "dashed")+
    scale_fill_viridis("Cases per 100k", "inferno")+
    scale_color_viridis("Cases per 100k", "inferno")+
    theme_bw(base_size = 12)+
    ylab("")+
    scale_x_continuous("", position = "top")


all_states_for_bar$forecast_per_100k <- (all_states_for_bar$current_inf / all_states_for_bar$population * 100000) *
    (1 + all_states_for_bar$forecast_increase / 100)
all_states_for_bar$region <- factor(all_states_for_bar$region) %>% fct_reorder(all_states_for_bar$forecast_per_100k)

forecast_inf_bar <- ggplot(all_states_for_bar)+
    geom_crossbar(aes(x = (current_inf / population) * 100000 * (1 + forecast_increase / 100),
                        xmin = inf_lower_50 / population * 100000,
                        xmax = inf_upper_50 / population * 100000,
                        y = region, text = forecast_inf_label((current_inf / population) * 100000 * (1 + forecast_increase / 100)),
                        color = current_inf / population * 100000 * forecast_increase,
    ), size = 0.1, fatten = 10, show.legend = FALSE)+
    # annotate("segment", x = 1, xend = 1, y = 0, yend = 60, linetype = "dashed")+
    scale_color_viridis("Cases per 100k", "inferno")+
    theme_bw(base_size = 12)+
    ylab("")+
    scale_x_continuous("", position = "top")


forecast_bar <- ggplot(all_states_for_bar)+
    geom_col(aes(x = forecast_increase, y = region, fill = forecast_increase, 
                 text = forecast_inf_label(forecast_increase)), show.legend = FALSE)+
    annotate("segment", x = 0, xend = 0, y = 0, yend = 50, linetype = "dashed")+
    scale_fill_gradient2("Rt", low = "#4daf4a", high = "#e41a1c", mid = "#377eb8",
                         midpoint = 0)+
    theme_bw(base_size = 12)+
    ylab("")+
    xlab("")

bar_list <- list(current_r_bar, current_inf_bar, forecast_inf_bar)


saveRDS(good_map_list, "docs/goodmaps.rds")
saveRDS(bar_list, "docs/bar_list.rds")
saveRDS(us_estimates, "docs/us_estimates.rds")
saveRDS(state_estimates_data, "docs/state_data.rds")
saveRDS(state_estimates_plot,"docs/state_plots.rds")
saveRDS(list(us_infection_plot, us_rt_plot), "docs/us_plots.rds")

