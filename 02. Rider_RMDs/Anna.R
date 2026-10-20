library(tidyverse)
library(ismev)
library(lubridate)
library(patchwork)

## Initial checks of the data
# Checking the coverage of the data
coverage_check <- river_data %>%
  arrange(Timestamp)
coverage_check

# Checking for duplicate timestamps or time gaps
data_checking <- coverage_check %>%
  summarise(
    n_observations = n(),
    first_timestamp = min(Timestamp),
    last_timestamp = max(Timestamp),
    missing_discharge = sum(is.na(Discharge)),
    duplicate_timestamps = sum(duplicated(Timestamp))
  )
data_checking

# Inspect time intervals to check they are consitent
interval_check <- river_data %>%
  mutate(
    interval_minutes = as.numeric(
      difftime(Timestamp, lag(Timestamp), units = "mins")
    )
  ) %>%
  count(interval_minutes, sort = TRUE)



## 3-5 peaks method
#Check the mean residual life plot
mrl.plot(river_data$Discharge) # Stable around 80-200 

# Checking the stability of parameter estimates
gpd.fitrange(river_data$Discharge, 80, 200)
gpd.fitrange(river_data$Discharge, 100)

# Trying 100 as a thresholdand visulasing this 
u <- 100

ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
  geom_line(colour = "steelblue", linewidth = 0.3) +
  geom_hline(yintercept = u,
             colour = "red",
             linetype = "dashed") +
  labs(
    title = "River Don discharge with POT threshold",
    x = "Date",
    y = "Discharge (m³/s)"
  ) +
  theme_minimal()

# trying multiple theresholds
ggplot(river_data, aes(Timestamp, Discharge)) +
  geom_line(colour = "steelblue", linewidth = 0.3) +
  geom_hline(
    yintercept = c(80, 100, 120, 150, 180),
    colour = c("pink", "red", "orange", "green", "purple"),
    linetype = "dashed"
  ) +
  labs(
    title = "Comparing POT thresholds",
    x = "Date",
    y = "Discharge (m³/s)"
  ) +
  theme_minimal()


# Looking at it by year
river_data_year <- river_data %>%
  mutate(year = year(Timestamp))

thresholds <- data.frame(
  threshold = c(80, 100, 120, 150, 180),
  colour = c("pink", "red", "orange", "green", "purple")
)

ggplot(river_data_year, aes(Timestamp, Discharge)) +
  geom_line(colour = "steelblue", linewidth = 0.3) +
  geom_hline(
    data = thresholds,
    aes(yintercept = threshold, colour = factor(threshold)),
    linetype = "dashed"
  ) +
  scale_colour_manual(
    name = "Threshold (m³/s)",
    values = thresholds$colour
  ) +
  facet_wrap(~ year, scales = "free_x") +
  labs(
    title = "Comparing POT thresholds by year",
    x = "Date",
    y = "Discharge (m³/s)"
  ) +
  theme_bw()


# Lowering the thresholds for more per year and examining 
thresholds <- data.frame(
  threshold = c(40, 60, 70, 80, 100),
  colour = c("blue", "green", "purple", "pink", "red")
)

ggplot(river_data_year, aes(Timestamp, Discharge)) +
  geom_line(colour = "steelblue", linewidth = 0.3) +
  geom_hline(
    data = thresholds,
    aes(yintercept = threshold, colour = factor(threshold)),
    linetype = "dashed"
  ) +
  scale_colour_manual(
    name = "Threshold (m³/s)",
    values = thresholds$colour
  ) +
  facet_wrap(~ year, scales = "free_x") +
  labs(
    title = "Comparing POT thresholds by year",
    x = "Date",
    y = "Discharge (m³/s)"
  ) +
  theme_bw()

# Calculating the number of candidate peaks at each threshold
thresholds <- c(40, 60, 70, 80, 100)

peak_counts <- data.frame(
  threshold = thresholds,
  candidate_peaks = sapply(thresholds, function(u) {
    sum(peaks$Discharge > u)
  })
)

peak_counts



# Intial data exploratory graph
background <- river_data_year %>%
  rename(background_year = year)

ggplot() +
  geom_line(
    data = background,
    aes(
      x = Timestamp,
      y = Discharge,
      group = background_year
    ),
    colour = "grey75",
    linewidth = 0.4,
    alpha = 0.4
  ) +
  geom_line(
    data = river_data_year,
    aes(
      x = Timestamp,
      y = Discharge,
      group = year
    ),
    colour = "#69b3a2",
    linewidth = 0.9
  ) +
  facet_wrap(~year, ncol = 3) +
  scale_x_continuous(
    breaks = c(1, 91, 182, 274, 365),
    labels = c("Jan", "Apr", "Jul", "Oct", "Dec")
  ) +
  labs(
    title = "Annual hydrograph comparison: Don at Parkhill",
    x = "Month",
    y = expression("Discharge (m"^3*"/s)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )



# Same graph but with the threshold and max peaks added

thresholds_final <- data.frame(
  threshold = 65)

ggplot() +
  geom_line(
    data = background,
    aes(
      x = Timestamp,
      y = Discharge,
      group = background_year
    ),
    colour = "grey75",
    linewidth = 0.4,
    alpha = 0.4
  ) +
  geom_line(
    data = river_data_year,
    aes(
      x = Timestamp,
      y = Discharge,
      group = year
    ),
    colour = "steelblue",
    linewidth = 0.9
  ) +
  geom_hline(
    data = thresholds_final,
    aes(yintercept = threshold),
    linetype = "dashed", colour = "red"
  ) + 
  facet_wrap(~year, ncol = 3) +
  scale_x_continuous(
    breaks = c(1, 91, 182, 274, 365),
    labels = c("Jan", "Apr", "Jul", "Oct", "Dec")
  ) +
  labs(
    title = "Annual hydrograph comparison: Don at Parkhill",
    x = "Month",
    y = expression("Discharge (m"^3*"/s)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )


# Threshold: the 95th percentile of flow
u95 <- quantile(river_data$Discharge, 0.95)
u95

# Peaks: the highest reading in each separate period above the threshold
peaks_95 <- river_data %>%
  mutate(above = Discharge > u95,
         start = above & !lag(above, default = FALSE),   # TRUE when the flow first goes above the threshold
         run   = cumsum(start)) %>%                       # numbers each period above the threshold
  filter(above) %>%
  group_by(run) %>%
  slice_max(Discharge, n = 1, with_ties = FALSE) %>%
  ungroup()

# Number of peaks, and average peaks per year
nrow(peaks_95)
nrow(peaks_95) / n_distinct(river_data$water_year)

# Plot
ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
  geom_line(color = "steelblue", linewidth = 0.4) +
  geom_hline(yintercept = u95, linetype = "dashed", color = "red") +
  geom_point(data = peaks_95, color = "red", size = 1) +
  labs(
    title = "Peaks above the 95th percentile threshold",
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()



# Check 95th percentile threshold
u95 <- as.numeric(
  quantile(river_data$Discharge, 0.95, na.rm = TRUE)
)

cat("95th percentile threshold:", u95, "m3/s\n")

# Check candidate thresholds against 3-5 events/year
# Note: use complete observation years where possible.

years_observed <- as.numeric(
  difftime(
    max(river_data$Timestamp),
    min(river_data$Timestamp),
    units = "days"
  )
) / 365.25

candidate_thresholds <- sort(unique(c(
  seq(50, 100, by = 5),
  u95
)))

threshold_comparison <- map_dfr(
  candidate_thresholds,
  function(u) {
    event_peaks <- independent_peaks(u)
    
    tibble(
      threshold = u,
      events = nrow(event_peaks),
      events_per_year = nrow(event_peaks) / years_observed,
      target_met = between(
        nrow(event_peaks) / years_observed, 3, 5
      )
    )
  }
)

print(threshold_comparison)

ggplot(
  threshold_comparison,
  aes(x = threshold, y = events_per_year)
) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 3, ymax = 5,
    fill = "grey80", alpha = 0.3
  ) +
  geom_line() +
  geom_point() +
  geom_vline(
    xintercept = u95,
    linetype = "dashed"
  ) +
  labs(
    title = "Independent extreme-flow events by threshold",
    x = "Threshold (m3/s)",
    y = "Independent events per year"
  ) +
  theme_minimal()

# ------------------------------------------------------------
# Choosing a final threshold
# ------------------------------------------------------------
#For example

final_u <- 65

# Recreate the time-series vectors in sorted order
q <- river_data$Discharge
t <- river_data$Timestamp


# independent_peaks() from Lewis's code.
# Re-run its function definition after sorting river_data.

events <- independent_peaks(final_u) %>%
  arrange(Timestamp)

stopifnot(nrow(events) >= 2)

# Calculating time to peak 

# Searching up to 7 days before the peak, without passing
# the previous selected peak.


peak_rows <- match(events$Timestamp, t)

previous_rows <- c(1L, head(peak_rows, -1))

find_start <- function(p, previous) {
  
  window <- which(
    t >= max(t[p] - 7 * 24 * 3600, t[previous]) &
      t <= t[p]
  )
  
  if (length(window) == 0) {
    return(NA_integer_)
  }
  
  low <- min(q[window], na.rm = TRUE)
  
  cutoff <- low + 0.05 * (q[p] - low)
  
  eligible <- window[q[window] <= cutoff]
  
  if (length(eligible) == 0) {
    return(NA_integer_)
  }
  
  max(eligible)
}

start_rows <- mapply(
  find_start,
  peak_rows,
  previous_rows
)

events <- events %>%
  mutate(
    start_time = t[start_rows],
    max_flow = Discharge,
    time_to_peak_h = as.numeric(
      difftime(Timestamp, start_time, units = "hours")
    )
  )

# Check for invalid results
print(
  events %>%
    filter(
      is.na(time_to_peak_h) |
        time_to_peak_h <= 0 |
        time_to_peak_h > 168
    )
)


#  Descriptive statistics
valid_events <- events %>%
  filter(
    is.finite(max_flow),
    is.finite(time_to_peak_h),
    time_to_peak_h > 0
  )

summarise_parameter <- function(x) {
  tibble(
    n = length(x),
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    IQR = IQR(x),
    minimum = min(x),
    p05 = as.numeric(quantile(x, 0.05)),
    p25 = as.numeric(quantile(x, 0.25)),
    p75 = as.numeric(quantile(x, 0.75)),
    p95 = as.numeric(quantile(x, 0.95)),
    maximum = max(x)
  )
}

max_flow_summary <- summarise_parameter(
  valid_events$max_flow
)

time_to_peak_summary <- summarise_parameter(
  valid_events$time_to_peak_h
)

cat("\nMAXIMUM FLOW SUMMARY\n")
print(max_flow_summary)

cat("\nTIME TO PEAK SUMMARY\n")
print(time_to_peak_summary)


#Distrbutions of these parameters
#Max flow
p_flow <- ggplot(
  valid_events,
  aes(x = max_flow)
) +
  geom_histogram(
    bins = 15,
    fill = "steelblue",
    colour = "white"
  ) +
  geom_vline(
    xintercept = median(valid_events$max_flow),
    linetype = "dashed"
  ) +
  labs(
    title = "Maximum flow across extreme events",
    x = "Maximum flow (m3/s)",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_flow)

# Time to peak
p_time <- ggplot(
  valid_events,
  aes(x = time_to_peak_h)
) +
  geom_histogram(
    bins = 15,
    fill = "red3",
    colour = "white"
  ) +
  geom_vline(
    xintercept = median(valid_events$time_to_peak_h),
    linetype = "dashed"
  ) +
  labs(
    title = "Time to peak across extreme events",
    x = "Time to peak (hours)",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_time)


# Investigating the relationship between the two parameters
p_relationship <- ggplot(
  valid_events,
  aes(x = max_flow, y = time_to_peak_h)
) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "loess",
    se = FALSE
  ) +
  labs(
    title = "Maximum flow and time to peak",
    x = "Maximum flow (m3/s)",
    y = "Time to peak (hours)"
  ) +
  theme_minimal()

print(p_relationship)

# Spearman rank correlation
association <- cor.test(
  valid_events$max_flow,
  valid_events$time_to_peak_h,
  method = "spearman",
  exact = FALSE
)

print(association)


# Varation on the time to peak within flow groups

# Splits events into three equal-sized rank groups to assess whether events with similar flow
# magnitudes have different times to peak.

flow_groups <- valid_events %>%
  mutate(
    flow_group = ntile(max_flow, 3),
    flow_group = factor(
      flow_group,
      levels = 1:3,
      labels = c(
        "Lower extreme flows",
        "Middle extreme flows",
        "Higher extreme flows"
      )
    )
  )

group_summary <- flow_groups %>%
  group_by(flow_group) %>%
  summarise(
    events = n(),
    median_flow = median(max_flow),
    median_time = median(time_to_peak_h),
    IQR_time = IQR(time_to_peak_h),
    .groups = "drop"
  )

print(group_summary)

ggplot(
  flow_groups,
  aes(x = flow_group, y = time_to_peak_h)
) +
  geom_boxplot(fill = "lightblue") +
  labs(
    title = "Time to peak at different flow magnitudes",
    x = "Maximum-flow group",
    y = "Time to peak (hours)"
  ) +
  theme_minimal()

# Seasonal effects
seasonal_events <- valid_events %>%
  mutate(
    month = lubridate::month(Timestamp),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      TRUE ~ "Autumn"
    ),
    season = factor(
      season,
      levels = c("Winter", "Spring", "Summer", "Autumn")
    )
  )

seasonal_summary <- seasonal_events %>%
  group_by(season) %>%
  summarise(
    events = n(),
    median_flow = median(max_flow),
    median_time_to_peak = median(time_to_peak_h),
    IQR_time_to_peak = IQR(time_to_peak_h),
    .groups = "drop"
  )

print(seasonal_summary)

ggplot(
  seasonal_events,
  aes(x = season, y = time_to_peak_h)
) +
  geom_boxplot(fill = "green4") +
  labs(
    title = "Seasonal variation in time to peak",
    x = "Season",
    y = "Time to peak (hours)"
  ) +
  theme_minimal()


# Sensitivity of our threshold

# (Similar table to Lewis's)
# This function repeats the event extraction and time-to-peak calculation for each candidate threshold.

analyse_threshold <- function(u) {
  
  ep <- independent_peaks(u) %>%
    arrange(Timestamp)
  
  if (nrow(ep) < 3) {
    return(tibble(
      threshold = u,
      events = nrow(ep),
      events_per_year = nrow(ep) / ep$water_year,
      median_flow = NA_real_,
      median_time = NA_real_,
      IQR_time = NA_real_,
      spearman_rho = NA_real_
    ))
  }
  
  rows <- match(ep$Timestamp, t)
  previous <- c(1L, head(rows, -1))
  
  starts <- mapply(find_start, rows, previous)
  
  ep <- ep %>%
    mutate(
      max_flow = Discharge,
      time_to_peak_h = as.numeric(
        difftime(Timestamp, t[starts], units = "hours")
      )
    ) %>%
    filter(
      is.finite(max_flow),
      is.finite(time_to_peak_h),
      time_to_peak_h > 0
    )
  
  tibble(
    threshold = u,
    events = nrow(ep),
    events_per_year = nrow(ep) / ep$water_year,
    median_flow = median(ep$max_flow),
    median_time = median(ep$time_to_peak_h),
    IQR_time = IQR(ep$time_to_peak_h),
    spearman_rho = cor(
      ep$max_flow,
      ep$time_to_peak_h,
      method = "spearman"
    )
  )
}

sensitivity <- map_dfr(
  c(55, 60, 65, 70, 75),
  analyse_threshold
)

print(sensitivity)

ggplot(
  sensitivity,
  aes(x = threshold, y = median_time)
) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Threshold sensitivity of time to peak",
    x = "Threshold (m3/s)",
    y = "Median time to peak (hours)"
  ) +
  theme_minimal()

ggplot(
  sensitivity,
  aes(x = threshold, y = spearman_rho)
) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Threshold sensitivity of association",
    x = "Threshold (m3/s)",
    y = "Spearman correlation"
  ) +
  theme_minimal()



#------------------------------------------------------
#Threshold Comparison 

# Plot 95
plot_95 <- ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
  geom_line(color = "steelblue", linewidth = 0.4) +
  geom_hline(yintercept = u95, linetype = "dashed", color = "red") +
  geom_point(data = peaks_95, color = "red", size = 1) +
  labs(
    title = "Peaks above the 95th percentile threshold",
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()


thresholds <- 0:max(river_data$Discharge)

# Peaks per year for each threshold
peaks_per_year <- sapply(thresholds, function(u) {
  above <- river_data$Discharge > u
  sum(above & !lag(above, default = FALSE))
}) / n_distinct(river_data$water_year)


# Find all thresholds giving 3-5 peaks per year
selected_thresholds <- thresholds[
  peaks_per_year >= 3 & peaks_per_year <= 5
]

# Identify the peak for each exceedance period at each threshold
peaks_35_all <- purrr::map_dfr(selected_thresholds, function(u) {
  
  river_data %>%
    arrange(Timestamp) %>%
    mutate(
      above = Discharge > u,
      start = above & !lag(above, default = FALSE),
      run = cumsum(start)
    ) %>%
    filter(above) %>%
    group_by(run) %>%
    slice_max(Discharge, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(threshold = u)
  
})

# Categorise the peaks by discharge
peaks_coloured <- peaks_35_all %>%
  distinct(Timestamp, Discharge) %>%
  mutate(
    Peak_category = if_else(
      Discharge >= min(selected_thresholds) &
        Discharge <= max(selected_thresholds),
      "Within threshold band",
      "Above threshold band"
    )
  )

# Plot all thresholds on the hydrograph
plot_35 <- ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
  geom_line(color = "steelblue", linewidth = 0.4) +
  geom_point(
    data = peaks_coloured,
    aes(color = Peak_category),
    size = 1.5
  )  +
  annotate(
    "rect",
    xmin = min(river_data$Timestamp),
    xmax = max(river_data$Timestamp),
    ymin = min(selected_thresholds),
    ymax = max(selected_thresholds),
    fill = "red",
    alpha = 0.15
  ) + scale_color_manual(
    values = c(
      "Within threshold band" = "red",
      "Above threshold band" = "red4"
    ),
    name = "Peak classification"
  ) +
  labs(
    title = "Peaks using the 3-5 peaks per year method",
    subtitle = "All qualifying thresholds (73-113 m³/s)",
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()

# GPD diagnostic threshold
# 80 was used for the exploratory GPD fit
u_gpd <- 80

# Identify peaks above 80 m³/s
peaks_gpd <- river_data %>%
  arrange(Timestamp) %>%
  mutate(
    above = Discharge > u_gpd,
    start = above & !lag(above, default = FALSE),
    run = cumsum(start)
  ) %>%
  filter(above) %>%
  group_by(run) %>%
  slice_max(Discharge, n = 1, with_ties = FALSE) %>%
  ungroup()

# Common plot settings
base_plot <- ggplot(
  river_data,
  aes(x = Timestamp, y = Discharge)
) +
  geom_line(color = "steelblue", linewidth = 0.35) +
  labs(
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal(base_size = 12)

# GPD exploratory threshold
plot_gpd <- base_plot +
  geom_hline(
    yintercept = u_gpd,
    color = "red",
    linetype = "dashed"
  ) +
  geom_point(
    data = peaks_gpd,
    color = "red",
    size = 1.3
  ) +
  labs(
    title = "Peaks using the Extreme-Value Theory (GPD)",
    subtitle = "Exploratory threshold = 80 m³/s"
  )

plot_35 / plot_95 / plot_gpd

#--------------------------------------------------------------------
# Scenario based analyse

ggplot(events,
       aes(x = time_to_peak_h,
           y = Discharge)) +
  geom_point(
    color = "steelblue",
    alpha = 0.7,
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "red"
  ) +
  labs(
    title = " Variation in Flood Hydrograph Parameters",
    subtitle = "River Don at Parkhill",
    x = "Time to peak (hours)",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()

ggplot(events,
       aes(x = max_flow,
           y = Discharge)) +
  geom_point(
    color = "steelblue",
    alpha = 0.7,
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "red"
  ) +
  labs(
    title = " Variation in Flood Hydrograph Parameters",
    subtitle = "River Don at Parkhill",
    x = "Max flow",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()



