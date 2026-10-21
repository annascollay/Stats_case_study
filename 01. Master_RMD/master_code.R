# ============================================================
# Code for Project 1
# ============================================================

# Load required libraries
library(tidyverse)
library(knitr)
library(kableExtra)
library(ismev)


#-----------------------------------------------------------
# LOAD THE DATA
#-----------------------------------------------------------
# Extend R's default HTTP timeout limit to 300 seconds
options(timeout = 300)

cache_file <- "river_don_parkhill_flow.rds"

if (file.exists(cache_file)) {
  # Load from local cache if file already exists
  river_data <- readRDS(cache_file)
  
} else {
  # If no local cache exists, fetch from SEPA in 1-year chunks
  start_years <- seq(2016, 2025, by = 1)
  
  fetch_chunk <- function(yr) {
    from_date <- paste0(yr, "-10-01T00:00:00")
    to_date   <- paste0(yr + 1, "-09-30T23:59:59")
    
    url <- paste0(
      "https://timeseries.sepa.org.uk/KiWIS/KiWIS?",
      "service=kisters&type=queryServices&datasource=0&",
      "request=getTimeseriesValues&ts_path=1/234292/Q/15m.Cmd&",
      "from=", from_date, "&to=", to_date, "&",
      "metadata=true&returnfields=Timestamp,Value,Quality%20Code&format=csv"
    )
    
    read_delim(
      url, 
      delim = ";", 
      comment = "#",
      col_names = c("Timestamp", "Discharge", "QualityCode"),
      show_col_types = FALSE
    )
  }
  
  river_data <- map_dfr(start_years, fetch_chunk) %>%
    mutate(Timestamp = as.POSIXct(Timestamp, format = "%Y-%m-%dT%H:%M:%S")) %>%
    filter(!is.na(Discharge)) %>%
    distinct(Timestamp, .keep_all = TRUE)
  
  # Save to local file so future runs don't need to use the API
  saveRDS(river_data, cache_file)
}


#-----------------------------------------------------------
# TIDY THE DATA
#-----------------------------------------------------------

# Data segmented by quality code
# Create the quality code table from SEPA guidance
quality_code_table <- tibble(
  `Code number` = c(50, 100, 140, 150, 200, 254),
  `Code name` = c(
    "Good",
    "Estimated",
    "Provisional",
    "Suspect",
    "Unchecked (imported from legacy database)",
    "Unchecked"
  ),
  `Abbreviated code name` = c("G", "E", "Prov", "S", "V", "U")
)

quality_combined <- quality_code_table %>%
  left_join(
    river_data %>%
      mutate(QualityCode = as.numeric(QualityCode)) %>%
      count(QualityCode, name = "Count"),
    by = c("Code number" = "QualityCode")
  ) %>%
  mutate(
    Count = coalesce(Count, 0L),
    Percentage = sprintf("%.2f%%", Count / sum(Count) * 100)
  )


#-----------------------------------------------------------
# EXPLORATORY ANALYSIS
#-----------------------------------------------------------

# Time variables added
river_data <- river_data %>%
  mutate(
    month      = month(Timestamp, label = TRUE),
    water_year = if_else(month(Timestamp) >= 10, year(Timestamp) + 1, year(Timestamp)),
    season     = case_when(
      month(Timestamp) %in% c(12, 1, 2) ~ "Winter",
      month(Timestamp) %in% c(3, 4, 5)  ~ "Spring",
      month(Timestamp) %in% c(6, 7, 8)  ~ "Summer",
      TRUE                              ~ "Autumn"
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))
  )

# Yearly analysis
yearly <- river_data %>%
  group_by(water_year) %>%
  summarise(
    Mean    = mean(Discharge),
    Median  = median(Discharge),
    Maximum = max(Discharge)
  )


# Mean discharge for each month
monthly <- river_data %>%
  mutate(month_num = as.numeric(format(Timestamp, "%m"))) %>%
  group_by(month_num) %>%
  summarise(value = mean(Discharge)) %>%
  mutate(
    month  = month.abb[month_num],
    season = case_when(
      month_num %in% c(9, 10, 11) ~ "Autumn",
      month_num %in% c(12, 1, 2)  ~ "Winter",
      month_num %in% c(3, 4, 5)   ~ "Spring",
      TRUE                        ~ "Summer"
    ),
    season = factor(season, levels = c("Autumn", "Winter", "Spring", "Summer")),
    month_order = (month_num - 10) %% 12   # Oct = 0, Nov = 1, ..., Sep = 11
  ) %>%
  arrange(month_order)
monthly$id <- 1:nrow(monthly)

# Angle and alignment of each month label
monthly$angle <- 90 - 360 * (monthly$id - 0.5) / nrow(monthly)
monthly$hjust <- ifelse(monthly$angle < -90, 1, 0)
monthly$angle <- ifelse(monthly$angle < -90, monthly$angle + 180, monthly$angle)

# Line under each season and position of its name
base_data <- monthly %>%
  group_by(season) %>%
  summarise(start = min(id), end = max(id)) %>%
  mutate(title = (start + end) / 2)

# Grid lines every 10 up to 40
grid_breaks <- c(10, 20, 30, 40)
y_max <- max(c(monthly$value, 40))


# Months in water year order (October to September)
wy_months <- month.abb[c(10:12, 1:9)]


# High flow defintion
high_flow <- quantile(river_data$Discharge, 0.95)

# Seasonality added
seasonal <- river_data %>%
  group_by(water_year, season) %>%
  summarise(
    Mean    = mean(Discharge),
    Maximum = max(Discharge),
    .groups = "drop"
  ) %>%
  mutate(season = factor(season, levels = c("Autumn", "Winter", "Spring", "Summer")))

season_colours <- c("Winter" = "#5298BC", "Spring" = "#BFDC80", "Summer" = "#FEE060", 
                    "Autumn" = "#AA1945" )

#-----------------------------------------------------------------------------------
# THRESHOLD ANALYSIS
#-----------------------------------------------------------------------------------
# 3-5 per year method
# Thresholds to try: 0, 1, 2, ... up to the maximum flow
thresholds <- 0:max(river_data$Discharge)

# Peaks per year for each threshold
peaks_per_year <- sapply(thresholds, function(u) {
  above <- river_data$Discharge > u
  sum(above & !lag(above, default = FALSE))
}) / n_distinct(river_data$water_year)

# Thresholds giving 3 to 5 peaks per year
thresholds[peaks_per_year >= 3 & peaks_per_year <= 5]


# Percentile method
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


# Count the peaks above a given percentile threshold
count_peaks <- function(p) {
  u <- quantile(river_data$Discharge, p)
  river_data %>%
    mutate(above = Discharge > u,
           start = above & !lag(above, default = FALSE)) %>%
    summarise(n = sum(start)) %>%
    pull(n)
}


# mrl
q <- river_data$Discharge

# Mean residual life plot: look for the lowest threshold above which the line is roughly straight
mrl.plot(q)

# How many readings are above a candidate threshold (change 100 to the value you read off the MRL plot)
length(which(q > 80))

# Stability of the parameter estimates across a range of thresholds
gpd.fitrange(q, 50, 200)

# Fit the GPD at the chosen threshold and check the parameter estimates (scale, shape)
gpd_fit <- gpd.fit(q, 80, npy = 96 * 365.25)
gpd_fit$mle

# Diagnostic plots
gpd.diag(gpd_fit)


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


#Threshold plots

# Common plot settings
base_plot <- ggplot(
  river_data,
  aes(x = Timestamp, y = Discharge)
) +
  geom_line(color = "#5298BC", linewidth = 0.35) +
  labs(
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal(base_size = 12)


# 95th percentile threshold plot
plot_95 <- base_plot +
  geom_hline(yintercept = u95, linetype = "dashed", color = "#ACD35A") +
  geom_point(data = peaks_95, color = "#AA1945", size = 1) +
  labs(
    title = "Peaks above the 95th percentile threshold",
    subtitle = expression("95th percentile threshold: 57.878 "*m^3*"/s")
  ) +
  theme_minimal()


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

# Plot all 3-5 thresholds on the hydrograph
plot_35 <- base_plot + geom_point(
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
    fill = "#BFDC80",
    alpha = 0.3
  ) + scale_color_manual(
    values = c(
      "Within threshold band" = "#6E0C2A",
      "Above threshold band" = "#AA1945"
    ),
    name = "Peak classification"
  ) +
  labs(
    title = "Peaks using the 3-5 peaks per year method",
    subtitle = "All qualifying thresholds (73-113 m³/s)"
  ) +
  theme_minimal()

# GPD exploratory threshold
plot_gpd <- base_plot +
  geom_hline(
    yintercept = u_gpd,
    color = "#ACD35A",
    linetype = "dashed"
  ) +
  geom_point(
    data = peaks_gpd,
    color = "#AA1945",
    size = 1.3
  ) +
  labs(
    title = "Peaks using the Extreme-Value Theory (GPD)",
    subtitle = "Exploratory threshold = 80 m³/s"
  )

#Graph of all thresholds together
# Find the threshold closest to 3, 4 and 5 peaks per year
# Use the highest threshold meeting or exceeding each target rate
get_threshold <- function(target) {
  eligible <- thresholds[peaks_per_year >= target]
  
  if (length(eligible) == 0) {
    stop(paste("No threshold gives at least", target,
               "peaks per year."))
  }
  
  max(eligible)
}

u3 <- get_threshold(3)
u4 <- get_threshold(4)
u5 <- get_threshold(5)

# Create a table of the four threshold lines
threshold_lines <- tibble(
  Threshold = c(
    "95th percentile",
    "3 peaks/year",
    "4 peaks/year",
    "5 peaks/year"
  ),
  Value = c(u95, u3, u4, u5)
)


# Identify the peak of each exceedance period above
# the LOWEST of the four thresholds
u_min <- min(threshold_lines$Value)

peaks_combined <- river_plot %>%
  mutate(
    above = Discharge > u_min,
    start = above &
      !dplyr::lag(above, default = FALSE),
    run = cumsum(start)
  ) %>%
  filter(above) %>%
  group_by(run) %>%
  slice_max(
    Discharge,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Categorise each peak by exactly which thresholds it exceeds
peaks_combined <- peaks_combined %>%
  mutate(
    Peak_category = case_when(
      Discharge > max(threshold_lines$Value) ~
        "Above all thresholds",
      
      Discharge > u4 & Discharge > u5 & Discharge > u95 ~
        "Above 95th, 4-peak and 5-peak thresholds",
      
      Discharge > u5 & Discharge > u95 ~
        "Above 95th and 5-peak thresholds",
      
      Discharge > u95 ~
        "Above 95th percentile only",
      
      TRUE ~
        "Other threshold combinations"
    )
  )

# Combined hydrograph
plot_thresholds <- ggplot(
  river_plot,
  aes(x = Timestamp, y = Discharge)
) +
  
  geom_line(
    color = "#5298BC",
    linewidth = 0.35
  ) +
  
  geom_hline(
    data = threshold_lines,
    aes(
      yintercept = Value,
      color = Threshold
    ),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_point(
    data = peaks_combined,
    aes(color = Peak_category),
    size = 1.6
  ) +
  
  scale_color_manual(
    values = c(
      "95th percentile" = "pink",
      "3 peaks/year" = "#AA1945",
      "4 peaks/year" = "red",
      "5 peaks/year" = "pink3",
      "Above all thresholds" = "#AA1945",
      "Above 95th, 4-peak and 5-peak thresholds" = "red",
      "Above 95th and 5-peak thresholds" = "pink3",
      "Above 95th percentile only" = "pink",
      "Other threshold combinations" = "#555555"
    )
  ) +
  
  labs(
    title = "Comparison of Flood Peak Thresholds",
    subtitle = "95th percentile and 3, 4 and 5 peaks per year",
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)"),
    color = "Threshold / Peak category"
  ) +
  
  theme_minimal(base_size = 12)

#-----------------------------------------------------------------------------------
# INDEPENDECE 
#-----------------------------------------------------------------------------------
q <- river_data$Discharge
t <- river_data$Timestamp

independent_peaks <- function(u) {
  
  # Step 1: candidate peaks - the highest reading in each period above the threshold
  candidates <- river_data %>%
    mutate(row   = row_number(),
           above = Discharge > u,
           run   = cumsum(above & !lag(above, default = FALSE))) %>%
    filter(above) %>%
    group_by(run) %>%
    slice_max(Discharge, n = 1, with_ties = FALSE) %>%
    pull(row)
  
  # Step 2: average time to rise, calculated using at least five clean peaks exceeding the threshold
  # A clean peak has no other candidate peak in the 7 days before it
  days_since_last <- as.numeric(difftime(t[candidates], lag(t[candidates]), units = "days"))
  clean <- candidates[is.na(days_since_last) | days_since_last > 7]
  if (length(clean) < 5) stop("Fewer than five clean peaks - lower the threshold")
  
  # Time to rise = duration between a trough and its corresponding peak
  # Trough = last point before the peak still within 5% of the rise above the lowest flow in the previous 7 days
  rise_hours <- sapply(clean, function(p) {
    window <- which(t >= t[p] - 7 * 24 * 3600 & t <= t[p])
    low    <- min(q[window])
    trough <- max(window[q[window] <= low + 0.05 * (q[p] - low)])
    as.numeric(difftime(t[p], t[trough], units = "hours"))
  })
  
  min_gap_hours <- 3 * mean(rise_hours)
  
  # Step 3: check each pair of adjacent peaks against both rules
  kept <- candidates[1]
  for (p in candidates[-1]) {
    last   <- kept[length(kept)]                                  # the first peak of the pair
    gap    <- as.numeric(difftime(t[p], t[last], units = "hours"))
    trough <- min(q[last:p])                                      # discharge at the trough between them
    
    # Rule 1: at least three times the average time to rise apart
    # Rule 2: trough at least two-thirds lower than the first peak
    if (gap >= min_gap_hours & trough <= q[last] / 3) {
      kept <- c(kept, p)                 # independent: a new event
    } else if (q[p] > q[last]) {
      kept[length(kept)] <- p            # same event: keep the larger peak
    }
  }
  
  cat("Candidate peaks:", length(candidates), "\n")
  cat("Clean peaks used for time to rise:", length(clean), "\n")
  cat("Average time to rise (hours):", round(mean(rise_hours), 1), "\n")
  cat("Independent peaks:", length(kept), "\n")
  
  river_data[kept, ]
}

## Rechecking the threshold
# Removing dependent peaks lowers the number of peaks per year, 
#so we recheck which thresholds give 3 to 5 independent peaks per year.

check <- tibble(threshold = seq(50, 200, by = 10)) %>%
  mutate(peaks_per_year = sapply(threshold, function(u) nrow(independent_peaks(u))) /
           n_distinct(river_data$water_year))
check

## Final threshold
u <- 65   # change to your final threshold
peaks <- independent_peaks(u)

# Row number of each peak, and of the peak before it
peak_rows <- match(peaks$Timestamp, t)
previous  <- c(1, head(peak_rows, -1))

# Start of each rise: same definition as the trough in the time to rise,
# searching back up to 7 days but never past the previous peak
start_rows <- mapply(function(p, prev) {
  window <- which(t >= max(t[p] - 7 * 24 * 3600, t[prev]) & t <= t[p])
  low    <- min(q[window])
  max(window[q[window] <= low + 0.05 * (q[p] - low)])
}, peak_rows, previous)

# The two quantities SEPA asked for
peaks <- peaks %>%
  mutate(
    start_time     = t[start_rows],
    max_flow       = Discharge,
    time_to_peak_h = as.numeric(difftime(Timestamp, start_time, units = "hours"))
  )

peaks %>% select(start_time, Timestamp, max_flow, time_to_peak_h)
summary(peaks$time_to_peak_h)


#---------------------------------------------------------------------------
# UNCERTANITY ANALYSIS
#--------------------------------------------------------------------------

# Max flow
Max_flow_plot <- ggplot(peaks, aes(x = max_flow)) +
  geom_histogram(bins = 15, fill = "#5298BC", colour = "#2A5267") +
  labs(
    title = "Distribution of max flow",
    x = expression("Max flow ("*m^3*"/s)"),
    y = "Number of events"
  ) +
  theme_minimal()

# Time to peak
Time_to_peak_plot <- ggplot(peaks, aes(x = time_to_peak_h)) +
  geom_histogram(bins = 20, fill = "#5298BC", colour = "#2A5267") +
  labs(
    title = "Distribution of time to peak",
    x = "Time to peak (hours)",
    y = "Number of events"
  ) +
  theme_minimal()

Max_vs_time <- ggplot(peaks, aes(x = time_to_peak_h, y = max_flow)) +
  geom_point(color = "#5298BC") +
  geom_smooth(method = "lm", color = "#AA1945", fill = "#FBB1BA", alpha = 0.05) +
  labs(
    title = "Max flow against time to peak",
    x = "Time to peak (hours)",
    y = expression("Max flow ("*m^3*"/s)")
  ) +
  theme_minimal()

cor_result <- cor.test(peaks$time_to_peak_h, peaks$max_flow, method = "spearman", exact = FALSE)


# Analysis of uncertanity
# Number of events and typical values in each season
events_season_table <- peaks %>%
  group_by(season) %>%
  summarise(
    Events = n(),
    `Median max flow` = median(max_flow),
    `Median time to peak (h)` = median(time_to_peak_h)
  )

events_season_table <- knitr::kable(
  events_season_table,
  booktabs = TRUE,
  digits = 3,
  caption = "Summary of Max Flow and Time to Peak for each Season"
)

# Relationship, coloured by season
Max_vs_time_season <- ggplot(peaks, aes(x = time_to_peak_h, y = max_flow )) +
  geom_point(size = 2, aes(colour = season)) +
  scale_color_manual(values = c("Autumn" = "#AA1945", "Winter" = "#5298BC",
                                "Spring" = "#BFDC80", "Summer" = "#FEE060")) +
  labs(
    title = "Max flow against time to peak, by season",
    x = "Time to peak (hours)",
    y = expression("Max flow ("*m^3*"/s)"),
    color = "Season"
  ) +
  theme_minimal()

# Max flow by season
max_season_plot <- ggplot(peaks, aes(x = season, y = max_flow)) +
  geom_boxplot(fill = season_colours, alpha = 0.7) +
  labs(
    title = "Max flow by season",
    x = "Season",
    y = expression("Max flow ("*m^3*"/s)")
  ) +
  theme_minimal()

# Time to peak by season
time_season_plot <- peaks %>% group_by(season) %>% ggplot(aes(x = season, y = time_to_peak_h)) +
  geom_boxplot(fill = season_colours, alpha = 0.7) +
  labs(
    title = "Time to peak by season",
    x = "Season",
    y = "Time to peak (hours)"
  ) +
  theme_minimal()

# Uncertanity final table
results <- tibble()

for (u in c(55, 60, 65, 70, 75)) {
  peaks <- independent_peaks(u)
  
  # Start of each rise and time to peak (same code as before)
  peak_rows <- match(peaks$Timestamp, t)
  previous  <- c(1, head(peak_rows, -1))
  start_rows <- mapply(function(p, prev) {
    window <- which(t >= max(t[p] - 7 * 24 * 3600, t[prev]) & t <= t[p])
    low    <- min(q[window])
    max(window[q[window] <= low + 0.05 * (q[p] - low)])
  }, peak_rows, previous)
  time_to_peak_h <- as.numeric(difftime(peaks$Timestamp, t[start_rows], units = "hours"))
  
  # One row of results for this threshold
  results <- bind_rows(results, tibble(
    threshold             = u,
    events                = nrow(peaks),
    events_per_year       = nrow(peaks) / n_distinct(river_data$water_year),
    median_max_flow       = median(peaks$Discharge),
    median_time_to_peak_h = median(time_to_peak_h),
    spearman_rho          = cor(peaks$Discharge, time_to_peak_h, method = "spearman")
  ))
}


results_table <- knitr::kable(
  results,
  col.names = c("threshold", "events", "events_per_year" = "Events per Year",
  "median_max_flow" = "Median Max Flow",
  "median_time_to_peak_h" = "Median Time to Peak (Hours)", 
  "spearman_rho" = "Spearman Rho Corelation"),
  booktabs = T,
  digits = 3,
  caption = "Final Uncertanity Table"
) %>%
  kable_styling(latex_options = "striped")