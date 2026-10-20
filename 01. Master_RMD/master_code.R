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
    from_date <- paste0(yr, "-08-31T00:00:00")
    to_date   <- paste0(yr + 1, "-08-31T00:00:00")
    
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
    water_year = if_else(month(Timestamp) >= 9, year(Timestamp) + 1, year(Timestamp)),
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
    month_order = (month_num - 9) %% 12   # Sep = 0, Oct = 1, ..., Aug = 11
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


# Months in water year order (September to August)
wy_months <- month.abb[c(9:12, 1:8)]


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







