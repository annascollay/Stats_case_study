library(tidyverse)
library(ismev)
library(lubridate)

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

## Implementing indpendence Criterion
