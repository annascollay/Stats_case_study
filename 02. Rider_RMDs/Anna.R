library(tidyverse)
library(ismev)

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
ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
  geom_line(color = "steelblue", linewidth = 0.4) +
  labs(
    title = "Hydrograph: Don at Parkhill (15-min Interval)",
    subtitle = "August 2016 – August 2026",
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()

#Check the mean residual life plot
mrl.plot(river_data$Discharge) # Stable around 80-200 

# Checking the stability of parameter estimates
gpd.fitrange(river_data$Discharge, 80, 200)

# Trying 100-200 the stability of parameter estimates
gpd.fitrange(river_data$Discharge, 100, 200)

# Fit GPD
gpd_fit <- gpd.fit(river_data$Discharge, 80)

#Check parameter estimates
gpd_fit$mle

# Plot diagnostics
gpd.diag(gpd_fit)

# Trying 100
# Fit GPD
gpd_fit <- gpd.fit(river_data$Discharge, 100)

#Check parameter estimates
gpd_fit$mle

# Plot diagnostics
gpd.diag(gpd_fit)

