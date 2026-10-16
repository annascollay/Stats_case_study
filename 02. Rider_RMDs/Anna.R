library(tidyverse)

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
river_data %>%
  mutate(
    interval_minutes = as.numeric(
      difftime(Timestamp, lag(Timestamp), units = "mins")
    )
  ) %>%
  count(interval_minutes, sort = TRUE)

