library(tidyverse)
library(lubridate)

# From Part 1's river_data and peaks tidying

# Checking the number of events
nrow(peaks)

# Extracting 5 days before and after each peak
hydrographs <- map_dfr(seq_len(nrow(peaks)), function(i) {
  
  peak_time <- peaks$Timestamp[i]
  peak_flow <- peaks$max_flow[i]
  
  river_data %>%
    filter(
      Timestamp >= peak_time - days(5),
      Timestamp <= peak_time + days(5)
    ) %>%
    transmute(
      event_id = i,
      Timestamp,
      hours_from_peak = as.numeric(
        difftime(Timestamp, peak_time, units = "hours")
      ),
      Discharge,
      peak_flow = peak_flow,
      scaled_flow = Discharge / peak_flow
    )
})

# Finding the number of observations for each of the 41 events
hydrographs %>%
  count(event_id)


# For the data (10-day window and 15-minute observations) 
# each complete hydrograph should contain 961 readings, including both endpoints.
#Checking this is the case for all 41 events

# The expected number of observations
expected_n <- 10 * 24 * 4 + 1

# Check number of readings and missing time steps
event_quality <- hydrographs %>%
  group_by(event_id) %>%
  summarise(
    n_readings = n(),
    missing_flow = sum(is.na(Discharge)),
    unique_times = n_distinct(hours_from_peak),
    .groups = "drop"
  ) %>%
  mutate(
    complete = n_readings == expected_n &
      missing_flow == 0 &
      unique_times == expected_n
  )

print(event_quality)


# Events requiring investigation
# Since Part 1 code removes readings with bad quality codes
event_quality %>%
  filter(!complete)


# Identifying the missing timestamps in incomplete events
incomplete_events <- event_quality %>%
  filter(!complete) %>%
  pull(event_id)

missing_readings <- map_dfr(incomplete_events, function(i) {
  
  peak_time <- peaks$Timestamp[i]
  
  # Expected timestamps over the 10-day window
  expected_times <- seq(
    from = peak_time - 5 * 24 * 3600,
    to = peak_time + 5 * 24 * 3600,
    by = "15 min"
  )
  
  # Actual timestamps in this event
  actual_times <- hydrographs %>%
    filter(event_id == i) %>%
    pull(Timestamp)
  
  # Find missing timestamps
  missing_times <- expected_times[
    !expected_times %in% actual_times
  ]
  
  tibble(
    event_id = i,
    peak_time = peak_time,
    missing_timestamp = missing_times,
    hours_from_peak = as.numeric(
      difftime(missing_times, peak_time, units = "hours")
    )
  )
})

print(missing_readings, n = Inf)

# Finding the quality codes of the missing obervations and reason for exclusion
original_data <- readRDS("river_don_parkhill_flow.rds")

missing_quality <- missing_readings %>%
  left_join(
    original_data %>%
      select(Timestamp, Discharge, QualityCode),
    by = c("missing_timestamp" = "Timestamp")
  ) %>%
  mutate(
    reason = case_when(
      is.na(QualityCode) ~ "Absent from original dataset",
      as.numeric(QualityCode) == 254 ~ "Excluded: Unchecked",
      as.numeric(QualityCode) %in% c(50, 100) ~
        "Acceptable quality - investigate",
      TRUE ~ "Excluded: Other quality code"
    )
  )

print(missing_quality, n = Inf)

# Visulasing the missing obervations in the data
ggplot(
  hydrographs %>%
    filter(event_id %in% incomplete_events),
  aes(x = hours_from_peak, y = scaled_flow)
) +
  geom_line(colour = "#5298BC", linewidth = 0.7) +
  geom_vline(
    data = missing_readings,
    aes(xintercept = hours_from_peak),
    colour = "red",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    colour = "black",
    linetype = "dotted"
  ) +
  facet_wrap(~event_id, ncol = 1) +
  labs(
    title = "Incomplete flood hydrographs",
    subtitle = "Red dashed lines indicate missing observations",
    x = "Hours relative to peak",
    y = "Normalised discharge"
  ) +
  theme_minimal()

#Checking to see how accurate it is if we estimate these gaps
# Compare original unchecked values with interpolation

interpolation_check <- missing_quality %>%
  group_by(event_id) %>%
  group_modify(~ {
    
    i <- .y$event_id
    
    # Clean observations surrounding this event
    event_data <- hydrographs %>%
      filter(event_id == i) %>%
      arrange(Timestamp)
    
    # Interpolate at missing timestamps
    estimated <- approx(
      x = as.numeric(event_data$Timestamp),
      y = event_data$Discharge,
      xout = as.numeric(.x$missing_timestamp),
      method = "linear"
    )$y
    
    .x %>%
      mutate(
        interpolated_flow = estimated,
        difference = Discharge - interpolated_flow,
        absolute_difference = abs(difference),
        percentage_difference =
          100 * absolute_difference / interpolated_flow
      )
  }) %>%
  ungroup()

# Display results
interpolation_check %>%
  select(
    event_id,
    missing_timestamp,
    Discharge,
    interpolated_flow,
    difference,
    percentage_difference
  ) %>%
  print(n = Inf)


#Visulasing the difference between our esrtimations and actual data
ggplot(interpolation_check,
       aes(x = missing_timestamp)) +
  geom_point(
    aes(y = Discharge, colour = "Original unchecked"),
    size = 3
  ) +
  geom_point(
    aes(y = interpolated_flow, colour = "Interpolated"),
    size = 3,
    shape = 17
  ) +
  facet_wrap(~event_id, scales = "free") +
  scale_colour_manual(
    values = c(
      "Original unchecked" = "#AA1945",
      "Interpolated" = "#5298BC"
    )
  ) +
  labs(
    title = "Original unchecked vs interpolated discharge",
    x = "Timestamp",
    y = "Discharge (m³/s)",
    colour = "Value"
  ) +
  theme_minimal()


# Exploratory analysis
#Plotting all hydrographs shapes
ggplot(
  hydrographs,
  aes(
    x = hours_from_peak,
    y = scaled_flow,
    group = event_id
  )
) +
  geom_line(alpha = 0.5, colour = "#5298BC") +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "#AA1945"
  ) +
  labs(
    title = "Variation in extreme flood hydrograph shapes",
    subtitle = "The 41 independent events aligned at peak flow",
    x = "Hours relative to peak",
    y = "Discharge / Maximum discharge"
  ) +
  theme_minimal()


# Mean hydrograph and variation plot 
complete_ids <- event_quality %>%
  filter(complete) %>%
  pull(event_id)

hydro_complete <- hydrographs %>%
  filter(event_id %in% complete_ids)

shape_summary <- hydro_complete %>%
  group_by(hours_from_peak) %>%
  summarise(
    mean_flow = mean(scaled_flow),
    lower = quantile(scaled_flow, 0.10),
    upper = quantile(scaled_flow, 0.90),
    .groups = "drop"
  )

ggplot(
  shape_summary,
  aes(x = hours_from_peak, y = mean_flow)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#5298BC",
    alpha = 0.25
  ) +
  geom_line(
    colour = "#174B67",
    linewidth = 1
  ) +
  labs(
    title = "Mean normalised hydrograph",
    subtitle = "Shaded area: 10th to 90th percentile",
    x = "Hours relative to peak",
    y = "Normalised discharge"
  ) +
  theme_minimal()


# Comparing three event peaks;
# the shortest, median and longest time to peak from part 1 estimates
selected_events <- peaks %>%
  mutate(event_id = row_number()) %>%
  arrange(time_to_peak_h) %>%
  slice(1, ceiling(n() / 2), n()) %>%
  pull(event_id)

hydrographs %>%
  filter(event_id %in% selected_events) %>%
  ggplot(
    aes(
      x = hours_from_peak,
      y = scaled_flow,
      colour = factor(event_id)
    )
  ) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Comparison of selected flood events",
    subtitle = "A comparision of the shortest, median and longest time to peak events",
    x = "Hours relative to peak",
    y = "Normalised discharge",
    colour = "Event ID"
  ) +
  scale_color_manual(values = c("#AA1945", "#5298BC", "#BFDC80"))+
  theme_minimal()

