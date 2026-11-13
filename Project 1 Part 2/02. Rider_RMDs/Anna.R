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
all_hydrograph_shapes <- ggplot(
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

Variation_hydrograph <- ggplot(
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

selected_event_min <- peaks %>% group_by(event_id) %>% slice_min(time_to_peak_h)

Comparison_hydrograph <- hydrographs %>%
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
    y = "Discharge",
    colour = "Event ID"
  ) +
  scale_color_manual(values = c("#AA1945", "#5298BC", "#BFDC80"))+
  theme_minimal()

df_events_long %>% filter(event_id %in% c(3,13, 24)) %>% 
  ggplot(aes(x = hours_from_peak,
             y = Discharge,
             colour = factor(event_id))) +
  geom_line(linewidth = 0.8) +
  labs(title = "Comparison of selected flood events",
       subtitle = "A comparision of the smallest, median and biggest discharge events",
       x = "Hours relative to peak",
       y = "Discharge",
       colour = "Event ID") +
  scale_color_manual(values = c("#AA1945", "#5298BC", "#BFDC80"))+
  theme_minimal()



# -------------------


#Heat map comparing flood events
heat_map_plot <- ggplot(df_events_long,
                        aes(x = hours_from_peak,
                            y = event_id,
                            fill = Discharge)) +
  geom_tile() +
  scale_fill_gradient(
    low = "#DCEEF5",
    high = "#174B67"
  ) +
  labs(
    title = "Heatmap of Flood Hydrograph Events",
    x = "Hours Relative to Peak",
    y = "Flood Event",
    fill = "Discharge"
  ) +
  theme_minimal()


# Comparasion of hydrographs with simialr peak flows shapes
# Identify the peak time (t = 0)
peak_col <- which.min(abs(rel_time_hours))

# Observed peak discharge for each event
peak_values <- Y_raw[, peak_col]

stopifnot(
  N_events >= 2,
  all(is.finite(Y_raw)),
  all(peak_values > 0),
  length(peak_values) == N_events
)

# Peak-normalized hydrographs
Y_norm <- sweep(Y_raw, 1, peak_values, "/")

# Long-format dataset
df_norm <- tibble(
  event_id = rep(seq_len(N_events), each = ncol(Y_norm)),
  hours_from_peak = rep(rel_time_hours, times = N_events),
  normalized_flow = as.vector(t(Y_norm))
)

# Select 4 events with similar peak discharge
n_compare <- min(4, N_events)

peak_order <- order(peak_values)
sorted_peaks <- peak_values[peak_order]

# Find the group of 4 consecutive sorted peaks
# with the smallest difference between largest and smallest
peak_ranges <- sapply(
  seq_len(N_events - n_compare + 1),
  function(i) {
    sorted_peaks[i + n_compare - 1] - sorted_peaks[i]
  }
)

best_start <- which.min(peak_ranges)

selected_events <- peak_order[
  best_start:(best_start + n_compare - 1)
]

# Prepare data
df_similar <- df_norm %>%
  filter(event_id %in% selected_events) %>%
  mutate(
    event = factor(
      event_id,
      levels = selected_events,
      labels = paste0(
        "Event ", selected_events,
        " (", round(peak_values[selected_events], 1),
        " m³/s)"
      )
    )
  )

similar_peak_plot <- ggplot(
  df_similar,
  aes(
    x = hours_from_peak,
    y = normalized_flow,
    colour = event
  )
) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  labs(
    title = "Hydrograph Shape Variation at Similar Peak Flows",
    subtitle = "Four events selected for similar maximum discharge",
    x = "Hours relative to peak",
    y = "Discharge / Peak discharge",
    colour = "Flood event"
  ) + scale_color_manual(values = c("#8ECFE0","#5298BC","#247C91", "#2A5267")) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(similar_peak_plot)