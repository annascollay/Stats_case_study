# ==============================================================================
# Statistical Case Studies - Project 1 (Part 2)
# Authors: Anna Scollay, Lewis Seath, Zacharias Onyejiaka
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup & Package Loading
# ------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(knitr)
library(patchwork) # For combining diagnostic plots if installed
library(wavelets)
library(splines)   # part of base R; provides bs()
library(kableExtra)

# Set HTTP timeout limit for SEPA API fetch if local RDS is absent
options(timeout = 300)
cache_file <- "river_don_parkhill_flow.rds"

# ------------------------------------------------------------------------------
# 2. Load Data, Keep Codes 50/100/200/254, Set Aside Suspect (150)
# ------------------------------------------------------------------------------
if (file.exists(cache_file)) {
  river_data <- readRDS(cache_file)
} else {
  start_years <- seq(2001, 2025, by = 1)
  
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
    read_delim(url, delim = ";", comment = "#",
               col_names = c("Timestamp", "Discharge", "QualityCode"),
               show_col_types = FALSE)
  }
  
  river_data <- map_dfr(start_years, fetch_chunk) %>%
    mutate(Timestamp = as.POSIXct(Timestamp, format = "%Y-%m-%dT%H:%M:%S")) %>%
    filter(!is.na(Discharge)) %>%
    distinct(Timestamp, .keep_all = TRUE)
  
  saveRDS(river_data, cache_file)
}

# Quality codes kept for modelling, and the Suspect code checked separately (section 5)
retained_codes <- c(50, 100, 200, 254)   # Good, Estimated, Unchecked (legacy), Unchecked
suspect_code   <- 150                    # Suspect: beyond the acceptable extrapolation limit

# Full record before quality filtering, kept so the excluded Suspect readings
# can be checked against the event windows (section 5b)
river_data_raw <- river_data %>%
  arrange(Timestamp) %>%
  mutate(QualityCode = as.numeric(QualityCode))

river_data <- river_data_raw %>%
  filter(QualityCode %in% retained_codes) %>%
  mutate(
    is_unchecked = QualityCode %in% c(200, 254),
    month = month(Timestamp, label = TRUE),
    water_year = if_else(month(Timestamp) >= 10, year(Timestamp) + 1, year(Timestamp)),
    season = case_when(
      month(Timestamp) %in% c(12, 1, 2) ~ "Winter",
      month(Timestamp) %in% c(3, 4, 5) ~ "Spring",
      month(Timestamp) %in% c(6, 7, 8) ~ "Summer",
      TRUE ~ "Autumn"
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))
  )

# Extract convenient vectors for peak extraction
q <- river_data$Discharge
t <- river_data$Timestamp
n_years <- n_distinct(river_data$water_year)

# Subsection 2.1
#-------------------------------------------------------------------------------

# Data Description

# Count the number of missing readings by comparing our river data against a
# complete 15-minute sequence
complete_sequence <- seq(min(river_data_raw$Timestamp), 
                         max(river_data$Timestamp), by = "15 min")
number_missing <- sum(!complete_sequence %in% river_data_raw$Timestamp)

# Create data description table
data_description_table <- tibble(
  Item = c("Station", "Source", "Measurement", "Catchment area",
           "Resolution", "Station type", "Period", "Readings",
           "Missing readings"),
  Value = c(
    "Don at Parkhill",
    "SEPA time series data service",
    "River discharge ($m^3/s$)",
    "1,273 $km^2$",
    "15 minutes",
    "Velocity-area",
    paste(format(min(river_data_raw$Timestamp), "%d %B %Y"), "–",
          format(max(river_data_raw$Timestamp), "%d %B %Y")),
    format(nrow(river_data_raw), big.mark = ","),
    format(number_missing, big.mark = ",")
  )
)

# Data Quality

# Create quality code table from SEPA guidance, recording which codes we use
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
  `Abbreviated code name` = c("G", "E", "Prov", "S", "V", "U"),
  Used = 
    c("retained", "retained", "excluded", "excluded", "excluded", "excluded")
)

# Join quality code table with the number of readings for each code
quality_combined <- quality_code_table %>%
  left_join(
    river_data_raw %>% count(QualityCode = as.numeric(QualityCode), name = "Count"),
    by = c("Code number" = "QualityCode")
  ) %>%
  mutate(
    Count = coalesce(Count, 0L),
    Percentage = sprintf("%.2f\\%%", Count / sum(Count) * 100)
  )

# Combine into one table, with one row per quality code below the data description
data_table <- bind_rows(
  data_description_table,
  quality_combined %>%
    transmute(
      Item  = paste0(
        `Code number`, " – ", `Code name`, " (", `Abbreviated code name`, ")"),
      Value = paste0(format(Count, big.mark = ",", trim = TRUE), " readings (",
                     Percentage, "), ", Used)
    )
)

# Combined table

# Create a table describing the data
tidy_data_table <- kable(
  data_table, col.names = c("", ""),
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Data used in this analysis, and the quality codes assigned to each reading",
  label = "tidyingthedata") %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    font_size = 9,
    full_width = FALSE,
    position = "left"
  ) %>%
  pack_rows("Data description", 1, 9) %>%
  pack_rows("Quality codes", 10, nrow(data_table))

# ------------------------------------------------------------------------------
# 3. Peak Identification (Robson & Reed 1999 Criteria, u = 65 m³/s)
# ------------------------------------------------------------------------------

# Helper: Find rising limb breakout (P1P1 logic)
rise_start_row <- function(p, earliest = t[p] - 7 * 24 * 3600) {
  window <- which(t >= earliest & t <= t[p])
  lowest <- min(q[window])
  max(window[q[window] <= 1.05 * lowest])
}

# Candidate peaks above threshold u
candidate_peaks <- function(u, min_duration_hours = 1) {
  min_readings <- min_duration_hours * 4 
  river_data %>%
    mutate(
      row = row_number(),
      above = Discharge > u,
      run = cumsum(above & !lag(above, default = FALSE))
    ) %>%
    filter(above) %>%
    group_by(run) %>%
    filter(n() >= min_readings) %>% 
    slice_max(Discharge, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull(row)
}

# Independence criteria application
independent_peaks <- function(u) {
  candidates <- candidate_peaks(u)
  
  clean <- candidates[sapply(candidates, function(p) {
    window_mask <- abs(as.numeric(difftime(t[candidates], t[p], units = "days"))) <= 7
    q[p] == max(q[candidates[window_mask]])
  })]
  
  rise_hours <- sapply(clean, function(p) {
    as.numeric(difftime(t[p], t[rise_start_row(p)], units = "hours"))
  })
  min_gap_hours <- 3 * mean(rise_hours)
  
  kept <- candidates[1]
  for (p in candidates[-1]) {
    last <- kept[length(kept)]
    gap <- as.numeric(difftime(t[p], t[last], units = "hours"))
    trough <- min(q[last:p])
    
    if (gap >= min_gap_hours & trough <= q[last] / 3) {
      kept <- c(kept, p)
    } else if (q[p] > q[last]) {
      kept[length(kept)] <- p
    }
  }
  
  peaks_u <- river_data[kept, ]
  attr(peaks_u, "mean_rise_hours") <- mean(rise_hours)
  peaks_u
}

# Isolate independent peaks at final threshold u = 65 m³/s
u_final <- 65
peaks <- independent_peaks(u_final) %>%
  mutate(event_id = row_number())

# ------------------------------------------------------------------------------
# 4. Matrix Construction & Diagnostic Output
# ------------------------------------------------------------------------------
win_before_h   <- 48
win_after_h    <- 72
dt_hours       <- 0.25 # 15-minute sampling interval
rel_time_hours <- seq(-win_before_h, win_after_h, by = dt_hours)

extract_event_hydrograph <- function(event_idx, peak_row) {
  t_pk <- t[peak_row]
  idx  <- which(t >= t_pk - win_before_h * 3600 & t <= t_pk + win_after_h * 3600)
  df_win <- river_data[idx, ]
  
  event_t_rel <- as.numeric(difftime(df_win$Timestamp, t_pk, units = "hours"))
  
  # Resample onto the standard 481-point relative time grid
  approx(x = event_t_rel, y = df_win$Discharge, xout = rel_time_hours, rule = 2)$y
}

N_events  <- nrow(peaks)
peak_rows <- match(peaks$Timestamp, t)

# Construct N x m matrix Y_raw (118 x 481)
Y_raw <- t(sapply(1:N_events, function(i) extract_event_hydrograph(i, peak_rows[i])))
rownames(Y_raw) <- paste0("Event_", 1:N_events)
colnames(Y_raw) <- paste0("t_", sprintf("%+.2fh", rel_time_hours))

# Calculate empirical mean hydrograph vector Q_bar(t)
mean_hydrograph <- colMeans(Y_raw)

# ------------------------------------------------------------------------------
# 5. Unchecked readings vs interpolation for peak windows
# ------------------------------------------------------------------------------

# Step 1: for every unchecked reading in the record, interpolate linearly
# between the checked readings immediately before and after its burst.
# Done on the full record so that a burst touching the edge of a peak window
# still has its checked anchor points.
compare_unchecked <- function(data, unchecked_codes = c(200, 254)) {
  d <- data %>%
    mutate(
      is_unch  = QualityCode %in% unchecked_codes,
      burst_id = cumsum(is_unch & !lag(is_unch, default = FALSE)),
      burst_id = if_else(is_unch, burst_id, NA_integer_)
    )
  tnum <- as.numeric(d$Timestamp)
  out  <- list()
  
  for (b in na.omit(unique(d$burst_id))) {
    bi         <- which(d$burst_id == b)
    idx_before <- min(bi) - 1
    idx_after  <- max(bi) + 1
    
    # No checked anchor on one side (start/end of record): cannot interpolate
    if (idx_before < 1 || idx_after > nrow(d)) next
    
    interp <- approx(
      x    = tnum[c(idx_before, idx_after)],
      y    = d$Discharge[c(idx_before, idx_after)],
      xout = tnum[bi]            # exact timestamps of the unchecked readings
    )$y
    
    out[[length(out) + 1]] <- tibble(
      Burst_ID            = b,
      Timestamp           = d$Timestamp[bi],
      Discharge_Recorded = d$Discharge[bi],
      Discharge_Interp    = interp,
      Diff                = d$Discharge[bi] - interp,
      Burst_Length_Mins   = length(bi) * 15
    )
  }
  bind_rows(out)
}

unchecked_all <- compare_unchecked(river_data)

# Flagged readings falling inside a set of windows. starts/ends/peak_times
# are POSIXct vectors with one element per event. readings defaults to the
# Unchecked readings; pass suspect_all for the Suspect readings.
unchecked_in_windows <- function(starts, ends, peak_times, readings = unchecked_all) {
  if (nrow(readings) == 0) return(tibble())
  map_dfr(seq_along(starts), function(i) {
    readings %>%
      filter(Timestamp >= starts[i], Timestamp <= ends[i]) %>%
      mutate(Event_ID = i,
             Hours_From_Peak = as.numeric(difftime(Timestamp, peak_times[i],
                                                   units = "hours")))
  })
}

# Statistics, combined table and scatter plot for one set of windows.
# label is used in the plot ("Unchecked" or "Suspect").
# Returns NULL if there are no flagged readings in the windows.
summarise_unchecked_readings <- function(in_windows, n_events, label = "Unchecked") {
  if (nrow(in_windows) == 0) return(NULL)
  
  # each reading counted once even if two event windows overlap
  unique_readings <- distinct(in_windows, Timestamp, .keep_all = TRUE)
  rmse   <- sqrt(mean(unique_readings$Diff^2))
  max_ad <- max(abs(unique_readings$Diff))
  max_q  <- max(unique_readings$Discharge_Recorded)
  
  events_lookup <- in_windows %>%
    group_by(Burst_ID) %>%
    summarise(Events = paste(sort(unique(Event_ID)), collapse = ", "),
              .groups = "drop")
  
  burst_tbl <- unique_readings %>%
    group_by(Burst_ID) %>%
    summarise(
      Start         = format(min(Timestamp), "%Y-%m-%d %H:%M"),
      Readings      = n(),
      Duration_Mins = first(Burst_Length_Mins),
      Max_Flow      = max(Discharge_Recorded),
      RMSE          = sqrt(mean(Diff^2)),
      Max_Abs_Diff  = max(abs(Diff)),
      .groups = "drop"
    ) %>%
    left_join(events_lookup, by = "Burst_ID") %>%
    transmute(
      Burst        = as.character(row_number()),
      Event        = Events,
      Start,
      Readings     = as.character(Readings),
      Duration_Min = as.character(Duration_Mins),
      Max_Flow     = sprintf("%.1f", Max_Flow),
      RMSE         = sprintf("%.3f", RMSE),
      Max_Abs_Diff = sprintf("%.3f", Max_Abs_Diff)
    )
  
  overall_row <- tibble(
    Burst        = "**All**",
    Event        = sprintf("**%d of %d events**", n_distinct(in_windows$Event_ID), n_events),
    Start        = "—",
    Readings     = sprintf("**%d**", nrow(unique_readings)),
    Duration_Min = "—",
    Max_Flow     = sprintf("**%.1f**", max_q),
    RMSE         = sprintf("**%.3f**", rmse),
    Max_Abs_Diff = sprintf("**%.3f**", max_ad)
  )
  
  lims <- range(c(unique_readings$Discharge_Interp, unique_readings$Discharge_Recorded))
  
  plot <- ggplot(unique_readings, aes(x = Discharge_Interp, y = Discharge_Recorded)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#AA1945") +
    geom_point(colour = "#5298BC", alpha = 0.7, size = 2) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(
      title = paste(label, "readings vs linear interpolation"),
      subtitle = sprintf("n = %d | RMSE = %.3f m³/s | dashed line = 1:1",
                         nrow(unique_readings), rmse),
      x = expression("Interpolated discharge (" * m^3 * "/s)"),
      y = paste(label, "discharge (m³/s)")
    ) +
    theme_minimal()
  
  list(in_windows = in_windows, unique_readings = unique_readings,
       rmse = rmse, max_ad = max_ad, max_q = max_q,
       burst_tbl = burst_tbl, overall_row = overall_row,
       table = bind_rows(burst_tbl, overall_row), plot = plot)
}

# Hydrographs of the events containing flagged readings. Line = readings used
# in the model (across a gap this is exactly the linear interpolation used);
# points = flagged readings. One panel per affected event.
plot_flagged_events <- function(flagged, starts, ends, peak_times, label,
                                colour = "#AA1945") {
  if (is.null(flagged)) return(NULL)
  affected <- sort(unique(flagged$in_windows$Event_ID))
  panel    <- paste0("Event ", affected, " (", format(peak_times[affected], "%d %b %Y"), ")")
  
  hydro <- map_dfr(seq_along(affected), function(k) {
    i <- affected[k]
    river_data %>%
      filter(Timestamp >= starts[i], Timestamp <= ends[i]) %>%
      mutate(Event = panel[k])
  }) %>% mutate(Event = factor(Event, levels = panel))
  
  pts <- flagged$in_windows %>%
    mutate(Event = factor(panel[match(Event_ID, affected)], levels = panel))
  
  ggplot(hydro, aes(Timestamp, Discharge)) +
    geom_line(colour = "#5298BC", linewidth = 0.5) +
    geom_point(data = pts, aes(y = Discharge_Recorded), colour = colour, size = 1.2) +
    geom_hline(yintercept = u_final, linetype = "dotted") +
    facet_wrap(~ Event, scales = "free") +
    labs(
      title    = paste(label, "readings within event hydrographs"),
      subtitle = sprintf("Line: readings used in the model; points: %s readings; dotted: u = %g m³/s",
                         tolower(label), u_final),
      x = NULL, y = expression("Discharge (" * m^3 * "/s)")
    ) +
    theme_minimal()
}

# Fixed [-48h, +72h] windows
unchecked_in_events <- unchecked_in_windows(
  starts     = peaks$Timestamp - win_before_h * 3600,
  ends       = peaks$Timestamp + win_after_h * 3600,
  peak_times = peaks$Timestamp
)
unchecked_fixed <- summarise_unchecked_readings(unchecked_in_events, N_events)

if (!is.null(unchecked_fixed)) {
  # Create object names for Rmd
  unique_readings        <- unchecked_fixed$unique_readings
  rmse                   <- unchecked_fixed$rmse
  max_ad                 <- unchecked_fixed$max_ad
  burst_tbl              <- unchecked_fixed$burst_tbl
  overall_row            <- unchecked_fixed$overall_row
  unchecked_table        <- unchecked_fixed$table
  unchecked_scatter_plot <- unchecked_fixed$plot
} else {
  cat("No unchecked quality codes detected within the 120-hour peak event windows.\n\n")
}

# Suspect readings are EXCLUDED from river_data. Wherever one falls inside an
# event, the model bridges the gap by linear interpolation between the
# retained readings either side, which is what compare_unchecked() computes.
# The same table and plots therefore show (a) whether any Suspect reading
# falls in an event window, (b) how large it is relative to u, and (c) how far
# the excluded value is from the interpolated value the model uses instead.

# Record-wide summary (for the report text, whatever the window results)
suspect_readings <- river_data_raw %>% filter(QualityCode == suspect_code)
n_suspect        <- nrow(suspect_readings)
pct_suspect      <- 100 * n_suspect / nrow(river_data_raw)
n_suspect_above  <- sum(suspect_readings$Discharge > u_final)
max_suspect      <- if (n_suspect > 0) max(suspect_readings$Discharge) else NA_real_

# Each Suspect reading against the interpolation between its retained neighbours
suspect_all <- compare_unchecked(
  river_data_raw %>% filter(QualityCode %in% c(retained_codes, suspect_code)),
  unchecked_codes = suspect_code
)
# Bursts at the very start/end of the record have no anchor and are skipped
cat(sprintf("Suspect readings: %d in record (%.2f%%), %d above u, max %.1f m3/s; %d compared with interpolation\n",
            n_suspect, pct_suspect, n_suspect_above, max_suspect, nrow(suspect_all)))

# Fixed [-48h, +72h] windows
suspect_in_events <- unchecked_in_windows(
  starts     = peaks$Timestamp - win_before_h * 3600,
  ends       = peaks$Timestamp + win_after_h * 3600,
  peak_times = peaks$Timestamp,
  readings   = suspect_all
)
suspect_fixed <- summarise_unchecked_readings(suspect_in_events, N_events, label = "Suspect")

if (is.null(suspect_fixed)) {
  cat("No Suspect readings within the 120-hour peak event windows.\n\n")
} else {
  cat(sprintf("Fixed windows: %d Suspect readings in %d events\n\n",
              nrow(suspect_fixed$unique_readings),
              n_distinct(suspect_fixed$in_windows$Event_ID)))
}

# ------------------------------------------------------------------------------
# 6. Master Hydrograph Overlay Plot Generation
# ------------------------------------------------------------------------------
df_events_long <- as_tibble(Y_raw) %>%
  mutate(event_id = factor(1:N_events)) %>%
  pivot_longer(
    cols = starts_with("t_"),
    names_to = "t_col",
    values_to = "Discharge"
  ) %>%
  mutate(hours_from_peak = rep(rel_time_hours, times = N_events))

df_mean <- tibble(
  hours_from_peak = rel_time_hours,
  Discharge = mean_hydrograph
)

# Plot: All 118 hydrographs centered at peak arrival (greyed out + bold black mean)
hydrograph_overlay_plot <- ggplot() +
  geom_line(
    data = df_events_long,
    aes(x = hours_from_peak, y = Discharge, group = event_id),
    colour = "#5298BC",
    linewidth = 0.4,
    alpha = 0.25
  ) +
  geom_line(
    data = df_mean,
    aes(x = hours_from_peak, y = Discharge),
    colour = "#174B67",
    linewidth = 1.3
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#AA1945", linewidth = 0.6) +
  geom_hline(yintercept = u_final, linetype = "dotted", colour = "black", linewidth = 0.7) +
  annotate("text", x = 0, y = max(Y_raw) * 0.98, label = "Peak Arrival (t = 0)", 
           colour = "#AA1945", hjust = -0.1, size = 3.5, fontface = "bold") +
  annotate("text", x = -45, y = u_final + 5, label = paste0("u = ", u_final, " m³/s"), 
           colour = "black", hjust = 0, size = 3.2) +
  scale_x_continuous(breaks = seq(-48, 72, by = 12)) +
  labs(
    title = "River Don Hydrograph Shape Variability (N = 118 Events)",
    subtitle = "Centered at peak arrival (t = 0) [-48h, +72h]. Bold black: Mean curve Q̄(t). Unchecked values retaine; Suspect (code 150) excluded.",
    x = "Hours from Peak Arrival (t = 0)",
    y = expression("Discharge (" * m^3 * "/s)")
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, colour = "grey30")
  )

# Display master plot
print(hydrograph_overlay_plot)

#Heat map comparing extreme flow events
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
    title = "Heatmap of Extreme Flow Event Hydrograph",
    x = "Hours Relative to Peak",
    y = "Extreme Flow Event",
    fill = "Discharge"
  ) +
  theme_minimal() +
  scale_y_discrete(breaks = levels(df_events_long$event_id)[seq(1, N_events, 10)])

# ------------------------------------------------------------------------------
# 7. Event start, end and width (trough to trough)
# ------------------------------------------------------------------------------

search_days   <- 7     # how far either side of a peak to search for a trough
trough_factor <- 1.05  # trough = reading within 5% of the local minimum

# End of the falling limb: the first reading after the peak within 5% of the
# lowest flow before the limit (mirrors rise_start_row)
fall_end_row <- function(p, limit) {
  window <- which(t >= t[p] & t <= limit)
  lowest <- min(q[window])
  min(window[q[window] <= trough_factor * lowest])
}

# Where to stop searching for each event's end: 7 days after its peak,
# or the next event's peak if that comes sooner
end_limit <- pmin(
  t[peak_rows] + search_days * 24 * 3600,
  lead(t[peak_rows], default = max(t))
)

# One row per event: start, peak and end, with the key quantities for each
events <- tibble(
  event_id  = seq_along(peak_rows),
  peak_row  = peak_rows,
  start_row = map_int(peak_rows, rise_start_row),
  end_row   = map2_int(peak_rows, end_limit, fall_end_row)
) %>%
  mutate(
    peak_time      = t[peak_row],
    start_time     = t[start_row],
    end_time       = t[end_row],
    peak_flow      = q[peak_row],
    time_to_peak_h = as.numeric(difftime(peak_time, start_time, units = "hours")),
    recession_h    = as.numeric(difftime(end_time, peak_time, units = "hours")),
    width_h        = time_to_peak_h + recession_h,
    peak_position  = time_to_peak_h / width_h   # where the peak falls in [0, 1]
  )

unchecked_model_in_events <- unchecked_in_windows(
  starts     = events$start_time,
  ends       = events$end_time,
  peak_times = events$peak_time
)
unchecked_model <- summarise_unchecked_readings(unchecked_model_in_events, nrow(events))

if (!is.null(unchecked_model)) {
  fixed_ts <- if (is.null(unchecked_fixed)) as.POSIXct(character(0)) else
    unchecked_fixed$unique_readings$Timestamp
  cat(sprintf(
    "Modelling events (trough to trough): %d unchecked readings, %d of them outside the fixed windows.\n",
    nrow(unchecked_model$unique_readings),
    sum(!unchecked_model$unique_readings$Timestamp %in% fixed_ts)))
  cat(sprintf("RMSE vs interpolation: %.3f m3/s | max abs diff: %.3f m3/s\n\n",
              unchecked_model$rmse, unchecked_model$max_ad))
  # unchecked_model$table  -> kable() for the report
  # unchecked_model$plot   -> scatter plot for the report
} else {
  cat("No unchecked quality codes within the modelling events.\n\n")
}

# Suspect readings within the modelling events (trough to trough)
suspect_model_in_events <- unchecked_in_windows(
  starts     = events$start_time,
  ends       = events$end_time,
  peak_times = events$peak_time,
  readings   = suspect_all
)
suspect_model      <- summarise_unchecked_readings(suspect_model_in_events, nrow(events),
                                                   label = "Suspect")
suspect_model_plot <- plot_flagged_events(suspect_model, events$start_time,
                                          events$end_time, events$peak_time, "Suspect")

# Quantities for the report text
n_suspect_model_readings <- if (is.null(suspect_model)) 0 else nrow(suspect_model$unique_readings)
n_suspect_model_events   <- if (is.null(suspect_model)) 0 else n_distinct(suspect_model$in_windows$Event_ID)
n_suspect_model_bursts   <- if (is.null(suspect_model)) 0 else nrow(suspect_model$burst_tbl)

if (is.null(suspect_model)) {
  cat("No Suspect readings within the modelling events.\n\n")
} else {
  cat(sprintf("Modelling events: %d Suspect readings (%d bursts) in %d events | max %.1f m3/s | RMSE vs interpolation %.3f m3/s\n\n",
              n_suspect_model_readings, n_suspect_model_bursts, n_suspect_model_events,
              suspect_model$max_q, suspect_model$rmse))
}

# ------------------------------------------------------------------------------
# 8. Standardise each event and make the Y matrix
# ------------------------------------------------------------------------------
# For each event, one at a time:
#   a. rescale time onto [0, 1]   (0 = start of event, 1 = end)
#   b. realign, so the peak moves to u_star (the average peak position),
#      keeping the start at 0 and the end at 1
#   c. divide flow by peak flow   (so every event peaks at 1)
#   d. read the curve off at 64 equally spaced points between 0 and 1
#      (linear interpolation between the 15-minute readings), giving one
#      row of Y

u_star <- mean(events$peak_position)
n_grid <- 64
grid_u <- seq(0, 1, length.out = n_grid)

standardise_event <- function(i) {
  rows <- events$start_row[i]:events$end_row[i]   # readings in this event
  p    <- events$peak_position[i]                 # where its peak falls in [0, 1]
  
  # a. rescale time onto [0, 1]
  u <- as.numeric(difftime(t[rows], events$start_time[i], units = "hours")) /
    events$width_h[i]
  
  # b. realign: rising part -> [0, u_star], falling part -> [u_star, 1]
  u_reg <- ifelse(u <= p,
                  u * u_star / p,
                  u_star + (u - p) * (1 - u_star) / (1 - p))
  
  # c. divide flow by peak flow, d. read off at the 64 grid points
  approx(u_reg, q[rows] / events$peak_flow[i], xout = grid_u, rule = 2)$y
}

Y <- t(sapply(seq_len(nrow(events)), standardise_event))

cat("Y:", nrow(Y), "events x", ncol(Y), "points; peaks aligned at u =", round(u_star, 3), "\n")

# ------------------------------------------------------------------------------
# 9. Reconstruct a real hydrograph
# ------------------------------------------------------------------------------
# Reverses section 5: turns a standardised curve back into a real extreme flow event.
#   1. undo the realignment: move the peak from u_star back to the event's
#      own peak position
#   2. undo the time rescaling: hours since start = u x width
#   3. undo the flow scaling: discharge (m^3/s) = y x peak flow

reconstruct_hydrograph <- function(y, peak_position, width_h, peak_flow) {
  
  # 1. undo the realignment
  u <- ifelse(grid_u <= u_star,
              grid_u * peak_position / u_star,
              peak_position + (grid_u - u_star) * (1 - peak_position) / (1 - u_star))
  
  # 2. and 3. back to real time and real flow
  data.frame(
    hours     = u * width_h,
    discharge = y * peak_flow
  )
}

# ------------------------------------------------------------------------------
# 10. Shared functions for checking every method
# ------------------------------------------------------------------------------
# Every method is checked in the same way, so they can be compared fairly:
# each fitted curve is turned back into a real extreme flow event and compared
# with the observed event in m^3/s and hours.

# Set use_cache <- FALSE to force everything to be recomputed.
use_cache <- TRUE

# name : label for the cache file
# expr : the expensive code (only evaluated when there is no usable cache)
# key  : any settings the result depends on besides Y (ranges, chosen values)
# The file name contains a hash of Y and key, so the cache refreshes
# automatically if the data or those settings change.
# If the CODE of a function used inside expr is edited, delete the cache /
# folder.
cached <- function(name, expr, key = NULL) {
  if (!use_cache) return(expr)
  file <- file.path("cache",
                    paste0(name, "_", substr(rlang::hash(list(Y, key)), 1, 10), ".rds"))
  if (file.exists(file)) return(readRDS(file))
  result <- expr
  tryCatch({
    dir.create("cache", showWarnings = FALSE)
    saveRDS(result, file)
  }, error = function(e) NULL, warning = function(w) NULL)   # never fail on saving
  result
}

# Trapezium rule used for extreme flow volume. Needed because
# after undoing the realignment the points are not evenly spaced in hours.
area_under <- function(hours, flow) {
  sum(diff(hours) * (head(flow, -1) + tail(flow, -1)) / 2)
}

# Errors for one extreme flow event:
#   rmse         : overall closeness of the curve (m^3/s)
#   peak_error   : fitted peak minus observed peak (m^3/s)
#   timing_error : fitted peak time minus observed peak time (hours)
#   volume_error : % difference in total volume
flood_errors <- function(hours, observed, fitted) {
  c(
    rmse         = sqrt(mean((fitted - observed)^2)),
    peak_error   = max(fitted) - max(observed),
    timing_error = hours[which.max(fitted)] - hours[which.max(observed)],
    volume_error = 100 * (area_under(hours, fitted) - area_under(hours, observed)) /
      area_under(hours, observed)
  )
}

# Errors for every extreme flow event, for a matrix of fitted curves laid out like Y
method_errors <- function(Y_fit) {
  errors <- matrix(NA, nrow = nrow(Y), ncol = 4)
  for (i in 1:nrow(Y)) {
    observed <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                       events$width_h[i], events$peak_flow[i])
    fitted   <- reconstruct_hydrograph(Y_fit[i, ], events$peak_position[i],
                                       events$width_h[i], events$peak_flow[i])
    errors[i, ] <- flood_errors(observed$hours, observed$discharge, fitted$discharge)
  }
  colnames(errors) <- c("rmse", "peak_error", "timing_error", "volume_error")
  errors
}

# Average size of each error across extreme flow events
average_errors <- function(errors) {
  colMeans(abs(errors))
}

# Leave-one-out cross-validation: each event is left out in turn,
# the method is fitted to the other events, and the left-out event is
# predicted. fit_left_out(i) must return the method's curve for event i
# using only the other events. u_star is the average over all events; its
# effect on the left-out event is negligible
loocv <- function(fit_left_out) {
  Y_pred <- Y
  for (i in 1:nrow(Y)) {
    Y_pred[i, ] <- fit_left_out(i)
  }
  Y_pred
}

# Fixed training and test events: about 20% held out, chosen once
# with a fixed seed so same split always used
set.seed(2026)
test_floods  <- sort(sample(1:nrow(Y), size = round(0.2 * nrow(Y))))
train_floods <- setdiff(1:nrow(Y), test_floods)

# Plot observed vs fitted events in real units for six kinds of event
sharpness <- rowMeans(Y > 0.5)   # share of the event spent above half the peak

example_events <- unique(c(which.min(events$peak_position), which.max(events$peak_position),
                     which.min(events$width_h),       which.max(events$width_h),
                     which.max(sharpness),            which.min(sharpness)))

plot_example_fits <- function(Y_fit, title) {
  par(mfrow = c(2, 3))
  for (i in example_events) {
    observed <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                       events$width_h[i], events$peak_flow[i])
    fitted   <- reconstruct_hydrograph(Y_fit[i, ], events$peak_position[i],
                                       events$width_h[i], events$peak_flow[i])
    plot(observed$hours, observed$discharge, type = "l", lwd = 2,
         ylim = range(observed$discharge, fitted$discharge),
         xlab = "Hours since start", ylab = "Discharge (m^3/s)",
         main = paste(title, "- event", i))
    lines(fitted$hours, fitted$discharge, col = "red", lwd = 2, lty = 2)
    legend("topright", c("Observed", "Fitted"), col = c("black", "red"),
           lty = c(1, 2), lwd = 2, cex = 0.8)
  }
  par(mfrow = c(1, 1))
}

# Plot residuals against time for every event. Red line is average residual.
plot_residuals <- function(Y_fit, title) {
  residuals <- Y - Y_fit
  matplot(grid_u, t(residuals), type = "l", lty = 1, col = "grey",
          xlab = "Rescaled time u", ylab = "Observed - fitted",
          main = paste(title, ": residuals against time"))
  lines(grid_u, colMeans(residuals), col = "red", lwd = 2)
  abline(h = 0, lty = 2)
}

# Plot the four average errors against a method's setting with one line for each type of error
plot_choice <- function(choice_table, setting_name) {
  types <- unique(choice_table$type)
  par(mfrow = c(2, 2))
  for (error in c("rmse", "peak_error", "timing_error", "volume_error")) {
    plot(range(choice_table$setting), range(choice_table[[error]]), type = "n",
         xlab = setting_name, ylab = "Average absolute error", main = error)
    for (j in seq_along(types)) {
      rows <- choice_table$type == types[j]
      lines(choice_table$setting[rows], choice_table[[error]][rows],
            type = "b", pch = 19, col = j)
    }
    legend("topright", types, col = seq_along(types), lty = 1, cex = 0.7)
  }
  par(mfrow = c(1, 1))
}

# Fitted curves for every event, for any per-event linear model
fit_all_with <- function(fit_model, setting) {
  t(sapply(1:nrow(Y), function(i) fitted(fit_model(Y[i, ], setting))))
}

# Train/test errors by setting
train_test_choice <- function(settings, fit_all) {
  map_dfr(settings, function(s) {
    errors <- method_errors(fit_all(s))
    rbind(
      data.frame(setting = s, type = "train", t(average_errors(errors[train_floods, ]))),
      data.frame(setting = s, type = "test",  t(average_errors(errors[test_floods, ])))
    )
  })
}

# In-sample vs leave-one-out errors by setting
insample_loo_choice <- function(settings, rebuild, left_out) {
  map_dfr(settings, function(s) {
    rbind(
      data.frame(setting = s, type = "in-sample",
                 t(average_errors(method_errors(rebuild(s))))),
      data.frame(setting = s, type = "leave-one-out",
                 t(average_errors(method_errors(loocv(function(i) left_out(i, s))))))
    )
  })
}

# Summary row for the final comparison table
cmp_row <- function(method, n_params, type, errors, negative = NA) {
  data.frame(method = method, n_params = n_params, type = type,
             t(average_errors(errors)), negative_points = negative)
}

# A standardised curve as an event in m^3/s and hours, using event i's own
# peak flow, width and peak time
as_flood <- function(curve, i) {
  reconstruct_hydrograph(curve, events$peak_position[i], events$width_h[i], events$peak_flow[i])
}

# ------------------------------------------------------------------------------
# 11. SEPA's current approach
# -----------------------------------------------------------------------------
# SEPA uses one average shape for every event, scaled by peak flow. SEPA does
# not realign peaks, so the average is taken over the events on their own
# timing. Every method is compared against this.

# Each event's curve on its own timing on the 64 points
unaligned_curve <- function(i) {
  real <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                 events$width_h[i], events$peak_flow[i])
  approx(real$hours / events$width_h[i], Y[i, ], xout = grid_u, rule = 2)$y
}

Y_unaligned <- t(sapply(1:nrow(Y), unaligned_curve))

# Average shape given to event i, laid out like Y so it can be checked in the same way as every method
baseline_curve <- function(average_shape, i) {
  real <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                 events$width_h[i], events$peak_flow[i])
  approx(grid_u, average_shape, xout = real$hours / events$width_h[i], rule = 2)$y
}

# In-sample: every event gets the average of all events
average_shape <- colMeans(Y_unaligned)
average_shape <- average_shape / max(average_shape)   # rescale so the peak equals the max flow, as SEPA does
Y_baseline    <- t(sapply(1:nrow(Y), function(i) baseline_curve(average_shape, i)))

baseline_errors <- method_errors(Y_baseline)

# Leave-one-out: every event gets the average of the other events
baseline_errors_loocv <- cached("baseline_loocv",
                                method_errors(loocv(function(i) {
                                  s <- colMeans(Y_unaligned[-i, ])
                                  baseline_curve(s / max(s), i)
                                })))

cat("\nSEPA's approach (average shape), average errors:\n")
print(rbind(in_sample     = average_errors(baseline_errors),
            leave_one_out = average_errors(baseline_errors_loocv)))


# ------------------------------------------------------------------------------
# 12. Method 1: polynomial fitting
# ------------------------------------------------------------------------------
# Each event's curve is fitted with a polynomial in time, separately for each event, using lm():
#   y(u) = theta_0 + theta_1 u + theta_2 u^2 + ... + theta_p u^p + error
# The fitted coefficients theta_0, ..., theta_p are the event's shape
# parameters. Each event is fitted on its own, so leave-one-out errors are the
# same as in-sample errors. degree is checked on training and test events.

# Formula y ~ I(u^1) + I(u^2) + ... + I(u^p), written as in the handout
poly_formula <- function(p) {
  as.formula(paste("y ~", paste0("I(u^", 1:p, ")", collapse = " + ")))
}

# Fit one event's curve with a degree-p polynomial
fit_poly_model <- function(y, p) {
  lm(poly_formula(p), data = data.frame(y = y, u = grid_u))
}

# Fitted curves for every event at degree p
poly_fit_all <- function(p) fit_all_with(fit_poly_model, p)

# Choosing the degree
degrees <- 2:15

# Significance of the highest term: proportion of events where the highest term has p < 0.05.
# Neighbouring points on a curve are correlated, so these p-values are too optimistic
prop_significant <- cached("poly_prop_sig", sapply(degrees, function(p) {
                    mean(sapply(1:nrow(Y), function(i) {
                    coefs <- summary(fit_poly_model(Y[i, ], p))$coefficients
                    coefs[nrow(coefs), 4] < 0.05
                    }))}),
                    key = degrees
                    )

print(data.frame(degree = degrees, prop_significant = prop_significant))

# Errors by degree on training and test events. the chosen degree
# should work as well on the test events as on the training events
poly_choice <- cached("poly_choice",
                      train_test_choice(degrees, poly_fit_all), key = degrees)

print(poly_choice, digits = 3)
plot_choice(poly_choice, "Polynomial degree")

# Fit at the chosen degree
p_chosen <- 7

poly_models <- lapply(1:nrow(Y), function(i) fit_poly_model(Y[i, ], p_chosen))

# Fitted curves, one row per event
Y_poly <- t(sapply(poly_models, fitted))

# Shape parameters, one row per event
poly_coefs <- t(sapply(poly_models, coef))
colnames(poly_coefs) <- paste0("theta_", 0:p_chosen)

poly_errors <- method_errors(Y_poly)

cat("\nPolynomial (degree", p_chosen, "), average errors:\n")
print(average_errors(poly_errors))

# Fitted flow should never be negative
cat("Grid points with negative fitted flow:", sum(Y_poly < 0), "\n")

# Checking the fit
plot_example_fits(Y_poly, "Polynomial")
plot_residuals(Y_poly, "Polynomial")

# Residuals against each polynomial term
residuals_poly <- as.vector(t(Y - Y_poly))
u_all          <- rep(grid_u, times = nrow(Y))

par(mfrow = c(2, 3))
for (k in 1:p_chosen) {
  plot(u_all^k, residuals_poly, pch = 20, cex = 0.4, col = "grey40",
       xlab = paste0("u^", k), ylab = "Residual", main = paste0("Residuals vs u^", k))
  abline(h = 0, col = "red", lty = 2)
}
par(mfrow = c(1, 1))

# Standard diagnostic plots for the worst-fitted event:
# residuals vs fitted, Q-Q, scale-location, leverage. Normality and
# independence matter less here: the aim is a good basis, not inference.
worst_poly <- which.max(poly_errors[, "rmse"])

par(mfrow = c(2, 2))
plot(poly_models[[worst_poly]])
par(mfrow = c(1, 1))

# Coefficient table for the same event
print(summary(poly_models[[worst_poly]]))

# Correlation between shape parameters across events: strong
# correlations make modelling their uncertainty harder later
cat("\nCorrelation between polynomial shape parameters:\n")
print(round(cor(poly_coefs), 2))


# ------------------------------------------------------------------------------
# 13. Method 2: principal component analysis (PCA)
# ------------------------------------------------------------------------------
# PCA finds the main shapes in the extreme flow events from the data itself.
# Each event's curve is rebuilt as the mean curve plus k component curves:
#   y_i(u) = mean curve(u) + C_i1 a_1(u) + ... + C_ik a_k(u)
#   a_j  : loadings of component j (a curve over the 64 grid points)
#   C_ij : score of event i on component j
# The scores C_i1, ..., C_ik are the event's shape parameters.
# Curves are centred but not scaled, as all points are on the same scale

make_pca <- function(data, log = FALSE) {
  X    <- if (log) log(data) else data      # log version needs strictly positive flow
  back <- if (log) exp else identity
  fit  <- prcomp(X, center = TRUE, scale. = FALSE)
  
  # Rebuilt curves for all events using the first k components
  rebuild <- function(k) {
    scores   <- fit$x[, 1:k, drop = FALSE]
    loadings <- fit$rotation[, 1:k, drop = FALSE]
    back(sweep(scores %*% t(loadings), 2, colMeans(X), FUN = "+"))
  }
  
  # event i rebuilt from components found without it
  left_out <- function(i, k) {
    fit_train  <- prcomp(X[-i, ], center = TRUE, scale. = FALSE)
    mean_train <- colMeans(X[-i, ])
    loadings   <- fit_train$rotation[, 1:k, drop = FALSE]
    scores     <- (X[i, ] - mean_train) %*% loadings
    back(as.vector(mean_train + scores %*% t(loadings)))
  }
  
  list(fit = fit, rebuild = rebuild, left_out = left_out)
}

pca             <- make_pca(Y)
pca_result      <- pca$fit
pca_rebuild     <- pca$rebuild
pca_left_out    <- pca$left_out

# Components and variance explained
# Loadings of the first three components as curves over time
plot(grid_u, pca_result$rotation[, 1], type = "l", col = "red",
     ylim = range(pca_result$rotation[, 1:3]),
     xlab = "Rescaled time u", ylab = "Loading",
     main = "Loadings for PC1, PC2 and PC3")
lines(grid_u, pca_result$rotation[, 2], col = "blue")
lines(grid_u, pca_result$rotation[, 3], col = "darkgreen")
abline(h = 0, lty = 3)
legend("bottomright", c("PC1", "PC2", "PC3"),
       col = c("red", "blue", "darkgreen"), lty = 1)

# Variance of each component
plot(pca_result, main = "Variance of each principal component")

# Cumulative proportion of variance explained: choose the smallest k that captures enough
cum_var <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)

plot(seq_along(cum_var), cum_var, type = "b", pch = 19, ylim = c(0, 1),
     xlab = "Number of components k", ylab = "Cumulative proportion of variance",
     main = "Cumulative variance explained")
abline(h = 0.95, lty = 2, col = "red")

cat("Smallest k explaining at least 95% of variance:", which(cum_var >= 0.95)[1], "\n")

# Choosing k: In-sample errors always fall as k grows; leave-one-out errors show whether
# extra components still help for a new event
k_values <- 1:min(10, nrow(Y) - 2)

pca_choice <- cached("pca_choice",
                     insample_loo_choice(k_values, pca_rebuild, pca_left_out),
                     key = k_values)

print(pca_choice, digits = 3)
plot_choice(pca_choice, "Number of components k")

# Fit at the chosen k and check the fit.
# Update after looking at sections 9a and 9b

k_chosen <- 6

# Rebuilt curves, one row per event
Y_pca <- pca_rebuild(k_chosen)

# Shape parameters (scores), one row per event
pca_scores <- pca_result$x[, 1:k_chosen, drop = FALSE]

pca_errors       <- method_errors(Y_pca)
pca_errors_loocv <- cached("pca_loocv",
                           method_errors(loocv(function(i) pca_left_out(i, k_chosen))),
                           key = k_chosen)

cat("\nPCA (k =", k_chosen, "), average errors:\n")
print(rbind(in_sample     = average_errors(pca_errors),
            leave_one_out = average_errors(pca_errors_loocv)))

# Rebuilt flow should never be negative
cat("Grid points with negative rebuilt flow:", sum(Y_pca < 0), "\n")

plot_example_fits(Y_pca, "PCA")
plot_residuals(Y_pca, "PCA")

# Outliers in the scores: events far from the rest are unusual
pairs(pca_scores, pch = 19, main = "Scores on the kept components")

# Variance and correlation of the scores. Scores are
# uncorrelated by construction, which helps when modelling their uncertainty.
cat("\nVariance of each PCA score:\n")
print(round(apply(pca_scores, 2, var), 4))
cat("\nCorrelation between PCA scores:\n")
print(round(cor(pca_scores), 2))


# ------------------------------------------------------------------------------
# 14. Method 3: wavelets
# ------------------------------------------------------------------------------
# Each event's curve is split by the discrete wavelet transform into a smooth
# trend plus detail at several scales. Most detail coefficients
# are tiny, so only the largest few are kept and the rest set to zero.
# The same detail coefficients are kept for every event, chosen
# as those that are largest on average across events.
# Daubechies "la8" filter, 3 levels, periodic boundaries. Y has 64 = 2^6 points, so no padding is needed.
# Shape parameters: the smooth coefficients plus the kept detail coefficients.

wav_filter <- "la8"   # wavelet filter
wav_levels <- 3      # number of levels

# Wavelet transform of one curve
wavelet_transform <- function(y, filter = wav_filter, levels = wav_levels) {
  dwt(as.numeric(y), filter = filter, n.levels = levels, boundary = "periodic")
}

# Average size of each detail coefficient across a set of events, by level
average_detail_size <- function(wav_list) {
  lapply(1:length(wav_list[[1]]@W), function(j) {
    rowMeans(sapply(wav_list, function(wt) abs(wt@W[[j]][, 1])))
  })
}

# Rebuild one curve keeping only the "keep" detail coefficients that are
# largest on average (avg_size), setting the rest to zero
wavelet_rebuild <- function(wt, avg_size, keep) {
  cutoff <- sort(unlist(avg_size), decreasing = TRUE)[keep]
  for (j in seq_along(avg_size)) {
    wt@W[[j]][avg_size[[j]] < cutoff, 1] <- 0
  }
  as.numeric(idwt(wt))
}

# Wavelet transform of every event
wav_all  <- lapply(1:nrow(Y), function(i) wavelet_transform(Y[i, ]))
avg_size <- average_detail_size(wav_all)

# Rebuilt curves for every event, keeping "keep" detail coefficients
wavelet_fit_all <- function(keep) {
  t(sapply(wav_all, function(wt) wavelet_rebuild(wt, avg_size, keep)))
}

# Left-out event i, keeping the coefficients chosen from the other events
wavelet_left_out <- function(i, keep) {
  wavelet_rebuild(wav_all[[i]], average_detail_size(wav_all[-i]), keep)
}

# Decomposition of one event
wt_example <- wav_all[[1]]

par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
plot(wt_example@W[[1]], type = "l", col = "red",       xlab = "", ylab = "", main = "Level 1 detail (fine)")
plot(wt_example@W[[2]], type = "l", col = "blue",      xlab = "", ylab = "", main = "Level 2 detail")
plot(wt_example@W[[3]], type = "l", col = "darkgreen", xlab = "", ylab = "", main = "Level 3 detail (coarse)")
plot(wt_example@V[[3]], type = "l", col = "black",     xlab = "", ylab = "", main = "Level 3 smooth (trend)")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Keeping every coefficient should rebuild the curve exactly
cat("Max error of full rebuild:", max(abs(Y[1, ] - as.numeric(idwt(wt_example)))), "\n")

# Choosing how many detail coefficients to keep
keep_values <- c(1, 2, 3, 5, 8, 10, 15, 20)

wav_choice <- cached("wav_choice",
                     insample_loo_choice(keep_values, wavelet_fit_all, wavelet_left_out),
                     key = keep_values)

print(wav_choice, digits = 3)
plot_choice(wav_choice, "Detail coefficients kept")

# Fit at the chosen number and check the fit
# Update after looking at section 10b
keep_chosen <- 5

# Rebuilt curves, one row per event
Y_wav <- wavelet_fit_all(keep_chosen)

# Shape parameters: smooth coefficients + kept detail coefficients
cutoff_chosen     <- sort(unlist(avg_size), decreasing = TRUE)[keep_chosen]
kept_coefs       <- unlist(avg_size) >= cutoff_chosen
wav_params <- t(sapply(wav_all, function(wt) {
  c(wt@V[[wav_levels]][, 1], unlist(wt@W)[kept_coefs])
}))

cat("Wavelet shape parameters per event:", ncol(wav_params), "\n")

wav_errors       <- method_errors(Y_wav)
wav_errors_loocv <- cached("wav_loocv",
                           method_errors(loocv(function(i) wavelet_left_out(i, keep_chosen))),
                           key = keep_chosen)

cat("\nWavelets (", keep_chosen, "detail coefficients ), average errors:\n")
print(rbind(in_sample     = average_errors(wav_errors),
            leave_one_out = average_errors(wav_errors_loocv)))

# Rebuilt flow should never be negative
cat("Grid points with negative rebuilt flow:", sum(Y_wav < 0), "\n")

plot_example_fits(Y_wav, "Wavelets")
plot_residuals(Y_wav, "Wavelets")

# Sensitivity to the wavelet filter and number of levels, keeping
# the same number of detail coefficients. More levels means fewer smooth
# coefficients (64 / 2^levels), so fewer parameters.
run_wav_sensitivity <- function() {
  grid <- expand.grid(wf = c("haar", "d4", "la8"), wl = 3:5, stringsAsFactors = FALSE)
  map_dfr(seq_len(nrow(grid)), function(r) {
    wf <- grid$wf[r]; wl <- grid$wl[r]
    wav_list <- lapply(1:nrow(Y), function(i) wavelet_transform(Y[i, ], wf, wl))
    sizes    <- average_detail_size(wav_list)
    Y_fit    <- t(sapply(wav_list, function(wt) wavelet_rebuild(wt, sizes, keep_chosen)))
    data.frame(filter = wf, levels = wl, n_params = 64 / 2^wl + keep_chosen,
               t(average_errors(method_errors(Y_fit))))
  })
}

wav_sensitivity <- cached("wav_sensitivity", run_wav_sensitivity(),
                          key = list(keep_chosen, wav_filter, wav_levels))

print(wav_sensitivity, digits = 3)

# Variance and correlation of the shape parameters
cat("\nVariance of each wavelet shape parameter:\n")
print(round(apply(wav_params, 2, var), 4))
cat("\nCorrelation between wavelet shape parameters:\n")
print(round(cor(wav_params), 2))

# ------------------------------------------------------------------------------
# 15. Method 4: known hydrograph function (gamma-type)
# ------------------------------------------------------------------------------
# Fit a known function that approximates hydrographs,
# as its parameters can be interpreted. A standard choice is the gamma-type
# hydrograph. Our curves start and end at the trough level rather than zero,
# so a baseflow term b is added:
#   y(u) = b + (1 - b) (u / u_p)^m exp(m (1 - u / u_p))
#   b   : baseflow level, as a fraction of peak flow
#   u_p : time of the peak on [0, 1]
#   m   : sharpness (larger m = narrower, sharper peak)
# The curve equals 1 at u_p. It is not linear in its parameters, so it is
# fitted by nonlinear least squares with nls(), separately for each event.
# Each event is fitted on its own, so leave-one-out errors are the same as
# in-sample errors.

gamma_curve <- function(u, b, up, m) {
  b + (1 - b) * (u / up)^m * exp(m * (1 - u / up))
}

# Fit one event's curve
fit_gamma_model <- function(y) {
  tryCatch(
    nls(y ~ gamma_curve(u, b, up, m),
        data      = data.frame(y = y, u = grid_u),
        start     = list(b = max(min(y), 0.01), up = u_star, m = 5),
        algorithm = "port",
        lower     = c(b = 0,    up = 0.01, m = 0.1),
        upper     = c(b = 0.99, up = 0.99, m = 500)),
    error = function(e) NULL
  )
}

gamma_models <- lapply(1:nrow(Y), function(i) fit_gamma_model(Y[i, ]))

# events where the fit worked
gamma_ok <- which(!sapply(gamma_models, is.null))
cat("Gamma fits that failed to converge:", nrow(Y) - length(gamma_ok), "\n")

# Fitted curves, one row per event (failed fits left as the observed curve, and excluded from the errors below)
Y_gamma <- Y
for (i in gamma_ok) {
  Y_gamma[i, ] <- fitted(gamma_models[[i]])
}

# Shape parameters, one row per event
gamma_params <- matrix(NA, nrow = nrow(Y), ncol = 3,
                       dimnames = list(NULL, c("b", "up", "m")))
for (i in gamma_ok) {
  gamma_params[i, ] <- coef(gamma_models[[i]])
}

gamma_errors <- method_errors(Y_gamma)[gamma_ok, ]

cat("\nGamma-type function, average errors (converged fits):\n")
print(average_errors(gamma_errors))

# Fitted flow should never be negative
cat("Grid points with negative fitted flow:", sum(Y_gamma[gamma_ok, ] < 0), "\n")

plot_example_fits(Y_gamma, "Gamma-type")
plot_residuals(Y_gamma, "Gamma-type")

# The parameters can be interpreted: their spread and relationships
cat("\nSummary of gamma-type parameters:\n")
print(summary(gamma_params))

pairs(gamma_params[gamma_ok, ], pch = 19, main = "Gamma-type parameters across events")

cat("\nCorrelation between gamma-type parameters:\n")
print(round(cor(gamma_params[gamma_ok, ]), 2))


# ------------------------------------------------------------------------------
# 16. Method 5: B-splines
# ------------------------------------------------------------------------------
# Ramsay & Silverman's standard alternative to polynomials. A B-spline basis
# joins short cubic pieces smoothly, so each basis function only affects one
# part of the curve. This avoids the polynomial's wiggles in the tail, while
# still being a linear model fitted with lm(), like the polynomial:
#   y(u) = theta_0 + theta_1 B_1(u) + ... + theta_df B_df(u) + error
# df (the number of basis functions) plays the role of the polynomial degree.
# Each event is fitted on its own, so the number of basis functions is
# checked on training and test events, as for the polynomial.

# Fit one event's curve with n_basis B-spline basis functions
fit_bspline_model <- function(y, n_basis) {
  lm(y ~ bs(u, df = n_basis), data = data.frame(y = y, u = grid_u))
}

# Fitted curves for every event with n_basis basis functions
bspline_fit_all <- function(n_basis) fit_all_with(fit_bspline_model, n_basis)

# Choosing the number of basis functions
n_basis_values <- 4:15

bspline_choice <- cached("bspline_choice",
                         train_test_choice(n_basis_values, bspline_fit_all),
                         key = n_basis_values)

print(bspline_choice, digits = 3)
plot_choice(bspline_choice, "Number of B-spline basis functions")

# Fit at the chosen number and check the fit
# Update after doing section 13a

n_basis_chosen <- 8

bspline_models <- lapply(1:nrow(Y), function(i) fit_bspline_model(Y[i, ], n_basis_chosen))

# Fitted curves, one row per event
Y_bspline <- t(sapply(bspline_models, fitted))

# Shape parameters, one row per event
bspline_coefs <- t(sapply(bspline_models, coef))
colnames(bspline_coefs) <- paste0("theta_", 0:n_basis_chosen)

bspline_errors <- method_errors(Y_bspline)

cat("\nB-splines (", n_basis_chosen, "basis functions ), average errors:\n")
print(average_errors(bspline_errors))

# Fitted flow should never be negative
cat("Grid points with negative fitted flow:", sum(Y_bspline < 0), "\n")

plot_example_fits(Y_bspline, "B-spline")
plot_residuals(Y_bspline, "B-spline")

# Standard diagnostic plots for the worst-fitted event
worst_bspline <- which.max(bspline_errors[, "rmse"])

par(mfrow = c(2, 2))
plot(bspline_models[[worst_bspline]])
par(mfrow = c(1, 1))

# Correlation between shape parameters across events
cat("\nCorrelation between B-spline shape parameters:\n")
print(round(cor(bspline_coefs), 2))

# ------------------------------------------------------------------------------
# 18. Validation of methods
# ------------------------------------------------------------------------------

# Fitted curves from every method on same six contrasting events

# Same six events as plot_example_fits() (section 10), with all five methods
# drawn on each event so they can be compared directly.

fits_to_compare <- list(Y_poly, Y_pca, Y_wav, Y_gamma, Y_bspline)
method_names    <- c("Polynomial", "PCA", "Wavelets", "Gamma-type", "B-spline")
method_colours <- c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7")

# Fitted curves from every method on the same six events as plot_example_fits()
# (section 10). One row per method, one column per event.
plot_all_methods <- function() {
  par(mfrow = c(5, length(example_events)), mar = c(2, 2, 2, 0.5), oma = c(3, 3, 0, 0))
  for (j in 1:5) {
    for (i in example_events) {
      observed <- as_flood(Y[i, ], i)
      plot(observed$hours, observed$discharge, type = "l", lwd = 1.5, col = "black",
           xlab = "", ylab = "", main = paste(method_names[j], "- event", i),
           cex.main = 0.9)
      # skip the gamma-type curve if its fit failed for this event
      if (j == 4 && !(i %in% gamma_ok)) next
      fitted <- as_flood(fits_to_compare[[j]][i, ], i)
      lines(fitted$hours, fitted$discharge, col = method_colours[j], lwd = 2)
    }
  }
  mtext("Hours since start",   side = 1, outer = TRUE, line = 1)
  mtext("Discharge (m^3/s)",   side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
}

# Residuals against time for every method on same figure

plot_all_residuals <- function() {
  # residuals for each method (only events where the gamma-type fit worked)
  residuals_by_method <- list()
  for (j in 1:5) {
    rows <- if (j == 4) gamma_ok else 1:nrow(Y)
    residuals_by_method[[j]] <- Y[rows, ] - fits_to_compare[[j]][rows, ]
  }
  y_range <- range(unlist(residuals_by_method))   # same y-axis on every panel
  
  par(mfrow = c(2, 3))
  # one panel per method: every event faint, average bold
  for (j in 1:5) {
    matplot(grid_u, t(residuals_by_method[[j]]), type = "l", lty = 1,
            col = adjustcolor("grey40", alpha.f = 0.3), ylim = y_range,
            xlab = "Rescaled time u", ylab = "Observed - fitted", main = method_names[j])
    lines(grid_u, colMeans(residuals_by_method[[j]]), col = method_colours[j], lwd = 2.5)
    abline(h = 0, lty = 2)
  }
  # sixth panel: the five averages together, on their own scale
  averages <- sapply(residuals_by_method, colMeans)
  matplot(grid_u, averages, type = "l", lty = 1, lwd = 2, col = method_colours,
          xlab = "Rescaled time u", ylab = "Average residual", main = "Averages compared")
  abline(h = 0, lty = 2)
  legend("bottomright", method_names, col = method_colours, lwd = 2, cex = 0.7)
  par(mfrow = c(1, 1))
}


# ------------------------------------------------------------------------------
# 18. Comparing all methods
# ------------------------------------------------------------------------------
# Accuracy (average errors), parsimony (number of shape parameters) and the
# physical check (negative flows). For methods fitted to each event on its
# own (polynomial, gamma-type, B-spline), leave-one-out errors are the same
# as in-sample errors.

all_comparison <- rbind(
  cmp_row("SEPA (average shape)", 0, "in-sample", baseline_errors, sum(Y_baseline < 0)),
  cmp_row("SEPA (average shape)", 0, "leave-one-out", baseline_errors_loocv),
  cmp_row("Polynomial", p_chosen + 1, "in-sample", poly_errors, sum(Y_poly < 0)),
  cmp_row("PCA", k_chosen, "in-sample", pca_errors, sum(Y_pca < 0)),
  cmp_row("PCA", k_chosen, "leave-one-out", pca_errors_loocv),
  cmp_row("Wavelets", ncol(wav_params), "in-sample", wav_errors, sum(Y_wav < 0)),
  cmp_row("Wavelets", ncol(wav_params), "leave-one-out", wav_errors_loocv),
  cmp_row("Gamma-type", 3, "in-sample", gamma_errors, sum(Y_gamma[gamma_ok, ] < 0)),
  cmp_row("B-spline", n_basis_chosen + 1, "in-sample", bspline_errors, sum(Y_bspline < 0))
)

n_smooth_cmp <- length(wav_all[[1]]@V[[wav_levels]])   # 64 / 2^levels smooth coefficients

# Methods that learn from other events (SEPA average, PCA basis, wavelet
# coefficient choice) are rebuilt by leave-one-out. Per-event methods
# (polynomial, B-spline, gamma-type) learn nothing from other events, so their
# held-out curves are their in-sample fits.
Y_sepa_loo <- cached("cmp_sepa_loo_curves", loocv(function(i) {
  s <- colMeans(Y_unaligned[-i, ])
  baseline_curve(s / max(s), i)
}))
Y_pca_loo <- cached("cmp_pca_loo_curves",
                    loocv(function(i) pca_left_out(i, k_chosen)), key = k_chosen)
Y_wav_loo <- cached("cmp_wav_loo_curves",
                    loocv(function(i) wavelet_left_out(i, keep_chosen)),
                    key = list(keep_chosen, wav_filter, wav_levels))

all_rows <- seq_len(nrow(Y))

# One entry per method
# fit_in / fit_out : in-sample and held-out curves, laid out like Y
# rows             : events used (gamma-type: converged fits only)
# K                : number of shape parameters
# inputs           : emulator inputs = K + peak flow, duration, peak position
#                    (SEPA: peak flow and duration only)
# params           : shape parameters, one row per event (for correlations)
cmp <- list(
  "SEPA (average shape)" = list(fit_in = Y_baseline, fit_out = Y_sepa_loo, rows = all_rows,
                                K = 0, inputs = 2, params = NULL),
  "Gamma-type" = list(fit_in = Y_gamma, fit_out = Y_gamma, rows = gamma_ok,
                      K = 3, inputs = 3 + 3, params = gamma_params[gamma_ok, ]),
  "Polynomial" = list(fit_in = Y_poly, fit_out = Y_poly, rows = all_rows,
                      K = p_chosen + 1, inputs = p_chosen + 1 + 3, params = poly_coefs),
  "B-spline"   = list(fit_in = Y_bspline, fit_out = Y_bspline, rows = all_rows,
                      K = n_basis_chosen + 1, inputs = n_basis_chosen + 1 + 3,
                      params = bspline_coefs),
  "PCA"        = list(fit_in = Y_pca, fit_out = Y_pca_loo, rows = all_rows,
                      K = k_chosen, inputs = k_chosen + 3, params = pca_scores),
  "Wavelets"   = list(fit_in = Y_wav, fit_out = Y_wav_loo, rows = all_rows,
                      K = ncol(wav_params), inputs = ncol(wav_params) + 3,
                      params = wav_params)
)
cmp_order   <- names(cmp)
cmp_colours <- c("SEPA (average shape)" = "grey40", setNames(method_colours, method_names))

# Per-event errors (section 10 definitions)
cmp_err <- map(cmp, function(m) list(
  in_s = method_errors(m$fit_in)[m$rows, , drop = FALSE],
  out  = method_errors(m$fit_out)[m$rows, , drop = FALSE]
))

# Extra event-relevant features (max, width, rise, fall)
# Both curves are put on a regular 15-minute grid in real hours, then:
#   duration_error : hours above u (fitted - observed)
#   excess_error   : volume above u, million m^3 (fitted - observed)
#   rise_rmse      : RMSE on the rising limb (up to the observed peak), m^3/s
#   fall_rmse      : RMSE on the falling limb, m^3/s
flood_features <- function(hours, observed, fitted, u = u_final, dt = 0.25) {
  grid_h <- seq(0, max(hours), by = dt)
  obs    <- approx(hours, observed, xout = grid_h, rule = 2)$y
  fit    <- approx(hours, fitted,   xout = grid_h, rule = 2)$y
  rising <- grid_h <= hours[which.max(observed)]
  c(
    duration_error = dt * (sum(fit > u) - sum(obs > u)),
    excess_error   = (sum(pmax(fit - u, 0)) - sum(pmax(obs - u, 0))) * dt * 3600 / 1e6,
    rise_rmse      = sqrt(mean((fit[rising]  - obs[rising])^2)),
    fall_rmse      = sqrt(mean((fit[!rising] - obs[!rising])^2))
  )
}

feature_errors <- function(Y_fit) {
  t(sapply(seq_len(nrow(Y)), function(i) {
    observed <- as_flood(Y[i, ], i)
    fitted   <- as_flood(Y_fit[i, ], i)
    flood_features(observed$hours, observed$discharge, fitted$discharge)
  }))
}

cmp_feat <- map(cmp, function(m) feature_errors(m$fit_out)[m$rows, , drop = FALSE])

# Summary: one row per method
max_abs_cor <- function(p) {
  if (is.null(p) || ncol(p) < 2) return(NA_real_)
  r <- cor(p)
  max(abs(r[upper.tri(r)]))
}

cmp_summary <- imap_dfr(cmp, function(m, name) {
  e_in  <- cmp_err[[name]]$in_s
  e_out <- cmp_err[[name]]$out
  f     <- cmp_feat[[name]]
  tibble(
    method      = name,
    K           = m$K,
    inputs      = m$inputs,
    n_events    = length(m$rows),
    rmse_in     = mean(e_in[, "rmse"]),
    rmse_out    = mean(e_out[, "rmse"]),
    nrmse_out   = mean(100 * e_out[, "rmse"] / events$peak_flow[m$rows]),
    peak_out    = mean(abs(e_out[, "peak_error"])),
    timing_out  = mean(abs(e_out[, "timing_error"])),
    volume_out  = mean(abs(e_out[, "volume_error"])),
    peak_bias   = mean(e_out[, "peak_error"]),
    timing_bias = mean(e_out[, "timing_error"]),
    volume_bias = mean(e_out[, "volume_error"]),
    rmse_p90    = unname(quantile(e_out[, "rmse"], 0.9)),
    rmse_max    = max(e_out[, "rmse"]),
    duration    = mean(abs(f[, "duration_error"])),
    excess      = mean(abs(f[, "excess_error"])),
    excess_bias = mean(f[, "excess_error"]),
    rise        = mean(f[, "rise_rmse"]),
    fall        = mean(f[, "fall_rmse"]),
    neg_events  = sum(apply(m$fit_out[m$rows, , drop = FALSE], 1, min) < 0),
    max_cor     = max_abs_cor(m$params)
  )
})

# Worst-fitted (held-out) event for each method, for the text on unusual events
worst_events <- imap_dfr(cmp_err, function(e, name) {
  rows <- cmp[[name]]$rows
  w    <- rows[which.max(e$out[, "rmse"])]
  tibble(method = name, event = w, rmse = max(e$out[, "rmse"]),
         peak_flow = events$peak_flow[w], width_h = events$width_h[w])
})

# Accuracy-parsimony trade-off:
# Every method is evaluated at every whole number of shape parameters in K_grid,
# so the curves can be compared at the same K. Errors are averaged over ALL
# events: leave-one-out for PCA and wavelets; for the per-event methods,
# held-out = in-sample (Section 3.2).
# Structural limits: polynomial K >= 2 (degree >= 1); B-spline K >= 4
# (bs() needs >= 3 basis functions); wavelets K >= n_smooth_cmp + 1 = 9
# (all smooth coefficients are kept); gamma-type is fixed at K = 3.
K_grid <- 1:16

avg_err <- function(Y_fit) average_errors(method_errors(Y_fit))

tradeoff <- cached("cmp_tradeoff_common_K", bind_rows(
  map_dfr(K_grid[K_grid >= 2], function(K)
    data.frame(method = "Polynomial", K = K, t(avg_err(poly_fit_all(K - 1))))),
  map_dfr(K_grid[K_grid >= 4], function(K)
    data.frame(method = "B-spline", K = K, t(avg_err(bspline_fit_all(K - 1))))),
  map_dfr(K_grid, function(K)
    data.frame(method = "PCA", K = K,
               t(avg_err(loocv(function(i) pca_left_out(i, K)))))),
  map_dfr(K_grid[K_grid > n_smooth_cmp], function(K)
    data.frame(method = "Wavelets", K = K,
               t(avg_err(loocv(function(i) wavelet_left_out(i, K - n_smooth_cmp)))))),
  data.frame(method = "Gamma-type", K = 3, t(average_errors(gamma_errors)))
), key = list(K_grid, wav_filter, wav_levels)) %>%
  select(method, K, rmse, peak_error, volume_error) %>%
  arrange(method, K)

chosen_K <- c("Polynomial" = p_chosen + 1, "B-spline" = n_basis_chosen + 1,
              "PCA" = k_chosen, "Wavelets" = ncol(wav_params), "Gamma-type" = 3)
tradeoff <- tradeoff %>% mutate(chosen = K == chosen_K[method])

metric_labels <- c(rmse = "RMSE (m³/s)", peak_error = "Peak error (m³/s)",
                   volume_error = "Volume error (%)")
tradeoff_long <- tradeoff %>%
  pivot_longer(c(rmse, peak_error, volume_error), names_to = "metric") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         method = factor(method, levels = cmp_order))
sepa_ref <- tibble(metric = names(metric_labels),
                   value  = average_errors(baseline_errors_loocv)[names(metric_labels)]) %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels))

# Are the differences real? Paired comparisons over events.
# Events are independent by construction (Robson & Reed criteria), so the
# per-event differences in absolute error are resampled (bootstrap 95% CI for
# the mean difference) and tested with a Wilcoxon signed-rank test.
pairs_to_test <- list(c("PCA", "SEPA (average shape)"),
                      c("Wavelets", "SEPA (average shape)"),
                      c("Wavelets", "PCA"),
                      c("PCA", "B-spline"))
paired_metrics <- c(rmse = "RMSE (m$^3$/s)", peak_error = "Peak error (m$^3$/s)",
                    volume_error = "Volume error (\\%)")

set.seed(2026)
paired_results <- map_dfr(pairs_to_test, function(p) {
  map_dfr(names(paired_metrics), function(mt) {
    a <- abs(cmp_err[[p[1]]]$out[, mt])
    b <- abs(cmp_err[[p[2]]]$out[, mt])
    d <- a - b
    boots <- replicate(2000, mean(sample(d, replace = TRUE)))
    tibble(comparison   = paste(p[1], "minus", p[2]),
           metric       = paired_metrics[[mt]],
           diff         = mean(d),
           lower        = unname(quantile(boots, 0.025)),
           upper        = unname(quantile(boots, 0.975)),
           first_better = 100 * mean(a < b),
           p_value      = suppressWarnings(
             wilcox.test(a, b, paired = TRUE, exact = FALSE)$p.value))
  })
})

# Numbers for the text
S <- function(m, col) cmp_summary[[col]][cmp_summary$method == m]
cmp_text <- list(
  rmse_sepa      = S("SEPA (average shape)", "rmse_out"),
  rmse_pca       = S("PCA", "rmse_out"),
  rmse_wav       = S("Wavelets", "rmse_out"),
  rmse_pca_in    = S("PCA", "rmse_in"),
  red_pca        = 100 * (1 - S("PCA", "rmse_out") / S("SEPA (average shape)", "rmse_out")),
  red_wav        = 100 * (1 - S("Wavelets", "rmse_out") / S("SEPA (average shape)", "rmse_out")),
  timing_sepa    = S("SEPA (average shape)", "timing_out"),
  volume_sepa    = S("SEPA (average shape)", "volume_bias"),   # signed: + = overestimates
  excess_sepa    = S("SEPA (average shape)", "excess_bias"),
  peak_poly      = S("Polynomial", "peak_out"),
  peak_bspline   = S("B-spline", "peak_out"),
  pca_at_wav_K   = tradeoff$rmse[tradeoff$method == "PCA" & tradeoff$K == ncol(wav_params)],
  gamma_failed   = nrow(Y) - length(gamma_ok),
  wav_minus_pca  = paired_results$diff[paired_results$comparison == "Wavelets minus PCA" &
                                         paired_results$metric == paired_metrics[["rmse"]]]
)

# Formatting helper: bold the best (lowest) value in a column
fmt_best <- function(x, digits = 2) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- "--"
  best <- which(x == min(x, na.rm = TRUE))
  out[best] <- paste0("\\textbf{", out[best], "}")
  out
}

cmp_main_tbl <- cmp_summary %>%
  transmute(
    Method   = method,
    K        = as.character(K),
    Inputs   = as.character(inputs),
    rmse_in  = fmt_best(rmse_in, 2),
    rmse_out = fmt_best(rmse_out, 2),
    nrmse    = fmt_best(nrmse_out, 1),
    peak     = fmt_best(peak_out, 2),
    timing   = fmt_best(timing_out, 2),
    volume   = fmt_best(volume_out, 2),
    neg      = as.character(neg_events),
    cor      = ifelse(is.na(max_cor), "--", formatC(max_cor, format = "f", digits = 2))
  )

cmp_feat_tbl <- cmp_summary %>%
  transmute(
    Method      = method,
    peak_bias   = sprintf("%+.2f", peak_bias),
    volume_bias = sprintf("%+.1f", volume_bias),
    p90         = fmt_best(rmse_p90, 1),
    worst       = fmt_best(rmse_max, 1),
    duration    = fmt_best(duration, 1),
    excess      = fmt_best(excess, 2),
    excess_bias = sprintf("%+.2f", excess_bias),
    rise        = fmt_best(rise, 2),
    fall        = fmt_best(fall, 2)
  )

per_event_long <- imap_dfr(cmp_err, function(e, name) {
  rows <- cmp[[name]]$rows
  tibble(method = name,
         `RMSE (% of peak flow)`                = 100 * e$out[, "rmse"] / events$peak_flow[rows],
         `Absolute peak error (% of peak flow)` = 100 * abs(e$out[, "peak_error"]) / events$peak_flow[rows],
         `Absolute volume error (%)`            = abs(e$out[, "volume_error"]))
}) %>%
  pivot_longer(-method, names_to = "metric") %>%
  mutate(method = factor(method, levels = cmp_order),
         metric = factor(metric, levels = unique(metric)))

per_event_long_plot <- ggplot(per_event_long, aes(method, value, fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8, linewidth = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = cmp_colours, guide = "none") +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank())

tradeoff_long_plot <- ggplot(tradeoff_long, aes(K, value, colour = method)) +
  geom_hline(data = sepa_ref, aes(yintercept = value), inherit.aes = FALSE,
             linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = n_smooth_cmp + 1, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  geom_point(data = filter(tradeoff_long, chosen), size = 3.5, shape = 21,
             fill = NA, stroke = 1.3) +
  scale_x_continuous(limits = c(1, max(K_grid)), breaks = seq(2, max(K_grid), by = 2)) +
  scale_y_log10() +
  scale_colour_manual(values = cmp_colours, name = NULL) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Number of shape parameters K", y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

paired_tbl <- paired_results %>%
  transmute(
    Metric     = metric,
    Difference = sprintf("%+.2f", diff),
    CI         = sprintf("(%+.2f, %+.2f)", lower, upper),
    Better     = sprintf("%.0f", first_better),
    p          = ifelse(p_value < 0.001, "$<$0.001", sprintf("%.3f", p_value))
  )

paired_kable <- kable(paired_tbl, format = "latex", booktabs = TRUE, escape = FALSE,
                      align = "lrrrr",
                      col.names = c("Metric", "Mean difference", "95\\% bootstrap CI",
                                    "Events where first is better (\\%)", "Wilcoxon $p$"),
                      caption = "Paired comparison of held-out absolute errors over the same extreme flow events. A negative difference means the first method has smaller errors.") %>%
  kable_styling(latex_options = c("hold_position"), font_size = 9)

groups <- unique(paired_results$comparison)
for (g in seq_along(groups)) {
  rows_g <- which(paired_results$comparison == groups[g])
  paired_kable <- pack_rows(paired_kable, groups[g], min(rows_g), max(rows_g))
}

f1 <- function(x, d = 1) formatC(x, format = "f", digits = d)

synthesis <- tribble(
  ~Method, ~Accuracy, ~Features, ~Unusual, ~Parsimony, ~Behaviour, ~Valid, ~Interpret,
  
  "SEPA (average shape)",
  sprintf("Poor: %s m$^3$/s", f1(S("SEPA (average shape)", "rmse_out"))),
  sprintf("Peak matched; timing %s h; volume %+.0f\\%%",
          f1(S("SEPA (average shape)", "timing_out"), 0), S("SEPA (average shape)", "volume_bias")),
  "One shape for every event",
  "0 shape parameters",
  "None to model",
  "Valid",
  "High",
  
  "Gamma-type",
  sprintf("Moderate: %s m$^3$/s", f1(S("Gamma-type", "rmse_out"))),
  sprintf("Peak %s m$^3$/s; volume %s\\%%", f1(S("Gamma-type", "peak_out"), 2),
          f1(S("Gamma-type", "volume_out"))),
  "Single peak only; same baseflow at start and end",
  "3",
  sprintf("Weakly correlated (max %s)", f1(S("Gamma-type", "max_cor"), 2)),
  sprintf("Valid; %d fits failed", cmp_text$gamma_failed),
  "High (baseflow, peak time, sharpness)",
  
  "Polynomial",
  sprintf("Moderate: %s m$^3$/s", f1(S("Polynomial", "rmse_out"))),
  sprintf("Peak too low by %s m$^3$/s", f1(S("Polynomial", "peak_out"))),
  "Fails multi-peak events; drifts at ends",
  as.character(S("Polynomial", "K")),
  sprintf("Near-collinear (max %s)", f1(S("Polynomial", "max_cor"), 2)),
  sprintf("%d events negative", S("Polynomial", "neg_events")),
  "Low",
  
  "B-spline",
  sprintf("Moderate: %s m$^3$/s", f1(S("B-spline", "rmse_out"))),
  sprintf("Peak too low by %s m$^3$/s", f1(S("B-spline", "peak_out"))),
  "Good tails; peaks flattened",
  as.character(S("B-spline", "K")),
  sprintf("Correlated (max %s)", f1(S("B-spline", "max_cor"), 2)),
  sprintf("%d events negative", S("B-spline", "neg_events")),
  "Low",
  
  "PCA",
  sprintf("Very good: %s m$^3$/s", f1(S("PCA", "rmse_out"), 2)),
  sprintf("Peak %s m$^3$/s; volume %s\\%%", f1(S("PCA", "peak_out"), 2),
          f1(S("PCA", "volume_out"), 2)),
  "Close on multi-peak; smaller peaks underestimated",
  as.character(S("PCA", "K")),
  "Uncorrelated by construction",
  sprintf("%d events negative", S("PCA", "neg_events")),
  "Moderate (breadth, rise vs recession)",
  
  "Wavelets",
  sprintf("Very good: %s m$^3$/s", f1(S("Wavelets", "rmse_out"), 2)),
  sprintf("Peak %s m$^3$/s; volume %s\\%%", f1(S("Wavelets", "peak_out"), 2),
          f1(S("Wavelets", "volume_out"), 2)),
  "Close, but jagged dips and end artefacts",
  as.character(S("Wavelets", "K")),
  sprintf("Smooth coefficients correlated (max %s)", f1(S("Wavelets", "max_cor"), 2)),
  sprintf("%d events negative", S("Wavelets", "neg_events")),
  "Low"
)
# Qualitative wording (Unusual, Interpret) summarises Sections 4.1-4.3:
# re-check it if the results change.

# ------------------------------------------------------------------------------
# 19. Uncertainty in the shape parameters
# ------------------------------------------------------------------------------
# Uncertainty analysis of the two chosen methods: PCA (k_chosen components)
# and wavelets (la8, 3 levels, 8 smooth + 5 detail = 13 shape parameters).
#
#
# The shape parameters themselves come from earlier sections:
#   PCA      : pca_scores  (section 13, one row per event, k_chosen columns)
#   Wavelets : wav_params  (section 14, one row per event, 13 columns)
# To simulate new events we need a probability distribution for these
# parameters.
#
#   19a. Curves from simulated parameters
#   19b. Looking at the shape parameters: distributions, dependence on the
#        size and timing of the event, and residuals (Section 6.1)
#   19c. Conditional models: a statistical model and a sampling scheme
#        (Section 6.2)
#   19d. Out-of-sample check (Section 6.3)
#   19e. Final simulations (Section 6.4)
#   19f. Sampling uncertainty (Section 6.5)

# Shape parameters from sections 13 and 14, labelled so plots and tables
# are readable:
#   PCA      : scores on the kept components (PC1, PC2, ...)
#   Wavelets : smooth coefficients (S1, S2, ...) followed by the kept
#              detail coefficients (D1, D2, ...). With la8 and 3 levels
#              there are 64 / 2^3 = 8 smooth coefficients.
# Both are stored in one list so steps can be run for both methods at once.
colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))
n_smooth <- length(wav_all[[1]]@V[[wav_levels]])
colnames(wav_params) <- c(paste0("S", 1:n_smooth),
                          paste0("D", 1:(ncol(wav_params) - n_smooth)))

params <- list(PCA = pca_scores, Wavelets = wav_params)


# ------------------------------------------------------------------------------
# 19a. Curves from simulated parameters
# ------------------------------------------------------------------------------
# Sections 13 and 14 rebuild the observed events only. These two functions
# turn any shape parameters (one row per event), e.g. simulated ones from
# 19c, into curves laid out like Y. By default they use the components / kept
# coefficients found from all events; 19d passes in ones found from training
# events only.

# PCA: mean curve + scores x loadings
pca_sim_curves <- function(scores, pca_fit = pca_result) {
  loadings <- pca_fit$rotation[, 1:ncol(scores)]
  sweep(scores %*% t(loadings), 2, pca_fit$center, FUN = "+")
}

# Wavelets: put the smooth and kept detail coefficients back into a wavelet
# object and invert the transform
wav_sim_curves <- function(params, kept = kept_coefs) {
  t(apply(params, 1, function(p) {
    wt <- wav_all[[1]]                                   # used as a template
    wt@V[[wav_levels]][, 1] <- p[1:n_smooth]
    detail       <- numeric(length(kept))
    detail[kept] <- p[-(1:n_smooth)]
    level <- rep(seq_along(wt@W), sapply(wt@W, nrow))   # level of each detail coefficient
    for (j in seq_along(wt@W)) wt@W[[j]][, 1] <- detail[level == j]
    as.numeric(idwt(wt))
  }))
}

# observed parameters should give back the curves from sections 13 and 14 (differences of 0)
cat("PCA check:",      max(abs(pca_sim_curves(pca_scores) - Y_pca)), "\n")
cat("Wavelets check:", max(abs(wav_sim_curves(wav_params) - Y_wav)), "\n")


# ------------------------------------------------------------------------------
# 19b. Looking at the shape parameters (Section 6.1)
# ------------------------------------------------------------------------------
# Normal Q-Q plot of each shape parameter: if a parameter is roughly Normal,
# the points lie close to the black line; curving away at the ends means skew
# or long tails. Mostly Normal with some heavy tails motivates fitting a
# Normal model and also a bootstrap, which does not assume Normality.
# Row 1: PCA scores (blue). Rows 2-3: wavelet parameters (orange), smooth
# coefficients S1-S8 then detail coefficients D1-D5.

pca_colour <- "#0072B2"
wav_colour <- "#E69F00"

# Normal Q-Q plot of one parameter
qq_panel <- function(values, name, colour) {
  qqnorm(values, pch = 19, cex = 0.7, col = colour, main = name, xlab = "", ylab = "")
  qqline(values, col = "black", lwd = 1.5)
}

plot_all_qq <- function() {
  par(mfrow = c(3, 7), mar = c(2.5, 2.5, 2, 0.5), oma = c(2.5, 2.5, 0, 0),
      mgp = c(2, 0.6, 0), cex.axis = 0.8)
  
  # row 1: PCA scores, then the legend
  for (j in 1:ncol(pca_scores)) {
    qq_panel(pca_scores[, j], colnames(pca_scores)[j], pca_colour)
  }
  
  plot.new()
  legend("center", c("PCA scores", "Wavelet parameters"), col = c(pca_colour, wav_colour),
         pch = 19, bty = "n")
  
  # rows 2-3: wavelet parameters
  for (j in 1:ncol(wav_params)) {
    qq_panel(wav_params[, j], colnames(wav_params)[j], wav_colour)
  }
  
  mtext("Theoretical Normal quantiles", side = 1, outer = TRUE, line = 1)
  mtext("Sample quantiles",             side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0),
      mgp = c(3, 1, 0), cex.axis = 1)
}

# Is shape related to the size and timing of the event? If so, shape
# should be simulated conditional on peak flow, duration and peak position
# (19c), rather than separately from them.
flood_size_timing <- cbind(peak_flow     = events$peak_flow,
                           duration      = events$width_h,
                           peak_position = events$peak_position)

pca_dependence <- round(cor(pca_scores, flood_size_timing, method = "spearman"), 2)
wav_dependence <- round(cor(wav_params, flood_size_timing, method = "spearman"), 2)

# Follow-up check: once size and timing are accounted for, is the remaining
# variation still roughly Normal? Normal Q-Q plot of the residuals from
# regressing each shape parameter on peak flow, duration and peak position
# (the regression used in 19c). PCA top row, wavelets rows 2-3.
pca_residuals <- resid(lm(pca_scores ~ flood_size_timing))
wav_residuals <- resid(lm(wav_params ~ flood_size_timing))

plot_residual_qq <- function() {
  par(mfrow = c(3, 7), mar = c(2.5, 2.5, 2, 0.5), oma = c(2.5, 2.5, 0, 0),
      mgp = c(2, 0.6, 0), cex.axis = 0.8)
  for (j in 1:ncol(pca_residuals)) {
    qq_panel(pca_residuals[, j], colnames(pca_residuals)[j], pca_colour)
  }
  plot.new()
  legend("center", c("PCA residuals", "Wavelet residuals"), col = c(pca_colour, wav_colour),
         pch = 19, bty = "n")
  for (j in 1:ncol(wav_residuals)) {
    qq_panel(wav_residuals[, j], colnames(wav_residuals)[j], wav_colour)
  }
  mtext("Theoretical Normal quantiles", side = 1, outer = TRUE, line = 1)
  mtext("Sample quantiles",             side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0),
      mgp = c(3, 1, 0), cex.axis = 1)
}

# ------------------------------------------------------------------------------
# 19c. Conditional models for the shape parameters (Section 6.2)
# ------------------------------------------------------------------------------
# Shape depends on peak flow, duration and peak position (19b), so new shape
# parameters are simulated given these three, as the sum of
#   expected : regression of the shape parameters on peak flow, duration and
#              peak position, predicted for the new event
#   variation: what the regression cannot explain, either
#     - conditional Normal (statistical model): drawn from a multivariate
#       Normal with the covariance of the regression residuals
#     - conditional bootstrap (sampling scheme): a real event's residual plus
#       a little Normal noise (size h from Silverman's rule of thumb), divided
#       by sqrt(1 + h^2) so the noise does not make the spread too wide
# Both take the same arguments:
#   shape_params    : observed shape parameters, one row per event
#   size_timing     : observed peak flow, duration and peak position
#   new_size_timing : size and timing of the events to simulate, one row each

silverman_h <- function(residuals) {
  (4 / ((ncol(residuals) + 2) * nrow(residuals)))^(1 / (ncol(residuals) + 4))
}

simulate_conditional_normal <- function(shape_params, size_timing, new_size_timing) {
  regression <- lm(shape_params ~ size_timing)
  expected   <- cbind(1, new_size_timing) %*% coef(regression)
  variation  <- MASS::mvrnorm(nrow(new_size_timing), rep(0, ncol(shape_params)),
                              cov(resid(regression)))
  expected + variation
}

simulate_conditional_bootstrap <- function(shape_params, size_timing, new_size_timing) {
  regression <- lm(shape_params ~ size_timing)
  residuals  <- resid(regression)
  expected   <- cbind(1, new_size_timing) %*% coef(regression)
  n_new      <- nrow(new_size_timing)
  h          <- silverman_h(residuals)
  variation  <- (residuals[sample(nrow(residuals), n_new, replace = TRUE), ] +
                   MASS::mvrnorm(n_new, rep(0, ncol(residuals)), h^2 * cov(residuals))) /
    sqrt(1 + h^2)
  expected + variation
}

simulators <- list("Conditional Normal"    = simulate_conditional_normal,
                   "Conditional bootstrap" = simulate_conditional_bootstrap)

# The four combinations, in panel order:
# PCA Normal, PCA bootstrap, wavelets Normal, wavelets bootstrap
combos <- expand.grid(scheme = names(simulators), method = names(params),
                      stringsAsFactors = FALSE)
combo_names   <- paste0(combos$method, ", ", combos$scheme)
method_colour <- c(PCA = pca_colour, Wavelets = wav_colour)
sim_to_curves <- list(PCA = pca_sim_curves, Wavelets = wav_sim_curves)


# ------------------------------------------------------------------------------
# 19d. Out-of-sample check (Section 6.3)
# ------------------------------------------------------------------------------
# Lecture 5: check whether the fitted distribution works well out of sample.
# Each event in turn is left out; the method (PCA components / wavelet
# coefficients kept) and the distribution are fitted to the other events, and
# n_each curves are simulated with the left-out event's own peak flow,
# duration and peak position, then compared with its observed curve in m^3/s:
#   coverage   : % of the observed curve inside the simulated 5-95% band
#                (should be close to 90)
#   crps       : average distance between the simulated and observed curve,
#                minus half the average distance between simulated curves
#                (m^3/s, lower is better; for a single curve, as for SEPA, it
#                is the mean absolute error)
#   negative   : % of simulated curves with negative flow
#   above_peak : % of simulated curves rising more than 1% above the peak

n_each <- 200

# Wavelet coefficients to keep, chosen from the training events only, and the
# shape parameters of a set of events using them
wav_kept_from <- function(rows) {
  sizes <- unlist(average_detail_size(wav_all[rows]))
  sizes >= sort(sizes, decreasing = TRUE)[keep_chosen]
}
wav_params_from <- function(rows, kept) {
  t(sapply(wav_all[rows], function(wt) c(wt@V[[wav_levels]][, 1], unlist(wt@W)[kept])))
}

# Simulated standardised curves for event i, fitted without event i
held_out_curves <- function(i, method, simulate) {
  train <- setdiff(1:nrow(Y), i)
  new_size_timing <- flood_size_timing[rep(i, n_each), ]
  if (method == "PCA") {
    pca_fit <- prcomp(Y[train, ], center = TRUE, scale. = FALSE)
    shapes  <- simulate(pca_fit$x[, 1:k_chosen], flood_size_timing[train, ], new_size_timing)
    pca_sim_curves(shapes, pca_fit)
  } else {
    kept   <- wav_kept_from(train)
    shapes <- simulate(wav_params_from(train, kept), flood_size_timing[train, ], new_size_timing)
    wav_sim_curves(shapes, kept)
  }
}

# CRPS of one observed value y against simulated values x
crps <- function(x, y) {
  mean(abs(x - y)) - 0.5 * mean(abs(outer(x, x, "-")))
}

# Scores for event i from its simulated standardised curves
score_event <- function(curves, i) {
  simulated <- curves * events$peak_flow[i]                     # m^3/s
  observed  <- Y[i, ] * events$peak_flow[i]
  lower <- apply(simulated, 2, quantile, 0.05)
  upper <- apply(simulated, 2, quantile, 0.95)
  c(coverage   = 100 * mean(observed >= lower & observed <= upper),
    crps       = mean(sapply(1:n_grid, function(j) crps(simulated[, j], observed[j]))),
    negative   = 100 * mean(apply(curves, 1, min) < 0),
    above_peak = 100 * mean(apply(curves, 1, max) > 1.01))
}

# SEPA's approach: one average shape of the other events, no uncertainty
score_sepa <- function(i) {
  average <- colMeans(Y_unaligned[-i, ])
  shape   <- baseline_curve(average / max(average), i)
  c(coverage = NA, crps = events$peak_flow[i] * mean(abs(shape - Y[i, ])),
    negative = NA, above_peak = NA)
}

# Average scores over all left-out events (takes a few minutes; cached)
out_of_sample <- cached("out_of_sample_check", {
  set.seed(2026)
  rbind(
    data.frame(method = "SEPA", scheme = "Average shape",
               t(colMeans(t(sapply(1:nrow(Y), score_sepa))))),
    do.call(rbind, lapply(names(params), function(method) {
      do.call(rbind, lapply(names(simulators), function(scheme) {
        scores <- t(sapply(1:nrow(Y), function(i)
          score_event(held_out_curves(i, method, simulators[[scheme]]), i)))
        data.frame(method = method, scheme = scheme, t(colMeans(scores)))
      }))
    }))
  )
}, key = list(n_each, k_chosen, wav_filter, wav_levels, keep_chosen, names(simulators)))

print(out_of_sample, digits = 3)

# Value from the out-of-sample table, for inline numbers in the report
oos <- function(method, scheme, column) {
  out_of_sample[[column]][out_of_sample$method == method & out_of_sample$scheme == scheme]
}

options(knitr.kable.NA = "--")
out_of_sample_kable <- kable(
  out_of_sample[, -1], format = "latex", booktabs = TRUE, digits = 2, row.names = FALSE,
  col.names = c("", "Coverage (\\%)", "CRPS (m$^3$/s)", "Negative flow (\\%)", "Above peak (\\%)"),
  escape = FALSE,
  caption = "Out-of-sample check: each event is left out of the fit, and 200 curves are simulated with its own peak flow, duration and peak position. Coverage: percentage of the observed curve inside the simulated 90\\% band (target 90\\%). CRPS: lower is better; for SEPA's single shape it is the mean absolute error. Last two columns: percentage of simulated curves that are physically invalid. Coverage, negative flow and above peak are measured on simulated curves, so do not apply to SEPA's single fixed shape.") %>%
  kable_styling(latex_options = "hold_position", font_size = 9) %>%
  pack_rows("SEPA", 1, 1) %>%
  pack_rows("PCA", 2, 3) %>%
  pack_rows("Wavelets", 4, 5)

# Figure: simulated 90% band against the observed curve, each event left out
# of the fit. One row per combination, one column per event. By default the
# latest-peak, longest and flattest of the Section 4.1 events.
plot_held_out_bands <- function(events_to_plot = example_events[c(2, 4, 5)]) {
  set.seed(2026)
  par(mfrow = c(4, length(events_to_plot)), mar = c(2, 2, 2, 0.5), oma = c(3, 3, 0, 0))
  for (r in 1:4) {
    colour <- method_colour[[combos$method[r]]]
    for (i in events_to_plot) {
      simulated <- held_out_curves(i, combos$method[r], simulators[[combos$scheme[r]]]) *
        events$peak_flow[i]
      observed  <- as_flood(Y[i, ], i)
      lower <- apply(simulated, 2, quantile, 0.05)
      upper <- apply(simulated, 2, quantile, 0.95)
      plot(observed$hours, observed$discharge, type = "n",
           ylim = range(lower, upper, observed$discharge),
           xlab = "", ylab = "", main = paste0(combo_names[r], " - event ", i), cex.main = 0.85)
      polygon(c(observed$hours, rev(observed$hours)), c(lower, rev(upper)),
              col = adjustcolor(colour, 0.3), border = NA)
      lines(observed$hours, observed$discharge, lwd = 2)
    }
  }
  mtext("Hours since start", side = 1, outer = TRUE, line = 1)
  mtext("Discharge (m^3/s)", side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
}


# ------------------------------------------------------------------------------
# 19e. Final simulations (Section 6.4)
# ------------------------------------------------------------------------------
# All four combinations, fitted to all events. Figures have one panel per
# combination.

# Simulated standardised curves for one combination, given size and timing
simulate_curves <- function(method, scheme, new_size_timing) {
  shapes <- simulators[[scheme]](params[[method]], flood_size_timing, new_size_timing)
  sim_to_curves[[method]](shapes)
}

# 1000 simulated events: size and timing taken from real events at random,
# then shape simulated given them
n_sim <- 1000
set.seed(2026)
sim_size_timing <- flood_size_timing[sample(nrow(Y), n_sim, replace = TRUE), ]
final_shapes <- lapply(1:4, function(r)
  simulate_curves(combos$method[r], combos$scheme[r], sim_size_timing))


# Fan charts: 5-95% and 25-75% bands and median of the simulated shapes, with
# the real events in grey
plot_fan_charts <- function() {
  bands <- lapply(final_shapes, function(s)
    apply(s, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
  y_range <- range(Y, unlist(bands))                     # same y-axis on every panel
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  for (r in 1:4) {
    colour <- method_colour[[combos$method[r]]]
    band   <- bands[[r]]
    plot(grid_u, band[3, ], type = "n", ylim = y_range, main = combo_names[r],
         xlab = "Rescaled time u", ylab = "Flow / peak flow")
    polygon(c(grid_u, rev(grid_u)), c(band[1, ], rev(band[5, ])),
            col = adjustcolor(colour, 0.2), border = NA)
    polygon(c(grid_u, rev(grid_u)), c(band[2, ], rev(band[4, ])),
            col = adjustcolor(colour, 0.4), border = NA)
    matlines(grid_u, t(Y), col = adjustcolor("grey30", 0.3), lty = 1)
    lines(grid_u, band[3, ], col = colour, lwd = 2)
  }
  legend("topright", c("5-95%", "25-75%", "Median", "Real events"), bty = "n", cex = 0.8,
         fill = c(adjustcolor("grey50", 0.2), adjustcolor("grey50", 0.4), NA, NA),
         border = NA, col = c(NA, NA, "grey20", "grey30"), lty = c(NA, NA, 1, 1))
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}


# What shape uncertainty means for flooding: for an event of typical peak
# flow, duration and peak position (the medians), the range of flood volumes
# from 1000 simulated shapes, against the single volume from SEPA's shape
typical <- apply(flood_size_timing, 2, median)
typical_size_timing <- matrix(typical, n_sim, 3, byrow = TRUE)

volume_of <- function(shape) {
  flood <- reconstruct_hydrograph(shape, typical["peak_position"], typical["duration"],
                                  typical["peak_flow"])
  area_under(flood$hours, flood$discharge) * 3600 / 1e6        # million m^3
}

sepa_shape  <- colMeans(Y_unaligned) / max(colMeans(Y_unaligned))
sepa_volume <- area_under(grid_u * typical["duration"], sepa_shape * typical["peak_flow"]) *
  3600 / 1e6

set.seed(2026)
volume_table <- rbind(
  data.frame(scheme = "SEPA (single shape)", lower = NA, median = unname(sepa_volume), upper = NA),
  do.call(rbind, lapply(1:4, function(r) {
    volumes <- apply(simulate_curves(combos$method[r], combos$scheme[r], typical_size_timing),
                     1, volume_of)
    data.frame(scheme = combo_names[r],
               lower  = unname(quantile(volumes, 0.05)),
               median = median(volumes),
               upper  = unname(quantile(volumes, 0.95)))
  }))
)
print(volume_table, digits = 3)

options(knitr.kable.NA = "--")
volume_kable <- kable(
  volume_table, format = "latex", booktabs = TRUE, digits = 2, row.names = FALSE,
  col.names = c("", "5\\%", "Median", "95\\%"), escape = FALSE,
  caption = paste0("Flood volume (million m$^3$) of an event with the median peak flow (",
                   round(typical["peak_flow"]), " m$^3$/s), duration and peak position: ",
                   "SEPA's single shape against 1000 simulated shapes from each scheme.")) %>%
  kable_styling(latex_options = "hold_position", font_size = 9)


# Simulating whole events in m^3/s and hours, as they would be fed to the
# flood inundation model or emulator: draw a peak flow, duration and peak
# position, then a shape given them, divided by its own maximum so the event
# peaks at exactly the drawn peak flow.
# PLACEHOLDER: draw_size_timing() resamples the real events; replace with the
# Part 1 distributions (e.g. the GPD for peak flow) if available.
draw_size_timing <- function(n) flood_size_timing[sample(nrow(Y), n, replace = TRUE), ]

simulate_events <- function(n, method, scheme) {
  size_timing <- draw_size_timing(n)
  shapes <- simulate_curves(method, scheme, size_timing)
  shapes <- shapes / apply(shapes, 1, max)
  lapply(1:n, function(j) reconstruct_hydrograph(shapes[j, ], size_timing[j, "peak_position"],
                                                 size_timing[j, "duration"],
                                                 size_timing[j, "peak_flow"]))
}

# 30 simulated events per combination (coloured) over the real events (grey)
plot_simulated_events <- function(n = 30) {
  real <- lapply(1:nrow(Y), function(i) as_flood(Y[i, ], i))
  set.seed(2026)
  simulated  <- lapply(1:4, function(r) simulate_events(n, combos$method[r], combos$scheme[r]))
  all_events <- c(real, unlist(simulated, recursive = FALSE))
  max_hours  <- max(sapply(all_events, function(f) max(f$hours)))
  flow_range <- range(sapply(all_events, function(f) range(f$discharge)))
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  for (r in 1:4) {
    plot(NA, xlim = c(0, max_hours), ylim = flow_range, main = combo_names[r],
         xlab = "Hours since start", ylab = "Discharge (m^3/s)")
    for (f in real) lines(f$hours, f$discharge, col = adjustcolor("grey30", 0.3))
    for (f in simulated[[r]]) lines(f$hours, f$discharge,
                                    col = adjustcolor(method_colour[[combos$method[r]]], 0.7))
  }
  legend("topright", c("Real events", "Simulated events"), lty = 1, bty = "n", cex = 0.8,
         col = c("grey30", "grey20"))
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}


# ------------------------------------------------------------------------------
# 19f. Sampling uncertainty (Section 6.5)
# ------------------------------------------------------------------------------
# Lecture 5: "your dataset is limited in size". For each combination, resample
# the events with replacement, refit the method and distribution, and record
# the simulated 90% band; if the bands from different refits are close, the
# results do not depend much on which events we happened to observe.
n_boot <- 100

# Simulated 5% and 95% curves from one refit on the events in rows
refitted_band <- function(rows, method, simulate) {
  new_size_timing <- flood_size_timing[sample(rows, n_sim, replace = TRUE), ]
  if (method == "PCA") {
    pca_fit <- prcomp(Y[rows, ], center = TRUE, scale. = FALSE)
    curves  <- pca_sim_curves(simulate(pca_fit$x[, 1:k_chosen], flood_size_timing[rows, ],
                                       new_size_timing), pca_fit)
  } else {
    kept   <- wav_kept_from(rows)
    curves <- wav_sim_curves(simulate(wav_params_from(rows, kept), flood_size_timing[rows, ],
                                      new_size_timing), kept)
  }
  apply(curves, 2, quantile, probs = c(0.05, 0.95))
}

# refitted_bands[[r]][[b]]: band from refit b of combination r (slow; cached)
refitted_bands <- cached("sampling_check_all", {
  set.seed(2026)
  lapply(1:4, function(r) {
    lapply(1:n_boot, function(b) {
      rows <- sample(nrow(Y), replace = TRUE)
      refitted_band(rows, combos$method[r], simulators[[combos$scheme[r]]])
    })
  })
}, key = list(n_boot, n_sim, k_chosen, wav_filter, wav_levels, keep_chosen, names(simulators)))

# Edges of the 90% band from every refit: lower and upper edges in the
# combination's colour
plot_sampling_bands <- function() {
  y_range <- range(unlist(refitted_bands))
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  for (r in 1:4) {
    colour <- adjustcolor(method_colour[[combos$method[r]]], 0.2)
    matplot(grid_u, sapply(refitted_bands[[r]], function(b) b[1, ]), type = "l", lty = 1,
            col = colour, ylim = y_range, main = combo_names[r],
            xlab = "Rescaled time u", ylab = "Flow / peak flow")
    matlines(grid_u, sapply(refitted_bands[[r]], function(b) b[2, ]), lty = 1, col = colour)
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

# Average width of the 90% band, and how much its edges move between refits
# (fractions of the peak flow)
sampling_table <- data.frame(
  scheme     = combo_names,
  band_width = sapply(refitted_bands, function(bands)
    mean(sapply(bands, function(b) mean(b[2, ] - b[1, ])))),
  lower_sd   = sapply(refitted_bands, function(bands)
    mean(apply(sapply(bands, function(b) b[1, ]), 1, sd))),
  upper_sd   = sapply(refitted_bands, function(bands)
    mean(apply(sapply(bands, function(b) b[2, ]), 1, sd)))
)
print(sampling_table, digits = 3)

sampling_kable <- kable(
  sampling_table, format = "latex", booktabs = TRUE, digits = 3, row.names = FALSE,
  col.names = c("", "Band width", "SD of lower edge", "SD of upper edge"),
  caption = paste0("Sampling uncertainty: average width of the simulated 90\\% band, and the ",
                   "average standard deviation of its edges across ", n_boot,
                   " refits to resampled events (fractions of the peak flow)."),
  escape = FALSE) %>%
  kable_styling(latex_options = "hold_position", font_size = 9)