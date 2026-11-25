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
library(patchwork)
library(wavelets)
library(splines)
library(kableExtra)
library(gridGraphics)

# Longer timeout for the SEPA download
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

# Codes kept for modelling. Suspect is checked separately in section 5
retained_codes <- c(50, 100, 200, 254)   # Good, Estimated, Unchecked (legacy), Unchecked
suspect_code   <- 150                    # Suspect: beyond the acceptable extrapolation limit

# Full record before filtering, used to check Suspect readings in section 5
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

# Vectors used for finding peaks
q <- river_data$Discharge
t <- river_data$Timestamp
n_years <- n_distinct(river_data$water_year)

# Subsection 2.1
#-------------------------------------------------------------------------------

# Data Description

# Missing readings, compared with a full 15-minute sequence
complete_sequence <- seq(min(river_data_raw$Timestamp), 
                         max(river_data$Timestamp), by = "15 min")
number_missing <- sum(!complete_sequence %in% river_data_raw$Timestamp)

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
    paste(format(min(river_data_raw$Timestamp), "%d %B %Y"), "-",
          format(max(river_data_raw$Timestamp), "%d %B %Y")),
    format(nrow(river_data_raw), big.mark = ","),
    format(number_missing, big.mark = ",")
  )
)

# Data Quality

# Quality codes from SEPA guidance, and which ones we use
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
    c("retained", "retained", "N/A", "excluded", "retained", "retained")
)

# Add the number of readings for each code
quality_combined <- quality_code_table %>%
  left_join(
    river_data_raw %>% count(QualityCode = as.numeric(QualityCode), name = "Count"),
    by = c("Code number" = "QualityCode")
  ) %>%
  mutate(
    Count = coalesce(Count, 0L),
    Percentage = sprintf("%.2f\\%%", Count / sum(Count) * 100)
  )

# One table: data description, then one row per quality code
data_table <- bind_rows(
  data_description_table,
  quality_combined %>%
    transmute(
      Item  = paste0(
        `Code number`, " - ", `Code name`, " (", `Abbreviated code name`, ")"),
      Value = paste0(format(Count, big.mark = ",", trim = TRUE), " readings (",
                     Percentage, "), ", Used)
    )
)

# Combined table
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
# 3. Peak Identification (Robson & Reed 1999 criteria, u = 65)
# ------------------------------------------------------------------------------

# Start of the rising limb (same as Part 1)
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

# Apply the independence criteria
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

# Independent peaks at u = 65
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
  
  # Resample onto the 481-point time grid
  approx(x = event_t_rel, y = df_win$Discharge, xout = rel_time_hours, rule = 2)$y
}

N_events  <- nrow(peaks)
peak_rows <- match(peaks$Timestamp, t)

# Y_raw: one row per event, 481 columns
Y_raw <- t(sapply(1:N_events, function(i) extract_event_hydrograph(i, peak_rows[i])))
rownames(Y_raw) <- paste0("Event_", 1:N_events)
colnames(Y_raw) <- paste0("t_", sprintf("%+.2fh", rel_time_hours))

# Mean hydrograph
mean_hydrograph <- colMeans(Y_raw)

# ------------------------------------------------------------------------------
# 5. Unchecked readings vs interpolation for peak windows
# ------------------------------------------------------------------------------

# Each unchecked burst vs linear interpolation between the readings either side.
# Uses the full record so bursts at a window edge still have anchors.
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
    
    # No anchor at the start or end of the record, so skip
    if (idx_before < 1 || idx_after > nrow(d)) next
    
    interp <- approx(
      x    = tnum[c(idx_before, idx_after)],
      y    = d$Discharge[c(idx_before, idx_after)],
      xout = tnum[bi]
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

# Flagged readings inside each event window. Defaults to Unchecked readings.
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

# Stats, table and scatter plot for flagged readings. NULL if none.
summarise_unchecked_readings <- function(in_windows, n_events, label = "Unchecked") {
  if (nrow(in_windows) == 0) return(NULL)
  
  # count each reading once if windows overlap
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
    Start        = "-",
    Readings     = sprintf("**%d**", nrow(unique_readings)),
    Duration_Min = "-",
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
      subtitle = bquote("n =" ~ .(nrow(unique_readings)) ~ "| RMSE =" ~
                          .(sprintf("%.3f", rmse)) ~ m^3/s ~ "| dashed line = 1:1"),
      x = expression("Interpolated discharge (" * m^3 * "/s)"),
      y = bquote(.(label) ~ "discharge (" * m^3 * "/s)")
    ) +
    theme_minimal()
  
  list(in_windows = in_windows, unique_readings = unique_readings,
       rmse = rmse, max_ad = max_ad, max_q = max_q,
       burst_tbl = burst_tbl, overall_row = overall_row,
       table = bind_rows(burst_tbl, overall_row), plot = plot)
}

# Hydrographs of events with flagged readings, one panel per event
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
      subtitle = bquote("Line: readings used in the model; points:" ~ .(tolower(label)) ~
                          "readings; dotted: u =" ~ .(u_final) ~ m^3/s),
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
  # Names used in the Rmd
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

# Suspect readings are excluded, so the model interpolates across them.
# Compare each one with that interpolation.

# Record-wide summary for the report text
suspect_readings <- river_data_raw %>% filter(QualityCode == suspect_code)
n_suspect        <- nrow(suspect_readings)
pct_suspect      <- 100 * n_suspect / nrow(river_data_raw)
n_suspect_above  <- sum(suspect_readings$Discharge > u_final)
max_suspect      <- if (n_suspect > 0) max(suspect_readings$Discharge) else NA_real_

suspect_all <- compare_unchecked(
  river_data_raw %>% filter(QualityCode %in% c(retained_codes, suspect_code)),
  unchecked_codes = suspect_code
)
# Bursts at either end of the record are skipped
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

# All event hydrographs centred on the peak, mean in bold
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
  annotate("text", x = -45, y = u_final + 5, label = paste0("u == ", u_final, " ~ m^3/s"), parse = TRUE,
           colour = "black", hjust = 0, size = 3.2) +
  scale_x_continuous(breaks = seq(-48, 72, by = 12)) +
  labs(
    title = paste0("River Don Hydrograph Shape Variability (N = ", N_events, " Events)"),
    subtitle = expression("Centred at peak arrival (t = 0) [-48h, +72h]. Bold line: mean curve" ~
                            bar(Q)(t) * ". Unchecked values retained; Suspect (code 150) excluded."),
    x = "Hours from Peak Arrival (t = 0)",
    y = expression("Discharge (" * m^3 * "/s)")
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, colour = "grey30")
  )

print(hydrograph_overlay_plot)

# ------------------------------------------------------------------------------
# 7. Event start, end and width (trough to trough)
# ------------------------------------------------------------------------------

search_days   <- 7     # how far either side of a peak to search for a trough
trough_factor <- 1.05  # trough = reading within 5% of the local minimum

# End of the falling limb. Mirrors rise_start_row.
fall_end_row <- function(p, limit) {
  window <- which(t >= t[p] & t <= limit)
  lowest <- min(q[window])
  min(window[q[window] <= trough_factor * lowest])
}

# Search up to 7 days after the peak, or to the next peak if sooner
end_limit <- pmin(
  t[peak_rows] + search_days * 24 * 3600,
  lead(t[peak_rows], default = max(t))
)

# One row per event: start, peak, end and timings
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
  # $table and $plot are used in the report
} else {
  cat("No unchecked quality codes within the modelling events.\n\n")
}

# Suspect readings within the trough-to-trough events
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
# Rescale time to [0, 1], move the peak to u_star, divide by peak flow,
# then read off at 64 points. One row of Y per event.

u_star <- mean(events$peak_position)
n_grid <- 64
grid_u <- seq(0, 1, length.out = n_grid)

standardise_event <- function(i) {
  rows <- events$start_row[i]:events$end_row[i]
  p    <- events$peak_position[i]

  # rescale time onto [0, 1]
  u <- as.numeric(difftime(t[rows], events$start_time[i], units = "hours")) /
    events$width_h[i]
  
  # realign so the rise is on [0, u_star] and the fall on [u_star, 1]
  u_reg <- ifelse(u <= p,
                  u * u_star / p,
                  u_star + (u - p) * (1 - u_star) / (1 - p))
  
  # divide by peak flow and read off at the grid points
  approx(u_reg, q[rows] / events$peak_flow[i], xout = grid_u, rule = 2)$y
}

Y <- t(sapply(seq_len(nrow(events)), standardise_event))

cat("Y:", nrow(Y), "events x", ncol(Y), "points; peaks aligned at u =", round(u_star, 3), "\n")

# ------------------------------------------------------------------------------
# 9. Reconstruct a real hydrograph
# ------------------------------------------------------------------------------
# Reverses section 8: standardised curve back to hours and m^3/s

reconstruct_hydrograph <- function(y, peak_position, width_h, peak_flow) {
  
  # undo the realignment
  u <- ifelse(grid_u <= u_star,
              grid_u * peak_position / u_star,
              peak_position + (grid_u - u_star) * (1 - peak_position) / (1 - u_star))
  
  # back to real time and flow
  data.frame(
    hours     = u * width_h,
    discharge = y * peak_flow
  )
}

# ------------------------------------------------------------------------------
# 10. Shared functions for checking every method
# ------------------------------------------------------------------------------
# Fitted curves are put back into m^3/s and hours and compared with the data

# Set to FALSE to recompute everything
use_cache <- TRUE

# Returns expr, saved in cache/ under a hash of Y and key
cached <- function(name, expr, key = NULL) {
  if (!use_cache) return(expr)
  file <- file.path("cache",
                    paste0(name, "_", substr(rlang::hash(list(Y, key)), 1, 10), ".rds"))
  if (file.exists(file)) return(readRDS(file))
  result <- expr
  tryCatch({
    dir.create("cache", showWarnings = FALSE)
    saveRDS(result, file)
  }, error = function(e) NULL, warning = function(w) NULL)
  result
}

# Volume by the trapezium rule (points are not evenly spaced in hours)
area_under <- function(hours, flow) {
  sum(diff(hours) * (head(flow, -1) + tail(flow, -1)) / 2)
}

# RMSE, peak, timing and % volume errors for one event
flood_errors <- function(hours, observed, fitted) {
  c(
    rmse         = sqrt(mean((fitted - observed)^2)),
    peak_error   = max(fitted) - max(observed),
    timing_error = hours[which.max(fitted)] - hours[which.max(observed)],
    volume_error = 100 * (area_under(hours, fitted) - area_under(hours, observed)) /
      area_under(hours, observed)
  )
}

# Errors for every event, given fitted curves laid out like Y
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

# Mean absolute value of each error
average_errors <- function(errors) {
  colMeans(abs(errors))
}

# Leave-one-out predictions. fit_left_out(i) fits without event i.
loocv <- function(fit_left_out) {
  Y_pred <- Y
  for (i in 1:nrow(Y)) {
    Y_pred[i, ] <- fit_left_out(i)
  }
  Y_pred
}

# Fixed train/test split, 20% held out
set.seed(2026)
test_floods  <- sort(sample(1:nrow(Y), size = round(0.2 * nrow(Y))))
train_floods <- setdiff(1:nrow(Y), test_floods)

# Observed vs fitted for six contrasting events
sharpness <- rowMeans(Y > 0.5)   # share of time above half the peak

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
    lines(fitted$hours, fitted$discharge, col = "#AA1945", lwd = 2, lty = 2)
    legend("topright", c("Observed", "Fitted"), col = c("black", "#AA1945"),
           lty = c(1, 2), lwd = 2, cex = 0.8)
  }
  par(mfrow = c(1, 1))
}

# Residuals against time, red line is the mean
plot_residuals <- function(Y_fit, title) {
  residuals <- Y - Y_fit
  matplot(grid_u, t(residuals), type = "l", lty = 1, col = "grey",
          xlab = "Rescaled time u", ylab = "Observed - fitted",
          main = paste(title, ": residuals against time"))
  lines(grid_u, colMeans(residuals), col = "#AA1945", lwd = 2)
  abline(h = 0, lty = 2)
}

# The four average errors against a method's setting
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

# Standardised curve back in real units, using event i's size and timing
as_flood <- function(curve, i) {
  reconstruct_hydrograph(curve, events$peak_position[i], events$width_h[i], events$peak_flow[i])
}

# ------------------------------------------------------------------------------
# 11. SEPA's current approach
# -----------------------------------------------------------------------------
# One average shape scaled by peak flow, with no peak realignment

# Each event's curve on its own timing
unaligned_curve <- function(i) {
  real <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                 events$width_h[i], events$peak_flow[i])
  approx(real$hours / events$width_h[i], Y[i, ], xout = grid_u, rule = 2)$y
}

Y_unaligned <- t(sapply(1:nrow(Y), unaligned_curve))

# Average shape for event i, laid out like Y
baseline_curve <- function(average_shape, i) {
  real <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                 events$width_h[i], events$peak_flow[i])
  approx(grid_u, average_shape, xout = real$hours / events$width_h[i], rule = 2)$y
}

# In-sample: every event gets the average of all events
average_shape <- colMeans(Y_unaligned)
average_shape <- average_shape / max(average_shape)   # peak of 1, as SEPA does
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
# Degree-p polynomial fitted to each event with lm(). Coefficients are the
# shape parameters. Degree chosen on the train/test split.

# y ~ I(u^1) + ... + I(u^p), as in the handout
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

# Share of events where the top term has p < 0.05. Too optimistic, as
# neighbouring points are correlated.
prop_significant <- cached("poly_prop_sig", sapply(degrees, function(p) {
                    mean(sapply(1:nrow(Y), function(i) {
                    coefs <- summary(fit_poly_model(Y[i, ], p))$coefficients
                    coefs[nrow(coefs), 4] < 0.05
                    }))}),
                    key = degrees
                    )

print(data.frame(degree = degrees, prop_significant = prop_significant))

# Train and test errors by degree
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

# Diagnostic plots for the worst-fitted event
worst_poly <- which.max(poly_errors[, "rmse"])

par(mfrow = c(2, 2))
plot(poly_models[[worst_poly]])
par(mfrow = c(1, 1))

# Coefficient table for the same event
print(summary(poly_models[[worst_poly]]))

# Correlation between shape parameters
cat("\nCorrelation between polynomial shape parameters:\n")
print(round(cor(poly_coefs), 2))


# ------------------------------------------------------------------------------
# 13. Method 2: principal component analysis (PCA)
# ------------------------------------------------------------------------------
# Each curve = mean curve + k components. Scores are the shape parameters.
# Centred but not scaled.

make_pca <- function(data, log = FALSE) {
  X    <- if (log) log(data) else data      # log needs positive flow
  back <- if (log) exp else identity
  fit  <- prcomp(X, center = TRUE, scale. = FALSE)
  
  # all events rebuilt from the first k components
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

# Loadings of the first three components
plot(grid_u, pca_result$rotation[, 1], type = "l", col = "red",
     ylim = range(pca_result$rotation[, 1:3]),
     xlab = "Rescaled time u", ylab = "Loading",
     main = "Loadings for PC1, PC2 and PC3")
lines(grid_u, pca_result$rotation[, 2], col = "blue")
lines(grid_u, pca_result$rotation[, 3], col = "darkgreen")
abline(h = 0, lty = 3)
legend("bottomright", c("PC1", "PC2", "PC3"),
       col = c("red", "blue", "darkgreen"), lty = 1)

plot(pca_result, main = "Variance of each principal component")

# Cumulative variance explained
cum_var <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)

plot(seq_along(cum_var), cum_var, type = "b", pch = 19, ylim = c(0, 1),
     xlab = "Number of components k", ylab = "Cumulative proportion of variance",
     main = "Cumulative variance explained")
abline(h = 0.95, lty = 2, col = "red")

cat("Smallest k explaining at least 95% of variance:", which(cum_var >= 0.95)[1], "\n")

# Choosing k from in-sample and leave-one-out errors
k_values <- 1:min(10, nrow(Y) - 2)

pca_choice <- cached("pca_choice",
                     insample_loo_choice(k_values, pca_rebuild, pca_left_out),
                     key = k_values)

print(pca_choice, digits = 3)
plot_choice(pca_choice, "Number of components k")

# Fit at the chosen k
k_chosen <- 6

# Rebuilt curves, one row per event
Y_pca <- pca_rebuild(k_chosen)

# Scores, one row per event
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

# Look for outlying events in the scores
pairs(pca_scores, pch = 19, main = "Scores on the kept components")

# Variance and correlation of the scores
cat("\nVariance of each PCA score:\n")
print(round(apply(pca_scores, 2, var), 4))
cat("\nCorrelation between PCA scores:\n")
print(round(cor(pca_scores), 2))


# ------------------------------------------------------------------------------
# 14. Method 3: wavelets
# ------------------------------------------------------------------------------
# DWT of each curve, keeping the smooth coefficients and the detail
# coefficients that are largest on average. la8 filter, 3 levels, periodic.

wav_filter <- "la8"
wav_levels <- 3

# Wavelet transform of one curve
wavelet_transform <- function(y, filter = wav_filter, levels = wav_levels) {
  dwt(as.numeric(y), filter = filter, n.levels = levels, boundary = "periodic")
}

# Mean size of each detail coefficient, by level
average_detail_size <- function(wav_list) {
  lapply(1:length(wav_list[[1]]@W), function(j) {
    rowMeans(sapply(wav_list, function(wt) abs(wt@W[[j]][, 1])))
  })
}

# Rebuild one curve from the largest "keep" detail coefficients
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

# Rebuilt curves for every event
wavelet_fit_all <- function(keep) {
  t(sapply(wav_all, function(wt) wavelet_rebuild(wt, avg_size, keep)))
}

# Event i, with coefficients chosen from the other events
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

# Full rebuild should be exact
cat("Max error of full rebuild:", max(abs(Y[1, ] - as.numeric(idwt(wt_example)))), "\n")

# Choosing how many detail coefficients to keep
keep_values <- c(1, 2, 3, 5, 8, 10, 15, 20)

wav_choice <- cached("wav_choice",
                     insample_loo_choice(keep_values, wavelet_fit_all, wavelet_left_out),
                     key = keep_values)

print(wav_choice, digits = 3)
plot_choice(wav_choice, "Detail coefficients kept")

# Fit at the chosen number
keep_chosen <- 5

# Rebuilt curves, one row per event
Y_wav <- wavelet_fit_all(keep_chosen)

# Shape parameters: smooth plus kept detail coefficients
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

# Sensitivity to filter and number of levels, same number of details kept
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
# y(u) = b + (1 - b) (u / u_p)^m exp(m (1 - u / u_p)), fitted with nls().
# b is baseflow, u_p the peak time and m the sharpness.

gamma_curve <- function(u, b, up, m) {
  b + (1 - b) * (u / up)^m * exp(m * (1 - u / up))
}

# Fit one event, NULL if nls fails
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

# Events where the fit worked
gamma_ok <- which(!sapply(gamma_models, is.null))
cat("Gamma fits that failed to converge:", nrow(Y) - length(gamma_ok), "\n")

# Fitted curves. Failed fits stay as the observed curve and are left out below.
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

# Spread of the parameters and how they relate
cat("\nSummary of gamma-type parameters:\n")
print(summary(gamma_params))

pairs(gamma_params[gamma_ok, ], pch = 19, main = "Gamma-type parameters across events")

cat("\nCorrelation between gamma-type parameters:\n")
print(round(cor(gamma_params[gamma_ok, ]), 2))


# ------------------------------------------------------------------------------
# 16. Method 5: B-splines
# ------------------------------------------------------------------------------
# Cubic B-spline basis fitted with lm() (Ramsay & Silverman). Number of basis
# functions chosen on the train/test split, like the polynomial degree.

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

# Fit at the chosen number
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

# Diagnostic plots for the worst-fitted event
worst_bspline <- which.max(bspline_errors[, "rmse"])

par(mfrow = c(2, 2))
plot(bspline_models[[worst_bspline]])
par(mfrow = c(1, 1))

# Correlation between shape parameters
cat("\nCorrelation between B-spline shape parameters:\n")
print(round(cor(bspline_coefs), 2))

# ------------------------------------------------------------------------------
# 17. Validation of methods
# ------------------------------------------------------------------------------

fits_to_compare <- list(Y_poly, Y_pca, Y_wav, Y_gamma, Y_bspline)
method_names    <- c("Polynomial", "PCA", "Wavelets", "Gamma-type", "B-spline")
method_colours <- c("#9C78BA", "#5298BC", "#CC6679", "#FDD835", "#B6D86F")

# Every method on the six example events from section 10
plot_all_methods <- function() {
  par(mfrow = c(5, length(example_events)), mar = c(2, 2, 2, 0.5), oma = c(3, 3, 0, 0))
  for (j in 1:5) {
    for (i in example_events) {
      observed <- as_flood(Y[i, ], i)
      plot(observed$hours, observed$discharge, type = "l", lwd = 1.5, col = "black",
           xlab = "", ylab = "", main = paste(method_names[j], "- event", i),
           cex.main = 0.9)
      # skip failed gamma fits
      if (j == 4 && !(i %in% gamma_ok)) next
      fitted <- as_flood(fits_to_compare[[j]][i, ], i)
      lines(fitted$hours, fitted$discharge, col = method_colours[j], lwd = 2)
    }
  }
  mtext("Hours since start",   side = 1, outer = TRUE, line = 1)
  mtext("Discharge (m^3/s)",   side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
}

# Residuals against time for every method
plot_all_residuals <- function() {
  # gamma-type uses converged fits only
  residuals_by_method <- list()
  for (j in 1:5) {
    rows <- if (j == 4) gamma_ok else 1:nrow(Y)
    residuals_by_method[[j]] <- Y[rows, ] - fits_to_compare[[j]][rows, ]
  }
  y_range <- range(unlist(residuals_by_method))   # same y-axis on every panel

  par(mfrow = c(2, 3))
  # one panel per method, mean in bold
  for (j in 1:5) {
    matplot(grid_u, t(residuals_by_method[[j]]), type = "l", lty = 1,
            col = adjustcolor("grey40", alpha.f = 0.3), ylim = y_range,
            xlab = "Rescaled time u", ylab = "Observed - fitted", main = method_names[j])
    lines(grid_u, colMeans(residuals_by_method[[j]]), col = method_colours[j], lwd = 2.5)
    abline(h = 0, lty = 2)
  }
  # last panel: all five means
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
# Accuracy, number of shape parameters and negative flows for each method

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

n_smooth_cmp <- length(wav_all[[1]]@V[[wav_levels]])

# Held-out curves. Per-event methods just use their in-sample fits.
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

# One entry per method. inputs is K plus peak flow, duration and peak position.
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
cmp_colours <- c("SEPA (average shape)" = "#174B67", setNames(method_colours, method_names))

# Per-event errors
cmp_err <- map(cmp, function(m) list(
  in_s = method_errors(m$fit_in)[m$rows, , drop = FALSE],
  out  = method_errors(m$fit_out)[m$rows, , drop = FALSE]
))

# Errors in time above u, volume above u (million m^3) and rise/fall RMSE
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

# Worst-fitted event for each method
worst_events <- imap_dfr(cmp_err, function(e, name) {
  rows <- cmp[[name]]$rows
  w    <- rows[which.max(e$out[, "rmse"])]
  tibble(method = name, event = w, rmse = max(e$out[, "rmse"]),
         peak_flow = events$peak_flow[w], width_h = events$width_h[w])
})

# Errors against number of shape parameters K. Smallest K is 2 for the
# polynomial, 4 for B-splines and 9 for wavelets. Gamma-type is fixed at 3.
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

# plotmath labels, drawn with label_parsed
metric_labels <- c(rmse = "RMSE~(m^3/s)", peak_error = "Peak~error~(m^3/s)",
                   volume_error = "Volume~error~('%')")
tradeoff_long <- tradeoff %>%
  pivot_longer(c(rmse, peak_error, volume_error), names_to = "metric") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         method = factor(method, levels = cmp_order))
sepa_ref <- tibble(metric = names(metric_labels),
                   value  = average_errors(baseline_errors_loocv)[names(metric_labels)]) %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels))

# Paired comparisons: bootstrap CI for the mean difference and a Wilcoxon
# signed-rank test. Events are independent by the Robson & Reed criteria.
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
  volume_sepa    = S("SEPA (average shape)", "volume_bias"),   # positive means overestimate
  excess_sepa    = S("SEPA (average shape)", "excess_bias"),
  peak_poly      = S("Polynomial", "peak_out"),
  peak_bspline   = S("B-spline", "peak_out"),
  pca_at_wav_K   = tradeoff$rmse[tradeoff$method == "PCA" & tradeoff$K == ncol(wav_params)],
  gamma_failed   = nrow(Y) - length(gamma_ok),
  wav_minus_pca  = paired_results$diff[paired_results$comparison == "Wavelets minus PCA" &
                                         paired_results$metric == paired_metrics[["rmse"]]]
)

# Bold the lowest value in a column
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
             linetype = "dashed", colour = "#174B67") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  geom_point(data = filter(tradeoff_long, chosen), size = 3.5, shape = 21,
             fill = NA, stroke = 1.3) +
  scale_x_continuous(limits = c(1, max(K_grid)), breaks = seq(2, max(K_grid), by = 2)) +
  scale_y_log10() +
  scale_colour_manual(values = cmp_colours, name = NULL) +
  facet_wrap(~ metric, scales = "free_y", labeller = label_parsed) +
  labs(x = "Number of shape parameters K", y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# Display order as in the text. Row groups are added in the Rmd.
paired_order <- c("PCA minus SEPA (average shape)",
                  "Wavelets minus SEPA (average shape)",
                  "PCA minus B-spline",
                  "Wavelets minus PCA")
paired_tbl <- paired_results %>%
  mutate(comparison = factor(comparison, levels = paired_order)) %>%
  arrange(comparison) %>%
  transmute(
    Comparison = factor(sub(" minus ", " vs ", comparison),
                        levels = sub(" minus ", " vs ", paired_order)),
    Metric     = metric,
    Difference = sprintf("%+.2f", diff),
    CI         = sprintf("(%+.2f, %+.2f)", lower, upper),
    Better     = sprintf("%.0f", first_better),
    p          = ifelse(p_value < 0.001, "$<$0.001", sprintf("%.3f", p_value))
  )

f1 <-function(x, d = 1) formatC(x, format = "f", digits = d)

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
# TODO: Unusual and Interpret columns are typed by hand from Sections 4.1-4.3

# ------------------------------------------------------------------------------
# 19. Uncertainty in the shape parameters
# ------------------------------------------------------------------------------
# Uncertainty for the two chosen methods, PCA and wavelets

# Label the shape parameters: PC1, PC2, ... and S1, ..., D1, ...
colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))
n_smooth <- length(wav_all[[1]]@V[[wav_levels]])
colnames(wav_params) <- c(paste0("S", 1:n_smooth),
                          paste0("D", 1:(ncol(wav_params) - n_smooth)))

params <- list(PCA = pca_scores, Wavelets = wav_params)


# ------------------------------------------------------------------------------
# 19a. Curves from simulated parameters
# ------------------------------------------------------------------------------
# Turn any shape parameters into curves laid out like Y

# PCA: mean curve plus scores times loadings
pca_sim_curves <- function(scores, pca_fit = pca_result) {
  loadings <- pca_fit$rotation[, 1:ncol(scores)]
  sweep(scores %*% t(loadings), 2, pca_fit$center, FUN = "+")
}

# Wavelets: put the coefficients back and invert the transform
wav_sim_curves <- function(params, kept = kept_coefs) {
  t(apply(params, 1, function(p) {
    wt <- wav_all[[1]]                                   # template
    wt@V[[wav_levels]][, 1] <- p[1:n_smooth]
    detail       <- numeric(length(kept))
    detail[kept] <- p[-(1:n_smooth)]
    level <- rep(seq_along(wt@W), sapply(wt@W, nrow))   # level of each detail coefficient
    for (j in seq_along(wt@W)) wt@W[[j]][, 1] <- detail[level == j]
    as.numeric(idwt(wt))
  }))
}

# Should give back the section 13 and 14 curves exactly
cat("PCA check:",      max(abs(pca_sim_curves(pca_scores) - Y_pca)), "\n")
cat("Wavelets check:", max(abs(wav_sim_curves(wav_params) - Y_wav)), "\n")


# ------------------------------------------------------------------------------
# 19b. Looking at the shape parameters (Section 6.1)
# ------------------------------------------------------------------------------
# Normal Q-Q plots of each shape parameter. PCA on row 1, wavelets on rows 2-3.

pca_colour <- "#5298BC"
wav_colour <- "#CC6679"

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

# Is shape related to peak flow, duration and peak position?
flood_size_timing <- cbind(peak_flow     = events$peak_flow,
                           duration      = events$width_h,
                           peak_position = events$peak_position)

pca_dependence <- round(cor(pca_scores, flood_size_timing, method = "spearman"), 2)
wav_dependence <- round(cor(wav_params, flood_size_timing, method = "spearman"), 2)

# Heatmap of the Spearman correlations
plot_dependence <- function() {
  dependence <- rbind(pca_dependence, wav_dependence)
  dependence <- dependence[nrow(dependence):1, ]          # PC1 at the top
  n_params   <- nrow(dependence)
  colours    <- colorRampPalette(c("#B2182B", "white", "#2166AC"))(101)
  values     <- round(dependence, 2)
  values[values == 0] <- 0                                # avoid -0.00
  
  par(mar = c(0.5, 3.5, 2, 0.5))
  image(1:3, 1:n_params, t(dependence), zlim = c(-1, 1), col = colours,
        axes = FALSE, xlab = "", ylab = "")
  abline(h = 0:n_params + 0.5, v = 0:3 + 0.5, col = "white", lwd = 2)   # cell borders
  abline(h = nrow(wav_dependence) + 0.5, col = "grey20", lwd = 1.5)     # PCA / wavelets divider
  text(rep(1:3, each = n_params), rep(1:n_params, 3), sprintf("%.2f", values), cex = 0.8,
       font = ifelse(abs(values) >= 0.5, 2, 1))                         # bold if |r| >= 0.5
  axis(3, at = 1:3, labels = c("Peak flow", "Duration", "Peak position"),
       tick = FALSE, line = -0.8, cex.axis = 0.85)
  axis(2, at = 1:n_params, labels = rownames(dependence), las = 1, tick = FALSE,
       line = -0.5, cex.axis = 0.85)
  par(mar = c(5, 4, 4, 2) + 0.1)
}

# Q-Q plots of the residuals after regressing on size and timing (as in 19c)
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
# Regression prediction plus either multivariate Normal residuals or
# resampled residuals with noise. h from Silverman's rule of thumb.

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

# Bandwidths and share of variance from the noise
h_pca <- silverman_h(pca_residuals)
h_wav <- silverman_h(wav_residuals)

noise_share <- function(h) h^2 / (1 + h^2)
noise_share_pca <- noise_share(h_pca)
noise_share_wav <- noise_share(h_wav)

# The four method and scheme combinations, in panel order
combos <- expand.grid(scheme = names(simulators), method = names(params),
                      stringsAsFactors = FALSE)
combo_names   <- paste0(combos$method, ", ", combos$scheme)
method_colour <- c(PCA = pca_colour, Wavelets = wav_colour)
sim_to_curves <- list(PCA = pca_sim_curves, Wavelets = wav_sim_curves)


# ------------------------------------------------------------------------------
# 19d. Out-of-sample check (Section 6.3)
# ------------------------------------------------------------------------------
# Lecture 5: leave each event out, refit, simulate n_each curves for it and
# score them by coverage, CRPS and negative or above-peak curves

n_each <- 200

# Kept wavelet coefficients and shape parameters from a subset of events
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
  simulated <- curves * events$peak_flow[i]
  observed  <- Y[i, ] * events$peak_flow[i]
  lower <- apply(simulated, 2, quantile, 0.05)
  upper <- apply(simulated, 2, quantile, 0.95)
  c(coverage   = 100 * mean(observed >= lower & observed <= upper),
    crps       = mean(sapply(1:n_grid, function(j) crps(simulated[, j], observed[j]))),
    negative   = 100 * mean(apply(curves, 1, min) < 0),
    above_peak = 100 * mean(apply(curves, 1, max) > 1.01))
}

# SEPA: average shape of the other events, no uncertainty
score_sepa <- function(i) {
  average <- colMeans(Y_unaligned[-i, ])
  shape   <- baseline_curve(average / max(average), i)
  c(coverage = NA, crps = events$peak_flow[i] * mean(abs(shape - Y[i, ])),
    negative = NA, above_peak = NA)
}

# Average scores over all left-out events
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

# One value from the table, for inline numbers in the report
oos <- function(method, scheme, column) {
  out_of_sample[[column]][out_of_sample$method == method & out_of_sample$scheme == scheme]
}

# Table is built in the report chunk "out-of-sample"
options(knitr.kable.NA = "--")

# Held-out 90% bands against the observed curve for three Section 4.1 events
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
# All four combinations fitted to all events

# Simulated standardised curves for one combination
simulate_curves <- function(method, scheme, new_size_timing) {
  shapes <- simulators[[scheme]](params[[method]], flood_size_timing, new_size_timing)
  sim_to_curves[[method]](shapes)
}

# 1000 events with size and timing resampled from the real events
n_sim <- 1000
set.seed(2026)
sim_size_timing <- flood_size_timing[sample(nrow(Y), n_sim, replace = TRUE), ]
final_shapes <- lapply(1:4, function(r)
  simulate_curves(combos$method[r], combos$scheme[r], sim_size_timing))


# Fan charts of the simulated shapes, real events in grey
plot_fan_charts <- function() {
  bands <- lapply(final_shapes, function(s)
    apply(s, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
  y_range <- range(Y, unlist(bands))
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


# Volumes for a typical event from simulated shapes vs SEPA's single shape
typical <- apply(flood_size_timing, 2, median)
typical_size_timing <- matrix(typical, n_sim, 3, byrow = TRUE)

volume_of <- function(shape) {
  flood <- reconstruct_hydrograph(shape, typical["peak_position"], typical["duration"],
                                  typical["peak_flow"])
  area_under(flood$hours, flood$discharge) * 3600 / 1e6   # million m^3
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

# Table is built in the report chunk "volumes"


# Simulate whole events in real units, rescaled to peak at the drawn flow
# TODO: draw peak flow from the Part 1 GPD instead of resampling
draw_size_timing <- function(n) flood_size_timing[sample(nrow(Y), n, replace = TRUE), ]

simulate_events <- function(n, method, scheme) {
  size_timing <- draw_size_timing(n)
  shapes <- simulate_curves(method, scheme, size_timing)
  shapes <- shapes / apply(shapes, 1, max)
  lapply(1:n, function(j) reconstruct_hydrograph(shapes[j, ], size_timing[j, "peak_position"],
                                                 size_timing[j, "duration"],
                                                 size_timing[j, "peak_flow"]))
}

# 30 simulated events per combination over the real events
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
# Lecture 5: resample the events, refit and see how much the 90% band moves
n_boot <- 100

# 5% and 95% curves from one refit on the given rows
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

# refitted_bands[[r]][[b]] is refit b of combination r
refitted_bands <- cached("sampling_check_all", {
  set.seed(2026)
  lapply(1:4, function(r) {
    lapply(1:n_boot, function(b) {
      rows <- sample(nrow(Y), replace = TRUE)
      refitted_band(rows, combos$method[r], simulators[[combos$scheme[r]]])
    })
  })
}, key = list(n_boot, n_sim, k_chosen, wav_filter, wav_levels, keep_chosen, names(simulators)))

# Edges of the 90% band from every refit
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

# Mean band width and how much its edges move, as fractions of peak flow
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
# Table is built in the report chunk "sampling"