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
# 2. Load Data & Retain Unchecked Quality Codes (50, 100, 200, 254)
# ------------------------------------------------------------------------------
if (file.exists(cache_file)) {
  river_data <- readRDS(cache_file)
} else {
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

# Key Modification: Include Unchecked codes (200 & 254) along with 50 & 100
river_data <- river_data %>%
  arrange(Timestamp) %>%
  mutate(QualityCode = as.numeric(QualityCode)) %>%
  filter(QualityCode %in% c(50, 100, 200, 254)) %>%
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

# Construct N x m matrix Y_raw (41 x 481)
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
      Discharge_Unchecked = d$Discharge[bi],
      Discharge_Interp    = interp,
      Diff                = d$Discharge[bi] - interp,
      Burst_Length_Mins   = length(bi) * 15
    )
  }
  bind_rows(out)
}

unchecked_all <- compare_unchecked(river_data)

# Unchecked readings falling inside a set of windows. starts/ends/peak_times
# are POSIXct vectors with one element per event.
unchecked_in_windows <- function(starts, ends, peak_times) {
  map_dfr(seq_along(starts), function(i) {
    unchecked_all %>%
      filter(Timestamp >= starts[i], Timestamp <= ends[i]) %>%
      mutate(Event_ID = i,
             Hours_From_Peak = as.numeric(difftime(Timestamp, peak_times[i],
                                                   units = "hours")))
  })
}

# Statistics, combined table and scatter plot for one set of windows.
# Returns NULL if there are no unchecked readings in the windows.
summarise_unchecked_readings <- function(in_windows, n_events) {
  if (nrow(in_windows) == 0) return(NULL)
  
  # each reading counted once even if two event windows overlap
  unique_readings <- distinct(in_windows, Timestamp, .keep_all = TRUE)
  rmse   <- sqrt(mean(unique_readings$Diff^2))
  max_ad <- max(abs(unique_readings$Diff))
  
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
      RMSE         = sprintf("%.3f", RMSE),
      Max_Abs_Diff = sprintf("%.3f", Max_Abs_Diff)
    )
  
  overall_row <- tibble(
    Burst        = "**All**",
    Event        = sprintf("**%d of %d events**", n_distinct(in_windows$Event_ID), n_events),
    Start        = "—",
    Readings     = sprintf("**%d**", nrow(unique_readings)),
    Duration_Min = "—",
    RMSE         = sprintf("**%.3f**", rmse),
    Max_Abs_Diff = sprintf("**%.3f**", max_ad)
  )
  
  lims <- range(c(unique_readings$Discharge_Interp, unique_readings$Discharge_Unchecked))
  
  plot <- ggplot(unique_readings, aes(x = Discharge_Interp, y = Discharge_Unchecked)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#AA1945") +
    geom_point(colour = "#5298BC", alpha = 0.7, size = 2) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(
      title = "Unchecked readings vs linear interpolation",
      subtitle = sprintf("n = %d | RMSE = %.3f m³/s | dashed line = 1:1",
                         nrow(unique_readings), rmse),
      x = expression("Interpolated discharge (" * m^3 * "/s)"),
      y = expression("Unchecked discharge (" * m^3 * "/s)")
    ) +
    theme_minimal()
  
  list(in_windows = in_windows, unique_readings = unique_readings,
       rmse = rmse, max_ad = max_ad,
       burst_tbl = burst_tbl, overall_row = overall_row,
       table = bind_rows(burst_tbl, overall_row), plot = plot)
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

# Plot: All 41 hydrographs centered at peak arrival (greyed out + bold black mean)
hydrograph_overlay_plot <- ggplot() +
  geom_line(
    data = df_events_long,
    aes(x = hours_from_peak, y = Discharge, group = event_id),
    colour = "#5298BC",
    linewidth = 0.4,
    alpha = 0.5
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
    title = "River Don Hydrograph Shape Variability (N = 41 Events)",
    subtitle = "Centered at peak arrival (t = 0) [-48h, +72h]. Bold black: Mean curve Q̄(t). Unchecked values retained.",
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
# Reverses section 5: turns a standardised curve back into a real flood.
#   1. undo the realignment: move the peak from u_star back to the flood's
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
# each fitted curve is turned back into a real flood and compared
# with the observed flood in m^3/s and hours.

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

# Trapezium rule used for flood volume. Needed because
# after undoing the realignment the points are not evenly spaced in hours.
area_under <- function(hours, flow) {
  sum(diff(hours) * (head(flow, -1) + tail(flow, -1)) / 2)
}

# Errors for one flood:
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

# Errors for every flood, for a matrix of fitted curves laid out like Y
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

# Average size of each error across floods
average_errors <- function(errors) {
  colMeans(abs(errors))
}

# Leave-one-out cross-validation: each flood is left out in turn,
# the method is fitted to the other floods, and the left-out flood is
# predicted. fit_left_out(i) must return the method's curve for flood i
# using only the other floods. u_star is the average over all floods; its
# effect on the left-out flood is negligible
loocv <- function(fit_left_out) {
  Y_pred <- Y
  for (i in 1:nrow(Y)) {
    Y_pred[i, ] <- fit_left_out(i)
  }
  Y_pred
}

# Fixed training and test floods: about 20% held out, chosen once
# with a fixed seed so same split always used
set.seed(2026)
test_floods  <- sort(sample(1:nrow(Y), size = round(0.2 * nrow(Y))))
train_floods <- setdiff(1:nrow(Y), test_floods)

# Plot observed vs fitted floods in real units for six kinds of flood
sharpness <- rowMeans(Y > 0.5)   # share of the event spent above half the peak

plot_example_fits <- function(Y_fit, title) {
  examples <- unique(c(which.min(events$peak_position), which.max(events$peak_position),
                       which.min(events$width_h),       which.max(events$width_h),
                       which.max(sharpness),            which.min(sharpness)))
  par(mfrow = c(2, 3))
  for (i in examples) {
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

# Plot residuals against time for every flood. Red line is average residual.
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

# Fitted curves for every flood, for any per-flood linear model
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

# ------------------------------------------------------------------------------
# 11. SEPA's current approach
# -----------------------------------------------------------------------------
# SEPA uses one average shape for every flood, scaled by peak flow. SEPA does
# not realign peaks, so the average is taken over the floods on their own
# timing. Every method is compared against this.

# Each flood's curve on its own timing on the 64 points
unaligned_curve <- function(i) {
  real <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                 events$width_h[i], events$peak_flow[i])
  approx(real$hours / events$width_h[i], Y[i, ], xout = grid_u, rule = 2)$y
}

Y_unaligned <- t(sapply(1:nrow(Y), unaligned_curve))

# Average shape given to flood i, laid out like Y so it can be checked in the same way as every method
baseline_curve <- function(average_shape, i) {
  real <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                 events$width_h[i], events$peak_flow[i])
  approx(grid_u, average_shape, xout = real$hours / events$width_h[i], rule = 2)$y
}

# In-sample: every flood gets the average of all floods
average_shape <- colMeans(Y_unaligned)
average_shape <- average_shape / max(average_shape)   # rescale so the peak equals the max flow, as SEPA does
Y_baseline    <- t(sapply(1:nrow(Y), function(i) baseline_curve(average_shape, i)))

baseline_errors <- method_errors(Y_baseline)

# Leave-one-out: every flood gets the average of the other floods
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
# Each flood's curve is fitted with a polynomial in time, separately for each flood, using lm():
#   y(u) = theta_0 + theta_1 u + theta_2 u^2 + ... + theta_p u^p + error
# The fitted coefficients theta_0, ..., theta_p are the flood's shape
# parameters. Each flood is fitted on its own, so leave-one-out errors are the
# same as in-sample errors. degree is checked on training and test floods.

# Formula y ~ I(u^1) + I(u^2) + ... + I(u^p), written as in the handout
poly_formula <- function(p) {
  as.formula(paste("y ~", paste0("I(u^", 1:p, ")", collapse = " + ")))
}

# Fit one flood's curve with a degree-p polynomial
fit_poly_model <- function(y, p) {
  lm(poly_formula(p), data = data.frame(y = y, u = grid_u))
}

# Fitted curves for every flood at degree p
poly_fit_all <- function(p) fit_all_with(fit_poly_model, p)

# Choosing the degree
degrees <- 2:10

# Significance of the highest term: proportion of floods where the highest term has p < 0.05.
# Neighbouring points on a curve are correlated, so these p-values are too optimistic
prop_significant <- cached("poly_prop_sig", sapply(degrees, function(p) {
                    mean(sapply(1:nrow(Y), function(i) {
                    coefs <- summary(fit_poly_model(Y[i, ], p))$coefficients
                    coefs[nrow(coefs), 4] < 0.05
                    }))}),
                    key = degrees
                    )

print(data.frame(degree = degrees, prop_significant = prop_significant))

# Errors by degree on training and test floods. the chosen degree
# should work as well on the test floods as on the training floods
poly_choice <- cached("poly_choice",
                      train_test_choice(degrees, poly_fit_all), key = degrees)

print(poly_choice, digits = 3)
plot_choice(poly_choice, "Polynomial degree")

# Fit at the chosen degree
p_chosen <- 7

poly_models <- lapply(1:nrow(Y), function(i) fit_poly_model(Y[i, ], p_chosen))

# Fitted curves, one row per flood
Y_poly <- t(sapply(poly_models, fitted))

# Shape parameters, one row per flood
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

# Standard diagnostic plots for the worst-fitted flood:
# residuals vs fitted, Q-Q, scale-location, leverage. Normality and
# independence matter less here: the aim is a good basis, not inference.
worst_poly <- which.max(poly_errors[, "rmse"])

par(mfrow = c(2, 2))
plot(poly_models[[worst_poly]])
par(mfrow = c(1, 1))

# Coefficient table for the same flood
print(summary(poly_models[[worst_poly]]))

# Correlation between shape parameters across floods: strong
# correlations make modelling their uncertainty harder later
cat("\nCorrelation between polynomial shape parameters:\n")
print(round(cor(poly_coefs), 2))


# ------------------------------------------------------------------------------
# 13. Method 2: principal component analysis (PCA)
# ------------------------------------------------------------------------------
# PCA finds the main shapes in the floods from the data itself.
# Each flood's curve is rebuilt as the mean curve plus k component curves:
#   y_i(u) = mean curve(u) + C_i1 a_1(u) + ... + C_ik a_k(u)
#   a_j  : loadings of component j (a curve over the 64 grid points)
#   C_ij : score of flood i on component j
# The scores C_i1, ..., C_ik are the flood's shape parameters.
# Curves are centred but not scaled, as all points are on the same scale

make_pca <- function(data, log = FALSE) {
  X    <- if (log) log(data) else data      # log version needs strictly positive flow
  back <- if (log) exp else identity
  fit  <- prcomp(X, center = TRUE, scale. = FALSE)
  
  # Rebuilt curves for all floods using the first k components
  rebuild <- function(k) {
    scores   <- fit$x[, 1:k, drop = FALSE]
    loadings <- fit$rotation[, 1:k, drop = FALSE]
    back(sweep(scores %*% t(loadings), 2, colMeans(X), FUN = "+"))
  }
  
  # Flood i rebuilt from components found without it
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
# extra components still help for a new flood
k_values <- 1:min(10, nrow(Y) - 2)

pca_choice <- cached("pca_choice",
                     insample_loo_choice(k_values, pca_rebuild, pca_left_out),
                     key = k_values)

print(pca_choice, digits = 3)
plot_choice(pca_choice, "Number of components k")

# Fit at the chosen k and check the fit.
# Update after looking at sections 9a and 9b

k_chosen <- 6

# Rebuilt curves, one row per flood
Y_pca <- pca_rebuild(k_chosen)

# Shape parameters (scores), one row per flood
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

# Outliers in the scores: floods far from the rest are unusual
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
# Each flood's curve is split by the discrete wavelet transform into a smooth
# trend plus detail at several scales. Most detail coefficients
# are tiny, so only the largest few are kept and the rest set to zero.
# The same detail coefficients are kept for every flood, chosen
# as those that are largest on average across floods.
# Daubechies "d4" filter, 3 levels, periodic boundaries. Y has 64 = 2^6 points, so no padding is needed.
# Shape parameters: the smooth coefficients plus the kept detail coefficients.

wav_filter <- "d4"   # wavelet filter
wav_levels <- 3      # number of levels

# Wavelet transform of one curve
wavelet_transform <- function(y, filter = wav_filter, levels = wav_levels) {
  dwt(as.numeric(y), filter = filter, n.levels = levels, boundary = "periodic")
}

# Average size of each detail coefficient across a set of floods, by level
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

# Wavelet transform of every flood
wav_all  <- lapply(1:nrow(Y), function(i) wavelet_transform(Y[i, ]))
avg_size <- average_detail_size(wav_all)

# Rebuilt curves for every flood, keeping "keep" detail coefficients
wavelet_fit_all <- function(keep) {
  t(sapply(wav_all, function(wt) wavelet_rebuild(wt, avg_size, keep)))
}

# Left-out flood i, keeping the coefficients chosen from the other floods
wavelet_left_out <- function(i, keep) {
  wavelet_rebuild(wav_all[[i]], average_detail_size(wav_all[-i]), keep)
}

# Decomposition of one flood
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

# Rebuilt curves, one row per flood
Y_wav <- wavelet_fit_all(keep_chosen)

# Shape parameters: smooth coefficients + kept detail coefficients
cutoff_chosen     <- sort(unlist(avg_size), decreasing = TRUE)[keep_chosen]
kept_coefs       <- unlist(avg_size) >= cutoff_chosen
wav_params <- t(sapply(wav_all, function(wt) {
  c(wt@V[[wav_levels]][, 1], unlist(wt@W)[kept_coefs])
}))

cat("Wavelet shape parameters per flood:", ncol(wav_params), "\n")

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
# fitted by nonlinear least squares with nls(), separately for each flood.
# Each flood is fitted on its own, so leave-one-out errors are the same as
# in-sample errors.

gamma_curve <- function(u, b, up, m) {
  b + (1 - b) * (u / up)^m * exp(m * (1 - u / up))
}

# Fit one flood's curve
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

# Floods where the fit worked
gamma_ok <- which(!sapply(gamma_models, is.null))
cat("Gamma fits that failed to converge:", nrow(Y) - length(gamma_ok), "\n")

# Fitted curves, one row per flood (failed fits left as the observed curve, and excluded from the errors below)
Y_gamma <- Y
for (i in gamma_ok) {
  Y_gamma[i, ] <- fitted(gamma_models[[i]])
}

# Shape parameters, one row per flood
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

pairs(gamma_params[gamma_ok, ], pch = 19, main = "Gamma-type parameters across floods")

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
# Each flood is fitted on its own, so the number of basis functions is
# checked on training and test floods, as for the polynomial.

# Fit one flood's curve with n_basis B-spline basis functions
fit_bspline_model <- function(y, n_basis) {
  lm(y ~ bs(u, df = n_basis), data = data.frame(y = y, u = grid_u))
}

# Fitted curves for every flood with n_basis basis functions
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

# Fitted curves, one row per flood
Y_bspline <- t(sapply(bspline_models, fitted))

# Shape parameters, one row per flood
bspline_coefs <- t(sapply(bspline_models, coef))
colnames(bspline_coefs) <- paste0("theta_", 0:n_basis_chosen)

bspline_errors <- method_errors(Y_bspline)

cat("\nB-splines (", n_basis_chosen, "basis functions ), average errors:\n")
print(average_errors(bspline_errors))

# Fitted flow should never be negative
cat("Grid points with negative fitted flow:", sum(Y_bspline < 0), "\n")

plot_example_fits(Y_bspline, "B-spline")
plot_residuals(Y_bspline, "B-spline")

# Standard diagnostic plots for the worst-fitted flood
worst_bspline <- which.max(bspline_errors[, "rmse"])

par(mfrow = c(2, 2))
plot(bspline_models[[worst_bspline]])
par(mfrow = c(1, 1))

# Correlation between shape parameters across floods
cat("\nCorrelation between B-spline shape parameters:\n")
print(round(cor(bspline_coefs), 2))

# ------------------------------------------------------------------------------
# 18. Validation of methods
# ------------------------------------------------------------------------------

# Fitted curves from every method on same six contrasting floods

# Same six floods as plot_example_fits() (section 10), with all five methods
# drawn on each flood so they can be compared directly.

fits_to_compare <- list(Y_poly, Y_pca, Y_wav, Y_gamma, Y_bspline)
method_names    <- c("Polynomial", "PCA", "Wavelets", "Gamma-type", "B-spline")
method_colours <- c("#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7")

# Fitted curves from every method on the same six floods as plot_example_fits()
# (section 10). One row per method, one column per flood.
plot_all_methods <- function() {
  examples <- unique(c(which.min(events$peak_position), which.max(events$peak_position),
                       which.min(events$width_h),       which.max(events$width_h),
                       which.max(sharpness),            which.min(sharpness)))
  par(mfrow = c(5, length(examples)), mar = c(2, 2, 2, 0.5), oma = c(3, 3, 0, 0))
  for (j in 1:5) {
    for (i in examples) {
      observed <- as_flood(Y[i, ], i)
      plot(observed$hours, observed$discharge, type = "l", lwd = 1.5, col = "black",
           xlab = "", ylab = "", main = paste(method_names[j], "- event", i),
           cex.main = 0.9)
      # skip the gamma-type curve if its fit failed for this flood
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
  # residuals for each method (only floods where the gamma-type fit worked)
  residuals_by_method <- list()
  for (j in 1:5) {
    rows <- if (j == 4) gamma_ok else 1:nrow(Y)
    residuals_by_method[[j]] <- Y[rows, ] - fits_to_compare[[j]][rows, ]
  }
  y_range <- range(unlist(residuals_by_method))   # same y-axis on every panel
  
  par(mfrow = c(2, 3))
  # one panel per method: every flood faint, average bold
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
# physical check (negative flows). For methods fitted to each flood on its
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

print(all_comparison, digits = 3, row.names = FALSE)

# ------------------------------------------------------------------------------
# 19. Uncertainty in the shape parameters
# ------------------------------------------------------------------------------
# Uncertainty analysis of the two chosen methods: PCA (k_chosen components)
# and wavelets (la8, 3 levels, 8 smooth + 5 detail = 13 shape parameters).
#
# BEFORE RUNNING: in section 14 set wav_filter <- "la8", wav_levels <- 3 and
# keep_chosen <- 5, so that wav_params holds the 13 la8 parameters.
#
# The shape parameters themselves come from earlier sections:
#   PCA      : pca_scores  (section 13, one row per flood, k_chosen columns)
#   Wavelets : wav_params  (section 14, one row per flood, 13 columns)
# To simulate new floods we need a probability distribution for these
# parameters.
#
#   19a. Curves from simulated parameters
#   19b. Looking at the shape parameters: distributions, correlations, and
#        dependence on the size and timing of the flood
#   19c. Simulating new shape parameters: a statistical model and a sampling 
#   scheme
#   19d. Out-of-sample check of the simulated shapes against held-out floods
#   19e. Plots of the out-of-sample check
#   19f. Final simulations (provisionally PCA with the Normal)
#   19g. Sampling uncertainty

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
# Sections 13 and 14 rebuild the observed floods only. These two functions
# turn any shape parameters (one row per flood), e.g. simulated ones from
# 19c, into curves laid out like Y. By default they use the components / kept coefficients found from all
# floods. 19d passes in ones found from training floods.

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
# 19b. Looking at the shape parameters
# ------------------------------------------------------------------------------

# Distributions (1): histogram of each parameter with a Normal curve (same mean
# and standard deviation) drawn on top. The multivariate Normal in 19c
# assumes each parameter roughly follows this curve: look for skew, long
# tails or outlying floods.
plot_param_histograms <- function(shape_params, method_name) {
  par(mfrow = c(ceiling(ncol(shape_params) / 4), 4), mar = c(4, 4, 2, 1))
  for (j in 1:ncol(shape_params)) {
    values <- shape_params[, j]
    hist(values, freq = FALSE, breaks = 10, col = "#DCEEF5", xlab = "",
         main = paste(method_name, colnames(shape_params)[j]))
    curve(dnorm(x, mean(values), sd(values)), add = TRUE, col = "red", lwd = 2)
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

# Distributions (2): Normal Q-Q plot of each parameter: if the parameter is
# Normal, the points lie close to the red line. Curving away at the ends
# means skew or long tails.
plot_param_qq <- function(shape_params, method_name) {
  par(mfrow = c(ceiling(ncol(shape_params) / 4), 4), mar = c(4, 4, 2, 1))
  for (j in 1:ncol(shape_params)) {
    values <- shape_params[, j]
    qqnorm(values, pch = 19, cex = 0.7,
           main = paste(method_name, colnames(shape_params)[j]))
    qqline(values, col = "red", lwd = 2)
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

plot_param_histograms(pca_scores, "PCA")
plot_param_qq(pca_scores, "PCA")
plot_param_histograms(wav_params, "Wavelet")
plot_param_qq(wav_params, "Wavelet")

# Correlations between the parameters (printed in sections 13 and 14). PCA
# scores are uncorrelated by construction; the wavelet parameters are not,
# which the multivariate Normal in 19c has to capture (13 x 14 / 2 = 91
# correlations estimated from 41 floods).
wav_correlations <- cor(wav_params)
image(1:13, 1:13, wav_correlations, zlim = c(-1, 1),
      col = colorRampPalette(c("#AA1945", "white", "#174B67"))(41),
      axes = FALSE, xlab = "", ylab = "", main = "Correlation between wavelet parameters")
axis(1, 1:13, colnames(wav_correlations), las = 2)
axis(2, 1:13, colnames(wav_correlations), las = 2)

# Dependence on the size and timing of the flood. We divided by peak flow and
# rescaled time, so shape is modelled separately from Part 1's peak flow and
# time to peak. That is only justified if the shape parameters are (roughly)
# unrelated to them.
flood_size_timing <- cbind(log_peak_flow = log(events$peak_flow),
                           log_width     = log(events$width_h),
                           peak_position = events$peak_position)

# Spearman correlation (based on ranks, so not driven by a few extreme
# floods). With 41 floods, a correlation above about 0.31 in size would arise
# by chance less than 1 time in 20; with this many correlations a few will
# pass by chance, so look for strong, consistent relationships.
size_timing_correlations <- lapply(params, function(shape_params)
  round(cor(shape_params, flood_size_timing, method = "spearman"), 2))
print(size_timing_correlations)

# Scatter plots of each PCA score against each size/timing variable, with a
# smooth trend line (red)
par(mfrow = c(3, ncol(pca_scores)), mar = c(4, 4, 1, 1))
for (v in colnames(flood_size_timing)) {
  for (pc in colnames(pca_scores)) {
    plot(flood_size_timing[, v], pca_scores[, pc], pch = 19, cex = 0.6, xlab = v, ylab = pc)
    lines(lowess(flood_size_timing[, v], pca_scores[, pc]), col = "red")
  }
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


# ------------------------------------------------------------------------------
# 19c. Simulating new shape parameters
# ------------------------------------------------------------------------------
# Two ways to simulate shape parameters for new floods.
# Each takes the parameters of a set of floods (one row per flood)
# and returns n_sim simulated parameter vectors (one per row).
# MASS::mvrnorm() draws from a multivariate Normal (MASS comes with R).

# 1. Statistical model - multivariate Normal: each parameter is Normal, with
#    the means and correlations estimated from the floods
simulate_normal <- function(shape_params, n_sim) {
  MASS::mvrnorm(n_sim, mu = colMeans(shape_params), Sigma = cov(shape_params))
}

# 2. Sampling scheme - smoothed bootstrap: pick one of the observed floods at
#    random and add a little Normal noise, so new floods are similar to, but
#    not copies of, observed ones. No Normality assumption.
#    noise_scale sets the size of the noise relative to the spread of the
#    floods. It uses Silverman's rule of thumb, the standard default choice
#    (the same rule R's density() uses): less noise when there are more
#    floods. Here it is about 0.65 for PCA and 0.75 for wavelets.
simulate_bootstrap <- function(shape_params, n_sim) {
  n_floods    <- nrow(shape_params)
  n_params    <- ncol(shape_params)
  noise_scale <- (4 / ((n_params + 2) * n_floods))^(1 / (n_params + 4))
  
  picked <- shape_params[sample(n_floods, n_sim, replace = TRUE), ]
  noise  <- MASS::mvrnorm(n_sim, mu = rep(0, n_params),
                          Sigma = noise_scale^2 * cov(shape_params))
  picked + noise
}

# Simulate 1000 parameter vectors for each method, fitted to all floods
n_sim <- 1000
set.seed(2026)
simulated_params <- list(
  PCA      = list(Normal    = simulate_normal(pca_scores, n_sim),
                  Bootstrap = simulate_bootstrap(pca_scores, n_sim)),
  Wavelets = list(Normal    = simulate_normal(wav_params, n_sim),
                  Bootstrap = simulate_bootstrap(wav_params, n_sim))
)

# Check 1: distribution of the first four parameters - observed (black)
# against simulated (blue = Normal, red = bootstrap)
compare_param_distributions <- function(observed, simulated, method_name) {
  par(mfrow = c(1, 4), mar = c(4, 4, 2.5, 1))
  for (j in 1:4) {
    obs_density  <- density(observed[, j])
    norm_density <- density(simulated$Normal[, j])
    boot_density <- density(simulated$Bootstrap[, j])
    plot(obs_density, lwd = 2, xlab = "", main = paste(method_name, colnames(observed)[j]),
         ylim = c(0, max(obs_density$y, norm_density$y, boot_density$y)))
    lines(norm_density, col = "#5298BC", lwd = 2, lty = 2)
    lines(boot_density, col = "#AA1945", lwd = 2, lty = 3)
  }
  legend("topright", c("Observed", "Normal", "Bootstrap"),
         col = c("black", "#5298BC", "#AA1945"), lty = 1:3, lwd = 2, cex = 0.7)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

compare_param_distributions(pca_scores, simulated_params$PCA, "PCA")
compare_param_distributions(wav_params, simulated_params$Wavelets, "Wavelet")

# Check 2: 20 simulated curves next to 20 observed ones. Do the simulated
# floods look like real ones (one peak, no negative flow, similar spread)?
plot_20_curves <- function(curves, title, colour) {
  matplot(grid_u, t(curves[1:20, ]), type = "l", lty = 1, col = colour,
          xlab = "Rescaled time u", ylab = "Flow / peak flow", main = title)
}

par(mfrow = c(2, 3))
plot_20_curves(Y, "Observed", "grey40")
plot_20_curves(pca_sim_curves(simulated_params$PCA$Normal),    "PCA - Normal",    "#5298BC")
plot_20_curves(pca_sim_curves(simulated_params$PCA$Bootstrap), "PCA - Bootstrap", "#5298BC")
plot_20_curves(Y, "Observed", "grey40")
plot_20_curves(wav_sim_curves(simulated_params$Wavelets$Normal),    "Wavelets - Normal",    "#AA1945")
plot_20_curves(wav_sim_curves(simulated_params$Wavelets$Bootstrap), "Wavelets - Bootstrap", "#AA1945")
par(mfrow = c(1, 1))


# ------------------------------------------------------------------------------
# 19d. Out-of-sample check
# ------------------------------------------------------------------------------
# Check whether the fitted distribution works well out
# of sample. For a held-out flood:
#   1. using the training floods only, find the PCA components / wavelet
#      coefficients to keep, and fit the Normal or bootstrap to their
#      shape parameters
#   2. simulate n_sim shapes and turn them into floods using the held-out
#      flood's own peak flow, width and peak time (so only SHAPE is tested)
#   3. compare them with the real held-out flood:
#        rmse, peak_error, timing_error, volume_error : the median simulated
#            flood against the real one (flood_errors(), section 10)
#        coverage90 : share of the real curve inside the simulated 90% band
#            (should be about 0.90; lower = band too narrow)
#        crps : average distance between the simulated and real flood, minus
#            half the average distance between two simulated floods (m^3/s,
#            lower is better). It rewards being close AND having the right
#            spread; for a single curve (SEPA) it equals the mean absolute
#            error.
#        negative : share of simulated curves with negative flow
#        peak_off : share whose highest point is not within 5% of the peak flow
# Done for leave-one-out (every flood held out in turn) and for the training
# and test floods from section 10.

# A standardised curve as a flood in m^3/s and hours, using flood i's own
# peak flow, width and peak time
as_flood <- function(curve, i) {
  reconstruct_hydrograph(curve, events$peak_position[i], events$width_h[i], events$peak_flow[i])
}

# Wavelet coefficients to keep, chosen from the training floods only, and the
# shape parameters of a set of floods using those coefficients
wav_kept_from <- function(train_rows) {
  sizes <- unlist(average_detail_size(wav_all[train_rows]))
  sizes >= sort(sizes, decreasing = TRUE)[keep_chosen]
}

wav_params_from <- function(rows, kept) {
  t(sapply(wav_all[rows], function(wt) c(wt@V[[wav_levels]][, 1], unlist(wt@W)[kept])))
}

# Simulated shapes for a held-out flood, using only the training floods
simulate_shapes <- function(train_rows, method, simulate) {
  if (method == "PCA") {
    pca_fit <- prcomp(Y[train_rows, ], center = TRUE, scale. = FALSE)
    pca_sim_curves(simulate(pca_fit$x[, 1:k_chosen], n_sim), pca_fit)
  } else {
    kept <- wav_kept_from(train_rows)
    wav_sim_curves(simulate(wav_params_from(train_rows, kept), n_sim), kept)
  }
}

# Scores for one held-out flood i
check_held_out <- function(i, train_rows, method, simulate) {
  sim_shapes <- simulate_shapes(train_rows, method, simulate)
  real_shape <- Y[i, ]
  
  real_flood   <- as_flood(real_shape, i)
  median_flood <- as_flood(apply(sim_shapes, 2, median), i)
  
  band_low  <- apply(sim_shapes, 2, quantile, probs = 0.05)
  band_high <- apply(sim_shapes, 2, quantile, probs = 0.95)
  
  half <- n_sim / 2   # compare simulations 1-500 with 501-1000 for the spread
  crps_by_time <- sapply(1:n_grid, function(j) {
    sims <- sim_shapes[, j]
    mean(abs(sims - real_shape[j])) -
      0.5 * mean(abs(sims[1:half] - sims[(half + 1):n_sim]))
  })
  
  c(flood_errors(real_flood$hours, real_flood$discharge, median_flood$discharge),
    crps       = events$peak_flow[i] * mean(crps_by_time),
    coverage90 = mean(real_shape >= band_low & real_shape <= band_high),
    negative   = mean(apply(sim_shapes, 1, min) < 0),
    peak_off   = mean(abs(apply(sim_shapes, 1, max) - 1) > 0.05))
}

# SEPA's approach for comparison: one average shape of the training floods,
# no shape uncertainty
check_held_out_sepa <- function(i, train_rows) {
  average <- colMeans(Y_unaligned[train_rows, ])
  shape   <- baseline_curve(average / max(average), i)
  real_flood <- as_flood(Y[i, ], i)
  sepa_flood <- as_flood(shape, i)
  c(flood_errors(real_flood$hours, real_flood$discharge, sepa_flood$discharge),
    crps = events$peak_flow[i] * mean(abs(shape - Y[i, ])),
    coverage90 = NA, negative = NA, peak_off = NA)
}

# Every model to check: check(i, train_rows) returns the scores for flood i
models_to_check <- list(
  "SEPA average shape" = check_held_out_sepa,
  "PCA - Normal"       = function(i, train_rows) check_held_out(i, train_rows, "PCA", simulate_normal),
  "PCA - Bootstrap"    = function(i, train_rows) check_held_out(i, train_rows, "PCA", simulate_bootstrap),
  "Wavelets - Normal"    = function(i, train_rows) check_held_out(i, train_rows, "Wavelets", simulate_normal),
  "Wavelets - Bootstrap" = function(i, train_rows) check_held_out(i, train_rows, "Wavelets", simulate_bootstrap)
)

# Scores for every held-out flood (one row per flood)
leave_one_out_scores <- function(check) {
  set.seed(2026)
  t(sapply(1:nrow(Y), function(i) check(i, setdiff(1:nrow(Y), i))))
}
test_flood_scores <- function(check) {
  set.seed(2026)
  t(sapply(test_floods, function(i) check(i, train_floods)))
}

# Takes a minute or two
held_out_scores <- cached("uncertainty_held_out",
                          lapply(models_to_check, function(check) list(
                            "leave-one-out" = leave_one_out_scores(check),
                            "test floods"   = test_flood_scores(check))),
                          key = list(k_chosen, wav_filter, wav_levels, keep_chosen, n_sim))

# Summary: average over held-out floods. Errors are averaged as absolute
# values, as in average_errors() (section 10).
uncertainty_summary <- do.call(rbind, lapply(c("leave-one-out", "test floods"), function(type) {
  do.call(rbind, lapply(names(held_out_scores), function(model) {
    scores <- held_out_scores[[model]][[type]]
    data.frame(model = model, check = type,
               t(colMeans(abs(scores[, c("rmse", "peak_error", "timing_error",
                                         "volume_error", "crps")]))),
               t(colMeans(scores[, c("coverage90", "negative", "peak_off")])))
  }))
}))

print(uncertainty_summary, digits = 3, row.names = FALSE)
# ------------------------------------------------------------------------------
# 19e. Plots of the out-of-sample check
# ------------------------------------------------------------------------------
# The simulated 90% band against the real flood, for the same six kinds of
# flood as plot_example_fits() (section 10): earliest and latest peak,
# shortest and longest, sharpest and flattest. Each flood is left out of the
# fit (as in 19d), so this shows how well the model covers a flood it has
# not seen. The real flood (black) should mostly sit inside the band.

plot_held_out_bands <- function(method, simulate, title) {
  examples <- unique(c(which.min(events$peak_position), which.max(events$peak_position),
                       which.min(events$width_h),       which.max(events$width_h),
                       which.max(sharpness),            which.min(sharpness)))
  set.seed(2026)
  par(mfrow = c(2, 3))
  for (i in examples) {
    sim_shapes <- simulate_shapes(setdiff(1:nrow(Y), i), method, simulate)
    
    real_flood   <- as_flood(Y[i, ], i)
    low_flood    <- as_flood(apply(sim_shapes, 2, quantile, probs = 0.05), i)
    median_flood <- as_flood(apply(sim_shapes, 2, median), i)
    high_flood   <- as_flood(apply(sim_shapes, 2, quantile, probs = 0.95), i)
    
    plot(real_flood$hours, real_flood$discharge, type = "n",
         ylim = range(low_flood$discharge, high_flood$discharge, real_flood$discharge),
         xlab = "Hours since start", ylab = "Discharge (m^3/s)",
         main = paste(title, "- event", i))
    polygon(c(low_flood$hours, rev(high_flood$hours)),
            c(low_flood$discharge, rev(high_flood$discharge)),
            col = "#DCEEF5", border = NA)
    lines(median_flood$hours, median_flood$discharge, col = "#5298BC", lwd = 2, lty = 2)
    lines(real_flood$hours, real_flood$discharge, lwd = 2)
    legend("topright", c("Real flood", "Median simulated", "90% band"),
           col = c("black", "#5298BC", NA), lty = c(1, 2, NA), lwd = 2,
           fill = c(NA, NA, "#DCEEF5"), border = NA, cex = 0.7)
  }
  par(mfrow = c(1, 1))
}

plot_held_out_bands("PCA",      simulate_normal,    "PCA - Normal")
plot_held_out_bands("PCA",      simulate_bootstrap, "PCA - Bootstrap")
plot_held_out_bands("Wavelets", simulate_normal,    "Wavelets - Normal")
plot_held_out_bands("Wavelets", simulate_bootstrap, "Wavelets - Bootstrap")


# ------------------------------------------------------------------------------
# 19f. Final simulations
# ------------------------------------------------------------------------------
# The chosen method and way of simulating, fitted to all 41 floods (19c).
# Provisionally PCA with the multivariate Normal.

chosen_method       <- "PCA"      # "PCA" or "Wavelets"   - UPDATE after 19d
chosen_distribution <- "Normal"   # "Normal" or "Bootstrap" - UPDATE after 19d

chosen_simulate <- list(Normal = simulate_normal,
                        Bootstrap = simulate_bootstrap)[[chosen_distribution]]
chosen_curves   <- if (chosen_method == "PCA") pca_sim_curves else wav_sim_curves
chosen_params   <- params[[chosen_method]]

# The n_sim simulated shapes from 19c for the chosen model
final_shapes <- chosen_curves(simulated_params[[chosen_method]][[chosen_distribution]])

# Fan chart: simulated 5-95% and 25-75% bands and median, with the real
# floods in grey
band <- apply(final_shapes, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

plot(grid_u, band[3, ], type = "n", ylim = range(band, Y),
     xlab = "Rescaled time u", ylab = "Flow / peak flow",
     main = paste("Simulated flood shapes:", chosen_method, "-", chosen_distribution))
polygon(c(grid_u, rev(grid_u)), c(band[1, ], rev(band[5, ])), col = "#DCEEF5", border = NA)
polygon(c(grid_u, rev(grid_u)), c(band[2, ], rev(band[4, ])), col = "#9CC7DD", border = NA)
matlines(grid_u, t(Y), col = adjustcolor("grey30", 0.3), lty = 1)
lines(grid_u, band[3, ], col = "#174B67", lwd = 2)
legend("topright", c("5-95%", "25-75%", "Median", "Real floods"),
       fill = c("#DCEEF5", "#9CC7DD", NA, NA), border = NA,
       col = c(NA, NA, "#174B67", "grey30"), lty = c(NA, NA, 1, 1), cex = 0.8)

# Physical checks on the simulated shapes (as in 19d)
cat("Share of simulated shapes with negative flow:",
    round(mean(apply(final_shapes, 1, min) < 0), 3), "\n")
cat("Share whose highest point is not within 5% of the peak flow:",
    round(mean(abs(apply(final_shapes, 1, max) - 1) > 0.05), 3), "\n")

# Simulating whole floods in m^3/s and hours: what would be fed to the flood
# inundation model or emulator. Each simulated flood combines
#   - a simulated shape, divided by its own maximum so the flood peaks at
#     exactly the drawn peak flow (as in SEPA's current approach)
#   - a peak flow, width and peak position for that flood
#
# PLACEHOLDER: draw_flood_size_timing() currently resamples the observed
# floods (which keeps the relationships between the three). REPLACE it with
# the Part 1 distributions (e.g. the GPD for peak flow) for the final analysis.
draw_flood_size_timing <- function(n) {
  events[sample(nrow(events), n, replace = TRUE),
         c("peak_flow", "width_h", "peak_position")]
}

simulate_floods <- function(n) {
  shapes      <- chosen_curves(chosen_simulate(chosen_params, n))
  shapes      <- shapes / apply(shapes, 1, max)
  size_timing <- draw_flood_size_timing(n)
  lapply(1:n, function(j) reconstruct_hydrograph(shapes[j, ], size_timing$peak_position[j],
                                                 size_timing$width_h[j], size_timing$peak_flow[j]))
}

# Plot: 30 simulated floods next to the real floods, in m^3/s and hours
set.seed(2026)
simulated_floods <- simulate_floods(30)
real_floods      <- lapply(1:nrow(Y), function(i) as_flood(Y[i, ], i))

max_hours <- max(sapply(c(real_floods, simulated_floods), function(f) max(f$hours)))
max_flow  <- max(sapply(c(real_floods, simulated_floods), function(f) max(f$discharge)))

par(mfrow = c(1, 2))
plot(NA, xlim = c(0, max_hours), ylim = c(0, max_flow),
     xlab = "Hours since start", ylab = "Discharge (m^3/s)", main = "Real floods")
for (f in real_floods) lines(f$hours, f$discharge, col = adjustcolor("grey30", 0.6))
plot(NA, xlim = c(0, max_hours), ylim = c(0, max_flow),
     xlab = "Hours since start", ylab = "Discharge (m^3/s)",
     main = paste("Simulated floods:", chosen_method, "-", chosen_distribution))
for (f in simulated_floods) lines(f$hours, f$discharge, col = adjustcolor("#5298BC", 0.6))
par(mfrow = c(1, 1))


# ------------------------------------------------------------------------------
# 19g. Sampling uncertainty
# ------------------------------------------------------------------------------
# Our simulations are based on only 41 floods; a different 41 floods would
# give a slightly different model. To see how much this matters:
#   1. resample the 41 floods with replacement (some floods appear more than
#      once, some not at all)
#   2. refit everything on the resampled floods (PCA components / wavelet
#      coefficients to keep, and the Normal or bootstrap) and simulate shapes
#   3. record the simulated 90% band
# Repeated n_boot times. If the bands from the refits are close together,
# the model does not depend much on which floods we happened to observe.

# Uses chosen_distribution from 19f; both methods are shown for comparison.

n_boot <- 100

# Lower (5%) and upper (95%) edges of the simulated band for each refit,
# one row per refit
refitted_bands <- function(method) {
  set.seed(2026)
  lower <- matrix(NA, n_boot, n_grid)
  upper <- matrix(NA, n_boot, n_grid)
  for (b in 1:n_boot) {
    resampled  <- sample(nrow(Y), replace = TRUE)
    sim_shapes <- simulate_shapes(resampled, method, chosen_simulate)
    lower[b, ] <- apply(sim_shapes, 2, quantile, probs = 0.05)
    upper[b, ] <- apply(sim_shapes, 2, quantile, probs = 0.95)
  }
  list(lower = lower, upper = upper)
}

# Takes a minute or two
bootstrap_bands <- cached("uncertainty_sampling",
                          list(PCA      = refitted_bands("PCA"),
                               Wavelets = refitted_bands("Wavelets")),
                          key = list(k_chosen, wav_filter, wav_levels, keep_chosen,
                                     n_sim, n_boot, chosen_distribution))

# Band from all 41 floods (19c), for comparison
full_band <- list(
  PCA      = pca_sim_curves(simulated_params$PCA[[chosen_distribution]]),
  Wavelets = wav_sim_curves(simulated_params$Wavelets[[chosen_distribution]])
)

# Plot: band edges from every refit (blue = lower, red = upper), with the
# band from all floods in black
par(mfrow = c(1, 2))
for (method in names(bootstrap_bands)) {
  bands <- bootstrap_bands[[method]]
  matplot(grid_u, t(rbind(bands$lower, bands$upper)), type = "l", lty = 1,
          col = c(rep(adjustcolor("#5298BC", 0.2), n_boot),
                  rep(adjustcolor("#AA1945", 0.2), n_boot)),
          xlab = "Rescaled time u", ylab = "Flow / peak flow",
          main = paste0(method, ": 90% band, ", n_boot, " refits"))
  lines(grid_u, apply(full_band[[method]], 2, quantile, probs = 0.05), lwd = 2)
  lines(grid_u, apply(full_band[[method]], 2, quantile, probs = 0.95), lwd = 2)
}
par(mfrow = c(1, 1))

# Summary (fractions of the peak flow, averaged over the curve):
#   band_width : average width of the 90% band
#   lower_sd, upper_sd : how much each edge of the band moves between refits
# If lower_sd and upper_sd are small compared with band_width, sampling
# error is a small part of the overall uncertainty.
sampling_summary <- sapply(bootstrap_bands, function(bands) c(
  band_width = mean(bands$upper - bands$lower),
  lower_sd   = mean(apply(bands$lower, 2, sd)),
  upper_sd   = mean(apply(bands$upper, 2, sd))
))
print(round(sampling_summary, 3))
