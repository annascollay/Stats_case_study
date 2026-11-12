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
# 4. Standardized Asymmetric Time Window Setup (Option B: [-48h, +72h])
# ------------------------------------------------------------------------------
win_before_h   <- 48
win_after_h    <- 72
dt_hours       <- 0.25 # 15-minute sampling interval
rel_time_hours <- seq(-win_before_h, win_after_h, by = dt_hours)
m_points       <- length(rel_time_hours)

extract_event_hydrograph <- function(event_idx, peak_row) {
  t_pk <- t[peak_row]
  idx  <- which(t >= t_pk - win_before_h * 3600 & t <= t_pk + win_after_h * 3600)
  df_win <- river_data[idx, ]
  
  event_t_rel <- as.numeric(difftime(df_win$Timestamp, t_pk, units = "hours"))
  
  # Resample onto the standard 481-point relative time grid
  approx(x = event_t_rel, y = df_win$Discharge, xout = rel_time_hours, rule = 2)$y
}

# ------------------------------------------------------------------------------
# 5. Matrix Construction & Diagnostic Output
# ------------------------------------------------------------------------------
N_events  <- nrow(peaks)
peak_rows <- match(peaks$Timestamp, t)

# Construct N x m matrix Y_raw (41 x 481)
Y_raw <- t(sapply(1:N_events, function(i) extract_event_hydrograph(i, peak_rows[i])))
rownames(Y_raw) <- paste0("Event_", 1:N_events)
colnames(Y_raw) <- paste0("t_", sprintf("%+.2fh", rel_time_hours))

# Calculate empirical mean hydrograph vector Q_bar(t)
mean_hydrograph <- colMeans(Y_raw)

# Construct mean-centered matrix Y_centered for Milestone 2 PCA/Basis models
Y_centered <- sweep(Y_raw, 2, mean_hydrograph, FUN = "-")

# Unchecked readings vs interpolation for peak windows
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

# Step 2: keep only unchecked readings that fall inside a peak window.
# These are the ones that end up in Y_raw
unchecked_in_events <- map_dfr(seq_len(N_events), function(i) {
  t_pk <- peaks$Timestamp[i]
  unchecked_all %>%
    filter(Timestamp >= t_pk - win_before_h * 3600,
           Timestamp <= t_pk + win_after_h * 3600) %>%
    mutate(Event_ID = i,
           Hours_From_Peak = as.numeric(difftime(Timestamp, t_pk, units = "hours")))
})

# Step 3: summaries and scatter plot
if (nrow(unchecked_in_events) > 0) {
  
  # Per-burst summary (readings inside windows only)
  df_unchecked_summary <- unchecked_in_events %>%
    group_by(Event_ID, Burst_ID) %>%
    summarise(
      Burst_Start      = format(min(Timestamp), "%Y-%m-%d %H:%M"),
      N_Readings       = n(),
      Duration_Mins    = first(Burst_Length_Mins),
      Max_Abs_Diff_m3s = max(abs(Diff)),
      .groups = "drop"
    )
  
  # Overall agreement statistics. each reading counted once even if two
  # event windows overlap
  unique_readings <- distinct(unchecked_in_events, Timestamp, .keep_all = TRUE)
  rmse   <- sqrt(mean(unique_readings$Diff^2))
  max_ad <- max(abs(unique_readings$Diff))
  
  # Table for agreement statistics
  events_lookup <- unchecked_in_events %>%
    group_by(Burst_ID) %>%
    summarise(Events = paste(sort(unique(Event_ID)), collapse = ", "), .groups = "drop")
  
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
    Event        = sprintf("**%d of %d events**", n_distinct(unchecked_in_events$Event_ID), N_events),
    Start        = "—",
    Readings     = sprintf("**%d**", nrow(unique_readings)),
    Duration_Min = "—",
    RMSE         = sprintf("**%.3f**", rmse),
    Max_Abs_Diff = sprintf("**%.3f**", max_ad)
  )
  
  # Scatter: interpolated (x) vs unchecked (y) at identical timestamps
  lims <- range(c(unique_readings$Discharge_Interp,
                  unique_readings$Discharge_Unchecked))
  
  unchecked_scatter_plot <- ggplot(unique_readings,
                                   aes(x = Discharge_Interp, y = Discharge_Unchecked)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#AA1945") +
    geom_point(colour = "#174B67", alpha = 0.7, size = 2) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(
      title = "Unchecked readings vs linear interpolation (peak windows)",
      subtitle = sprintf("n = %d | RMSE = %.3f m³/s | dashed line = 1:1",
                         nrow(unique_readings), rmse),
      x = expression("Interpolated discharge (" * m^3 * "/s)"),
      y = expression("Unchecked discharge (" * m^3 * "/s)")
    ) +
    theme_minimal()
  
  print(unchecked_scatter_plot)
  
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


# Mean hydrograph and variation plot 
shape_summary <- df_events_long %>%
  group_by(hours_from_peak) %>%
  summarise(
    mean_flow = mean(Discharge),
    lower = quantile(Discharge, 0.10),
    upper = quantile(Discharge, 0.90),
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