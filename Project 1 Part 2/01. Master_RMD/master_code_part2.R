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

# Storage for diagnostic plots comparing unchecked values against interpolation
unchecked_diagnostic_plots <- list()
unchecked_summary_list <- list()

extract_event_hydrograph <- function(event_idx, peak_row) {
  t_pk <- t[peak_row]
  t_start_win <- t_pk - win_before_h * 3600
  t_end_win   <- t_pk + win_after_h * 3600
  
  idx <- which(t >= t_start_win & t <= t_end_win)
  df_win <- river_data[idx, ]
  
  event_t_rel <- as.numeric(difftime(df_win$Timestamp, t_pk, units = "hours"))
  event_q     <- df_win$Discharge
  
  # Identify distinct continuous bursts of unchecked quality codes (200/254)
  df_win <- df_win %>%
    mutate(
      is_unchecked = QualityCode %in% c(200, 254),
      # Assign a unique ID to each distinct contiguous burst
      burst_id = cumsum(is_unchecked & !lag(is_unchecked, default = FALSE)),
      burst_id = if_else(is_unchecked, burst_id, NA_integer_)
    )
  
  # Initialize gap-only interpolation column
  df_win$Discharge_Interp <- NA_real_
  
  unique_bursts <- unique(na.omit(df_win$burst_id))
  
  if (length(unique_bursts) > 0) {
    for (b in unique_bursts) {
      burst_indices <- which(df_win$burst_id == b)
      idx_start <- min(burst_indices)
      idx_end   <- max(burst_indices)
      
      # Identify bounding checked points immediately before and after this burst
      idx_before <- max(1, idx_start - 1)
      idx_after  <- min(nrow(df_win), idx_end + 1)
      gap_idx    <- idx_before:idx_after
      
      # Interpolate strictly across the gap [idx_before, idx_after]
      interp_vals <- approx(
        x = as.numeric(df_win$Timestamp[c(idx_before, idx_after)]),
        y = df_win$Discharge[c(idx_before, idx_after)],
        xout = as.numeric(df_win$Timestamp[gap_idx]),
        rule = 2
      )$y
      
      # Assign interpolated values ONLY across this gap
      df_win$Discharge_Interp[gap_idx] <- interp_vals
      
      # Metrics for this specific burst
      burst_rows <- df_win[burst_indices, ]
      burst_duration_mins <- nrow(burst_rows) * 15
      
      # Difference between raw unchecked readings and interpolation
      diffs <- abs(df_win$Discharge[gap_idx] - df_win$Discharge_Interp[gap_idx])
      max_diff <- max(diffs[gap_idx %in% burst_indices], na.rm = TRUE)
      
      # Store burst summary
      unchecked_summary_list[[length(unchecked_summary_list) + 1]] <<- tibble(
        Event_ID = event_idx,
        Burst_ID = b,
        Burst_Start = format(df_win$Timestamp[idx_start], "%Y-%m-%d %H:%M"),
        Duration_Mins = burst_duration_mins,
        Max_Diff_m3s = max_diff
      )
      
      # Dedicated plot zoomed tightly around THIS specific burst (+/- 3 hours)
      t_start_b  <- df_win$Timestamp[idx_start]
      t_end_b    <- df_win$Timestamp[idx_end]
      zoom_start <- t_start_b - hours(3)
      zoom_end   <- t_end_b + hours(3)
      
      p_diag <- ggplot(df_win, aes(x = Timestamp)) +
        # Raw continuous river flow series
        geom_line(aes(y = Discharge, colour = "Raw Series"), linewidth = 0.6) +
        # Dashed line drawn strictly across the gap
        geom_line(aes(y = Discharge_Interp, colour = "Interpolated (Checked Boundary)"), 
                  linetype = "dashed", linewidth = 1.0, na.rm = TRUE) +
        # Highlight unchecked points for this burst
        geom_point(data = burst_rows, 
                   aes(y = Discharge, fill = "Unchecked Reading (Code 200/254)"), 
                   shape = 21, size = 3, colour = "black") +
        coord_cartesian(xlim = c(zoom_start, zoom_end)) +
        scale_colour_manual(values = c("Raw Series" = "#5298BC", 
                                       "Interpolated (Checked Boundary)" = "#E58FA6")) +
        scale_fill_manual(values = c("Unchecked Reading (Code 200/254)" = "#AA1945")) +
        labs(
          title = sprintf("Event %d (%s): Burst %d Diagnostic", 
                          event_idx, format(t_pk, "%d %b %Y"), b),
          subtitle = sprintf("Burst Duration: %d mins | Max Diff vs Linear Interp: %.3f m³/s", 
                             burst_duration_mins, max_diff),
          x = "Timestamp",
          y = expression("Discharge (" * m^3 * "/s)"),
          colour = "Series",
          fill = "Point Status"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      unchecked_diagnostic_plots[[paste0("Event_", event_idx, "_Burst_", b)]] <<- p_diag
    }
  }
  
  # Return interpolated curve on standard 481-point relative time grid
  approx(x = event_t_rel, y = event_q, xout = rel_time_hours, rule = 2)$y
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

# Print Summary of Unchecked Values found within Peak Windows
if (length(unchecked_summary_list) > 0) {
  df_unchecked_summary <- bind_rows(unchecked_summary_list)
  cat("Unchecked Value Events Summary:\n")
  print(kable(df_unchecked_summary, digits = 3))
  cat("\nDiagnostic comparison plots saved in list `unchecked_diagnostic_plots`.\n\n")
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

# Example: Display diagnostic plot for the first event with unchecked readings (if any exist)
if (length(unchecked_diagnostic_plots) > 0) {
  print(unchecked_diagnostic_plots[[1]])
}


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