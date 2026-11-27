# ==============================================================================
# Statistical Case Studies - Project 1 (Part 2)
# Authors: Anna Scollay, Lewis Seath, Zacharias Onyejiaka
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup & Package Loading
# ------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(knitr)
library(patchwork)
library(wavelets)
library(splines)
library(kableExtra)

if (!isTRUE(getOption("knitr.in.progress"))) pdf(NULL)

# Longer timeout for the SEPA download
options(timeout = 300)
cache_file <- "river_don_parkhill_flow.rds"

# ------------------------------------------------------------------------------
# Section 2: Data and exploratory data analysis
# ------------------------------------------------------------------------------

# Load the data, and keep codes 50/100/200/254, but exclude Suspect (150)
# ------------------------------------------------------------------------------

# Load 25 years of data from the River Don at Parkhill
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
retained_codes <- c(50, 100, 200, 254)   
suspect_code   <- 150                    

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


# Data Description
# ------------------------------------------------------------------------------

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


# Data quality
# ------------------------------------------------------------------------------

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
    latex_options = c("striped", "HOLD_position"),
    font_size = 9,
    full_width = FALSE,
    position = "left"
  ) %>%
  pack_rows("Data description", 1, 9) %>%
  pack_rows("Quality codes", 10, nrow(data_table))


# Peak Identification as in Part 1, using Robson and Reed 1999
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

# Independent peaks at final threshold of 65 m^3/s
u_final <- 65
peaks <- independent_peaks(u_final) %>%
  mutate(event_id = row_number())


# Matrix Construction before standardisation
# ------------------------------------------------------------------------------

win_before_h   <- 48
win_after_h    <- 72
dt_hours       <- 0.25 
rel_time_hours <- seq(-win_before_h, win_after_h, by = dt_hours)

extract_event_hydrograph <- function(event_idx, peak_row) {
  t_pk <- t[peak_row]
  idx <- which(t >= t_pk - win_before_h * 3600 & t <= t_pk + win_after_h * 3600)
  df_win <- river_data[idx, ]
  
  event_t_rel <- as.numeric(difftime(df_win$Timestamp, t_pk, units = "hours"))
  
  # Resample onto the 481-point time grid
  approx(x = event_t_rel, y = df_win$Discharge, xout = rel_time_hours, rule = 2)$y
}

N_events <- nrow(peaks)
peak_rows <- match(peaks$Timestamp, t)

# Y_raw: one row per event, 481 columns
Y_raw <- t(sapply(1:N_events, function(i) extract_event_hydrograph(i, peak_rows[i])))
rownames(Y_raw) <- paste0("Event_", 1:N_events)
colnames(Y_raw) <- paste0("t_", sprintf("%+.2fh", rel_time_hours))

# Unchecked readings analysis
# ------------------------------------------------------------------------------

# Each unchecked burst vs linear interpolation between the readings either side.
# Uses the full record so bursts at a window edge still have anchors.
compare_unchecked <- function(data, unchecked_codes = c(200, 254)) {
  d <- data %>%
    mutate(
      is_unch = QualityCode %in% unchecked_codes,
      burst_id = cumsum(is_unch & !lag(is_unch, default = FALSE)),
      burst_id = if_else(is_unch, burst_id, NA_integer_)
    )
  tnum <- as.numeric(d$Timestamp)
  out <- list()
  
  for (b in na.omit(unique(d$burst_id))) {
    bi <- which(d$burst_id == b)
    idx_before <- min(bi) - 1
    idx_after  <- max(bi) + 1
    
    # No anchor at the start or end of the record, so skip
    if (idx_before < 1 || idx_after > nrow(d)) next
    
    interp <- approx(
      x = tnum[c(idx_before, idx_after)],
      y = d$Discharge[c(idx_before, idx_after)],
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

# Stats, table and scatter plot for flagged readings, where NULL if none
summarise_unchecked_readings <- function(in_windows, n_events, label = "Unchecked") {
  if (nrow(in_windows) == 0) return(NULL)
  
  # Count each reading once if windows overlap
  unique_readings <- distinct(in_windows, Timestamp, .keep_all = TRUE)
  rmse <- sqrt(mean(unique_readings$Diff^2))
  max_ad <- max(abs(unique_readings$Diff))
  max_q  <- max(unique_readings$Discharge_Recorded)
  
  events_lookup <- in_windows %>%
    group_by(Burst_ID) %>%
    summarise(Events = paste(sort(unique(Event_ID)), collapse = ", "),
              .groups = "drop")
  
  burst_tbl <- unique_readings %>%
    group_by(Burst_ID) %>%
    summarise(
      Start  = format(min(Timestamp), "%Y-%m-%d %H:%M"),
      Readings = n(),
      Duration_Mins = first(Burst_Length_Mins),
      Max_Flow = max(Discharge_Recorded),
      RMSE = sqrt(mean(Diff^2)),
      Max_Abs_Diff = max(abs(Diff)),
      .groups = "drop"
    ) %>%
    left_join(events_lookup, by = "Burst_ID") %>%
    transmute(
      Burst = as.character(row_number()),
      Event = Events,
      Start,
      Readings = as.character(Readings),
      Duration_Min = as.character(Duration_Mins),
      Max_Flow = sprintf("%.1f", Max_Flow),
      RMSE = sprintf("%.3f", RMSE),
      Max_Abs_Diff = sprintf("%.3f", Max_Abs_Diff)
    )
  
  overall_row <- tibble(
    Burst = "**All**",
    Event = sprintf("**%d of %d events**", n_distinct(in_windows$Event_ID), n_events),
    Start = "-",
    Readings  = sprintf("**%d**", nrow(unique_readings)),
    Duration_Min = "-",
    Max_Flow = sprintf("**%.1f**", max_q),
    RMSE = sprintf("**%.3f**", rmse),
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
                          .(sprintf("%.3f", rmse)) ~ m^3 * "/s" ~ "| dashed line = 1:1"),
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


# Suspect readings analysis
# ------------------------------------------------------------------------------

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


# Defining event start, end and width, and storing extreme flow event information
# ------------------------------------------------------------------------------

search_days   <- 7     
trough_factor <- 1.05  

# Define end of the falling limb, which closely matches rise_start_row from Part 1
fall_end_row <- function(p, limit) {
  window <- which(t >= t[p] & t <= limit)
  lowest <- min(q[window])
  min(window[q[window] <= trough_factor * lowest])
}

# Search up to 7 days after the peak, or to the next peak if it is sooner
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
    peak_time = t[peak_row],
    start_time = t[start_row],
    end_time = t[end_row],
    peak_flow = q[peak_row],
    time_to_peak_h = as.numeric(difftime(peak_time, start_time, units = "hours")),
    recession_h = as.numeric(difftime(end_time, peak_time, units = "hours")),
    width_h = time_to_peak_h + recession_h,
    peak_position = time_to_peak_h / width_h   
  )


# Plots for unchecked readings analysis
# ------------------------------------------------------------------------------

unchecked_model_in_events <- unchecked_in_windows(
  starts = events$start_time,
  ends = events$end_time,
  peak_times = events$peak_time
)
unchecked_model <- summarise_unchecked_readings(unchecked_model_in_events, nrow(events))

# Unchecked readings table and plot for the report, shown side by side
unchecked_side_tbl <- bind_rows(
  unchecked_model$burst_tbl %>%
    transmute(Event, Start = substr(Start, 1, 10), Readings, Max_Flow, RMSE, Max_Abs_Diff),
  tibble(
    Event = "\\textbf{All}",
    Start = "-",
    Readings = sprintf("\\textbf{%d}", nrow(unchecked_model$unique_readings)),
    Max_Flow = sprintf("\\textbf{%.1f}", unchecked_model$max_q),
    RMSE = sprintf("\\textbf{%.3f}", unchecked_model$rmse),
    Max_Abs_Diff = sprintf("\\textbf{%.3f}", unchecked_model$max_ad)
  )
)

unchecked_side_plot <- unchecked_model$plot + labs(title = NULL, subtitle = NULL)


# Plots for suspect reading analysis
# ------------------------------------------------------------------------------

# Suspect readings within the trough-to-trough events
suspect_model_in_events <- unchecked_in_windows(
  starts = events$start_time,
  ends = events$end_time,
  peak_times = events$peak_time,
  readings = suspect_all
)
suspect_model <- summarise_unchecked_readings(suspect_model_in_events, nrow(events),
                                                   label = "Suspect")
suspect_model_plot <- plot_flagged_events(suspect_model, events$start_time,
                                          events$end_time, events$peak_time, "Suspect")

# Quantities for the report text
n_suspect_model_readings <- if (is.null(suspect_model)) 0 else nrow(suspect_model$unique_readings)
n_suspect_model_events <- if (is.null(suspect_model)) 0 else n_distinct(suspect_model$in_windows$Event_ID)
n_suspect_model_bursts <- if (is.null(suspect_model)) 0 else nrow(suspect_model$burst_tbl)

if (is.null(suspect_model)) {
  cat("No Suspect readings within the modelling events.\n\n")
} else {
  cat(sprintf("Modelling events: %d Suspect readings (%d bursts) in %d events | max %.1f m3/s | RMSE vs interpolation %.3f m3/s\n\n",
              n_suspect_model_readings, n_suspect_model_bursts, n_suspect_model_events,
              suspect_model$max_q, suspect_model$rmse))
}


# Standardising each extreme flow event
# ------------------------------------------------------------------------------

# Define necessary parameters for standardisation
u_star <- mean(events$peak_position)
n_grid <- 64
grid_u <- seq(0, 1, length.out = n_grid)

# To standardise, we rescale time to [0, 1], move the peak to u_star, divide by peak flow,
# then read off at 64 points.
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

# Construct the Y matrix after standardisation
# ------------------------------------------------------------------------------

# Each standardised event is a row of Y
Y <- t(sapply(seq_len(nrow(events)), standardise_event))

cat("Y:", nrow(Y), "events x", ncol(Y), "points; peaks aligned at u =", round(u_star, 3), "\n")


# Exploratory data analysis plots
# ------------------------------------------------------------------------------

# Plotting the events aligned at their peak, from 48 hours before to 72 hours after
overlay_panel <- function(curves, ylab, title) {
  long <- tibble(event = rep(seq_len(nrow(curves)), each = ncol(curves)),
                 hours = rep(rel_time_hours, times = nrow(curves)),
                 flow  = as.vector(t(curves)))
  ggplot() +
    geom_line(data = long, aes(hours, flow, group = event),
              colour = "#5298BC", linewidth = 0.4, alpha = 0.25) +
    geom_line(data = tibble(hours = rel_time_hours, flow = colMeans(curves)),
              aes(hours, flow), colour = "#174B67", linewidth = 1.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#AA1945", linewidth = 0.6) +
    scale_x_continuous(breaks = seq(-48, 72, by = 24)) +
    labs(title = title, x = "Hours from peak", y = ylab) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11))
}

hydrograph_overlay_eda_plot <-
  overlay_panel(Y_raw, expression("Discharge (" * m^3 * "/s)"), NULL) +
  geom_hline(yintercept = u_final, linetype = "dotted", linewidth = 0.7)


# Plotting the distributions and correlations of the standardisation parameters
size_timing_df <- tibble(peak_flow = events$peak_flow, duration = events$width_h,
                         peak_position = events$peak_position)

hist_panel <- function(var, xlab) {
  ggplot(size_timing_df, aes(.data[[var]])) +
    geom_histogram(bins = 20, fill = "#5298BC", colour = "white", alpha = 0.8) +
    labs(x = xlab, y = "Events") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11))
}
scatter_panel <- function(x, y, xlab, ylab) {
  r <- cor(size_timing_df[[x]], size_timing_df[[y]], method = "spearman")
  ggplot(size_timing_df, aes(.data[[x]], .data[[y]])) +
    geom_point(colour = "#174B67", alpha = 0.6, size = 1.6) +
    labs(x = xlab, y = ylab, subtitle = sprintf("Spearman r = %.2f", r)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11))
}

size_timing_plot <-
  (hist_panel("peak_flow", expression("Peak flow (" * m^3 * "/s)")) |
     hist_panel("duration", "Duration (h)") |
     hist_panel("peak_position", "Relative peak position")) /
  (scatter_panel("peak_flow", "duration", expression("Peak flow (" * m^3 * "/s)"),
                 "Duration (h)") |
     scatter_panel("peak_flow", "peak_position", expression("Peak flow (" * m^3 * "/s)"),
                   "Relative peak position") |
     scatter_panel("duration", "peak_position", "Duration (h)", "Relative peak position"))


# Plotting each step of the standardisation
standardise_steps <- map_dfr(seq_len(nrow(events)), function(i) {
  rows <- events$start_row[i]:events$end_row[i]
  p    <- events$peak_position[i]
  u    <- as.numeric(difftime(t[rows], events$start_time[i], units = "hours")) /
    events$width_h[i]
  tibble(event     = i,
         u         = u,
         u_aligned = ifelse(u <= p, u * u_star / p, u_star + (u - p) * (1 - u_star) / (1 - p)),
         flow      = q[rows])
})

standardised_long <- tibble(event = rep(seq_len(nrow(Y)), each = n_grid),
                            u     = rep(grid_u, times = nrow(Y)),
                            flow  = as.vector(t(Y)))

step_panel <- function(data, x, ylab, title) {
  ggplot(data, aes(.data[[x]], flow, group = event)) +
    geom_line(colour = "#5298BC", linewidth = 0.4, alpha = 0.25) +
    labs(title = title, x = "Standardised time u", y = ylab) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11))
}

standardised_curves_plot <-
  step_panel(standardise_steps, "u", expression("Discharge (" * m^3 * "/s)"),
             "1. Time rescaled to [0, 1]") |
  (step_panel(standardise_steps, "u_aligned", expression("Discharge (" * m^3 * "/s)"),
              "2. Peaks aligned at u*") +
     geom_vline(xintercept = u_star, linetype = "dashed", colour = "#AA1945",
                linewidth = 0.6)) |
  (step_panel(standardised_long, "u", "Flow / peak flow", "3. Divided by peak flow") +
     geom_vline(xintercept = u_star, linetype = "dashed", colour = "#AA1945",
                linewidth = 0.6) +
     geom_line(data = tibble(u = grid_u, flow = colMeans(Y)), aes(u, flow),
               inherit.aes = FALSE, colour = "#174B67", linewidth = 1.3))


# Reversing the standardised curve back to a hydrograph in real units
# ------------------------------------------------------------------------------

reconstruct_hydrograph <- function(y, peak_position, width_h, peak_flow) {
  
  # Undo the realignment
  u <- ifelse(grid_u <= u_star,
              grid_u * peak_position / u_star,
              peak_position + (grid_u - u_star) * (1 - peak_position) / (1 - u_star))
  
  # Reverse back to m^3/s and hours
  data.frame(
    hours     = u * width_h,
    discharge = y * peak_flow
  )
}

# ------------------------------------------------------------------------------
# Section 3: Methods for representing hydrograph shape
# ------------------------------------------------------------------------------

# Shared functions for errors, plotting, comparing and validating across each method
# ------------------------------------------------------------------------------

# We cached these functions to speed up the run time of the script
use_cache <- TRUE

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


# Calculate volume by the trapezium rule 
area_under <- function(hours, flow) {
  sum(diff(hours) * (head(flow, -1) + tail(flow, -1)) / 2)
}


# Computes the RMSE, peak, timing and % volume errors for one event
flood_errors <- function(hours, observed, fitted) {
  c(
    rmse         = sqrt(mean((fitted - observed)^2)),
    peak_error   = max(fitted) - max(observed),
    timing_error = hours[which.max(fitted)] - hours[which.max(observed)],
    volume_error = 100 * (area_under(hours, fitted) - area_under(hours, observed)) /
      area_under(hours, observed)
  )
}


# Computes the errors of the fitted against observed, for each event
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


# Compute the mean absolute value of each error
average_errors <- function(errors) {
  colMeans(abs(errors))
}


# Compute the Leave-one-out predictions, where fit_left_out(i) fits without event i included
loocv <- function(fit_left_out) {
  Y_pred <- Y
  for (i in 1:nrow(Y)) {
    Y_pred[i, ] <- fit_left_out(i)
  }
  Y_pred
}


# Plots observed against fitted for six contrasting events (the earliest and latest, the
# shortest and longest, and the sharpest and flattest)
sharpness <- rowMeans(Y > 0.5)   

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


# Plots residuals against time for all extreme flow events
plot_residuals <- function(Y_fit, title) {
  residuals <- Y - Y_fit
  matplot(grid_u, t(residuals), type = "l", lty = 1, col = "grey",
          xlab = "Standardised time u", ylab = "Observed - fitted",
          main = paste(title, ": residuals against time"))
  lines(grid_u, colMeans(residuals), col = "#AA1945", lwd = 2)
  abline(h = 0, lty = 2)
}


# Plots the RMSE, peak error, timing error, and volume error for different numbers of parameters
plot_choice <- function(choice_table, setting_name) {
  types <- unique(choice_table$type)
  error_titles <- list(rmse         = expression("RMSE (" * m^3 * "/s)"),
                       peak_error   = expression("Peak error (" * m^3 * "/s)"),
                       timing_error = "Timing error (h)",
                       volume_error = "Volume error (%)")
  par(mfrow = c(2, 2))
  for (error in names(error_titles)) {
    plot(range(choice_table$setting), range(choice_table[[error]]), type = "n",
         xlab = setting_name, ylab = "Average absolute error", main = error_titles[[error]])
    for (j in seq_along(types)) {
      rows <- choice_table$type == types[j]
      lines(choice_table$setting[rows], choice_table[[error]][rows],
            type = "b", pch = 19, col = j)
    }
    legend("topright", types, col = seq_along(types), lty = 1, cex = 0.7)
  }
  par(mfrow = c(1, 1))
}


# For per-event methods that are a linear model, this computes fitted curves for every event 
fit_all_with <- function(fit_model, setting) {
  t(sapply(1:nrow(Y), function(i) fitted(fit_model(Y[i, ], setting))))
}


# For per-event methods, this computes the errors over all events, as leaving events out doesn't change the fit
all_events_choice <- function(settings, fit_all) {
  map_dfr(settings, function(s) {
    data.frame(setting = s, type = "all events",
               t(average_errors(method_errors(fit_all(s)))))
  })
}


# The in-sample vs leave-one-out errors by number of parameters for the methods that require it
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


# Reconstruct event i, using its own size and timing
as_flood <- function(curve, i) {
  reconstruct_hydrograph(curve, events$peak_position[i], events$width_h[i], events$peak_flow[i])
}


# SEPA's current approach
# -----------------------------------------------------------------------------

# Hours from the peak, up to 7 days either side
sepa_hours <- seq(-search_days * 24, search_days * 24, by = dt_hours)

# Event i's flow divided by its peak flow, against hours from the peak, and
# NA outside the event's own start and end.
scaled_event <- function(i) {
  rows  <- events$start_row[i]:events$end_row[i]
  hours <- as.numeric(difftime(t[rows], events$peak_time[i], units = "hours"))
  approx(hours, q[rows] / events$peak_flow[i], xout = sepa_hours)$y
}

# One row per event
Y_sepa_raw <- t(sapply(1:nrow(Y), scaled_event))

# Average hydrograph over the chosen events, with a peak of 1
sepa_average <- function(rows) {
  average <- colMeans(Y_sepa_raw[rows, , drop = FALSE], na.rm = TRUE)
  average / max(average, na.rm = TRUE)
}

# The average hydrograph placed at event i's peak, laid out like Y
baseline_curve <- function(average_shape, i) {
  real <- as_flood(Y[i, ], i)
  hours_from_peak <- real$hours - events$time_to_peak_h[i]
  approx(sepa_hours, average_shape, xout = hours_from_peak, rule = 2)$y
}

# Compute in-sample, every event gets the average of all events
average_shape <- sepa_average(1:nrow(Y))
Y_baseline <- t(sapply(1:nrow(Y), function(i) baseline_curve(average_shape, i)))
baseline_errors <- method_errors(Y_baseline)

# Compute leave-one-out, every event gets the average of the other events
baseline_errors_loocv <- cached("baseline_loocv_fixed",
                                method_errors(loocv(function(i) {
                                  baseline_curve(sepa_average(-i), i)
                                })))

cat("\nSEPA's approach (average shape), average errors:\n")
print(rbind(in_sample     = average_errors(baseline_errors),
            leave_one_out = average_errors(baseline_errors_loocv)))


# Polynomial regression
# ------------------------------------------------------------------------------

# y ~ I(u^1) + ... + I(u^p), as in the lecture handout
poly_formula <- function(p) {
  as.formula(paste("y ~", paste0("I(u^", 1:p, ")", collapse = " + ")))
}

# Fit one event's curve with a degree-p polynomial
fit_poly_model <- function(y, p) {
  lm(poly_formula(p), data = data.frame(y = y, u = grid_u))
}

# Fitted curves for every event at degree p
poly_fit_all <- function(p) fit_all_with(fit_poly_model, p)

# Range of degree considered
degrees <- 2:15

# Share of events where the highest-order term has p < 0.05
prop_significant <- cached("poly_prop_sig", sapply(degrees, function(p) {
                    mean(sapply(1:nrow(Y), function(i) {
                    coefs <- summary(fit_poly_model(Y[i, ], p))$coefficients
                    coefs[nrow(coefs), 4] < 0.05
                    }))}),
                    key = degrees
                    )

print(data.frame(degree = degrees, prop_significant = prop_significant))

# Errors over all events by degree
poly_choice <- cached("poly_choice",
                      all_events_choice(degrees, poly_fit_all), key = degrees)

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


# Principal component analysis (PCA)
# ------------------------------------------------------------------------------

# Defining a function that fits PCA to the standardised curves and rebuilds events
# from the first K components
# This contains a log option which we previously considered, but have not included
# in our final report
make_pca <- function(data, log = FALSE) {
  X    <- if (log) log(data) else data      
  back <- if (log) exp else identity
  fit  <- prcomp(X, center = TRUE, scale. = FALSE)
  
  # All events rebuilt from the first k components, and the mean curve
  rebuild <- function(k) {
    scores   <- fit$x[, 1:k, drop = FALSE]
    loadings <- fit$rotation[, 1:k, drop = FALSE]
    back(sweep(scores %*% t(loadings), 2, colMeans(X), FUN = "+"))
  }
  
  # Event i rebuilt from components found without it
  left_out <- function(i, k) {
    fit_train  <- prcomp(X[-i, ], center = TRUE, scale. = FALSE)
    mean_train <- colMeans(X[-i, ])
    loadings   <- fit_train$rotation[, 1:k, drop = FALSE]
    scores     <- (X[i, ] - mean_train) %*% loadings
    back(as.vector(mean_train + scores %*% t(loadings)))
  }
  
  list(fit = fit, rebuild = rebuild, left_out = left_out)
}

# Apply our PCA function to Y, and store the results
pca <- make_pca(Y)
pca_result <- pca$fit
pca_rebuild <- pca$rebuild
pca_left_out <- pca$left_out

# Plot the loadings of the first three components
plot(grid_u, pca_result$rotation[, 1], type = "l", col = "red",
     ylim = range(pca_result$rotation[, 1:3]),
     xlab = "Standardised time u", ylab = "Loading",
     main = "Loadings for PC1, PC2 and PC3")
lines(grid_u, pca_result$rotation[, 2], col = "blue")
lines(grid_u, pca_result$rotation[, 3], col = "darkgreen")
abline(h = 0, lty = 3)
legend("bottomright", c("PC1", "PC2", "PC3"),
       col = c("red", "blue", "darkgreen"), lty = 1)

# Plots of the variance of each component
plot(pca_result, main = "Variance of each principal component")

# Calculate the cumulative variance explained by each component
cum_var <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)

# Plot the cumulative variance explained against number of components, with 95% marked
plot(seq_along(cum_var), cum_var, type = "b", pch = 19, ylim = c(0, 1),
     xlab = "Number of components k", ylab = "Cumulative proportion of variance",
     main = "Cumulative variance explained")
abline(h = 0.95, lty = 2, col = "red")

cat("Smallest k explaining at least 95% of variance:", which(cum_var >= 0.95)[1], "\n")

# The range of components to be considered
k_values <- 1:min(10, nrow(Y) - 2)

# Choosing our number of components K, comparing in-sample and leave-one-out errors
pca_choice <- cached("pca_choice",
                     insample_loo_choice(k_values, pca_rebuild, pca_left_out),
                     key = k_values)

print(pca_choice, digits = 3)
plot_choice(pca_choice, "Number of components k")

# Choose the number of components from the plot above
k_chosen <- 6

# Store the rebuilt curves at our chosen number of parameters, with one row per event
Y_pca <- pca_rebuild(k_chosen)

# Store the PCA scores, the shape parameters, with one row per event
pca_scores <- pca_result$x[, 1:k_chosen, drop = FALSE]

# In-sample and leave-one-out errors at the chosen K
pca_errors <- method_errors(Y_pca)
pca_errors_loocv <- cached("pca_loocv",
                           method_errors(loocv(function(i) pca_left_out(i, k_chosen))),
                           key = k_chosen)

cat("\nPCA (k =", k_chosen, "), average errors:\n")
print(rbind(in_sample     = average_errors(pca_errors),
            leave_one_out = average_errors(pca_errors_loocv)))

# Rebuilt flow should never be negative
cat("Grid points with negative rebuilt flow:", sum(Y_pca < 0), "\n")

# Checking the fit and residuals
plot_example_fits(Y_pca, "PCA")
plot_residuals(Y_pca, "PCA")

# Look for outlying events in the scores
pairs(pca_scores, pch = 19, main = "Scores on the kept components")

# Variance and correlation of the scores
cat("\nVariance of each PCA score:\n")
print(round(apply(pca_scores, 2, var), 4))
cat("\nCorrelation between PCA scores:\n")
print(round(cor(pca_scores), 2))


# Wavelets
# ------------------------------------------------------------------------------

# DWT of each curve, where we choose the la8 filter, 3 levels, and periodic
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

# Rebuild one curve from the largest kept detail coefficients
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
plot(wt_example@W[[1]], type = "l", col = "red", xlab = "", ylab = "", main = "Level 1 detail (fine)")
plot(wt_example@W[[2]], type = "l", col = "blue", xlab = "", ylab = "", main = "Level 2 detail")
plot(wt_example@W[[3]], type = "l", col = "darkgreen", xlab = "", ylab = "", main = "Level 3 detail (coarse)")
plot(wt_example@V[[3]], type = "l", col = "black", xlab = "", ylab = "", main = "Level 3 smooth (trend)")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

cat("Max error of full rebuild:", max(abs(Y[1, ] - as.numeric(idwt(wt_example)))), "\n")

# The range of detail coefficients we consider
keep_values <- c(1, 2, 3, 5, 8, 10, 15, 20)

# Choosing the number of detail coefficients to keep, comparing in-sample and leave-one-out erros
wav_choice <- cached("wav_choice",
                     insample_loo_choice(keep_values, wavelet_fit_all, wavelet_left_out),
                     key = keep_values)

print(wav_choice, digits = 3)
plot_choice(wav_choice, "Detail coefficients kept")

# Our choice of detail coefficients to keep
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

# In-sample and leave-one-out errors at the chosen number of detail coefficients
wav_errors <- method_errors(Y_wav)
wav_errors_loocv <- cached("wav_loocv",
                           method_errors(loocv(function(i) wavelet_left_out(i, keep_chosen))),
                           key = keep_chosen)

cat("\nWavelets (", keep_chosen, "detail coefficients ), average errors:\n")
print(rbind(in_sample     = average_errors(wav_errors),
            leave_one_out = average_errors(wav_errors_loocv)))

# Rebuilt flow should never be negative
cat("Grid points with negative rebuilt flow:", sum(Y_wav < 0), "\n")

# Check the fit and residuals
plot_example_fits(Y_wav, "Wavelets")
plot_residuals(Y_wav, "Wavelets")

# Sensitivity table of filter and number of levels, where the same number of details kept
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


# Gamma-type function
# ------------------------------------------------------------------------------

# y(u) = b + (1 - b) (u / u_p)^m exp(m (1 - u / u_p)), fitted with nls().
# b is baseflow, u_p the peak time and m the sharpness.
gamma_curve <- function(u, b, up, m) {
  b + (1 - b) * (u / up)^m * exp(m * (1 - u / up))
}

# Fit one event's curve, which returns NULL is nonlinear least squares fails
fit_gamma_model <- function(y) {
  tryCatch(
    nls(y ~ gamma_curve(u, b, up, m),
        data = data.frame(y = y, u = grid_u),
        start = list(b = max(min(y), 0.01), up = u_star, m = 5),
        algorithm = "port",
        lower = c(b = 0,    up = 0.01, m = 0.1),
        upper = c(b = 0.99, up = 0.99, m = 500)),
    error = function(e) NULL
  )
}

# Fit the gamma models to our standardised curves
gamma_models <- lapply(1:nrow(Y), function(i) fit_gamma_model(Y[i, ]))

# Check how many gamma fits failed to converge
gamma_ok <- which(!sapply(gamma_models, is.null))
cat("Gamma fits that failed to converge:", nrow(Y) - length(gamma_ok), "\n")

# The fitted curves, with one row per event. Failed fits stay as the observed curve 
# and are left our of the errors and parameters below
Y_gamma <- Y
for (i in gamma_ok) {
  Y_gamma[i, ] <- fitted(gamma_models[[i]])
}

# Shape parameters, one row per event, NA for failed fits
gamma_params <- matrix(NA, nrow = nrow(Y), ncol = 3,
                       dimnames = list(NULL, c("b", "up", "m")))
for (i in gamma_ok) {
  gamma_params[i, ] <- coef(gamma_models[[i]])
}

# The gamma errors for the converged fits only
gamma_errors <- method_errors(Y_gamma)[gamma_ok, ]

cat("\nGamma-type function, average errors (converged fits):\n")
print(average_errors(gamma_errors))

# The fitted flow should never be negative
cat("Grid points with negative fitted flow:", sum(Y_gamma[gamma_ok, ] < 0), "\n")

# Check the fit and residuals
plot_example_fits(Y_gamma, "Gamma-type")
plot_residuals(Y_gamma, "Gamma-type")

# Spread of the parameters and their correlations
cat("\nSummary of gamma-type parameters:\n")
print(summary(gamma_params))

pairs(gamma_params[gamma_ok, ], pch = 19, main = "Gamma-type parameters across events")

cat("\nCorrelation between gamma-type parameters:\n")
print(round(cor(gamma_params[gamma_ok, ]), 2))


# B-splines
# ------------------------------------------------------------------------------

# Fit one event's curve with n_basis B-spline basis functions
fit_bspline_model <- function(y, n_basis) {
  lm(y ~ bs(u, df = n_basis), data = data.frame(y = y, u = grid_u))
}

# Fitted curves for every event with n_basis basis functions
bspline_fit_all <- function(n_basis) fit_all_with(fit_bspline_model, n_basis)

# The range of basis functions we consider
n_basis_values <- 4:15

# Choosing the number of basis functions for B-spline
bspline_choice <- cached("bspline_choice",
                         all_events_choice(n_basis_values, bspline_fit_all),
                         key = n_basis_values)

print(bspline_choice, digits = 3)
plot_choice(bspline_choice, "Number of B-spline basis functions")

# Our chosen number of basis functions
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
# Validation and comparison of methods
# ------------------------------------------------------------------------------

# Validation of each method on observed events
# ------------------------------------------------------------------------------

# Store the fits of the five methods
fits_to_compare <- list(Y_poly, Y_pca, Y_wav, Y_gamma, Y_bspline)
method_names <- c("Polynomial", "PCA", "Wavelets", "Gamma-type", "B-spline")
method_colours <- c("#9C78BA", "#5298BC", "#CC6679", "#FDD835", "#B6D86F")

# Plot every method on the six contrasting events we previously defined
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
  mtext(expression("Discharge (" * m^3 * "/s)"), side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
}

# Plot residuals against time for every method
plot_all_residuals <- function() {
  
  # The gamma-type uses converged fits only
  residuals_by_method <- list()
  for (j in 1:5) {
    rows <- if (j == 4) gamma_ok else 1:nrow(Y)
    residuals_by_method[[j]] <- Y[rows, ] - fits_to_compare[[j]][rows, ]
  }
  y_range <- range(unlist(residuals_by_method))   # same y-axis on every panel

  par(mfrow = c(2, 3))

  for (j in 1:5) {
    matplot(grid_u, t(residuals_by_method[[j]]), type = "l", lty = 1,
            col = adjustcolor("grey40", alpha.f = 0.3), ylim = y_range,
            xlab = "Standardised time u", ylab = "Observed - fitted", main = method_names[j])
    lines(grid_u, colMeans(residuals_by_method[[j]]), col = method_colours[j], lwd = 2.5)
    abline(h = 0, lty = 2)
  }
  
  averages <- sapply(residuals_by_method, colMeans)
  matplot(grid_u, averages, type = "l", lty = 1, lwd = 2, col = method_colours,
          xlab = "Standardised time u", ylab = "Average residual", main = "Averages compared")
  abline(h = 0, lty = 2)
  legend("bottomright", method_names, col = method_colours, lwd = 2, cex = 0.7)
  par(mfrow = c(1, 1))
}

# Comparison of all five methods
# ------------------------------------------------------------------------------

# Conducting the comparison, extending the evaluation criteria from Section 3

# Store the number of smooth coefficients at the coarsets wavelet level, as these are always kept
n_smooth_cmp <- length(wav_all[[1]]@V[[wav_levels]])

# The held-out curves for the appropriate methods, while per-event methods just use 
# their in-sample fits.
Y_sepa_loo <- cached("cmp_sepa_loo_curves_fixed", loocv(function(i) {
  baseline_curve(sepa_average(-i), i)
}))
Y_pca_loo <- cached("cmp_pca_loo_curves",
                    loocv(function(i) pca_left_out(i, k_chosen)), key = k_chosen)
Y_wav_loo <- cached("cmp_wav_loo_curves",
                    loocv(function(i) wavelet_left_out(i, keep_chosen)),
                    key = list(keep_chosen, wav_filter, wav_levels))

all_rows <- seq_len(nrow(Y))

# One entry per method, storing the in-sample and held-out fits, the events used, 
# the number of shape parameters K, and the shape parameters themselves
cmp <- list(
  "SEPA (average shape)" = list(fit_in = Y_baseline, fit_out = Y_sepa_loo, rows = all_rows,
                                K = 0, inputs = 1, params = NULL),
  "Gamma-type" = list(fit_in = Y_gamma, fit_out = Y_gamma, rows = gamma_ok,
                      K = 3, inputs = 3 + 3, params = gamma_params[gamma_ok, ]),
  "Polynomial" = list(fit_in = Y_poly, fit_out = Y_poly, rows = all_rows,
                      K = p_chosen + 1, inputs = p_chosen + 1 + 3, params = poly_coefs),
  "B-spline" = list(fit_in = Y_bspline, fit_out = Y_bspline, rows = all_rows,
                      K = n_basis_chosen + 1, inputs = n_basis_chosen + 1 + 3,
                      params = bspline_coefs),
  "PCA" = list(fit_in = Y_pca, fit_out = Y_pca_loo, rows = all_rows,
                      K = k_chosen, inputs = k_chosen + 3, params = pca_scores),
  "Wavelets" = list(fit_in = Y_wav, fit_out = Y_wav_loo, rows = all_rows,
                      K = ncol(wav_params), inputs = ncol(wav_params) + 3,
                      params = wav_params)
)

# Use the method order and colours defined previously in every table and plot below
cmp_order <- names(cmp)
cmp_colours <- c("SEPA (average shape)" = "#174B67", setNames(method_colours, method_names))

# Per-event errors, in-sample and held-out, for each method
cmp_err <- map(cmp, function(m) list(
  in_s = method_errors(m$fit_in)[m$rows, , drop = FALSE],
  out  = method_errors(m$fit_out)[m$rows, , drop = FALSE]
))

# Errors in time above u, volume above u (million m^3) and rise/fall RMSE
flood_features <- function(hours, observed, fitted, u = u_final, dt = 0.25) {
  
  # Put both curves on a common 15-minute grid
  grid_h <- seq(0, max(hours), by = dt)
  obs <- approx(hours, observed, xout = grid_h, rule = 2)$y
  fit <- approx(hours, fitted,   xout = grid_h, rule = 2)$y
  
  # The rising limb runs up to the observed peak
  rising <- grid_h <= hours[which.max(observed)]
  c(
    duration_error = dt * (sum(fit > u) - sum(obs > u)),
    excess_error = (sum(pmax(fit - u, 0)) - sum(pmax(obs - u, 0))) * dt * 3600 / 1e6,
    rise_rmse = sqrt(mean((fit[rising]  - obs[rising])^2)),
    fall_rmse = sqrt(mean((fit[!rising] - obs[!rising])^2))
  )
}

# Computes the flood feature errors for each event, in real units
feature_errors <- function(Y_fit) {
  t(sapply(seq_len(nrow(Y)), function(i) {
    observed <- as_flood(Y[i, ], i)
    fitted   <- as_flood(Y_fit[i, ], i)
    flood_features(observed$hours, observed$discharge, fitted$discharge)
  }))
}

# Feature errors for each method, using the held-out fits
cmp_feat <- map(cmp, function(m) feature_errors(m$fit_out)[m$rows, , drop = FALSE])

# Calculates the largest absolute correlation between a method's shape parameters,
# and returns NA if less than two
max_abs_cor <- function(p) {
  if (is.null(p) || ncol(p) < 2) return(NA_real_)
  r <- cor(p)
  max(abs(r[upper.tri(r)]))
}

# Summary table with one row per method. The absolute errors measure size, and the bias
# columns keep the sign
cmp_summary <- imap_dfr(cmp, function(m, name) {
  e_in <- cmp_err[[name]]$in_s
  e_out <- cmp_err[[name]]$out
  f <- cmp_feat[[name]]
  tibble(
    method = name,
    K  = m$K,
    inputs = m$inputs,
    n_events = length(m$rows),
    rmse_in = mean(e_in[, "rmse"]),
    rmse_out = mean(e_out[, "rmse"]),
    nrmse_out = mean(100 * e_out[, "rmse"] / events$peak_flow[m$rows]),
    peak_out = mean(abs(e_out[, "peak_error"])),
    timing_out = mean(abs(e_out[, "timing_error"])),
    volume_out = mean(abs(e_out[, "volume_error"])),
    peak_bias = mean(e_out[, "peak_error"]),
    timing_bias = mean(e_out[, "timing_error"]),
    volume_bias = mean(e_out[, "volume_error"]),
    rmse_p90 = unname(quantile(e_out[, "rmse"], 0.9)),
    rmse_max = max(e_out[, "rmse"]),
    duration = mean(abs(f[, "duration_error"])),
    excess = mean(abs(f[, "excess_error"])),
    excess_bias = mean(f[, "excess_error"]),
    rise = mean(f[, "rise_rmse"]),
    fall = mean(f[, "fall_rmse"]),
    neg_events = sum(apply(m$fit_out[m$rows, , drop = FALSE], 1, min) < 0),
    max_cor = max_abs_cor(m$params)
  )
})

# Errors against number of shape parameters K. PCA and wavelets use leave-one-out, while
# the per-event methods use their in-sample fits
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

# Circle the chosen K for each method
chosen_K <- c("Polynomial" = p_chosen + 1, "B-spline" = n_basis_chosen + 1,
              "PCA" = k_chosen, "Wavelets" = ncol(wav_params), "Gamma-type" = 3)
tradeoff <- tradeoff %>% mutate(chosen = K == chosen_K[method])

metric_labels <- c(rmse = "RMSE~(m^3*'/s')", peak_error = "Peak~error~(m^3*'/s')",
                   volume_error = "Volume~error~('%')")
tradeoff_long <- tradeoff %>%
  pivot_longer(c(rmse, peak_error, volume_error), names_to = "metric") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         method = factor(method, levels = cmp_order))

# SEPA's leave-one-out errors, which are drawn as a dashed reference line in each panel
sepa_ref <- tibble(metric = names(metric_labels),
                   value  = average_errors(baseline_errors_loocv)[names(metric_labels)]) %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels))

# Paired comparisons, using bootstrap confidence intervals for the mean difference and a Wilcoxon
# signed-rank test. The events are independent by the Robson & Reed criteria in Section 2.
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
    tibble(comparison = paste(p[1], "minus", p[2]),
           metric = paired_metrics[[mt]],
           diff = mean(d),
           lower = unname(quantile(boots, 0.025)),
           upper = unname(quantile(boots, 0.975)),
           first_better = 100 * mean(a < b),
           p_value = suppressWarnings(
             wilcox.test(a, b, paired = TRUE, exact = FALSE)$p.value))
  })
})

# Defining a function to extract numbers for the text
S <- function(m, col) cmp_summary[[col]][cmp_summary$method == m]
cmp_text <- list(
  rmse_sepa = S("SEPA (average shape)", "rmse_out"),
  rmse_pca = S("PCA", "rmse_out"),
  rmse_wav = S("Wavelets", "rmse_out"),
  rmse_pca_in = S("PCA", "rmse_in"),
  red_pca  = 100 * (1 - S("PCA", "rmse_out") / S("SEPA (average shape)", "rmse_out")),
  red_wav = 100 * (1 - S("Wavelets", "rmse_out") / S("SEPA (average shape)", "rmse_out")),
  timing_sepa = S("SEPA (average shape)", "timing_out"),
  volume_sepa = S("SEPA (average shape)", "volume_bias"),   # positive means overestimate
  excess_sepa = S("SEPA (average shape)", "excess_bias"),
  peak_poly = S("Polynomial", "peak_out"),
  peak_bspline = S("B-spline", "peak_out"),
  pca_at_wav_K = tradeoff$rmse[tradeoff$method == "PCA" & tradeoff$K == ncol(wav_params)],
  gamma_failed = nrow(Y) - length(gamma_ok),
  wav_minus_pca = paired_results$diff[paired_results$comparison == "Wavelets minus PCA" &
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


# Displaying the comparisons carried out in tables and figures for our rmd

# Constructing the main comparison table, containing accuracy and fitness for purpose
cmp_main_tbl <- cmp_summary %>%
  transmute(
    Method  = method,
    K = as.character(K),
    Inputs = as.character(inputs),
    rmse_in  = fmt_best(rmse_in, 2),
    rmse_out = fmt_best(rmse_out, 2),
    nrmse = fmt_best(nrmse_out, 1),
    peak = fmt_best(peak_out, 2),
    timing = fmt_best(timing_out, 2),
    volume = fmt_best(volume_out, 2),
    neg = as.character(neg_events),
    cor = ifelse(is.na(max_cor), "--", formatC(max_cor, format = "f", digits = 2))
  )


# Constructing the event-feature table, containing bias, worst cases and flood feature events
cmp_feat_tbl <- cmp_summary %>%
  transmute(
    Method = method,
    peak_bias = sprintf("%+.2f", peak_bias),
    volume_bias = sprintf("%+.1f", ifelse(abs(volume_bias) < 0.05, 0, volume_bias)),   # no "-0.0"
    p90 = fmt_best(rmse_p90, 1),
    worst = fmt_best(rmse_max, 1),
    duration = fmt_best(duration, 1),
    excess = fmt_best(excess, 2),
    excess_bias = sprintf("%+.2f", excess_bias),
    rise = fmt_best(rise, 2),
    fall = fmt_best(fall, 2)
  )


# The per-event held-out errors for the boxplots, scaled by peak flow so events of
# different sizes can be compared
per_event_long <- imap_dfr(cmp_err, function(e, name) {
  rows <- cmp[[name]]$rows
  tibble(method = name,
         `RMSE (% of peak flow)` = 100 * e$out[, "rmse"] / events$peak_flow[rows],
         `Absolute peak error (% of peak flow)` = 100 * abs(e$out[, "peak_error"]) / events$peak_flow[rows],
         `Absolute volume error (%)` = abs(e$out[, "volume_error"]))
}) %>%
  pivot_longer(-method, names_to = "metric") %>%
  mutate(method = factor(method, levels = cmp_order),
         metric = factor(metric, levels = unique(metric)))


# Boxplots of the per-event errors for each method, on a log scale
per_event_long_plot <- ggplot(per_event_long, aes(method, value, fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8, linewidth = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = cmp_colours, guide = "none") +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank())


# Plot of errors against K for each method, with the chosen K circled and SEPA a dashed line
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

# The paired comparison table, where row groups are added in the rmd
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

# Decimal place formatting
f1 <-function(x, d = 1) formatC(x, format = "f", digits = d)

# Final synthesis table, which summarises each method against each criterion
synthesis <- tribble(
  ~Method, ~Accuracy, ~Features, ~Unusual, ~Parsimony, ~Behaviour, ~Valid, ~Interpret,
  
  "SEPA (average shape)",
  sprintf("Poor: %s m$^3$/s", f1(S("SEPA (average shape)", "rmse_out"))),
  sprintf("Peak matched; timing %s h; volume %s\\%%",
          f1(S("SEPA (average shape)", "timing_out")), f1(S("SEPA (average shape)", "volume_out"))),
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
  sprintf("Valid; %d %s failed", cmp_text$gamma_failed, ifelse(cmp_text$gamma_failed == 1, "fit", "fits")),
  "High (baseflow, peak time, sharpness)",
  
  "Polynomial",
  sprintf("Moderate: %s m$^3$/s", f1(S("Polynomial", "rmse_out"))),
  sprintf("Peak error %s m$^3$/s, too low", f1(S("Polynomial", "peak_out"))),
  "Fails multi-peak events; drifts at ends",
  as.character(S("Polynomial", "K")),
  sprintf("Near-collinear (max %s)", f1(S("Polynomial", "max_cor"), 2)),
  sprintf("%d events negative", S("Polynomial", "neg_events")),
  "Low",
  
  "B-spline",
  sprintf("Moderate: %s m$^3$/s", f1(S("B-spline", "rmse_out"))),
  sprintf("Peak error %s m$^3$/s, too low", f1(S("B-spline", "peak_out"))),
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


# ------------------------------------------------------------------------------
# Uncertainty analysis
# ------------------------------------------------------------------------------

# Set-up
# ------------------------------------------------------------------------------

# Label the shape parameters for the PCA and wavelets methods taken forward: 
# PC1, PC2, ... and S1, ..., D1, ...
colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))
n_smooth <- length(wav_all[[1]]@V[[wav_levels]])
colnames(wav_params) <- c(paste0("S", 1:n_smooth),
                          paste0("D", 1:(ncol(wav_params) - n_smooth)))

params <- list(PCA = pca_scores, Wavelets = wav_params)

# Turn any shape parameters into curves laid out like Y

# PCA: mean curve plus scores times loadings
pca_sim_curves <- function(scores, pca_fit = pca_result) {
  loadings <- pca_fit$rotation[, 1:ncol(scores)]
  sweep(scores %*% t(loadings), 2, pca_fit$center, FUN = "+")
}

# Wavelets: put the coefficients back and invert the transform
wav_sim_curves <- function(params, kept = kept_coefs) {
  t(apply(params, 1, function(p) {
    wt <- wav_all[[1]]                                  
    wt@V[[wav_levels]][, 1] <- p[1:n_smooth]
    detail       <- numeric(length(kept))
    detail[kept] <- p[-(1:n_smooth)]
    level <- rep(seq_along(wt@W), sapply(wt@W, nrow))   
    for (j in seq_along(wt@W)) wt@W[[j]][, 1] <- detail[level == j]
    as.numeric(idwt(wt))
  }))
}

cat("PCA check:",      max(abs(pca_sim_curves(pca_scores) - Y_pca)), "\n")
cat("Wavelets check:", max(abs(wav_sim_curves(wav_params) - Y_wav)), "\n")


# Distribution and dependencies of the shape parameters
# ------------------------------------------------------------------------------

pca_colour <- "#5298BC"
wav_colour <- "#CC6679"

# Defining a function for the Normal Q-Q plot of one parameter
qq_panel <- function(values, name, colour) {
  qqnorm(values, pch = 19, cex = 0.7, col = colour, main = name, xlab = "", ylab = "")
  qqline(values, col = "black", lwd = 1.5)
}

# Storing the size and timing for each event
flood_size_timing <- cbind(peak_flow = events$peak_flow,
                           duration = events$width_h,
                           peak_position = events$peak_position)

# Computing the dependencies of the shape parameters of each method with the size and timing,
# using Spearman correlations
pca_dependence <- round(cor(pca_scores, flood_size_timing, method = "spearman"), 2)
wav_dependence <- round(cor(wav_params, flood_size_timing, method = "spearman"), 2)

# Heatmap of the Spearman correlations, where positive is blue and negative is red,
# and correlations with absolute value greater than 0.5 are in bold
plot_dependence <- function() {
  dependence <- rbind(pca_dependence, wav_dependence)
  dependence <- dependence[nrow(dependence):1, ]          
  n_params <- nrow(dependence)
  colours <- colorRampPalette(c("#B2182B", "white", "#2166AC"))(101)
  values <- round(dependence, 2)
  values[values == 0] <- 0                                
  
  par(mar = c(0.5, 3.5, 2, 0.5))
  image(1:3, 1:n_params, t(dependence), zlim = c(-1, 1), col = colours,
        axes = FALSE, xlab = "", ylab = "")
  abline(h = 0:n_params + 0.5, v = 0:3 + 0.5, col = "white", lwd = 2)   
  abline(h = nrow(wav_dependence) + 0.5, col = "grey20", lwd = 1.5)     
  text(rep(1:3, each = n_params), rep(1:n_params, 3), sprintf("%.2f", values), cex = 0.8,
       font = ifelse(abs(values) >= 0.5, 2, 1))                         
  axis(3, at = 1:3, labels = c("Peak flow", "Duration", "Peak position"),
       tick = FALSE, line = -0.8, cex.axis = 0.85)
  axis(2, at = 1:n_params, labels = rownames(dependence), las = 1, tick = FALSE,
       line = -0.5, cex.axis = 0.85)
  par(mar = c(5, 4, 4, 2) + 0.1)
}

# Q-Q plots of the residuals after regressing on size and timing
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
  mtext("Sample quantiles", side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0),
      mgp = c(3, 1, 0), cex.axis = 1)
}

# Residuals against fitted values from the same regressions, to check the
# constant-variance assumption, displayed the same way as the Q-Q plots
pca_fitted <- fitted(lm(pca_scores ~ flood_size_timing))
wav_fitted <- fitted(lm(wav_params ~ flood_size_timing))

plot_residual_fitted <- function() {
  par(mfrow = c(3, 7), mar = c(2.5, 2.5, 2, 0.5), oma = c(2.5, 2.5, 0, 0),
      mgp = c(2, 0.6, 0), cex.axis = 0.8)
  fitted_panel <- function(fitted_values, residuals, name, colour) {
    plot(fitted_values, residuals, pch = 19, cex = 0.6, col = adjustcolor(colour, 0.6),
         main = name, xlab = "", ylab = "")
    abline(h = 0, lty = 2)
  }
  for (j in 1:ncol(pca_residuals)) {
    fitted_panel(pca_fitted[, j], pca_residuals[, j], colnames(pca_residuals)[j], pca_colour)
  }
  plot.new()
  legend("center", c("PCA residuals", "Wavelet residuals"), col = c(pca_colour, wav_colour),
         pch = 19, bty = "n")
  for (j in 1:ncol(wav_residuals)) {
    fitted_panel(wav_fitted[, j], wav_residuals[, j], colnames(wav_residuals)[j], wav_colour)
  }
  mtext("Fitted values", side = 1, outer = TRUE, line = 1)
  mtext("Residuals",     side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0),
      mgp = c(3, 1, 0), cex.axis = 1)
}

# Fitting the statistical model and sampling scheme
# ------------------------------------------------------------------------------

# Defining Silverman's rule of thumb
silverman_h <- function(residuals) {
  (4 / ((ncol(residuals) + 2) * nrow(residuals)))^(1 / (ncol(residuals) + 4))
}

# Fitting the conditional multivariate normal distribution
simulate_conditional_normal <- function(shape_params, size_timing, new_size_timing) {
  regression <- lm(shape_params ~ size_timing)
  expected  <- cbind(1, new_size_timing) %*% coef(regression)
  variation  <- MASS::mvrnorm(nrow(new_size_timing), rep(0, ncol(shape_params)),
                              cov(resid(regression)))
  expected + variation
}

# Fitting the conditional smoothed bootstrap sampling scheme
simulate_conditional_bootstrap <- function(shape_params, size_timing, new_size_timing) {
  regression <- lm(shape_params ~ size_timing)
  residuals <- resid(regression)
  expected <- cbind(1, new_size_timing) %*% coef(regression)
  n_new <- nrow(new_size_timing)
  h <- silverman_h(residuals)
  variation <- (residuals[sample(nrow(residuals), n_new, replace = TRUE), ] +
                   MASS::mvrnorm(n_new, rep(0, ncol(residuals)), h^2 * cov(residuals))) /
    sqrt(1 + h^2)
  expected + variation
}

simulators <- list("Conditional Normal" = simulate_conditional_normal,
                   "Conditional bootstrap" = simulate_conditional_bootstrap)

# Calculate h for each method, using Silverman's rule of thumb
h_pca <- silverman_h(pca_residuals)
h_wav <- silverman_h(wav_residuals)

# Calculate the noise share for each method
noise_share <- function(h) h^2 / (1 + h^2)
noise_share_pca <- noise_share(h_pca)
noise_share_wav <- noise_share(h_wav)

# The four method and scheme combinations, in panel order
combos <- expand.grid(scheme = names(simulators), method = names(params),
                      stringsAsFactors = FALSE)
combo_names <- paste0(combos$method, ", ", combos$scheme)
method_colour <- c(PCA = pca_colour, Wavelets = wav_colour)
sim_to_curves <- list(PCA = pca_sim_curves, Wavelets = wav_sim_curves)


# Evaluation criteria and validating the simulations
# ------------------------------------------------------------------------------

# Number of simulations for the leave-one-out validation
n_each <- 200

# The wavelet coefficients and shape parameters kept from a subset of events
wav_kept_from <- function(rows) {
  sizes <- unlist(average_detail_size(wav_all[rows]))
  sizes >= sort(sizes, decreasing = TRUE)[keep_chosen]
}
wav_params_from <- function(rows, kept) {
  t(sapply(wav_all[rows], function(wt) c(wt@V[[wav_levels]][, 1], unlist(wt@W)[kept])))
}

# Define a function to conduct leave-one-out validation, where each event i is left out
# of the fit and 200 curves are simulated  with its own standardisation parameters
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

# Define the CRPS measure, as explained in the report
crps <- function(x, y) {
  mean(abs(x - y)) - 0.5 * mean(abs(outer(x, x, "-")))
}

# Scores for event i from its simulated standardised curves
score_event <- function(curves, i) {
  simulated <- curves * events$peak_flow[i]
  observed  <- Y[i, ] * events$peak_flow[i]
  lower <- apply(simulated, 2, quantile, 0.05)
  upper <- apply(simulated, 2, quantile, 0.95)
  
  # Defines the remaining evaluation criteria
  c(coverage   = 100 * mean(observed >= lower & observed <= upper),
    crps = mean(sapply(1:n_grid, function(j) crps(simulated[, j], observed[j]))),
    negative = 100 * mean(apply(curves, 1, min) < 0),
    above_peak = 100 * mean(apply(curves, 1, max) > 1.01))
}

# Scores for SEPA's current approach
score_sepa <- function(i) {
  shape <- baseline_curve(sepa_average(-i), i)
  c(coverage = NA, crps = events$peak_flow[i] * mean(abs(shape - Y[i, ])),
    negative = NA, above_peak = NA)
}

# The average scores over all left-out events
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

# Extracts values for the inline code in rmd
oos <- function(method, scheme, column) {
  out_of_sample[[column]][out_of_sample$method == method & out_of_sample$scheme == scheme]
}

options(knitr.kable.NA = "--")

# The held-out 90% bands against the observed curve for three of the contrasting events
# that we previously explored
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
  mtext(expression("Discharge (" * m^3 * "/s)"), side = 2, outer = TRUE, line = 1)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
}


# Simulating new extreme flow events
# ------------------------------------------------------------------------------

# Define a function that simulates the standardised curves for one combination
simulate_curves <- function(method, scheme, new_size_timing) {
  shapes <- simulators[[scheme]](params[[method]], flood_size_timing, new_size_timing)
  sim_to_curves[[method]](shapes)
}

# Simulate 1000 events, where size and timing are resampled from the real events
n_sim <- 1000
set.seed(2026)
sim_size_timing <- flood_size_timing[sample(nrow(Y), n_sim, replace = TRUE), ]
final_shapes <- lapply(1:4, function(r)
  simulate_curves(combos$method[r], combos$scheme[r], sim_size_timing))


# A fan chart for each combination, displaying the median, 25-75% and 5-95% bands,
# alongside the observed events in grey
plot_fan_charts <- function() {
  bands <- lapply(final_shapes, function(s)
    apply(s, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
  y_range <- range(Y, unlist(bands))
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  for (r in 1:4) {
    colour <- method_colour[[combos$method[r]]]
    band   <- bands[[r]]
    plot(grid_u, band[3, ], type = "n", ylim = y_range, main = combo_names[r],
         xlab = "Standardised time u", ylab = "Flow / peak flow")
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


# The volumes for a typical event from simulated shapes, using the median of the
# size and timing parameters 
typical <- apply(flood_size_timing, 2, median)
typical_size_timing <- matrix(typical, n_sim, 3, byrow = TRUE)

volume_of <- function(shape) {
  flood <- reconstruct_hydrograph(shape, typical["peak_position"], typical["duration"],
                                  typical["peak_flow"])
  area_under(flood$hours, flood$discharge) * 3600 / 1e6   # million m^3
}

# SEPA's average hydrograph for the typical event, which is placed at its peak
typical_hours <- reconstruct_hydrograph(rep(1, n_grid), typical["peak_position"],
                                        typical["duration"], typical["peak_flow"])$hours
typical_from_peak <- typical_hours - typical["peak_position"] * typical["duration"]
sepa_flow   <- approx(sepa_hours, average_shape, xout = typical_from_peak, rule = 2)$y *
  typical["peak_flow"]
sepa_volume <- area_under(typical_hours, sepa_flow) * 3600 / 1e6

# Construct the table for total volume of a typical event, for each of the combinations
# and SEPA's single shape
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

# The flood volumes from the recommended combination, which we use in the rmd
pca_boot_volume <- volume_table[volume_table$scheme == "PCA, Conditional bootstrap", ]

# Carrying out the same analysis but keeping only physically valid curves
set.seed(2026)
pca_boot_curves <- simulate_curves("PCA", "Conditional bootstrap", typical_size_timing)
valid_curve     <- apply(pca_boot_curves, 1, min) >= 0 & apply(pca_boot_curves, 1, max) <= 1.01
valid_volumes   <- apply(pca_boot_curves[valid_curve, ], 1, volume_of)
pca_boot_valid_volume <- list(share = 100 * mean(valid_curve),
                              lower = unname(quantile(valid_volumes, 0.05)),
                              upper = unname(quantile(valid_volumes, 0.95)))

# Other uncertainties
# ------------------------------------------------------------------------------

# Sampling error

# Define the number of refits
n_boot <- 100

# Simulate the 5% and 95% curves from one refit on the rows of a chosen method
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

# Conduct the refit for all combinations of methods and schemes
# We cached this, as it was slow to run
refitted_bands <- cached("sampling_check_all", {
  set.seed(2026)
  lapply(1:4, function(r) {
    lapply(1:n_boot, function(b) {
      rows <- sample(nrow(Y), replace = TRUE)
      refitted_band(rows, combos$method[r], simulators[[combos$scheme[r]]])
    })
  })
}, key = list(n_boot, n_sim, k_chosen, wav_filter, wav_levels, keep_chosen, names(simulators)))

# The simulated edges of the 90% band from every refit
plot_sampling_bands <- function() {
  y_range <- range(unlist(refitted_bands))
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  for (r in 1:4) {
    colour <- adjustcolor(method_colour[[combos$method[r]]], 0.2)
    matplot(grid_u, sapply(refitted_bands[[r]], function(b) b[1, ]), type = "l", lty = 1,
            col = colour, ylim = y_range, main = combo_names[r],
            xlab = "Standardised time u", ylab = "Flow / peak flow")
    matlines(grid_u, sapply(refitted_bands[[r]], function(b) b[2, ]), lty = 1, col = colour)
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

# Compute the mean band width for each combination as well as the standard deviations
# of the edges, expressing each as a fraction of peak flow
sampling_table <- data.frame(
  scheme = combo_names,
  band_width = sapply(refitted_bands, function(bands)
    mean(sapply(bands, function(b) mean(b[2, ] - b[1, ])))),
  lower_sd = sapply(refitted_bands, function(bands)
    mean(apply(sapply(bands, function(b) b[1, ]), 1, sd))),
  upper_sd = sapply(refitted_bands, function(bands)
    mean(apply(sapply(bands, function(b) b[2, ]), 1, sd)))
)
print(sampling_table, digits = 3)