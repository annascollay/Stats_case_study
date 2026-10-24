# ============================================================
# Statistical Case Studies - Project 1 (Part 1)
# Authors: Anna Scollay, Lewis Seath, Zacharias Onyejiaka
# ============================================================

# --------------------------------------------------------------
# Set-up
# --------------------------------------------------------------
  
# Load packages
library(tidyverse)
library(knitr)
library(kableExtra)
library(ismev)

# --------------------------------------------------------------
# Section 2 - Data and Exploratory Analysis
# --------------------------------------------------------------

# Load data
# --------------------------------------------------------------

# Extend R's default HTTP timeout limit to 300 seconds
options(timeout = 300)

cache_file <- "river_don_parkhill_flow.rds"

if (file.exists(cache_file)) {
  # Load from local cache if file already exists
  river_data <- readRDS(cache_file)
  
} else {
  # If no local cache exists, fetch from SEPA in 1-year chunks
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
    
    read_delim(
      url, 
      delim = ";", 
      comment = "#",
      col_names = c("Timestamp", "Discharge", "QualityCode"),
      show_col_types = FALSE
    )
  }
  
  river_data <- map_dfr(start_years, fetch_chunk) %>%
    mutate(Timestamp = as.POSIXct(Timestamp, format = "%Y-%m-%dT%H:%M:%S")) %>%
    filter(!is.na(Discharge)) %>%
    distinct(Timestamp, .keep_all = TRUE)
  
  # Save to local file so future runs don't need to use the API
  saveRDS(river_data, cache_file)
}

# Subsection 2.1 - Filtering and Cleaning
# --------------------------------------------------------------

# Data Description

# Subsection 2.1 - Filtering and Cleaning
# --------------------------------------------------------------

# Data Description

# Count the number of missing readings by comparing our river data against a
# complete 15-minute sequence
complete_sequence <- seq(min(river_data$Timestamp), max(river_data$Timestamp), by = "15 min")
number_missing <- sum(!complete_sequence %in% river_data$Timestamp)

# Create data description table
data_description_table <- tibble(
  Item = c("Station", "Source", "Measurement", "Resolution",
           "Period", "Readings", "Missing readings"),
  Value = c(
    "Don at Parkhill",
    "SEPA time series service (KiWIS)",
    "River discharge (m³/s)",
    "15 minutes",
    paste(format(min(river_data$Timestamp), "%d %B %Y"), "–",
          format(max(river_data$Timestamp), "%d %B %Y")),
    format(nrow(river_data), big.mark = ","),
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
  Used = c("retained", "retained", "excluded", "excluded", "excluded", "excluded")
)

# Join quality code table with the number of readings for each code
quality_combined <- quality_code_table %>%
  left_join(
    river_data %>% count(QualityCode = as.numeric(QualityCode), name = "Count"),
    by = c("Code number" = "QualityCode")
  ) %>%
  mutate(
    Count = coalesce(Count, 0L),
    Percentage = sprintf("%.2f%%", Count / sum(Count) * 100)
  )

# Combine into one table, with one row per quality code below the data description
data_table <- bind_rows(
  data_description_table,
  quality_combined %>%
    transmute(
      Item  = paste0(`Code number`, " – ", `Code name`, " (", `Abbreviated code name`, ")"),
      Value = paste0(format(Count, big.mark = ",", trim = TRUE), " readings (",
                     Percentage, "), ", Used)
    )
)

# Create a table describing the data
tidy_data_table <- kable(data_table, col.names = c("", ""), caption = "Data used in this analysis, and the quality codes assigned to each reading") %>%
  kable_styling(
    bootstrap_options = c("bordered", "striped", "hover"),
    font_size = 9,
    full_width = FALSE,
    position = "left"
  ) %>%
  pack_rows("Data description", 1, 7) %>%
  pack_rows("Quality codes", 8, nrow(data_table))
# Filter for Quality Codes 50 and 100, joining with quality_code_table labels
qc_plot_data <- river_data %>%
  mutate(QualityCode = as.numeric(QualityCode)) %>%
  filter(QualityCode %in% c(50, 100)) %>%
  left_join(
    quality_code_table %>% select(`Code number`, `Code name`),
    by = c("QualityCode" = "Code number")
  ) %>%
  mutate(
    Quality_Label = factor(
      paste0("Code ", QualityCode, " (", `Code name`, ")"),
      levels = c("Code 50 (Good)", "Code 100 (Estimated)")
    )
  )

# Generate discharge by quality code box-plot
qc_plot <- ggplot(qc_plot_data, aes(x = Quality_Label, y = Discharge, 
                                    fill = Quality_Label, colour = Quality_Label)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1, outlier.alpha = 0.3) +
  scale_fill_manual(values = c(
    "Code 50 (Good)" = "#5298BC",
    "Code 100 (Estimated)" = "#D0436A"
  )) +
  scale_colour_manual(values = c(
    "Code 50 (Good)" = "#2A5267",
    "Code 100 (Estimated)" = "#AA1945"
  )) +
  labs(
    title = "Discharge distribution by quality code",
    x = NULL,
    y = expression("Discharge (" * m^3 * "/s)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

# Keep only good (50) and estimated (100) readings, matching the Used column above
river_data <- river_data %>%
  filter(as.numeric(QualityCode) %in% c(50, 100))

# Subsection 2.2 - Exploratory Analysis
# --------------------------------------------------------------

# Sort by time and create month, season and water year variables to support our analysis
river_data <- river_data %>%
  arrange(Timestamp) %>%
  mutate(
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

# Shorter names for discharge and time, used in the peak analysis
q <- river_data$Discharge
t <- river_data$Timestamp

# Number of water years
n_years <- n_distinct(river_data$water_year)

# Find where the rise before a peak starts
# Looks back over the previous 7 days (or back to 'earliest' if given) and returns the
# last reading within 5% of the rise above the lowest flow in that window
rise_start_row <- function(p, earliest = t[p] - 7 * 24 * 3600) {
  window <- which(t >= earliest & t <= t[p])
  lowest <- min(q[window])
  max(window[q[window] <= lowest + 0.05 * (q[p] - lowest)])
}

# Months in water year order (October to September)
wy_months <- month.abb[c(10:12, 1:9)]

# Colours for each season, used across all seasonal plots
season_colours <- c("Winter" = "#5298BC", "Spring" = "#BFDC80", "Summer" = "#FEE060",
                    "Autumn" = "#AA1945")

# High flow definition: the 95th percentile of discharge
high_flow <- quantile(river_data$Discharge, 0.95)

# Hydrograph
hydrograph <- ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
  geom_line(colour = "#5298BC", linewidth = 0.4) +
  labs(
    title = "Hydrograph: Don at Parkhill (15-min Interval)",
    subtitle = "1 October 2016 - 30 September 2026",
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()

# Flow Record and Distribution

# Summary statistics for the 15-minute discharge readings
summary_stats <- river_data %>%
  summarise(
    `Mean`               = mean(Discharge),
    `Standard deviation` = sd(Discharge),
    `Minimum`            = min(Discharge),
    `25th percentile`    = quantile(Discharge, 0.25),
    `Median`             = median(Discharge),
    `75th percentile`    = quantile(Discharge, 0.75),
    `95th percentile`    = quantile(Discharge, 0.95),
    `99th percentile`    = quantile(Discharge, 0.99),
    `Maximum`            = max(Discharge)
  ) %>%
  pivot_longer(everything(), names_to = "Statistic", values_to = "Value (m³/s)") %>%
  mutate(`Value (m³/s)` = sprintf("%.2f", `Value (m³/s)`))

# Distrubtion of the data
distribution_plot <- ggplot(river_data, aes(x = Discharge)) +
  geom_histogram(bins = 60, fill = "#5298BC", colour = "white") +
  geom_vline(xintercept = quantile(river_data$Discharge, c(0.5, 0.95, 0.99)),
             linetype = "dashed", colour = "#AA1945") +
  scale_x_log10() +
  labs(
    title = "Distribution of 15-minute discharge",
    subtitle = "Dashed lines: median, 95th and 99th percentiles (log scale)",
    x = expression("Discharge ("*m^3*"/s)"),
    y = "Number of readings"
  ) +
  theme_minimal()

# Seasonality of High Flows

# Median discharge for each month, for the polar plot
monthly <- river_data %>%
  mutate(month_num = as.numeric(format(Timestamp, "%m"))) %>%
  group_by(month_num) %>%
  summarise(value = median(Discharge)) %>%
  mutate(
    month  = month.abb[month_num],
    season = case_when(
      month_num %in% c(9, 10, 11) ~ "Autumn",
      month_num %in% c(12, 1, 2)  ~ "Winter",
      month_num %in% c(3, 4, 5)   ~ "Spring",
      TRUE                        ~ "Summer"
    ),
    season = factor(season, levels = c("Autumn", "Winter", "Spring", "Summer")),
    month_order = (month_num - 10) %% 12   # Oct = 0, Nov = 1, ..., Sep = 11
  ) %>%
  arrange(month_order)
monthly$id <- 1:nrow(monthly)

# Angle and alignment of each month label
monthly$angle <- 90 - 360 * (monthly$id - 0.5) / nrow(monthly)
monthly$hjust <- ifelse(monthly$angle < -90, 1, 0)
monthly$angle <- ifelse(monthly$angle < -90, monthly$angle + 180, monthly$angle)

# Grid lines every 10 up to 40
grid_breaks <- c(10, 20, 30, 40)
y_max <- max(c(monthly$value, 40))

# Days spent above the 95th percentile in each month
days_above <- river_data %>%
  filter(Discharge > high_flow) %>%
  mutate(month = factor(month, levels = wy_months)) %>%
  count(month, season) %>%
  mutate(days = n / 96)   # 96 readings per day

#Plot of the distubtion of mean discharges by month
mean_dis_plot <- ggplot(monthly, aes(x = as.factor(id), y = value, fill = season)) +
  geom_blank() +
  geom_col(alpha = 0.6, colour = "black", linewidth = 0.3) +
  geom_hline(yintercept = grid_breaks, colour = "grey50", linewidth = 0.3, alpha = 0.5) +
  annotate("text", x = 0.5, y = grid_breaks, label = grid_breaks,
           colour = "grey30", size = 3, fontface = "bold", vjust = -0.3) +
  scale_fill_manual(name = "Seasons", values = season_colours) +
  ylim(-y_max, y_max * 1.4) +
  coord_polar() +
  labs(title = "Median discharge by month") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  geom_text(aes(x = id, y = value + y_max * 0.05, label = month,
                hjust = hjust, angle = angle),
            colour = "black", alpha = 0.7, size = 3.5, fontface = "bold",
            inherit.aes = FALSE)

# Plot showing time above the 95th percentile by month
time_above_95_plot <- ggplot(days_above, aes(x = month, y = days, fill = season)) +
  geom_col(alpha = 0.6, colour = "black") +
  scale_fill_manual(values = season_colours) +
  labs(
    title = "Time spent above the 95th percentile, by month",
    x = "Month",
    y = "Total days above the 95th percentile"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


# Anatomy of a Flood Event

# The largest recorded peak
flood_peak <- river_data %>% slice_max(Discharge, n = 1, with_ties = FALSE)
peak_time  <- flood_peak$Timestamp
peak_flow  <- flood_peak$Discharge

# Start of the rise: the last reading in the week before the peak that is within
# 5% of that week's lowest flow (the same definition as in the independence analysis)
rise_start <- river_data %>%
  filter(Timestamp >= peak_time - 7 * 24 * 3600, Timestamp < peak_time) %>%
  filter(Discharge <= 1.05 * min(Discharge)) %>%
  slice_max(Timestamp, n = 1, with_ties = FALSE)
start_time <- rise_start$Timestamp

# Time to peak in hours
time_to_peak <- as.numeric(difftime(peak_time, start_time, units = "hours"))

# Readings from two days before the rise to five days after the peak
flood_event <- river_data %>%
  filter(Timestamp >= start_time - 2 * 24 * 3600,
         Timestamp <= peak_time + 5 * 24 * 3600)

#Plot showing the anatomy of a flood year
flood_anatomy_plot <- ggplot(flood_event, aes(x = Timestamp, y = Discharge)) +
  geom_line(colour = "#5298BC", linewidth = 0.6) +
  geom_vline(xintercept = c(start_time, peak_time), linetype = "dashed", colour = "gray60") +
  geom_point(data = flood_peak, colour = "#AA1945", size = 2.5) +
  annotate("segment", x = start_time, xend = peak_time,
           y = peak_flow * 1.1, yend = peak_flow * 1.1,
           arrow = arrow(ends = "both", length = unit(0.15, "cm"))) +
  annotate("text", x = start_time + (peak_time - start_time) / 2, y = peak_flow * 1.16,
           label = paste("Time to peak:", round(time_to_peak), "hours"), size = 3.5) +
  annotate("text", x = peak_time, y = peak_flow, hjust = -0.1,
           label = paste("Max flow:", round(peak_flow, 1), "m³/s"), size = 3.5) +
  labs(
    title = "Anatomy of a flood event",
    subtitle = format(peak_time, "Peak on %d %B %Y"),
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()

# Variation in Flood Hydrograph Shape

# The largest peak in each water year
annual_peaks <- river_data %>%
  group_by(water_year) %>%
  slice_max(Discharge, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(water_year, peak_time = Timestamp, peak_flow = Discharge)

# Readings from five days before to five days after each annual peak, with time
# measured from the peak and discharge scaled by the peak flow
annual_peak_shapes <- map_dfr(1:nrow(annual_peaks), function(i) {
  peak <- annual_peaks[i, ]
  river_data %>%
    filter(Timestamp >= peak$peak_time - 5 * 24 * 3600,
           Timestamp <= peak$peak_time + 5 * 24 * 3600) %>%
    mutate(
      water_year      = peak$water_year,
      hours_from_peak = as.numeric(difftime(Timestamp, peak$peak_time, units = "hours")),
      scaled_flow     = Discharge / peak$peak_flow
    )
})

# Plot visualising largest peaks by water year
year_colours <- c(
  "2017" = "#6BBED1",
  "2018" = "#347F89",
  "2019" = "#82B7E0",
  "2020" = "#4285C5",
  "2021" = "#5298BC",
  "2022" = "#4169E1",
  "2023" = "#2764A5",
  "2024" = "#2445A2",
  "2025" = "#174B67",
  "2026" = "#142D59"
)
peak_wateryear_plot <- ggplot(annual_peak_shapes, aes(x = hours_from_peak, y = scaled_flow,
                                                      colour = factor(water_year), group = water_year)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#AA1945") +
  scale_colour_manual(values = year_colours) +
  labs(
    title = "Largest peak in each water year",
    subtitle = "Aligned at the peak and scaled by peak flow",
    x = "Hours from peak",
    y = "Discharge as a proportion of peak flow",
    colour = "Water year"
  ) +
  guides(colour = guide_legend(
    override.aes = list(linewidth = 1.5)
  ))+
  theme_minimal()
#-----------------------------------------------------------------------------------
# Section 3 - Threshold Selection and Independence
#-----------------------------------------------------------------------------------

# Candidate peaks: the highest reading in each separate period above threshold u,
# returned as row numbers of river_data
candidate_peaks <- function(u) {
  river_data %>%
    mutate(row = row_number(),
           above = Discharge > u,
           run = cumsum(above & !lag(above, default = FALSE))) %>%
    filter(above) %>%
    group_by(run) %>%
    slice_max(Discharge, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull(row)
}


# Subsection 3.1 - Initial Threshold
# --------------------------------------------------------------

# Candidate peaks per year at every whole-number threshold, before independence
crossing_thresholds <- 0:floor(max(q))
crossings_per_year <- sapply(crossing_thresholds, function(u) {
  above <- q > u
  sum(above & !lag(above, default = FALSE))
}) / n_years

# Highest threshold giving at least a target number of candidate peaks per year
crossing_threshold <- function(target) {
  max(crossing_thresholds[crossings_per_year >= target])
}

# The four candidate thresholds
threshold_lines <- tibble(
  Method = c("95th percentile", "5 peaks/year", "4 peaks/year", "3 peaks/year"),
  Value  = c(as.numeric(quantile(q, 0.95)),
             crossing_threshold(5), crossing_threshold(4), crossing_threshold(3))
)

# Initial threshold: the 95th percentile, the lowest of the four
u_initial <- crossing_threshold(5)

# Hydrograph with the four candidate thresholds
threshold_plot <- ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
  geom_line(colour = "#5298BC", linewidth = 0.35) +
  geom_hline(data = threshold_lines, aes(yintercept = Value, colour = Method),
             linetype = "dashed", linewidth = 0.7) +
  geom_text(data = threshold_lines,
            aes(x = max(river_data$Timestamp), y = Value, label = round(Value), colour = Method),
            hjust = 1, vjust = -0.4, size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = c("95th percentile" = "#7FA83A",
                                 "5 peaks/year"    = "#E58FA6",
                                 "4 peaks/year"    = "#D0436A",
                                 "3 peaks/year"    = "#AA1945")) +
  labs(
    title = "Candidate thresholds",
    subtitle = "Peaks per year counted before independence",
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)"),
    colour = "Threshold"
  ) +
  theme_minimal()


# Subsection 3.2 - Independence of Peaks
# --------------------------------------------------------------

# Apply the Robson and Reed (1999) independence rules to the candidate peaks above u
independent_peaks <- function(u) {
  
  candidates <- candidate_peaks(u)
  
  # Clean peaks: no other candidate peak in the previous 7 days
  days_since_last <- as.numeric(difftime(t[candidates], lag(t[candidates]), units = "days"))
  clean <- candidates[is.na(days_since_last) | days_since_last > 7]
  if (length(clean) < 5) stop("Fewer than five clean peaks above ", u)
  
  # Average time to rise over the clean peaks
  rise_hours <- sapply(clean, function(p) {
    as.numeric(difftime(t[p], t[rise_start_row(p)], units = "hours"))
  })
  min_gap_hours <- 3 * mean(rise_hours)
  
  # Check each candidate against the last peak kept
  kept <- candidates[1]
  for (p in candidates[-1]) {
    last <- kept[length(kept)]                                  # the first peak of the pair
    gap <- as.numeric(difftime(t[p], t[last], units = "hours"))
    trough <- min(q[last:p])                                      # flow at the trough between them
    
    # Rule 1: at least three times the average time to rise apart
    # Rule 2: trough at least two-thirds lower than the first peak
    if (gap >= min_gap_hours & trough <= q[last] / 3) {
      kept <- c(kept, p)                 # independent: a new event
    } else if (q[p] > q[last]) {
      kept[length(kept)] <- p            # same event: keep the larger peak
    }
  }
  
  # Return the independent peaks, recording the average time to rise
  peaks_u <- river_data[kept, ]
  attr(peaks_u, "mean_rise_hours") <- mean(rise_hours)
  peaks_u
}

# Candidate and independent peaks at the initial threshold
candidates_initial <- river_data[candidate_peaks(u_initial), ]
peaks_initial <- independent_peaks(u_initial)
mean_rise_hours <- attr(peaks_initial, "mean_rise_hours")
min_gap_hours <- 3 * mean_rise_hours

# Hydrograph with a set of peaks marked
plot_peaks <- function(peak_data, u, title) {
  ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
    geom_line(colour = "#5298BC", linewidth = 0.35) +
    geom_hline(yintercept = u, linetype = "dashed", colour = "#7FA83A") +
    geom_point(data = peak_data, colour = "#AA1945", size = 1.2) +
    labs(
      title = title,
      subtitle = paste0(nrow(peak_data), " peaks above ", round(u, 1), " m³/s"),
      x = "Date",
      y = expression("Discharge ("*m^3*"/s)")
    ) +
    theme_minimal()
}

# Before and after independence plots, shown side by side in the report
plot_before <- plot_peaks(candidates_initial, u_initial, "Before independence")
plot_after <- plot_peaks(peaks_initial, u_initial, "After independence")


# Subsection 3.3 - Final Threshold
# --------------------------------------------------------------

# Independent peaks per year across a range of thresholds
independent_rates <- tibble(threshold = 50:120) %>%
  mutate(peaks_per_year = sapply(threshold, function(u) {
    peaks_u <- tryCatch(independent_peaks(u), error = function(e) NULL)
    if (is.null(peaks_u)) NA_real_ else nrow(peaks_u) / n_years
  })) %>%
  filter(!is.na(peaks_per_year))

# Thresholds giving 3 to 5 independent peaks per year
thresholds_3to5 <- independent_rates %>%
  filter(peaks_per_year >= 3, peaks_per_year <= 5) %>%
  pull(threshold)

# Final threshold, chosen from the range above using the GPD checks
u_final <- 65

# Independent peaks per year against threshold, with the 3-5 band and final threshold
rates_plot <- ggplot(independent_rates, aes(x = threshold, y = peaks_per_year)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3, ymax = 5,
           fill = "#BFDC80", alpha = 0.3) +
  geom_line(colour = "#5298BC") +
  geom_point(colour = "#5298BC", size = 1) +
  geom_vline(xintercept = u_final, linetype = "dashed", colour = "#AA1945") +
  labs(
    title = "Independent peaks per year by threshold",
    subtitle = "Shaded band: 3 to 5 peaks per year | Dashed line: chosen threshold",
    x = expression("Threshold ("*m^3*"/s)"),
    y = "Independent peaks per year"
  ) +
  theme_minimal()

# Independent peaks at the final threshold
peaks <- independent_peaks(u_final)

# GPD fitted to the excesses of the independent peak flows above the final threshold
gpd_fit <- gpd.fit(peaks$Discharge, u_final, npy = nrow(peaks) / n_years, show = FALSE)

# Start of each rise: searching back up to 7 days, but never past the previous peak
peak_rows <- match(peaks$Timestamp, t)
previous_rows <- c(1, head(peak_rows, -1))
start_rows <- mapply(function(p, prev) {
  rise_start_row(p, earliest = max(t[p] - 7 * 24 * 3600, t[prev]))
}, peak_rows, previous_rows)

# The two quantities SEPA asked for
peaks <- peaks %>%
  mutate(
    start_time = t[start_rows],
    max_flow = Discharge,
    time_to_peak_h = as.numeric(difftime(Timestamp, start_time, units = "hours"))
  )

#---------------------------------------------------------------------------
# Section 4 - Uncertainty Analysis
#--------------------------------------------------------------------------

# Max flow
Max_flow_plot <- ggplot(peaks, aes(x = max_flow)) +
  geom_histogram(bins = 15, fill = "#5298BC", colour = "#2A5267") +
  labs(
    title = "Distribution of max flow",
    x = expression("Max flow ("*m^3*"/s)"),
    y = "Number of events"
  ) +
  theme_minimal()

# Time to peak
Time_to_peak_plot <- ggplot(peaks, aes(x = time_to_peak_h)) +
  geom_histogram(bins = 20, fill = "#5298BC", colour = "#2A5267") +
  labs(
    title = "Distribution of time to peak",
    x = "Time to peak (hours)",
    y = "Number of events"
  ) +
  theme_minimal()

Max_vs_time <- ggplot(peaks, aes(x = time_to_peak_h, y = max_flow)) +
  geom_point(color = "#5298BC") +
  geom_smooth(method = "lm", color = "#AA1945", fill = "#FBB1BA", alpha = 0.5) +
  labs(
    title = "Max flow against time to peak",
    x = "Time to peak (hours)",
    y = expression("Max flow ("*m^3*"/s)")
  ) +
  theme_minimal()

cor_result <- cor.test(peaks$time_to_peak_h, peaks$max_flow, method = "spearman", exact = FALSE)


# Correlation table
cor_table <- tibble(
  Test = "Spearman's rank correlation",
  Rho = unname(cor_result$estimate),
  `P-value` = cor_result$p.value,
  `S statistic` = unname(cor_result$statistic)
)

cor_result_table <- knitr::kable(
  cor_table,
  digits = 3,
  booktabs = TRUE,
  caption = "Spearman's Rank Correlation Between Peak Discharge and Time to Peak"
)


# Analysis of uncertanity
# Number of events and typical values in each season
events_season_table <- peaks %>%
  group_by(season) %>%
  summarise(
    Events = n(),
    `Median max flow` = median(max_flow),
    `Median time to peak (h)` = median(time_to_peak_h)
  )

events_season_table <- knitr::kable(
  events_season_table,
  booktabs = TRUE,
  digits = 3,
  caption = "Summary of Max Flow and Time to Peak for each Season"
)

# Relationship, coloured by season
Max_vs_time_season <- ggplot(peaks, aes(x = time_to_peak_h, y = max_flow )) +
  geom_point(size = 2, aes(colour = season)) +
  scale_color_manual(values = c("Autumn" = "#AA1945", "Winter" = "#5298BC",
                                "Spring" = "#BFDC80", "Summer" = "#FEE060")) +
  labs(
    title = "Max flow against time to peak, by season",
    x = "Time to peak (hours)",
    y = expression("Max flow ("*m^3*"/s)"),
    color = "Season"
  ) +
  theme_minimal()

# Max flow by season
max_season_plot <- ggplot(peaks, aes(x = season, y = max_flow)) +
  geom_boxplot(fill = season_colours, alpha = 0.7) +
  labs(
    title = "Max flow by season",
    x = "Season",
    y = expression("Max flow ("*m^3*"/s)")
  ) +
  theme_minimal()

# Time to peak by season
time_season_plot <- peaks %>% group_by(season) %>% ggplot(aes(x = season, y = time_to_peak_h)) +
  geom_boxplot(fill = season_colours, alpha = 0.7) +
  labs(
    title = "Time to peak by season",
    x = "Season",
    y = "Time to peak (hours)"
  ) +
  theme_minimal()

# Uncertanity final table
results <- tibble()

for (u in c(55, 60, 65, 70, 75)) {
  peaks_test <- independent_peaks(u)
  
  # Start of each rise and time to peak (same code as before)
  peak_rows <- match(peaks_test$Timestamp, t)
  previous  <- c(1, head(peak_rows, -1))
  start_rows <- mapply(function(p, prev) {
    window <- which(t >= max(t[p] - 7 * 24 * 3600, t[prev]) & t <= t[p])
    low    <- min(q[window])
    max(window[q[window] <= low + 0.05 * (q[p] - low)])
  }, peak_rows, previous)
  time_to_peak_h <- as.numeric(difftime(peaks_test$Timestamp, t[start_rows], units = "hours"))
  
  # One row of results for this threshold
  results <- bind_rows(results, tibble(
    threshold             = u,
    events                = nrow(peaks_test),
    events_per_year       = nrow(peaks_test) / n_distinct(river_data$water_year),
    median_max_flow       = median(peaks_test$Discharge),
    median_time_to_peak_h = median(time_to_peak_h),
    spearman_rho          = cor(peaks_test$Discharge, time_to_peak_h, method = "spearman")
  ))
}


results_table <- knitr::kable(
  results,
  col.names = c("threshold", "events", "events_per_year" = "Events per Year",
                "median_max_flow" = "Median Max Flow",
                "median_time_to_peak_h" = "Median Time to Peak (Hours)", 
                "spearman_rho" = "Spearman Rho Corelation"),
  booktabs = T,
  digits = 3,
  caption = "Final Uncertanity Table"
) %>%
  kable_styling(latex_options = "striped")