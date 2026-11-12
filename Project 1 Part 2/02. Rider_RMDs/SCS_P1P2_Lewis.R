# ==============================================================================
# Statistical Case Studies - Project 1 (Part 2): Data preparation
# Authors: Anna Scollay, Lewis Seath, Zacharias Onyejiaka
# ==============================================================================
#
# Builds the standardised hydrograph dataset that every shape method uses:
#
#   1. Data cleaned exactly as in Part 1: only quality codes 50 and 100 kept,
#      with every exclusion recorded (Lecture 1: record any deletion).
#   2. Independent peaks extracted exactly as in Part 1 (u = 65 m^3/s,
#      Robson & Reed 1999 independence criteria; Lecture 2).
#   3. Each event runs from the trough before the peak to the trough after it.
#      The start uses Part 1's definition, so the rising part of each event
#      equals the time to peak modelled in Part 1.
#   4. Discharge is divided by the event's peak flow, separating shape from
#      magnitude (consistent with SEPA's current peak-scaling approach).
#   5. Time is rescaled onto [0, 1] by the event width (Lecture 4), and the
#      width is kept as a separate parameter.
#
# Main outputs:
#   events     : one row per event (peak flow, width, time to peak, etc.)
#   Y          : N x n_grid matrix of standardised curves, one row per event
#   mean_curve : mean standardised hydrograph
#   Y_centred  : Y minus the mean curve (input for PCA)
#   gap_check  : unchecked readings vs linear interpolation (justifies filling)
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages and settings
# ------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)

cache_file <- "river_don_parkhill_flow.rds"
options(timeout = 300)   # SEPA API can be slow if the cache is absent

u_final       <- 65    # POT threshold from Part 1 (m^3/s)
search_days   <- 7     # Part 1 search span for the trough either side of a peak
trough_factor <- 1.05  # trough = last/first reading within 5% of the local minimum
n_grid        <- 64    # points on the common [0, 1] grid


# ------------------------------------------------------------------------------
# 2. Load raw data (cached locally, otherwise fetched from the SEPA API)
# ------------------------------------------------------------------------------
if (file.exists(cache_file)) {
  river_raw <- readRDS(cache_file)
} else {
  fetch_chunk <- function(yr) {
    url <- paste0(
      "https://timeseries.sepa.org.uk/KiWIS/KiWIS?",
      "service=kisters&type=queryServices&datasource=0&",
      "request=getTimeseriesValues&ts_path=1/234292/Q/15m.Cmd&",
      "from=", yr, "-10-01T00:00:00&to=", yr + 1, "-09-30T23:59:59&",
      "metadata=true&returnfields=Timestamp,Value,Quality%20Code&format=csv"
    )
    read_delim(url, delim = ";", comment = "#",
               col_names = c("Timestamp", "Discharge", "QualityCode"),
               show_col_types = FALSE)
  }
  
  river_raw <- map_dfr(2016:2025, fetch_chunk) %>%
    mutate(Timestamp = as.POSIXct(Timestamp, format = "%Y-%m-%dT%H:%M:%S")) %>%
    filter(!is.na(Discharge)) %>%
    distinct(Timestamp, .keep_all = TRUE)
  
  saveRDS(river_raw, cache_file)
}

# as.character() first, so the true code is kept even if stored as a factor
river_raw <- river_raw %>%
  mutate(QualityCode = as.numeric(as.character(QualityCode))) %>%
  arrange(Timestamp)


# ------------------------------------------------------------------------------
# 3. Clean as in Part 1: keep quality codes 50 and 100 only
# ------------------------------------------------------------------------------
river_data <- river_raw %>%
  filter(QualityCode %in% c(50, 100))

# Record of what was excluded, by quality code (for the data section)
exclusion_summary <- river_raw %>%
  filter(!QualityCode %in% c(50, 100)) %>%
  count(QualityCode, name = "readings_excluded") %>%
  mutate(percent_of_all = 100 * readings_excluded / nrow(river_raw))

print(exclusion_summary)

# Vectors used by the peak-extraction functions
q     <- river_data$Discharge
times <- river_data$Timestamp


# ------------------------------------------------------------------------------
# 4. Independent peaks, as in Part 1 (Robson & Reed 1999)
# ------------------------------------------------------------------------------

# Start of the rising limb (Part 1 definition): the last reading before the
# peak that lies within 5% of the minimum flow over the preceding 7 days
rise_start_row <- function(p, earliest = times[p] - search_days * 24 * 3600) {
  window <- which(times >= earliest & times <= times[p])
  lowest <- min(q[window])
  max(window[q[window] <= trough_factor * lowest])
}

# Candidate peaks: the maximum of each run above u lasting at least 1 hour
candidate_peaks <- function(u, min_duration_hours = 1) {
  min_readings <- min_duration_hours * 4
  river_data %>%
    mutate(
      row   = row_number(),
      above = Discharge > u,
      run   = cumsum(above & !lag(above, default = FALSE))
    ) %>%
    filter(above) %>%
    group_by(run) %>%
    filter(n() >= min_readings) %>%
    slice_max(Discharge, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull(row)
}

# Independence: peaks at least 3 x mean time to rise apart, and the trough
# between them at least two-thirds below the first peak
independent_peaks <- function(u) {
  candidates <- candidate_peaks(u)
  
  # "Clean" peaks (largest within +/- 7 days) used to estimate mean rise time
  clean <- candidates[sapply(candidates, function(p) {
    near <- abs(as.numeric(difftime(times[candidates], times[p], units = "days"))) <= 7
    q[p] == max(q[candidates[near]])
  })]
  
  rise_hours <- sapply(clean, function(p) {
    as.numeric(difftime(times[p], times[rise_start_row(p)], units = "hours"))
  })
  min_gap_hours <- 3 * mean(rise_hours)
  
  kept <- candidates[1]
  for (p in candidates[-1]) {
    last   <- kept[length(kept)]
    gap    <- as.numeric(difftime(times[p], times[last], units = "hours"))
    trough <- min(q[last:p])
    
    if (gap >= min_gap_hours & trough <= q[last] / 3) {
      kept <- c(kept, p)
    } else if (q[p] > q[last]) {
      kept[length(kept)] <- p
    }
  }
  kept
}

peak_rows <- independent_peaks(u_final)
cat("Number of independent events:", length(peak_rows), "\n")


# ------------------------------------------------------------------------------
# 5. Event start, end and width (trough to trough)
# ------------------------------------------------------------------------------

# End of the recession, mirroring the Part 1 start rule: the first reading
# after the peak within 5% of the minimum flow over the following 7 days.
# The search stops at the next event's peak, so the end can never fall
# inside the next flood.
recession_end_row <- function(p, latest) {
  window <- which(times >= times[p] & times <= latest)
  lowest <- min(q[window])
  min(window[q[window] <= trough_factor * lowest])
}

search_limit <- pmin(
  times[peak_rows] + search_days * 24 * 3600,
  lead(times[peak_rows], default = max(times))
)

events <- tibble(
  event_id  = seq_along(peak_rows),
  peak_row  = peak_rows,
  start_row = map_int(peak_rows, rise_start_row),
  end_row   = map2_int(peak_rows, search_limit, recession_end_row)
) %>%
  mutate(
    peak_time      = times[peak_row],
    start_time     = times[start_row],
    end_time       = times[end_row],
    peak_flow      = q[peak_row],
    time_to_peak_h = as.numeric(difftime(peak_time, start_time, units = "hours")),
    recession_h    = as.numeric(difftime(end_time, peak_time, units = "hours")),
    width_h        = time_to_peak_h + recession_h,
    peak_position  = time_to_peak_h / width_h   # where the peak falls in [0, 1]
  )

print(events %>% select(event_id, peak_time, peak_flow, time_to_peak_h,
                        recession_h, width_h, peak_position))


# ------------------------------------------------------------------------------
# 6. Gap check: are the excluded unchecked readings close to interpolation?
# ------------------------------------------------------------------------------
# Excluding unchecked readings leaves short gaps inside some events. These are
# filled by linear interpolation when curves are put on the common grid. This
# check compares each excluded unchecked reading (codes 200/254) inside an event
# with the value linear interpolation gives at that time, to justify the filling.

gap_check <- events %>%
  select(event_id, start_time, end_time) %>%
  pmap_dfr(function(event_id, start_time, end_time) {
    removed <- river_raw %>%
      filter(Timestamp >= start_time, Timestamp <= end_time,
             QualityCode %in% c(200, 254))
    if (nrow(removed) == 0) return(NULL)
    
    kept <- river_data %>%
      filter(Timestamp >= start_time, Timestamp <= end_time)
    
    removed %>%
      mutate(
        event_id     = event_id,
        interpolated = approx(as.numeric(kept$Timestamp), kept$Discharge,
                              xout = as.numeric(Timestamp))$y,
        abs_diff     = abs(Discharge - interpolated),
        pct_diff     = 100 * abs_diff / interpolated
      )
  })

if (nrow(gap_check) > 0) {
  gap_check %>%
    group_by(event_id) %>%
    summarise(
      readings_filled = n(),
      max_abs_diff    = max(abs_diff, na.rm = TRUE),
      max_pct_diff    = max(pct_diff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    print()
} else {
  cat("No unchecked readings fall inside any event.\n")
}


# ------------------------------------------------------------------------------
# 7. Standardise: divide by peak flow, rescale time onto [0, 1]
# ------------------------------------------------------------------------------
# For event i with start s_i and width w_i:
#   u = (t - s_i) / w_i  in [0, 1],   y(u) = Q(t) / Q_peak,i
# Each curve is then linearly interpolated onto the same grid of n_grid
# points, which also fills the gaps checked in section 6.

grid_u <- seq(0, 1, length.out = n_grid)

standardise_event <- function(start_row, end_row, start_time, width_h, peak_flow) {
  idx <- start_row:end_row
  u   <- as.numeric(difftime(times[idx], start_time, units = "hours")) / width_h
  approx(u, q[idx] / peak_flow, xout = grid_u)$y
}

Y <- events %>%
  select(start_row, end_row, start_time, width_h, peak_flow) %>%
  pmap(standardise_event) %>%
  do.call(rbind, .)

rownames(Y) <- paste0("event_", events$event_id)
colnames(Y) <- paste0("u_", sprintf("%.3f", grid_u))

stopifnot(!anyNA(Y))

mean_curve <- colMeans(Y)
Y_centred  <- sweep(Y, 2, mean_curve)

cat("Standardised data:", nrow(Y), "events x", ncol(Y), "grid points\n")


# ------------------------------------------------------------------------------
# 8. Exploratory plots
# ------------------------------------------------------------------------------

# 8a. Descriptive: events in real time and m^3/s, aligned at the peak
#     (48 h before to 72 h after, for the data section of the report)
real_time_long <- events %>%
  select(event_id, peak_time) %>%
  pmap_dfr(function(event_id, peak_time) {
    river_data %>%
      filter(Timestamp >= peak_time - hours(48),
             Timestamp <= peak_time + hours(72)) %>%
      transmute(
        event_id,
        hours_from_peak = as.numeric(difftime(Timestamp, peak_time, units = "hours")),
        Discharge
      )
  })

plot_real_time <- ggplot(real_time_long,
                         aes(hours_from_peak, Discharge, group = event_id)) +
  geom_line(colour = "#5298BC", linewidth = 0.4, alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#AA1945") +
  geom_hline(yintercept = u_final, linetype = "dotted") +
  scale_x_continuous(breaks = seq(-48, 72, by = 12)) +
  labs(
    title = sprintf("Extreme flow events on the Don at Parkhill (N = %d)",
                    nrow(events)),
    subtitle = sprintf("Aligned at the peak; dotted line: threshold u = %d m³/s",
                       u_final),
    x = "Hours from peak",
    y = expression("Discharge (" * m^3 * "/s)")
  ) +
  theme_minimal()

# 8b. Standardised curves with mean and 10th-90th percentile band
standardised_long <- as_tibble(Y) %>%
  mutate(event_id = events$event_id) %>%
  pivot_longer(-event_id, values_to = "scaled_flow") %>%
  mutate(u = rep(grid_u, times = nrow(Y)))

band <- standardised_long %>%
  group_by(u) %>%
  summarise(
    mean_flow = mean(scaled_flow),
    lower     = quantile(scaled_flow, 0.10),
    upper     = quantile(scaled_flow, 0.90),
    .groups   = "drop"
  )

plot_standardised <- ggplot() +
  geom_line(data = standardised_long,
            aes(u, scaled_flow, group = event_id),
            colour = "#5298BC", linewidth = 0.4, alpha = 0.4) +
  geom_ribbon(data = band, aes(u, ymin = lower, ymax = upper),
              fill = "#5298BC", alpha = 0.25) +
  geom_line(data = band, aes(u, mean_flow),
            colour = "#174B67", linewidth = 1.1) +
  labs(
    title = "Standardised hydrographs",
    subtitle = "Bold: mean curve; band: 10th-90th percentile",
    x = "Rescaled time (0 = start of event, 1 = end of event)",
    y = "Discharge / peak discharge"
  ) +
  theme_minimal()

# 8c. Event width against peak flow (do larger floods last longer?)
plot_width_peak <- ggplot(events, aes(peak_flow, width_h)) +
  geom_point(colour = "#174B67", size = 2) +
  labs(
    title = "Event width against peak flow",
    x = expression("Peak flow (" * m^3 * "/s)"),
    y = "Event width (hours)"
  ) +
  theme_minimal()

print(plot_real_time)
print(plot_standardised)
print(plot_width_peak)

cat("Rank correlation, width vs peak flow:",
    round(cor(events$width_h, events$peak_flow, method = "spearman"), 3), "\n")
cat("Peak position in [0, 1]: range",
    paste(round(range(events$peak_position), 3), collapse = " to "), "\n")

# ==============================================================================
# SHARED EVALUATION FUNCTIONS
# ==============================================================================
# Every method is checked with these same functions, so that the methods are
# compared fairly (same conversion, same metrics, same unseen events).
#
#   9.  to_real()      : standardised curve -> hours and m^3/s
#   10. fit_metrics()  : errors for one event; evaluate_fits() for all events
#   11. loocv()        : leave-one-out cross-validation for any method
#   12. test_ids       : fixed held-out events for testing the uncertainty model
#   13. Baseline       : SEPA's current approach (mean shape), as an example
# ==============================================================================


# ------------------------------------------------------------------------------
# 9. Reconstruction: back to real time and discharge
# ------------------------------------------------------------------------------
# A standardised curve y (on grid_u) for an event with width w and peak flow
# Q_peak becomes:  hours since start = u * w,   discharge = y * Q_peak

to_real <- function(y, width_h, peak_flow) {
  data.frame(
    hours     = grid_u * width_h,
    discharge = y * peak_flow
  )
}


# ------------------------------------------------------------------------------
# 10. Metrics: how close is a fitted curve to the observed one?
# ------------------------------------------------------------------------------
# Both curves are converted to real units first, so errors are in m^3/s and
# hours. The event's own width and peak flow are used for both curves, so
# only the error in SHAPE is measured.
#
#   rmse         : overall closeness of the curve (m^3/s)
#   peak_error   : fitted max minus observed max (m^3/s)
#   timing_error : fitted time of max minus observed time of max (hours)
#   volume_error : % difference in total volume (affects flood duration)
#
# Peak, timing and volume errors keep their sign (+ = too high / too late /
# too much water); summarise_metrics() averages their absolute values.

fit_metrics <- function(y_obs, y_fit, width_h, peak_flow) {
  obs <- to_real(y_obs, width_h, peak_flow)
  fit <- to_real(y_fit, width_h, peak_flow)
  
  c(
    rmse         = sqrt(mean((fit$discharge - obs$discharge)^2)),
    peak_error   = max(fit$discharge) - max(obs$discharge),
    timing_error = fit$hours[which.max(fit$discharge)] -
      obs$hours[which.max(obs$discharge)],
    volume_error = 100 * (sum(fit$discharge) - sum(obs$discharge)) /
      sum(obs$discharge)
  )
}

# Metrics for every event. Y_fit has the same layout as Y (one row per event).
evaluate_fits <- function(Y_fit) {
  results <- t(sapply(seq_len(nrow(Y)), function(i) {
    fit_metrics(Y[i, ], Y_fit[i, ], events$width_h[i], events$peak_flow[i])
  }))
  data.frame(event_id = events$event_id, results)
}

# One-line summary of a method: average absolute error for each metric
summarise_metrics <- function(results) {
  colMeans(abs(results[, c("rmse", "peak_error", "timing_error", "volume_error")]))
}


# ------------------------------------------------------------------------------
# 11. Leave-one-out cross-validation (Lecture 5)
# ------------------------------------------------------------------------------
# Each method supplies ONE function of its own:
#
#   fit_and_predict(Y_train, y_new)
#     - fits the method using only the curves in Y_train
#     - returns its approximation of the new curve y_new (length 64)
#
# loocv() leaves out each event in turn, fits on the other events, predicts
# the left-out event, and returns all predictions in the same layout as Y.
# Pass the result to evaluate_fits() to get out-of-sample errors.

loocv <- function(fit_and_predict) {
  Y_pred <- Y
  for (i in seq_len(nrow(Y))) {
    Y_pred[i, ] <- fit_and_predict(Y_train = Y[-i, , drop = FALSE],
                                   y_new   = Y[i, ])
  }
  Y_pred
}


# ------------------------------------------------------------------------------
# 12. Fixed held-out events (for testing the uncertainty model in Phase 3)
# ------------------------------------------------------------------------------
# About 20% of events, chosen once with a fixed seed so everyone uses the same.

set.seed(2026)
test_ids  <- sort(sample(events$event_id, size = round(0.2 * nrow(events))))
train_ids <- setdiff(events$event_id, test_ids)

cat("Held-out events:", test_ids, "\n")


# ------------------------------------------------------------------------------
# 13. Baseline: SEPA's current approach (one mean shape for every event)
# ------------------------------------------------------------------------------
# Every method's errors should be compared against this. It also shows how
# to use the functions above.

baseline_fit_and_predict <- function(Y_train, y_new) {
  colMeans(Y_train)   # ignores y_new: every event gets the same shape
}

# In-sample: every event gets the mean of all events
Y_baseline <- matrix(mean_curve, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
baseline_in_sample <- evaluate_fits(Y_baseline)

# Out-of-sample: every event gets the mean of the OTHER events
baseline_loocv <- evaluate_fits(loocv(baseline_fit_and_predict))

cat("\nBaseline (mean shape), average absolute errors:\n")
print(rbind(
  in_sample  = summarise_metrics(baseline_in_sample),
  out_sample = summarise_metrics(baseline_loocv)
))

# ==============================================================================
# METHOD 1: POLYNOMIAL FITTING (Lecture 4, slides 9-12)
# ==============================================================================
# Each standardised curve is approximated by a degree-p polynomial in time:
#
#   y_i(u) = theta_i0 + theta_i1 u + theta_i2 u^2 + ... + theta_ip u^p + error
#
# fitted by least squares separately for each event, using lm() with the
# polynomial terms as explanatory variables and the standardised discharge
# as the response (as in Lecture 4 and the linear models handout). The
# intercept is kept: the handout notes that beta_0 = 0 "will not be true for
# your data". The estimated coefficients theta_i are each event's shape
# parameters.
#
# Raw powers u, u^2, ... are used, as in the lecture design matrix. These are
# strongly correlated with each other, so degrees are capped at 10 here, and
# the correlation between fitted coefficients is checked in section 16
# (Lecture 5, slide 12).
# ==============================================================================


# ------------------------------------------------------------------------------
# 14. Choosing the degree p (Lecture 5, slide 12)
# ------------------------------------------------------------------------------

# Build the formula y ~ I(u^1) + I(u^2) + ... + I(u^p), as in the handout
poly_formula <- function(p) {
  as.formula(paste("y ~", paste0("I(u^", 1:p, ")", collapse = " + ")))
}

fit_poly_model <- function(y, p) {
  lm(poly_formula(p), data = data.frame(y = y, u = grid_u))
}

degrees <- 2:10

# 14a. Statistical significance of the highest-order term.
#      As in the handout (quadratic vs cubic), a degree is supported if its
#      highest-order term is significant. For each degree, report the
#      proportion of events where that term has p < 0.05.
#      Caution (Lecture 5, slides 16-17): residuals within a curve are
#      autocorrelated, so standard errors are too small and p-values too
#      optimistic. Use this alongside 14b and 14c, not on its own.
top_term_significant <- sapply(degrees, function(p) {
  mean(sapply(seq_len(nrow(Y)), function(i) {
    coefs <- summary(fit_poly_model(Y[i, ], p))$coefficients
    coefs[nrow(coefs), "Pr(>|t|)"] < 0.05
  }))
})

# 14b. Prediction error, checked separately on training and test events.
#      If the degree chosen on the training events also works well on the
#      held-out test events, the choice generalises (Lecture 5, slide 12).
#      Errors are reported for the curve overall and for the features that
#      matter for flooding: peak, timing and volume.
degree_errors <- function(p, ids) {
  rows  <- which(events$event_id %in% ids)
  Y_fit <- t(sapply(seq_len(nrow(Y)), function(i) fitted(fit_poly_model(Y[i, ], p))))
  summarise_metrics(evaluate_fits(Y_fit)[rows, ])
}

degree_table <- do.call(rbind, lapply(degrees, function(p) {
  rbind(
    data.frame(degree = p, set = "train", t(degree_errors(p, train_ids))),
    data.frame(degree = p, set = "test",  t(degree_errors(p, test_ids)))
  )
}))

print(data.frame(degree = degrees, prop_top_term_significant = top_term_significant))
print(degree_table, digits = 3)

# 14c. Plot errors against degree for training and test events.
#      Small changes in degree should not change the answers much.
plot_degree <- degree_table %>%
  pivot_longer(c(rmse, peak_error, timing_error, volume_error),
               names_to = "metric", values_to = "error") %>%
  ggplot(aes(degree, error, colour = set)) +
  geom_line() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_continuous(breaks = degrees) +
  scale_colour_manual(values = c(train = "#174B67", test = "#AA1945")) +
  labs(
    title = "Polynomial fit: average absolute error by degree",
    x = "Polynomial degree p",
    y = "Average absolute error",
    colour = "Events"
  ) +
  theme_minimal()

print(plot_degree)

# Note: volume error is always zero. With an intercept and equally spaced
# points, least-squares residuals sum to zero, so each fitted curve has the
# same total volume as the observed curve. This is a property of the method.


# ------------------------------------------------------------------------------
# 15. Fit at the chosen degree
# ------------------------------------------------------------------------------
# UPDATE p_chosen after looking at section 14.

p_chosen <- 6

poly_models <- lapply(seq_len(nrow(Y)), function(i) fit_poly_model(Y[i, ], p_chosen))

Y_poly <- t(sapply(poly_models, fitted))

# Shape parameters: one row per event, p_chosen + 1 coefficients
poly_coefs <- t(sapply(poly_models, coef))
colnames(poly_coefs) <- paste0("theta_", 0:p_chosen)
rownames(poly_coefs) <- rownames(Y)

poly_results <- evaluate_fits(Y_poly)

cat("\nPolynomial (degree", p_chosen, ") vs baseline, average absolute errors:\n")
print(rbind(
  baseline   = summarise_metrics(baseline_in_sample),
  polynomial = summarise_metrics(poly_results)
))

# Physical check: fitted curves should not go below zero
cat("Grid points with negative fitted flow:", sum(Y_poly < 0), "\n")


# ------------------------------------------------------------------------------
# 16. Validation (Lecture 5, slide 11; linear models handout)
# ------------------------------------------------------------------------------

# 16a. Fitted vs observed curves for different types of peak: earliest peak,
#      latest peak, longest event and the worst-fitted event.
example_ids <- unique(c(
  events$event_id[which.min(events$peak_position)],
  events$event_id[which.max(events$peak_position)],
  events$event_id[which.max(events$width_h)],
  poly_results$event_id[which.max(poly_results$rmse)]
))

examples_long <- data.frame(
  event_id = rep(events$event_id, each = ncol(Y)),
  u        = rep(grid_u, times = nrow(Y)),
  observed = as.vector(t(Y)),
  fitted   = as.vector(t(Y_poly))
) %>%
  filter(event_id %in% example_ids) %>%
  pivot_longer(c(observed, fitted), names_to = "curve", values_to = "scaled_flow")

plot_examples <- ggplot(examples_long, aes(u, scaled_flow, colour = curve)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~event_id, labeller = label_both) +
  scale_colour_manual(values = c(observed = "#174B67", fitted = "#AA1945")) +
  labs(
    title = sprintf("Polynomial (degree %d): example fits", p_chosen),
    x = "Rescaled time u",
    y = "Discharge / peak discharge",
    colour = NULL
  ) +
  theme_minimal()

print(plot_examples)

# 16b. Residuals against time, all events together. Look for systematic
#      structure (e.g. consistent under-prediction on the rising limb), not
#      just small residual size. Bold line: average residual.
residuals_long <- data.frame(
  event_id = rep(events$event_id, each = ncol(Y)),
  u        = rep(grid_u, times = nrow(Y)),
  residual = as.vector(t(Y - Y_poly))
)

plot_residuals_time <- ggplot(residuals_long, aes(u, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(group = event_id), colour = "#5298BC", alpha = 0.3) +
  stat_summary(fun = mean, geom = "line", colour = "#AA1945", linewidth = 1) +
  labs(
    title = sprintf("Polynomial (degree %d): residuals against time", p_chosen),
    x = "Rescaled time u",
    y = "Observed - fitted (standardised)"
  ) +
  theme_minimal()

print(plot_residuals_time)

# 16c. Residuals against each polynomial term (Lecture 5, slide 11; the
#      handout suggests plotting residuals against the predictor variables).
residual_terms_long <- do.call(rbind, lapply(1:p_chosen, function(k) {
  data.frame(term = paste0("u^", k), value = residuals_long$u^k,
             residual = residuals_long$residual)
}))

plot_residuals_terms <- ggplot(residual_terms_long, aes(value, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(colour = "#5298BC", alpha = 0.2, size = 0.6) +
  facet_wrap(~term, scales = "free_x") +
  labs(
    title = "Residuals against each polynomial term",
    x = "Value of term",
    y = "Residual"
  ) +
  theme_minimal()

print(plot_residuals_terms)

# 16d. Standard linear model diagnostic plots, as in the handout, for the
#      worst-fitted event. Residuals vs fitted (missing structure), Q-Q
#      (normality), scale-location (constant variance), leverage.
#      Normality and independence are not essential here (Lecture 5,
#      slides 15-17): the aim is a good basis, not inference.
worst <- which.max(poly_results$rmse)
par(mfrow = c(2, 2))
plot(poly_models[[worst]])
par(mfrow = c(1, 1))

# 16e. Coefficient table for the same event, as in the handout.
print(summary(poly_models[[worst]]))

# 16f. Correlation between fitted coefficients across events (Lecture 5,
#      slide 12). Strong collinearity will complicate modelling the
#      uncertainty in these coefficients between hydrographs.
cat("\nCorrelation between shape parameters across events:\n")
print(round(cor(poly_coefs), 2))

# ==============================================================================
# METHOD 2: PRINCIPAL COMPONENT ANALYSIS (Lecture 4, slides 14-22;
#           lecture PCA example in R)
# ==============================================================================
# Instead of fixing the basis in advance, PCA estimates one from the curves.
# Each standardised curve is reconstructed as the mean curve plus a weighted
# sum of k principal component curves:
#
#   y_i(u) ~ ybar(u) + sum_{j=1}^{k} C_ij a_j(u)
#
#   a_j  : loadings of principal component j (a curve over the 64 grid points)
#   C_ij : score of event i on component j (how much of component j it has)
#
# The scores C_i1, ..., C_ik are each event's shape parameters.
# Following the lecture example, curves are centred but not scaled, since
# every grid point is already on the same scale.
# ==============================================================================


# ------------------------------------------------------------------------------
# 17. Fit PCA (as in the lecture example)
# ------------------------------------------------------------------------------

# PCA algorithm: centre but don't scale, as already on the same scale
pca_result <- prcomp(Y, center = TRUE, scale. = FALSE)

# pca_result$rotation: column j holds the loadings a_j over the 64 grid points
# pca_result$x       : row i holds the scores (C_i1, C_i2, ...) for event i
dim(pca_result$rotation)
dim(pca_result$x)

# 17a. Plot the loadings of the first principal components as functions of
#      time (functional view, as in the lecture example)
plot(grid_u, pca_result$rotation[, 1], type = "l", col = "red",
     ylim = range(pca_result$rotation[, 1:3]),
     xlab = "Rescaled time u", ylab = "Loading",
     main = "Loadings for PC1, PC2 and PC3")
lines(grid_u, pca_result$rotation[, 2], col = "blue")
lines(grid_u, pca_result$rotation[, 3], col = "darkgreen")
abline(h = 0, lty = 3)
legend("bottomright", legend = c("PC1", "PC2", "PC3"),
       col = c("red", "blue", "darkgreen"), lty = 1)

# 17b. Variance of each component (as in the lecture example)
plot(pca_result, main = "Variance of each principal component")

# 17c. Cumulative proportion of variance explained (Lecture 5, slide 18):
#      choose the smallest k that captures an acceptable proportion
#      (Lecture 4 suggests, e.g., 95%).
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cum_var       <- cumsum(var_explained)

print(round(data.frame(k = seq_along(cum_var), proportion = var_explained,
                       cumulative = cum_var)[1:10, ], 3))

plot(seq_along(cum_var), cum_var, type = "b", pch = 19,
     xlab = "Number of components k", ylab = "Cumulative proportion of variance",
     main = "Cumulative variance explained", ylim = c(0, 1))
abline(h = 0.95, lty = 2, col = "red")

cat("Smallest k explaining at least 95% of variance:", which(cum_var >= 0.95)[1], "\n")


# ------------------------------------------------------------------------------
# 18. Choosing k: in-sample and leave-one-out errors (Lecture 5, slide 18)
# ------------------------------------------------------------------------------
# Reconstruction, as in the lecture example:
#   recon = mean curve + scores[, 1:k] %*% t(loadings[, 1:k])

pca_reconstruct <- function(k) {
  mean_hydro <- colMeans(Y)
  recon <- pca_result$x[, 1:k, drop = FALSE] %*%
    t(pca_result$rotation[, 1:k, drop = FALSE])
  sweep(recon, 2, mean_hydro, FUN = "+")
}

# Leave-one-out cross-validation (Lecture 5, slide 18): remove a hydrograph,
# recompute the basis from the rest, then project and reconstruct the
# held-out curve. The held-out curve is centred with the TRAINING mean, and
# its scores are found by projecting it onto the training loadings.
pca_fit_and_predict <- function(k) {
  function(Y_train, y_new) {
    pca_train  <- prcomp(Y_train, center = TRUE, scale. = FALSE)
    mean_train <- colMeans(Y_train)
    loadings   <- pca_train$rotation[, 1:k, drop = FALSE]
    scores_new <- (y_new - mean_train) %*% loadings
    as.vector(mean_train + scores_new %*% t(loadings))
  }
}

# Compare k = 1 to 10 (or fewer if there are few events)
k_values <- 1:min(10, nrow(Y) - 2)

k_table <- do.call(rbind, lapply(k_values, function(k) {
  rbind(
    data.frame(k = k, type = "in-sample",
               t(summarise_metrics(evaluate_fits(pca_reconstruct(k))))),
    data.frame(k = k, type = "leave-one-out",
               t(summarise_metrics(evaluate_fits(loocv(pca_fit_and_predict(k))))))
  )
}))

print(k_table, digits = 3)

# In-sample errors always fall as k grows; leave-one-out errors show whether
# extra components still help for a NEW event. Choose k where leave-one-out
# errors stop improving meaningfully, alongside the variance threshold.
plot_k <- k_table %>%
  pivot_longer(c(rmse, peak_error, timing_error, volume_error),
               names_to = "metric", values_to = "error") %>%
  ggplot(aes(k, error, colour = type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_continuous(breaks = k_values) +
  scale_colour_manual(values = c("in-sample" = "#174B67", "leave-one-out" = "#AA1945")) +
  labs(
    title = "PCA: average absolute error by number of components",
    x = "Number of components k",
    y = "Average absolute error",
    colour = NULL
  ) +
  theme_minimal()

print(plot_k)


# ------------------------------------------------------------------------------
# 19. Fit at the chosen k and validate (Lecture 5, slide 18)
# ------------------------------------------------------------------------------
# UPDATE k_chosen after looking at sections 17c and 18.

k_chosen <- 3

Y_pca <- pca_reconstruct(k_chosen)

# Shape parameters: one row per event, k_chosen scores
pca_scores <- pca_result$x[, 1:k_chosen, drop = FALSE]
rownames(pca_scores) <- rownames(Y)

pca_results       <- evaluate_fits(Y_pca)
pca_results_loocv <- evaluate_fits(loocv(pca_fit_and_predict(k_chosen)))

cat("\nAverage absolute errors (in-sample and leave-one-out):\n")
print(rbind(
  baseline_in_sample   = summarise_metrics(baseline_in_sample),
  baseline_loocv       = summarise_metrics(baseline_loocv),
  polynomial_in_sample = summarise_metrics(poly_results),
  pca_in_sample        = summarise_metrics(pca_results),
  pca_loocv            = summarise_metrics(pca_results_loocv)
))

# Physical check: reconstructed curves should not go below zero
cat("Grid points with negative reconstructed flow:", sum(Y_pca < 0), "\n")

# 19a. Reconstructed vs original curves (as in the lecture example), for
#      different types of event: earliest peak, latest peak, longest event
#      and the worst-reconstructed event.
pca_example_ids <- unique(c(
  events$event_id[which.min(events$peak_position)],
  events$event_id[which.max(events$peak_position)],
  events$event_id[which.max(events$width_h)],
  pca_results$event_id[which.max(pca_results$rmse)]
))

par(mfrow = c(2, 2))
for (id in pca_example_ids) {
  i <- which(events$event_id == id)
  plot(grid_u, Y[i, ], type = "l", col = "black", lwd = 2,
       ylim = range(c(Y[i, ], Y_pca[i, ])),
       xlab = "Rescaled time u", ylab = "Discharge / peak discharge",
       main = sprintf("Event %d: original vs PCA (k = %d)", id, k_chosen))
  lines(grid_u, Y_pca[i, ], col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Original", "Recon"),
         col = c("black", "red"), lty = c(1, 2), lwd = 2, cex = 0.8)
}
par(mfrow = c(1, 1))

# 19b. Residuals against time, all events: where is the fit deficient?
#      Bold line: average residual.
pca_residuals_long <- data.frame(
  event_id = rep(events$event_id, each = ncol(Y)),
  u        = rep(grid_u, times = nrow(Y)),
  residual = as.vector(t(Y - Y_pca))
)

plot_pca_residuals <- ggplot(pca_residuals_long, aes(u, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(group = event_id), colour = "#5298BC", alpha = 0.3) +
  stat_summary(fun = mean, geom = "line", colour = "#AA1945", linewidth = 1) +
  labs(
    title = sprintf("PCA (k = %d): residuals against time", k_chosen),
    x = "Rescaled time u",
    y = "Original - reconstructed (standardised)"
  ) +
  theme_minimal()

print(plot_pca_residuals)

# 19c. Outliers in the scores (Lecture 5, slide 18): pairs plot of the
#      retained scores. Points far from the rest are unusual hydrographs.
pairs(pca_scores, pch = 19, col = "#174B67",
      main = sprintf("Scores on the first %d principal components", k_chosen))

# 19d. Variance and correlation of the scores. Scores are uncorrelated by
#      construction (Lecture 4, slide 16), which simplifies modelling their
#      uncertainty later; the variances show how much each component varies.
cat("\nVariance of each retained score:\n")
print(round(apply(pca_scores, 2, var), 4))
cat("\nCorrelation between retained scores:\n")
print(round(cor(pca_scores), 2))

# ==============================================================================
# METHOD 3: WAVELETS (Lecture 4, slides 24-27; lecture wavelet example in R)
# ==============================================================================
# Each standardised curve is decomposed with the Discrete Wavelet Transform
# (DWT) into a smooth trend plus detail at each scale (Lecture 4, slide 26):
#
#   y(u) = sum_k a_{J,k} phi_{J,k}(u)  +  sum_{j<=J} sum_k d_{j,k} psi_{j,k}(u)
#           smooth (scaling) coeffs        detail (wavelet) coeffs
#
# Most detail coefficients are very small, so only the largest few are kept.
# As in the lecture example, the Daubechies "d4" filter, 3 levels and
# periodic boundaries are used. The grid already has 64 = 2^6 points, so no
# padding is needed (the lecture example pads its 50 points up to 64).
#
# One change from the single-curve lecture example: the example keeps the
# largest coefficients of ONE curve. With many curves, Lecture 4 (slide 26)
# says the same basis functions must be kept across all hydrographs, so the
# kept positions are chosen ONCE, as those with the largest average size
# across events. Every event then has the same parameters: the smooth
# coefficients plus the kept detail coefficients.
# ==============================================================================

#install.packages("wavelets")
library(wavelets)

wav_filter <- "d4"   # as in the lecture example
wav_levels <- 3      # as in the lecture example

stopifnot(log2(ncol(Y)) == round(log2(ncol(Y))))   # DWT needs length 2^n


# ------------------------------------------------------------------------------
# 20. Decompose every curve (as in the lecture example)
# ------------------------------------------------------------------------------

# 20a. DWT of one curve, returning detail and smooth coefficients
wavelet_coefs <- function(y, filter = wav_filter, levels = wav_levels) {
  wt <- dwt(as.numeric(y), filter = filter, n.levels = levels,
            boundary = "periodic")
  list(detail = unlist(wt@W),        # all detail coefficients, fine to coarse
       smooth = wt@V[[levels]])      # only the last smooth vector is needed
}

# 20b. Rebuild a curve from detail and smooth coefficients with the inverse
#      DWT, by putting them into a template decomposition
wavelet_curve <- function(detail, smooth, filter = wav_filter, levels = wav_levels) {
  wt   <- dwt(as.numeric(Y[1, ]), filter = filter, n.levels = levels,
              boundary = "periodic")
  lens <- lengths(wt@W)
  wt@W <- lapply(split(detail, rep(seq_along(lens), lens)), matrix, ncol = 1)
  wt@V[[levels]] <- matrix(smooth, ncol = 1)
  as.numeric(idwt(wt))
}

# 20c. Plot the decomposition of one event (as in the lecture example)
wt_example <- dwt(as.numeric(Y[1, ]), filter = wav_filter,
                  n.levels = wav_levels, boundary = "periodic")

par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
plot(seq(0, 1, length.out = length(wt_example@W[[1]])), wt_example@W[[1]],
     type = "l", col = "red", xlab = "", ylab = "", main = "L1 detail (fine)")
plot(seq(0, 1, length.out = length(wt_example@W[[2]])), wt_example@W[[2]],
     type = "l", col = "blue", xlab = "", ylab = "", main = "L2 detail")
plot(seq(0, 1, length.out = length(wt_example@W[[3]])), wt_example@W[[3]],
     type = "l", col = "darkgreen", xlab = "", ylab = "", main = "L3 detail (coarse)")
plot(seq(0, 1, length.out = length(wt_example@V[[3]])), wt_example@V[[3]],
     type = "l", col = "black", xlab = "", ylab = "", main = "L3 smooth (trend)")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# 20d. Check the full decomposition reconstructs the curve (as in the example)
recon_full <- as.numeric(idwt(wt_example))
cat("Max error of full reconstruction:", max(abs(Y[1, ] - recon_full)), "\n")

# 20e. Coefficients for every event: one row per event
all_coefs  <- lapply(seq_len(nrow(Y)), function(i) wavelet_coefs(Y[i, ]))
D_all      <- t(sapply(all_coefs, `[[`, "detail"))   # N x 56 detail coeffs
S_all      <- t(sapply(all_coefs, `[[`, "smooth"))   # N x 8 smooth coeffs

cat("Detail coefficients per event:", ncol(D_all),
    "| smooth coefficients per event:", ncol(S_all), "\n")


# ------------------------------------------------------------------------------
# 21. Choosing how many detail coefficients to keep (Lecture 5, slide 19)
# ------------------------------------------------------------------------------
# Kept positions = the keep_k positions with the largest average absolute
# detail coefficient across the events used for fitting. All other detail
# coefficients are set to zero (thresholding), as in the lecture example.

choose_positions <- function(D, keep_k) {
  order(colMeans(abs(D)), decreasing = TRUE)[seq_len(keep_k)]
}

reconstruct_kept <- function(detail, smooth, positions,
                             filter = wav_filter, levels = wav_levels) {
  kept <- rep(0, length(detail))
  kept[positions] <- detail[positions]
  wavelet_curve(kept, smooth, filter, levels)
}

# In-sample: positions chosen from all events
wavelet_reconstruct <- function(keep_k) {
  positions <- choose_positions(D_all, keep_k)
  t(sapply(seq_len(nrow(Y)), function(i) {
    reconstruct_kept(D_all[i, ], S_all[i, ], positions)
  }))
}

# Held-out check (Lecture 5, slide 19): positions chosen WITHOUT the event,
# then applied to it, using the shared leave-one-out function
wavelet_fit_and_predict <- function(keep_k) {
  function(Y_train, y_new) {
    D_train   <- t(apply(Y_train, 1, function(y) wavelet_coefs(y)$detail))
    positions <- choose_positions(D_train, keep_k)
    new       <- wavelet_coefs(y_new)
    reconstruct_kept(new$detail, new$smooth, positions)
  }
}

# Reconstruction error as a function of the number of retained coefficients
keep_values <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)

keep_table <- do.call(rbind, lapply(keep_values, function(k) {
  rbind(
    data.frame(keep_k = k, type = "in-sample",
               t(summarise_metrics(evaluate_fits(wavelet_reconstruct(k))))),
    data.frame(keep_k = k, type = "leave-one-out",
               t(summarise_metrics(evaluate_fits(loocv(wavelet_fit_and_predict(k))))))
  )
}))

print(keep_table, digits = 3)

plot_keep <- keep_table %>%
  pivot_longer(c(rmse, peak_error, timing_error, volume_error),
               names_to = "metric", values_to = "error") %>%
  ggplot(aes(keep_k, error, colour = type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y") +
  scale_colour_manual(values = c("in-sample" = "#174B67", "leave-one-out" = "#AA1945")) +
  labs(
    title = sprintf("Wavelets (%s, %d levels): error by number of kept detail coefficients",
                    wav_filter, wav_levels),
    x = "Number of detail coefficients kept",
    y = "Average absolute error",
    colour = NULL
  ) +
  theme_minimal()

print(plot_keep)


# ------------------------------------------------------------------------------
# 22. Fit at the chosen number of coefficients and validate (Lecture 5, slide 19)
# ------------------------------------------------------------------------------
# UPDATE keep_chosen after looking at section 21.

keep_chosen <- 5

kept_positions <- choose_positions(D_all, keep_chosen)
Y_wav          <- wavelet_reconstruct(keep_chosen)

# Shape parameters: smooth coefficients plus the kept detail coefficients
wav_params <- cbind(S_all, D_all[, kept_positions, drop = FALSE])
colnames(wav_params) <- c(paste0("smooth_", seq_len(ncol(S_all))),
                          paste0("detail_", kept_positions))
rownames(wav_params) <- rownames(Y)

cat("Parameters per event:", ncol(wav_params),
    "(", ncol(S_all), "smooth +", keep_chosen, "detail )\n")

wav_results       <- evaluate_fits(Y_wav)
wav_results_loocv <- evaluate_fits(loocv(wavelet_fit_and_predict(keep_chosen)))

cat("\nAverage absolute errors (in-sample and leave-one-out):\n")
print(rbind(
  baseline_in_sample   = summarise_metrics(baseline_in_sample),
  polynomial_in_sample = summarise_metrics(poly_results),
  pca_in_sample        = summarise_metrics(pca_results),
  pca_loocv            = summarise_metrics(pca_results_loocv),
  wavelet_in_sample    = summarise_metrics(wav_results),
  wavelet_loocv        = summarise_metrics(wav_results_loocv)
))

cat("Grid points with negative reconstructed flow:", sum(Y_wav < 0), "\n")

# 22a. Reconstructions for a range of hydrograph shapes (as in the lecture
#      example): earliest peak, latest peak, longest event, worst fit.
wav_example_ids <- unique(c(
  events$event_id[which.min(events$peak_position)],
  events$event_id[which.max(events$peak_position)],
  events$event_id[which.max(events$width_h)],
  wav_results$event_id[which.max(wav_results$rmse)]
))

par(mfrow = c(2, 2))
for (id in wav_example_ids) {
  i <- which(events$event_id == id)
  plot(grid_u, Y[i, ], type = "l", col = "black", lwd = 2,
       ylim = range(c(Y[i, ], Y_wav[i, ])),
       xlab = "Rescaled time u", ylab = "Discharge / peak discharge",
       main = sprintf("Event %d: %d detail coeffs", id, keep_chosen))
  lines(grid_u, Y_wav[i, ], col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Original", "Recon"),
         col = c("black", "red"), lty = c(1, 2), lwd = 2, cex = 0.8)
}
par(mfrow = c(1, 1))

# 22b. Residuals against time, all events (as for regression)
wav_residuals_long <- data.frame(
  event_id = rep(events$event_id, each = ncol(Y)),
  u        = rep(grid_u, times = nrow(Y)),
  residual = as.vector(t(Y - Y_wav))
)

plot_wav_residuals <- ggplot(wav_residuals_long, aes(u, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(group = event_id), colour = "#5298BC", alpha = 0.3) +
  stat_summary(fun = mean, geom = "line", colour = "#AA1945", linewidth = 1) +
  labs(
    title = sprintf("Wavelets (%d detail coeffs): residuals against time", keep_chosen),
    x = "Rescaled time u",
    y = "Original - reconstructed (standardised)"
  ) +
  theme_minimal()

print(plot_wav_residuals)

# 22c. Sensitivity to assumptions: wavelet choice and resolution levels
#      (Lecture 5, slide 19). Same number of kept detail coefficients
#      throughout. Note the number of smooth coefficients (and so parameters)
#      changes with the number of levels: 64 / 2^levels.
sensitivity <- expand.grid(filter = c("haar", "d4", "la8"), levels = 3:5,
                           stringsAsFactors = FALSE)

sensitivity_table <- do.call(rbind, lapply(seq_len(nrow(sensitivity)), function(r) {
  f <- sensitivity$filter[r]
  l <- sensitivity$levels[r]
  coefs     <- lapply(seq_len(nrow(Y)), function(i) wavelet_coefs(Y[i, ], f, l))
  D         <- t(sapply(coefs, `[[`, "detail"))
  S         <- t(sapply(coefs, `[[`, "smooth"))
  positions <- choose_positions(D, keep_chosen)
  Y_fit     <- t(sapply(seq_len(nrow(Y)), function(i) {
    reconstruct_kept(D[i, ], S[i, ], positions, f, l)
  }))
  data.frame(filter = f, levels = l, n_params = ncol(S) + keep_chosen,
             t(summarise_metrics(evaluate_fits(Y_fit))))
}))

print(sensitivity_table, digits = 3)

# 22d. Correlation and variance of the retained coefficients across events
#      (Lecture 5, slide 19): informs the later uncertainty modelling.
cat("\nVariance of each retained coefficient:\n")
print(round(apply(wav_params, 2, var), 4))
cat("\nCorrelation between retained coefficients:\n")
print(round(cor(wav_params), 2))

# ==============================================================================
# ADDITIONAL METHODS BEYOND THE LECTURE EXAMPLES
# ==============================================================================
#   23. Known hydrograph function (gamma-type)   - Lecture 4, slide 12 tip
#   24. B-spline basis                           - Ramsay & Silverman, FDA
#   25. Landmark registration, then PCA          - Ramsay & Silverman, FDA
#   26. PCA on log flow (positive functions)     - Ramsay & Silverman, FDA
#   27. Comparison of all methods                - Lecture 5, slides 21-25
# ==============================================================================

library(splines)   # part of base R; provides bs() for B-splines

# Helpers used by every method below ------------------------------------------

# Example reconstructions for different types of event (as in sections 16a,
# 19a and 22a): earliest peak, latest peak, longest event, worst fit
plot_method_examples <- function(Y_fit, results, title) {
  ids <- unique(c(
    events$event_id[which.min(events$peak_position)],
    events$event_id[which.max(events$peak_position)],
    events$event_id[which.max(events$width_h)],
    results$event_id[which.max(results$rmse)]
  ))
  par(mfrow = c(2, 2))
  for (id in ids) {
    i <- which(events$event_id == id)
    plot(grid_u, Y[i, ], type = "l", col = "black", lwd = 2,
         ylim = range(c(Y[i, ], Y_fit[i, ])),
         xlab = "Rescaled time u", ylab = "Discharge / peak discharge",
         main = sprintf("%s: event %d", title, id))
    lines(grid_u, Y_fit[i, ], col = "red", lwd = 2, lty = 2)
    legend("topright", legend = c("Original", "Fitted"),
           col = c("black", "red"), lty = c(1, 2), lwd = 2, cex = 0.8)
  }
  par(mfrow = c(1, 1))
}

# Residuals against time for all events, with the average residual in bold
plot_method_residuals <- function(Y_fit, title) {
  df <- data.frame(
    event_id = rep(events$event_id, each = ncol(Y)),
    u        = rep(grid_u, times = nrow(Y)),
    residual = as.vector(t(Y - Y_fit))
  )
  p <- ggplot(df, aes(u, residual)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(aes(group = event_id), colour = "#5298BC", alpha = 0.3) +
    stat_summary(fun = mean, geom = "line", colour = "#AA1945", linewidth = 1) +
    labs(title = paste(title, ": residuals against time"),
         x = "Rescaled time u", y = "Observed - fitted (standardised)") +
    theme_minimal()
  print(p)
}


# ==============================================================================
# 23. KNOWN HYDROGRAPH FUNCTION (gamma-type)
# ==============================================================================
# Lecture 4 (slide 12) suggests fitting a known function that approximates
# hydrographs, because its parameters may be interpretable. A standard
# choice is the gamma-type hydrograph. Our curves start and end at the
# trough level rather than zero, so a baseflow term b is added:
#
#   y(u) = b + (1 - b) * (u / u_p)^m * exp( m * (1 - u / u_p) )
#
#   b   : baseflow level, as a fraction of peak flow
#   u_p : time of the peak on [0, 1]
#   m   : sharpness (larger m = narrower, sharper peak)
#
# The curve equals 1 exactly at u_p. It is nonlinear in its parameters, so
# it is fitted by nonlinear least squares with nls(), with bounds keeping
# the parameters physically sensible. Being a fixed function, leave-one-out
# errors equal in-sample errors (as for the polynomial).

gamma_curve <- function(u, b, up, m) {
  b + (1 - b) * (u / up)^m * exp(m * (1 - u / up))
}

fit_gamma <- function(y) {
  start <- list(b = max(min(y), 0.01), up = grid_u[which.max(y)], m = 5)
  fit <- tryCatch(
    nls(y ~ gamma_curve(u, b, up, m),
        data = data.frame(y = y, u = grid_u), start = start,
        algorithm = "port",
        lower = c(b = 0, up = 0.01, m = 0.1),
        upper = c(b = 0.99, up = 0.99, m = 500)),
    error = function(e) NULL
  )
  fit
}

gamma_models <- lapply(seq_len(nrow(Y)), function(i) fit_gamma(Y[i, ]))
gamma_failed <- sapply(gamma_models, is.null)
cat("Gamma fits that failed to converge:", sum(gamma_failed), "\n")

# Shape parameters: one row per event (NA where the fit failed)
gamma_params <- t(sapply(gamma_models, function(fit) {
  if (is.null(fit)) c(b = NA, up = NA, m = NA) else coef(fit)
}))
rownames(gamma_params) <- rownames(Y)

Y_gamma <- t(sapply(seq_len(nrow(Y)), function(i) {
  if (gamma_failed[i]) rep(NA, ncol(Y)) else fitted(gamma_models[[i]])
}))

# Errors (events with failed fits are reported separately, not hidden)
gamma_ok      <- which(!gamma_failed)
gamma_results <- evaluate_fits(replace(Y_gamma, is.na(Y_gamma), 0))[gamma_ok, ]

cat("\nGamma-type function, average absolute errors (converged fits):\n")
print(summarise_metrics(gamma_results))

# Validation: example fits, residuals, parameters
Y_gamma_plot <- Y_gamma
Y_gamma_plot[gamma_failed, ] <- Y[gamma_failed, ]   # failed fits not plotted as fits
plot_method_examples(Y_gamma_plot, evaluate_fits(Y_gamma_plot), "Gamma-type")
plot_method_residuals(Y_gamma_plot, "Gamma-type")   # failed fits show zero residual

# The parameters are interpretable: look at their spread and relationships
cat("\nSummary of gamma-type parameters:\n")
print(summary(gamma_params))
pairs(gamma_params[gamma_ok, ], pch = 19, col = "#174B67",
      main = "Gamma-type parameters across events")
cat("\nCorrelation between gamma-type parameters:\n")
print(round(cor(gamma_params[gamma_ok, ]), 2))


# ==============================================================================
# 24. B-SPLINE BASIS
# ==============================================================================
# Ramsay & Silverman's standard alternative to polynomials. A B-spline basis
# joins short cubic pieces smoothly at "knots", so each basis function only
# affects one part of the curve. This avoids the polynomial's problems
# (wiggles in the tail, ill-conditioning) while still being an ordinary
# linear model fitted with lm(), exactly like the polynomial:
#
#   y_i(u) = theta_i0 + sum_{k=1}^{df} theta_ik B_k(u) + error
#
# df (number of basis functions) plays the role of the polynomial degree.
# Knots are equally spaced over [0, 1]. As a fixed basis, leave-one-out
# errors equal in-sample errors; the degree check uses train/test events.

fit_bspline_model <- function(y, df) {
  lm(y ~ bs(u, df = df), data = data.frame(y = y, u = grid_u))
}

df_values <- 4:15

bspline_table <- do.call(rbind, lapply(df_values, function(d) {
  Y_fit <- t(sapply(seq_len(nrow(Y)), function(i) fitted(fit_bspline_model(Y[i, ], d))))
  res   <- evaluate_fits(Y_fit)
  rbind(
    data.frame(df = d, set = "train",
               t(summarise_metrics(res[events$event_id %in% train_ids, ]))),
    data.frame(df = d, set = "test",
               t(summarise_metrics(res[events$event_id %in% test_ids, ])))
  )
}))

print(bspline_table, digits = 3)

plot_bspline_df <- bspline_table %>%
  pivot_longer(c(rmse, peak_error, timing_error, volume_error),
               names_to = "metric", values_to = "error") %>%
  ggplot(aes(df, error, colour = set)) +
  geom_line() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_continuous(breaks = df_values) +
  scale_colour_manual(values = c(train = "#174B67", test = "#AA1945")) +
  labs(title = "B-splines: average absolute error by number of basis functions",
       x = "Number of basis functions (df)", y = "Average absolute error",
       colour = "Events") +
  theme_minimal()

print(plot_bspline_df)

# UPDATE df_chosen after looking at the table and plot above.
df_chosen <- 8

bspline_models <- lapply(seq_len(nrow(Y)), function(i) fit_bspline_model(Y[i, ], df_chosen))
Y_bspline      <- t(sapply(bspline_models, fitted))

bspline_coefs <- t(sapply(bspline_models, coef))
colnames(bspline_coefs) <- paste0("theta_", 0:df_chosen)
rownames(bspline_coefs) <- rownames(Y)

bspline_results <- evaluate_fits(Y_bspline)

cat("\nB-splines (df =", df_chosen, "), average absolute errors:\n")
print(summarise_metrics(bspline_results))
cat("Grid points with negative fitted flow:", sum(Y_bspline < 0), "\n")

plot_method_examples(Y_bspline, bspline_results, sprintf("B-spline (df = %d)", df_chosen))
plot_method_residuals(Y_bspline, sprintf("B-spline (df = %d)", df_chosen))

# Standard linear model diagnostic plots for the worst-fitted event
par(mfrow = c(2, 2))
plot(bspline_models[[which.max(bspline_results$rmse)]])
par(mfrow = c(1, 1))

cat("\nCorrelation between B-spline coefficients across events:\n")
print(round(cor(bspline_coefs), 2))


# ==============================================================================
# 25. LANDMARK REGISTRATION, THEN PCA
# ==============================================================================
# Ramsay & Silverman separate timing (phase) variation from shape
# (amplitude) variation by registering curves before analysing them. Here
# the landmark is the peak: each curve is warped piecewise-linearly so that
# its peak moves to a common position u_star (the average peak position),
# keeping the start (0) and end (1) fixed:
#
#   rising part  [0, p_i] -> [0, u_star]
#   falling part [p_i, 1] -> [u_star, 1]
#
# PCA is then applied to the registered curves. Each event's parameters are
# its peak position p_i plus its k scores. To evaluate, reconstructed curves
# are warped back to each event's own timing before computing errors.

peak_pos <- function(y) grid_u[which.max(y)]

# Registered curve: value of y at the time that maps onto each grid point
register_curve <- function(y, p, u_star) {
  u <- ifelse(grid_u <= u_star,
              grid_u * p / u_star,
              p + (grid_u - u_star) * (1 - p) / (1 - u_star))
  approx(grid_u, y, xout = u, rule = 2)$y
}

# Inverse: warp a registered curve back to an event with peak position p
unregister_curve <- function(y_reg, p, u_star) {
  v <- ifelse(grid_u <= p,
              grid_u * u_star / p,
              u_star + (grid_u - p) * (1 - u_star) / (1 - p))
  approx(grid_u, y_reg, xout = v, rule = 2)$y
}

reg_peaks <- apply(Y, 1, peak_pos)
u_star    <- mean(reg_peaks)

Y_reg <- t(sapply(seq_len(nrow(Y)), function(i) register_curve(Y[i, ], reg_peaks[i], u_star)))

# Registered curves: peaks now line up at u_star
matplot(grid_u, t(Y_reg), type = "l", lty = 1, col = adjustcolor("#5298BC", 0.4),
        xlab = "Registered time", ylab = "Discharge / peak discharge",
        main = "Curves after landmark registration (peaks aligned)")
lines(grid_u, colMeans(Y_reg), col = "#174B67", lwd = 3)
abline(v = u_star, lty = 2, col = "#AA1945")

# PCA on the registered curves (as in section 17)
reg_pca <- prcomp(Y_reg, center = TRUE, scale. = FALSE)
reg_cum_var <- cumsum(reg_pca$sdev^2) / sum(reg_pca$sdev^2)

cat("\nCumulative variance explained, unregistered vs registered PCA:\n")
print(round(data.frame(k = 1:10, unregistered = cum_var[1:10],
                       registered = reg_cum_var[1:10]), 3))

# Reconstruction on each event's own timing, using k components
reg_reconstruct <- function(k) {
  recon_reg <- sweep(reg_pca$x[, 1:k, drop = FALSE] %*%
                       t(reg_pca$rotation[, 1:k, drop = FALSE]),
                     2, colMeans(Y_reg), FUN = "+")
  t(sapply(seq_len(nrow(Y)), function(i) unregister_curve(recon_reg[i, ], reg_peaks[i], u_star)))
}

# Leave-one-out: registration target, mean and components from the other
# events only; the left-out curve is registered, projected, then unwarped
reg_fit_and_predict <- function(k) {
  function(Y_train, y_new) {
    p_train  <- apply(Y_train, 1, peak_pos)
    u_s      <- mean(p_train)
    Y_tr_reg <- t(sapply(seq_len(nrow(Y_train)), function(i) {
      register_curve(Y_train[i, ], p_train[i], u_s)
    }))
    pca_tr   <- prcomp(Y_tr_reg, center = TRUE, scale. = FALSE)
    mean_tr  <- colMeans(Y_tr_reg)
    loadings <- pca_tr$rotation[, 1:k, drop = FALSE]
    p_new    <- peak_pos(y_new)
    y_reg    <- register_curve(y_new, p_new, u_s)
    recon    <- as.vector(mean_tr + ((y_reg - mean_tr) %*% loadings) %*% t(loadings))
    unregister_curve(recon, p_new, u_s)
  }
}

reg_k_table <- do.call(rbind, lapply(k_values, function(k) {
  rbind(
    data.frame(k = k, type = "in-sample",
               t(summarise_metrics(evaluate_fits(reg_reconstruct(k))))),
    data.frame(k = k, type = "leave-one-out",
               t(summarise_metrics(evaluate_fits(loocv(reg_fit_and_predict(k))))))
  )
}))

print(reg_k_table, digits = 3)

# UPDATE reg_k_chosen after looking at the table above.
reg_k_chosen <- 3

Y_regpca            <- reg_reconstruct(reg_k_chosen)
regpca_results      <- evaluate_fits(Y_regpca)
regpca_results_loocv <- evaluate_fits(loocv(reg_fit_and_predict(reg_k_chosen)))

# Parameters: peak position plus the k scores
regpca_params <- cbind(peak_position = reg_peaks,
                       reg_pca$x[, 1:reg_k_chosen, drop = FALSE])

plot_method_examples(Y_regpca, regpca_results,
                     sprintf("Registered PCA (k = %d)", reg_k_chosen))
plot_method_residuals(Y_regpca, sprintf("Registered PCA (k = %d)", reg_k_chosen))

cat("\nCorrelation between registered-PCA parameters:\n")
print(round(cor(regpca_params), 2))


# ==============================================================================
# 26. PCA ON LOG FLOW (positive functions)
# ==============================================================================
# Ramsay & Silverman represent functions that must stay positive through
# their logarithm. PCA is fitted to log(y) and reconstructions are
# exponentiated, so reconstructed (and later simulated) flows can never be
# negative. The same idea could be applied to any method above.

Y_log <- log(Y)   # standardised flows are always > 0

log_pca <- prcomp(Y_log, center = TRUE, scale. = FALSE)

log_reconstruct <- function(k) {
  exp(sweep(log_pca$x[, 1:k, drop = FALSE] %*%
              t(log_pca$rotation[, 1:k, drop = FALSE]),
            2, colMeans(Y_log), FUN = "+"))
}

log_fit_and_predict <- function(k) {
  function(Y_train, y_new) {
    L_train  <- log(Y_train)
    pca_tr   <- prcomp(L_train, center = TRUE, scale. = FALSE)
    mean_tr  <- colMeans(L_train)
    loadings <- pca_tr$rotation[, 1:k, drop = FALSE]
    exp(as.vector(mean_tr + ((log(y_new) - mean_tr) %*% loadings) %*% t(loadings)))
  }
}

log_k_table <- do.call(rbind, lapply(k_values, function(k) {
  rbind(
    data.frame(k = k, type = "in-sample",
               t(summarise_metrics(evaluate_fits(log_reconstruct(k))))),
    data.frame(k = k, type = "leave-one-out",
               t(summarise_metrics(evaluate_fits(loocv(log_fit_and_predict(k))))))
  )
}))

print(log_k_table, digits = 3)

# UPDATE log_k_chosen after looking at the table above.
log_k_chosen <- 3

Y_logpca             <- log_reconstruct(log_k_chosen)
logpca_results       <- evaluate_fits(Y_logpca)
logpca_results_loocv <- evaluate_fits(loocv(log_fit_and_predict(log_k_chosen)))

plot_method_examples(Y_logpca, logpca_results, sprintf("Log-flow PCA (k = %d)", log_k_chosen))
plot_method_residuals(Y_logpca, sprintf("Log-flow PCA (k = %d)", log_k_chosen))


# ==============================================================================
# 27. COMPARISON OF ALL METHODS (Lecture 5, slides 21-25)
# ==============================================================================
# Accuracy (in-sample and, where the basis is learnt from the data,
# leave-one-out), parsimony (number of shape parameters per event) and a
# physical check (negative flows). For fixed-basis methods (polynomial,
# B-spline, gamma), leave-one-out errors equal in-sample errors.

comparison <- rbind(
  data.frame(method = "Baseline (mean shape)", n_params = 0, type = "in-sample",
             t(summarise_metrics(baseline_in_sample)), negative_points = 0),
  data.frame(method = "Baseline (mean shape)", n_params = 0, type = "leave-one-out",
             t(summarise_metrics(baseline_loocv)), negative_points = 0),
  data.frame(method = sprintf("Polynomial (p = %d)", p_chosen), n_params = p_chosen + 1,
             type = "in-sample", t(summarise_metrics(poly_results)),
             negative_points = sum(Y_poly < 0)),
  data.frame(method = "Gamma-type function", n_params = 3, type = "in-sample",
             t(summarise_metrics(gamma_results)),
             negative_points = sum(Y_gamma < 0, na.rm = TRUE)),
  data.frame(method = sprintf("B-spline (df = %d)", df_chosen), n_params = df_chosen + 1,
             type = "in-sample", t(summarise_metrics(bspline_results)),
             negative_points = sum(Y_bspline < 0)),
  data.frame(method = sprintf("PCA (k = %d)", k_chosen), n_params = k_chosen,
             type = "in-sample", t(summarise_metrics(pca_results)),
             negative_points = sum(Y_pca < 0)),
  data.frame(method = sprintf("PCA (k = %d)", k_chosen), n_params = k_chosen,
             type = "leave-one-out", t(summarise_metrics(pca_results_loocv)),
             negative_points = NA),
  data.frame(method = sprintf("Registered PCA (k = %d)", reg_k_chosen),
             n_params = reg_k_chosen + 1, type = "in-sample",
             t(summarise_metrics(regpca_results)), negative_points = sum(Y_regpca < 0)),
  data.frame(method = sprintf("Registered PCA (k = %d)", reg_k_chosen),
             n_params = reg_k_chosen + 1, type = "leave-one-out",
             t(summarise_metrics(regpca_results_loocv)), negative_points = NA),
  data.frame(method = sprintf("Log-flow PCA (k = %d)", log_k_chosen), n_params = log_k_chosen,
             type = "in-sample", t(summarise_metrics(logpca_results)),
             negative_points = sum(Y_logpca < 0)),
  data.frame(method = sprintf("Log-flow PCA (k = %d)", log_k_chosen), n_params = log_k_chosen,
             type = "leave-one-out", t(summarise_metrics(logpca_results_loocv)),
             negative_points = NA),
  data.frame(method = sprintf("Wavelets (%d detail)", keep_chosen),
             n_params = ncol(wav_params), type = "in-sample",
             t(summarise_metrics(wav_results)), negative_points = sum(Y_wav < 0)),
  data.frame(method = sprintf("Wavelets (%d detail)", keep_chosen),
             n_params = ncol(wav_params), type = "leave-one-out",
             t(summarise_metrics(wav_results_loocv)), negative_points = NA)
)

print(comparison, digits = 3, row.names = FALSE)
