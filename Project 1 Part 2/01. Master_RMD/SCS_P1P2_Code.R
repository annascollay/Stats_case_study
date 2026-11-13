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
  mutate(QualityCode = as.numeric(as.character(QualityCode))) %>%
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
q     <- river_data$Discharge
times <- river_data$Timestamp
n_years <- n_distinct(river_data$water_year)

# ------------------------------------------------------------------------------
# 3. Peak Identification (Robson & Reed 1999 Criteria, u = 65 m³/s)
# ------------------------------------------------------------------------------

# Helper: Find rising limb breakout (P1P1 logic)
rise_start_row <- function(p, earliest = times[p] - 7 * 24 * 3600) {
  window <- which(times >= earliest & times <= times[p])
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
    window_mask <- abs(as.numeric(difftime(times[candidates], times[p], units = "days"))) <= 7
    q[p] == max(q[candidates[window_mask]])
  })]
  
  rise_hours <- sapply(clean, function(p) {
    as.numeric(difftime(times[p], times[rise_start_row(p)], units = "hours"))
  })
  min_gap_hours <- 3 * mean(rise_hours)
  
  kept <- candidates[1]
  for (p in candidates[-1]) {
    last <- kept[length(kept)]
    gap <- as.numeric(difftime(times[p], times[last], units = "hours"))
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

# Row number of each peak in river_data (used by the later sections)
peak_rows <- match(peaks$Timestamp, times)

cat("Number of independent events:", nrow(peaks), "\n")

# ------------------------------------------------------------------------------
# 4. Event start, end and width (trough to trough)
# ------------------------------------------------------------------------------

search_days   <- 7     # how far either side of a peak to search for a trough
trough_factor <- 1.05  # trough = reading within 5% of the local minimum

# End of the falling limb: the first reading after the peak within 5% of the
# lowest flow before the limit (mirrors rise_start_row)
fall_end_row <- function(p, limit) {
  window <- which(times >= times[p] & times <= limit)
  lowest <- min(q[window])
  min(window[q[window] <= trough_factor * lowest])
}

# Where to stop searching for each event's end: 7 days after its peak,
# or the next event's peak if that comes sooner
end_limit <- pmin(
  times[peak_rows] + search_days * 24 * 3600,
  lead(times[peak_rows], default = max(times))
)

# One row per event: start, peak and end, with the key quantities for each
events <- tibble(
  event_id  = seq_along(peak_rows),
  peak_row  = peak_rows,
  start_row = map_int(peak_rows, rise_start_row),
  end_row   = map2_int(peak_rows, end_limit, fall_end_row)
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

# ------------------------------------------------------------------------------
# 5. Standardise each event and make the Y matrix
# ------------------------------------------------------------------------------
# For each event, one at a time:
#   a. rescale time onto [0, 1]   (0 = start of event, 1 = end)
#   b. realign, so the peak moves to u_star (the average peak position),
#      keeping the start at 0 and the end at 1
#   c. divide flow by peak flow   (so every event peaks at 1)
#   d. read the curve off at 64 equally spaced points between 0 and 1
#      (linear interpolation between the 15-minute readings), giving one
#      row of Y (Lecture 4: stack the hydrographs on a common grid of times)

u_star <- mean(events$peak_position)
n_grid <- 64
grid_u <- seq(0, 1, length.out = n_grid)

Y <- matrix(NA, nrow = nrow(events), ncol = n_grid)

for (i in 1:nrow(events)) {
  rows <- events$start_row[i]:events$end_row[i]   # readings in this event
  p    <- events$peak_position[i]                 # where its peak falls in [0, 1]
  
  # a. rescale time onto [0, 1]
  u <- as.numeric(difftime(times[rows], events$start_time[i], units = "hours")) /
    events$width_h[i]
  
  # b. realign: stretch the rising part to [0, u_star], the falling part to [u_star, 1]
  u_reg <- ifelse(u <= p,
                  u * u_star / p,
                  u_star + (u - p) * (1 - u_star) / (1 - p))
  
  # c. divide flow by peak flow
  y <- q[rows] / events$peak_flow[i]
  
  # d. read off at the 64 grid points
  Y[i, ] <- approx(u_reg, y, xout = grid_u, rule = 2)$y
}

cat("Y:", nrow(Y), "events x", ncol(Y), "points; peaks aligned at u =", round(u_star, 3), "\n")

# ------------------------------------------------------------------------------
# 6. Reconstruct a real hydrograph
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
# 7. Shared functions for checking every method
# ------------------------------------------------------------------------------
# Every method is checked in the same way, so they can be compared fairly:
# each fitted curve is turned back into a real flood (section 6) and compared
# with the observed flood in m^3/s and hours.

# Area under a curve (trapezium rule), used for flood volume. Needed because
# after undoing the realignment the points are not evenly spaced in hours.
area_under <- function(hours, flow) {
  sum(diff(hours) * (head(flow, -1) + tail(flow, -1)) / 2)
}

# Errors for one flood (observed and fitted flow at the same hours):
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

# Average size of each error across floods (one line per method)
average_errors <- function(errors) {
  colMeans(abs(errors))
}

# Leave-one-out cross-validation (Lecture 5): each flood is left out in turn,
# the method is fitted to the other floods, and the left-out flood is
# predicted. fit_left_out(i) must return the method's curve for flood i
# using only the other floods. (u_star is the average over all floods; its
# effect on the left-out flood is negligible.)
loocv <- function(fit_left_out) {
  Y_pred <- Y
  for (i in 1:nrow(Y)) {
    Y_pred[i, ] <- fit_left_out(i)
  }
  Y_pred
}

# Fixed training and test floods (Lecture 5): about 20% held out, chosen once
# with a fixed seed so everyone uses the same split
set.seed(2026)
test_floods  <- sort(sample(1:nrow(Y), size = round(0.2 * nrow(Y))))
train_floods <- setdiff(1:nrow(Y), test_floods)

# Plot observed vs fitted floods in real units for four kinds of flood:
# earliest peak, latest peak, shortest and longest
plot_example_fits <- function(Y_fit, title) {
  examples <- unique(c(which.min(events$peak_position), which.max(events$peak_position),
                       which.min(events$width_h),       which.max(events$width_h)))
  par(mfrow = c(2, 2))
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

# Plot residuals against time for every flood. Look for patterns (e.g. always
# too low at the peak), not just small residuals. Red line: average residual.
plot_residuals <- function(Y_fit, title) {
  residuals <- Y - Y_fit
  matplot(grid_u, t(residuals), type = "l", lty = 1, col = "grey",
          xlab = "Rescaled time u", ylab = "Observed - fitted",
          main = paste(title, ": residuals against time"))
  lines(grid_u, colMeans(residuals), col = "red", lwd = 2)
  abline(h = 0, lty = 2)
}

# Plot the four average errors against a method's setting (e.g. degree),
# with one line for each type of error (e.g. in-sample and leave-one-out)
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


# ------------------------------------------------------------------------------
# 7a. SEPA's current approach (the baseline)
# ------------------------------------------------------------------------------
# SEPA uses one average shape for every flood, scaled by peak flow. SEPA does
# not realign peaks, so the average is taken over the floods on their own
# timing. Every method is compared against this.

# Each flood's curve on its own timing (not realigned), on the 64 points
unaligned_curve <- function(i) {
  real <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                 events$width_h[i], events$peak_flow[i])
  approx(real$hours / events$width_h[i], Y[i, ], xout = grid_u, rule = 2)$y
}

Y_unaligned <- t(sapply(1:nrow(Y), unaligned_curve))

# Average shape given to flood i, laid out like Y (so it can be checked in
# the same way as every method)
baseline_curve <- function(average_shape, i) {
  real <- reconstruct_hydrograph(Y[i, ], events$peak_position[i],
                                 events$width_h[i], events$peak_flow[i])
  approx(grid_u, average_shape, xout = real$hours / events$width_h[i], rule = 2)$y
}

# In-sample: every flood gets the average of all floods
average_shape <- colMeans(Y_unaligned)
Y_baseline    <- t(sapply(1:nrow(Y), function(i) baseline_curve(average_shape, i)))

baseline_errors <- method_errors(Y_baseline)

# Leave-one-out: every flood gets the average of the other floods
baseline_errors_loocv <- method_errors(loocv(function(i) {
  baseline_curve(colMeans(Y_unaligned[-i, ]), i)
}))

cat("\nSEPA's approach (average shape), average errors:\n")
print(rbind(in_sample     = average_errors(baseline_errors),
            leave_one_out = average_errors(baseline_errors_loocv)))


# ------------------------------------------------------------------------------
# 8. Method 1: polynomial fitting
# ------------------------------------------------------------------------------
# Each flood's curve (one row of Y) is fitted with a polynomial in time,
# separately for each flood, using lm() (Lecture 4):
#   y(u) = theta_0 + theta_1 u + theta_2 u^2 + ... + theta_p u^p + error
# The fitted coefficients theta_0, ..., theta_p are the flood's shape
# parameters. Each flood is fitted on its own, so leave-one-out errors are the
# same as in-sample errors; the degree is checked on training and test floods.

# Formula y ~ I(u^1) + I(u^2) + ... + I(u^p), written as in the handout
poly_formula <- function(p) {
  as.formula(paste("y ~", paste0("I(u^", 1:p, ")", collapse = " + ")))
}

# Fit one flood's curve with a degree-p polynomial
fit_poly_model <- function(y, p) {
  lm(poly_formula(p), data = data.frame(y = y, u = grid_u))
}

# Fitted curves for every flood at degree p
poly_fit_all <- function(p) {
  t(sapply(1:nrow(Y), function(i) fitted(fit_poly_model(Y[i, ], p))))
}


# ------------------------------------------------------------------------------
# 8a. Choosing the degree
# ------------------------------------------------------------------------------

degrees <- 2:10

# Significance of the highest term (as in the handout's quadratic vs cubic
# comparison): proportion of floods where the highest term has p < 0.05.
# Neighbouring points on a curve are correlated, so these p-values are too
# optimistic (Lecture 5); use alongside the errors below.
prop_significant <- sapply(degrees, function(p) {
  mean(sapply(1:nrow(Y), function(i) {
    coefs <- summary(fit_poly_model(Y[i, ], p))$coefficients
    coefs[nrow(coefs), 4] < 0.05
  }))
})

print(data.frame(degree = degrees, prop_significant = prop_significant))

# Errors by degree on training and test floods (Lecture 5): the chosen degree
# should work as well on the test floods as on the training floods
poly_choice <- data.frame()
for (p in degrees) {
  errors <- method_errors(poly_fit_all(p))
  poly_choice <- rbind(poly_choice,
                       data.frame(setting = p, type = "train", t(average_errors(errors[train_floods, ]))),
                       data.frame(setting = p, type = "test",  t(average_errors(errors[test_floods, ]))))
}

print(poly_choice, digits = 3)
plot_choice(poly_choice, "Polynomial degree")


# ------------------------------------------------------------------------------
# 8b. Fit at the chosen degree
# ------------------------------------------------------------------------------
# UPDATE after looking at section 8a

p_chosen <- 6

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


# ------------------------------------------------------------------------------
# 8c. Checking the fit (Lecture 5; the handout)
# ------------------------------------------------------------------------------

plot_example_fits(Y_poly, "Polynomial")
plot_residuals(Y_poly, "Polynomial")

# Residuals against each polynomial term (Lecture 5; the handout)
residuals_poly <- as.vector(t(Y - Y_poly))
u_all          <- rep(grid_u, times = nrow(Y))

par(mfrow = c(2, 3))
for (k in 1:p_chosen) {
  plot(u_all^k, residuals_poly, pch = 20, cex = 0.4, col = "grey40",
       xlab = paste0("u^", k), ylab = "Residual", main = paste0("Residuals vs u^", k))
  abline(h = 0, col = "red", lty = 2)
}
par(mfrow = c(1, 1))

# Standard diagnostic plots for the worst-fitted flood (as in the handout):
# residuals vs fitted, Q-Q, scale-location, leverage. Normality and
# independence matter less here: the aim is a good basis, not inference.
worst <- which.max(poly_errors[, "rmse"])

par(mfrow = c(2, 2))
plot(poly_models[[worst]])
par(mfrow = c(1, 1))

# Coefficient table for the same flood (as in the handout)
print(summary(poly_models[[worst]]))

# Correlation between shape parameters across floods (Lecture 5): strong
# correlations make modelling their uncertainty harder later
cat("\nCorrelation between polynomial shape parameters:\n")
print(round(cor(poly_coefs), 2))


# ------------------------------------------------------------------------------
# 9. Method 2: principal component analysis (PCA)
# ------------------------------------------------------------------------------
# PCA finds the main shapes in the floods from the data itself (Lecture 4).
# Each flood's curve is rebuilt as the mean curve plus k component curves:
#   y_i(u) = mean curve(u) + C_i1 a_1(u) + ... + C_ik a_k(u)
#   a_j  : loadings of component j (a curve over the 64 grid points)
#   C_ij : score of flood i on component j
# The scores C_i1, ..., C_ik are the flood's shape parameters.
# Curves are centred but not scaled, as all points are on the same scale
# (as in the lecture example).

pca_result <- prcomp(Y, center = TRUE, scale. = FALSE)

# Rebuilt curves using the first k components
pca_rebuild <- function(k) {
  scores   <- pca_result$x[, 1:k, drop = FALSE]
  loadings <- pca_result$rotation[, 1:k, drop = FALSE]
  sweep(scores %*% t(loadings), 2, colMeans(Y), FUN = "+")
}

# Left-out flood i, rebuilt from components found without it (Lecture 5):
# centre it with the other floods' mean, project onto their loadings, rebuild
pca_left_out <- function(i, k) {
  pca_train  <- prcomp(Y[-i, ], center = TRUE, scale. = FALSE)
  mean_train <- colMeans(Y[-i, ])
  loadings   <- pca_train$rotation[, 1:k, drop = FALSE]
  scores     <- (Y[i, ] - mean_train) %*% loadings
  as.vector(mean_train + scores %*% t(loadings))
}


# ------------------------------------------------------------------------------
# 9a. Components and variance explained (as in the lecture example)
# ------------------------------------------------------------------------------

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

# Cumulative proportion of variance explained (Lecture 5): choose the
# smallest k that captures enough (Lecture 4 suggests e.g. 95%)
cum_var <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)

plot(seq_along(cum_var), cum_var, type = "b", pch = 19, ylim = c(0, 1),
     xlab = "Number of components k", ylab = "Cumulative proportion of variance",
     main = "Cumulative variance explained")
abline(h = 0.95, lty = 2, col = "red")

cat("Smallest k explaining at least 95% of variance:", which(cum_var >= 0.95)[1], "\n")


# ------------------------------------------------------------------------------
# 9b. Choosing k: in-sample and leave-one-out errors (Lecture 5)
# ------------------------------------------------------------------------------
# In-sample errors always fall as k grows; leave-one-out errors show whether
# extra components still help for a new flood

k_values <- 1:min(10, nrow(Y) - 2)

pca_choice <- data.frame()
for (k in k_values) {
  pca_choice <- rbind(pca_choice,
                      data.frame(setting = k, type = "in-sample",
                                 t(average_errors(method_errors(pca_rebuild(k))))),
                      data.frame(setting = k, type = "leave-one-out",
                                 t(average_errors(method_errors(loocv(function(i) pca_left_out(i, k)))))))
}

print(pca_choice, digits = 3)
plot_choice(pca_choice, "Number of components k")


# ------------------------------------------------------------------------------
# 9c. Fit at the chosen k and check the fit (Lecture 5)
# ------------------------------------------------------------------------------
# UPDATE after looking at sections 9a and 9b

k_chosen <- 3

# Rebuilt curves, one row per flood
Y_pca <- pca_rebuild(k_chosen)

# Shape parameters (scores), one row per flood
pca_scores <- pca_result$x[, 1:k_chosen, drop = FALSE]

pca_errors       <- method_errors(Y_pca)
pca_errors_loocv <- method_errors(loocv(function(i) pca_left_out(i, k_chosen)))

cat("\nPCA (k =", k_chosen, "), average errors:\n")
print(rbind(in_sample     = average_errors(pca_errors),
            leave_one_out = average_errors(pca_errors_loocv)))

# Rebuilt flow should never be negative
cat("Grid points with negative rebuilt flow:", sum(Y_pca < 0), "\n")

plot_example_fits(Y_pca, "PCA")
plot_residuals(Y_pca, "PCA")

# Outliers in the scores (Lecture 5): floods far from the rest are unusual
pairs(pca_scores, pch = 19, main = "Scores on the kept components")

# Variance and correlation of the scores (Lecture 5). Scores are
# uncorrelated by construction, which helps when modelling their uncertainty.
cat("\nVariance of each PCA score:\n")
print(round(apply(pca_scores, 2, var), 4))
cat("\nCorrelation between PCA scores:\n")
print(round(cor(pca_scores), 2))


# ------------------------------------------------------------------------------
# 10. Method 3: wavelets
# ------------------------------------------------------------------------------
# Each flood's curve is split by the discrete wavelet transform into a smooth
# trend plus detail at several scales (Lecture 4). Most detail coefficients
# are tiny, so only the largest few are kept and the rest set to zero.
# The same detail coefficients are kept for every flood (Lecture 4), chosen
# as those that are largest on average across floods.
# As in the lecture example: Daubechies "d4" filter, 3 levels, periodic
# boundaries. Y has 64 = 2^6 points, so no padding is needed.
# Shape parameters: the smooth coefficients plus the kept detail coefficients.

library(wavelets)

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

# Rebuild one curve keeping only the `keep` detail coefficients that are
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

# Rebuilt curves for every flood, keeping `keep` detail coefficients
wavelet_fit_all <- function(keep) {
  t(sapply(wav_all, function(wt) wavelet_rebuild(wt, avg_size, keep)))
}

# Left-out flood i, keeping the coefficients chosen from the other floods
wavelet_left_out <- function(i, keep) {
  wavelet_rebuild(wav_all[[i]], average_detail_size(wav_all[-i]), keep)
}


# ------------------------------------------------------------------------------
# 10a. Decomposition of one flood (as in the lecture example)
# ------------------------------------------------------------------------------

wt_example <- wav_all[[1]]

par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
plot(wt_example@W[[1]], type = "l", col = "red",       xlab = "", ylab = "", main = "Level 1 detail (fine)")
plot(wt_example@W[[2]], type = "l", col = "blue",      xlab = "", ylab = "", main = "Level 2 detail")
plot(wt_example@W[[3]], type = "l", col = "darkgreen", xlab = "", ylab = "", main = "Level 3 detail (coarse)")
plot(wt_example@V[[3]], type = "l", col = "black",     xlab = "", ylab = "", main = "Level 3 smooth (trend)")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Keeping every coefficient should rebuild the curve exactly
cat("Max error of full rebuild:", max(abs(Y[1, ] - as.numeric(idwt(wt_example)))), "\n")


# ------------------------------------------------------------------------------
# 10b. Choosing how many detail coefficients to keep (Lecture 5)
# ------------------------------------------------------------------------------

keep_values <- c(1, 2, 3, 5, 8, 10, 15, 20)

wav_choice <- data.frame()
for (keep in keep_values) {
  wav_choice <- rbind(wav_choice,
                      data.frame(setting = keep, type = "in-sample",
                                 t(average_errors(method_errors(wavelet_fit_all(keep))))),
                      data.frame(setting = keep, type = "leave-one-out",
                                 t(average_errors(method_errors(loocv(function(i) wavelet_left_out(i, keep)))))))
}

print(wav_choice, digits = 3)
plot_choice(wav_choice, "Detail coefficients kept")


# ------------------------------------------------------------------------------
# 10c. Fit at the chosen number and check the fit (Lecture 5)
# ------------------------------------------------------------------------------
# UPDATE after looking at section 10b

keep_chosen <- 5

# Rebuilt curves, one row per flood
Y_wav <- wavelet_fit_all(keep_chosen)

# Shape parameters: smooth coefficients + kept detail coefficients
cutoff     <- sort(unlist(avg_size), decreasing = TRUE)[keep_chosen]
kept       <- unlist(avg_size) >= cutoff
wav_params <- t(sapply(wav_all, function(wt) {
  c(wt@V[[wav_levels]][, 1], unlist(wt@W)[kept])
}))

cat("Wavelet shape parameters per flood:", ncol(wav_params), "\n")

wav_errors       <- method_errors(Y_wav)
wav_errors_loocv <- method_errors(loocv(function(i) wavelet_left_out(i, keep_chosen)))

cat("\nWavelets (", keep_chosen, "detail coefficients ), average errors:\n")
print(rbind(in_sample     = average_errors(wav_errors),
            leave_one_out = average_errors(wav_errors_loocv)))

# Rebuilt flow should never be negative
cat("Grid points with negative rebuilt flow:", sum(Y_wav < 0), "\n")

plot_example_fits(Y_wav, "Wavelets")
plot_residuals(Y_wav, "Wavelets")

# Sensitivity to the wavelet filter and number of levels (Lecture 5), keeping
# the same number of detail coefficients. More levels means fewer smooth
# coefficients (64 / 2^levels), so fewer parameters.
wav_sensitivity <- data.frame()
for (filter in c("haar", "d4", "la8")) {
  for (levels in 3:5) {
    wav_list <- lapply(1:nrow(Y), function(i) wavelet_transform(Y[i, ], filter, levels))
    sizes    <- average_detail_size(wav_list)
    Y_fit    <- t(sapply(wav_list, function(wt) wavelet_rebuild(wt, sizes, keep_chosen)))
    wav_sensitivity <- rbind(wav_sensitivity,
                             data.frame(filter = filter, levels = levels,
                                        n_params = 64 / 2^levels + keep_chosen,
                                        t(average_errors(method_errors(Y_fit)))))
  }
}

print(wav_sensitivity, digits = 3)

# Variance and correlation of the shape parameters (Lecture 5)
cat("\nVariance of each wavelet shape parameter:\n")
print(round(apply(wav_params, 2, var), 4))
cat("\nCorrelation between wavelet shape parameters:\n")
print(round(cor(wav_params), 2))


# ------------------------------------------------------------------------------
# 11. Comparing the methods (Lecture 5)
# ------------------------------------------------------------------------------
# Accuracy (average errors), parsimony (number of shape parameters) and the
# physical check (negative flows). For the polynomial, leave-one-out errors
# are the same as in-sample errors, as each flood is fitted on its own.

comparison <- rbind(
  data.frame(method = "SEPA (average shape)", n_params = 0, type = "in-sample",
             t(average_errors(baseline_errors)), negative_points = sum(Y_baseline < 0)),
  data.frame(method = "SEPA (average shape)", n_params = 0, type = "leave-one-out",
             t(average_errors(baseline_errors_loocv)), negative_points = NA),
  data.frame(method = "Polynomial", n_params = p_chosen + 1, type = "in-sample",
             t(average_errors(poly_errors)), negative_points = sum(Y_poly < 0)),
  data.frame(method = "PCA", n_params = k_chosen, type = "in-sample",
             t(average_errors(pca_errors)), negative_points = sum(Y_pca < 0)),
  data.frame(method = "PCA", n_params = k_chosen, type = "leave-one-out",
             t(average_errors(pca_errors_loocv)), negative_points = NA),
  data.frame(method = "Wavelets", n_params = ncol(wav_params), type = "in-sample",
             t(average_errors(wav_errors)), negative_points = sum(Y_wav < 0)),
  data.frame(method = "Wavelets", n_params = ncol(wav_params), type = "leave-one-out",
             t(average_errors(wav_errors_loocv)), negative_points = NA)
)

print(comparison, digits = 3, row.names = FALSE)

# ------------------------------------------------------------------------------
# 12. Method 4: known hydrograph function (gamma-type)
# ------------------------------------------------------------------------------
# Lecture 4 suggests fitting a known function that approximates hydrographs,
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

# Fit one flood's curve; returns NULL if the fit fails to converge
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

# Fitted curves, one row per flood (failed fits left as the observed curve,
# and excluded from the errors below)
Y_gamma <- Y
for (i in gamma_ok) {
  Y_gamma[i, ] <- fitted(gamma_models[[i]])
}

# Shape parameters, one row per flood (NA where the fit failed)
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
# 13. Method 5: B-splines
# ------------------------------------------------------------------------------
# Ramsay & Silverman's standard alternative to polynomials. A B-spline basis
# joins short cubic pieces smoothly, so each basis function only affects one
# part of the curve. This avoids the polynomial's wiggles in the tail, while
# still being a linear model fitted with lm(), like the polynomial:
#   y(u) = theta_0 + theta_1 B_1(u) + ... + theta_df B_df(u) + error
# df (the number of basis functions) plays the role of the polynomial degree.
# Each flood is fitted on its own, so the number of basis functions is
# checked on training and test floods, as for the polynomial.

library(splines)   # part of base R; provides bs()

# Fit one flood's curve with df B-spline basis functions
fit_bspline_model <- function(y, df) {
  lm(y ~ bs(u, df = df), data = data.frame(y = y, u = grid_u))
}

# Fitted curves for every flood with df basis functions
bspline_fit_all <- function(df) {
  t(sapply(1:nrow(Y), function(i) fitted(fit_bspline_model(Y[i, ], df))))
}


# ------------------------------------------------------------------------------
# 13a. Choosing the number of basis functions
# ------------------------------------------------------------------------------

df_values <- 4:15

bspline_choice <- data.frame()
for (df in df_values) {
  errors <- method_errors(bspline_fit_all(df))
  bspline_choice <- rbind(bspline_choice,
                          data.frame(setting = df, type = "train", t(average_errors(errors[train_floods, ]))),
                          data.frame(setting = df, type = "test",  t(average_errors(errors[test_floods, ]))))
}

print(bspline_choice, digits = 3)
plot_choice(bspline_choice, "Number of B-spline basis functions")


# ------------------------------------------------------------------------------
# 13b. Fit at the chosen number and check the fit
# ------------------------------------------------------------------------------
# UPDATE after looking at section 13a

df_chosen <- 8

bspline_models <- lapply(1:nrow(Y), function(i) fit_bspline_model(Y[i, ], df_chosen))

# Fitted curves, one row per flood
Y_bspline <- t(sapply(bspline_models, fitted))

# Shape parameters, one row per flood
bspline_coefs <- t(sapply(bspline_models, coef))
colnames(bspline_coefs) <- paste0("theta_", 0:df_chosen)

bspline_errors <- method_errors(Y_bspline)

cat("\nB-splines (", df_chosen, "basis functions ), average errors:\n")
print(average_errors(bspline_errors))

# Fitted flow should never be negative
cat("Grid points with negative fitted flow:", sum(Y_bspline < 0), "\n")

plot_example_fits(Y_bspline, "B-spline")
plot_residuals(Y_bspline, "B-spline")

# Standard diagnostic plots for the worst-fitted flood (as in the handout)
worst <- which.max(bspline_errors[, "rmse"])

par(mfrow = c(2, 2))
plot(bspline_models[[worst]])
par(mfrow = c(1, 1))

# Correlation between shape parameters across floods
cat("\nCorrelation between B-spline shape parameters:\n")
print(round(cor(bspline_coefs), 2))


# ------------------------------------------------------------------------------
# 14. Method 6: PCA on log flow
# ------------------------------------------------------------------------------
# Ramsay & Silverman represent curves that must stay positive through their
# logarithm. PCA (as in section 9) is fitted to log(y), and rebuilt curves
# are exponentiated, so rebuilt (and later simulated) flows can never be
# negative.

Y_log <- log(Y)   # scaled flows are always above zero

log_pca_result <- prcomp(Y_log, center = TRUE, scale. = FALSE)

# Rebuilt curves using the first k components
log_pca_rebuild <- function(k) {
  scores   <- log_pca_result$x[, 1:k, drop = FALSE]
  loadings <- log_pca_result$rotation[, 1:k, drop = FALSE]
  exp(sweep(scores %*% t(loadings), 2, colMeans(Y_log), FUN = "+"))
}

# Left-out flood i, rebuilt from components found without it
log_pca_left_out <- function(i, k) {
  pca_train  <- prcomp(Y_log[-i, ], center = TRUE, scale. = FALSE)
  mean_train <- colMeans(Y_log[-i, ])
  loadings   <- pca_train$rotation[, 1:k, drop = FALSE]
  scores     <- (Y_log[i, ] - mean_train) %*% loadings
  exp(as.vector(mean_train + scores %*% t(loadings)))
}


# ------------------------------------------------------------------------------
# 14a. Choosing k: in-sample and leave-one-out errors
# ------------------------------------------------------------------------------

log_pca_choice <- data.frame()
for (k in k_values) {
  log_pca_choice <- rbind(log_pca_choice,
                          data.frame(setting = k, type = "in-sample",
                                     t(average_errors(method_errors(log_pca_rebuild(k))))),
                          data.frame(setting = k, type = "leave-one-out",
                                     t(average_errors(method_errors(loocv(function(i) log_pca_left_out(i, k)))))))
}

print(log_pca_choice, digits = 3)
plot_choice(log_pca_choice, "Number of components k (log flow)")


# ------------------------------------------------------------------------------
# 14b. Fit at the chosen k and check the fit
# ------------------------------------------------------------------------------
# UPDATE after looking at section 14a

log_k_chosen <- 3

# Rebuilt curves, one row per flood
Y_log_pca <- log_pca_rebuild(log_k_chosen)

# Shape parameters (scores), one row per flood
log_pca_scores <- log_pca_result$x[, 1:log_k_chosen, drop = FALSE]

log_pca_errors       <- method_errors(Y_log_pca)
log_pca_errors_loocv <- method_errors(loocv(function(i) log_pca_left_out(i, log_k_chosen)))

cat("\nPCA on log flow (k =", log_k_chosen, "), average errors:\n")
print(rbind(in_sample     = average_errors(log_pca_errors),
            leave_one_out = average_errors(log_pca_errors_loocv)))

plot_example_fits(Y_log_pca, "Log-flow PCA")
plot_residuals(Y_log_pca, "Log-flow PCA")

# Variance of each score (scores are uncorrelated by construction)
cat("\nVariance of each log-flow PCA score:\n")
print(round(apply(log_pca_scores, 2, var), 4))


# ------------------------------------------------------------------------------
# 15. Comparing all methods (Lecture 5)
# ------------------------------------------------------------------------------
# Accuracy (average errors), parsimony (number of shape parameters) and the
# physical check (negative flows). For methods fitted to each flood on its
# own (polynomial, gamma-type, B-spline), leave-one-out errors are the same
# as in-sample errors.

all_comparison <- rbind(
  comparison,
  data.frame(method = "Gamma-type", n_params = 3, type = "in-sample",
             t(average_errors(gamma_errors)),
             negative_points = sum(Y_gamma[gamma_ok, ] < 0)),
  data.frame(method = "B-spline", n_params = df_chosen + 1, type = "in-sample",
             t(average_errors(bspline_errors)), negative_points = sum(Y_bspline < 0)),
  data.frame(method = "Log-flow PCA", n_params = log_k_chosen, type = "in-sample",
             t(average_errors(log_pca_errors)), negative_points = sum(Y_log_pca < 0)),
  data.frame(method = "Log-flow PCA", n_params = log_k_chosen, type = "leave-one-out",
             t(average_errors(log_pca_errors_loocv)), negative_points = NA)
)

print(all_comparison, digits = 3, row.names = FALSE)

# ------------------------------------------------------------------------------
# Next steps
#
# - Compare models and choose (looking like PCA?)
# - Model uncertainty in its parameters. I believe by fitting a distribution for
# shape parameters, peak flow, width, and peak position. Check how they depend on
# each other
# - Simulate hydrographs and validate them using this
# - Show the emulator workflow (don't actually need to fit)
# ------------------------------------------------------------------------------