# Load libraries
library(tidyverse)

# Fetch 10 years of 15-minute flow data for River Don at Parkhill via SEPA
# KiWIS API
# Station parameters and time range: Aug-2016 to Aug-2026
url <- paste0(
  "https://timeseries.sepa.org.uk/KiWIS/KiWIS?",
  "service=kisters&type=queryServices&datasource=0&",
  "request=getTimeseriesValues&ts_path=1/234292/Q/15m.Cmd&",
  "from=2016-08-31T00:00:00&to=2026-08-31T00:00:00&",
  "metadata=true&returnfields=Timestamp,%20Value,%20Quality%20Code&format=csv"
)

# Download and parse CSV
river_data <- read_delim(url, delim = ";", comment = "#",
                         col_names = c("Timestamp",
                                       "Discharge",
                                       "QualityCode")) %>%
  mutate(Timestamp = as.POSIXct(Timestamp, format="%Y-%m-%dT%H:%M:%S")) %>%
  filter(!is.na(Discharge))

# Summary & Hydrograph Plot
summary(river_data$Discharge)

ggplot(river_data, aes(x = Timestamp, y = Discharge)) +
  geom_line(color = "steelblue", linewidth = 0.4) +
  labs(
    title = "Hydrograph: River Don at Parkhill (15-min interval)",
    x = "Date",
    y = expression("Discharge ("*m^3*"/s)")
  ) +
  theme_minimal()