# -----------------------------
# This script is for reproducing my codes for Shanika sanity checking.
# -----------------------------

library(Matrix)
library(tidyr)
library(ggplot2)

library(vctrs)
library(fabletools)
library(fable)
library(feasts)
library(tsibble)
library(lubridate)
library(dplyr)

# Load data and wrangling -----------------

visnights_raw <- readr::read_csv(
  ... # <your path to the file>
  ) |>
  mutate(Month = yearmonth(Month)) |>
  group_by(Month, Region) |>
  summarise(Nights = sum(Nights), .groups = "drop")

## Define a mapping of regions to states -----------
state_map <- list(
  NSW = c(
    "Blue Mountains", "Central Coast", "Central NSW", "Goulburn",
    "Hunter", "New England North West", "North Coast NSW", "Outback NSW",
    "Riverina", "Snowy Mountains", "South Coast", "Sydney",
    "Capital Country"
  ),
  VIC = c(
    "Ballarat", "Bendigo Loddon", "Central Murray", "Geelong and the Bellarine",
    "Gippsland", "Great Ocean Road", "High Country", "Macedon", "Mallee",
    "Melbourne", "Melbourne East", "Murray East", "Peninsula", "Phillip Island",
    "Spa Country", "The Murray", "Upper Yarra", "Western Grampians", "Wimmera"
  ),
  QLD = c(
    "Brisbane", "Bundaberg", "Capricorn", "Fraser Coast", "Gladstone",
    "Gold Coast", "Mackay", "Outback Queensland", "Southern Queensland Country",
    "Sunshine Coast", "Townsville", "Tropical North Queensland",
    "Whitsundays"
  ),
  SA = c(
    "Adelaide", "Adelaide Hills", "Barossa", "Clare Valley", "Eyre Peninsula",
    "Fleurieu Peninsula", "Flinders Ranges and Outback", "Kangaroo Island",
    "Limestone Coast", "Murray River, Lakes and Coorong", "Riverland",
    "Yorke Peninsula"
  ),
  WA = c(
    "Australia's Coral Coast", "Australia's Golden Outback",
    "Australia's North West", "Australia's South West", "Destination Perth"
  ),
  TAS = c(
    "Central Highlands", "East Coast", "Hobart and the South",
    "Launceston and the North", "Lakes", "North West", "West Coast"
  ),
  NT = c(
    "Alice Springs", "Barkly", "Darwin", "Katherine Daly",
    "Litchfield Kakadu Arnhem", "Lasseter", "MacDonnell"
  ),
  ACT = c(
    "Canberra"
  )
)

region_state_lookup <- tibble(
  State  = rep(names(state_map), lengths(state_map)),
  Region = unlist(state_map)
)

visnights <- visnights_raw %>%
  left_join(region_state_lookup, by = "Region") %>%
  select(Month, State, Region, Nights)

# aggregate
visnights_full <- visnights %>%
  as_tsibble(index = Month, key = c(State, Region)) %>%
  aggregate_key(State/Region, Nights = sum(Nights))


# Modelling --------------------

fit <- visnights_full |>
  model(base = ETS(Nights))

# base forecasts
h <- 12
fc <- fit %>%
  forecast(h = h)


# # # # # # # # # # # #
# Obtain S matrix ------
# # # # # # # # # # # #

key_data <- key_data(fit)
agg_data <- fabletools:::build_key_data_smat(key_data)

S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vec_c(!!!agg_data$agg)))
S[length(agg_data$agg)*(vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L

# set colnames
colnames(S) <- key_data %>%
  filter(!is_aggregated(Region)) %>%
  pull(Region) %>%
  as.character()

# set rownames
row_names <- key_data %>%
  distinct(State, Region) %>%
  mutate(
    name = case_when(
      Region == "<aggregated>" & State == "<aggregated>" ~ "Total",
      Region == "<aggregated>" ~ as.character(State),
      TRUE ~ as.character(Region)
    )
  )
rownames(S) <- row_names$name %>%
  as.character()


# Get data from fit and convert to matrices --------
fit_augment <- fit %>%
  augment() %>%
  filter(.model == "base") %>%
  left_join(row_names, by = c("State", "Region")) %>%
  select(Month, State, Region, name, .fitted, Nights)

y_hat <- fit_augment %>%
  as_tibble() %>%
  select(.fitted, name, Month) %>%
  pivot_wider(names_from = name, values_from = .fitted) %>%
  select(-Month) %>%
  as.matrix()

y <- fit_augment %>%
  as_tibble() %>%
  select(Nights, name, Month) %>%
  pivot_wider(names_from = name, values_from = Nights) %>%
  select(-Month) %>%
  as.matrix()

base_fc <- fc %>%
  as_tibble() %>%
  filter(.model == "base") %>%
  left_join(row_names, by = c("State", "Region")) %>%
  select(name, .mean, Month) %>%
  pivot_wider(names_from = name, values_from = .mean) %>%
  select(-Month) %>%
  as.matrix()

actual <- visnights_full %>%
  filter(year(Month) > year_filter & year(Month) <= year_filter +1) %>%
  as_tibble() %>%
  arrange(State, Region) %>%
  left_join(row_names, by = c("State", "Region")) %>%
  select(name, Nights, Month) %>%
  pivot_wider(names_from = name, values_from = Nights) %>%
  select(-Month) %>%
  as.matrix()



# Check the order of series between y/y_hat/base_fc/actual and S matrix
all.equal(
  rownames(S),
  colnames(y) # same for y_hat, base_fc, actual
)
