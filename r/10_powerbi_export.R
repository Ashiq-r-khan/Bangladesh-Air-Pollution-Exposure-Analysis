# 10 Power BI export: dimension, fact, statistics and model tables, plus relationships, dictionary, glossary and findings
# Bangladesh Air Pollution Project | Md. Ashiqur Rahman Khan

library(arrow)
library(dplyr)
library(tidyr)
library(sf)
library(terra)
setwd("E:/Bangladesh Air Pollution Project")
options(scipen = 0)
out <- "outputs/powerbi"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tabs <- list()
put <- function(d, name) {d <- as_tibble(d) |> mutate(across(where(\(x) is.double(x) && !inherits(x, "Date")), \(x) if (any(abs(x) >= 1e5, na.rm = TRUE)) round(x) else signif(x, 6))); tabs[[name]] <<- d
  write.csv(d, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "", fileEncoding = "UTF-8"); cat(sprintf("%-28s %7d rows  %2d columns\n", name, nrow(d), ncol(d)))}
rt <- \(f) read.csv(file.path("outputs/R/tables", f), check.names = FALSE)
gpkg <- "data/processed/bgd_admin.gpkg"
lv <- c("national", "division", "district", "upazila")

# ---- 1. Dimensions ----
div <- st_read(gpkg, layer = "division", quiet = TRUE) |> st_drop_geometry()
dis <- st_read(gpkg, layer = "district", quiet = TRUE) |> st_drop_geometry()
upz <- st_read(gpkg, layer = "upazila", quiet = TRUE)
uy <- read_parquet("data/processed/pm25_unit_year.parquet")
pop <- uy |> filter(year == 2024) |> select(pcode, people = pop)
dim_geo <- bind_rows(
  tibble(pcode = "BD", level = "national", name = "Bangladesh", district_pcode = NA, division_pcode = NA, area_km2 = sum(div$area_km2)),
  div |> transmute(pcode, level = "division", name, district_pcode = NA, division_pcode = pcode, area_km2),
  dis |> transmute(pcode, level = "district", name, district_pcode = pcode, division_pcode, area_km2),
  st_drop_geometry(upz) |> transmute(pcode, level = "upazila", name, district_pcode, division_pcode, area_km2)) |>
  left_join(div |> select(division_pcode = pcode, division = name), by = "division_pcode") |>
  left_join(dis |> select(district_pcode = pcode, district = name), by = "district_pcode") |>
  group_by(level, name) |> mutate(label = if (n() > 1) paste0(name, " (", district, ")") else name) |> ungroup() |>
  left_join(pop, by = "pcode") |> left_join(st_drop_geometry(upz) |> select(pcode) |> bind_cols(as_tibble(st_coordinates(st_transform(st_point_on_surface(st_geometry(st_transform(upz, 32646))), 4326))) |> rename(lon = X, lat = Y)), by = "pcode") |>
  mutate(level_sort = match(level, lv), density_per_km2 = people / area_km2) |> arrange(level_sort, division, district, name) |>
  mutate(name_sort = row_number()) |> relocate(pcode, level, level_sort, name, label, name_sort)
stopifnot(!anyDuplicated(dim_geo$pcode), !anyDuplicated(paste(dim_geo$level, dim_geo$label)), nrow(dim_geo) == 580)
put(dim_geo, "dim_geo")
put(tibble(year = 1998:2024, period = if_else(year <= 2010, "1998-2010", "2011-2024"), five_year = paste0(year - (year - 1998) %% 5, "-", pmin(year - (year - 1998) %% 5 + 4, 2024)),
           is_2018_spike = year == 2018, year_label = as.character(year)), "dim_year")
put(tibble(month = 1:12, month_name = month.name, month_short = month.abb,
           season = c("Winter", "Winter", "Pre-monsoon", "Pre-monsoon", "Pre-monsoon", "Monsoon", "Monsoon", "Monsoon", "Monsoon", "Post-monsoon", "Post-monsoon", "Winter"),
           dry_or_monsoon = if_else(month %in% 6:9, "Monsoon (Jun-Sep)", if_else(month %in% c(11, 12, 1, 2), "Dry (Nov-Feb)", "Transition"))), "dim_month")
bands <- c(0, 5, 10, 15, 25, 35, 50, 75, 100, Inf)
dim_band <- tibble(band = c("<=5 guideline", "5-10 (IT-4)", "10-15 (IT-3)", "15-25 (IT-2)", "25-35 (IT-1)", "35-50", "50-75", "75-100", ">100"), band_sort = 1:9,
                   lower = head(bands, -1), upper = tail(bands, -1), meets = c("WHO guideline", "IT-4", "IT-3", "IT-2", "IT-1 and national standard", rep("None", 4)),
                   above_35 = band_sort >= 6, color = c("#1a9850", "#66bd63", "#a6d96a", "#d9ef8b", "#fee08b", "#fdae61", "#f46d43", "#d73027", "#7f0000"))
put(dim_band, "dim_band")
preds <- c(log_pop_density = "Population density (log)", urban_share = "Urban share", dist_dhaka_km = "Distance to Dhaka", log_elevation = "Elevation (log)", elevation_m = "Elevation",
           log_roads = "Major road density (log)", road_km_per_km2 = "Major road density", dry_wind_ms = "Dry-season wind speed", rain_mm_day = "Rainfall")
put(tibble(predictor = names(preds), predictor_label = unname(preds), predictor_sort = seq_along(preds),
           unit = c("log people per km2", "share 0-1", "km", "log m", "m", "log km per km2", "km per km2", "m/s", "mm per day")), "dim_predictor")

# ---- 2. Facts: yearly and monthly PM2.5 for every unit, STL components ----
put(uy |> transmute(pcode, year, pm25_popw, pm25_area, weighting_gap = pm25_popw - pm25_area, people = pop), "fact_pm25_year")
mo <- read_parquet("data/processed/pm25_unit_month.parquet") |> mutate(date = as.Date(date))
put(mo |> transmute(pcode, date, year = as.integer(format(date, "%Y")), month = as.integer(format(date, "%m")), pm25_popw), "fact_pm25_month")
put(mo |> filter(pcode %in% dim_geo$pcode[dim_geo$level %in% lv[1:2]]) |> arrange(date) |> group_by(pcode) |> group_modify(\(d, k) {
  f <- stl(ts(log(d$pm25_popw), start = c(1998, 1), frequency = 12), s.window = 13, robust = TRUE)$time.series
  tibble(date = d$date, observed = d$pm25_popw, trend = exp(f[, "trend"]), seasonal_factor = exp(f[, "seasonal"]), remainder_log = f[, "remainder"])}) |> ungroup(), "fact_stl_month")

# ---- 3. Facts: fine exposure bands, percentiles and recent exposure for every unit, from the 0.01 deg grid ----
grid <- sds("data/processed/pm25_annual_grid.nc")
pm <- grid[["pm25"]]
names(pm) <- 1998:2024
cells <- as.data.frame(c(grid[["pop_2022"]], pm), xy = TRUE) |> rename(pop = pop_2022) |> filter(pop > 0)
pts <- st_as_sf(cells[c("x", "y")], coords = c("x", "y"), crs = 4326) |> st_transform(32646)
upz_m <- st_transform(upz, 32646)
idx <- sapply(st_intersects(pts, upz_m), \(i) if (length(i)) i[1] else NA_integer_)
idx[is.na(idx)] <- st_nearest_feature(pts[is.na(idx), ], upz_m)
cells <- cells |> mutate(upazila = upz$pcode[idx], district = upz$district_pcode[idx], division = upz$division_pcode[idx], national = "BD", recent = rowMeans(across(all_of(as.character(2020:2024)))))
long <- cells |> select(pop, all_of(lv), all_of(as.character(1998:2024))) |> pivot_longer(all_of(as.character(1998:2024)), names_to = "year", values_to = "pm25") |>
  mutate(year = as.integer(year), band = as.character(cut(pm25, bands, dim_band$band)))
by_levels <- function(d, f) bind_rows(lapply(lv, \(l) d |> rename(pcode = all_of(l)) |> f()))
wq <- function(x, w, p) {o <- order(x); cw <- cumsum(w[o]) / sum(w); sapply(p, \(q) x[o][which(cw >= q)[1]])}
band_year <- by_levels(long, \(d) d |> group_by(pcode, year, band) |> summarise(people = sum(pop), .groups = "drop")) |>
  complete(nesting(pcode), year, band = dim_band$band, fill = list(people = 0)) |> group_by(pcode, year) |> mutate(share = people / sum(people)) |> ungroup() |>
  left_join(pop, by = "pcode", suffix = c("", "_unit")) |> mutate(people = share * people_unit) |> select(-people_unit)
put(band_year, "fact_band_year")
put(by_levels(long, \(d) d |> group_by(pcode, year) |> summarise(tibble::as_tibble_row(setNames(wq(pm25, pop, c(0.1, 0.25, 0.5, 0.75, 0.9)), c("p10", "p25", "p50", "p75", "p90"))),
                                                                  min = min(pm25), max = max(pm25), .groups = "drop")), "stat_exposure_percentiles")
exposure <- by_levels(cells, \(d) d |> group_by(pcode) |> summarise(pm25_2020_2024 = weighted.mean(recent, pop),
  cut_to_35_pct = 100 * weighted.mean(pmax(0, 1 - 35 / recent), pop), cut_to_25_pct = 100 * weighted.mean(pmax(0, 1 - 25 / recent), pop), cut_to_5_pct = 100 * weighted.mean(pmax(0, 1 - 5 / recent), pop),
  years_lost_who = weighted.mean(0.098 * pmax(0, recent - 5), pop), years_gain_35 = weighted.mean(0.098 * pmax(0, recent - 35), pop), years_gain_25 = weighted.mean(0.098 * pmax(0, recent - 25), pop),
  life_years_lost_m = sum(pop * 0.098 * pmax(0, recent - 5)) / 1e6, .groups = "drop"))
put(exposure, "stat_exposure_recent")

# ---- 4. Statistics from scripts 06 to 08, keyed to pcode ----
key_of <- \(d) d |> left_join(dim_geo |> select(level, name, pcode) |> filter(level %in% lv[1:2]), by = c("level", "name"))
put(rt("06_trend_all_units.csv") |> select(-level, -name) |>
      left_join(key_of(rt("06_trend_periods.csv")) |> select(pcode, sen_1998_2010, sen_2011_2024, p_1998_2010, p_2011_2024), by = "pcode") |>
      left_join(key_of(rt("06_seasonal_mann_kendall.csv") |> select(-pcode)) |> select(pcode, smk_z, smk_p, sea_sen), by = "pcode") |>
      mutate(significant_rise = q_hr < 0.05 & sen > 0), "stat_trend")
put(rt("06_month_slopes.csv") |> select(pcode, month, base_1998_2002 = base, mean_2020_2024 = last5, sen, p, pct_decade), "stat_month_slope")
put(key_of(rt("06_stl_summary.csv")) |> select(-level, -name) |> relocate(pcode), "stat_stl_summary")
put(key_of(rt("06_stl_anomalies.csv")) |> transmute(pcode, date = as.Date(date), observed, trend, remainder_log = remainder, robust_z = rz, direction), "stat_stl_anomaly")
put(rt("07_upazila_clusters.csv") |> select(pcode, starts_with("lisa"), starts_with("gi")), "stat_spatial_cluster")
put(bind_rows(rt("07_global_moran.csv") |> mutate(variable = var, test = "Global Moran's I") |> select(-var),
              rt("07_ml_residual_moran.csv") |> mutate(variable = paste("Residuals", target, cv, model), test = "Moran's I of ML residuals", weights = "queen") |> select(-target, -cv, -model)) |>
      relocate(test, variable, weights), "stat_moran")
put(rt("07_correlogram.csv") |> rename(variable = var, lag_order = order, mean_distance_km = km), "stat_correlogram")
put(rt("08_model_comparison.csv"), "model_spatial_comparison")
put(rt("08_effects.csv") |> rename(predictor = var), "model_spatial_effects")
put(read.csv("outputs/models/05_cv_results.csv") |> rename(r2 = R2, mae = MAE, rmse = RMSE), "model_ml_cv")
shap <- read.csv("outputs/models/05_shap_values.csv") |> pivot_longer(-pcode, names_to = "predictor", values_to = "shap")
put(shap |> mutate(target = "pm25_recent"), "model_shap_upazila")
put(shap |> group_by(predictor) |> summarise(mean_abs_shap = mean(abs(shap))) |> arrange(desc(mean_abs_shap)) |> mutate(rank = row_number(), target = "pm25_recent"), "model_shap_importance")

# ---- 5. Findings, glossary, relationships, dictionary ----
put(bind_rows(read.csv("outputs/tables/04_key_numbers.csv") |> mutate(source = "Notebook 04", target = NA),
              lapply(c("06", "07", "08", "09"), \(s) rt(paste0(s, "_key_numbers.csv")) |> mutate(source = paste("R script", s)))) |>
      mutate(finding_id = row_number(), section = case_when(source == "Notebook 04" ~ "Exposure overview", source == "R script 06" ~ "Trends and seasonality",
             source == "R script 07" ~ "Spatial clustering", source == "R script 08" ~ "Spatial regression", TRUE ~ "WHO exposure and health")) |>
      relocate(finding_id, section, source, target), "findings")
put(tribble(~Term, ~Meaning, ~`How calculated`, ~`How to read`, ~Source,
  "Population-weighted PM2.5", "The annual fine particle level the average person breathes", "Mean of 0.01 degree cells weighted by WorldPop 2022 people", "Higher is worse; compare with WHO 5 and national 35", "ACAG SatPM2.5 V6.GL.03, WorldPop",
  "Area-weighted PM2.5", "The average over land, ignoring where people live", "Mean of cells weighted by cell area", "Gap to the population-weighted value shows where people sit", "ACAG SatPM2.5 V6.GL.03",
  "Weighting gap", "Population-weighted minus area-weighted PM2.5", "pm25_popw - pm25_area", "Positive: people live in the dirtier parts of the area", "Notebook 02",
  "WHO guideline", "WHO 2021 annual PM2.5 guideline", "5 ug/m3", "Nobody in Bangladesh meets it", "WHO 2021",
  "Interim targets", "WHO stepping stones IT-4 to IT-1", "10, 15, 25 and 35 ug/m3", "IT-1 (35) is also Bangladesh's annual standard", "WHO 2021",
  "National standard", "Bangladesh annual PM2.5 standard", "35 ug/m3", "Share above it is the main exceedance KPI", "Air Pollution (Control) Rules 2022",
  "Sen's slope", "Typical yearly change in PM2.5", "Median of all pairwise slopes, 1998-2024", "ug/m3 per year; robust to the 2018 spike", "R script 06",
  "Mann-Kendall test", "Test for a steady upward or downward trend", "Hamed-Rao corrected for autocorrelation, BH-adjusted across units", "q < 0.05: significant trend", "R script 06",
  "Relative rise", "Trend as a share of the starting level", "10 x Sen's slope / 1998-2002 mean", "% per decade", "R script 06",
  "Seasonal swing", "Worst month divided by cleanest month", "exp(max - min) of the STL seasonal component", "Lower: winter and monsoon closer together", "R script 06",
  "Moran's I", "How much an area resembles its neighbours", "Queen contiguity, row standardised, 999 permutations", "0 random, 1 perfectly clustered", "R script 07",
  "LISA cluster", "Local cluster type for one upazila", "Local Moran's I, 999 permutations, p < 0.05", "High-High: hot cluster; Low-Low: clean cluster", "R script 07",
  "Gi* hot spot", "Area where high values concentrate", "Getis-Ord Gi* z-score with the upazila and its neighbours", "90, 95, 99% confidence classes", "R script 07",
  "Spatial error model", "Regression that allows errors to cluster in space", "spatialreg::errorsarlm, queen weights", "Effect per 1 SD of a predictor", "R script 08",
  "Cut needed", "Reduction needed to reach a target", "1 - target / PM2.5, averaged over people, 2020-2024", "% of today's level", "R script 09",
  "Years of life lost (AQLI)", "Life expectancy lost relative to the WHO guideline", "0.098 years per ug/m3 above 5, averaged over people, 2020-2024", "An estimate from a published relationship, not a measured outcome", "AQLI, Ebenstein et al. 2017",
  "Spatial cross-validation", "Model test on whole regions held out", "Folds by region, notebook 05", "Lower than random CV when models lean on neighbours", "Notebook 05"), "glossary")
rel <- tribble(~from_table, ~from_column, ~to_table, ~to_column,
  "fact_pm25_year", "pcode", "dim_geo", "pcode", "fact_pm25_year", "year", "dim_year", "year",
  "fact_pm25_month", "pcode", "dim_geo", "pcode", "fact_pm25_month", "year", "dim_year", "year", "fact_pm25_month", "month", "dim_month", "month",
  "fact_stl_month", "pcode", "dim_geo", "pcode",
  "fact_band_year", "pcode", "dim_geo", "pcode", "fact_band_year", "year", "dim_year", "year", "fact_band_year", "band", "dim_band", "band",
  "stat_exposure_percentiles", "pcode", "dim_geo", "pcode", "stat_exposure_percentiles", "year", "dim_year", "year",
  "stat_exposure_recent", "pcode", "dim_geo", "pcode", "stat_trend", "pcode", "dim_geo", "pcode",
  "stat_month_slope", "pcode", "dim_geo", "pcode", "stat_month_slope", "month", "dim_month", "month",
  "stat_stl_summary", "pcode", "dim_geo", "pcode", "stat_stl_anomaly", "pcode", "dim_geo", "pcode", "stat_spatial_cluster", "pcode", "dim_geo", "pcode",
  "model_spatial_effects", "predictor", "dim_predictor", "predictor", "model_shap_upazila", "pcode", "dim_geo", "pcode",
  "model_shap_upazila", "predictor", "dim_predictor", "predictor", "model_shap_importance", "predictor", "dim_predictor", "predictor") |>
  mutate(cardinality = "many-to-one", cross_filter = "single")
for (i in seq_len(nrow(rel))) with(rel[i, ], {f <- tabs[[from_table]][[from_column]]; k <- tabs[[to_table]][[to_column]]
  if (anyDuplicated(k) || !all(na.omit(f) %in% k)) stop(sprintf("Key check failed: %s.%s -> %s.%s", from_table, from_column, to_table, to_column))})
put(rel, "relationships")
desc <- c(pcode = "HDX pcode of the area (BD = national)", level = "Admin level: national, division, district or upazila", level_sort = "Sort order of level", name = "Area name", label = "Unique display name (district added to repeated upazila names)",
  name_sort = "Sort order: level, division, district, name", district_pcode = "Parent district pcode", division_pcode = "Parent division pcode", area_km2 = "Area in km2", division = "Division name", district = "District name",
  people = "People, WorldPop 2022", lon = "Longitude of a point inside the upazila", lat = "Latitude of a point inside the upazila", density_per_km2 = "People per km2",
  year = "Calendar year", period = "1998-2010 or 2011-2024", five_year = "Five-year group", is_2018_spike = "TRUE for 2018, the spike year in V6.GL.03", year_label = "Year as text",
  month = "Month number", month_name = "Month name", month_short = "Short month name", season = "Bangladesh season", dry_or_monsoon = "Dry (Nov-Feb), monsoon (Jun-Sep) or transition",
  band = "Annual PM2.5 band", band_sort = "Band sort order", lower = "Lower bound, ug/m3", upper = "Upper bound, ug/m3", meets = "Highest WHO level met in this band", above_35 = "TRUE if above the national standard", color = "Hex colour for the band",
  predictor = "Predictor code", predictor_label = "Predictor name", predictor_sort = "Predictor sort order", unit = "Predictor unit",
  pm25_popw = "Population-weighted PM2.5, ug/m3", pm25_area = "Area-weighted PM2.5, ug/m3", weighting_gap = "pm25_popw minus pm25_area", date = "First day of the month",
  observed = "Observed PM2.5, ug/m3", trend = "STL trend, ug/m3", seasonal_factor = "STL seasonal factor (1 = average month)", remainder_log = "STL remainder on the log scale",
  share = "Share of the area's people", p10 = "10th percentile of exposure across people", p25 = "25th percentile", p50 = "Median person's exposure", p75 = "75th percentile", p90 = "90th percentile", min = "Lowest cell", max = "Highest cell",
  pm25_2020_2024 = "Population-weighted PM2.5, 2020-2024 mean", cut_to_35_pct = "Average % cut needed to reach 35", cut_to_25_pct = "Average % cut needed to reach 25", cut_to_5_pct = "Average % cut needed to reach 5",
  years_lost_who = "Life expectancy lost vs the WHO guideline, years (AQLI)", years_gain_35 = "Years gained by meeting 35", years_gain_25 = "Years gained by meeting 25", life_years_lost_m = "Total life-years lost, millions",
  tau = "Kendall tau", z = "Mann-Kendall z", p = "p-value", z_hr = "Hamed-Rao corrected z", p_hr = "Hamed-Rao corrected p-value", sen = "Sen's slope, ug/m3 per year", sen_lo = "Sen's slope 95% CI low", sen_hi = "Sen's slope 95% CI high",
  sen_no2018 = "Sen's slope without 2018", fit0 = "Sen line value in 1998", first5 = "1998-2002 mean", last5 = "2020-2024 mean", pct_decade = "Rise as % of the 1998-2002 level per decade", q_hr = "BH-adjusted Hamed-Rao p-value",
  sen_1998_2010 = "Sen's slope 1998-2010", sen_2011_2024 = "Sen's slope 2011-2024", p_1998_2010 = "p-value 1998-2010", p_2011_2024 = "p-value 2011-2024", smk_z = "Seasonal Mann-Kendall z", smk_p = "Seasonal Mann-Kendall p-value",
  sea_sen = "Seasonal Sen's slope, ug/m3 per year", significant_rise = "TRUE if q < 0.05 and the slope is positive", base_1998_2002 = "Month mean 1998-2002", mean_2020_2024 = "Month mean 2020-2024",
  f_trend = "STL strength of trend (0-1)", f_season = "STL strength of seasonality (0-1)", trend_1998 = "STL trend, 1998 mean", trend_2024 = "STL trend, 2024 mean", swing_1998 = "Seasonal swing 1998", swing_2024 = "Seasonal swing 2024",
  swing_sen = "Sen's slope of the seasonal swing", swing_p = "p-value of the swing trend", peak = "Usual peak month", trough = "Usual cleanest month", robust_z = "Robust z of the STL remainder", direction = "High or low anomaly",
  lisa_I_pm25_recent = "Local Moran's I, level", lisa_p_pm25_recent = "LISA pseudo p, level", lisa_q_pm25_recent = "LISA BH-adjusted p, level", lisa_pm25_recent = "LISA cluster, level",
  lisa_I_sen = "Local Moran's I, trend", lisa_p_sen = "LISA pseudo p, trend", lisa_q_sen = "LISA BH-adjusted p, trend", lisa_sen = "LISA cluster, trend",
  gi_z_pm25_recent = "Gi* z, level", gi_pm25_recent = "Gi* class, level", gi_q_pm25_recent = "Gi* BH-adjusted p, level", gi_z_sen = "Gi* z, trend", gi_sen = "Gi* class, trend", gi_q_sen = "Gi* BH-adjusted p, trend",
  test = "Test", variable = "Variable tested", weights = "Spatial weights", I = "Moran's I", expected = "Expected I under no clustering", p_perm = "Permutation p-value",
  lag_order = "Neighbour order", mean_distance_km = "Mean centroid distance at this order, km", se = "Standard error",
  target = "Model target", model = "Model", k = "Parameters", logLik = "Log-likelihood", AIC = "Akaike information criterion", rho = "Spatial lag parameter", lambda = "Spatial error parameter",
  pseudo_r2 = "Squared correlation of fitted and observed (includes the spatial term)", rmse = "Root mean squared error", resid_moran = "Moran's I of residuals", resid_moran_p = "p-value of residual Moran's I", dAIC = "AIC minus the best AIC", best = "TRUE for the lowest AIC",
  effect = "Direct, Spillover or Total", est = "Effect per 1 SD", lo = "95% CI low", hi = "95% CI high", sd_raw = "1 SD of the predictor in raw units",
  cv = "Random or spatial cross-validation", r2 = "Out-of-fold R2", mae = "Mean absolute error", shap = "SHAP value, ug/m3", mean_abs_shap = "Mean absolute SHAP, ug/m3", rank = "Importance rank",
  finding_id = "Finding number", section = "Dashboard section", source = "Notebook or script", measure = "Finding", value = "Value as text",
  Term = "Glossary term", Meaning = "Plain meaning", `How calculated` = "Method", `How to read` = "Reading guide", Source = "Source",
  from_table = "Many side table", from_column = "Many side column", to_table = "One side table", to_column = "One side column", cardinality = "Relationship type", cross_filter = "Filter direction",
  table = "Table name", column = "Column name", type = "R type when exported", description = "What the column holds")
dict <- bind_rows(lapply(names(tabs), \(t) tibble(table = t, column = names(tabs[[t]]), type = sapply(tabs[[t]], \(x) class(x)[1]))))
dict <- bind_rows(dict, tibble(table = "table_dictionary", column = c("table", "column", "type", "description"), type = "character")) |> mutate(description = desc[column])
if (anyNA(dict$description)) stop("Columns without a description: ", paste(unique(dict$column[is.na(dict$description)]), collapse = ", "))
put(dict |> mutate(description = coalesce(description, "")), "table_dictionary")

# ---- 6. Checks against the earlier outputs ----
nat <- function(t) filter(tabs[[t]], pcode == "BD")
rec <- uy |> filter(year >= 2020) |> group_by(pcode) |> summarise(nb = mean(pm25_popw)) |> inner_join(tabs$stat_exposure_recent, by = "pcode")
chk <- tibble(check = c("Share above 35 in 2024 vs notebook 02", "National 2020-2024 mean vs script 06", "Years lost, national, vs script 09", "Years lost, Savar, vs script 09",
                        "Band people vs unit people, largest relative gap", "Upazilas in the trend, exposure and cluster tables", "2020-2024 mean, grid method vs notebook 02, largest upazila gap"),
              here = c(sum(filter(nat("fact_band_year"), year == 2024, band %in% dim_band$band[6:9])$share), nat("stat_exposure_recent")$pm25_2020_2024, nat("stat_exposure_recent")$years_lost_who,
                       filter(tabs$stat_exposure_recent, pcode == filter(dim_geo, name == "Savar")$pcode)$years_lost_who,
                       max(abs(tabs$fact_band_year |> filter(year == 2024) |> group_by(pcode) |> summarise(s = sum(people)) |> inner_join(pop, by = "pcode") |> mutate(g = (s - people) / people) |> pull(g))),
                       min(sapply(c("stat_trend", "stat_exposure_recent", "stat_spatial_cluster"), \(t) sum(tabs[[t]]$pcode %in% dim_geo$pcode[dim_geo$level == "upazila"]))),
                       max(abs(rec$pm25_2020_2024 - rec$nb))),
              earlier = c(sum(filter(read_parquet("data/processed/pm25_who_bands.parquet"), level == "national", year == 2024, band == ">35 above IT-1")$share),
                          nat("stat_trend")$last5, filter(rt("09_life_expectancy.csv"), division == "Bangladesh")$years_lost_who, filter(rt("09_upazila_exposure.csv"), name == "Savar")$years_lost_who, 0, 507, 0),
              tolerance = c(0.001, 0.01, 0.01, 0.01, 1e-4, 0, 0.5)) |>
  mutate(ok = abs(here - earlier) <= tolerance)
print(chk, width = Inf)
if (!all(chk$ok)) warning("Some checks differ from the earlier outputs; see the table above")
cat(sprintf("\n%d tables written to %s\n", length(tabs), normalizePath(out)))
