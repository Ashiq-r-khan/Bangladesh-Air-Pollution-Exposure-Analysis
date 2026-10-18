# 06 Trends and seasonality: Mann-Kendall, Sen's slope, seasonal Mann-Kendall, STL
# Bangladesh Air Pollution Project | Md. Ashiqur Rahman Khan

library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(trend)
library(modifiedmk)
library(patchwork)
setwd("E:/Bangladesh Air Pollution Project")
fig_dir <- "outputs/R/figures"
tab_dir <- "outputs/R/tables"
invisible(lapply(c(fig_dir, tab_dir), dir.create, recursive = TRUE, showWarnings = FALSE))
save_fig <- function(p, name, w = 11, h = 7) {print(p); ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 300, bg = "white")}
save_tab <- function(d, name, n = 20) {write.csv(d, file.path(tab_dir, name), row.names = FALSE); print(as_tibble(d), n = n, width = Inf)}
theme_set(theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), plot.title.position = "plot", panel.grid.minor = element_blank(), legend.position = "bottom"))
ug <- "\u00b5g/m\u00b3"

# ---- 1. Load data ----
yr <- read_parquet("data/processed/pm25_unit_year.parquet")
mo <- read_parquet("data/processed/pm25_unit_month.parquet") |> mutate(date = as.Date(date), year = as.integer(format(date, "%Y")), month = as.integer(format(date, "%m")))
upz <- st_read("data/processed/bgd_admin.gpkg", layer = "upazila", quiet = TRUE)
div <- st_read("data/processed/bgd_admin.gpkg", layer = "division", quiet = TRUE)
unit_order <- c("Bangladesh", yr |> filter(level == "division", year >= 2020) |> group_by(name) |> summarise(v = mean(pm25_popw)) |> arrange(desc(v)) |> pull(name))
div_cols <- setNames(c("black", RColorBrewer::brewer.pal(8, "Dark2")), c("Bangladesh", sort(unit_order[-1])))
cat("Units by level:\n"); print(table(distinct(yr, level, pcode)$level))

# ---- 2. Annual trends for every unit (Mann-Kendall with Hamed-Rao correction, Sen's slope) ----
trend_row <- function(x, year) {
  s <- sens.slope(x)
  h <- mmkh(x)
  b <- s$estimates[[1]]
  tibble(tau = h[["Tau"]], z = h[["Original Z"]], p = h[["old P.value"]], z_hr = h[["Corrected Zc"]], p_hr = h[["new P-value"]],
         sen = b, sen_lo = s$conf.int[1], sen_hi = s$conf.int[2], sen_no2018 = sens.slope(x[year != 2018])$estimates[[1]],
         fit0 = median(x - b * (year - min(year))), first5 = mean(x[year <= min(year) + 4]), last5 = mean(x[year >= max(year) - 4]),
         pct_decade = 1000 * b / mean(x[year <= min(year) + 4]))
}
trends <- yr |> arrange(year) |> group_by(level, pcode, name) |> summarise(trend_row(pm25_popw, year), .groups = "drop") |>
  group_by(level) |> mutate(q_hr = p.adjust(p_hr, "BH")) |> ungroup() |>
  mutate(level = factor(level, c("national", "division", "district", "upazila"))) |> arrange(level, desc(sen))
save_tab(trends, "06_trend_all_units.csv", n = 9)
cat("\nSignificant increase (Hamed-Rao p, BH-adjusted q < 0.05) by level:\n")
print(trends |> group_by(level) |> summarise(units = n(), sig_up = sum(q_hr < 0.05 & sen > 0), median_sen = median(sen), min_sen = min(sen), max_sen = max(sen)))

# ---- 3. Has the rise sped up? Sen's slope by period ----
pfit <- yr |> filter(level %in% c("national", "division")) |> arrange(year) |>
  mutate(period = if_else(year <= 2010, "1998_2010", "2011_2024")) |>
  group_by(level, pcode, name, period) |> summarise(s = list(sens.slope(pm25_popw)), y0 = min(year), a = median(pm25_popw - s[[1]]$estimates[[1]] * (year - y0)), .groups = "drop") |>
  mutate(sen = sapply(s, \(t) t$estimates[[1]]), p = sapply(s, \(t) t$p.value)) |> select(-s)
periods <- pfit |> select(level, name, period, sen, p) |> pivot_wider(names_from = period, values_from = c(sen, p)) |> arrange(match(name, unit_order))
save_tab(periods, "06_trend_periods.csv")

# ---- 4. Figure: national and division annual series with Sen's slope lines ----
fit_lab <- trends |> filter(level %in% c("national", "division")) |>
  mutate(panel = factor(sprintf("%s\n+%.2f %s per year", name, sen, ug), levels = sprintf("%s\n+%.2f %s per year", name, sen, ug)[match(unit_order, name)]))
ann <- yr |> filter(level %in% c("national", "division")) |> mutate(period = if_else(year <= 2010, "1998_2010", "2011_2024")) |>
  inner_join(select(fit_lab, pcode, panel, sen, fit0), by = "pcode") |> inner_join(select(pfit, pcode, period, y0, a, sen_p = sen), by = c("pcode", "period")) |>
  mutate(fit = fit0 + sen * (year - 1998), fit_p = a + sen_p * (year - y0))
p_trend <- ggplot(ann, aes(year)) +
  geom_line(aes(y = pm25_popw), colour = "grey55") + geom_point(aes(y = pm25_popw), colour = "grey35", size = 0.9) +
  geom_line(aes(y = fit), colour = "#1967D2", linewidth = 1) + geom_line(aes(y = fit_p, group = period), colour = "#6b2d5c", linetype = "dashed", linewidth = 0.8) +
  geom_point(data = filter(ann, year == 2018), aes(y = pm25_popw), colour = "#d08c60", size = 2.4) +
  facet_wrap(~panel, ncol = 3) +
  labs(title = "PM2.5 climbed until about 2010, then held at a high plateau",
       subtitle = "Grey: population-weighted annual PM2.5. Blue: Sen's slope 1998-2024. Dashed: Sen's slope for 1998-2010 and 2011-2024.\nOrange: 2018, a spike in the V6.GL.03 release; Sen's slope barely moves without it",
       x = NULL, y = paste0("PM2.5 (", ug, ")"), caption = "Source: ACAG SatPM2.5 V6.GL.03, WorldPop 2022")
save_fig(p_trend, "06_annual_trends_sen.png", h = 9)

# ---- 5. Upazila trend map and baseline vs slope ----
upz_tr <- upz |> select(pcode, division_pcode) |> left_join(filter(trends, level == "upazila"), by = "pcode")
sig_n <- sum(upz_tr$q_hr < 0.05 & upz_tr$sen > 0, na.rm = TRUE)
map_panel <- function(var, title, lab) ggplot() + geom_sf(data = upz_tr, aes(fill = .data[[var]]), colour = NA) +
  geom_sf(data = filter(upz_tr, q_hr >= 0.05), fill = "grey80", colour = NA) + geom_sf(data = div, fill = NA, colour = "white", linewidth = 0.4) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = lab) + labs(title = title) +
  theme_void() + theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"), plot.title = element_text(face = "bold", hjust = 0.5))
p_map <- (map_panel("sen", "Absolute rise", paste0(ug, " per year")) | map_panel("pct_decade", "Relative rise", "% of 1998-2002 level per decade")) +
  plot_annotation(title = "Where PM2.5 has risen fastest, by upazila (1998-2024)",
                  subtitle = sprintf("Sen's slope; %d of %d upazilas show a significant rise (Hamed-Rao Mann-Kendall, BH-adjusted q < 0.05). Grey: not significant", sig_n, nrow(upz_tr)),
                  theme = theme(plot.title = element_text(face = "bold")))
save_fig(p_map, "06_upazila_sen_map.png", h = 8)
upz_sc <- trends |> filter(level == "upazila") |> left_join(st_drop_geometry(select(upz, pcode, division_pcode)), by = "pcode") |>
  left_join(st_drop_geometry(select(div, division_pcode = pcode, division = name)), by = "division_pcode")
rho <- cor.test(upz_sc$first5, upz_sc$sen, method = "spearman", exact = FALSE)
rho_pct <- cor.test(upz_sc$first5, upz_sc$pct_decade, method = "spearman", exact = FALSE)
p_sc <- ggplot(upz_sc, aes(first5, sen, colour = division)) + geom_point(alpha = 0.7, size = 1.6) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, colour = "black", linetype = "dashed", se = FALSE, linewidth = 0.6) +
  scale_colour_manual(values = div_cols, name = NULL) +
  labs(title = "Dirtier upazilas rose slightly faster in absolute terms, cleaner ones faster in relative terms",
       subtitle = sprintf("Spearman rho: baseline vs absolute slope %.2f (p = %.2g); baseline vs relative slope %.2f (p = %.2g)", rho$estimate, rho$p.value, rho_pct$estimate, rho_pct$p.value),
       x = paste0("1998-2002 mean PM2.5 (", ug, ")"), y = paste0("Sen's slope (", ug, " per year)"))
save_fig(p_sc, "06_baseline_vs_slope.png", w = 10, h = 6.5)

# ---- 6. Seasonal Mann-Kendall and Sen's slope by calendar month ----
mo_main <- mo |> filter(level %in% c("national", "division")) |> arrange(date)
smk <- mo_main |> group_by(level, pcode, name) |> summarise(x = list(ts(pm25_popw, start = c(1998, 1), frequency = 12)), .groups = "drop") |>
  mutate(smk_z = sapply(x, \(t) unname(smk.test(t)$statistic)), smk_p = sapply(x, \(t) smk.test(t)$p.value), sea_sen = sapply(x, sea.sens.slope)) |>
  select(-x) |> arrange(match(name, unit_order))
save_tab(smk, "06_seasonal_mann_kendall.csv")
month_slopes <- mo_main |> group_by(level, pcode, name, month) |>
  summarise(s = list(sens.slope(pm25_popw)), base = mean(pm25_popw[year <= 2002]), last5 = mean(pm25_popw[year >= 2020]), .groups = "drop") |>
  mutate(sen = sapply(s, \(t) t$estimates[[1]]), p = sapply(s, \(t) t$p.value), pct_decade = 1000 * sen / base) |> select(-s) |> arrange(match(name, unit_order), month)
save_tab(month_slopes, "06_month_slopes.csv", n = 12)
hm <- month_slopes |> mutate(name = factor(name, rev(unit_order)), month = factor(month.abb[month], month.abb))
heat <- function(var, title, fmt) ggplot(hm, aes(month, name, fill = .data[[var]])) + geom_tile(colour = "white") +
  geom_text(aes(label = if_else(p < 0.05, sprintf(fmt, .data[[var]]), paste0("(", sprintf(fmt, .data[[var]]), ")"))), size = 2.8) +
  scale_fill_distiller(palette = "RdPu", direction = 1, name = NULL) + labs(title = title, x = NULL, y = NULL) + theme(panel.grid = element_blank())
p_heat <- (heat("sen", paste0("Absolute: ", ug, " per year"), "%.1f") / heat("pct_decade", "Relative: % of 1998-2002 level per decade", "%.0f")) +
  plot_annotation(title = "Winter rose most in absolute terms, the monsoon most in relative terms",
                  subtitle = "Sen's slope of each calendar month across 1998-2024. Values in brackets: Mann-Kendall p >= 0.05",
                  theme = theme(plot.title = element_text(face = "bold")))
save_fig(p_heat, "06_month_slope_heatmap.png", w = 11, h = 9)

# ---- 7. STL decomposition (log scale, seasonal window 13 so the seasonal shape can drift) ----
stl_fit <- mo_main |> group_by(level, pcode, name) |> group_modify(\(d, k) {
  f <- stl(ts(log(d$pm25_popw), start = c(1998, 1), frequency = 12), s.window = 13, robust = TRUE)$time.series
  tibble(date = d$date, year = d$year, month = d$month, observed = d$pm25_popw, trend = exp(f[, "trend"]), seasonal = f[, "seasonal"], remainder = f[, "remainder"])
}) |> ungroup()
swing <- stl_fit |> group_by(level, name, year) |> summarise(swing = exp(max(seasonal) - min(seasonal)), peak = month.abb[month[which.max(seasonal)]], trough = month.abb[month[which.min(seasonal)]], .groups = "drop")
stl_summary <- stl_fit |> group_by(level, name) |>
  summarise(f_trend = max(0, 1 - var(remainder) / var(log(trend) + remainder)), f_season = max(0, 1 - var(remainder) / var(seasonal + remainder)),
            trend_1998 = mean(trend[year == 1998]), trend_2024 = mean(trend[year == 2024]), .groups = "drop") |>
  left_join(swing |> arrange(year) |> group_by(level, name) |> summarise(swing_1998 = first(swing), swing_2024 = last(swing), swing_sen = sens.slope(swing)$estimates[[1]], swing_p = sens.slope(swing)$p.value, peak = names(which.max(table(peak))), trough = names(which.max(table(trough))), .groups = "drop"), by = c("level", "name")) |>
  arrange(match(name, unit_order))
save_tab(stl_summary, "06_stl_summary.csv")
anomalies <- stl_fit |> group_by(level, name) |> mutate(rz = (remainder - median(remainder)) / mad(remainder)) |> ungroup() |>
  filter(abs(rz) > 3) |> mutate(direction = if_else(rz > 0, "high", "low")) |> select(level, name, date, observed, trend, remainder, rz, direction) |> arrange(date)
save_tab(anomalies, "06_stl_anomalies.csv")
cat("\nSTL anomaly months (|robust z| > 3) by year, national and divisions:\n"); print(table(format(anomalies$date, "%Y"), anomalies$direction))
nat_stl <- stl_fit |> filter(level == "national") |> select(date, observed, trend, seasonal, remainder) |> pivot_longer(-date) |>
  mutate(name = factor(name, c("observed", "trend", "seasonal", "remainder"), c(paste0(c("Observed (", "Trend ("), ug, ")"), "Seasonal (log)", "Remainder (log)")))
p_stl <- ggplot(nat_stl, aes(date, value)) + geom_line(colour = "#1967D2", linewidth = 0.5) +
  geom_point(data = filter(nat_stl, name == "Remainder (log)", date %in% filter(anomalies, level == "national")$date), colour = "#d08c60", size = 2) +
  facet_wrap(~name, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(title = "STL decomposition of national monthly PM2.5, 1998-2024", subtitle = "Robust STL on log PM2.5 (seasonal window 13); trend back-transformed to \u00b5g/m\u00b3. Orange: anomaly months (|robust z| > 3)", x = NULL, y = NULL) +
  theme(strip.placement = "outside")
save_fig(p_stl, "06_stl_national.png", h = 9)
p_div <- (ggplot(stl_fit |> mutate(name = factor(name, unit_order)), aes(date, trend, colour = name, linewidth = name == "Bangladesh")) + geom_line() +
            scale_linewidth_manual(values = c(0.6, 1.3), guide = "none") + scale_colour_manual(values = div_cols, name = NULL) +
            labs(title = "STL trend component", x = NULL, y = paste0("PM2.5 trend (", ug, ")"))) |
  (ggplot(swing |> mutate(name = factor(name, unit_order)), aes(year, swing, colour = name, linewidth = name == "Bangladesh")) + geom_line() +
     scale_linewidth_manual(values = c(0.6, 1.3), guide = "none") + scale_colour_manual(values = div_cols, name = NULL) +
     labs(title = "Seasonal swing (peak month / trough month)", x = NULL, y = "Ratio"))
p_div <- p_div + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p_div <- p_div + plot_annotation(title = "The trend rose to a plateau while the winter-monsoon gap narrowed everywhere except Chattogram", theme = theme(plot.title = element_text(face = "bold")))
save_fig(p_div, "06_stl_divisions.png", w = 12, h = 6.5)

# ---- 8. Key numbers ----
nat <- filter(trends, level == "national")
dv <- filter(trends, level == "division")
up <- filter(trends, level == "upazila")
nm <- filter(month_slopes, level == "national")
ns <- filter(stl_summary, level == "national")
np <- filter(periods, level == "national")
key <- tibble(measure = c("National Sen's slope, 1998-2024", "National Sen's slope without 2018", "National relative rise",
                          "Mann-Kendall (Hamed-Rao corrected)", "National slope 1998-2010 vs 2011-2024 (p)", "Fastest-rising division", "Slowest-rising division",
                          "Upazilas with a significant rise (BH q < 0.05)", "Fastest-rising upazila", "Slowest-rising upazila",
                          "Baseline vs absolute slope (upazila Spearman rho)", "Month with the largest absolute rise", "Month with the largest relative rise",
                          "National seasonal swing (peak / trough), 1998 vs 2024", "STL strength of trend / seasonality"),
              value = c(sprintf("+%.2f ug/m3 per year (95%% CI %.2f to %.2f)", nat$sen, nat$sen_lo, nat$sen_hi), sprintf("+%.2f ug/m3 per year", nat$sen_no2018),
                        sprintf("+%.0f%% of the 1998-2002 level per decade", nat$pct_decade), sprintf("tau %.2f, Zc %.2f, p %.2g", nat$tau, nat$z_hr, nat$p_hr),
                        sprintf("%.2f (%.2g) vs %.2f (%.2g) ug/m3 per year", np$sen_1998_2010, np$p_1998_2010, np$sen_2011_2024, np$p_2011_2024),
                        with(slice_max(dv, sen), sprintf("%s (+%.2f)", name, sen)), with(slice_min(dv, sen), sprintf("%s (+%.2f)", name, sen)),
                        sprintf("%d of %d", sum(up$q_hr < 0.05 & up$sen > 0), nrow(up)),
                        with(slice_max(up, sen), sprintf("%s (+%.2f)", name, sen)), with(slice_min(up, sen), sprintf("%s (+%.2f)", name, sen)),
                        sprintf("%.2f (p %.2g)", rho$estimate, rho$p.value),
                        with(slice_max(nm, sen), sprintf("%s (+%.2f ug/m3 per year)", month.name[month], sen)),
                        with(slice_max(nm, pct_decade), sprintf("%s (+%.0f%% per decade)", month.name[month], pct_decade)),
                        sprintf("%.1fx vs %.1fx", ns$swing_1998, ns$swing_2024), sprintf("%.2f / %.2f", ns$f_trend, ns$f_season)))
save_tab(key, "06_key_numbers.csv")
