# 09 WHO exposure bands, exposure distribution, reductions needed and AQLI-style life expectancy
# Bangladesh Air Pollution Project | Md. Ashiqur Rahman Khan

library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)
library(patchwork)
setwd("E:/Bangladesh Air Pollution Project")
fig_dir <- "outputs/R/figures"
tab_dir <- "outputs/R/tables"
invisible(lapply(c(fig_dir, tab_dir), dir.create, recursive = TRUE, showWarnings = FALSE))
save_fig <- function(p, name, w = 11, h = 7) {print(p); ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 300, bg = "white")}
save_tab <- function(d, name, n = 20) {write.csv(d, file.path(tab_dir, name), row.names = FALSE); print(as_tibble(d), n = n, width = Inf)}
theme_set(theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), plot.title.position = "plot", panel.grid.minor = element_blank(), legend.position = "bottom"))
ug <- "\u00b5g/m\u00b3"
who <- c("WHO guideline" = 5, "IT-4" = 10, "IT-3" = 15, "IT-2" = 25, "IT-1 / national standard" = 35)
aqli <- 0.098
bands <- c(0, 5, 10, 15, 25, 35, 50, 75, 100, Inf)
band_labs <- c("<=5 guideline", "5-10 (IT-4)", "10-15 (IT-3)", "15-25 (IT-2)", "25-35 (IT-1)", "35-50", "50-75", "75-100", ">100")
band_cols <- setNames(c("#1a9850", "#66bd63", "#a6d96a", "#d9ef8b", "#fee08b", "#fdae61", "#f46d43", "#d73027", "#7f0000"), band_labs)
wq <- function(x, w, p) {o <- order(x); cw <- cumsum(w[o]) / sum(w); sapply(p, \(q) x[o][which(cw >= q)[1]])}

# ---- 1. Grid cells with people, each assigned to an upazila (centre inside, else nearest) ----
gpkg <- "data/processed/bgd_admin.gpkg"
upz <- st_read(gpkg, layer = "upazila", quiet = TRUE)
div <- st_read(gpkg, layer = "division", quiet = TRUE)
grid <- sds("data/processed/pm25_annual_grid.nc")
pm <- grid[["pm25"]]
names(pm) <- 1998:2024
cells <- as.data.frame(c(grid[["pop_2022"]], pm), xy = TRUE) |> rename(pop = pop_2022) |> filter(pop > 0)
pts <- st_as_sf(cells[c("x", "y")], coords = c("x", "y"), crs = 4326) |> st_transform(32646)
upz_m <- st_transform(upz, 32646)
idx <- sapply(st_intersects(pts, upz_m), \(i) if (length(i)) i[1] else NA_integer_)
cat(sprintf("Cells with people: %d; centre outside every upazila: %d (%.3f%% of people), assigned to the nearest\n", nrow(cells), sum(is.na(idx)), 100 * sum(cells$pop[is.na(idx)]) / sum(cells$pop)))
idx[is.na(idx)] <- st_nearest_feature(pts[is.na(idx), ], upz_m)
cells <- cells |> mutate(pcode = upz$pcode[idx], division = div$name[match(upz$division_pcode[idx], div$pcode)], recent = rowMeans(across(all_of(as.character(2020:2024)))))
long <- cells |> select(pop, pcode, division, all_of(as.character(1998:2024))) |> pivot_longer(all_of(as.character(1998:2024)), names_to = "year", values_to = "pm25") |>
  mutate(year = as.integer(year), band = cut(pm25, bands, band_labs, right = TRUE))

# ---- 2. QA against the notebook 02 outputs ----
uy <- read_parquet("data/processed/pm25_unit_year.parquet")
wb <- read_parquet("data/processed/pm25_who_bands.parquet")
nat_check <- long |> group_by(year) |> summarise(popw = weighted.mean(pm25, pop)) |> inner_join(filter(uy, level == "national") |> select(year, nb02 = pm25_popw), by = "year")
upz_check <- long |> filter(year == 2024) |> group_by(pcode) |> summarise(popw = weighted.mean(pm25, pop)) |> inner_join(filter(uy, level == "upazila", year == 2024) |> select(pcode, nb02 = pm25_popw), by = "pcode")
qa <- tibble(check = c("People on the grid", "Max gap vs notebook 02, national yearly mean", "Max gap vs notebook 02, upazila 2024 mean", "Share above 35 in 2024: here vs notebook 02"),
             value = c(format(round(sum(cells$pop)), big.mark = ","), sprintf("%.3f %s", max(abs(nat_check$popw - nat_check$nb02)), ug), sprintf("%.2f %s", max(abs(upz_check$popw - upz_check$nb02)), ug),
                       sprintf("%.2f%% vs %.2f%%", 100 * with(filter(long, year == 2024), sum(pop[pm25 > 35]) / sum(pop)), 100 * sum(filter(wb, level == "national", year == 2024, band == ">35 above IT-1")$share))))
save_tab(qa, "09_qa.csv")

# ---- 3. National bands by year and the population-weighted exposure distribution ----
nat_bands <- long |> group_by(year, band) |> summarise(people = sum(pop), .groups = "drop") |> complete(year, band, fill = list(people = 0)) |>
  group_by(year) |> mutate(share = people / sum(people)) |> ungroup()
save_tab(nat_bands |> filter(year %in% c(1998, 2010, 2024)) |> mutate(people_m = people / 1e6), "09_national_bands_by_year.csv", n = 27)
write.csv(nat_bands, file.path(tab_dir, "09_national_bands_by_year.csv"), row.names = FALSE)
pct <- bind_rows(long |> mutate(division = "Bangladesh"), long) |> group_by(division, year) |>
  summarise(mean = weighted.mean(pm25, pop), tibble::as_tibble_row(setNames(wq(pm25, pop, c(0.1, 0.25, 0.5, 0.75, 0.9)), c("p10", "p25", "p50", "p75", "p90"))), min = min(pm25), max = max(pm25), .groups = "drop")
save_tab(pct |> filter(year %in% c(1998, 2024), division == "Bangladesh"), "09_exposure_percentiles.csv")
write.csv(pct, file.path(tab_dir, "09_exposure_percentiles.csv"), row.names = FALSE)
p_area <- ggplot(nat_bands, aes(year, share, fill = band)) + geom_area(colour = "white", linewidth = 0.2) +
  geom_line(data = filter(pct, division == "Bangladesh"), aes(year, mean / 130), inherit.aes = FALSE, colour = "black", linewidth = 0.9) +
  annotate("text", x = 2018, y = 1.02, label = "2018 spike", size = 3, vjust = 0) +
  scale_fill_manual(values = band_cols, drop = FALSE, name = paste0("Annual PM2.5 (", ug, ")")) + guides(fill = guide_legend(nrow = 1, reverse = FALSE)) +
  scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(~ . * 130, name = paste0("Population-weighted mean (", ug, ", black line)")), expand = expansion(mult = c(0, 0.06))) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = with(nat_bands |> group_by(year) |> summarise(s = sum(share[band %in% band_labs[6:9]])), sprintf("People living above 35 \u00b5g/m\u00b3 went from %.0f%% in 1998 to %.1f%% in 2024", 100 * s[year == 1998], 100 * s[year == 2024])), subtitle = "Share of the population by annual PM2.5 band, 1998-2024. 2022 WorldPop population held fixed, so the change is all air quality",
       x = NULL, y = "Share of population")
save_fig(p_area, "09_band_shares_over_time.png", w = 12, h = 6.5)
exc <- long |> filter(year %in% c(1998, 2010, 2024)) |> group_by(year) |> arrange(pm25, .by_group = TRUE) |>
  mutate(above = 1 - cumsum(pop) / sum(pop)) |> slice(unique(round(seq(1, n(), length.out = 600)))) |> ungroup()
p_exc <- ggplot(exc, aes(pm25, above, colour = factor(year))) + geom_vline(xintercept = who, colour = "grey70", linetype = "dashed") +
  annotate("text", x = who, y = 0.03, label = names(who), angle = 90, hjust = 0, vjust = -0.4, size = 3, colour = "grey40") + geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("1998" = "#fdae61", "2010" = "#d73027", "2024" = "#7f0000"), name = NULL) +
  scale_y_continuous(labels = scales::percent) + scale_x_continuous(breaks = seq(0, 140, 20)) +
  labs(title = "The exposure curve jumped right after 1998; 2024 sits a little left of 2010", subtitle = "Share of people living where annual PM2.5 is above each value. The upper tail was heavier in 2010 than in 2024", x = paste0("Annual PM2.5 (", ug, ")"), y = "Share of people above")
save_fig(p_exc, "09_exceedance_curve.png", w = 10, h = 6)

# ---- 4. Division bands, 2024 ----
div_bands <- long |> filter(year == 2024) |> group_by(division, band) |> summarise(people = sum(pop), .groups = "drop") |>
  group_by(division) |> mutate(share = people / sum(people)) |> ungroup()
save_tab(div_bands |> filter(people > 0) |> mutate(people_m = people / 1e6), "09_division_bands_2024.csv", n = 40)
div_order <- pct |> filter(year == 2024, division != "Bangladesh") |> arrange(mean) |> pull(division)
p_div <- ggplot(div_bands |> mutate(division = factor(division, div_order)), aes(people / 1e6, division, fill = band)) + geom_col(width = 0.7) +
  geom_text(data = pct |> filter(year == 2024, division != "Bangladesh") |> left_join(div_bands |> group_by(division) |> summarise(tot = sum(people)), by = "division") |> mutate(division = factor(division, div_order)),
            aes(tot / 1e6, division, label = sprintf("mean %.0f", mean)), inherit.aes = FALSE, hjust = -0.15, size = 3.2) +
  scale_fill_manual(values = band_cols, drop = TRUE, name = paste0("Annual PM2.5 (", ug, ")")) + scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = sprintf("In 2024, %.0f%% of Dhaka division's people lived above 75 \u00b5g/m\u00b3", 100 * sum(filter(div_bands, division == "Dhaka", band %in% band_labs[8:9])$share)), subtitle = "People by annual PM2.5 band and division, 2024", x = "People (millions)", y = NULL)
save_fig(p_div, "09_division_bands_2024.png", w = 11, h = 6)

# ---- 5. Cuts needed and AQLI-style life expectancy, 2020-2024 mean ----
upz_exp <- cells |> group_by(pcode) |> summarise(pm25 = weighted.mean(recent, pop),
  cut_to_35 = 100 * weighted.mean(pmax(0, 1 - 35 / recent), pop), cut_to_25 = 100 * weighted.mean(pmax(0, 1 - 25 / recent), pop), cut_to_5 = 100 * weighted.mean(pmax(0, 1 - 5 / recent), pop),
  years_lost_who = weighted.mean(aqli * pmax(0, recent - 5), pop), years_gain_35 = weighted.mean(aqli * pmax(0, recent - 35), pop), pop = sum(pop)) |>
  left_join(st_drop_geometry(upz) |> select(pcode, name, division_pcode), by = "pcode") |> left_join(st_drop_geometry(div) |> select(division_pcode = pcode, division = name), by = "division_pcode") |>
  relocate(pcode, name, division) |> arrange(desc(years_lost_who))
save_tab(upz_exp, "09_upazila_exposure.csv", n = 10)
le <- bind_rows(cells |> mutate(division = "Bangladesh"), cells) |> group_by(division) |>
  summarise(people_m = sum(pop) / 1e6, pm25 = weighted.mean(recent, pop), years_lost_who = weighted.mean(aqli * pmax(0, recent - 5), pop), years_gain_35 = weighted.mean(aqli * pmax(0, recent - 35), pop),
            years_gain_25 = weighted.mean(aqli * pmax(0, recent - 25), pop), life_years_lost_m = sum(pop * aqli * pmax(0, recent - 5)) / 1e6) |> arrange(desc(years_lost_who))
save_tab(le, "09_life_expectancy.csv")
upz_map <- upz |> select(pcode) |> inner_join(upz_exp, by = "pcode")
map_panel <- function(v, title, lab, pal) ggplot() + geom_sf(data = upz_map, aes(fill = .data[[v]]), colour = NA) + geom_sf(data = div, fill = NA, colour = "white", linewidth = 0.4) +
  scale_fill_distiller(palette = pal, direction = 1, name = lab) + labs(title = title) + theme_void() +
  theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"), plot.title = element_text(face = "bold", hjust = 0.5))
p_map <- (map_panel("cut_to_35", "Cut needed to meet the national standard (35)", "% reduction", "OrRd") | map_panel("years_lost_who", "Life expectancy lost vs the WHO guideline", "Years (AQLI method)", "RdPu")) +
  plot_annotation(title = "How far each upazila is from clean air, 2020-2024", caption = "AQLI method: 0.98 years of life expectancy per 10 \u00b5g/m\u00b3 above 5 (Ebenstein et al. 2017). An estimate from a published relationship, not a measured outcome",
                  theme = theme(plot.title = element_text(face = "bold")))
save_fig(p_map, "09_reduction_and_life_expectancy_map.png", h = 8)
p_le <- le |> filter(division != "Bangladesh") |> select(division, `Lost vs WHO guideline (5)` = years_lost_who, `Gained by meeting the national standard (35)` = years_gain_35) |>
  pivot_longer(-division) |> mutate(division = factor(division, rev(le$division[le$division != "Bangladesh"])), name = factor(name, c("Lost vs WHO guideline (5)", "Gained by meeting the national standard (35)"))) |>
  ggplot(aes(value, division, fill = name)) + geom_col(position = position_dodge(width = 0.75), width = 0.7) + geom_text(aes(label = sprintf("%.1f", value)), position = position_dodge(width = 0.75), hjust = -0.2, size = 3.2) +
  geom_vline(xintercept = le$years_lost_who[le$division == "Bangladesh"], linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("#6b2d5c", "#1967D2"), name = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = sprintf("The average person loses about %.1f years of life expectancy to PM2.5", le$years_lost_who[le$division == "Bangladesh"]),
       subtitle = sprintf("AQLI-style estimate from the 2020-2024 mean (dashed line: national average).\nMeeting Bangladesh's own annual standard of 35 would win back %.1f years on average", le$years_gain_35[le$division == "Bangladesh"]),
       x = "Years of life expectancy", y = NULL)
save_fig(p_le, "09_life_expectancy_by_division.png", w = 11, h = 6)

# ---- 6. Key numbers ----
n24 <- filter(nat_bands, year == 2024)
p24 <- filter(pct, division == "Bangladesh", year == 2024)
p98 <- filter(pct, division == "Bangladesh", year == 1998)
bd <- filter(le, division == "Bangladesh")
yrs99 <- nat_bands |> group_by(year) |> summarise(s = sum(share[band %in% band_labs[6:9]])) |> filter(s >= 0.99) |> pull(year)
key <- tibble(measure = c("People above the WHO guideline (5), every year", "People above 35 (IT-1 and national standard), 2024", "People above 50 / 75, 2024",
                          "Years with at least 99% of people above 35", "Median person's exposure, 1998 vs 2024", "Cleanest-air 10% vs dirtiest 10% of people, 2024 (P10 / P90)",
                          "Average cut needed to meet the national standard (35), 2020-2024", "Average cut needed to meet the WHO guideline (5)",
                          "Life expectancy lost vs WHO guideline (AQLI method), national", "Worst and best division, years lost", "Worst upazila, years lost",
                          "Gain from meeting the national standard (35)", "Total life-years at stake"),
              value = c(sprintf("%.1f%%", 100 * min(nat_bands |> group_by(year) |> summarise(s = 1 - sum(share[band == band_labs[1]])) |> pull(s))),
                        sprintf("%.1f million (%.1f%%)", sum(n24$people[n24$band %in% band_labs[6:9]]) / 1e6, 100 * sum(n24$share[n24$band %in% band_labs[6:9]])),
                        sprintf("%.0f%% / %.0f%%", 100 * sum(n24$share[n24$band %in% band_labs[7:9]]), 100 * sum(n24$share[n24$band %in% band_labs[8:9]])),
                        sprintf("%d of 27 (first %d; %s)", length(yrs99), min(yrs99), paste(yrs99, collapse = ", ")), sprintf("%.0f vs %.0f %s", p98$p50, p24$p50, "ug/m3"), sprintf("%.0f / %.0f ug/m3", p24$p10, p24$p90),
                        sprintf("%.0f%% (population-weighted)", weighted.mean(upz_exp$cut_to_35, upz_exp$pop)), sprintf("%.0f%%", weighted.mean(upz_exp$cut_to_5, upz_exp$pop)),
                        sprintf("%.1f years per person", bd$years_lost_who),
                        sprintf("%s %.1f / %s %.1f", le$division[le$division != "Bangladesh"][1], le$years_lost_who[le$division != "Bangladesh"][1], tail(le$division, 1), tail(le$years_lost_who, 1)),
                        sprintf("%s (%s) %.1f", upz_exp$name[1], upz_exp$division[1], upz_exp$years_lost_who[1]),
                        sprintf("%.1f years per person", bd$years_gain_35), sprintf("%.0f million life-years across %.0f million people", bd$life_years_lost_m, bd$people_m)))
save_tab(key, "09_key_numbers.csv")
