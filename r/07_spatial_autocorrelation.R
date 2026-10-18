# 07 Spatial autocorrelation: global Moran's I, correlogram, LISA, Getis-Ord Gi*
# Bangladesh Air Pollution Project | Md. Ashiqur Rahman Khan

library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(spdep)
library(patchwork)
setwd("E:/Bangladesh Air Pollution Project")
fig_dir <- "outputs/R/figures"
tab_dir <- "outputs/R/tables"
invisible(lapply(c(fig_dir, tab_dir), dir.create, recursive = TRUE, showWarnings = FALSE))
save_fig <- function(p, name, w = 11, h = 7) {print(p); ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 300, bg = "white")}
save_tab <- function(d, name, n = 20) {write.csv(d, file.path(tab_dir, name), row.names = FALSE); print(as_tibble(d), n = n, width = Inf)}
theme_set(theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), plot.title.position = "plot", panel.grid.minor = element_blank(), legend.position = "bottom"))
map_theme <- theme_void() + theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))
ug <- "\u00b5g/m\u00b3"
set.seed(2026)

# ---- 1. Load data (projected to UTM 46N so distances are in metres) ----
gpkg <- "data/processed/bgd_admin.gpkg"
div <- st_read(gpkg, layer = "division", quiet = TRUE) |> st_transform(32646)
d <- st_read(gpkg, layer = "upazila", quiet = TRUE) |> st_transform(32646) |> select(pcode) |>
  inner_join(read_parquet("data/processed/upazila_covariates.parquet"), by = "pcode") |>
  left_join(read.csv("outputs/R/tables/06_trend_all_units.csv") |> filter(level == "upazila") |> select(pcode, sen, pct_decade), by = "pcode") |>
  left_join(read_parquet("data/processed/pm25_unit_year.parquet") |> filter(level == "upazila", year == 2024) |> select(pcode, pop), by = "pcode") |>
  left_join(st_drop_geometry(div) |> select(division_pcode = pcode, division = name), by = "division_pcode")
stopifnot(nrow(d) == 507, !anyNA(d$sen), !anyNA(d$pop))
vars <- c(pm25_recent = "PM2.5 level, 2020-2024", sen = "PM2.5 trend, 1998-2024 (Sen's slope)")

# ---- 2. Spatial weights: queen contiguity (main), 6 nearest neighbours (robustness) ----
cent <- st_coordinates(st_point_on_surface(st_geometry(d)))
nb_q <- poly2nb(d, queen = TRUE)
cat("Upazilas with no touching neighbour:", d$name[card(nb_q) == 0], "-> linked to the nearest upazila\n")
k1 <- knn2nb(knearneigh(cent, k = 1))
for (i in which(card(nb_q) == 0)) {nb_q[[i]] <- k1[[i]]; nb_q[[k1[[i]]]] <- sort(c(nb_q[[k1[[i]]]], i))}
nb_k <- make.sym.nb(knn2nb(knearneigh(cent, k = 6)))
lws <- list(queen = nb2listw(nb_q, style = "W"), knn6 = nb2listw(nb_k, style = "W"))
saveRDS(list(queen = nb_q, knn6 = nb_k), "data/processed/upazila_neighbours.rds")
print(tibble(weights = names(lws), links = sapply(lws, \(w) sum(card(w$neighbours))), mean_neighbours = sapply(lws, \(w) mean(card(w$neighbours))),
             min = sapply(lws, \(w) min(card(w$neighbours))), max = sapply(lws, \(w) max(card(w$neighbours))),
             mean_link_km = sapply(lws, \(w) mean(unlist(nbdists(w$neighbours, cent))) / 1000)))

# ---- 3. Global Moran's I (analytical and 999 permutations) ----
moran_row <- function(x, lw) {t <- moran.test(x, lw); m <- moran.mc(x, lw, nsim = 999); tibble(I = t$estimate[[1]], expected = t$estimate[[2]], z = t$statistic[[1]], p = t$p.value, p_perm = m$p.value)}
global <- expand_grid(var = names(vars), weights = names(lws)) |> group_by(var, weights) |> summarise(moran_row(d[[var]], lws[[weights]]), .groups = "drop")
save_tab(global, "07_global_moran.csv")
res_moran <- read_parquet("outputs/models/05_oof_predictions.parquet") |>
  pivot_longer(contains("|"), names_to = c("target", "cv", "model"), names_sep = "\\|", values_to = "pred") |>
  mutate(resid = if_else(target == "pm25_recent", pm25_recent, pm25_change) - pred) |>
  group_by(target, cv, model) |> summarise(moran_row(resid[match(d$pcode, pcode)], lws$queen), .groups = "drop") |> arrange(target, cv, desc(I))
cat("\nMoran's I of the notebook 05 out-of-fold residuals (queen weights):\n")
save_tab(res_moran, "07_ml_residual_moran.csv", n = 12)

# ---- 4. Moran scatterplots ----
ms <- bind_rows(lapply(names(vars), \(v) tibble(var = vars[[v]], z = as.numeric(scale(d[[v]])), lag = lag.listw(lws$queen, as.numeric(scale(d[[v]])))))) |>
  mutate(quad = case_when(z >= 0 & lag >= 0 ~ "High-High", z < 0 & lag < 0 ~ "Low-Low", z >= 0 ~ "High-Low", TRUE ~ "Low-High"))
quad_cols <- c("High-High" = "#b2182b", "Low-Low" = "#2166ac", "High-Low" = "#ef8a62", "Low-High" = "#67a9cf", "Not significant" = "grey88")
ms_lab <- global |> filter(weights == "queen") |> mutate(var = vars[var], lab = sprintf("Moran's I = %.2f", I))
p_ms <- ggplot(ms, aes(z, lag)) + geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(aes(colour = quad), alpha = 0.7, size = 1.5) + geom_smooth(method = "lm", formula = y ~ x, colour = "black", se = FALSE, linewidth = 0.7) +
  geom_text(data = ms_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.5, fontface = "bold") +
  facet_wrap(~var) + scale_colour_manual(values = quad_cols, name = NULL) +
  labs(title = "Upazilas look like their neighbours, for both the level and the trend", subtitle = "Moran scatterplot: standardised value against the average of its queen-contiguity neighbours (slope = Moran's I).\nThe very high I partly reflects how smooth the satellite PM2.5 surface is, so read it as strong regional structure, not 507 independent readings",
       x = "Upazila value (z-score)", y = "Neighbour average (spatial lag of z)")
save_fig(p_ms, "07_moran_scatter.png", h = 6)

# ---- 5. Spatial correlogram: how far does the clustering reach? ----
lag_km <- sapply(nblag(nb_q, 10), \(n) mean(unlist(nbdists(n, cent)))) / 1000
corr <- bind_rows(lapply(names(vars), \(v) {r <- sp.correlogram(nb_q, d[[v]], order = 10, method = "I", style = "W")$res
  tibble(var = vars[[v]], order = 1:10, km = lag_km, I = r[, 1], expected = r[, 2], se = sqrt(r[, 3]), p = 2 * pnorm(-abs((r[, 1] - r[, 2]) / sqrt(r[, 3]))))}))
save_tab(corr, "07_correlogram.csv")
p_corr <- ggplot(corr, aes(km, I, colour = var)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = I - 1.96 * se, ymax = I + 1.96 * se, fill = var), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = c("#1967D2", "#6b2d5c"), name = NULL) + scale_fill_manual(values = c("#1967D2", "#6b2d5c"), name = NULL) +
  scale_x_continuous(breaks = round(lag_km), labels = \(x) paste0(x, " km\n(lag ", match(x, round(lag_km)), ")")) +
  labs(title = "Clustering is strong at short range and fades with distance", subtitle = "Moran's I by neighbour order (1 = touching, 2 = neighbours of neighbours...), plotted at the mean centroid distance of each order; band = 95% CI",
       x = NULL, y = "Moran's I")
save_fig(p_corr, "07_correlogram.png", w = 11, h = 6)

# ---- 6. Local Moran's I (LISA, 999 permutations) and Getis-Ord Gi* ----
w_gi <- nb2listw(include.self(nb_q), style = "B")
gi_breaks <- c(-Inf, -2.58, -1.96, -1.65, 1.65, 1.96, 2.58, Inf)
gi_labs <- c("Cold spot, 99%", "Cold spot, 95%", "Cold spot, 90%", "Not significant", "Hot spot, 90%", "Hot spot, 95%", "Hot spot, 99%")
for (v in names(vars)) {
  lm <- localmoran_perm(d[[v]], lws$queen, nsim = 999)
  p <- lm[, "Pr(folded) Sim"]
  gz <- as.numeric(localG(d[[v]], w_gi))
  d[[paste0("lisa_I_", v)]] <- lm[, "Ii"]
  d[[paste0("lisa_p_", v)]] <- p
  d[[paste0("lisa_q_", v)]] <- p.adjust(p, "BH")
  d[[paste0("lisa_", v)]] <- if_else(p < 0.05, as.character(attr(lm, "quadr")$mean), "Not significant")
  d[[paste0("gi_z_", v)]] <- gz
  d[[paste0("gi_", v)]] <- as.character(cut(gz, gi_breaks, gi_labs))
  d[[paste0("gi_q_", v)]] <- p.adjust(2 * pnorm(-abs(gz)), "BH")
}
clusters <- st_drop_geometry(d) |> select(pcode, name, district_pcode, division, pop, pm25_recent, sen, pct_decade, starts_with("lisa"), starts_with("gi"))
save_tab(clusters, "07_upazila_clusters.csv", n = 5)
summ <- bind_rows(lapply(names(vars), \(v) clusters |> mutate(cluster = .data[[paste0("lisa_", v)]]) |> group_by(cluster) |>
  summarise(upazilas = n(), fdr_survive = sum(.data[[paste0("lisa_q_", v)]] < 0.05), people_m = sum(pop) / 1e6, pm25_recent = weighted.mean(pm25_recent, pop), sen = weighted.mean(sen, pop),
            main_division = names(which.max(table(division)))) |> mutate(var = vars[[v]], share_people = people_m / sum(people_m)))) |>
  mutate(cluster = factor(cluster, names(quad_cols))) |> relocate(var, cluster, upazilas, fdr_survive, people_m, share_people) |> arrange(var, cluster)
cat("\nLISA clusters (p < 0.05; fdr_survive = still significant after BH adjustment):\n")
save_tab(summ, "07_cluster_summary.csv")
cat("\nGi* classes:\n"); print(lapply(names(vars), \(v) table(factor(d[[paste0("gi_", v)]], gi_labs))) |> setNames(names(vars)))

# ---- 7. Cluster maps ----
lisa_map <- function(v) ggplot() + geom_sf(data = d, aes(fill = factor(.data[[paste0("lisa_", v)]], names(quad_cols))), colour = "white", linewidth = 0.05) +
  geom_sf(data = div, fill = NA, colour = "grey25", linewidth = 0.35) + scale_fill_manual(values = quad_cols, drop = FALSE, name = NULL) +
  labs(title = vars[[v]]) + map_theme
p_lisa <- (lisa_map("pm25_recent") | lisa_map("sen")) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p_lisa <- p_lisa + plot_annotation(title = "LISA clusters: where high and low values group together",
                                   subtitle = "Local Moran's I, queen contiguity, 999 permutations, p < 0.05.\nHigh-High: high value surrounded by high values; High-Low and Low-High: spatial outliers",
                                   theme = theme(plot.title = element_text(face = "bold")))
save_fig(p_lisa, "07_lisa_maps.png", h = 8)
gi_cols <- setNames(c("#2166ac", "#67a9cf", "#d1e5f0", "grey92", "#fddbc7", "#ef8a62", "#b2182b"), gi_labs)
gi_map <- function(v) ggplot() + geom_sf(data = d, aes(fill = factor(.data[[paste0("gi_", v)]], gi_labs)), colour = "white", linewidth = 0.05) +
  geom_sf(data = div, fill = NA, colour = "grey25", linewidth = 0.35) + scale_fill_manual(values = gi_cols, drop = FALSE, name = NULL) +
  labs(title = vars[[v]]) + map_theme
p_gi <- (gi_map("pm25_recent") | gi_map("sen")) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p_gi <- p_gi + plot_annotation(title = "Getis-Ord Gi* hot and cold spots", subtitle = "Gi* z-score with the upazila and its queen-contiguity neighbours (binary weights); confidence from the normal approximation",
                               theme = theme(plot.title = element_text(face = "bold")))
save_fig(p_gi, "07_gistar_maps.png", h = 8)

# ---- 8. Key numbers ----
gq <- filter(global, weights == "queen")
gk <- filter(global, weights == "knn6")
hh <- filter(summ, var == vars[["pm25_recent"]], cluster == "High-High")
ll <- filter(summ, var == vars[["pm25_recent"]], cluster == "Low-Low")
hh_s <- filter(summ, var == vars[["sen"]], cluster == "High-High")
rx <- filter(res_moran, target == "pm25_recent", model == "XGBoost")
fade <- \(v) with(filter(corr, var == vars[[v]]), {i <- which(I < 0.1 | p >= 0.05)[1]; if (is.na(i)) "beyond lag 10" else sprintf("lag %d (about %.0f km)", order[i], km[i])})
key <- tibble(measure = c("Global Moran's I, PM2.5 level (queen / knn6)", "Global Moran's I, PM2.5 trend (queen / knn6)",
                          "Clustering of the level fades below I = 0.1 at", "Clustering of the trend fades below I = 0.1 at",
                          "High-High cluster, PM2.5 level", "Low-Low cluster, PM2.5 level", "High-High cluster, PM2.5 trend",
                          "Gi* hot spots at 99%, level / trend", "Moran's I of XGBoost residuals, random vs spatial CV (level)"),
              value = c(sprintf("%.2f / %.2f (permutation p %.3f)", gq$I[gq$var == "pm25_recent"], gk$I[gk$var == "pm25_recent"], gq$p_perm[gq$var == "pm25_recent"]),
                        sprintf("%.2f / %.2f (permutation p %.3f)", gq$I[gq$var == "sen"], gk$I[gk$var == "sen"], gq$p_perm[gq$var == "sen"]),
                        fade("pm25_recent"), fade("sen"),
                        sprintf("%d upazilas, %.1f million people (%.0f%%), mean %.1f ug/m3, mostly %s", hh$upazilas, hh$people_m, 100 * hh$share_people, hh$pm25_recent, hh$main_division),
                        sprintf("%d upazilas, %.1f million people (%.0f%%), mean %.1f ug/m3, mostly %s", ll$upazilas, ll$people_m, 100 * ll$share_people, ll$pm25_recent, ll$main_division),
                        sprintf("%d upazilas, %.1f million people, mean slope +%.2f ug/m3 per year, mostly %s", hh_s$upazilas, hh_s$people_m, hh_s$sen, hh_s$main_division),
                        sprintf("%d / %d upazilas", sum(d$gi_pm25_recent == "Hot spot, 99%"), sum(d$gi_sen == "Hot spot, 99%")),
                        sprintf("%.2f vs %.2f", rx$I[rx$cv == "random"], rx$I[rx$cv == "spatial"])))
save_tab(key, "07_key_numbers.csv")
