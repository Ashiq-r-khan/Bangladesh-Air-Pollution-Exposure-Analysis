# 08 Spatial regression: OLS, SLX, spatial lag (SAR), spatial error (SEM), spatial Durbin (SDM, SDEM)
# Bangladesh Air Pollution Project | Md. Ashiqur Rahman Khan

library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(spdep)
library(spatialreg)
library(car)
library(lmtest)
library(patchwork)
setwd("E:/Bangladesh Air Pollution Project")
fig_dir <- "outputs/R/figures"
tab_dir <- "outputs/R/tables"
invisible(lapply(c(fig_dir, tab_dir), dir.create, recursive = TRUE, showWarnings = FALSE))
save_fig <- function(p, name, w = 11, h = 7) {print(p); ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 300, bg = "white")}
save_tab <- function(d, name, n = 20) {write.csv(d, file.path(tab_dir, name), row.names = FALSE); print(as_tibble(d), n = n, width = Inf)}
theme_set(theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), plot.title.position = "plot", panel.grid.minor = element_blank(), legend.position = "bottom"))
ug <- "\u00b5g/m\u00b3"
set.seed(2026)

# ---- 1. Data and predictors (skewed ones log-transformed, all standardised so effects are per 1 SD) ----
gpkg <- "data/processed/bgd_admin.gpkg"
div <- st_read(gpkg, layer = "division", quiet = TRUE) |> st_transform(32646)
d <- st_read(gpkg, layer = "upazila", quiet = TRUE) |> st_transform(32646) |> select(pcode) |>
  inner_join(read_parquet("data/processed/upazila_covariates.parquet"), by = "pcode") |>
  left_join(read.csv("outputs/R/tables/06_trend_all_units.csv") |> filter(level == "upazila") |> select(pcode, sen), by = "pcode") |>
  mutate(log_elevation = log1p(pmax(elevation_m, 0)), log_roads = log1p(road_km_per_km2))
nb <- readRDS("data/processed/upazila_neighbours.rds")$queen
stopifnot(nrow(d) == 507, length(nb) == 507, !anyNA(d$sen))
lw <- nb2listw(nb, style = "W")
tr <- trW(as(lw, "CsparseMatrix"), type = "mult")
preds <- c(log_pop_density = "Population density (log)", urban_share = "Urban share", dist_dhaka_km = "Distance to Dhaka", log_elevation = "Elevation (log)",
           log_roads = "Major road density (log)", dry_wind_ms = "Dry-season wind speed", rain_mm_day = "Rainfall")
targets <- c(pm25_recent = paste0("PM2.5 level, 2020-2024 (", ug, ")"), sen = paste0("PM2.5 trend, 1998-2024 (", ug, " per year)"))
dm <- st_drop_geometry(d) |> mutate(across(all_of(names(preds)), \(x) as.numeric(scale(x))))
cat("Variance inflation factors (OLS, level model):\n"); print(round(vif(lm(reformulate(names(preds), "pm25_recent"), dm)), 2))

# ---- 2. Fit the six models for each target ----
rs_tests <- get0("lm.RStests", ifnotfound = lm.LMtests)
fit_models <- function(y) {
  f <- reformulate(names(preds), y)
  list(OLS = lm(f, dm), SLX = lmSLX(f, dm, lw), SAR = lagsarlm(f, dm, lw, trs = tr), SEM = errorsarlm(f, dm, lw, trs = tr),
       SDM = lagsarlm(f, dm, lw, type = "mixed", trs = tr), SDEM = errorsarlm(f, dm, lw, Durbin = TRUE, trs = tr))
}
fits <- lapply(names(targets), fit_models) |> setNames(names(targets))
rs <- bind_rows(lapply(names(targets), \(y) {t <- rs_tests(fits[[y]]$OLS, lw, test = "all")
  tibble(target = y, test = names(t), statistic = sapply(t, \(x) unname(x$statistic)), p = sapply(t, \(x) x$p.value))}))
cat("\nRao score (Lagrange multiplier) tests on the OLS residuals:\n")
save_tab(rs, "08_lm_tests.csv")

# ---- 3. Model comparison ----
res_of <- function(m) if (inherits(m, "Sarlm")) residuals(m) else resid(m)
comparison <- bind_rows(lapply(names(targets), \(y) bind_rows(lapply(names(fits[[y]]), \(m) {f <- fits[[y]][[m]]; mi <- moran.test(res_of(f), lw)
  tibble(target = y, model = m, k = attr(logLik(f), "df"), logLik = as.numeric(logLik(f)), AIC = AIC(f), rho = if (!is.null(f$rho)) unname(f$rho) else NA_real_,
         lambda = if (!is.null(f$lambda)) unname(f$lambda) else NA_real_, pseudo_r2 = cor(dm[[y]] - res_of(f), dm[[y]])^2,
         rmse = sqrt(mean(res_of(f)^2)), resid_moran = mi$estimate[[1]], resid_moran_p = mi$p.value)})))) |>
  group_by(target) |> mutate(dAIC = AIC - min(AIC), best = dAIC == 0) |> ungroup() |> mutate(model = factor(model, names(fits[[1]])))
save_tab(comparison, "08_model_comparison.csv", n = 12)
lr <- bind_rows(lapply(names(targets), \(y) {m <- fits[[y]]
  bind_rows(lapply(list(c("SDM", "SAR"), c("SDM", "SEM"), c("SDM", "SLX"), c("SDEM", "SEM"), c("SDEM", "SLX")), \(p) {t <- LR.Sarlm(m[[p[1]]], m[[p[2]]])
    tibble(target = y, full = p[1], restricted = p[2], LR = as.numeric(t$statistic), df = as.numeric(t$parameter), p = as.numeric(t$p.value))}))}))
cat("\nLikelihood-ratio tests (full model against the restricted one):\n")
save_tab(lr, "08_lr_tests.csv")
bp <- bind_rows(lapply(names(targets), \(y) tibble(target = y, test = c("Breusch-Pagan, OLS", "Spatial Breusch-Pagan, SEM"),
  BP = c(unname(bptest(fits[[y]]$OLS)$statistic), unname(bptest.Sarlm(fits[[y]]$SEM)$statistic)), p = c(bptest(fits[[y]]$OLS)$p.value, bptest.Sarlm(fits[[y]]$SEM)$p.value))))
cat("\nHeteroskedasticity checks:\n"); print(bp)
cat(sprintf("\nSpatial multiplier 1 / (1 - rho) in the SDM: %.0f (level), %.0f (trend). Global impacts this large are unstable, so SEM is the model interpreted below; SDEM adds local spillovers.\n",
            1 / (1 - fits$pm25_recent$SDM$rho), 1 / (1 - fits$sen$SDM$rho)))

# ---- 4. Effects per 1 SD of each predictor: OLS, SEM, SEM with 6-nearest-neighbour weights, SDEM direct and spillover ----
sds <- sapply(st_drop_geometry(d)[names(preds)], sd)
tidy_ci <- function(est, se) tibble(est = as.numeric(est), lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)))
beta_tbl <- function(s) bind_cols(tibble(var = names(preds), effect = "Direct"), tidy_ci(s[names(preds), 1], s[names(preds), 2]))
sdem_tbl <- function(m) {im <- summary(impacts(m))
  bind_rows(lapply(c(Direct = "direct", Spillover = "indirect", Total = "total"), \(k) bind_cols(tibble(var = names(preds)), tidy_ci(im$impacts[[k]], im$se[[k]]))), .id = "effect")}
lw_k <- nb2listw(readRDS("data/processed/upazila_neighbours.rds")$knn6, style = "W")
effects <- bind_rows(lapply(names(targets), \(y) bind_rows(
  beta_tbl(summary(fits[[y]]$OLS)$coefficients) |> mutate(model = "OLS"),
  beta_tbl(summary(fits[[y]]$SEM)$Coef) |> mutate(model = "SEM"),
  beta_tbl(summary(errorsarlm(reformulate(names(preds), y), dm, lw_k))$Coef) |> mutate(model = "SEM (knn6)"),
  sdem_tbl(fits[[y]]$SDEM) |> mutate(model = "SDEM")) |> mutate(target = y))) |>
  mutate(sd_raw = sds[var]) |> relocate(target, model, var, effect)
save_tab(effects, "08_effects.csv", n = 28)

# ---- 5. Figures ----
p_cmp <- ggplot(comparison, aes(model, dAIC, fill = best)) + geom_col(width = 0.65) + geom_text(aes(label = sprintf("I = %.2f", resid_moran)), vjust = -0.4, size = 3.2) +
  facet_wrap(~target, scales = "free_y", labeller = labeller(target = targets)) + scale_fill_manual(values = c("grey70", "#1967D2"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Modelling spatial dependence removes most of the clustering OLS leaves behind",
       subtitle = "AIC above the best model (lower is better); label = Moran's I of the model residuals (0 = no leftover clustering)", x = NULL, y = "AIC minus the best AIC")
save_fig(p_cmp, "08_model_comparison.png", h = 5.5)
show <- effects |> filter(model %in% c("OLS", "SEM", "SDEM"), effect != "Total") |>
  mutate(var = factor(preds[var], rev(preds)), model = factor(model, c("OLS", "SEM", "SDEM")),
         effect = factor(if_else(effect == "Direct", "Own upazila (direct)", "Neighbouring upazilas (spillover, SDEM only)"), c("Own upazila (direct)", "Neighbouring upazilas (spillover, SDEM only)")),
         sig = if_else(p < 0.05, "p < 0.05", "p >= 0.05"))
coef_plot <- function(y) ggplot(filter(show, target == y), aes(est, var, colour = model)) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi, shape = sig), position = position_dodge(width = 0.6), size = 0.35) +
  facet_wrap(~effect, nrow = 1) + scale_colour_manual(values = c(OLS = "grey55", SEM = "#1967D2", SDEM = "#6b2d5c"), name = NULL) +
  scale_shape_manual(values = c("p < 0.05" = 16, "p >= 0.05" = 1), name = NULL) + labs(title = targets[[y]], x = "Change per 1 SD increase in the predictor", y = NULL)
p_coef <- (coef_plot("pm25_recent") / coef_plot("sen")) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p_coef <- p_coef + plot_annotation(title = "Once spatial dependence is modelled, distance to Dhaka stays; wind and rainfall shrink",
                                   subtitle = "OLS against the spatial error model (SEM) and the spatial Durbin error model (SDEM); 95% CI; hollow = not significant.\nAssociations across upazilas, not causal effects",
                                   theme = theme(plot.title = element_text(face = "bold")))
save_fig(p_coef, "08_effects.png", w = 12, h = 10)
rmap <- bind_rows(lapply(names(targets), \(y) bind_rows(lapply(c("OLS", "SEM"), \(m) tibble(pcode = d$pcode, target = y, model = m, resid = res_of(fits[[y]][[m]])))))) |>
  left_join(comparison |> select(target, model, resid_moran) |> mutate(model = as.character(model)), by = c("target", "model")) |>
  group_by(target) |> mutate(z = resid / sd(resid[model == "OLS"]), panel = sprintf("%s\n%s residuals (Moran's I = %.2f)", sub(" \\(.*", "", targets[target]), model, resid_moran)) |> ungroup()
rmap <- d |> select(pcode) |> inner_join(rmap, by = "pcode", relationship = "one-to-many")
p_res <- ggplot() + geom_sf(data = rmap, aes(fill = pmax(pmin(z, 3), -3)), colour = NA) + geom_sf(data = div, fill = NA, colour = "grey25", linewidth = 0.3) +
  facet_wrap(~panel, nrow = 1) + scale_fill_distiller(palette = "RdBu", limits = c(-3, 3), name = "Residual (SD of OLS residuals)") +
  labs(title = "OLS leaves regional patches of error; the spatial error model removes most of them", subtitle = "Red: model under-predicts; blue: model over-predicts. Capped at \u00b13 SD") +
  theme_void() + theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"), plot.title = element_text(face = "bold"), strip.text = element_text(size = 9))
save_fig(p_res, "08_residual_maps.png", w = 14, h = 7)

# ---- 6. Key numbers ----
key <- bind_rows(lapply(names(targets), \(y) {cm <- filter(comparison, target == y); o <- filter(cm, model == "OLS"); s <- filter(cm, model == "SEM")
  ef <- filter(effects, target == y, effect == "Direct") |> select(model, var, est, p) |> pivot_wider(names_from = model, values_from = c(est, p))
  top <- ef |> arrange(desc(abs(est_SEM))) |> slice(1:3)
  tibble(target = y, measure = c("AIC ranking", "Rao score tests on OLS (robust lag / robust error)", "Spatial parameters (SDM rho / SEM lambda)", "Residual Moran's I, OLS vs SEM",
                                 "Largest SEM effects per 1 SD (OLS in brackets)", "Coefficients that shrink most from OLS to SEM"),
         value = c(paste(arrange(cm, AIC)$model, collapse = " < "),
                   with(filter(rs, target == y), sprintf("%.1f (p %.2g) / %.1f (p %.2g)", statistic[grepl("^(adjRSlag|RLMlag)$", test)], p[grepl("^(adjRSlag|RLMlag)$", test)], statistic[grepl("^(adjRSerr|RLMerr)$", test)], p[grepl("^(adjRSerr|RLMerr)$", test)])),
                   sprintf("%.2f / %.2f", filter(cm, model == "SDM")$rho, s$lambda), sprintf("%.2f vs %.2f", o$resid_moran, s$resid_moran),
                   paste(sprintf("%s %+.2f (%+.2f)", preds[top$var], top$est_SEM, top$est_OLS), collapse = "; "),
                   paste(sprintf("%s %+.2f to %+.2f", preds[(ef |> arrange(desc(abs(est_OLS - est_SEM))) |> slice(1:2))$var], (ef |> arrange(desc(abs(est_OLS - est_SEM))) |> slice(1:2))$est_OLS, (ef |> arrange(desc(abs(est_OLS - est_SEM))) |> slice(1:2))$est_SEM), collapse = "; ")))}))
cat("\nPredictor SDs in raw units (log variables in log units):\n"); print(round(sds, 3))
save_tab(key, "08_key_numbers.csv")
