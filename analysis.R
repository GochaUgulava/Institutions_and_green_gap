
# =====================================================================
# analysis.R 
# =====================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(zoo); library(sandwich); library(lmtest)
  library(cluster); library(ggrepel)
})
set.seed(2026)

# ---------------------------------------------------------------------
# 0. PARAMETERS
# ---------------------------------------------------------------------
PROJECT_DIR <- "G:/My Drive/Working Papers/Institutional Quality and Multidimensional Sustainable Development in EU Candidate Countries"
PROCESSED_DIR <- "data/processed"
OUT_DIR   <- "output"
setwd(PROJECT_DIR)
INPUT <- file.path(PROCESSED_DIR, "master_panel_wide.csv")


EST_WIN <- 2012:2023           
REF_YEAR <- 2023              
INST_WIN <- 2009:2025          
P0 <- 2012:2014; P1 <- 2021:2023    
MAX_GAP <- 2                  
USE_EXTENDED <- FALSE          
N_MC <- 5000; N_BOOT <- 2000   

 
events <- tribble(
  ~iso3, ~year, ~label,
  "GEO", 2023, "candidate status",
  "GEO", 2024, "elections; accession paused",
  "MDA", 2022, "candidate status",
  "MDA", 2024, "talks opened; referendum"
)

inst_core <- c("wgi_va","wgi_pv","wgi_ge","wgi_rq","wgi_rl","wgi_cc","vdem_libdem","vdem_legcon")
cap_vars  <- c("wgi_ge","wgi_rq","wgi_rl","wgi_cc","vdem_rule","vdem_corr_inv")    
dem_vars  <- c("wgi_va","vdem_libdem","vdem_polyarchy","vdem_frefair","vdem_altinf",
               "vdem_frassoc","vdem_cspart","vdem_legcon","vdem_jucon")            

 
ind <- tribble(
  ~var,        ~dim,     ~sub,          ~tier,      ~dir, ~log,  ~scored,
  "gdp_pc",    "econ",   "econ",        "core",      1,   TRUE,  FALSE,
  "emp",       "econ",   "econ",        "core",      1,   FALSE, FALSE,
  "p0p50",     "social", "social",      "core",      1,   FALSE, FALSE,
  "life_exp",  "social", "social",      "core",      1,   FALSE, FALSE,
  "mys",       "social", "social",      "core",      1,   FALSE, FALSE,
  "epi_PMD",   "env",    "env_health",  "core",      1,   FALSE, TRUE,
  "epi_HFD",   "env",    "env_health",  "core",      1,   FALSE, TRUE,
  "epi_NOD",   "env",    "env_health",  "core",      1,   FALSE, TRUE,
  "epi_OZD",   "env",    "env_health",  "core",      1,   FALSE, TRUE,
  "epi_WASH",  "env",    "env_health",  "core",      1,   FALSE, TRUE,
  "epi_PBD",   "env",    "env_health",  "extended",  1,   FALSE, TRUE,
  "epi_TBN",   "env",    "env_eco",     "core",      1,   FALSE, TRUE,
  "epi_TKP",   "env",    "env_eco",     "core",      1,   FALSE, TRUE,
  "epi_TCL",   "env",    "env_eco",     "core",      1,   FALSE, TRUE,
  "epi_RLI",   "env",    "env_eco",     "core",      1,   FALSE, TRUE,
  "epi_PAR",   "env",    "env_eco",     "extended",  1,   FALSE, TRUE,
  "epi_PAC",   "env",    "env_eco",     "extended",  1,   FALSE, TRUE,
  "epi_KFL",   "env",    "env_eco",     "extended",  1,   FALSE, TRUE,
  "epi_SNM",   "env",    "env_agri",    "core",      1,   FALSE, TRUE,
  "epi_PSU",   "env",    "env_agri",    "core",      1,   FALSE, TRUE,
  "co2int",    "env",    "env_climate", "core",     -1,   TRUE,  FALSE,
  "re_elec",   "env",    "env_climate", "core",      1,   FALSE, FALSE
)
ENV_SUBS <- c("sub_env_health", "sub_env_eco", "sub_env_agri", "sub_env_climate")

dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------
save_tbl <- function(x, name) write_csv(x, file.path(OUT_DIR, "tables", paste0(name, ".csv")))
save_fig <- function(p, name, w = 8, h = 6) ggsave(file.path(OUT_DIR, "figures", paste0(name, ".png")), p, width = w, height = h, dpi = 200)
sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))

fill_short <- function(x, maxgap = MAX_GAP) {
  if (sum(!is.na(x)) < 2) return(x)
  x <- zoo::na.approx(x, maxgap = maxgap, na.rm = FALSE)
  zoo::na.locf(x, maxgap = maxgap, na.rm = FALSE)
}

 
dim_scores <- function(df, ind_u, prefix = "n_") {
  out <- df[, c("iso3", "year")]
  for (sb in unique(ind_u$sub)) {
    cols <- paste0(prefix, ind_u$var[ind_u$sub == sb])
    v <- rowMeans(df[cols], na.rm = TRUE); v[is.nan(v)] <- NA_real_
    out[[paste0("sub_", sb)]] <- v
  }
  for (d in unique(ind_u$dim)) {
    scols <- paste0("sub_", unique(ind_u$sub[ind_u$dim == d]))
    out[[d]] <- rowMeans(out[scols])
  }
  out
}
compose <- function(M, w = rep(1/3, 3), type = c("arith", "geo")) {
  type <- match.arg(type)
  if (type == "arith") as.vector(M %*% w) else exp(as.vector(log(pmax(M, 1)) %*% w))
}
tidy_hc3 <- function(m) {
  v <- sandwich::vcovHC(m, type = "HC3"); ct <- lmtest::coeftest(m, vcov. = v); ci <- lmtest::coefci(m, vcov. = v)
  tibble(term = rownames(ct), est = ct[, 1], se = ct[, 2], p = ct[, 4], lo = ci[, 1], hi = ci[, 2],
         n = nobs(m), r2 = summary(m)$r.squared)
}
 
fit_pc <- function(vars, data_fit, data_all) {
  X <- as.matrix(data_fit[vars]); ok <- complete.cases(X)
  pc <- prcomp(X[ok, , drop = FALSE], center = TRUE, scale. = TRUE)
  rot <- pc$rotation[, 1]; sgn <- if (sum(rot) < 0) -1 else 1
  sc_fit <- sgn * pc$x[, 1]; m <- mean(sc_fit); s <- sd(sc_fit)
  Xa <- as.matrix(data_all[vars]); oka <- complete.cases(Xa)
  out <- rep(NA_real_, nrow(Xa))
  out[oka] <- (sgn * predict(pc, Xa[oka, , drop = FALSE])[, 1] - m) / s
  list(score = out, loadings = tibble(var = vars, loading = sgn * rot), var_expl = summary(pc)$importance[2, 1])
}

# ---------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------
raw <- read_csv(INPUT, show_col_types = FALSE)
base_need <- c("iso3","name","group","year","pop","gdp_pc","gdp_ppp","pm25","epi","co2_mt","primary_energy",
               setdiff(unique(c(inst_core, cap_vars, dem_vars)), "vdem_corr_inv"), "vdem_corr")
need <- unique(c(base_need, ind$var[ind$tier == "core"]))
miss <- setdiff(need, names(raw))
if (length(miss)) stop("Missing columns in the panel: ", paste(miss, collapse = ", "))
message("Panel: ", nrow(raw), " rows, ", n_distinct(raw$iso3), " countries, years ", min(raw$year), "–", max(raw$year))

panel <- raw %>% arrange(iso3, year) %>% mutate(vdem_corr_inv = 1 - vdem_corr)
ind_all <- ind %>% filter(var %in% names(panel)) %>%
  filter(map_lgl(var, ~ any(!is.na(panel[[.x]][panel$year %in% EST_WIN]))))
ind_use <- ind_all %>% filter(tier == "core" | USE_EXTENDED)

# Fill gaps <= MAX_GAP within the estimation window (SD indicators only)
INPUT0 <- sub("\\.csv$", "_unfilled.csv", INPUT)
raw0 <- if (file.exists(INPUT0)) read_csv(INPUT0, show_col_types = FALSE) else NULL
if (is.null(raw0) || !all(ind_all$var %in% names(raw0))) {
  warning("Missing file: ", INPUT0, " (or it lacks the required columns): D1 will be calculated from the already-filled panel")
  raw0 <- raw
}
na_before <- map_int(ind_all$var, ~ sum(is.na(raw0[[.x]][raw0$year %in% EST_WIN])))
panel <- panel %>% group_by(iso3) %>%
  mutate(across(all_of(ind_all$var), ~ { x <- .x; w <- year %in% EST_WIN; x[w] <- fill_short(x[w]); x })) %>% ungroup()
na_after <- map_int(ind_all$var, ~ sum(is.na(panel[[.x]][panel$year %in% EST_WIN])))
miss_by <- map_chr(ind_all$var, function(v) {
  w <- panel[panel$year %in% EST_WIN, c("iso3", v)]
  k <- tapply(is.na(w[[v]]), w$iso3, sum); k <- k[k > 0]
  if (length(k)) paste0(names(k), "(", k, ")", collapse = " ") else "" })
dq <- tibble(var = ind_all$var, tier = ind_all$tier, na_before_fill = na_before, na_after_fill = na_after,
             n_cells_filled = na_before - na_after, countries_still_missing = miss_by)
save_tbl(dq, "D1_data_quality_sd_indicators")

# ---------------------------------------------------------------------
# 2. INSTITUTIONS: overall index and two factors (capacity / accountability), 2009–2025
# ---------------------------------------------------------------------
win_rows <- panel$year %in% EST_WIN
f_all <- fit_pc(inst_core, panel[win_rows, ], panel); panel$inst_z <- f_all$score
f_cap <- fit_pc(cap_vars,  panel[win_rows, ], panel); panel$cap_z  <- f_cap$score
f_dem <- fit_pc(dem_vars,  panel[win_rows, ], panel); panel$dem_z  <- f_dem$score
lg <- log(panel$gdp_pc)
panel$lgdp_z <- (lg - mean(lg[win_rows], na.rm = TRUE)) / sd(lg[win_rows], na.rm = TRUE)

save_tbl(bind_rows(mutate(f_all$loadings, factor = "inst_overall", var_expl_PC1 = f_all$var_expl),
                   mutate(f_cap$loadings,  factor = "capacity",      var_expl_PC1 = f_cap$var_expl),
                   mutate(f_dem$loadings,  factor = "accountability", var_expl_PC1 = f_dem$var_expl)),
         "I0_institutional_pca_loadings")
corr_inst <- cor(panel[win_rows, unique(c(inst_core, cap_vars, dem_vars))], use = "pairwise.complete.obs")
save_tbl(as_tibble(corr_inst, rownames = "var"), "I0b_institutional_corr")
cd_2023 <- panel %>% filter(year == REF_YEAR)
message(sprintf("Institutions: overall PC1 %.0f%%, capacity %.0f%%, accountability %.0f%%; corr(capacity, accountability) in %d = %.2f",
                100 * f_all$var_expl, 100 * f_cap$var_expl, 100 * f_dem$var_expl, REF_YEAR, cor(cd_2023$cap_z, cd_2023$dem_z)))

# ---------------------------------------------------------------------
# 3. SD INDEX (estimation window): normalization, sub-blocks, dimensions, composites
# ---------------------------------------------------------------------
idx <- panel %>% filter(year %in% EST_WIN)
gp <- map_dfr(ind_all$var[!ind_all$scored], function(v) {
  x <- idx[[v]]; if (ind_all$log[ind_all$var == v]) x <- log(x)
  tibble(var = v, lo = min(x, na.rm = TRUE), hi = max(x, na.rm = TRUE)) })
for (i in seq_len(nrow(ind_all))) {
  v <- ind_all$var[i]; x <- idx[[v]]; if (ind_all$log[i]) x <- log(x)
  if (ind_all$scored[i]) { n <- x } else { g <- gp[gp$var == v, ]; n <- (x - g$lo) / (g$hi - g$lo) * 100 }
  idx[[paste0("n_", v)]] <- if (ind_all$dir[i] == -1) 100 - n else n
  idx[[paste0("z_", v)]] <- as.numeric(scale(x)) * ind_all$dir[i]            # alternative z-normalization
}
idx <- idx %>% left_join(dim_scores(idx, ind_use), by = c("iso3", "year"))
M <- as.matrix(idx[c("econ", "social", "env")])
idx$sd_arith <- compose(M, type = "arith")
idx$sd_geo   <- compose(M, type = "geo")
idx$sd_nogdp <- rowMeans(idx[c("social", "env")])
cc <- complete.cases(M); pc_sd <- prcomp(M[cc, ], center = TRUE, scale. = TRUE)
s <- pc_sd$x[, 1]; if (cor(s, rowMeans(M[cc, ])) < 0) s <- -s
idx$sd_pca <- NA_real_; idx$sd_pca[cc] <- (s - min(s)) / (max(s) - min(s)) * 100
dz <- dim_scores(idx, ind_use, "z_"); idx$sd_z <- rowMeans(dz[c("econ", "social", "env")])
idx$co2int_log <- log(idx$co2int); idx$pm25_log <- log(idx$pm25)
# Variant using the overall EPI score instead of green sub-blocks (REF_YEAR only; for validation)
ref <- idx$year == REF_YEAR & !is.na(idx$epi)
idx$env_epi <- NA_real_
idx$env_epi[ref] <- (idx$epi[ref] - min(idx$epi[ref])) / (max(idx$epi[ref]) - min(idx$epi[ref])) * 100
idx$sd_epi <- rowMeans(cbind(idx$econ, idx$social, idx$env_epi))
save_tbl(idx %>% select(iso3, name, group, year, econ, social, env, all_of(ENV_SUBS), sd_arith, sd_geo, sd_nogdp, sd_pca, sd_z,
                        inst_z, cap_z, dem_z), "S0_indices_panel_2012_2023")

# ---------------------------------------------------------------------
# 4. CROSS-SECTION (REF_YEAR): H1 gaps, H2 institutions ~ SD, H3 EU benchmark
# ---------------------------------------------------------------------
xs <- idx %>% filter(year == REF_YEAR) %>% mutate(cand = as.integer(group == "candidate")) %>%
  drop_na(sd_arith, sd_nogdp, inst_z, cap_z, dem_z, lgdp_z)
eu <- filter(xs, group == "EU")
gap_vars <- c("econ", "social", "env", ENV_SUBS, "sd_arith", "inst_z", "cap_z", "dem_z")
eu_w <- eu %>% summarise(across(all_of(gap_vars), ~ weighted.mean(.x, pop, na.rm = TRUE)))
q <- quantile(eu$sd_arith, c(.25, .5, .75))
T1 <- xs %>% arrange(desc(sd_arith)) %>% mutate(rank = row_number(),
        eu_quartile = cut(sd_arith, c(-Inf, q, Inf), labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))) %>%
  select(rank, iso3, name, group, econ, social, env, all_of(ENV_SUBS), sd_arith, sd_geo, inst_z, cap_z, dem_z, eu_quartile)
save_tbl(T1, "T1_ranking")
cands <- filter(xs, group == "candidate")
T2 <- cands %>% transmute(iso3, name,
        across(all_of(c("econ", "social", "env", ENV_SUBS)), ~ .x - eu_w[[cur_column()]], .names = "gap_vs_EUmean_{.col}"),
        across(all_of(c("econ", "social", "env", ENV_SUBS)), ~ .x - min(eu[[cur_column()]]), .names = "gap_vs_EUmin_{.col}"))
dimgaps <- as.matrix(T2[paste0("gap_vs_EUmean_", c("econ", "social", "env"))])
T2$binding_dimension <- c("econ", "social", "env")[max.col(-dimgaps, "first")]
save_tbl(T2, "T2_gaps_candidates")

m1  <- lm(sd_nogdp ~ inst_z + lgdp_z, data = xs)
m2  <- lm(sd_nogdp ~ inst_z * cand + lgdp_z, data = xs)
m3  <- lm(sd_nogdp ~ cap_z + dem_z + lgdp_z, data = xs)
m4  <- lm(sd_nogdp ~ (cap_z + dem_z) * cand + lgdp_z, data = xs)
mods <- list(M1_inst = m1, M2_inst_x_cand = m2, M3_two_factors = m3, M4_two_factors_x_cand = m4,
             M5_env = lm(env ~ inst_z + lgdp_z, xs), M6_social = lm(social ~ inst_z + lgdp_z, xs),
             M7_env_2f = lm(env ~ cap_z + dem_z + lgdp_z, xs), M8_social_2f = lm(social ~ cap_z + dem_z + lgdp_z, xs))
for (sb in ENV_SUBS) mods[[paste0("M9_", sb)]] <- lm(as.formula(paste(sb, "~ inst_z + lgdp_z")), xs)
for (sb in ENV_SUBS) mods[[paste0("M10_", sb, "_2f")]] <- lm(as.formula(paste(sb, "~ cap_z + dem_z + lgdp_z")), xs)
if (any(!is.na(xs$sd_epi))) mods$M11_sd_nogdp_EPIvariant <- lm(rowMeans(cbind(social, env_epi)) ~ inst_z + lgdp_z, xs)
T3 <- bind_rows(map(mods, tidy_hc3), .id = "model"); save_tbl(T3, "T3_cross_section_models")
pcor <- function(y, x, z) cor(resid(lm(y ~ z)), resid(lm(x ~ z)), method = "spearman")
save_tbl(tibble(pair = c("sd_nogdp~inst|lgdp", "sd_nogdp~cap|lgdp", "sd_nogdp~dem|lgdp", "env~inst|lgdp", "social~inst|lgdp"),
                partial_spearman = c(pcor(xs$sd_nogdp, xs$inst_z, xs$lgdp_z), pcor(xs$sd_nogdp, xs$cap_z, xs$lgdp_z),
                                     pcor(xs$sd_nogdp, xs$dem_z, xs$lgdp_z), pcor(xs$env, xs$inst_z, xs$lgdp_z),
                                     pcor(xs$social, xs$inst_z, xs$lgdp_z))), "T3b_partial_spearman")

# Bootstrap (by country) and jackknife
boot <- t(replicate(N_BOOT, { d <- xs[sample(nrow(xs), replace = TRUE), ]
  c(inst = unname(coef(lm(sd_nogdp ~ inst_z + lgdp_z, d))["inst_z"]),
    cap = unname(coef(lm(sd_nogdp ~ cap_z + dem_z + lgdp_z, d))["cap_z"]),
    dem = unname(coef(lm(sd_nogdp ~ cap_z + dem_z + lgdp_z, d))["dem_z"])) }))
T3c <- tibble(coef = colnames(boot), est = c(coef(m1)["inst_z"], coef(m3)["cap_z"], coef(m3)["dem_z"]),
              boot_lo = apply(boot, 2, quantile, .025, na.rm = TRUE), boot_hi = apply(boot, 2, quantile, .975, na.rm = TRUE),
              share_positive = apply(boot, 2, function(z) mean(z > 0, na.rm = TRUE)))
save_tbl(T3c, "T3c_bootstrap_ci")
jk <- map_dfr(xs$iso3, function(i) { d <- filter(xs, iso3 != i)
  tibble(dropped = i, inst = coef(lm(sd_nogdp ~ inst_z + lgdp_z, d))["inst_z"],
         cap = coef(lm(sd_nogdp ~ cap_z + dem_z + lgdp_z, d))["cap_z"], dem = coef(lm(sd_nogdp ~ cap_z + dem_z + lgdp_z, d))["dem_z"]) })
jk$cooks_M1 <- as.numeric(cooks.distance(m1))[match(jk$dropped, xs$iso3)]
save_tbl(arrange(jk, desc(cooks_M1)), "T3d_jackknife")

# H3: EU-27 benchmark -> predictions for candidate countries
bench <- function(y, xv) {
  m <- lm(as.formula(paste(y, "~", paste(c(xv, "lgdp_z"), collapse = "+"))), data = eu)
  pr <- predict(m, newdata = cands, interval = "prediction")
  tibble(outcome = y, spec = paste(xv, collapse = "+"), iso3 = cands$iso3, name = cands$name, actual = cands[[y]],
         fit = pr[, 1], lwr = pr[, 2], upr = pr[, 3], gap = actual - fit, outside_95PI = actual < lwr | actual > upr, r2_EU = summary(m)$r.squared)
}
T4 <- bind_rows(map_dfr(c("sd_nogdp", "env", "social"), bench, xv = "inst_z"),
                map_dfr(c("sd_nogdp", "env", "social"), bench, xv = c("cap_z", "dem_z")))
save_tbl(T4, "T4_eu_benchmark_gaps")

# ---------------------------------------------------------------------
# 5. PANEL 2012–2023: institutions (1- and 3-year lags, available from 2009) -> SD and sub-blocks
# ---------------------------------------------------------------------
lagcols <- c("inst_z", "cap_z", "dem_z", "lgdp_z")
pl <- panel %>% arrange(iso3, year) %>% group_by(iso3) %>%
  mutate(across(all_of(lagcols), list(l1 = ~ lag(.x, 1), l3 = ~ lag(.x, 3)), .names = "{.col}_{.fn}")) %>% ungroup()
pan <- pl %>% filter(year %in% EST_WIN) %>%
  select(iso3, year, ends_with("_l1"), ends_with("_l3")) %>%
  left_join(select(idx, iso3, year, group, sd_arith, sd_nogdp, social, env, all_of(ENV_SUBS)), by = c("iso3", "year")) %>%
  mutate(cand = as.integer(group == "candidate"))
run_panel <- function(y, xs_names, L, fe = c("twoway", "year"), interact = FALSE) {
  fe <- match.arg(fe)
  xv <- paste0(xs_names, "_l", L); gv <- paste0("lgdp_z_l", L)
  rhs <- c(xv, gv); if (interact) rhs <- c(rhs, paste0(xv, ":cand"))
  fe_part <- if (fe == "twoway") "factor(iso3) + factor(year)" else "factor(year) + cand"
  d <- pan %>% drop_na(all_of(c(y, xv, gv)))
  m <- lm(as.formula(paste(y, "~", paste(rhs, collapse = " + "), "+", fe_part)), data = d)
  ct <- coeftest(m, vcov. = vcovCL(m, cluster = d$iso3))
  keep <- which(grepl(paste0("^(", paste(xv, collapse = "|"), ")"), rownames(ct)))
  tibble(outcome = y, x = paste(xs_names, collapse = "+"), lag = L, fe = fe, interact = interact, term = rownames(ct)[keep],
         est = ct[keep, 1], se = ct[keep, 2], p = ct[keep, 4], mde80 = 2.8 * ct[keep, 2], n_obs = nobs(m), n_countries = n_distinct(d$iso3))
}
outs <- c("sd_nogdp", "social", "env", ENV_SUBS)
specs <- list(inst = "inst_z", two_factor = c("cap_z", "dem_z"))
T5 <- map_dfr(outs, function(y) map_dfr(specs, function(xn) map_dfr(c(1, 3), function(L)
        map_dfr(c("twoway", "year"), function(fe) map_dfr(c(FALSE, TRUE), function(it) run_panel(y, xn, L, fe, it))))))
save_tbl(T5, "T5_panel_models")

# ---------------------------------------------------------------------
# 6. DYNAMICS: decoupling, Kaya decomposition (LMDI), convergence
# ---------------------------------------------------------------------
per_avg <- function(df, yrs, vars) df %>% filter(year %in% yrs) %>% group_by(iso3, name, group) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
dec_calc <- function(y0, y1, label) {
  a <- per_avg(panel, y0, c("gdp_ppp", "co2_mt")); b <- per_avg(panel, y1, c("gdp_ppp", "co2_mt"))
  left_join(a, b, by = c("iso3", "name", "group"), suffix = c("_0", "_1")) %>%
    mutate(g_gdp = gdp_ppp_1 / gdp_ppp_0 - 1, g_co2 = co2_mt_1 / co2_mt_0 - 1,
           decoupling_factor = 1 - (co2_mt_1 / co2_mt_0) / (gdp_ppp_1 / gdp_ppp_0),
           class = case_when(g_gdp > 0 & g_co2 < 0 ~ "absolute decoupling",
                             g_gdp > 0 & g_co2 >= 0 & g_co2 < g_gdp ~ "relative decoupling",
                             g_gdp > 0 & g_co2 >= g_gdp ~ "no decoupling",
                             g_gdp <= 0 & g_co2 < 0 ~ "recession (both fall)", TRUE ~ "other"), window = label)
}
DEC <- bind_rows(dec_calc(P0, P1, "main 2012-14 -> 2021-23"), dec_calc(P0, 2017:2019, "pre-COVID 2012-14 -> 2017-19"),
                 dec_calc(P0, 2022:2024, "to 2024: 2012-14 -> 2022-24"))
save_tbl(DEC, "B1_decoupling")
DECs <- DEC %>% group_by(window, group, class) %>% summarise(n = n(), .groups = "drop") %>% arrange(window, group, desc(n))
save_tbl(DECs, "B1b_decoupling_class_counts")
save_tbl(DEC %>% filter(iso3 != "UKR") %>% group_by(window, group) %>%
           summarise(mean_decoupling_factor_exUKR = mean(decoupling_factor), median_g_gdp = median(g_gdp), median_g_co2 = median(g_co2), .groups = "drop"),
         "B1c_decoupling_means_exUKR")

# LMDI-I (additive) Kaya: CO2 = P * (GDP/P) * (E/GDP) * (CO2/E)
Lm <- function(a, b) ifelse(abs(a - b) < 1e-9, a, (a - b) / (log(a) - log(b)))
ka <- per_avg(panel, P0, c("pop", "gdp_ppp", "primary_energy", "co2_mt"))
kb <- per_avg(panel, P1, c("pop", "gdp_ppp", "primary_energy", "co2_mt"))
KAYA <- left_join(ka, kb, by = c("iso3", "name", "group"), suffix = c("_0", "_1")) %>%
  mutate(L = Lm(co2_mt_1, co2_mt_0), d_co2_Mt = co2_mt_1 - co2_mt_0,
         eff_population = L * log(pop_1 / pop_0),
         eff_affluence  = L * log((gdp_ppp_1 / pop_1) / (gdp_ppp_0 / pop_0)),
         eff_energy_intensity = L * log((primary_energy_1 / gdp_ppp_1) / (primary_energy_0 / gdp_ppp_0)),
         eff_carbon_intensity_of_energy = L * log((co2_mt_1 / primary_energy_1) / (co2_mt_0 / primary_energy_0)),
         check_sum_error = d_co2_Mt - (eff_population + eff_affluence + eff_energy_intensity + eff_carbon_intensity_of_energy)) %>%
  mutate(across(starts_with("eff_"), ~ 100 * .x / co2_mt_0, .names = "pct_{.col}"), pct_total = 100 * d_co2_Mt / co2_mt_0) %>%
  select(iso3, name, group, co2_Mt_0 = co2_mt_0, co2_Mt_1 = co2_mt_1, d_co2_Mt, starts_with("eff_"), starts_with("pct_"), check_sum_error)
save_tbl(KAYA, "B2_kaya_lmdi")
save_tbl(KAYA %>% group_by(group) %>% summarise(across(starts_with("pct_"), list(median = median)), .groups = "drop"), "B2b_kaya_group_medians")

# Convergence: beta (conditional and unconditional) and sigma
idx$inst_ref <- idx$inst_z
conv_vars <- c("sd_arith", "econ", "social", "env", ENV_SUBS, "co2int_log", "pm25_log", "re_elec", "cap_z", "dem_z")
Tspan <- mean(P1) - mean(P0)
conv_one <- function(v) {
  d <- idx %>% mutate(x = .data[[v]]) %>% filter(year %in% c(P0, P1)) %>%
    group_by(iso3, name, group) %>%
    summarise(x0 = mean(x[year %in% P0], na.rm = TRUE), x1 = mean(x[year %in% P1], na.rm = TRUE),
              inst0 = mean(inst_z[year %in% P0], na.rm = TRUE), lg0 = mean(lgdp_z[year %in% P0], na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = (x1 - x0) / Tspan, cand = as.integer(group == "candidate")) %>% drop_na(rate, x0, inst0, lg0)
  f <- list(unconditional = lm(rate ~ x0, d), conditional = lm(rate ~ x0 + inst0 + lg0, d), by_group = lm(rate ~ x0 * cand, d))
  res <- imap_dfr(f, function(m, nm) tidy_hc3(m) %>% filter(grepl("^x0", term)) %>% mutate(spec = nm))
  res %>% mutate(var = v, T_years = Tspan,
                 lambda = ifelse(term == "x0" & est < 0 & 1 + est * Tspan > 0, -log(1 + est * Tspan) / Tspan, NA_real_),
                 half_life_years = log(2) / lambda)
}
CONV <- map_dfr(conv_vars, conv_one) %>% select(var, spec, term, est, se, p, lo, hi, lambda, half_life_years, n, r2, T_years)
save_tbl(CONV, "B3_beta_convergence")
SIG <- idx %>% select(iso3, group, year, all_of(conv_vars)) %>% pivot_longer(all_of(conv_vars), names_to = "var", values_to = "x") %>%
  group_by(var, year) %>% summarise(sd_all = sd(x, na.rm = TRUE), sd_EU = sd(x[group == "EU"], na.rm = TRUE),
                                    sd_cand = sd(x[group == "candidate"], na.rm = TRUE), mean_EU = mean(x[group == "EU"], na.rm = TRUE),
                                    mean_cand = mean(x[group == "candidate"], na.rm = TRUE), gap_cand_minus_EU = mean_cand - mean_EU, .groups = "drop")
save_tbl(SIG, "B4_sigma_convergence_and_gap")
SIGt <- SIG %>% group_by(var) %>% summarise(
  slope_sd_all = coef(lm(sd_all ~ year))[2], p_sd_all = summary(lm(sd_all ~ year))$coefficients[2, 4],
  slope_gap = coef(lm(gap_cand_minus_EU ~ year))[2], p_gap = summary(lm(gap_cand_minus_EU ~ year))$coefficients[2, 4], .groups = "drop")
save_tbl(SIGt, "B4b_sigma_trends")

# ---------------------------------------------------------------------
# 7. INSTITUTIONAL MONITORING 2009–2025 + source triangulation (Georgia, Moldova)
# ---------------------------------------------------------------------
mon <- panel %>% filter(year %in% INST_WIN) %>% select(iso3, name, group, year, inst_z, cap_z, dem_z)
save_tbl(mon, "I1_institutional_trajectories_2009_2025")
chg <- mon %>% select(iso3, name, group, year, cap_z, dem_z) %>% filter(year %in% c(2023, 2024, 2025)) %>%
  pivot_wider(names_from = year, values_from = c(cap_z, dem_z), names_sep = "_") %>%
  mutate(d_cap_2023_2025 = cap_z_2025 - cap_z_2023, d_dem_2023_2025 = dem_z_2025 - dem_z_2023) %>% arrange(d_dem_2023_2025)
save_tbl(chg, "I2_changes_2023_2025")
brk <- panel %>% filter(year %in% 2015:2025) %>% group_by(iso3, name, group) %>% group_modify(function(d, k) {
  map_dfr(c("cap_z", "dem_z"), function(v) { tr <- filter(d, year <= 2023)
    m <- lm(reformulate("year", v), tr)
    tibble(factor = v, resid_2024 = d[[v]][d$year == 2024] - predict(m, data.frame(year = 2024)),
           resid_2025 = d[[v]][d$year == 2025] - predict(m, data.frame(year = 2025))) }) }) %>% ungroup() %>%
  group_by(factor) %>% mutate(pctile_resid_2025 = percent_rank(resid_2025)) %>% ungroup() %>% arrange(factor, resid_2025)
save_tbl(brk, "I3_break_vs_2015_2023_trend")

tri <- tribble(~var, ~panel_type, ~source,
  "wgi_va", "Democracy / accountability", "WGI voice & accountability", "vdem_libdem", "Democracy / accountability", "V-Dem liberal democracy",
  "vdem_polyarchy", "Democracy / accountability", "V-Dem polyarchy", "fh_total", "Democracy / accountability", "Freedom House total",
  "EIU_Index", "Democracy / accountability", "EIU Democracy Index",
  "wgi_cc", "State capacity", "WGI control of corruption", "wgi_ge", "State capacity", "WGI government effectiveness",
  "vdem_rule", "State capacity", "V-Dem rule of law") %>% filter(var %in% names(panel))
# Guard against carried-forward values: if >= 90% of countries have exactly the same value as in the previous year, the year is treated as NOT genuine
CARRY_SHARE <- 0.9
carry_chk <- panel %>% arrange(iso3, year) %>% group_by(iso3) %>%
  mutate(across(all_of(tri$var), ~ .x == lag(.x), .names = "same_{.col}")) %>% ungroup() %>%
  filter(year %in% INST_WIN) %>% group_by(year) %>%
  summarise(across(starts_with("same_"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-year, names_to = "var", names_prefix = "same_", values_to = "share_identical_to_prev_year") %>%
  mutate(flagged_as_carried = !is.na(share_identical_to_prev_year) & share_identical_to_prev_year >= CARRY_SHARE)
save_tbl(carry_chk, "D2_carry_forward_check_triangulation_sources")
if (any(carry_chk$flagged_as_carried)) message("WARNING: carried-forward values (excluded from the triangulation figure and table): ",
  paste(paste(carry_chk$var[carry_chk$flagged_as_carried], carry_chk$year[carry_chk$flagged_as_carried]), collapse = "; "))
tri_long <- panel %>% filter(year %in% INST_WIN) %>% select(iso3, name, group, year, all_of(tri$var)) %>%
  pivot_longer(all_of(tri$var), names_to = "var", values_to = "value") %>%
  left_join(select(carry_chk, year, var, flagged_as_carried), by = c("year", "var")) %>%
  mutate(value = ifelse(coalesce(flagged_as_carried, FALSE), NA_real_, value)) %>% left_join(tri, by = "var") %>%
  group_by(var) %>% mutate(lo = min(value[year %in% EST_WIN], na.rm = TRUE), hi = max(value[year %in% EST_WIN], na.rm = TRUE),
                           value01 = (value - lo) / (hi - lo) * 100) %>% ungroup()
tri_tab <- tri_long %>% filter(iso3 %in% c("GEO", "MDA"), year %in% c(2022, 2023, 2024, 2025)) %>%
  select(iso3, source, panel_type, year, value01) %>% pivot_wider(names_from = year, values_from = value01, names_prefix = "y") %>%
  mutate(last_real_year = ifelse(!is.na(y2025), 2025, ifelse(!is.na(y2024), 2024, 2023)), change_2023_to_last = coalesce(y2025, y2024) - y2023, direction = sign(change_2023_to_last))
save_tbl(tri_tab, "I4_triangulation_GEO_MDA")
save_tbl(tri_long %>% filter(iso3 %in% c("GEO", "MDA")) %>% select(iso3, source, panel_type, year, value, value01), "I4b_triangulation_series_GEO_MDA")

# ---------------------------------------------------------------------
# 8. ROBUSTNESS
# ---------------------------------------------------------------------
# R1: Aggregation method
aggr <- xs %>% select(any_of(c("sd_arith", "sd_geo", "sd_pca", "sd_z", "sd_epi")))
save_tbl(as_tibble(cor(aggr, method = "spearman", use = "pairwise.complete.obs"), rownames = "method"), "R1_aggregation_spearman")
# R2: Random weights (3 dimensions and 4 green sub-blocks), ranks
D2 <- as.matrix(xs[c("econ", "social")]); S4 <- as.matrix(xs[ENV_SUBS])
W3 <- matrix(rgamma(3 * N_MC, 1), ncol = 3); W3 <- W3 / rowSums(W3)
W4 <- matrix(rgamma(4 * N_MC, 1), ncol = 4); W4 <- W4 / rowSums(W4)
envmc <- S4 %*% t(W4)
sc <- sweep(matrix(D2[, 1], nrow(D2), N_MC), 2, W3[, 1], "*") + sweep(matrix(D2[, 2], nrow(D2), N_MC), 2, W3[, 2], "*") + sweep(envmc, 2, W3[, 3], "*")
Rk <- apply(-sc, 2, rank, ties.method = "min")
R2 <- tibble(iso3 = xs$iso3, name = xs$name, group = xs$group, rank_base = rank(-xs$sd_arith, ties.method = "min"),
             rank_p05 = apply(Rk, 1, quantile, .05), rank_median = apply(Rk, 1, median), rank_p95 = apply(Rk, 1, quantile, .95)) %>% arrange(rank_base)
save_tbl(R2, "R2_rank_intervals_random_weights")
# R3: Indicator / sub-block composition scenarios
sd_variant <- function(df_ref, drop = NULL, drop_sub = NULL, add_ext = FALSE) {
  iv <- (if (add_ext) ind_all else ind_use) %>% filter(!var %in% drop, !sub %in% drop_sub)
  ds <- dim_scores(df_ref, iv)
  tibble(iso3 = ds$iso3, sd_arith = compose(as.matrix(ds[c("econ", "social", "env")])), sd_nogdp = rowMeans(ds[c("social", "env")]))
}
ref_df <- filter(idx, year == REF_YEAR)
scen <- c(setNames(map(ind_use$var, ~ list(drop = .x)), paste0("drop_", ind_use$var)),
          setNames(map(c("env_health", "env_eco", "env_agri", "env_climate"), ~ list(drop_sub = .x)), paste0("dropsub_", c("env_health", "env_eco", "env_agri", "env_climate"))),
          list(only_OWID_climate_env = list(drop_sub = c("env_health", "env_eco", "env_agri")),
               only_EPI_env = list(drop_sub = "env_climate"),
               with_extended_indicators = list(add_ext = TRUE)))
base_rank <- rank(-filter(xs)$sd_arith)
R3 <- imap_dfr(scen, function(a, nm) {
  v <- do.call(sd_variant, c(list(df_ref = ref_df), a)) %>% inner_join(select(xs, iso3, inst_z, lgdp_z), by = "iso3")
  m <- lm(sd_nogdp ~ inst_z + lgdp_z, v); ci <- suppressWarnings(coefci(m, vcov. = vcovHC(m, "HC3")))["inst_z", ]
  tibble(scenario = nm, spearman_vs_base = sp(rank(-v$sd_arith), base_rank[match(v$iso3, xs$iso3)]),
         beta_inst = coef(m)["inst_z"], lo = ci[1], hi = ci[2], r2 = summary(m)$r.squared) })
save_tbl(R3, "R3_composition_scenarios")
cand_ranks <- imap_dfr(scen, function(a, nm) { v <- do.call(sd_variant, c(list(df_ref = ref_df), a)); v$rank <- rank(-v$sd_arith, ties.method = "min")
  v %>% filter(iso3 %in% cands$iso3) %>% transmute(scenario = nm, iso3, rank) }) %>% pivot_wider(names_from = iso3, values_from = rank)
save_tbl(cand_ranks, "R3b_candidate_ranks_by_scenario")
# R4: Excluding countries
excl <- c("UKR", "TUR", "GEO", "ALB", "MDA", "MNE", "SRB", "MKD", "BIH", "IRL", "LUX", "CYP", "MLT")
R4 <- map_dfr(excl, function(x) { d <- filter(xs, iso3 != x); m <- lm(sd_nogdp ~ inst_z + lgdp_z, d); ci <- coefci(m, vcov. = vcovHC(m, "HC3"))["inst_z", ]
  m2 <- lm(sd_nogdp ~ cap_z + dem_z + lgdp_z, d)
  tibble(excluded = x, beta_inst = coef(m)["inst_z"], lo = ci[1], hi = ci[2], beta_cap = coef(m2)["cap_z"], beta_dem = coef(m2)["dem_z"]) }) %>%
  bind_rows(tibble(excluded = "IRL+LUX+CYP+MLT",
                   beta_inst = coef(lm(sd_nogdp ~ inst_z + lgdp_z, filter(xs, !iso3 %in% c("IRL", "LUX", "CYP", "MLT"))))["inst_z"]),
            tibble(excluded = "UKR+ALB", beta_inst = coef(lm(sd_nogdp ~ inst_z + lgdp_z, filter(xs, !iso3 %in% c("UKR", "ALB"))))["inst_z"]))
save_tbl(R4, "R4_exclude_countries")
# R5: Cross-sectional year
R5 <- map_dfr(c(2018, 2020, 2021, 2022, 2023), function(y) { d <- idx %>% filter(year == y) %>% drop_na(sd_arith, sd_nogdp, inst_z, lgdp_z)
  m <- lm(sd_nogdp ~ inst_z + lgdp_z, d); ci <- coefci(m, vcov. = vcovHC(m, "HC3"))["inst_z", ]
  tibble(ref_year = y, beta_inst = coef(m)["inst_z"], lo = ci[1], hi = ci[2], r2 = summary(m)$r.squared,
         spearman_rank_vs_2023 = sp(d$sd_arith, xs$sd_arith[match(d$iso3, xs$iso3)])) })
save_tbl(R5, "R5_reference_year_sensitivity")
# R6: Validation of the green dimension against the overall EPI 2026 score
val <- xs %>% filter(!is.na(epi)) %>% summarise(across(c(env, all_of(ENV_SUBS), sd_arith), ~ sp(.x, epi), .names = "rho_{.col}_vs_EPI"))
save_tbl(val, "R6_validation_vs_EPI_overall")

# ---------------------------------------------------------------------
# 9. CLUSTERS (36 countries; Ward) and dimension correlations
# ---------------------------------------------------------------------
feat <- xs %>% select(iso3, econ, social, env, cap_z, dem_z) %>% column_to_rownames("iso3")
Z <- scale(feat); hc <- hclust(dist(Z), method = "ward.D2")
silh <- map_dbl(2:8, ~ mean(silhouette(cutree(hc, .x), dist(Z))[, 3])); k_best <- which.max(silh) + 1
save_tbl(tibble(k = 2:8, mean_silhouette = silh), "C1_silhouette")
CL <- tibble(iso3 = rownames(feat), cluster_best = cutree(hc, k_best), cluster_k3 = cutree(hc, 3), cluster_k4 = cutree(hc, 4)) %>%
  left_join(select(xs, iso3, name, group), by = "iso3") %>% arrange(cluster_best, group)
save_tbl(CL, "C2_clusters")
save_tbl(as_tibble(cor(xs[c("econ", "social", "env", ENV_SUBS, "inst_z", "cap_z", "dem_z")], method = "spearman"), rownames = "var"), "C3_dimension_correlations")

# ---------------------------------------------------------------------
# 10. FIGURES
# ---------------------------------------------------------------------
th <- theme_minimal(base_size = 11)
col_g <- c(EU = "#4C78A8", candidate = "#E45756")
F1 <- ggplot(xs, aes(inst_z, sd_nogdp, colour = group)) +
  geom_smooth(data = eu, method = "lm", se = TRUE, colour = "grey40", fill = "grey85", formula = y ~ x) +
  geom_point(size = 2) + geom_text_repel(aes(label = iso3), size = 3, max.overlaps = 40, show.legend = FALSE) +
  scale_colour_manual(values = col_g) + th +
  labs(x = "Institutional index (z, PC1)", y = "Social + green SD score (0-100)", colour = NULL,
       title = "Institutions and sustainable development", subtitle = "Grey line: model fitted on EU-27 only")
save_fig(F1, "F1_institutions_vs_SD")

pos <- mon %>% select(iso3, name, group, year, cap_z, dem_z) %>% filter(year %in% c(2023, 2025)) %>% pivot_wider(names_from = year, values_from = c(cap_z, dem_z), names_sep = "_")
F2 <- ggplot() +
  geom_point(data = filter(pos, group == "EU"), aes(dem_z_2023, cap_z_2023), colour = "grey70", size = 1.8) +
  geom_segment(data = filter(pos, group == "candidate"), aes(x = dem_z_2023, y = cap_z_2023, xend = dem_z_2025, yend = cap_z_2025, colour = iso3),
               arrow = arrow(length = unit(.18, "cm")), linewidth = .7, show.legend = FALSE) +
  geom_text_repel(data = filter(pos, group == "candidate"), aes(dem_z_2025, cap_z_2025, label = iso3, colour = iso3), size = 3.2, show.legend = FALSE) +
  th + labs(x = "Democratic accountability factor (z)", y = "State capacity factor (z)",
            title = "Institutional position of candidate countries, 2023 -> 2025", subtitle = "Grey: EU-27 in 2023; arrows: candidates from 2023 to 2025")
save_fig(F2, "F2_capacity_vs_accountability_2023_2025")

ev <- events %>% left_join(distinct(panel, iso3, name), by = "iso3")
cap_txt <- paste0("Dashed lines: ", paste(paste0(ev$iso3, " ", ev$year, ": ", ev$label), collapse = "; "))
F3 <- ggplot(filter(tri_long, iso3 %in% c("GEO", "MDA")), aes(year, value01, colour = source)) +
  geom_vline(data = ev, aes(xintercept = year), linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = .8, na.rm = TRUE) + geom_point(size = 1.2, na.rm = TRUE) +
  facet_grid(panel_type ~ name) + th + theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 3)) +
  labs(x = NULL, y = "Score rescaled to 0-100 (2012-2023 pooled min-max)", colour = NULL,
       title = "Institutional trajectories of Georgia and Moldova, 2009-2025", subtitle = "Triangulation across independent sources", caption = cap_txt)
save_fig(F3, "F3_trajectories_GEO_MDA_triangulation", 10, 7)

F4 <- ggplot(R2, aes(reorder(iso3, -rank_base), rank_median, colour = group)) + geom_pointrange(aes(ymin = rank_p05, ymax = rank_p95), size = .3) +
  coord_flip() + scale_y_reverse() + scale_colour_manual(values = col_g) + th +
  labs(x = NULL, y = "Rank (median and 5-95% under random weights)", colour = NULL, title = "Robustness of SD index ranks")
save_fig(F4, "F4_rank_intervals", 7, 9)

gl <- T2 %>% select(iso3, starts_with("gap_vs_EUmean_")) %>% pivot_longer(-iso3, names_to = "component", values_to = "gap", names_prefix = "gap_vs_EUmean_") %>%
  filter(component %in% c("econ", "social", ENV_SUBS))
F5 <- ggplot(gl, aes(gap, reorder(iso3, gap), fill = gap > 0)) + geom_col(show.legend = FALSE) + facet_wrap(~ component, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "#54A24B", `FALSE` = "#E45756")) + th +
  labs(x = "Gap vs EU population-weighted mean (index points)", y = NULL, title = "Profile of candidate-country gaps: economy, social, green sub-blocks")
save_fig(F5, "F5_gap_profile_candidates", 12, 5)

F6 <- ggplot(filter(DEC, window == "main 2012-14 -> 2021-23"), aes(g_gdp, g_co2, colour = group)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") + geom_point(size = 2) + geom_text_repel(aes(label = iso3), size = 3, max.overlaps = 40, show.legend = FALSE) +
  scale_colour_manual(values = col_g) + th + labs(x = "GDP growth (PPP), 2012-14 -> 2021-23", y = "CO2 change", colour = NULL,
                                                  title = "Decoupling of CO2 from GDP", subtitle = "Below the dashed line: relative decoupling; below zero: absolute decoupling")
save_fig(F6, "F6_decoupling_scatter")

cv <- idx %>% filter(year %in% c(P0, P1)) %>% group_by(iso3, group) %>%
  summarise(across(c(sd_arith, co2int_log), list(x0 = ~ mean(.x[year %in% P0], na.rm = TRUE), x1 = ~ mean(.x[year %in% P1], na.rm = TRUE))), .groups = "drop")
cvl <- bind_rows(transmute(cv, iso3, group, var = "SD index (0-100)", x0 = sd_arith_x0, rate = (sd_arith_x1 - sd_arith_x0) / Tspan),
                 transmute(cv, iso3, group, var = "Carbon intensity (log kg CO2/$)", x0 = co2int_log_x0, rate = (co2int_log_x1 - co2int_log_x0) / Tspan))
F7 <- ggplot(cvl, aes(x0, rate, colour = group)) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey40", formula = y ~ x) +
  geom_point(size = 2) + geom_text_repel(aes(label = iso3), size = 2.8, max.overlaps = 30, show.legend = FALSE) +
  facet_wrap(~ var, scales = "free") + scale_colour_manual(values = col_g) + th +
  labs(x = "Initial level (mean 2012-14)", y = "Annual change 2012-14 -> 2021-23", colour = NULL, title = "Beta convergence")
save_fig(F7, "F7_beta_convergence", 10, 5)

F8 <- ggplot(filter(SIG, var %in% c("sd_arith", "econ", "social", "env")), aes(year, gap_cand_minus_EU)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(colour = "#E45756", linewidth = .9) + facet_wrap(~ var, scales = "free_y") + th +
  labs(x = NULL, y = "Candidates minus EU mean (index points)", title = "Catching up: gap between candidate and EU means")
save_fig(F8, "F8_gap_over_time", 9, 6)
png(file.path(OUT_DIR, "figures", "F9_dendrogram.png"), 1800, 1000, res = 150)
plot(hc, main = paste0("Ward clustering, EU-27 + candidates (k = ", k_best, ")"), xlab = "", sub = ""); rect.hclust(hc, k = k_best, border = 2:6); invisible(dev.off())


# ---------------------------------------------------------------------
# 10-a. FIGURES (BLACK & WHITE / MONOCHROME EDITION)
# ---------------------------------------------------------------------

# Setting up a strict black and white theme
th <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey80"), # Clear panel borders for facets
    legend.key = element_blank()
  )

# Universal monochrome palettes for groups
col_bw   <- c(EU = "grey30", candidate = "grey70")
shape_bw <- c(EU = 16, candidate = 17) # Circle for EU, Triangle for candidate countries

# F1: Institutions and sustainable development
F1 <- ggplot(xs, aes(inst_z, sd_nogdp, shape = group, colour = group)) +
  geom_smooth(data = eu, method = "lm", se = TRUE, colour = "black", fill = "grey90", formula = y ~ x, linetype = "solid") +
  geom_point(size = 2.5) + 
  geom_text_repel(aes(label = iso3), size = 3, max.overlaps = 40, show.legend = FALSE, colour = "black") +
  scale_colour_manual(values = col_bw) + 
  scale_shape_manual(values = shape_bw) + 
  th +
  labs(x = "Institutional index (z, PC1)", y = "Social + green SD score (0-100)", colour = NULL, shape = NULL,
       title = "Institutions and sustainable development", subtitle = "Black line: model fitted on EU-27 only")
save_fig(F1, "F1_institutions_vs_SD")


# F2: Institutional trajectories of candidate countries
# Using distinct linetypes for country arrows instead of colors to keep it readable in B&W
pos <- mon %>% select(iso3, name, group, year, cap_z, dem_z) %>% filter(year %in% c(2023, 2025)) %>% pivot_wider(names_from = year, values_from = c(cap_z, dem_z), names_sep = "_")
F2 <- ggplot() +
  geom_point(data = filter(pos, group == "EU"), aes(dem_z_2023, cap_z_2023), colour = "grey75", shape = 16, size = 1.8) +
  geom_segment(data = filter(pos, group == "candidate"), 
               aes(x = dem_z_2023, y = cap_z_2023, xend = dem_z_2025, yend = cap_z_2025, linetype = iso3),
               arrow = arrow(length = unit(.18, "cm")), linewidth = .7, colour = "black", show.legend = TRUE) +
  geom_text_repel(data = filter(pos, group == "candidate"), aes(dem_z_2025, cap_z_2025, label = iso3), 
                  size = 3.2, show.legend = FALSE, colour = "black") +
  th + 
  theme(legend.position = "right") +
  labs(x = "Democratic accountability factor (z)", y = "State capacity factor (z)", linetype = "Country",
       title = "Institutional position of candidate countries, 2023 -> 2025", 
       subtitle = "Grey: EU-27 in 2023; Arrows & Types: Candidates from 2023 to 2025")
save_fig(F2, "F2_capacity_vs_accountability_2023_2025")


# F3: Institutional trajectories of Georgia and Moldova
# Differentiating data sources using varying linetypes and point shapes
ev <- events %>% left_join(distinct(panel, iso3, name), by = "iso3")
cap_txt <- paste0("Dashed lines: ", paste(paste0(ev$iso3, " ", ev$year, ": ", ev$label), collapse = "; "))
F3 <- ggplot(filter(tri_long, iso3 %in% c("GEO", "MDA")), aes(year, value01, linetype = source, shape = source)) +
  geom_vline(data = ev, aes(xintercept = year), linetype = "dotted", colour = "grey40", linewidth = 0.5) +
  geom_line(linewidth = .8, na.rm = TRUE, colour = "black") + 
  geom_point(size = 1.5, na.rm = TRUE, colour = "black") +
  facet_grid(panel_type ~ name) + 
  th + 
  theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.box = "horizontal") +
  guides(linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "Score rescaled to 0-100 (2012-2023 pooled min-max)", linetype = NULL, shape = NULL,
       title = "Institutional trajectories of Georgia and Moldova, 2009-2025", 
       subtitle = "Triangulation across independent sources", caption = cap_txt)
save_fig(F3, "F3_trajectories_GEO_MDA_triangulation", 10, 7)


# F4: Robustness of SD index ranks
F4 <- ggplot(R2, aes(reorder(iso3, -rank_base), rank_median, shape = group)) + 
  geom_pointrange(aes(ymin = rank_p05, ymax = rank_p95), size = .4, colour = "black") +
  coord_flip() + 
  scale_y_reverse() + 
  scale_shape_manual(values = shape_bw) + 
  th +
  labs(x = NULL, y = "Rank (median and 5-95% under random weights)", shape = NULL, title = "Robustness of SD index ranks")
save_fig(F4, "F4_rank_intervals", 7, 9)


# F5: Profile of candidate-country gaps
# Replacing color with contrast fills: dark grey for negative gaps, light grey for positive gaps
gl <- T2 %>% select(iso3, starts_with("gap_vs_EUmean_")) %>% pivot_longer(-iso3, names_to = "component", values_to = "gap", names_prefix = "gap_vs_EUmean_") %>%
  filter(component %in% c("econ", "social", ENV_SUBS))
F5 <- ggplot(gl, aes(gap, reorder(iso3, gap), fill = gap > 0)) + 
  geom_col(show.legend = FALSE, colour = "black", linewidth = 0.3) + 
  facet_wrap(~ component, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "grey85", `FALSE` = "grey40")) + 
  th +
  labs(x = "Gap vs EU population-weighted mean (index points)", y = NULL, 
       title = "Profile of candidate-country gaps: economy, social, green sub-blocks",
       subtitle = "Dark grey: negative gap (below EU) | Light grey: positive gap (above EU)")
save_fig(F5, "F5_gap_profile_candidates", 12, 5)


# F6: Decoupling of CO2 from GDP
F6 <- ggplot(filter(DEC, window == "main 2012-14 -> 2021-23"), aes(g_gdp, g_co2, shape = group, colour = group)) + 
  geom_hline(yintercept = 0, colour = "grey50", linetype = "solid") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotdash", colour = "grey40") + 
  geom_point(size = 2.5) + 
  geom_text_repel(aes(label = iso3), size = 3, max.overlaps = 40, show.legend = FALSE, colour = "black") +
  scale_colour_manual(values = col_bw) + 
  scale_shape_manual(values = shape_bw) + 
  th + 
  labs(x = "GDP growth (PPP), 2012-14 -> 2021-23", y = "CO2 change", colour = NULL, shape = NULL,
       title = "Decoupling of CO2 from GDP", 
       subtitle = "Below dot-dash line: relative decoupling; below zero line: absolute decoupling")
save_fig(F6, "F6_decoupling_scatter")


# F7: Beta convergence (Cleaned version without warnings)
cv <- idx %>% filter(year %in% c(P0, P1)) %>% group_by(iso3, group) %>%
  summarise(across(c(sd_arith, co2int_log), list(x0 = ~ mean(.x[year %in% P0], na.rm = TRUE), x1 = ~ mean(.x[year %in% P1], na.rm = TRUE))), .groups = "drop")
cvl <- bind_rows(transmute(cv, iso3, group, var = "SD index (0-100)", x0 = sd_arith_x0, rate = (sd_arith_x1 - sd_arith_x0) / Tspan),
                 transmute(cv, iso3, group, var = "Carbon intensity (log kg CO2/$)", x0 = co2int_log_x0, rate = (co2int_log_x1 - co2int_log_x0) / Tspan))
F7 <- ggplot(cvl, aes(x0, rate)) + 
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", linetype = "dashed", formula = y ~ x) +
  geom_point(aes(shape = group, colour = group), size = 2.5) + 
  geom_text_repel(aes(label = iso3), size = 2.8, max.overlaps = 30, show.legend = FALSE, colour = "black") +
  facet_wrap(~ var, scales = "free") + 
  scale_colour_manual(values = col_bw) + 
  scale_shape_manual(values = shape_bw) + 
  th +
  labs(x = "Initial level (mean 2012-14)", y = "Annual change 2012-14 -> 2021-23", colour = NULL, shape = NULL, title = "Beta convergence")

save_fig(F7, "F7_beta_convergence", 10, 5)


# F8: Catching up (Gap over time)
# Increased line thickness and pure black color for enhanced contrast against background grids
F8 <- ggplot(filter(SIG, var %in% c("sd_arith", "econ", "social", "env")), aes(year, gap_cand_minus_EU)) + 
  geom_hline(yintercept = 0, colour = "grey50", linetype = "solid") +
  geom_line(colour = "black", linewidth = 1.1) + 
  facet_wrap(~ var, scales = "free_y") + 
  th +
  labs(x = NULL, y = "Candidates minus EU mean (index points)", title = "Catching up: gap between candidate and EU means")
save_fig(F8, "F8_gap_over_time", 9, 6)


# F9: Dendrogram
# Setting monochrome cluster boundaries using distinct linetypes (lty 1 to k_best)
 
png(file.path(OUT_DIR, "figures", "F9_dendrogram.png"), 1800, 1000, res = 150)
plot(hc, main = paste0("Ward clustering, EU-27 + candidates (k = ", k_best, ")"), xlab = "", sub = "", col = "black")
lty_sequence <- 1:k_best
for (i in 1:k_best) {
  rect.hclust(hc, k = k_best, which = i, border = "grey30")
}
invisible(dev.off())



