# =============================================================================
# SCRIPT 2b — SENSITIVITY ANALYSIS AND TABLE 4.1
# =============================================================================
#  Purpose
#    Run baseline, group-specific, and sensitivity DST configurations, then
#    compile the full Table 4.1 summary statistics. Equal-weight rows (f, g, h)
#    are loaded from Script 2a outputs rather than re-run.
#
#  Sections
#    0  Environment         libraries, paths, directories, log connection
#    1  Staleness check     hash-verify inputs against Script 1 and 2a metadata
#    2  Load data           d1, dt.m, Script 2a best/duo/trio outputs
#    3  Helpers             add_gap_columns(); table_row() — one Table 4.1 row
#                           per DST output; columns match table41_stats.csv
#    4  Baseline (row a)    reference DST run (parameters under DST configurations)
#    5  Group runs (b-e)    one management group restricted at a time
#    6  Sensitivity runs    user weights (i-t), timeframe (u-x), cumulative
#                           effects (y-z); results collected into sens_stats
#    7  Compile Table 4.1   bind all rows; log summary
#    8  Run metadata        provenance record saved to run_meta_script2b.rds
#    9  Save outputs        table41_stats.csv; re-save RDS files with metadata
#   10  Figure              Fig_sensitivity.png — three-panel bar+line chart (Figure 4.2)
#                           across all sensitivity configurations
#
#  Inputs
#    results/d1.rds                from Script 1
#    results/dt_m.rds              from Script 1
#    results/run_meta_script1.rds  upstream staleness check
#    results/run_meta_script2a.rds chain-of-evidence check
#    results/out_best_n1.rds       from Script 2a (Table 4.1 row f)
#    results/out_duo_n2.rds        from Script 2a (Table 4.1 row g)
#    results/out_trio_n3.rds       from Script 2a (Table 4.1 row h)
#
#  Outputs
#    results/table41_stats.csv     compiled Table 4.1 (26 rows)
#    results/out_group_*.rds       group-specific best-measure outputs (b-e)
#    results/out_sens_*.rds        sensitivity best-measure outputs (i-z)
#    results/run_meta_script2b.rds run metadata
#    maps/Fig_sensitivity.png      three-panel sensitivity figure (report Figure 4.2)
#    logs/script_2b_report.txt
#
#  DST configurations
#    Baseline  (a)   nmax = 1, uw = c(1,1,1), simyear = 5
#    Groups    (b-e) nmax = 1, uw = c(1,1,1), one group at a time
#      b: Nutrient type      CF-MF, OF-MF
#      c: Nutrient eff.      EE, RFT, RFP
#      d: Soil conservation  NT-CT, RT-CT
#      e: Cropping           CC, ROT, RES
#    User-weight sensitivity (i-t) — nmax = 1/2/3 per row
#      i-k: uw = c(2,1,1)   yield-prioritised
#      l-n: uw = c(1,2,2)   environment-prioritised
#      o-q: uw = c(4,3,3)   multi-stakeholder, yield-leaning
#      r-t: uw = c(3,4,4)   multi-stakeholder, environment-leaning
#    Timeframe sensitivity (u-x) — uw = c(1,1,1), nmax = 1
#      u: simyear = 10       v: simyear = 15
#      w: simyear = 20       x: simyear = 100
#    Cumulative-effects sensitivity (y-z) — uw = c(1,1,1), simyear = 5, nmax = 1
#      y: yield_factor = 3, nsu_factor = 1   cumulative yield
#      z: yield_factor = 3, nsu_factor = 2   cumulative yield + N surplus
#
#  References
#    Young et al. (2021, 2022); dst_functions_fert.R, dst_main_fert.R;
#    phd_maddy GitHub repository.
# =============================================================================


# 0. Environment --------------------------------------------------------------

setwd('C:/Users/STMP2/Documents/Thesis/Github/Run/1')
rm(list = ls())

library(data.table)
library(digest)
library(readxl)

source('scripts/paths_original.R')
source(P$fn_fert)

if (!dir.exists('results')) dir.create('results')
if (!dir.exists('logs'))    dir.create('logs')
if (!dir.exists('maps'))    dir.create('maps')

run_ts <- Sys.time()

require_input(P$d1)
require_input(P$dt_m)
require_input(P$out_best)
require_input(P$out_duo)
require_input(P$out_trio)
require_input(P$meta_s1)
require_input(P$meta_s2a)

report_con <- file(P$log_s2b, open = 'wt')
log <- function(...) {
  msg <- paste0(...)
  cat(msg, '\n')
  tryCatch(
    cat(msg, '\n', file = report_con, append = TRUE),
    error = function(e) invisible(NULL)
  )
}
log('=== Script 2b run: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')


# 1. Staleness verification ---------------------------------------------------

script1_meta  <- readRDS(P$meta_s1)
script2a_meta <- readRDS(P$meta_s2a)
log('Script 1 finished  : ', format(script1_meta$finished_ts,  '%Y-%m-%d %H:%M:%S'))
log('Script 2a finished : ', format(script2a_meta$finished_ts, '%Y-%m-%d %H:%M:%S'))

db_hash_now <- digest::digest(file = P$db, algo = 'md5')
if (db_hash_now != script1_meta$input_hashes[['db']])
  stop('Database has changed since Script 1 ran. Re-run Script 1.')

d1_hash_now <- digest::digest(file = P$d1, algo = 'md5')
if (!is.null(script2a_meta$input_hashes[['d1']]) &&
    d1_hash_now != script2a_meta$input_hashes[['d1']])
  stop('d1 has changed since Script 2a ran. Re-run Script 2a.')


# 2. Load data ----------------------------------------------------------------

d1   <- readRDS(P$d1)
dt.m <- readRDS(P$dt_m)
dt.m <- dt.m[man_code != 'RFR']

ee_ha <- dt.m[man_code == 'EE' & indicator == 'Y', mean(ha_m1, na.rm = TRUE)]
if (ee_ha <= 0)
  stop('EE applicable area is zero or negative. Re-run Script 1.')

log('Measures for optimisation: ',
    paste(sort(unique(dt.m$man_code)), collapse = ', '))

out.best <- readRDS(P$out_best)
out.duo  <- readRDS(P$out_duo)
out.trio <- readRDS(P$out_trio)


# 3. Helpers ------------------------------------------------------------------

# Reference gap columns (dst_main_fert.R).
add_gap_columns <- function(out) {
  out[, soc_ref_s    := bd * soc_ref * 30]
  out[, D_Y          := ((1 + dY)   * yield_ref) - yield_ref]
  out[, D_SOC_s      := ((1 + dSOC) * soc_ref_s) - soc_ref_s]
  out[, D_Nsu        := ((1 + dNsu) * n_sp_ref)  - n_sp_ref]
  out[, yield_new    := yield_ref + D_Y]
  out[, yield_new_t  := yield_new / 1000]
  out[, soc_new      := soc_ref_s + D_SOC_s]
  out[, soc_new_t    := soc_new   / 1000]
  out[, n_sp_new     := n_sp_ref  + D_Nsu]
  out[, soc_targ_t   := (bd * soc_target * 30) / 1000]
  out[, soc_ref_t    := soc_ref_s    / 1000]
  out[, yield_targ_t := yield_target / 1000]
  out[, yield_ref_t  := yield_ref    / 1000]
  out[, yield_gap_t     := yield_targ_t - yield_ref_t]
  out[, soc_gap_t       := soc_targ_t   - soc_ref_t]
  out[, n_sp_gap        := n_sp_ref     - n_sp_crit]
  out[, yield_gap_fin_t := yield_targ_t - yield_new_t]
  out[, soc_gap_fin_t   := soc_targ_t   - soc_new_t]
  out[, n_sp_gap_fin    := n_sp_new     - n_sp_crit]
  out[, yield_gap_diff_p := (yield_gap_fin_t - yield_gap_t) / abs(yield_gap_t) * 100]
  out[, soc_gap_diff_p   := (soc_gap_fin_t   - soc_gap_t)   / abs(soc_gap_t)   * 100]
  out[, n_sp_gap_diff_p  := (n_sp_gap_fin    - n_sp_gap)    / abs(n_sp_gap)    * 100]
  # Inf/NaN -> 0: NCUs already at target have a zero denominator
  for (col in c('yield_gap_diff_p', 'soc_gap_diff_p', 'n_sp_gap_diff_p'))
    out[is.infinite(get(col)) | is.nan(get(col)) | is.na(get(col)), (col) := 0]
  out
}

# Compute a single Table 4.1 row from a DST output object.
# Yield and SOC means are geometric (log-transform then back-transform)
# per the reference Table 7 footnote of Young et al. (2022).
# Column names match table41_stats.csv: pct_area_Y / pct_area_C / pct_area_N.
# pct_area_* = area-weighted share; pct_ncu_* = unweighted NCU-count share.
# Area share of a per-NCU 0/1 target-met flag, weighted by NCU cropland area.
# The report expresses indicator attainment as a share of agricultural AREA;
# a plain sum(flag)/n_total would instead give the unweighted NCU-count share.
aw_pct <- function(flag, area)
  round(sum(flag * area, na.rm = TRUE) / sum(area, na.rm = TRUE) * 100, 1)

# Geometric mean in Mg ha-1 (log-transform, back-transform), guarding against
# non-positive / non-finite values that would otherwise collapse it to NaN.
geomean_Mg <- function(x) {
  v <- x / 1000
  v <- v[is.finite(v) & v > 0]
  if (length(v) == 0) return(NA_real_)
  exp(mean(log(v)))
}

table_row <- function(out, label) {
  if (!'soc_ref_s' %in% colnames(out)) {
    out <- copy(out)
    out[, soc_ref_s   := bd * soc_ref * 30]
    out[, n_sp_new    := (1 + dNsu) * n_sp_ref]
    out[, yield_new_t := (1 + dY)   * yield_ref / 1000]
    out[, soc_new_t   := (1 + dSOC) * soc_ref_s / 1000]
  }
  n_total   <- uniqueN(out$ncu)
  has_bipmc <- 'bipmc' %in% names(out) && var(out[['bipmc']], na.rm = TRUE) > 0
  data.table(
    scenario   = label,
    Y_Mg       = round(geomean_Mg(out$yield_ref), 2),
    Y_Mg_sd    = round(sd(out$yield_ref / 1000,              na.rm = TRUE), 2),
    dY_pct     = round(mean(out$dY,   na.rm = TRUE) * 100, 1),
    dY_pct_sd  = round(sd(out$dY,     na.rm = TRUE) * 100, 1),
    C_Mg       = round(geomean_Mg(out$soc_ref_s), 1),
    C_Mg_sd    = round(sd(out$soc_ref_s / 1000,              na.rm = TRUE), 1),
    dC_pct     = round(mean(out$dSOC,  na.rm = TRUE) * 100, 2),
    dC_pct_sd  = round(sd(out$dSOC,    na.rm = TRUE) * 100, 2),
    N_kg       = round(mean(out$n_sp_ref, na.rm = TRUE), 1),
    N_kg_sd    = round(sd(out$n_sp_ref,   na.rm = TRUE), 1),
    dN_pct     = round(mean(out$dNsu,  na.rm = TRUE) * 100, 1),
    dN_pct_sd  = round(sd(out$dNsu,    na.rm = TRUE) * 100, 1),
    MCA_index  = if (has_bipmc) round(mean(out[['bipmc']], na.rm = TRUE), 3)
    else NA_real_,
    MCA_sd     = if (has_bipmc) round(sd(out[['bipmc']],   na.rm = TRUE), 3)
    else NA_real_,
    pct_area_Y = aw_pct(out$tm_Y, out$area_ncu_ha_tot),
    pct_area_C = aw_pct(out$tm_C, out$area_ncu_ha_tot),
    pct_area_N = aw_pct(out$tm_N, out$area_ncu_ha_tot),
    pct_ncu_Y  = round(sum(out$tm_Y, na.rm = TRUE) / n_total * 100, 1),
    pct_ncu_C  = round(sum(out$tm_C, na.rm = TRUE) / n_total * 100, 1),
    pct_ncu_N  = round(sum(out$tm_N, na.rm = TRUE) / n_total * 100, 1)
  )
}


# 4. Baseline (Table 4.1 row a) ------------------------------------------------

log('Running baseline (row a)...')
sim_base <- runDST(db = d1, dt.m = dt.m, output = 'best_impact',
                   uw = c(1, 1, 1), simyear = 5, quiet = TRUE,
                   nmax = 1, nopt = FALSE)
base_out <- add_gap_columns(sim_base$impact_best)

n_total    <- uniqueN(base_out$ncu)
pct_Y_base <- sum(base_out$ti_Y, na.rm = TRUE) / n_total * 100
pct_C_base <- sum(base_out$ti_C, na.rm = TRUE) / n_total * 100
pct_N_base <- sum(base_out$ti_N, na.rm = TRUE) / n_total * 100
pct_Y_base_aw <- aw_pct(base_out$ti_Y, base_out$area_ncu_ha_tot)
pct_C_base_aw <- aw_pct(base_out$ti_C, base_out$area_ncu_ha_tot)
pct_N_base_aw <- aw_pct(base_out$ti_N, base_out$area_ncu_ha_tot)
log(sprintf('  Baseline targets met (area share): Y=%.1f%%  SOC=%.1f%%  N=%.1f%%',
            pct_Y_base_aw, pct_C_base_aw, pct_N_base_aw))

baseline_stats <- data.table(
  scenario   = 'Baseline (a)',
  Y_Mg       = round(geomean_Mg(base_out$yield_ref), 2),
  Y_Mg_sd    = round(sd(base_out$yield_ref / 1000,            na.rm = TRUE), 2),
  dY_pct     = 0, dY_pct_sd = 0,
  C_Mg       = round(geomean_Mg(base_out$soc_ref_s), 1),
  C_Mg_sd    = round(sd(base_out$soc_ref_s / 1000,            na.rm = TRUE), 1),
  dC_pct     = 0, dC_pct_sd = 0,
  N_kg       = round(mean(base_out$n_sp_ref, na.rm = TRUE), 1),
  N_kg_sd    = round(sd(base_out$n_sp_ref,   na.rm = TRUE), 1),
  dN_pct     = 0, dN_pct_sd = 0,
  MCA_index  = 0, MCA_sd    = 0,
  pct_area_Y = pct_Y_base_aw,
  pct_area_C = pct_C_base_aw,
  pct_area_N = pct_N_base_aw,
  pct_ncu_Y  = round(pct_Y_base, 1),
  pct_ncu_C  = round(pct_C_base, 1),
  pct_ncu_N  = round(pct_N_base, 1)
)
rm(sim_base)


# 5. Group runs (rows b-e) ---------------------------------------------------

groups <- list(
  list(label = 'Nutrient_type_b',       measures = c('CF-MF', 'OF-MF')),
  list(label = 'Nutrient_efficiency_c', measures = c('EE', 'RFT', 'RFP')),
  list(label = 'Soil_conservation_d',   measures = c('NT-CT', 'RT-CT')),
  list(label = 'Cropping_e',            measures = c('CC', 'ROT', 'RES'))
)

group_results <- list()
group_stats   <- list()
for (g in groups) {
  log('Running group: ', g$label)
  dt.m_grp <- dt.m[man_code %in% g$measures]
  sim_grp  <- runDST(db = d1, dt.m = dt.m_grp, output = 'best_impact',
                     uw = c(1, 1, 1), simyear = 5, quiet = TRUE,
                     nmax = 1, nopt = FALSE)
  out_grp <- add_gap_columns(sim_grp$impact_best)
  saveRDS(out_grp, file.path('results', paste0('out_group_', g$label, '.rds')))
  group_results[[g$label]] <- out_grp
  group_stats[[g$label]]   <- table_row(out_grp, g$label)
  log('  Measure frequency: ',
      paste(names(sort(table(out_grp$man_code), decreasing = TRUE)),
            collapse = ', '))
}


# 6. Sensitivity runs (rows i-z) ---------------------------------------------

# User-weight sensitivity (rows i-t) — nmax = 1/2/3 per row.
# Under unequal priorities the composite index uses within-site score
# normalisation; values are not directly comparable with the baseline.
uw_configs <- list(
  list(label = 'uw211_i', uw = c(2,1,1), nmax = 1, output = 'best_impact'),
  list(label = 'uw211_j', uw = c(2,1,1), nmax = 2, output = 'score_duo'),
  list(label = 'uw211_k', uw = c(2,1,1), nmax = 3, output = 'score_trio'),
  list(label = 'uw122_l', uw = c(1,2,2), nmax = 1, output = 'best_impact'),
  list(label = 'uw122_m', uw = c(1,2,2), nmax = 2, output = 'score_duo'),
  list(label = 'uw122_n', uw = c(1,2,2), nmax = 3, output = 'score_trio'),
  list(label = 'uw433_o', uw = c(4,3,3), nmax = 1, output = 'best_impact'),
  list(label = 'uw433_p', uw = c(4,3,3), nmax = 2, output = 'score_duo'),
  list(label = 'uw433_q', uw = c(4,3,3), nmax = 3, output = 'score_trio'),
  list(label = 'uw344_r', uw = c(3,4,4), nmax = 1, output = 'best_impact'),
  list(label = 'uw344_s', uw = c(3,4,4), nmax = 2, output = 'score_duo'),
  list(label = 'uw344_t', uw = c(3,4,4), nmax = 3, output = 'score_trio')
)

# Timeframe sensitivity (rows u-x) — uw = c(1,1,1), nmax = 1.
# Only SOC accumulates over simyear; yield and N-surplus are single-season.
simyear_configs <- list(
  list(label = 'sim_u_10y',  simyear = 10),
  list(label = 'sim_v_15y',  simyear = 15),
  list(label = 'sim_w_20y',  simyear = 20),
  list(label = 'sim_x_100y', simyear = 100)
)

# Cumulative-effects sensitivity (rows y-z) — uw = c(1,1,1), simyear = 5.
# Replicates Young et al. (2022) Table 11. Requires the modified runDST with
# yield_factor and nsu_factor arguments (default 1 = reference behaviour).
cumeff_configs <- list(
  list(label = 'cum_y_y',  yield_factor = 3, nsu_factor = 1),
  list(label = 'cum_z_yn', yield_factor = 3, nsu_factor = 2)
)

sens_results <- list()
sens_stats   <- list()

for (s in uw_configs) {
  log('Running sensitivity: ', s$label)
  sim_s <- runDST(db = d1, dt.m = dt.m, output = s$output,
                  uw = s$uw, simyear = 5,
                  quiet = TRUE, nmax = s$nmax, nopt = FALSE)
  out_s <- switch(s$output,
                  'best_impact' = sim_s$impact_best,
                  'score_duo'   = sim_s$score_duo,
                  'score_trio'  = sim_s$score_trio)
  if (s$output == 'best_impact') out_s <- add_gap_columns(out_s)
  saveRDS(out_s, file.path('results', paste0('out_sens_', s$label, '.rds')))
  sens_results[[s$label]] <- out_s
  sens_stats[[s$label]]   <- table_row(out_s, s$label)
}

for (s in simyear_configs) {
  log('Running sensitivity: ', s$label, '  (simyear = ', s$simyear, ')')
  sim_s <- runDST(db = d1, dt.m = dt.m, output = 'best_impact',
                  uw = c(1,1,1), simyear = s$simyear,
                  quiet = TRUE, nmax = 1, nopt = FALSE)
  out_s <- add_gap_columns(sim_s$impact_best)
  saveRDS(out_s, file.path('results', paste0('out_sens_', s$label, '.rds')))
  sens_results[[s$label]] <- out_s
  sens_stats[[s$label]]   <- table_row(out_s, s$label)
}

for (s in cumeff_configs) {
  log('Running sensitivity: ', s$label,
      '  (yield_factor = ', s$yield_factor, ', nsu_factor = ', s$nsu_factor, ')')
  sim_s <- runDST(db = d1, dt.m = dt.m, output = 'best_impact',
                  uw = c(1,1,1), simyear = 5,
                  quiet = TRUE, nmax = 1, nopt = FALSE,
                  yield_factor = s$yield_factor,
                  nsu_factor   = s$nsu_factor)
  out_s <- add_gap_columns(sim_s$impact_best)
  saveRDS(out_s, file.path('results', paste0('out_sens_', s$label, '.rds')))
  sens_results[[s$label]] <- out_s
  sens_stats[[s$label]]   <- table_row(out_s, s$label)
}


# 7. Compile Table 4.1 ----------------------------------------------------------

main_stats <- rbind(
  table_row(out.best, 'One measure (f)'),
  table_row(out.duo,  'One or two combined (g)'),
  table_row(out.trio, 'One two or three combined (h)')
)

table41 <- rbind(
  baseline_stats,
  rbindlist(group_stats, fill = TRUE),
  main_stats,
  rbindlist(sens_stats,  fill = TRUE),
  fill = TRUE
)

log('Table 4.1: ', nrow(table41), ' rows compiled.')
log(paste(capture.output(
  print(table41[, .(scenario, dY_pct, dN_pct, MCA_index,
                   pct_area_Y, pct_area_C, pct_area_N)])
), collapse = '\n'))


# 8. Run metadata -------------------------------------------------------------

run_meta <- list(
  script            = 'Script_2b_Sensitivity_Analyses.R',
  started_ts        = run_ts,
  finished_ts       = Sys.time(),
  r_version         = R.version.string,
  session_info      = sessionInfo(),
  input_hashes      = list(
    d1       = digest::digest(file = P$d1,       algo = 'md5'),
    dt_m     = digest::digest(file = P$dt_m,     algo = 'md5'),
    out_best = digest::digest(file = P$out_best, algo = 'md5'),
    out_duo  = digest::digest(file = P$out_duo,  algo = 'md5'),
    out_trio = digest::digest(file = P$out_trio, algo = 'md5')
  ),
  upstream_meta     = list(
    script1_finished  = script1_meta$finished_ts,
    script2a_finished = script2a_meta$finished_ts
  ),
  baseline_uw       = c(1, 1, 1),
  baseline_simyear  = 5,
  group_runs        = sapply(groups,          `[[`, 'label'),
  uw_runs           = sapply(uw_configs,      `[[`, 'label'),
  simyear_runs      = sapply(simyear_configs, `[[`, 'label'),
  cumeff_runs       = sapply(cumeff_configs,  `[[`, 'label'),
  n_table_rows      = nrow(table41)
)

attr(table41, 'run_meta') <- run_meta

# 9. Save outputs -------------------------------------------------------------

fwrite(table41, 'results/table41_stats.csv')
log('Saved: results/table41_stats.csv')

# Re-save group and sensitivity RDS files with run_meta attribute for provenance
for (g_label in names(group_results)) {
  out_g <- group_results[[g_label]]
  attr(out_g, 'run_meta') <- run_meta
  saveRDS(out_g, file.path('results', paste0('out_group_', g_label, '.rds')))
}
for (s_label in names(sens_results)) {
  out_s <- sens_results[[s_label]]
  attr(out_s, 'run_meta') <- run_meta
  saveRDS(out_s, file.path('results', paste0('out_sens_', s_label, '.rds')))
}
log('Re-saved group and sensitivity RDS files with run_meta attribute.')


# 10. Figure 4.2 — sensitivity of best-measure outcomes to model parameters -
# Produces maps/Fig_sensitivity.png = report Figure 4.2.
# Three-panel layout: (a) user priority settings, (b) simulation timeframe,
# (c) cumulative effects. Conventions follow Script 3 / dst_plotfigures.R:
# ggplot2, theme_bw(), ggsave 36 x 12 cm / 400 dpi.

if (requireNamespace('ggplot2',   quietly = TRUE) &&
    requireNamespace('patchwork', quietly = TRUE)) {

  library(ggplot2)
  library(patchwork)
  theme_set(theme_bw())

  tx <- fread('results/table41_stats.csv')

  C_Y   <- '#4393c3'   # blue       — yield target bar
  C_C   <- '#92c5de'   # mid-blue   — SOC target bar
  C_N   <- '#d1e5f0'   # pale-blue  — N surplus target bar
  C_MCA <- '#d73027'   # red        — composite index line
  C_DY  <- '#2166ac'   # dark blue  — ΔYield line  (panel c)
  C_DN  <- '#4dac26'   # green      — ΔN line    (panel c)

  # ── panels (a) and (b) ──────
  make_ab_panel <- function(scenarios, x_labels, panel_title, footnote) {
    ds <- tx[scenario %in% scenarios]
    ds[, x_pos := match(scenario, scenarios) - 1L]
    ww <- 0.22

    ggplot(ds, aes(x = x_pos)) +
      geom_col(aes(y = pct_area_Y, fill = 'Yield target met (%)'),
               width = ww, position = position_nudge(x = -ww)) +
      geom_col(aes(y = pct_area_C, fill = 'SOC target met (%)'),
               width = ww, position = position_nudge(x = 0)) +
      geom_col(aes(y = pct_area_N, fill = 'N surplus target met (%)'),
               width = ww, position = position_nudge(x =  ww)) +
      geom_line(aes(y = MCA_index * 100, colour = 'MCA index', group = 1),
                linewidth = 1.4, linetype = 'dashed') +
      geom_point(aes(y = MCA_index * 100, colour = 'MCA index'),
                 size = 3, shape = 16) +
      scale_fill_manual(
        name   = NULL,
        values = c('Yield target met (%)'    = C_Y,
                   'SOC target met (%)'      = C_C,
                   'N surplus target met (%)' = C_N),
        breaks = c('Yield target met (%)','SOC target met (%)','N surplus target met (%)')) +
      scale_colour_manual(
        name   = NULL,
        values = c('MCA index' = C_MCA)) +
      scale_x_continuous(breaks = seq_along(x_labels) - 1, labels = x_labels) +
      scale_y_continuous(
        name     = 'Area with target met (%)',
        limits   = c(0, 100), expand = c(0, 0),
        breaks   = seq(0, 100, by = 20),
        sec.axis = sec_axis(~ . / 100, name = 'MCA index',
                            breaks = seq(0, 1, by = 0.2))) +
      xlab(NULL) +
      ggtitle(panel_title) +
      labs(caption = footnote) +
      guides(fill   = guide_legend(order = 1,
                                   override.aes = list(linetype = 0, shape = NA)),
             colour = guide_legend(order = 2)) +
      theme(axis.text.x        = element_text(size = 8),
            axis.title.y.right = element_text(colour = C_MCA, size = 8),
            axis.text.y.right  = element_text(colour = C_MCA),
            plot.title         = element_text(size = 10, hjust = 0.5),
            plot.caption       = element_text(size = 7, hjust = 0,
                                              colour = 'grey40', face = 'italic'),
            legend.position    = c(0.01, 0.99),
            legend.justification = c(0, 1),
            legend.background  = element_rect(fill = 'white', colour = NA),
            legend.key.size    = unit(0.4, 'cm'),
            legend.text        = element_text(size = 7),
            panel.border       = element_rect(colour = 'grey60'))
  }

  # ── panel (c) ─────────────────────
  make_c_panel <- function(scenarios, x_labels, panel_title, footnote) {
    ds <- tx[scenario %in% scenarios]
    ds[, x_pos  := match(scenario, scenarios) - 1L]
    ds[, abs_dN := abs(dN_pct)]
    ww <- 0.30

    ggplot(ds, aes(x = x_pos)) +
      geom_col(aes(y = pct_area_Y, fill = 'Yield target met (%)'),
               width = ww, position = position_nudge(x = -ww / 2)) +
      geom_col(aes(y = pct_area_N, fill = 'N surplus target met (%)'),
               width = ww, position = position_nudge(x =  ww / 2)) +
      geom_line(aes(y = dY_pct,          colour = 'ΔYield (%)',    group = 1),
                linewidth = 1.4, linetype = 'solid') +
      geom_point(aes(y = dY_pct,         colour = 'ΔYield (%)'),
                 size = 3, shape = 15) +
      geom_line(aes(y = abs_dN,          colour = 'ΔN surplus (%)', group = 1),
                linewidth = 1.4, linetype = 'solid') +
      geom_point(aes(y = abs_dN,         colour = 'ΔN surplus (%)'),
                 size = 3, shape = 17) +
      geom_line(aes(y = MCA_index * 100, colour = 'MCA index',          group = 1),
                linewidth = 1.4, linetype = 'dashed') +
      geom_point(aes(y = MCA_index * 100, colour = 'MCA index'),
                 size = 3, shape = 16) +
      scale_fill_manual(
        name   = NULL,
        values = c('Yield target met (%)'     = C_Y,
                   'N surplus target met (%)' = C_N),
        breaks = c('Yield target met (%)','N surplus target met (%)')) +
      scale_colour_manual(
        name   = NULL,
        values = c('ΔYield (%)'        = C_DY,
                   'ΔN surplus (%)'  = C_DN,
                   'MCA index'              = C_MCA),
        breaks = c('ΔYield (%)','ΔN surplus (%)','MCA index')) +
      scale_x_continuous(breaks = seq_along(x_labels) - 1, labels = x_labels) +
      scale_y_continuous(
        name     = 'Area with target met (%)',
        limits   = c(0, 100), expand = c(0, 0),
        breaks   = seq(0, 100, by = 20),
        sec.axis = sec_axis(~ ., name = 'Change (%) and MCA Index (×100)',
                            breaks = seq(0, 100, by = 20))) +
      xlab(NULL) +
      ggtitle(panel_title) +
      labs(caption = footnote) +
      guides(fill   = guide_legend(order = 1,
                                   override.aes = list(linetype = 0, shape = NA)),
             colour = guide_legend(order = 2)) +
      theme(axis.text.x          = element_text(size = 8),
            axis.title.y.right   = element_text(size = 8),
            plot.title           = element_text(size = 10, hjust = 0.5),
            plot.caption         = element_text(size = 7, hjust = 0,
                                                colour = 'grey40', face = 'italic'),
            legend.position      = c(0.01, 0.99),
            legend.justification = c(0, 1),
            legend.background    = element_rect(fill = 'white', colour = NA),
            legend.key.size      = unit(0.4, 'cm'),
            legend.text          = element_text(size = 7),
            panel.border         = element_rect(colour = 'grey60'))
  }

  p_wt <- make_ab_panel(
    c('One measure (f)', 'uw211_i', 'uw122_l', 'uw433_o', 'uw344_r'),
    c('1:1:1\n(baseline)', '2:1:1\n(yield)', '1:2:2\n(SOC+N)',
      '4:3:3\n(yield)', '3:4:4\n(SOC+N)'),
    '(a) User priority settings',
    'Index values are not comparable across weighting schemes (\u00a72.2.3).')

  p_tf <- make_ab_panel(
    c('One measure (f)', 'sim_u_10y', 'sim_v_15y', 'sim_w_20y', 'sim_x_100y'),
    c('5y\n(baseline)', '10y', '15y', '20y', '100y'),
    '(b) Simulation timeframe',
    'Only SOC accumulates over simyear; yield and N surplus are\nsingle-season quantities by design.')

  p_cu <- make_c_panel(
    c('One measure (f)', 'cum_y_y', 'cum_z_yn'),
    c('Single-season\n(baseline)', 'Yield \u00d73\n(cumulative yield)',
      'Yield \u00d73 + N\u00d72\n(cumulative yield+N)'),
    '(c) Cumulative effects',
    'yield: yield_factor = 3; yield+N: yield_factor = 3 + nsu_factor = 2.\nIndex values not comparable across cumulative scenarios.')

  fig <- p_wt + p_tf + p_cu +
    plot_layout(ncol = 3) +
    plot_annotation(
      title    = 'Sensitivity of best measure outcomes to model parameters',
      subtitle = 'Bars: share of agricultural area with indicator target met. Line: MCA index.',
      theme    = theme(
        plot.title    = element_text(size = 13, hjust = 0.5, face = 'bold'),
        plot.subtitle = element_text(size =  9, hjust = 0.5)))

  ggsave('maps/Fig_sensitivity.png', fig,
         width = 36, height = 12, units = 'cm', dpi = 400)
  log('10: Fig_sensitivity.png saved to maps/')

} else {
  log('10: ggplot2 and/or patchwork not installed — sensitivity figure skipped.')
  log('    Run install.packages(c("ggplot2","patchwork")) once to enable.')
}

# =============================================================================
saveRDS(run_meta, P$meta_s2b)
log('Saved: ', P$meta_s2b)

log('=== Script 2b complete: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')
log('Duration: ', sprintf('%.1f min',
                          as.numeric(difftime(Sys.time(), run_ts, units = 'mins'))))
close(report_con)
# End of Script 2b ------------------------------------------------------------
