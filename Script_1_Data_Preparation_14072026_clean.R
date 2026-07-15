# =============================================================================
# SCRIPT 1 — DATA PREPARATION
# =============================================================================
#  Purpose
#    Build the d1, ma_models and dt.m objects required by Script 2a/2b for the
#    10-measure DSF optimisation across 29,476 NCUs covering EU-24 + UK.
#    (Croatia (HR), Cyprus (CY), and Malta (MT) absent due to missing NUTS2 codes in 2016 EUROSTAT
#    FSS data.)
#
#  Sections
#    0  Environment                   libraries, paths, log connection
#    1  Load database                 db_final_europe.csv via dst_loaddb.R conventions
#    2  Yield swap                    restore Young et al. (2022) §2.1 column order
#    3  Column checks                 verify required columns are present
#    4  4R partial areas              compute parea.ee / .rft / .rfp from FSS data
#    5  fertilizer cross-check        validate partial-area sums against total
#    6  Meta-analytical models        load and patch MA models; Class A/B/Cat-2 fills
#    7  Build dt.m                    pivot to long format; attach run_meta; save
#
#  Inputs
#    data/db_final_europe.csv
#    data/MA models template AGEE.xlsx
#
#  Outputs
#    results/d1.rds
#    results/ma_models.rds
#    results/dt_m.rds                (with run_meta attribute)
#    results/run_meta_script1.rds
#    logs/script_1_report.txt
#
#  References
#    Young et al. (2021)  doi:10.1016/j.agee.2021.107551
#    Young et al. (2022)  doi:10.1016/j.agee.2022.108229
#    Reference scripts: dst_loaddb.R, dst_prepare_input_fert.R,
#                       dst_functions_fert.R, dst_main_fert.R
# =============================================================================


# 0. Environment -------------------------------------------------------------

setwd('C:/Users/STMP2/Documents/Thesis/Github/Run/1')
rm(list = ls())

library(data.table)
library(readxl)
library(digest)

source('scripts/paths_original.R')              # P$* paths + require_input()
source(P$fn_fert)                      # Reference lmam, cIMAm, runDST

require_input(P$db)
require_input(P$ma_template)

run_ts     <- Sys.time()
report_con <- file(P$log_s1, open = 'wt')
log <- function(...) {
  msg <- paste0(...)
  cat(msg, '\n')
  cat(msg, '\n', file = report_con, append = TRUE)
}

log('=== Script 1 run: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')


# 1. Load db_final_europe.csv ------------------------------------------------
# The CSV is produced upstream by dst_loaddb.R and dst_prepare_input_fert.R
# and consolidates NCU-level inputs from INTEGRATOR, CAPRI, NCU covariates,
# NUTS spatial links, EUROSTAT 2016 FSS arable / management / manure data,
# and fertilizer-class downscaling (Young et al., 2022, Table 5).
#
# nup, nue and the 4R partial areas (parea.ee/.rft/.rfp/.rfr) are NOT in the
# CSV; they are derived in Step 4 below.

d1 <- fread(P$db)

log('')
log('--- Step 1: db_final_europe.csv loaded ---')
log('  Total rows:    ', nrow(d1))
log('  Unique NCUs:   ', length(unique(d1$ncu)))
log('  Total columns: ', ncol(d1))


# 2. Restore Young 2022 §2.1 yield convention --------------------------------
# yield_ref = current actual yield, yield_target = 80% water-limited potential.
# The raw CSV has the two columns transposed at source; this swap restores
# the convention required by Young et al. (2022) §2.1.

setnames(d1, c('yield_ref', 'yield_target'), c('yield_target', 'yield_ref'))

log('')
log('--- Step 2: Yield columns restored to Young 2022 §2.1 convention ---')
log('  mean(yield_ref)    = ',
    round(mean(d1$yield_ref,    na.rm = TRUE)), ' kg/ha (current actual)')
log('  mean(yield_target) = ',
    round(mean(d1$yield_target, na.rm = TRUE)), ' kg/ha (80% WLP)')


# 3. Verify required columns -------------------------------------------------

required_cols <- c(
  'ncu', 'country', 'crop', 'crop_code', 'area_ncu',
  'yield_ref', 'yield_target', 'soc_ref', 'soc_target',
  'n_sp_ref', 'n_sp_sw_crit', 'n_sp_gw_crit', 'density',
  'n_fert', 'n_man', 'n_fix', 'n_dep',
  'cov_soil', 'cov_clim', 'cov_crop', 'cov_fert', 'cov_soc',
  'area_ncu_ha',
  'parea.cc', 'parea.cres', 'parea.cr', 'parea.ntct', 'parea.rtct',
  'c_man_ncu',
  'parea.nifnof', 'parea.mifnof', 'parea.hifnof',
  'parea.mifhof', 'parea.hifhof', 'fert_type'
)
missing_cols <- setdiff(required_cols, colnames(d1))

log('')
log('--- Step 3: Column completeness check ---')
if (length(missing_cols) == 0) {
  log('  All ', length(required_cols), ' required columns present in d1.')
} else {
  log('  WARNING — missing columns: ', paste(missing_cols, collapse = ', '))
}

if (!'fr_rtnt' %in% colnames(d1)) {
  d1[, fr_rtnt := fifelse(
    (parea.ntct + parea.rtct) > 0,
    parea.ntct / (parea.ntct + parea.rtct),
    0
  )]
  log('  fr_rtnt absent — derived as parea.ntct / (parea.ntct + parea.rtct).')
}


# 4. Compute 4R partial areas ------------------------------------------------
# The reference regime is asymmetric: EE/RFT/RFP use the unclamped form
# (negative values pass through harmlessly to cIMAm's weighted.mean
# aggregation), while RFR alone clamps with pmax(0, .) because nue can
# exceed 1 in some NCUs.

d1[, nup       := n_fert + n_man + n_fix + n_dep - n_sp_ref]
d1[, nue       := nup / (n_fert + n_man + n_fix + n_dep)]
d1[, parea.ee  := pmin(1, yield_target / yield_ref - 1)        * area_ncu_ha]
d1[, parea.rft := pmin(1, yield_target / yield_ref - 1)        * area_ncu_ha]
d1[, parea.rfp := pmin(1, yield_target / yield_ref - 1)        * area_ncu_ha]
d1[, parea.rfr := pmin(1, pmax(0, 1 - nue))                    * area_ncu_ha]

log('')
log('--- Step 4: 4R partial areas computed ---')
for (col in c('parea.ee', 'parea.rft', 'parea.rfp', 'parea.rfr')) {
  log(sprintf(
    '  %-10s  min=%8.4f  mean=%9.2f  max=%9.2f  n_pos=%d  n_neg=%d',
    col,
    min(d1[[col]],  na.rm = TRUE),
    mean(d1[[col]], na.rm = TRUE),
    max(d1[[col]],  na.rm = TRUE),
    sum(d1[[col]] >  0, na.rm = TRUE),
    sum(d1[[col]] <  0, na.rm = TRUE)
  ))
}


# 5. Cross-check fertilizer-class partial areas ------------------------------

d1[, parea_fert_sum := parea.nifnof + parea.mifnof + parea.hifnof +
     parea.mifhof + parea.hifhof]
fert_check <- d1[, .(parea_fert_sum_ncu = parea_fert_sum[1],
                     area_ncu_ha_ncu    = area_ncu_ha[1]),
                 by = ncu]
fert_check[, ratio := parea_fert_sum_ncu / area_ncu_ha_ncu]

log('')
log('--- Step 5: fertilizer-class partial-area sum vs area_ncu_ha ---')
log('  Median ratio:           ', round(median(fert_check$ratio, na.rm = TRUE), 4))
log('  Mean ratio:             ', round(mean(  fert_check$ratio, na.rm = TRUE), 4))
log('  Min  ratio:             ', round(min(   fert_check$ratio, na.rm = TRUE), 4))
log('  Max  ratio:             ', round(max(   fert_check$ratio, na.rm = TRUE), 4))
log('  NCUs with ratio < 0.95: ', sum(fert_check$ratio < 0.95, na.rm = TRUE))
log('  NCUs with ratio > 1.05: ', sum(fert_check$ratio > 1.05, na.rm = TRUE))
log('  NCUs with NA   ratio:   ', sum(is.na(fert_check$ratio)))

d1[, parea_fert_sum := NULL]


# 6. Load and patch meta-analytical models -----------------------------------
# ma_cov_mean layout: (man_code, mods, group, mean_Y, mean_SOC, mean_Nsu)
# where `group` holds the covariate level (e.g. 'cereals', 'maize') and
# `mods` holds the covariate type (e.g. 'cov_crop', 'cov_soil').
#
# Two classes of patches are applied below:
#   6.1  Class A — replace lmam aggregates dominated by single non-European
#                  references with European LTE evidence (Zavattaro et al.,
#                  2014). Criterion: Zav n >= 10, |Δ lmam − Zav| >= 10 pp,
#                  and a single non-Zavattaro reference dominates the lmam
#                  contribution by inverse-variance weight.
#   6.1b Cat-2  — populate ma_cov_mean cells currently NA where Zavattaro
#                  rows exist with combined n >= 10 (the SEyr filter in
#                  lmam discards these because SE is missing).
# Plus published-value Class B fills for OF-MF and 4R Nsu, and the
# half-of-other rule for ROT/RT-CT N surplus (no meta-study model
# available; values are analyst imputations following the reference run).

ma_models <- lmam(fname = P$ma_template)

log('')
log('--- Step 6: ma_models loaded and patched ---')
log('  ma_mean man_codes:   ',
    paste(sort(unique(ma_models$ma_mean$man_code)), collapse = ', '))
log('  ma_cov_mean columns: ',
    paste(names(ma_models$ma_cov_mean), collapse = ', '))


# 6.1 Class A — Zavattaro European LTE overrides
ma_models$ma_cov_mean[man_code == 'ROT' & group == 'cereals', mean_Y :=  3.71]
ma_models$ma_cov_mean[man_code == 'RES' & group == 'cereals', mean_Y := -3.43]
ma_models$ma_cov_mean[man_code == 'RES' & group == 'maize',   mean_Y := -1.54]
ma_models$ma_cov_mean[man_code == 'CC'  & group == 'maize',   mean_Y :=  4.71]


# 6.1b Cat-2 fills — Zavattaro European LTE cells currently NA
# Five cells already exist in ma_cov_mean with another indicator populated;
# they are updated in place. Three cells are absent and appended as new rows.
# Treating both as appends would create duplicate (man_code, mods, group)
# keys and a cartesian-explosion error at cIMAm's left-join.

ma_models$ma_cov_mean[man_code == 'CC'  & mods == 'cov_crop' & group == 'cereals', mean_Y :=  8.52]
ma_models$ma_cov_mean[man_code == 'CC'  & mods == 'cov_soil' & group == 'coarse',  mean_Y :=  7.00]
ma_models$ma_cov_mean[man_code == 'CC'  & mods == 'cov_soil' & group == 'medium',  mean_Y := -0.82]
ma_models$ma_cov_mean[man_code == 'RES' & mods == 'cov_soil' & group == 'coarse',  mean_Y :=  6.50]
ma_models$ma_cov_mean[man_code == 'ROT' & mods == 'cov_soil' & group == 'coarse',  mean_Y := 11.00]

new_rows <- data.table(
  man_code = c('CF-MF',    'CF-MF',    'RES'    ),
  mods     = c('cov_soil', 'cov_soil', 'cov_soil'),
  group    = c('coarse',   'medium',   'medium' ),
  mean_Y   = c( 2.55,      -16.87,    -6.17    ),
  mean_SOC = c( 3.73,        2.55,    NA_real_ ),
  mean_Nsu = NA_real_
)
ma_models$ma_cov_mean <- rbind(ma_models$ma_cov_mean, new_rows, fill = TRUE)

stopifnot(
  "ma_cov_mean has duplicate (man_code, mods, group) keys" =
    !any(duplicated(ma_models$ma_cov_mean[, .(man_code, mods, group)]))
)


# 6.2 Class B — OF-MF N surplus (Young et al., 2022, Table 3)
# The Y row for OF-MF exists in the xlsx; the Nsu row is added here from
# the published global mean. SD set per Young et al. (2022) Table 3.
ma_models$ma_mean[man_code == 'OF-MF', mean_Nsu := -17.0]
ma_models$ma_sd  [man_code == 'OF-MF', sd_Nsu   :=   1.0]


# 6.2b Class B — 4R published global means (Young et al., 2022, Table 3)
# The xlsx `main` sheet has no rows for EE, RFT, RFP (covar = FALSE measures).
# These three rows are appended from the published Y22 Table 3 globals.
# Without this step, cIMAm() returns mmean = 0 for the 4R measures and
# they contribute nothing to the optimisation.
#
#              Yield (%)    SOC (% yr-1)    N surplus (%)
#   EE        +6.1 ± 0.3       n/a          -39 ± 1
#   RFT       +4.6 ± 0.4       n/a          -31 ± 1
#   RFP       +4.3 ± 0.3       n/a            n/a
fourR_main_rows <- data.table(
  man_code = c('EE',  'RFT',  'RFP'),
  mean_Y   = c( 6.1,    4.6,    4.3),
  mean_SOC = NA_real_,
  mean_Nsu = c(-39.0, -31.0,  NA_real_)
)
fourR_sd_rows <- data.table(
  man_code = c('EE',  'RFT',  'RFP'),
  sd_Y     = c( 0.3,    0.4,    0.3),
  sd_SOC   = NA_real_,
  sd_Nsu   = c( 1.0,    1.0,  NA_real_)
)
for (m in c('EE', 'RFT', 'RFP')) {
  if (!m %in% ma_models$ma_mean$man_code) {
    ma_models$ma_mean <- rbind(ma_models$ma_mean,
                               fourR_main_rows[man_code == m], fill = TRUE)
    ma_models$ma_sd   <- rbind(ma_models$ma_sd,
                               fourR_sd_rows[man_code == m], fill = TRUE)
  }
}
setkey(ma_models$ma_mean, man_code)
setkey(ma_models$ma_sd,   man_code)

log('  Class B 4R fills applied (EE, RFT, RFP).')
log('  ma_mean man_codes after patch: ',
    paste(sort(unique(ma_models$ma_mean$man_code)), collapse = ', '))


# 6.3 Half-of-other rule for ROT / RT-CT N surplus
# No meta-analytic model is available for these two cells (verified against
# both `main` and `covariates` sheets of the xlsx). The half-of-other rule
# (ROT Nsu = ½ RES Nsu; RT-CT Nsu = ½ NT-CT Nsu) is an analyst imputation
# applied in the reference Results compilation. We retain it here for
# numerical compatibility with the reference run. SD is set to NA because
# (i) the rule does not specify an SD treatment and (ii) msd_Nsu is unused
# downstream under nopt = FALSE / no Monte Carlo.
ma_models$ma_mean[man_code == 'ROT',
                  mean_Nsu := 0.5 * ma_models$ma_mean[man_code == 'RES',   mean_Nsu]]
ma_models$ma_sd  [man_code == 'ROT',   sd_Nsu := NA_real_]

ma_models$ma_mean[man_code == 'RT-CT',
                  mean_Nsu := 0.5 * ma_models$ma_mean[man_code == 'NT-CT', mean_Nsu]]
ma_models$ma_sd  [man_code == 'RT-CT', sd_Nsu := NA_real_]


# 7. Build dt.m via reference cIMAm() ----------------------------------------
# 10 measures, 3 indicators (Y, SOC, Nsu) -> 10 x 3 x 29,476 = 884,280 rows.
# EE/RFT/RFP use covar = FALSE (global 4R means); the other seven measures
# use covar = TRUE for covariate-disaggregated effects.
# RFR is intentionally excluded — matches the reference configuration in
# dst_main_fert.R.

dt.m1  <- cIMAm('EE',    d1, ma_models, covar = FALSE)
dt.m2  <- cIMAm('RFP',   d1, ma_models, covar = FALSE)
dt.m3  <- cIMAm('RFT',   d1, ma_models, covar = FALSE)
dt.m4  <- cIMAm('CF-MF', d1, ma_models, covar = TRUE )
dt.m5  <- cIMAm('OF-MF', d1, ma_models, covar = TRUE )
dt.m6  <- cIMAm('CC',    d1, ma_models, covar = TRUE )
dt.m7  <- cIMAm('NT-CT', d1, ma_models, covar = TRUE )
dt.m8  <- cIMAm('RES',   d1, ma_models, covar = TRUE )
dt.m9  <- cIMAm('ROT',   d1, ma_models, covar = TRUE )
dt.m10 <- cIMAm('RT-CT', d1, ma_models, covar = TRUE )

dt.m <- rbind(dt.m1, dt.m2, dt.m3, dt.m4, dt.m5,
              dt.m6, dt.m7, dt.m8, dt.m9, dt.m10)
rm(dt.m1, dt.m2, dt.m3, dt.m4, dt.m5,
   dt.m6, dt.m7, dt.m8, dt.m9, dt.m10)

n_neg  <- sum(dt.m$ha_m1 <  0, na.rm = TRUE)
n_na   <- sum(is.na(dt.m$ha_m1))
n_zero <- sum(dt.m$ha_m1 == 0, na.rm = TRUE)
neg_by_measure <- dt.m[ha_m1 < 0, .N, by = man_code]

log('')
log('--- Step 7: dt.m built ---')
log('  Total dt.m rows:     ', nrow(dt.m))
log('  Measures included:   ', paste(sort(unique(dt.m$man_code)), collapse = ', '))
log('  Negative ha_m1 rows: ', n_neg)
log('  Zero ha_m1 rows:     ', n_zero)
log('  NA   ha_m1 rows:     ', n_na)

# Verify 4R measures now have non-zero mmean values
fourR_check <- dt.m[man_code %in% c('EE', 'RFT', 'RFP'),
                    .(mean_mmean = mean(mmean, na.rm = TRUE),
                      n_NA       = sum(is.na(mmean))),
                    by = .(man_code, indicator)]
log('  4R mmean check:')
for (i in seq_len(nrow(fourR_check))) {
  log(sprintf('    %-4s %-4s  mean=%7.4f  n_NA=%d',
              fourR_check$man_code[i], fourR_check$indicator[i],
              fourR_check$mean_mmean[i], fourR_check$n_NA[i]))
}

log('  Negative ha_m1 per measure:')
if (nrow(neg_by_measure) > 0) {
  for (i in seq_len(nrow(neg_by_measure))) {
    log(sprintf('    %-6s  n_neg = %d',
                neg_by_measure$man_code[i], neg_by_measure$N[i]))
  }
} else {
  log('    none')
}


# 8. Save outputs with reproducibility metadata ------------------------------

finished_ts <- Sys.time()
run_meta <- list(
  script         = 'Script_1_Data_Preparation.R',
  finished_ts    = finished_ts,
  run_ts         = finished_ts,            # alias used by Script 2a
  measures       = sort(unique(dt.m$man_code)),
  n_ncu          = length(unique(d1$ncu)),
  n_dt_m_rows    = nrow(dt.m),
  negative_ha_m1 = list(
    n_total    = n_neg,
    by_measure = neg_by_measure,
    min_value  = min(dt.m$ha_m1, na.rm = TRUE),
    note       = paste('Negative ha_m1 in EE/RFT/RFP is reference',
                       '(NCUs above 80% WLP target).')
  ),
  input_hashes   = list(
    db          = digest(file = P$db,          algo = 'md5'),
    ma_template = digest(file = P$ma_template, algo = 'md5')
  )
)

attr(dt.m, 'run_meta') <- run_meta

saveRDS(d1,        P$d1)
saveRDS(ma_models, P$ma_models)
saveRDS(dt.m,      P$dt_m)
saveRDS(run_meta,  P$meta_s1)

log('')
log('--- Step 8: Outputs saved ---')
log('  ', P$d1,        '  (', nrow(d1), ' rows, ', ncol(d1), ' columns)')
log('  ', P$ma_models, '  (patched ma_mean + ma_cov_mean)')
log('  ', P$dt_m,      '  (', nrow(dt.m), ' rows, ',
    length(unique(dt.m$man_code)), ' measures, run_meta attached)')
log('  ', P$meta_s1,   '  (staleness verification metadata)')
log('  db_final_europe.csv MD5: ', substr(run_meta$input_hashes$db, 1, 16), '...')

log('')
log('=== Script 1 complete: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')
close(report_con)


# End of Script 1 ------------------------------------------------------------
