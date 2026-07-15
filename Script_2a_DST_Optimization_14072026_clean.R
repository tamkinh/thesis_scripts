# =============================================================================
# SCRIPT 2a — DST OPTIMISATION (best single, duo, trio)
# =============================================================================
#  Purpose
#    Run the three DST optimisations (nmax = 1, 2, 3) under the reference
#    configuration. Outputs feed Script 3 (Figs 7-9) and Script 2b (Table 4.1
#    rows f, g, h).
#
#  Sections
#    0  Environment         libraries, paths, log connection
#    1  Staleness check     hash-verify db against Script 1 metadata
#    2  Load data           d1 and dt.m from Script 1
#    3  Helpers             add_gap_columns(); enrich_for_gap_cols();
#                           d1_ncu area-weighted NCU aggregation
#    4  DST runs            best single (nmax=1), duo (nmax=2), trio (nmax=3)
#                           under reference uw = c(1,1,1), simyear = 5
#    5  Save outputs        RDS files with run_meta attribute; log summary
#
#  Inputs
#    results/d1.rds                from Script 1
#    results/dt_m.rds              from Script 1
#    results/run_meta_script1.rds  staleness verification
#
#  Outputs
#    results/out_best_n1.rds       Fig 7 input
#    results/out_runnerup_n1.rds   2nd-best single measure scores (Script 4, SRQ2)
#    results/out_allscores_n1.rds  all single-measure scores (optional; only if
#                                  runDST returns impact_allscores)
#    results/out_duo_n2.rds        Fig 8 input
#    results/out_trio_n3.rds       Fig 9 input
#    results/run_meta_script2a.rds run metadata
#    logs/script_2a_report.txt
#
#  Uses dst_main_fert.R (June 2024 run):
#    runDST(db = d1, dt.m = dt.m, output = 'best_impact',
#           uw = c(1,1,1), simyear = 5, nmax = 1, nopt = FALSE)
#
#  Gap-metric formulas (dst_main_fert.R, lines 380-429):
#    yield_gap_diff_p = (yield_gap_fin_t - yield_gap_t) / abs(yield_gap_t) * 100
#    soc_gap_diff_p   = (soc_gap_fin_t   - soc_gap_t)   / abs(soc_gap_t)   * 100
#    n_sp_gap_diff_p  = (n_sp_gap_fin    - n_sp_gap)    / abs(n_sp_gap)    * 100
#  ±Inf from zero-denominator NCUs collapses to 0 (dst_plotfigures.R L555-556).
#
#  References
#    Young et al. (2021, 2022); dst_main_fert.R, dst_functions_fert.R
# =============================================================================


# 0. Environment -------------------------------------------------------------

setwd('C:/Users/STMP2/Documents/Thesis/Github/Run/1')
rm(list = ls())

library(data.table)
library(digest)

source('scripts/paths_original.R')
source(P$fn_fert)

run_ts <- Sys.time()

require_input(P$d1)
require_input(P$dt_m)
require_input(P$meta_s1)

report_con <- file(P$log_s2a, open = 'wt')
log <- function(...) {
  msg <- paste0(...)
  cat(msg, '\n')
  tryCatch(
    cat(msg, '\n', file = report_con, append = TRUE),
    error = function(e) invisible(NULL)
  )
}

log('=== Script 2a run: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')


# 1. Staleness verification --------------------------------------------------

script1_meta <- readRDS(P$meta_s1)
log('Script 1 finished : ',
    format(script1_meta$finished_ts, '%Y-%m-%d %H:%M:%S'))

db_hash_now <- digest::digest(file = P$db, algo = 'md5')
if (db_hash_now != script1_meta$input_hashes[['db']])
  stop('Database has changed since Script 1 ran. Re-run Script 1.')


# 2. Load Script 1 outputs ---------------------------------------------------

d1        <- readRDS(P$d1)
dt.m      <- readRDS(P$dt_m)
dt_m_meta <- attr(dt.m, 'run_meta')

log('  dt.m origin run  : ',
    format(dt_m_meta$run_ts, '%Y-%m-%d %H:%M:%S'))
log('  dt.m rows        : ', nrow(dt.m))
log('  dt.m measures    : ', paste(dt_m_meta$measures, collapse = ', '))


# 3. Helpers -----------------------------------------------------------------

# Add reference gap columns to a DST output (dst_main_fert.R, L380-429).
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

  for (col in c('yield_gap_diff_p', 'soc_gap_diff_p', 'n_sp_gap_diff_p')) {
    out[is.infinite(get(col)), (col) := 0]
    out[is.nan(get(col)),      (col) := 0]
    out[is.na(get(col)),       (col) := 0]
  }

  out
}

# Aggregate d1 to one row per NCU using area-weighted means within the NCU,
# matching how runDST aggregates internally (dst_functions_fert.R, L105).
# Used to enrich the leaner duo/trio outputs with target/density columns.
d1_ncu <- d1[, .(
  yield_ref    = weighted.mean(yield_ref,    area_ncu, na.rm = TRUE),
  yield_target = weighted.mean(yield_target, area_ncu, na.rm = TRUE),
  soc_ref      = weighted.mean(soc_ref,      area_ncu, na.rm = TRUE),
  soc_target   = weighted.mean(soc_target,   area_ncu, na.rm = TRUE),
  n_sp_ref     = weighted.mean(n_sp_ref,     area_ncu, na.rm = TRUE),
  n_sp_crit    = pmin(weighted.mean(n_sp_sw_crit, area_ncu, na.rm = TRUE),
                      weighted.mean(n_sp_gw_crit, area_ncu, na.rm = TRUE)),
  bd           = weighted.mean(density,      area_ncu, na.rm = TRUE)
), by = ncu]

enrich_for_gap_cols <- function(out) {
  need <- c('yield_target', 'soc_target', 'n_sp_crit', 'bd',
            'yield_ref', 'soc_ref', 'n_sp_ref')
  miss <- setdiff(need, names(out))
  if (length(miss) > 0) {
    out <- merge(out,
                 d1_ncu[, c('ncu', miss), with = FALSE],
                 by = 'ncu', all.x = TRUE)
  }
  out
}


# 4. DST runs ----------------------------------------------------------------
# quiet = FALSE so the long-running optimisation shows progress to console.

log('')
log('--- Running DST nmax = 1 (best single) ---')
sim.n1 <- runDST(db = d1, dt.m = dt.m, output = 'best_impact',
                 uw = c(1, 1, 1), simyear = 5, quiet = FALSE,
                 nmax = 1, nopt = FALSE)
out.best <- add_gap_columns(sim.n1$impact_best)
# Runner-up (2nd-best single) scores for the SRQ2 selection-margin decomposition,
# and, if exported, the full single-measure score table for the extended variants
# (margin against the next distinct alternative; leave-one-out). Neither needs the
# gap columns added by add_gap_columns().
out.best.ru  <- sim.n1$impact_runnerup
out.best.all <- sim.n1$impact_allscores

log('Best single frequency:')
print(sort(table(out.best$man_code), decreasing = TRUE))

log('')
log('--- Running DST nmax = 2 (best duo) ---')
sim.n2 <- runDST(db = d1, dt.m = dt.m, output = 'score_duo',
                 uw = c(1, 1, 1), simyear = 5, quiet = FALSE,
                 nmax = 2, nopt = FALSE)
out.duo <- add_gap_columns(enrich_for_gap_cols(sim.n2$score_duo))

log('Top 10 duos:')
print(sort(table(out.duo$man_code), decreasing = TRUE)[1:10])


log('')
log('--- Running DST nmax = 3 (best trio) ---')
sim.n3 <- runDST(db = d1, dt.m = dt.m, output = 'score_trio',
                 uw = c(1, 1, 1), simyear = 5, quiet = FALSE,
                 nmax = 3, nopt = FALSE)
out.trio <- add_gap_columns(enrich_for_gap_cols(sim.n3$score_trio))

log('Top 10 trios:')
print(sort(table(out.trio$man_code), decreasing = TRUE)[1:10])


# 5. Save outputs ------------------------------------------------------------

# Derive the SRQ2 export paths from P$out_best when they are not defined in the
# paths file, so the script runs against either version of it.
if (is.null(P$out_runnerup))
  P$out_runnerup  <- sub('out_best_n1', 'out_runnerup_n1',  P$out_best)
if (is.null(P$out_allscores))
  P$out_allscores <- sub('out_best_n1', 'out_allscores_n1', P$out_best)

run_meta <- list(
  script              = 'Script_2a_DST_Optimization.R',
  started_ts          = run_ts,
  finished_ts         = Sys.time(),
  r_version           = R.version.string,
  session_info        = sessionInfo(),
  script1_finished_ts = script1_meta$finished_ts,
  input_hashes        = list(
    d1   = digest::digest(file = P$d1,   algo = 'md5'),
    dt_m = digest::digest(file = P$dt_m, algo = 'md5')
  ),
  out_best_n          = nrow(out.best),
  out_duo_n           = nrow(out.duo),
  out_trio_n          = nrow(out.trio),
  out_runnerup_n      = nrow(out.best.ru),
  out_allscores_n     = if (is.null(out.best.all)) NA_integer_ else nrow(out.best.all),
  uw                  = c(1, 1, 1),
  simyear             = 5,
  nopt                = FALSE,
  functions_src       = 'dst_functions_fert.R'
)

attr(out.best, 'run_meta') <- run_meta
attr(out.duo,  'run_meta') <- run_meta
attr(out.trio, 'run_meta') <- run_meta

saveRDS(out.best, P$out_best)
saveRDS(out.best.ru, P$out_runnerup)
if (!is.null(out.best.all)) {
  saveRDS(out.best.all, P$out_allscores)
} else {
  log('NOTE: runDST did not return impact_allscores - all-scores export skipped ',
      '(only needed for the extended SRQ2 variants).')
}
saveRDS(out.duo,  P$out_duo)
saveRDS(out.trio, P$out_trio)
saveRDS(run_meta, P$meta_s2a)

log('Saved: out_best_n1 (', nrow(out.best), ' rows), out_runnerup_n1 (',
    nrow(out.best.ru), ' rows), out_duo_n2 (', nrow(out.duo), ' rows), out_trio_n3 (',
    nrow(out.trio), ' rows)')

log('')
log('=== Script 2a complete: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')
log('Duration: ', sprintf('%.1f min',
                          as.numeric(difftime(Sys.time(), run_ts, units = 'mins'))))
close(report_con)


# End of Script 2a -----------------------------------------------------------
