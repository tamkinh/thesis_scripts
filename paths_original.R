# =============================================================================
#  scripts/paths_original.R — REFERENCE PROJECT PATHS
# =============================================================================
#
#  Purpose
#  -------
#  Single source of truth for all input, output, and log paths used across
#  the DSF pipeline (Scripts 1, 2a, 2b, 3, 4, 5). Every analysis script begins
#  with `source('scripts/paths.R')` and references paths via the `P` list
#  (e.g. `readRDS(P$dt_m)` rather than `readRDS('results/dt_m.rds')`).
#
#  Why
#  ---
#  Hard-coded paths scattered across scripts created confusion (May 2026:
#  Script 1 wrote to outputs/, Scripts 2a/2b/3 read from results/, dual-
#  write workarounds proliferated). This file eliminates that class of
#  bug. Changing a folder name now requires editing exactly one line here.
#
#  Convention
#  ----------
#    data/      raw inputs                (read-only, never written to)
#    scripts/   all R source files
#    results/   all .rds outputs from any script (single canonical location)
#    maps/      figures (PNG/PDF) from Script 3
#    logs/      text reports from each script's run
#
#  Usage
#  -----
#    setwd('C:/Users/STMP2/Documents/Thesis/Github/Run/1')
#    source('scripts/paths_original.R')
#    d1 <- readRDS(P$d1)
#    saveRDS(out.best, P$out_best)
# =============================================================================


# ---- Directory roots ------------------------------------------------------

DATA_DIR    <- 'data/'
RESULTS_DIR <- 'results/'
MAPS_DIR    <- 'maps/'
LOG_DIR     <- 'logs/'
SCRIPTS_DIR <- 'scripts/'


# ---- Ensure output directories exist at runtime ---------------------------
# data/ and scripts/ must exist before any script runs (raw inputs +
# source files). The three output-side directories are created on demand.

for (.d in c(RESULTS_DIR, MAPS_DIR, LOG_DIR)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}
rm(.d)


# ---- Canonical artifact paths (P list) ------------------------------------
# Reference paths symbolically via P$<name> in every script. NEVER hard-
# code 'results/x.rds' or 'data/y.csv' in any analysis script.

P <- list(

  # ---- Raw inputs (data/) --------------------------------------------
  db          = paste0(DATA_DIR,    'db_final_europe.csv'),
  ma_template = paste0(DATA_DIR,    'MA models template AGEE.xlsx'),
  ncu_raster  = paste0(DATA_DIR,    'gncu2010_ext.asc'),

  # ---- Script 1 outputs (results/) -----------------------------------
  d1          = paste0(RESULTS_DIR, 'd1.rds'),
  ma_models   = paste0(RESULTS_DIR, 'ma_models.rds'),
  dt_m        = paste0(RESULTS_DIR, 'dt_m.rds'),
  meta_s1     = paste0(RESULTS_DIR, 'run_meta_script1.rds'),

  # ---- Script 2a outputs ---------------------------------------------
  out_best      = paste0(RESULTS_DIR, 'out_best_n1.rds'),      # Fig 7 input; best single per NCU
  out_runnerup  = paste0(RESULTS_DIR, 'out_runnerup_n1.rds'),  # 2nd-ranked single (bipmcs == 2); SRQ2 selection margin
  out_allscores = paste0(RESULTS_DIR, 'out_allscores_n1.rds'), # full single-measure score ladder; SRQ2 'none'-class diagnostics
  out_duo       = paste0(RESULTS_DIR, 'out_duo_n2.rds'),       # Fig 8 input
  out_trio      = paste0(RESULTS_DIR, 'out_trio_n3.rds'),      # Fig 9 input
  meta_s2a      = paste0(RESULTS_DIR, 'run_meta_script2a.rds'),

  # ---- Script 2b outputs ---------------------------------------------
  table41_stats = paste0(RESULTS_DIR, 'table41_stats.csv'),    # compiled Table 4.1 (26 rows)
  meta_s2b      = paste0(RESULTS_DIR, 'run_meta_script2b.rds'),

  # ---- Script 3 outputs ----------------------------------------------
  meta_s3       = paste0(RESULTS_DIR, 'run_meta_script3.rds'),

  # ---- Script 4 outputs ----------------------------------------------
  driver_ncu    = paste0(RESULTS_DIR, 'srq2_driver_by_ncu.csv'),
  meta_s4       = paste0(RESULTS_DIR, 'run_meta_script4.rds'),

  # ---- Script 5 outputs ----------------------------------------------
  nuts2_agg     = paste0(RESULTS_DIR, 'nuts2_aggregated.csv'),
  country_agg   = paste0(RESULTS_DIR, 'country_aggregated.csv'),
  nd_wfd        = paste0(RESULTS_DIR, 'nd_wfd_assessment.csv'),
  meta_s5       = paste0(RESULTS_DIR, 'run_meta_script5.rds'),

  # ---- Log files (logs/) ---------------------------------------------
  log_s1      = paste0(LOG_DIR,     'script_1_report.txt'),
  log_s2a     = paste0(LOG_DIR,     'script_2a_report.txt'),
  log_s2b     = paste0(LOG_DIR,     'script_2b_report.txt'),
  log_s3      = paste0(LOG_DIR,     'script_3_report.txt'),
  log_s4      = paste0(LOG_DIR,     'script_4_report.txt'),
  log_s5      = paste0(LOG_DIR,     'script_5_report.txt'),

  # ---- Sourced function files (scripts/) -----------------------------
  fn_fert      = paste0(SCRIPTS_DIR, 'dst_functions_fert.R')
)

# ---- Helper: confirm an input file exists or stop ------------------------
# Use at the top of any script after sourcing paths.R:
#   require_input(P$db); require_input(P$ma_template)

require_input <- function(path) {
  if (!file.exists(path)) {
    stop('Required input file not found: ', path,
         '\nCheck paths.R and file location.', call. = FALSE)
  }
  invisible(TRUE)
}


# End of paths.R -----------------------------------------------------------
