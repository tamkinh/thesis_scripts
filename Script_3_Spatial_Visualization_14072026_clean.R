# =============================================================================
# SCRIPT 3 — SPATIAL VISUALISATION (reference maps)
# =============================================================================
#  Purpose
#    Generate the nine publication maps from Script 2a outputs. Palette,
#    breaks and projection are those found in dst_plotfigures.R.
#
#  Sections
#    0  Environment         libraries, paths, paletteer, theme
#    1  Staleness check     hash-verify db and Script 2a chain
#    2  Measure definitions  man_num map, combo parser, Fig 8/9 top-15 lists
#    3  Plotting helpers     visualize_discrete2() and visualize_discrete()
#                            following dst_plotfigures.R conventions
#    4  Load data            out_best/duo/trio joined to NCU raster
#
#    Reproduced reference maps (Appendix B, Fig B.1 and B.3):
#
#    5  Figs 1–3             actual indicator values (YlGn, Greys, RdPu)
#    6  Figs 4–6             relative change in indicator gap (Spectral)
#    7  Fig 7                best single measure (Tableau_10)
#    8  Fig 8                best two combined measures (Classic_20, top 15)
#    9  Fig 9                best three combined measures (Classic_20, top 15)
#
#  Inputs
#    results/out_best_n1.rds         from Script 2a
#    results/out_duo_n2.rds          from Script 2a
#    results/out_trio_n3.rds         from Script 2a
#    results/run_meta_script2a.rds   staleness verification
#    data/gncu2010_ext.asc           NCU raster (EPSG:3035)
#
#  Outputs (in maps/)
#    Fig1_yield_ref.png              Actual crop yield
#    Fig2_soc_ref.png                Actual SOC stock
#    Fig3_n_surplus_ref.png          Actual N surplus
#    Fig4_yield_gap_change.png       Relative change in yield gap
#    Fig5_soc_gap_change.png         Relative change in SOC gap
#    Fig6_n_gap_change.png           Relative change in N surplus gap
#    Fig7_best_single.png            Best management measure
#    Fig8_best_duo.png               Best two combined measures
#    Fig9_best_trio.png              Best three combined measures
#
#  Conventions
#    Figs 1-6   visualize_discrete2() with RColorBrewer palettes
#    Fig 7      visualize_discrete()  with Tableau_10. man_num 1-10 follows
#               the canonical assignment in dst_plotfigures.R (L271-281).
#    Figs 8-9   build_combo_raster() ranks combinations by frequency in the
#               current run and assigns Classic_20 colours to the top 15;
#               the remainder collapse to 'other' (bin 16). man_code strings
#               are canonicalised (components sorted alphabetically) before
#               counting so that 'EE-CC' and 'CC-EE' are treated as identical.
#
#  References
#    Young et al. (2022); dst_plotfigures.R
# =============================================================================


# 0. Environment -------------------------------------------------------------

setwd('C:/Users/STMP2/Documents/Thesis/Github/Run/1')
rm(list = ls())

library(data.table)
library(digest)
library(terra)
library(ggplot2)
library(paletteer)

source('scripts/paths_original.R')

require_input(P$out_best)
require_input(P$out_duo)
require_input(P$out_trio)
require_input(P$meta_s2a)
require_input(P$meta_s1)

if (!dir.exists('maps')) dir.create('maps')
theme_set(theme_bw())

run_ts <- Sys.time()
cat('=== Script 3 run: ', format(run_ts, '%Y-%m-%d %H:%M:%S'), ' ===\n')


# 1. Staleness verification --------------------------------------------------

script2a_meta <- readRDS(P$meta_s2a)
script1_meta  <- readRDS(P$meta_s1)

if (script2a_meta$script1_finished_ts != script1_meta$finished_ts)
  stop('Script 2a was not based on the current Script 1 output. Re-run Script 2a.')

db_hash_now <- digest::digest(file = P$db, algo = 'md5')
if (db_hash_now != script2a_meta$db_hash)
  stop('Database has changed since Script 2a ran. Re-run Scripts 1 and 2a.')


# 2. Measure definitions, combination parser, and Fig 8/9 constants ----------

MEASURES <- c('CF-MF', 'NT-CT', 'OF-MF', 'RT-CT',
              'CC', 'EE', 'RES', 'RFP', 'RFT', 'ROT')

ABBREV <- c('CC'    = 'CC',  'CF-MF' = 'CF',  'EE'  = 'EE',
            'NT-CT' = 'NT',  'OF-MF' = 'OF',  'RES' = 'RES',
            'RFP'   = 'RFP', 'RFT'   = 'RFT', 'ROT' = 'ROT',
            'RT-CT' = 'RT')

# Parse hyphen-concatenated man_code into components. Longest-first matching
# handles internal hyphens in CF-MF / OF-MF / NT-CT / RT-CT.
parse_combo <- function(combo) {
  ms    <- MEASURES[order(-nchar(MEASURES))]
  parts <- character(0)
  rem   <- combo
  while (nchar(rem) > 0) {
    hit <- FALSE
    for (m in ms) {
      pat <- paste0('^', gsub('-', '\\-', m, fixed = TRUE))
      if (grepl(pat, rem)) {
        parts <- c(parts, m)
        rem   <- sub(paste0(pat, '-?'), '', rem)
        hit   <- TRUE
        break
      }
    }
    if (!hit) break
  }
  parts
}

canonicalise_code <- function(code) {
  parts <- parse_combo(code)
  if (length(parts) == 0) return(code)
  paste(sort(parts), collapse = '-')
}

make_label <- function(combo, sep) {
  parts <- parse_combo(combo)
  if (length(parts) == 0) return(combo)
  paste(ABBREV[parts], collapse = sep)
}

# Vectorised wrapper around canonicalise_code() — alphabetises components
# of each man_code so 'EE-CC' and 'CC-EE' are treated identically.
# Used by build_combo_raster() to deduplicate combinations before counting.
normalise_combo <- function(codes) {
  vapply(codes, canonicalise_code, character(1), USE.NAMES = FALSE)
}

# Hardcoded top-15 combination lists from canonical dst_plotfigures.R L368-400.
# Order = man_num assignment for Reference 8 / Reference 9 legend matching.
# Used when mode = 'hardcoded' in build_combo_raster().
HARDCODED_DUO_15 <- c(
  'CC-EE',    'EE-RES',   'EE-RFT',     'CC-RES',  'RES-ROT',
  'EE-ROT',   'CC-ROT',   'NT-CT-RES',  'RFT-ROT', 'OF-MF-ROT',
  'CC-RFT',   'RFP-ROT',  'RES-RFT',    'EE-RFP',  'RES-RT-CT'
)

HARDCODED_TRIO_15 <- c(
  'CC-EE-RES',    'EE-RES-RFT',     'EE-RFP-RFT',  'CC-EE-RFT',
  'EE-RES-ROT',   'CC-EE-ROT',      'EE-RFT-ROT',  'CC-NT-CT-RES',
  'NT-CT-RES-ROT','CC-RES-ROT',     'EE-NT-CT-RES','RES-ROT-RT-CT',
  'CC-CF-MF-RES', 'EE-RFP-ROT',     'EE-RES-RT-CT'
)

# 3. Plotting helpers --------------------------------------------------------

visualize_discrete2 <- function(raster, layer, mapcolor, direction,
                                name, breaks, labels, ftitle) {
  df <- as.data.frame(raster[[layer]], xy = TRUE)
  colnames(df) <- c('x', 'y', 'variable')
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y,
                             fill = cut(variable, breaks, labels = labels))) +
    coord_sf(crs = 3035, lims_method = 'box') +
    scale_fill_brewer(palette      = mapcolor,
                      na.translate = TRUE,
                      na.value     = 'white',
                      drop         = FALSE,
                      direction    = direction) +
    xlab('') + ylab('') + labs(fill = name) +
    theme(text              = element_text(size = 28),
          legend.text       = element_text(size = 20),
          legend.position   = c(0.19, 0.85),
          legend.background = element_rect(fill = 'white', color = 'white'),
          panel.border      = element_blank(),
          plot.title        = element_text(hjust = 0.5)) +
    ggtitle(ftitle)
}

visualize_discrete <- function(raster, layer, mapcolor, name,
                               breaks, labels, leg_cols, ftitle) {
  df <- as.data.frame(raster[[layer]], xy = TRUE)
  colnames(df) <- c('x', 'y', 'variable')
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y,
                             fill = cut(variable, breaks, labels = labels))) +
    coord_sf(crs = 3035, lims_method = 'box') +
    mapcolor +
    xlab('') + ylab('') + labs(fill = name) +
    theme(text              = element_text(size = 24),
          legend.text       = element_text(size = 16),
          legend.position   = c(0.28, 0.88),
          legend.background = element_rect(fill = 'white', color = 'white'),
          panel.border      = element_blank(),
          plot.title        = element_text(hjust = 0.5)) +
    guides(fill = guide_legend(ncol = leg_cols)) +
    ggtitle(ftitle)
}

# -----------------------------------------------------------------------------
# build_combo_raster — supports two label-set modes for Figs 8 and 9
#
# 'hardcoded': uses Young (2022) / dst_plotfigures.R L368-400 published
#              top-15 lists. Each NCU's actual best combination is matched
#              to this list; non-matches collapse to bin 16 ('other').
#              Use this mode to reproduce Reference 8 / Reference 9 legends.
#
# 'dynamic':   ranks combinations by frequency in the current run and
#              picks the top 15. Use this for the most faithful empirical
#              representation of this run's distribution.
# -----------------------------------------------------------------------------
build_combo_raster <- function(r1.p, out.combo, sep,
                               mode = 'hardcoded', nmax) {

  stopifnot(mode %in% c('hardcoded', 'dynamic'))
  stopifnot(nmax %in% c(2L, 3L))

  out.combo <- copy(out.combo)
  out.combo[, man_code_norm := normalise_combo(man_code)]

  if (mode == 'hardcoded') {
    hardlist      <- if (nmax == 2L) HARDCODED_DUO_15 else HARDCODED_TRIO_15
    hardlist_norm <- normalise_combo(hardlist)
    labels16      <- c(sapply(hardlist, make_label, sep = sep), 'other')
    top_or_hard   <- hardlist
  } else {  # dynamic
    freq          <- out.combo[, .N, by = man_code_norm][order(-N)]
    hardlist_norm <- freq[seq_len(min(15L, nrow(freq))), man_code_norm]
    lk            <- unique(out.combo[, .(man_code_norm, man_code)])
    top_raw       <- lk[match(hardlist_norm, man_code_norm), man_code]
    labels16      <- c(sapply(top_raw, make_label, sep = sep), 'other')
    top_or_hard   <- top_raw
  }

  r.combo <- merge(r1.p,
                   out.combo[, .(ncu, man_code_norm)],
                   by.x = 'gncu2010_ext', by.y = 'ncu')

  r.combo[, man_num := 16L]
  for (i in seq_along(hardlist_norm)) {
    r.combo[man_code_norm == hardlist_norm[i], man_num := i]
  }

  r.vals <- r.combo[, .(x, y, man_num)]
  r.fin  <- terra::rast(r.vals, type = 'xyz')
  terra::crs(r.fin) <- 'epsg:3035'

  list(raster = r.fin, labels = labels16, top15 = top_or_hard, mode = mode)
}


# 4. Load data and build reference raster ------------------------------------

out.best <- readRDS(P$out_best)
out.duo  <- readRDS(P$out_duo)
out.trio <- readRDS(P$out_trio)

cat(sprintf('Inputs loaded: out.best %d rows, out.duo %d rows, out.trio %d rows\n',
            nrow(out.best), nrow(out.duo), nrow(out.trio)))

r1   <- terra::rast('data/gncu2010_ext.asc')
terra::crs(r1) <- 'epsg:3035'
r1.p <- as.data.table(as.data.frame(r1, xy = TRUE))

ref.cols <- out.best[, .(ncu,
                         yield_ref_t,    yield_targ_t,   yield_new_t,
                         soc_ref_t,      soc_targ_t,     soc_new_t,
                         n_sp_ref,       n_sp_crit,      n_sp_new,
                         yield_gap_diff_p, soc_gap_diff_p, n_sp_gap_diff_p)]

r.ref <- merge(r1.p, ref.cols, by.x = 'gncu2010_ext', by.y = 'ncu')
setcolorder(r.ref, c('x', 'y', 'gncu2010_ext',
                     'yield_ref_t',    'yield_targ_t',   'yield_new_t',
                     'soc_ref_t',      'soc_targ_t',     'soc_new_t',
                     'n_sp_ref',       'n_sp_crit',      'n_sp_new',
                     'yield_gap_diff_p', 'soc_gap_diff_p', 'n_sp_gap_diff_p'))

r.fin.ref <- terra::rast(r.ref, type = 'xyz')
terra::crs(r.fin.ref) <- 'epsg:3035'


# 5. Figures 1-3 — Current indicator values ---------------------------------

p <- visualize_discrete2(r.fin.ref, 'yield_ref_t', 'YlGn', 1,
                         expression('Mg ha'^-1),
                         c(0, 4, 6, 8, 10, 15, 55),
                         c('< 4', '4 - 6', '6 - 8', '8 - 10', '10 - 15', '> 15'),
                         'Actual crop yield')
ggsave('maps/Fig1_yield_ref.png', p, width = 25, height = 25,
       units = 'cm', dpi = 400)

p <- visualize_discrete2(r.fin.ref, 'soc_ref_t', 'Greys', 1,
                         expression('Mg ha'^-1),
                         c(0, 30, 45, 55, 65, 75, 2000),
                         c('< 30', '30 - 45', '45 - 55', '55 - 65', '65 - 75', '> 75'),
                         'Actual SOC stock')
ggsave('maps/Fig2_soc_ref.png', p, width = 25, height = 25,
       units = 'cm', dpi = 400)

p <- visualize_discrete2(r.fin.ref, 'n_sp_ref', 'RdPu', 1,
                         expression('kg ha'^-1),
                         c(-500, 24, 31, 38, 45, 55, 500),
                         c('< 24', '24 - 31', '31 - 38', '38 - 45', '45 - 55', '> 55'),
                         'Actual N surplus')
ggsave('maps/Fig3_n_surplus_ref.png', p, width = 25, height = 25,
       units = 'cm', dpi = 400)

cat('Figs 1-3 saved\n')


# 6. Figures 4-6 — Relative change in indicator gap -------------------------

p <- visualize_discrete2(r.fin.ref, 'yield_gap_diff_p', 'Spectral', -1,
                         '% of gap',
                         c(-2.7e+07, -150, -100, -50, -25, -5, 0, 5, 1.5e+07),
                         c('< -150', '-150 to -100', '-100 to -50', '-50 to -25',
                           '-25 to -5', '-5 to 0', '0 to 5', '> 5'),
                         'Relative change in crop yield gap')
ggsave('maps/Fig4_yield_gap_change.png', p, width = 25, height = 25,
       units = 'cm', dpi = 400)

p <- visualize_discrete2(r.fin.ref, 'soc_gap_diff_p', 'Spectral', -1,
                         '% of gap',
                         c(-200000, -5, -2, -1, 0, 1, 2, 5, 100000),
                         c('< -5', '-5 to -2', '-2 to -1', '-1 to 0',
                           '0 to 1', '1 to 2', '2 to 5', '> 5'),
                         'Relative change in SOC gap')
ggsave('maps/Fig5_soc_gap_change.png', p, width = 25, height = 25,
       units = 'cm', dpi = 400)

p <- visualize_discrete2(r.fin.ref, 'n_sp_gap_diff_p', 'Spectral', -1,
                         '% of gap',
                         c(-2.5e+06, -300, -150, -100, -50, 0, 5, 10, 100),
                         c('< -300', '-300 to -150', '-150 to -100', '-100 to -50',
                           '-50 to 0', '0 to 5', '5 to 10', '> 10'),
                         'Relative change in N surplus gap')
ggsave('maps/Fig6_n_gap_change.png', p, width = 25, height = 25,
       units = 'cm', dpi = 400)

cat('Figs 4-6 saved\n')


# 7. Figure 7 — Best single management measure ------------------------------
# man_num assignment per dst_plotfigures.R L271-281.

r.man <- out.best[, .(ncu, area_ncu_ha_tot, man_code, bipmc)]
r.man <- merge(r1.p, r.man, by.x = 'gncu2010_ext', by.y = 'ncu')

man_num_map <- c('CF-MF' = 1, 'OF-MF' = 2, 'EE'  = 3, 'RFT' = 4, 'RFP' = 5,
                 'RT-CT' = 6, 'NT-CT' = 7, 'ROT' = 8, 'CC'  = 9, 'RES' = 10)
r.man[, man_num := man_num_map[man_code]]

setcolorder(r.man, c('x', 'y', 'gncu2010_ext', 'man_code'))
r.fin.man <- terra::rast(r.man, type = 'xyz')
terra::crs(r.fin.man) <- 'epsg:3035'

p <- visualize_discrete(r.fin.man, 'man_num',
                        paletteer::scale_fill_paletteer_d('ggthemes::Tableau_10',
                                                          direction = -1),
                        '',
                        seq(0.5, 10.5, by = 1),
                        c('combined fert. (CF)',
                          'organic fert. (OF)',
                          'enhanced efficiency fert. (EE)',
                          'right fert. timing (RFT)',
                          'right fert. placement (RFP)',
                          'reduced tillage (RT)',
                          'no tillage (NT)',
                          'crop rotation (ROT)',
                          'cover cropping (CC)',
                          'residue retention (RES)'),
                        1, 'Best management measures')
ggsave('maps/Fig7_best_single.png', p, width = 25, height = 25,
       units = 'cm', dpi = 400)
cat('Fig 7 saved\n')


# 8 + 9. Figures 8 and 9 — Best two and three combined measures -------------
# Runs in both labelling modes for direct visual comparison:
#   - 'hardcoded': uses Young (2022) / dst_plotfigures.R L368-400 published
#                  top-15 lists. Output: maps/FigN_best_*_hardcoded.png
#   - 'dynamic'  : ranks combinations by frequency in current run.
#                  Output: maps/FigN_best_*_dynamic.png
# Both files are produced each run; pick the preferred one for the thesis
# after side-by-side comparison with Reference 8 and Reference 9.
# -----------------------------------------------------------------------------

for (mode in c('hardcoded', 'dynamic')) {

  # --- Figure 8 — Best two combined measures ---
  duo <- build_combo_raster(r1.p      = r1.p,
                            out.combo = out.duo,
                            sep       = ' + ',
                            mode      = mode,
                            nmax      = 2L)

  cat(sprintf('\nTop 15 duo combinations (%s mode):\n', mode))
  for (i in 1:15) cat(sprintf('  %2d  %s\n', i, duo$labels[i]))

  p <- visualize_discrete(duo$raster, 'man_num',
                          paletteer::scale_fill_paletteer_d('ggthemes::Classic_20',
                                                            direction = 1),
                          '',
                          c(seq(0.5, 15.5, by = 1), 16.5),
                          duo$labels, 2, 'Best two combined measures')
  ggsave(sprintf('maps/Fig8_best_duo_%s.png', mode), p,
         width = 25, height = 25, units = 'cm', dpi = 400)
  cat(sprintf('Fig 8 (%s) saved\n', mode))

  # --- Figure 9 — Best three combined measures ---
  trio <- build_combo_raster(r1.p      = r1.p,
                             out.combo = out.trio,
                             sep       = '+',
                             mode      = mode,
                             nmax      = 3L)

  cat(sprintf('\nTop 15 trio combinations (%s mode):\n', mode))
  for (i in 1:15) cat(sprintf('  %2d  %s\n', i, trio$labels[i]))

  p <- visualize_discrete(trio$raster, 'man_num',
                          paletteer::scale_fill_paletteer_d('ggthemes::Classic_20',
                                                            direction = 1),
                          '',
                          c(seq(0.5, 15.5, by = 1), 16.5),
                          trio$labels, 2, 'Best three combined measures')
  ggsave(sprintf('maps/Fig9_best_trio_%s.png', mode), p,
         width = 25, height = 25, units = 'cm', dpi = 400)
  cat(sprintf('Fig 9 (%s) saved\n', mode))
}

# 10. Save run metadata ------------------------------------------------------

run_meta <- list(
  script               = 'Script_3_Spatial_Visualisation.R',
  run_ts               = run_ts,
  finished_ts          = Sys.time(),
  script1_finished_ts  = script1_meta$finished_ts,
  script2a_finished_ts = script2a_meta$finished_ts,
  db_hash              = db_hash_now
)
saveRDS(run_meta, P$meta_s3)

cat('=== Script 3 complete ===\n')
cat('  Duration: ',
    sprintf('%.1f min', as.numeric(difftime(Sys.time(), run_ts, units = 'mins'))),
    '\n')


# End of Script 3 ------------------------------------------------------------
