# =============================================================================
# SCRIPT 5 — NUTS-2 AGGREGATION AND ENVIRONMENTAL POLICY ASSESSMENT
# =============================================================================
#  Purpose
#    Translate the per-NCU best-measure result into policy-relevant spatial
#    units, benchmarked against two N-management policy frameworks.
#
#  Benchmarks
#    Absolute:  nd_ok = 1 where n_sp_new <= n_sp_crit =
#              min(n_sp_sw_crit, n_sp_gw_crit), derived from the groundwater
#              nitrate standard of 11.3 mg NO3-N L-1 (50 mg NO3 L-1) via a
#              critical leaching rate divided by a soil-specific leaching
#              fraction (De Vries et al., 2020; Velthof et al., 2009).
#              Reported as pct_area_nd_wfd (% of NCU area meeting the limit).
#
#    Relative: Farm to Fork Strategy 2030 — caveated N-surplus proxy.
#              Three scope limitations apply:
#                (a) F2F covers N and P; this analysis is N-only.
#                (b) N surplus is a pressure indicator, not a loss measure;
#                    the surplus-loss relationship is non-linear and transfer
#                    coefficients vary by an order of magnitude across regions
#                    (De Vries et al., 2020).
#                (c) The modelled N-surplus change cannot be decomposed into
#                    fertiliser-use reduction; the F2F 20%-fertiliser target
#                    is not assessable from these outputs.
#              Reported as pct_area_nsu50 to distinguish it from ND/WFD.
#
#  Sections
#     0  Environment              libraries, paths, log connection
#     1  Staleness verification   hash-verify out_best against Script 2a
#     2  Load inputs              out_best_n1.rds and d1.rds
#     3  Apportionment map        NCU to NUTS-2 fractional area weights
#     4  Join and area            merge best-measure result; effective ha
#     5  Policy indicators        nd_ok (ND/WFD) and nsu50 (F2F proxy)
#     6  Aggregation helpers      area-weighted aggregate_impacts() and
#                                 modal_measure() functions
#     7  NUTS-2 aggregation       mean outcomes per NUTS-2 region
#     8  Country aggregation      country totals + EU24 bloc row
#     9  Assessment table         nd_wfd_assessment.csv with benchmarks (Table 4.3)
#    10  Run metadata             provenance record saved to
#                                 run_meta_script5.rds
#    11  Figures                  Fig_nuts2_ecdf.png   ECDF N-surplus (Fig 4.5a)
#                                 Fig_nuts2_ecdf_nd.png ECDF ND/WFD compliance (Fig 4.5b)
#                                 Fig_country_bars.png country bar chart (Fig 4.8)
#                                 Fig_nd_vs_f2f.png    country scatter (Fig 4.7)
#                                 Fig_nuts2_nsu_red.png NUTS-2 choropleth (Fig 4.6a)
#                                 Fig_nuts2_nd_wfd.png  NUTS-2 choropleth (Fig 4.6b)
#                                 Fig_benchmark_maps.png two-map divergence panel (report) (Fig 4.6)
#                                 Fig_benchmark_ecdf.png two-ECDF panel (report) (Fig 4.5)
#                                 Spectral palette, coord_sf(crs=3035),
#                                 ggsave 25x25 cm/400 dpi (§11d-e);
#                                 requires ggplot2 + terra
#
#  Aggregation
#    Cross-boundary NCUs apportioned by
#    frac = area_int_nuts / sum(area_int_nuts) within each NCU.
#    All NUTS-2 and country statistics area-weighted by effective area.
#
#  Inputs
#    results/out_best_n1.rds         uw = c(1,1,1); simyear = 5; nmax = 1
#    results/d1.rds
#    results/run_meta_script2a.rds   staleness verification
#    data/gncu2010_ext.asc           NCU raster (§11d-e; optional)
#
#  Outputs (in results/ unless noted)
#    nuts2_aggregated.csv
#    country_aggregated.csv
#    nd_wfd_assessment.csv
#    maps/Fig_nuts2_ecdf.png
#    maps/Fig_nuts2_ecdf_nd.png
#    maps/Fig_country_bars.png       (or Fig_country_bars_a/b.png if patchwork absent)
#    maps/Fig_nd_vs_f2f.png
#    maps/Fig_nuts2_nsu_red.png      (§11d; requires terra)
#    maps/Fig_nuts2_nd_wfd.png       (§11e; requires terra)
#    maps/Fig_benchmark_maps.png     (report figure; requires patchwork)
#    maps/Fig_benchmark_ecdf.png     (report figure; requires patchwork)
#    run_meta_script5.rds
#    logs/script_5_report.txt
#
#  References
#    Young et al. (2021, 2022); dst_plotfigures.R, Script 3.
#    De Vries et al. (2020); Velthof et al. (2009).
#    European Commission (2020). Farm to Fork Strategy.
# =============================================================================


# 0. Environment --------------------------------------------------------------
setwd('C:/Users/STMP2/Documents/Thesis/Github/Run/1')
rm(list = ls())

library(data.table)
library(digest)
library(ggplot2)

theme_set(theme_bw())
source('scripts/paths_original.R')
run_ts <- Sys.time()

P5 <- list(
  nuts2   = 'results/nuts2_aggregated.csv',
  country = 'results/country_aggregated.csv',
  nd_wfd  = 'results/nd_wfd_assessment.csv',
  meta    = 'results/run_meta_script5.rds',
  log     = 'logs/script_5_report.txt'
)
NSU50_THRESHOLD <- 0.50

require_input(P$out_best)
require_input(P$d1)
require_input(P$meta_s2a)

if (!dir.exists('results')) dir.create('results')
if (!dir.exists('logs'))    dir.create('logs')
if (!dir.exists('maps'))    dir.create('maps')

report_con <- file(P5$log, open = 'wt')
log <- function(...) {
  msg <- paste0(...)
  cat(msg, '\n')
  tryCatch(
    cat(msg, '\n', file = report_con, append = TRUE),
    error = function(e) invisible(NULL)
  )
}
log('=== Script 5 run: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')


# 1. Staleness verification ------------------------------------------------
script2a_meta <- readRDS(P$meta_s2a)
out_best_hash <- digest::digest(file = P$out_best, algo = 'md5')
if (!is.null(script2a_meta$input_hashes[['out_best']]) &&
    out_best_hash != script2a_meta$input_hashes[['out_best']])
  log('WARN: out_best_n1 hash differs from Script 2a record. Verify baseline.')
log('Script 2a finished: ',
    format(script2a_meta$finished_ts, '%Y-%m-%d %H:%M:%S'))


# 2. Load inputs -----------------------------------------------------------
best <- as.data.table(readRDS(P$out_best))
d1   <- as.data.table(readRDS(P$d1))
log('NCUs in best-measure result: ', uniqueN(best$ncu))


# 3. NCU to NUTS-2 apportionment map ---------------------------------------
nuts_map <- unique(d1[, .(ncu, NUTS2, country, area_int_nuts)])
nuts_map  <- nuts_map[!is.na(area_int_nuts) & area_int_nuts > 0]
nuts_map[, frac := area_int_nuts / sum(area_int_nuts), by = ncu]

n_multi <- nuts_map[, .N, by = ncu][N > 1, .N]
log('NCUs spanning >1 NUTS-2: ', n_multi,
    ' (', round(100 * n_multi / uniqueN(nuts_map$ncu), 1), '%)')
log('NUTS-2 regions: ', uniqueN(nuts_map$NUTS2),
    ' | countries: ',   uniqueN(nuts_map$country))


# 4. Join and compute effective area ---------------------------------------
keep <- intersect(
  c('ncu', 'area_ncu_ha_tot', 'man_code',
    'dY', 'dSOC', 'dNsu', 'tm_Y', 'tm_C', 'tm_N',
    'n_sp_ref', 'n_sp_new', 'n_sp_crit', 'n_fert'),
  names(best))
bn <- merge(best[, ..keep], nuts_map, by = 'ncu', allow.cartesian = TRUE)
bn[, aw := area_ncu_ha_tot * frac]

area_in  <- best[, sum(area_ncu_ha_tot, na.rm = TRUE)]
area_out <- bn[,   sum(aw,              na.rm = TRUE)]
log('Area check: total NCU ha = ', round(area_in),
    ' | apportioned ha = ', round(area_out),
    ' | ratio = ', round(area_out / area_in, 4))

na_crit <- bn[is.na(n_sp_crit), sum(aw, na.rm = TRUE)]
log('Area with undefined critical limit (n_sp_crit is NA): ',
    round(na_crit / 1000, 1), ' kha (',
    round(100 * na_crit / area_out, 2), '% of total). These NCUs have no ND/WFD ',
    'limit defined and are excluded from the absolute benchmark by na.rm.')


# 5. Policy indicators per NCU-NUTS2 piece ---------------------------------
# Absolute (ND/WFD): an NCU complies where its post-measure N surplus lies at or
# below its site-specific critical limit, as specified in Section 3.5.3 of the
# report. The comparison is made directly on the levels, which is also the form
# used by Young et al. (2022): their weighting function sN is built on
# (1 + dN) * N_ref / N_critical, i.e. on n_sp_new / n_sp_crit, and reaches its
# floor exactly where n_sp_new <= n_sp_crit.
#
# The framework's own flag tm_N is NOT used here. tm_N is dNsu <= 1 - dist_N (the
# tm_* block of runDST in dst_functions_fert.R). That form is exact only if dNsu
# were scaled by n_sp_crit, whereas dNsu is a meta-analytical fractional change of
# the reference (dNsu := mmean_Nsu * 0.01; n_sp_new := n_sp_ref * (1 + dNsu)), so
# the two conditions coincide only at dist_N == 1. tm_N is retained for the
# reproduction of the reference target-met shares (Script 2b, Table 4.1), where
# matching the framework is the objective; it is not the right instrument for
# benchmarking an external policy limit.
#
# n_sp_crit is NA where no critical limit is defined; nd_ok is left NA there and
# those NCUs drop out of the benchmark via na.rm, their area reported above.
bn[, nd_ok := as.integer(n_sp_new <= n_sp_crit)]

# Relative: nsu50 = 1 where N-surplus reduction >= NSU50_THRESHOLD
bn[, nsu_red := fifelse(n_sp_ref > 0, (n_sp_ref - n_sp_new) / n_sp_ref, NA_real_)]
bn[, nsu50   := fifelse(n_sp_ref > 0,
                        as.integer(nsu_red >= NSU50_THRESHOLD),
                        as.integer(n_sp_new <= 0))]


# 6. Aggregation helpers ---------------------------------------------------
aggregate_impacts <- function(dt, by_cols) {
  dt[, .(
    area_kha        = sum(aw, na.rm = TRUE) / 1000,
    dY_pct          = weighted.mean(dY,   aw, na.rm = TRUE) * 100,
    dSOC_pct        = weighted.mean(dSOC, aw, na.rm = TRUE) * 100,
    base_surplus    = weighted.mean(n_sp_ref, aw, na.rm = TRUE),   # <-- ADD: area-weighted baseline N surplus (kg N/ha)
    crit_limit      = weighted.mean(n_sp_crit, aw, na.rm = TRUE),  # <-- ADD: area-weighted critical limit (kg N/ha)
    nsu_red_kgha    = weighted.mean(n_sp_ref - n_sp_new, aw, na.rm = TRUE),
    nsu_red_pct     = weighted.mean(nsu_red, aw, na.rm = TRUE) * 100,
    pct_area_tmY    = weighted.mean(tm_Y,  aw, na.rm = TRUE) * 100,
    pct_area_tmC    = weighted.mean(tm_C,  aw, na.rm = TRUE) * 100,
    pct_area_nd_wfd = weighted.mean(nd_ok, aw, na.rm = TRUE) * 100,
    # Diagnostic only, not a benchmark: the same area share under the DSF's own
    # target-met flag. Carried so that the divergence between the framework's
    # internal flag and the policy condition is readable from the output rather
    # than asserted (report Table 3.4). Table 4.1 reports on this basis; the
    # ND/WFD benchmark of Section 4.4 reports on pct_area_nd_wfd above.
    pct_area_tmN_dsf = weighted.mean(tm_N, aw, na.rm = TRUE) * 100,
    pct_area_nsu50  = weighted.mean(nsu50, aw, na.rm = TRUE) * 100
  ), by = by_cols]
}

modal_measure <- function(dt, by_cols) {
  rk <- dt[, .(a = sum(aw, na.rm = TRUE)), by = c(by_cols, 'man_code')]
  rk[order(-a), .(modal_measure = man_code[1]), by = by_cols]
}


# 7. NUTS-2 aggregation ----------------------------------------------------
nuts2_agg <- aggregate_impacts(bn, c('NUTS2', 'country'))
nuts2_agg <- merge(nuts2_agg, modal_measure(bn, 'NUTS2'), by = 'NUTS2')
setorder(nuts2_agg, country, NUTS2)
num_cols <- setdiff(names(nuts2_agg), c('NUTS2', 'country', 'modal_measure'))
nuts2_agg[, (num_cols) := lapply(.SD, round, 2), .SDcols = num_cols]
fwrite(nuts2_agg, P5$nuts2)
log('Saved: ', P5$nuts2, '  (', nrow(nuts2_agg), ' NUTS-2 regions)')


# 8. Country aggregation + EU24 bloc ---------------------------------------
country_agg <- aggregate_impacts(bn, 'country')
country_agg <- merge(country_agg, modal_measure(bn, 'country'), by = 'country')

bn[, bloc := fifelse(country == 'UK', 'UK', 'EU24')]
bloc_agg <- aggregate_impacts(bn, 'bloc')
bloc_agg <- merge(bloc_agg, modal_measure(bn, 'bloc'), by = 'bloc')
setnames(bloc_agg, 'bloc', 'country')
bloc_agg <- bloc_agg[country == 'EU24']   # UK row already in country_agg

country_out <- rbind(country_agg, bloc_agg, fill = TRUE)
setorder(country_out, country)
num_cols <- setdiff(names(country_out), c('country', 'modal_measure'))
country_out[, (num_cols) := lapply(.SD, round, 2), .SDcols = num_cols]
fwrite(country_out, P5$country)
log('Saved: ', P5$country,
    '  (', nrow(country_out), ' rows: ',
    uniqueN(country_agg$country), ' countries + EU24 bloc)')


# 9. ND/WFD + F2F assessment table -----------------------------------------
nd_cols <- c('country', 'area_kha', 'nsu_red_kgha', 'nsu_red_pct',
             'pct_area_nd_wfd', 'pct_area_nsu50', 'pct_area_tmN_dsf')
nd_tab <- country_out[, ..nd_cols]
fwrite(nd_tab, P5$nd_wfd)
log('Saved: ', P5$nd_wfd)

eu24 <- nd_tab[country == 'EU24']
uk   <- nd_tab[country == 'UK']
log('EU24  N-surplus red = ', round(eu24$nsu_red_pct, 1),
    '%  | ND/WFD area = ', round(eu24$pct_area_nd_wfd, 1),
    '%  | F2F-N proxy = ', round(eu24$pct_area_nsu50,  1), '%')
log('UK    N-surplus red = ', round(uk$nsu_red_pct,    1),
    '%  | ND/WFD area = ', round(uk$pct_area_nd_wfd,  1),
    '%  | F2F-N proxy = ', round(uk$pct_area_nsu50,   1), '%')

# Divergence between the ND/WFD policy condition (n_sp_new <= n_sp_crit, Section
# 3.5.3) and the DSF's internal target-met flag tm_N, which Table 4.1 reports.
# Both are area shares on the same apportioned basis. Recorded for Table 3.4.
log('ND/WFD (n_sp_new <= n_sp_crit) vs DSF flag tm_N -- EU24: ',
    round(eu24$pct_area_nd_wfd, 1), '% vs ', round(eu24$pct_area_tmN_dsf, 1),
    '%  | UK: ', round(uk$pct_area_nd_wfd, 1), '% vs ',
    round(uk$pct_area_tmN_dsf, 1), '%')


# 10. Run metadata ---------------------------------------------------------
run_meta <- list(
  script              = 'Script_5_nuts_2_N_assessment.R',
  started_ts          = run_ts,
  finished_ts         = Sys.time(),
  r_version           = R.version.string,
  session_info        = sessionInfo(),
  scope               = '24 EU member states + UK (n = 25)',
  nsu50_threshold     = NSU50_THRESHOLD,
  benchmark_absolute  = paste('ND/WFD groundwater-nitrate limit (n_sp_crit; 11.3 mg NO3-N/L;',
                              'site-specific). Compliance = n_sp_new <= n_sp_crit, evaluated on',
                              'levels per report Section 3.5.3; the DSF flag tm_N is not used here'),
  benchmark_relative  = 'F2F N-surplus proxy >=50% reduction (N-only; caveated; uniform)',
  n_nuts2             = uniqueN(nuts2_agg$NUTS2),
  n_countries         = uniqueN(country_agg$country),
  area_ratio          = area_out / area_in,
  na_crit_kha         = na_crit / 1000,
  input_hashes        = list(
    out_best = digest::digest(file = P$out_best, algo = 'md5'),
    d1       = digest::digest(file = P$d1,       algo = 'md5')),
  upstream_meta       = list(script2a_finished = script2a_meta$finished_ts)
)
saveRDS(run_meta, P5$meta)
log('Saved: ', P5$meta)


# 11. Figures — all Chapter 4 SRQ3 figures ---------------------------------
# Colour conventions:
#   EU-24 = '#1a6faf' (blue), UK = '#e08214' (orange)
#   ND/WFD = '#1a9641' (green), F2F proxy = '#d73027' (red)
# ggplot2 required for §11a-c; terra additionally required for §11d-e.

C_EU  <- '#1a6faf'; C_UK <- '#e08214'
C_ND  <- '#1a9641'; C_F2F <- '#d73027'

eu_nuts <- nuts2_agg[country != 'UK']
uk_nuts <- nuts2_agg[country == 'UK']
eu_cnt  <- country_out[country != 'UK' & country != 'EU24']

# §11a — ECDF of NUTS-2 N-surplus reduction and ND/WFD compliance
# Two complementary benchmarks, one ECDF each: the relative F2F >=50%
# N-surplus-reduction proxy, and the absolute ND/WFD compliance limit.

# §11a-i: N-surplus reduction ECDF (Relative F2F 50% proxy)
ecdf_eu <- data.frame(
  x = sort(eu_nuts$nsu_red_pct),
  y = seq_along(eu_nuts$nsu_red_pct) / nrow(eu_nuts) * 100)
ecdf_uk <- data.frame(
  x = sort(uk_nuts$nsu_red_pct),
  y = seq_along(uk_nuts$nsu_red_pct) / nrow(uk_nuts) * 100)

eu_above50 <- round(mean(eu_nuts$nsu_red_pct >= 50) * 100)
uk_above50 <- round(mean(uk_nuts$nsu_red_pct >= 50) * 100)

p_ecdf_nsu <- ggplot() +
  geom_line(data = ecdf_eu, aes(x = x, y = y),
            colour = C_EU, linewidth = 1.4) +
  geom_line(data = ecdf_uk, aes(x = x, y = y),
            colour = C_UK, linewidth = 1.4, linetype = 'dashed') +
  geom_vline(xintercept = 50, colour = C_F2F,
             linewidth = 1, linetype = 'dotted') +
  annotate('text', x = 50.5, y = 10, hjust = 0, size = 5,
           colour = C_F2F, label = 'F2F 50% reference') +
  annotate('text', x = 60, y = 100 - eu_above50 + 4, hjust = 1, size = 5,
           colour = C_EU,
           label = paste0('EU-24: ', eu_above50, '% of regions \u2265 50%')) +
  annotate('text', x = 55, y = 100 - uk_above50 - 6, hjust = 1, size = 5,
           colour = C_UK,
           label = paste0('UK: ', uk_above50, '% of regions \u2265 50%')) +
  scale_x_continuous(limits = c(15, 65), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  xlab('Mean N-surplus reduction (%)') +
  ylab('Cumulative share of NUTS-2 regions (%)') +
  ggtitle('N-surplus reduction by NUTS-2 region') +
  theme(text         = element_text(size = 18),
        plot.title   = element_text(hjust = 0.5),
        panel.grid   = element_line(colour = 'grey90'),
        panel.border = element_rect(colour = 'grey60'))
ggsave('maps/Fig_nuts2_ecdf.png', p_ecdf_nsu,
       width = 18, height = 14, units = 'cm', dpi = 400)
log('11a-i: Fig_nuts2_ecdf.png saved')

# §11a-ii: ND/WFD compliance ECDF (absolute limit)
ecdf_nd_eu <- data.frame(
  x = sort(eu_nuts$pct_area_nd_wfd),
  y = seq_along(eu_nuts$pct_area_nd_wfd) / nrow(eu_nuts) * 100)
ecdf_nd_uk <- data.frame(
  x = sort(uk_nuts$pct_area_nd_wfd),
  y = seq_along(uk_nuts$pct_area_nd_wfd) / nrow(uk_nuts) * 100)

eu_nd_above50 <- round(mean(eu_nuts$pct_area_nd_wfd >= 50) * 100)
uk_nd_above50 <- round(mean(uk_nuts$pct_area_nd_wfd >= 50) * 100)

p_ecdf_nd <- ggplot() +
  geom_line(data = ecdf_nd_eu, aes(x = x, y = y),
            colour = C_EU, linewidth = 1.4) +
  geom_line(data = ecdf_nd_uk, aes(x = x, y = y),
            colour = C_UK, linewidth = 1.4, linetype = 'dashed') +
  annotate('text', x = max(ecdf_nd_eu$x) * 0.95,
           y = 100 - eu_nd_above50 + 4, hjust = 1, size = 5,
           colour = C_EU,
           label = paste0('EU-24: ', eu_nd_above50,
                          '% of regions \u2265 50% compliance')) +
  annotate('text', x = max(ecdf_nd_uk$x) * 0.95,
           y = 100 - uk_nd_above50 - 6, hjust = 1, size = 5,
           colour = C_UK,
           label = paste0('UK: ', uk_nd_above50,
                          '% of regions \u2265 50% compliance')) +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  xlab('Area meeting ND/WFD groundwater-nitrate limit (%)') +
  ylab('Cumulative share of NUTS-2 regions (%)') +
  ggtitle('ND/WFD groundwater-nitrate compliance by NUTS-2 region') +
  theme(text         = element_text(size = 18),
        plot.title   = element_text(hjust = 0.5, size = 15),
        panel.grid   = element_line(colour = 'grey90'),
        panel.border = element_rect(colour = 'grey60'))
ggsave('maps/Fig_nuts2_ecdf_nd.png', p_ecdf_nd,
       width = 18, height = 14, units = 'cm', dpi = 400)
log('11a-ii: Fig_nuts2_ecdf_nd.png saved')

# §11b — Country bar chart (N-surplus reduction left; ND/WFD compliance right)
cnt_plot <- eu_cnt[order(nsu_red_pct)]
cnt_plot[, country := factor(country, levels = country)]

p_bar1 <- ggplot(cnt_plot, aes(x = country, y = nsu_red_pct)) +
  geom_col(fill = C_EU, width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 70), expand = c(0, 0)) +
  xlab('') + ylab('Mean N-surplus reduction (%)') +
  ggtitle('(a) N-surplus reduction') +
  theme(text = element_text(size = 14), plot.title = element_text(size = 13),
        panel.border = element_blank(), axis.line.x = element_line())

p_bar2 <- ggplot(cnt_plot, aes(x = country, y = pct_area_nd_wfd)) +
  geom_col(fill = C_ND, width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  xlab('') + ylab('Area meeting ND/WFD limit (%)') +
  ggtitle('(b) ND/WFD compliance') +
  theme(text = element_text(size = 14), plot.title = element_text(size = 13),
        axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
        panel.border = element_blank(), axis.line.x  = element_line())

# side-by-side using cowplot if available, else patchwork, else separate files
if (requireNamespace('patchwork', quietly = TRUE)) {
  library(patchwork)
  p_bars <- p_bar1 + p_bar2 +
    plot_annotation(title = 'N-surplus reduction and ND/WFD compliance by country',
                    theme = theme(plot.title = element_text(size = 15, hjust = 0.5)))
  ggsave('maps/Fig_country_bars.png', p_bars,
         width = 28, height = 20, units = 'cm', dpi = 400)
  log('11b: Fig_country_bars.png saved (patchwork)')
} else {
  ggsave('maps/Fig_country_bars_a.png', p_bar1,
         width = 16, height = 20, units = 'cm', dpi = 400)
  ggsave('maps/Fig_country_bars_b.png', p_bar2,
         width = 14, height = 20, units = 'cm', dpi = 400)
  log('11b: Fig_country_bars_a/b.png saved (patchwork not installed; panels split)')
  log('     Run install.packages("patchwork") for a combined figure.')
}

# §11c — Country scatter: F2F proxy (x) vs ND/WFD compliance (y)
cnt_s <- country_out[!country %in% c('EU24','UK')]
uk_s  <- country_out[country == 'UK']

p_scat <- ggplot() +
  geom_abline(slope = 1, intercept = 0,
              colour = 'grey60', linetype = 'dashed', linewidth = 0.8) +
  geom_point(data = cnt_s,
             aes(x = pct_area_nsu50, y = pct_area_nd_wfd,
                 size = area_kha),
             colour = C_EU, alpha = 0.75, stroke = 0.3) +
  geom_point(data = uk_s,
             aes(x = pct_area_nsu50, y = pct_area_nd_wfd,
                 size = area_kha),
             colour = C_UK, shape = 18, alpha = 0.9) +
  geom_text(data = rbind(cnt_s[, .(country, pct_area_nsu50, pct_area_nd_wfd)],
                         uk_s[,  .(country, pct_area_nsu50, pct_area_nd_wfd)]),
            aes(x = pct_area_nsu50 + 1.5, y = pct_area_nd_wfd + 1.5,
                label = country),
            size = 3.5, colour = 'grey30') +
  geom_vline(xintercept = 50, colour = C_F2F,
             linewidth = 0.8, linetype = 'dotted') +
  geom_hline(yintercept = 50, colour = C_ND,
             linewidth = 0.8, linetype = 'dotted') +
  scale_size_continuous(range = c(2, 10), guide = 'none') +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  xlab('Area achieving \u2265 50% N-surplus reduction — F2F N proxy (%)') +
  ylab('Area meeting ND/WFD groundwater-nitrate compliance (%)') +
  ggtitle('Relative vs absolute N benchmarks by country') +
  theme(text         = element_text(size = 15),
        plot.title   = element_text(hjust = 0.5),
        panel.border = element_rect(colour = 'grey60'))
ggsave('maps/Fig_nd_vs_f2f.png', p_scat,
       width = 20, height = 18, units = 'cm', dpi = 400)
log('11c: Fig_nd_vs_f2f.png saved')

# §11d-e — NUTS-2 choropleths via NCU raster (terra required)
# Palette: Spectral direction = -1 (red = low/bad, green = high/good)
# Same convention as Figs 4-6 in Script 3 / dst_plotfigures.R.
# N-surplus reduction is a mean % change per region (not % of area).
# ND/WFD compliance is % of NCU area meeting the limit.
if (requireNamespace('terra', quietly = TRUE)) {

  library(terra)

  ncu_raster_path <- if (file.exists('data/gncu2010_ext.asc')) 'data/gncu2010_ext.asc'
  else if (file.exists('gncu2010_ext.asc')) 'gncu2010_ext.asc'
  else NA_character_

  if (!is.na(ncu_raster_path)) {

    r <- rast(ncu_raster_path)

    # assign each NCU the value of its dominant NUTS-2 (largest area_int_nuts)
    ncu_nuts2 <- nuts_map[, .(area_int_nuts = sum(area_int_nuts)), by = .(ncu, NUTS2)]
    ncu_nuts2 <- ncu_nuts2[, .SD[which.max(area_int_nuts)], by = ncu]
    ncu_nuts2 <- merge(ncu_nuts2[, .(ncu, NUTS2)],
                       nuts2_agg[, .(NUTS2, nsu_red_pct, pct_area_nd_wfd)],
                       by = 'NUTS2', all.x = TRUE)

    make_ncu_layer <- function(ncu_dt, value_col) {
      mat <- as.matrix(ncu_dt[!is.na(get(value_col)), .(ncu, get(value_col))])
      classify(r, mat, others = NA)
    }

    visualize_discrete2_raster <- function(rst, layer_name, mapcolor, direction,
                                           unit_label, breaks, labels, ftitle) {
      names(rst) <- layer_name
      df <- as.data.frame(rst, xy = TRUE)
      colnames(df) <- c('x', 'y', 'variable')
      ggplot() +
        geom_tile(data = df,
                  aes(x = x, y = y,
                      fill = cut(variable, breaks, labels = labels,
                                 include.lowest = TRUE))) +
        coord_sf(crs = 3035, lims_method = 'box') +
        scale_fill_brewer(palette      = mapcolor,
                          direction    = direction,
                          na.translate = FALSE,
                          na.value     = 'white',
                          drop         = FALSE) +
        xlab('') + ylab('') + labs(fill = unit_label) +
        theme(text              = element_text(size = 28),
              legend.text       = element_text(size = 20),
              legend.position      = c(0.28, 0.88),
              legend.background    = element_rect(fill = 'white', color = 'white'),
              panel.border         = element_blank(),
              plot.title           = element_text(hjust = 0.5)) +
        ggtitle(ftitle)
    }

    # §11d  N-surplus reduction (mean % per NUTS-2 region)
    r_nsu <- make_ncu_layer(ncu_nuts2, 'nsu_red_pct')
    p_nsu <- visualize_discrete2_raster(
      r_nsu, 'nsu_red_pct', 'Spectral', 1,
      'Mean reduction (%)',
      breaks = c(0, 30, 40, 45, 50, 55, 60),
      labels = c('< 30','30 - 40','40 - 45','45 - 50','50 - 55','55 - 60'),
      ftitle = 'N-surplus reduction by NUTS-2 region')
    ggsave('maps/Fig_nuts2_nsu_red.png', p_nsu,
           width = 25, height = 25, units = 'cm', dpi = 400)
    log('11d: Fig_nuts2_nsu_red.png saved')

    # §11e  ND/WFD compliance (% of NCU area meeting limit)
    r_nd <- make_ncu_layer(ncu_nuts2, 'pct_area_nd_wfd')
    p_nd <- visualize_discrete2_raster(
      r_nd, 'pct_area_nd_wfd', 'Spectral', 1,
      'Area meeting limit (%)',
      breaks = c(0, 10, 20, 30, 50, 70, 90, 100),
      labels = c('< 10','10 - 20','20 - 30','30 - 50',
                 '50 - 70','70 - 90','> 90'),
      ftitle = 'ND/WFD compliance by NUTS-2 region')
    ggsave('maps/Fig_nuts2_nd_wfd.png', p_nd,
           width = 30, height = 30, units = 'cm', dpi = 400)
    log('11e: Fig_nuts2_nd_wfd.png saved')

  } else {
    log('11d-e: gncu2010_ext.asc not found — NUTS-2 choropleth maps skipped.')
  }

} else {
  log('11d-e: terra not installed — NUTS-2 choropleth maps skipped.')
}

# =======================================================================================
# Fig — geographic divergence: two NUTS-2 maps (separate legends; different quantities)
# Requires patchwork; p_nsu (§11d) and p_nd (§11e) must exist.
if (exists('p_nsu') && exists('p_nd') && requireNamespace('patchwork', quietly = TRUE)) {
  library(patchwork)
  mtheme <- theme(text            = element_text(size = 12),
                  plot.title      = element_text(size = 13, hjust = 0.5),
                  legend.title    = element_text(size = 9),
                  legend.text     = element_text(size = 8),
                  legend.key.size = unit(0.35, 'cm'),
                  legend.position = c(0.24, 0.85))
  a <- p_nsu + ggtitle('(a) Relative: N-surplus reduction (%)') + mtheme
  b <- p_nd  + ggtitle('(b) Absolute: area meeting ND/WFD limit (%)') + mtheme
  fig_maps <- (a | b) +
    plot_annotation(
      title = 'Geographic divergence of the relative and absolute nitrogen benchmarks',
      theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))
  ggsave('maps/Fig_benchmark_maps.png', fig_maps,
         width = 34, height = 18, units = 'cm', dpi = 400)
  log('Fig: Fig_benchmark_maps.png saved (two-map divergence panel)')
}

# Fig — threshold attainment: two ECDFs (separate x-axes; different quantities)
# Requires patchwork; p_ecdf_nsu (§11a-i) and p_ecdf_nd (§11a-ii) must exist.
if (exists('p_ecdf_nsu') && exists('p_ecdf_nd') && requireNamespace('patchwork', quietly = TRUE)) {
  library(patchwork)
  etheme <- theme(text       = element_text(size = 12),
                  plot.title = element_text(size = 13, hjust = 0.5))
  a <- p_ecdf_nsu + ggtitle('(a) Relative benchmark (F2F 50% proxy)') + etheme
  b <- p_ecdf_nd  + ggtitle('(b) Absolute benchmark (ND/WFD limit)')  + etheme
  fig_ecdf <- (a | b) +
    plot_annotation(
      title = 'Cumulative attainment of the two nitrogen benchmarks across NUTS-2 regions',
      theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))
  ggsave('maps/Fig_benchmark_ecdf.png', fig_ecdf,
         width = 30, height = 14, units = 'cm', dpi = 400)
  log('Fig: Fig_benchmark_ecdf.png saved (two-ECDF attainment panel)')
}

# =============================================================================
log('=== Script 5 complete: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')
log('Duration: ', sprintf('%.1f min',
                          as.numeric(difftime(Sys.time(), run_ts, units = 'mins'))))
close(report_con)
# End of Script 5 -------------------------------------------------------------
