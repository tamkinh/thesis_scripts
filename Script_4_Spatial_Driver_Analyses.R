# =============================================================================
# SCRIPT 4 — SPATIAL DRIVERS AND INDICATOR-TARGET DECOMPOSITION
# =============================================================================
#  Purpose
#    Characterise the spatial drivers and indicator-target basis of best-
#    measure selection from the Script 2a equal-priority output. Extends the
#    analytical scope of the DSF (Young et al., 2021, 2022) with covariate
#    association analysis and an indicator-contribution decomposition.
#
#  Sections
#    0  Environment          libraries, paths, log connection
#    1  Staleness check      hash-verify d1 against Script 2a metadata
#    2  Load data            d1, out_best_n1, duo, trio; merge site covariates
#    3  Covariate cross-tabs best-measure share by site-property covariate;
#                            yield-gap descriptive table
#    4  Covariate signature  Cramér's V + multinomial logit on five site
#                            covariates; ranked bar chart (Fig_cramersV)
#    5  Indicator decomp.    which indicator target (yield / SOC / N surplus)
#                            predominantly drove selection; stacked bar
#                            (Fig_driver_by_measure)
#    6  Mediation            site covariates predicting the §5 driver class
#    7  Portfolios           frequency, country concentration, indicator
#                            drivers of the top-10 two- and three-measure
#                            portfolios; top-15 single/duo/trio combination
#                            leaderboard figure (Fig_combo_leaderboard, §7b)
#    8  Run metadata         provenance record saved to run_meta_script4.rds
#    9  Map                  NCU-level indicator-target driver map;
#                            conventions follow dst_plotfigures.R and
#                            Script 3 (geom_tile, coord_sf crs = 3035,
#                            scale_fill_manual, ggsave 25 x 25 cm / 400 dpi)
#                            plus a combined report panel: (a) driver bar +
#                            (b) driver map (Fig_driver_panel)
#
#  Site covariates (§2, §4) are the five exogenous site-property categories
#  from the INTEGRATOR dataset used by Young et al. (2022) to stratify meta-
#  analytical effect sizes. The yield-gap class is excluded from §2 and §4
#  and appears only as a descriptive distance-to-target table (§3).
#
#  PREREQUISITE for §5–§9
#  sY_combi / sSOC_combi / sNsu_combi must be present in out_best_n1.rds.
#  Add these to the column lists at dst_functions_fert.R l.268 and l.326,
#  then re-run Script 2a (nmax = 1, uw = c(1,1,1)). Until then §5/§6/§9
#  self-skip; §3/§4/§7 run normally.
#
#  Inputs
#    results/d1.rds
#    results/out_best_n1.rds         from Script 2a
#    results/out_duo_n2.rds          from Script 2a
#    results/out_trio_n3.rds         from Script 2a
#    results/run_meta_script2a.rds   staleness verification
#    data/gncu2010_ext.asc           NCU raster (§9; optional)
#
#  Outputs (in results/ unless noted)
#    sxt_best_by_{soil,clim,crop,fert,soc}.csv
#    sxt_best_by_ygap_DESCRIPTIVE.csv
#    sxt_combined_summary.csv
#    sxt_effect_sizes.csv
#    mlogit_common_coefs.csv
#    srq2_driver_by_ncu.csv
#    srq2_driver_by_measure.csv
#    srq2_mediation_cramersV.csv
#    top10_duos.csv, top10_trios.csv
#    portfolio_country_concentration.csv
#    maps/Fig_cramersV.png
#    maps/Fig_driver_by_measure.png
#    maps/Fig_combo_leaderboard.png  (§7b; requires patchwork)
#    maps/Fig_driver_map.png         (§9; requires terra + ggplot2)
#    maps/Fig_driver_panel.png       (report figure; requires patchwork)
#    run_meta_script4.rds
#    logs/script_4_report.txt
#
#  References
#    Young et al. (2021, 2022); dst_plotfigures.R, Script 3
# =============================================================================


# 0. Environment --------------------------------------------------------------
setwd('C:/Users/STMP2/Documents/Thesis/Github/Run/1')
rm(list = ls())

library(data.table)
library(digest)
library(nnet)
library(ggplot2)

theme_set(theme_bw())
source('scripts/paths_original.R')
run_ts <- Sys.time()

require_input(P$d1)
require_input(P$out_best)
require_input(P$meta_s2a)
have_portfolios <- !is.null(P$out_duo) && !is.null(P$out_trio) &&
  file.exists(P$out_duo) && file.exists(P$out_trio)

if (!dir.exists('results')) dir.create('results')
if (!dir.exists('logs'))    dir.create('logs')
if (!dir.exists('maps'))    dir.create('maps')
if (is.null(P$log_s4))  P$log_s4  <- 'logs/script_4_report.txt'
if (is.null(P$meta_s4)) P$meta_s4 <- 'results/run_meta_script4.rds'

report_con <- file(P$log_s4, open = 'wt')
log <- function(...) {
  msg <- paste0(...)
  cat(msg, '\n')
  tryCatch(
    cat(msg, '\n', file = report_con, append = TRUE),
    error = function(e) invisible(NULL)
  )
}
log('=== Script 4 run: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')


# helpers ---------------------------------------------------------------------
reference_levels <- c('EE','CC','ROT','RES','CF-MF','OF-MF','NT-CT','RT-CT','RFT','RFP')
COV   <- c('cov_soil','cov_clim','cov_crop','cov_fert','cov_soc')
COV_labels <- c('Soil texture','Climate zone','Crop type',
                'Fertilisation intensity','SOC level')

cramersV <- function(x, y) {
  if (is.factor(x)) x <- droplevels(x)
  tab <- table(x, y)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
  k   <- min(nrow(tab), ncol(tab)) - 1
  if (k <= 0) return(NA_real_)
  as.numeric(sqrt(chi / (sum(tab) * k)))
}

parse_portfolio <- function(pid_str, valid = reference_levels) {
  tk <- strsplit(pid_str, '-', fixed = TRUE)[[1]]
  out <- character(0); i <- 1
  while (i <= length(tk)) {
    two <- if (i < length(tk)) paste(tk[i:(i + 1)], collapse = '-') else ''
    if (two %in% valid) { out <- c(out, two); i <- i + 2; next }
    if (tk[i] %in% valid) out <- c(out, tk[i])
    i <- i + 1
  }
  out
}


# 1. Staleness verification ---------------------------------------------------
script2a_meta <- readRDS(P$meta_s2a)
d1_hash_now   <- digest::digest(file = P$d1, algo = 'md5')
if (!is.null(script2a_meta$input_hashes[['d1']]) &&
    d1_hash_now != script2a_meta$input_hashes[['d1']])
  log('WARN: d1 hash differs from Script 2a record. Check provenance.')
log('Script 2a finished: ',
    format(script2a_meta$finished_ts, '%Y-%m-%d %H:%M:%S'))


# 2. Load inputs and merge ----------------------------------------------------
d1 <- readRDS(P$d1)
ob <- readRDS(P$out_best)
ob[, man_code := factor(man_code, levels = reference_levels)]

if (length(setdiff(COV, names(d1))))
  stop('d1 missing covariates: ', paste(setdiff(COV, names(d1)), collapse = ', '))
ncu_cov <- d1[, .SD[1], by = ncu, .SDcols = c(COV, 'area_ncu_ha')]

has_scores <- all(c('sY_combi','sSOC_combi','sNsu_combi') %in% names(ob))
if (!has_scores)
  log('WARN: per-indicator scores absent in out_best_n1 — §3/§4/§7 skipped. ',
      'Add sY_combi/sSOC_combi/sNsu_combi to dst_functions_fert.R l.268 and l.326 ',
      'and re-run Script 2a (nmax = 1, uw = c(1,1,1)).')

keep <- c('ncu', 'man_code',
          if ('yield_gap_t' %in% names(ob)) 'yield_gap_t',
          if (has_scores) c('sY_combi','sSOC_combi','sNsu_combi'))
bwc <- merge(ob[, ..keep], ncu_cov, by = 'ncu')

if ('yield_gap_t' %in% names(bwc)) {
  bwc[, ygap_class := cut(yield_gap_t, c(-Inf,-2,-1,0,1,2,Inf),
                          labels = c('strong_below','mod_below','mild_below',
                                     'mild_above','mod_above','strong_above'))]
  # drop empty yield-gap bins (descriptive table only; some ranges hold no NCUs)
  bwc[, ygap_class := droplevels(ygap_class)]
}
log('Merged: ', nrow(bwc), ' NCUs')


# 3. Covariate cross-tabs --------------------------------------------------
cross_tab <- function(dt, covar, outfile) {
  tabN <- table(dt[[covar]], dt$man_code)
  tabP <- 100 * prop.table(tabN, 1)
  out  <- data.table(level = rownames(tabN))
  for (m in colnames(tabN)) {
    out[, (paste0(m, '_n'))   := as.integer(tabN[, m])]
    out[, (paste0(m, '_pct')) := round(tabP[, m], 1)]
  }
  out[, total_n := rowSums(tabN)]
  fwrite(out, outfile)
  invisible(out)
}

covar_files <- setNames(
  sprintf('results/sxt_best_by_%s.csv', c('soil','clim','crop','fert','soc')), COV)
for (cv in COV) cross_tab(bwc, cv, covar_files[[cv]])
log('3: site-covariate cross-tabs written (', length(COV), ').')

if ('ygap_class' %in% names(bwc)) {
  cross_tab(bwc, 'ygap_class', 'results/sxt_best_by_ygap_DESCRIPTIVE.csv')
  log('3: yield-gap descriptive table written (distance-to-target context only).')
}

comb <- rbindlist(lapply(COV, function(cv) {
  tab <- table(bwc[[cv]], bwc$man_code)
  rbindlist(lapply(rownames(tab), function(lv) {
    cs <- tab[lv, ]
    if (sum(cs) == 0) return(NULL)
    data.table(covariate = cv, level = lv, n_ncus = as.integer(sum(cs)),
               modal_measure = names(cs)[which.max(cs)],
               modal_pct     = round(100 * max(cs) / sum(cs), 1))
  }))
}))
setorder(comb, covariate, -n_ncus)
fwrite(comb, 'results/sxt_combined_summary.csv')


# 4. Covariate signature — Cramér's V + multinomial logit + bar chart ------
eff <- data.table(
  covariate = COV,
  label     = COV_labels,
  cramers_v = round(sapply(COV, function(v) cramersV(bwc[[v]], bwc$man_code)), 3))
setorder(eff, -cramers_v)
fwrite(eff[, .(covariate, cramers_v)], 'results/sxt_effect_sizes.csv')
log('4: Cramér\'s V (site covariates), ranked:')
log(paste(capture.output(print(eff[, .(label, cramers_v)])), collapse = '\n'))

# Cramér's V bar chart — Fig_cramersV
p_cv <- ggplot(eff, aes(x = reorder(label, cramers_v), y = cramers_v)) +
  geom_col(fill = '#4393c3', width = 0.65) +
  geom_text(aes(label = sprintf('%.3f', cramers_v)),
            hjust = -0.15, size = 5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(eff$cramers_v) * 1.18),
                     expand  = c(0, 0)) +
  xlab('') + ylab("Cramér's V") +
  ggtitle("Site-covariate association with best measure selection") +
  theme(text         = element_text(size = 18),
        plot.title   = element_text(hjust = 0.5, size = 16),
        panel.border = element_blank(),
        axis.line.x  = element_line())
ggsave('maps/Fig_cramersV.png', p_cv,
       width = 20, height = 14, units = 'cm', dpi = 400)
log('4: Fig_cramersV.png saved')

common_set <- c('EE','CC','ROT','RES')
bc <- bwc[man_code %in% common_set]
bc[, man_code := relevel(factor(man_code, levels = common_set), ref = 'EE')]
for (v in COV) bc[, (v) := factor(get(v))]

m_common <- multinom(as.formula(paste('man_code ~', paste(COV, collapse = ' + '))),
                     data = bc, maxit = 500, trace = FALSE)
m_null   <- multinom(man_code ~ 1, data = bc, maxit = 500, trace = FALSE)
mcf      <- 1 - deviance(m_common) / deviance(m_null)

s  <- summary(m_common)
co <- melt(as.data.table(s$coefficients,    keep.rownames = 'measure'),
           id.vars = 'measure', variable.name = 'covariate', value.name = 'beta')
se <- melt(as.data.table(s$standard.errors, keep.rownames = 'measure'),
           id.vars = 'measure', variable.name = 'covariate', value.name = 'se')
tidy <- merge(co, se, by = c('measure', 'covariate'))
tidy[, z     := beta / se]
tidy[, p     := 2 * (1 - pnorm(abs(z)))]
tidy[, abeta := abs(beta)]
setorder(tidy, measure, -abeta)
tidy[, abeta := NULL]
fwrite(tidy, 'results/mlogit_common_coefs.csv')
log('4: multinomial McFadden R2 = ', round(mcf, 3),
    '  (n = ', nrow(bc), ', converged = ', m_common$convergence == 0, ')')

m_int <- multinom(
  as.formula(paste('man_code ~',
                   paste(setdiff(COV, c('cov_fert','cov_soc')), collapse = ' + '),
                   '+ cov_fert * cov_soc')),
  data = bc, maxit = 500, trace = FALSE)
delta_aic <- AIC(m_int) - AIC(m_common)
log('4: cov_fert x cov_soc interaction: dAIC = ', round(delta_aic, 1),
    ' (negative favours interaction)')


# 5. Indicator-target decomposition + stacked bar chart --------------------
driver_modal <- NULL

if (has_scores) {
  bwc[, s_tot := sY_combi + sSOC_combi + sNsu_combi]
  bwc[, `:=`(
    c_Y   = fifelse(s_tot > 0, sY_combi   / s_tot, NA_real_),
    c_SOC = fifelse(s_tot > 0, sSOC_combi / s_tot, NA_real_),
    c_N   = fifelse(s_tot > 0, sNsu_combi / s_tot, NA_real_))]
  bwc[, driver := fcase(
    is.na(s_tot) | s_tot == 0,      'none',
    c_Y   >= c_SOC & c_Y  >= c_N,   'yield',
    c_SOC >= c_Y   & c_SOC >= c_N,  'soc',
    default =                        'nsurplus')]

  fwrite(bwc[, .(ncu, man_code, sY_combi, sSOC_combi, sNsu_combi,
                 c_Y, c_SOC, c_N, driver)],
         'results/srq2_driver_by_ncu.csv')

  dd <- bwc[, .N, by = .(man_code, driver)]
  dd[, pct := round(100 * N / sum(N), 1), by = man_code]
  fwrite(dcast(dd, man_code ~ driver, value.var = 'pct', fill = 0),
         'results/srq2_driver_by_measure.csv')

  driver_modal <- bwc[, .N, by = .(man_code, driver)][order(man_code, -N)][,
                                                                           .(driver = driver[1]), by = man_code]

  log('5: indicator-target driver split (% of all NCUs):')
  log(paste(capture.output(print(
    bwc[, .(pct = round(100 * .N / nrow(bwc), 1)), by = driver][order(-pct)]
  )), collapse = '\n'))

  n_none <- bwc[driver == 'none', .N]
  if (n_none > 0)
    log('NOTE: ', n_none, ' NCUs have all-zero indicator scores (driver = none).')

  # stacked bar: % of NCUs by driver x measure — Fig_driver_by_measure
  # dd_plot: drop the rare all-zero-score 'none' rows; the plotted factor has
  # exactly the three indicator-gap classes, in the order used by the §9 map.
  dd_plot <- dd[driver %in% c('yield','nsurplus','soc')]
  dd_plot[, driver_label := factor(driver,
                                   levels = c('yield','nsurplus','soc'),
                                   labels = c('Yield','N surplus','SOC'))]

  driver_cols <- c('Yield'     = '#E69F00',   # orange
                   'N surplus' = '#0072B2',   # blue
                   'SOC'       = '#009E73')    # green

  p_drv <- ggplot(dd_plot,
                  aes(x = man_code, y = pct, fill = driver_label)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = driver_cols) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) +
    xlab('') + ylab('Share of NCUs (%)') +
    labs(fill = 'Dominant gap') +
    ggtitle('Indicator target predominantly driving best-measure selection') +
    theme(text            = element_text(size = 18),
          plot.title      = element_text(hjust = 0.5, size = 14),
          legend.position = 'right',
          panel.border    = element_blank(),
          axis.line.x     = element_line())
  ggsave('maps/Fig_driver_by_measure.png', p_drv,
         width = 22, height = 18, units = 'cm', dpi = 400)
  log('5: Fig_driver_by_measure.png saved')

} else {
  log('5: skipped (per-indicator scores absent).')
}


# 6. Mediation — site covariates predicting indicator-target driver class --
if (has_scores) {
  dc <- bwc[driver != 'none']
  med <- data.table(
    covariate = COV,
    label     = COV_labels,
    cramers_v = round(sapply(COV, function(v) cramersV(dc[[v]], dc$driver)), 3))
  setorder(med, -cramers_v)
  fwrite(med[, .(covariate, cramers_v)], 'results/srq2_mediation_cramersV.csv')
  log('6: mediation Cramér\'s V saved')
  log(paste(capture.output(print(med[, .(label, cramers_v)])), collapse = '\n'))

  for (v in COV) dc[, (v) := factor(get(v))]
  dc[, driver := relevel(factor(driver), ref = 'nsurplus')]
  md  <- multinom(as.formula(paste('driver ~', paste(COV, collapse = ' + '))),
                  data = dc, maxit = 500, trace = FALSE)
  md0 <- multinom(driver ~ 1, data = dc, maxit = 500, trace = FALSE)
  log('6: driver ~ site covariates McFadden R2 = ',
      round(1 - deviance(md) / deviance(md0), 3))

} else {
  log('6: skipped (per-indicator scores absent).')
}


# 7. Portfolios — frequency, country concentration, indicator drivers ------
if (have_portfolios) {
  od <- readRDS(P$out_duo)
  ot <- readRDS(P$out_trio)

  freq <- function(x) {
    f <- x[, .(ncu, portfolio_id = man_code)][, .N, by = portfolio_id][order(-N)]
    f[, rank    := .I]
    f[, pct_ncu := round(100 * N / sum(N), 1)]
    f[]
  }
  duo_freq  <- freq(od)
  trio_freq <- freq(ot)
  top10_duos  <- duo_freq[rank  <= 10]
  top10_trios <- trio_freq[rank <= 10]

  ncu_country <- unique(d1[, .(ncu, country)], by = 'ncu')
  conc <- function(x, top) {
    m  <- merge(x[, .(ncu, portfolio_id = man_code)], ncu_country, by = 'ncu')
    m  <- m[portfolio_id %in% top$portfolio_id]
    cc <- m[, .N, by = .(portfolio_id, country)][order(portfolio_id, -N)]
    cc[, share := round(100 * N / sum(N), 1), by = portfolio_id]
    cc[, head(.SD, 3), by = portfolio_id]
  }
  country_conc <- rbind(
    conc(od, top10_duos)[,  type := 'duo'],
    conc(ot, top10_trios)[, type := 'trio'])
  fwrite(country_conc, 'results/portfolio_country_concentration.csv')

  if (!is.null(driver_modal)) {
    drv_code <- c(yield = 'Y', soc = 'C', nsurplus = 'N', none = '')
    port_mech <- function(pid) {
      mem <- parse_portfolio(pid)
      paste(intersect(c('Y','C','N'),
                      unique(drv_code[driver_modal[man_code %in% mem, driver]])),
            collapse = '')
    }
    top10_duos[,  mechanism := sapply(portfolio_id, port_mech)]
    top10_trios[, mechanism := sapply(portfolio_id, port_mech)]
  }
  fwrite(top10_duos,  'results/top10_duos.csv')
  fwrite(top10_trios, 'results/top10_trios.csv')
  log('7: top duo  ', top10_duos[1,  portfolio_id], ' (',
      top10_duos[1,  pct_ncu], '%);  top trio ', top10_trios[1, portfolio_id],
      ' (', top10_trios[1, pct_ncu], '%).')

  # 7b. Combination leaderboard — Fig 4.4 -----------------------------------
  # Horizontal bar charts of the most frequent best-measure selections, as a
  # share of NCUs: (a) single measures, (b) top-15 two-measure combinations,
  # (c) top-15 three-measure combinations. Combinations beyond the top 15 are
  # grouped as 'Other'. Shares use the NCU-count denominator quoted in 4.1.2.
  if (requireNamespace('patchwork', quietly = TRUE)) {
    library(patchwork)

    single_freq <- ob[, .(N = .N), by = man_code][order(-N)]
    single_freq[, pct := round(100 * N / sum(N), 1)]
    single_freq[, label := factor(as.character(man_code),
                                  levels = as.character(man_code)[order(pct)])]

    top_with_other <- function(fr, k = 15) {
      fr  <- copy(fr)[order(-N)]
      top <- fr[seq_len(min(k, nrow(fr)))]
      other_n <- fr[-seq_len(min(k, nrow(fr))), sum(N)]
      if (other_n > 0)
        top <- rbind(top,
                     data.table(portfolio_id = 'Other', N = other_n,
                                rank = NA_integer_, pct_ncu = NA_real_),
                     fill = TRUE)
      top[, pct := round(100 * N / sum(fr$N), 1)]
      lev <- top[portfolio_id != 'Other'][order(pct), portfolio_id]
      top[, label := factor(portfolio_id, levels = c('Other', lev))]  # 'Other' pinned to base
      top[]
    }
    duo_top  <- top_with_other(duo_freq)
    trio_top <- top_with_other(trio_freq)

    bar_panel <- function(dt, ttl, fillcol = '#4393c3') {
      ggplot(dt, aes(x = label, y = pct)) +
        geom_col(fill = fillcol, width = 0.72) +
        geom_text(aes(label = sprintf('%.1f', pct)), hjust = -0.15, size = 3) +
        coord_flip() +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max(dt$pct) * 1.15)) +
        xlab('') + ylab('Share of NCUs (%)') + ggtitle(ttl) +
        theme(text = element_text(size = 12),
              plot.title   = element_text(size = 12, hjust = 0.5),
              panel.border = element_blank(), axis.line.x = element_line())
    }

    fig44 <- (bar_panel(single_freq, '(a) Single measure') |
                bar_panel(duo_top,     '(b) Two-measure combinations') |
                bar_panel(trio_top,    '(c) Three-measure combinations')) +
      plot_annotation(
        title = 'Most frequent best single, two, and three measure selections',
        theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))
    ggsave('maps/Fig_combo_leaderboard.png', fig44,
           width = 34, height = 16, units = 'cm', dpi = 400)
    log('7b: Fig_combo_leaderboard.png saved (single/duo/trio frequency; Fig 4.4)')
  }

} else {
  log('7: skipped (out_duo / out_trio not available).')
}


# 8. Run metadata ----------------------------------------------------------
run_meta <- list(
  script            = 'Script_4_Analyses_23062026_clean.R',
  started_ts        = run_ts,
  finished_ts       = Sys.time(),
  r_version         = R.version.string,
  session_info      = sessionInfo(),
  input_hashes      = list(
    d1       = d1_hash_now,
    out_best = digest::digest(file = P$out_best, algo = 'md5')),
  n_ncus            = nrow(bwc),
  covariates        = COV,
  has_scores        = has_scores,
  multinom_mcfadden = mcf,
  interaction_dAIC  = delta_aic,
  upstream_meta     = list(script2a_finished = script2a_meta$finished_ts)
)
saveRDS(run_meta, P$meta_s4)
log('8: Saved: ', P$meta_s4)


# 9. Map — indicator-target driver by NCU ----------------------------------
# Conventions follow Script 3 / dst_plotfigures.R: geom_tile on the NCU
# raster, coord_sf(crs = 3035), scale_fill_manual, ggsave 25 x 25 cm / 400
# dpi. Requires terra (geom_tile reads raster via as.data.frame).
# Category colours use Okabe-Ito (perceptually equidistant, colorblind-safe):
#   Yield gap dominant  #E69F00 (orange)
#   N surplus dominant  #0072B2 (blue)
#   SOC gap dominant    #009E73 (green)
# These avoid confusion with the Tableau_10 measure colours in Fig 7.
if (has_scores && requireNamespace('terra', quietly = TRUE)) {

  library(terra)

  ncu_raster_path <- if (file.exists('data/gncu2010_ext.asc')) 'data/gncu2010_ext.asc'
  else if (file.exists('gncu2010_ext.asc')) 'gncu2010_ext.asc'
  else NA_character_

  if (!is.na(ncu_raster_path)) {

    driver_ncu <- fread('results/srq2_driver_by_ncu.csv')[, .(ncu, driver)]
    driver_ncu[, driver_int := fcase(
      driver == 'yield',    1L,
      driver == 'nsurplus', 2L,
      driver == 'soc',      3L,
      default =             NA_integer_)]

    # 2-column classify matrix: [NCU value, becomes driver_int]
    rcl      <- as.matrix(driver_ncu[!is.na(driver_int), .(ncu, driver_int)])
    r        <- rast(ncu_raster_path)
    r_driver <- classify(r, rcl, others = NA)
    names(r_driver) <- 'driver'

    df <- as.data.frame(r_driver, xy = TRUE)
    colnames(df) <- c('x', 'y', 'driver_int')
    df <- df[!is.na(df$driver_int), ]
    df$driver_label <- factor(df$driver_int, levels = c(1,2,3),
                              labels = c('Yield','N surplus','SOC'))

    driver_cols <- c('Yield'     = '#E69F00',
                     'N surplus' = '#0072B2',
                     'SOC'       = '#009E73')

    p_map <- ggplot() +
      geom_tile(data = df, aes(x = x, y = y, fill = driver_label)) +
      coord_sf(crs = 3035, lims_method = 'box') +
      scale_fill_manual(values      = driver_cols,
                        na.translate = FALSE,
                        na.value    = 'white',
                        drop        = FALSE) +
      xlab('') + ylab('') +
      labs(fill = 'Dominant gap') +
      ggtitle('Best-measure selection driver') +
      theme(text              = element_text(size = 28),
            legend.text       = element_text(size = 20),
            legend.position   = c(0.28, 0.88),
            legend.background = element_rect(fill = 'white', color = 'white'),
            panel.border      = element_blank(),
            plot.title        = element_text(hjust = 0.5))
    ggsave('maps/Fig_driver_map.png', p_map,
           width = 25, height = 25, units = 'cm', dpi = 400)
    log('9: indicator-target driver map saved to maps/Fig_driver_map.png')

  } else {
    log('9: gncu2010_ext.asc not found in data/ or working directory — map skipped.')
  }

} else if (!has_scores) {
  log('9: skipped (per-indicator scores absent — §3 must complete first).')
} else {
  log('9: terra not installed — map skipped.')
}

# ============================================================================
# Fig 4.3 — combined two-panel figure: (a) driver bar + (b) driver map
# Requires patchwork; p_drv (§5) and p_map (§9) must exist.
# Produces maps/Fig_driver_panel.png (the version for the report).
# ============================================================================
if (exists('p_drv') && exists('p_map') &&
    requireNamespace('patchwork', quietly = TRUE)) {

  library(patchwork)
  BASE <- 12; TITLE <- 13; LEG <- 9          # harmonised; legend reduced

  shared <- theme(text            = element_text(size = BASE),
                  plot.title      = element_text(size = TITLE, hjust = 0.5),
                  legend.title    = element_text(size = LEG + 1),
                  legend.text     = element_text(size = LEG),
                  legend.key.size = unit(0.40, 'cm'),
                  legend.position = 'right')   # reset map's inline c() so collect works

  a <- p_drv + ggtitle('(a) By measure') +
    shared + theme(axis.title.x = element_text(size = BASE),
                   panel.border = element_blank(), axis.line.x = element_line())

  b <- p_map + ggtitle('(b) By location') +
    shared + theme(panel.border = element_blank()) +
    guides(fill = 'none')

  fig43 <- (a | b) +
    plot_layout(widths = c(0.65, 1.35), guides = 'collect') +
    plot_annotation(
      title = 'Dominant indicator gap behind best-measure selection',
      theme = theme(plot.title = element_text(size = TITLE + 2, hjust = 0.5))) &
    theme(legend.position = 'bottom')

  ggsave('maps/Fig_driver_panel.png', fig43,
         width = 32, height = 18, units = 'cm', dpi = 400)
  log('Fig 4.3: Fig_driver_panel.png saved (combined a+b)')
}

# =============================================================================
log('=== Script 4 complete: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')
log('Duration: ', sprintf('%.1f min',
                          as.numeric(difftime(Sys.time(), run_ts, units = 'mins'))))
close(report_con)
# End of Script 4 -------------------------------------------------------------
