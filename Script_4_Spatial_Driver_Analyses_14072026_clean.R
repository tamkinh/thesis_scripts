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
#    2  Load data            d1, out_best_n1, duo, trio; build the NCU-level
#                            covariate table
#    3  Covariate cross-tabs area share of each best measure by covariate level;
#                            measure availability mask (§3b)
#    4  Covariate signature  Cramér's V + multinomial logit on the five
#                            stratification covariates; ranked bar chart
#                            (Fig_cramersV)
#    5  Indicator decomp.    which indicator (yield / SOC / N surplus) gives the
#                            best measure its margin over the runner-up
#                            (selection margin); stacked bar
#                            (Fig_driver_by_measure)
#    6  Portfolios           frequency and country concentration of the top-10
#                            two- and three-measure portfolios; top-15
#                            single/duo/trio leaderboard (Fig_combo_leaderboard,
#                            report Figure 4.1; §6b)
#    7  Run metadata         provenance record saved to run_meta_script4.rds
#    8  Map                  NCU-level indicator-target driver map;
#                            conventions follow dst_plotfigures.R and
#                            Script 3 (geom_tile, coord_sf crs = 3035,
#                            scale_fill_manual, ggsave 25 x 25 cm / 400 dpi)
#                            plus a combined report panel: (a) driver bar +
#                            (b) driver map (Fig_driver_panel)
#
#  Conventions (rationale in report Section 3.5.2)
#    Covariates  cov_soil, cov_clim, cov_crop, cov_fert, cov_soc: the five axes
#                on which the meta-analytical effect sizes are stratified inside
#                cIMAm(). Referred to as stratification covariates, not site
#                covariates: crop type is a land use, and none is exogenous to
#                selection.
#    Unit        the NCU. runDST assigns one best measure per NCU.
#    Weight      area_ncu_ha_tot, matching the reference implementation, which
#                reports every share as sum(area where flag) / sum(area).
#    Covariate   per covariate, the NCU's area-dominant class.
#    Prerequisite for Sections 5 and 8: sY_combi / sSOC_combi / sNsu_combi in
#                out_best_n1.rds and out_runnerup_n1.rds, both written by
#                Script 2a at nmax = 1, uw = c(1,1,1). Absent either, Sections 5
#                and 8 self-skip and 3/4/6 run normally.
#
#  Inputs
#    results/d1.rds
#    results/out_best_n1.rds         from Script 2a
#    results/out_runnerup_n1.rds     from Script 2a (§5, §8; optional)
#    results/out_duo_n2.rds          from Script 2a
#    results/out_trio_n3.rds         from Script 2a
#    results/run_meta_script2a.rds   staleness verification
#    data/gncu2010_ext.asc           NCU raster (§8; optional)
#
#  Outputs (in results/ unless noted)
#    sxt_best_by_{soil,clim,crop,fert,soc}.csv
#    sxt_combined_summary.csv
#    measure_availability.csv
#    sxt_cramersV.csv
#    mlogit_common_coefs.csv
#    srq2_driver_by_ncu.csv
#    srq2_driver_by_measure.csv
#    srq2_driver_by_country.csv
#    top10_duos.csv, top10_trios.csv
#    portfolio_country_concentration.csv
#    maps/Fig_cramersV.png           (report Figure 4.3)
#    maps/Fig_driver_by_measure.png
#    maps/Fig_combo_leaderboard.png  (report Figure 4.1; §6b; requires patchwork)
#    maps/Fig_driver_map.png         (§8; requires terra + ggplot2)
#    maps/Fig_driver_panel.png       (report Figure 4.4; requires patchwork)
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

# Area-weighted Cramér's V between a covariate and a categorical outcome, one
# row per NCU, weighted by NCU agricultural area. V = sqrt(chi2 / (N (k-1))) is
# a normalised association measure and is invariant to rescaling the weights, so
# the raw hectare weights are passed through unchanged.
cramersV_w <- function(x, y, w) {
  d   <- data.frame(x = droplevels(as.factor(x)), y = droplevels(as.factor(y)),
                    w = as.numeric(w))
  d   <- d[is.finite(d$w) & d$w > 0, ]
  tab <- stats::xtabs(w ~ x + y, data = d)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (min(dim(tab)) < 2) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1))))
}

# Normalised inference weights for multinom(), which treats `weights` as case
# weights. Rescaling so sum(w) == nobs keeps the relative area weighting while
# returning the deviance to the observation scale (Section 3.5.2).
norm_weights <- function(area) area * (length(area) / sum(area))


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
ru <- if (file.exists(P$out_runnerup)) readRDS(P$out_runnerup) else NULL
if (is.null(ru))
  log('WARNING: ', P$out_runnerup, ' not found - SRQ2 driver decomposition will be ',
      'skipped. Re-run Script 2a to export the runner-up scores.')
ob[, man_code := factor(man_code, levels = reference_levels)]

if (length(setdiff(COV, names(d1))))
  stop('d1 missing covariates: ', paste(setdiff(COV, names(d1)), collapse = ', '))
if (!'area_ncu_ha_tot' %in% names(ob))
  stop('out_best_n1 lacks area_ncu_ha_tot; it is required as the NCU area weight.')

# NCU-level covariate table: each covariate reduced independently to the NCU's
# area-dominant class (Section 3.5.2). Climate is constant within an NCU and
# soil texture nearly so, so for those two the reduction is a no-op.
dom_class <- function(cv) {
  a <- d1[!is.na(get(cv)), .(w = sum(area_ncu, na.rm = TRUE)), by = c('ncu', cv)]
  setorderv(a, c('ncu', 'w'), c(1L, -1L))
  a[, .SD[1], by = ncu][, w := NULL][]
}
ncu_cov <- Reduce(function(x, y) merge(x, y, by = 'ncu'), lapply(COV, dom_class))

has_scores <- all(c('sY_combi','sSOC_combi','sNsu_combi') %in% names(ob))
if (!has_scores)
  log('WARN: per-indicator scores absent in out_best_n1 - sections 5 and 8 skipped. ',
      'sY_combi/sSOC_combi/sNsu_combi are exported at dst_functions_fert.R l.273 ',
      'and l.331; re-run Script 2a (nmax = 1, uw = c(1,1,1)).')

# ti_N is the framework's own measure-independent flag for "the N surplus is at
# or below its critical limit": ti_N := fifelse(dist_N <= 1, 1, 0), with
# dist_N := n_sp_ref / pmin(n_sp_sw_crit, n_sp_gw_crit) (dst_functions_fert.R
# l.102, l.292). Used in §5 to separate the two reasons an NCU can be yield- or
# SOC-driven.
keep <- c('ncu', 'man_code', 'area_ncu_ha_tot',
          if ('ti_N' %in% names(ob)) 'ti_N',
          if (has_scores) c('sY_combi','sSOC_combi','sNsu_combi'))
bwc <- merge(ob[, ..keep], ncu_cov, by = 'ncu')
if (nrow(bwc) != uniqueN(ob$ncu))
  stop('bwc is not one row per NCU (', nrow(bwc), ' rows / ',
       uniqueN(ob$ncu), ' NCUs).')
log('Merged: ', nrow(bwc), ' NCUs, ',
    round(sum(bwc$area_ncu_ha_tot) / 1e6, 2), ' Mha')


# 3. Covariate cross-tabs --------------------------------------------------
# Area share of each best measure within each covariate level, on the same
# NCU-level, area-weighted basis as §4. NCU counts are retained alongside
# the shares so that thinly populated levels are visible.
cross_tab <- function(dt, covar, outfile) {
  tabA <- xtabs(as.formula(paste('area_ncu_ha_tot ~', covar, '+ man_code')),
                data = dt)
  tabN <- table(dt[[covar]], dt$man_code)[rownames(tabA), colnames(tabA), drop = FALSE]
  tabP <- 100 * prop.table(tabA, 1)
  out  <- data.table(level = rownames(tabA))
  for (m in colnames(tabA)) {
    out[, (paste0(m, '_n'))       := as.integer(tabN[, m])]
    out[, (paste0(m, '_pct_area')) := round(tabP[, m], 1)]
  }
  out[, total_n    := rowSums(tabN)]
  out[, total_kha  := round(rowSums(tabA) / 1000, 1)]
  fwrite(out, outfile)
  invisible(out)
}

covar_files <- setNames(
  sprintf('results/sxt_best_by_%s.csv', c('soil','clim','crop','fert','soc')), COV)
for (cv in COV) cross_tab(bwc, cv, covar_files[[cv]])
log('3: covariate cross-tabs written (', length(COV), '; area shares).')

comb <- rbindlist(lapply(COV, function(cv) {
  tab <- xtabs(as.formula(paste('area_ncu_ha_tot ~', cv, '+ man_code')), data = bwc)
  rbindlist(lapply(rownames(tab), function(lv) {
    cs <- tab[lv, ]
    if (sum(cs) == 0) return(NULL)
    data.table(covariate = cv, level = lv,
               n_ncus        = bwc[get(cv) == lv, .N],
               area_kha      = round(sum(cs) / 1000, 1),
               modal_measure = names(cs)[which.max(cs)],
               modal_pct_area = round(100 * max(cs) / sum(cs), 1))
  }))
}))
setorder(comb, covariate, -area_kha)
fwrite(comb, 'results/sxt_combined_summary.csv')


# 3b. Measure availability ---------------------------------------------------
# A measure whose potential area ha_m1 is zero has its effect size zeroed inside
# cIMAm() and cannot be selected, so selection is gated by availability as well
# as by the covariates of Section 4. The mask is reported here as context for
# that association; the mechanism is described in report Section 3.5.2.
# ha_m1 is summed over an NCU's crop rows before the test, as cIMAm() does. The
# measure-to-parea mapping below reproduces cIMAm() exactly.
PAREA <- c(
  'CC'    = 'parea.cc',
  'RES'   = 'parea.cres',
  'ROT'   = 'parea.cr',
  'NT-CT' = 'parea.ntct',
  'RT-CT' = 'parea.rtct',
  'EE'    = 'parea.ee',
  'RFT'   = 'parea.rft',
  'RFP'   = 'parea.rfp',
  'CF-MF' = 'parea.nifnof + parea.hifnof + 0.75*parea.mifnof + 0.75*parea.hifhof + 0.25*parea.mifhof',
  'OF-MF' = 'parea.nifnof + parea.hifnof + 0.75*parea.mifnof + 0.75*parea.hifhof + 0.25*parea.mifhof')

avail_cols <- unique(unlist(regmatches(PAREA, gregexpr('parea[.][a-z]+', PAREA))))
if (all(avail_cols %in% names(d1))) {
  tot_area <- sum(bwc$area_ncu_ha_tot)
  avail <- rbindlist(lapply(names(PAREA), function(m) {
    ha <- d1[, .(ha_m1 = sum(eval(parse(text = PAREA[[m]])), na.rm = TRUE)), by = ncu]
    ha <- merge(ha, bwc[, .(ncu, man_code, area_ncu_ha_tot)], by = 'ncu')
    data.table(
      man_code          = m,
      n_unavailable     = ha[ha_m1 == 0, .N],
      pct_ncu_unavail   = round(100 * ha[ha_m1 == 0, .N] / nrow(ha), 1),
      pct_area_unavail  = round(100 * ha[ha_m1 == 0, sum(area_ncu_ha_tot)] / tot_area, 1),
      pct_area_selected = round(100 * ha[man_code == m, sum(area_ncu_ha_tot)] / tot_area, 1))
  }))
  setorder(avail, -pct_area_unavail)
  fwrite(avail, 'results/measure_availability.csv')
  log('3b: measure availability (ha_m1 == 0 -> effect zeroed, l.645):')
  log(paste(capture.output(print(avail)), collapse = '\n'))
} else {
  avail <- NULL
  log('3b: availability mask not computed - d1 lacks: ',
      paste(setdiff(avail_cols, names(d1)), collapse = ', '))
}


# 4. Covariate signature - area-weighted Cramér's V + multinomial + bar chart
# Report §4.3.1, Figure 4.3. One row per NCU, weighted by NCU agricultural area,
# across all ten measures. Bivariate: each V reports covariation with selection,
# not an independent contribution.
eff <- data.table(
  covariate = COV,
  label     = COV_labels,
  cramers_v = round(sapply(COV, function(v)
    cramersV_w(bwc[[v]], bwc$man_code, bwc$area_ncu_ha_tot)), 3))
setorder(eff, -cramers_v)
fwrite(eff[, .(covariate, cramers_v)], 'results/sxt_cramersV.csv')
log('4: Cramér\'s V (stratification covariates), ranked:')
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
  ggtitle("Covariate association with best measure selection") +
  theme(text         = element_text(size = 18),
        plot.title   = element_text(hjust = 0.5, size = 16),
        panel.border = element_blank(),
        axis.line.x  = element_line())
ggsave('maps/Fig_cramersV.png', p_cv,
       width = 20, height = 14, units = 'cm', dpi = 400)
log('4: Fig_cramersV.png saved')

# Estimation scope: the four measures that between them hold the great majority
# of NCUs and of agricultural area. The six minor measures are too sparse for
# stable multinomial coefficients and are excluded here; they remain in §3, §5.
common_set <- c('EE','CC','ROT','RES')
bc <- bwc[man_code %in% common_set]
bc[, man_code := relevel(factor(man_code, levels = common_set), ref = 'EE')]
for (v in COV) bc[, (v) := factor(get(v))]
bc[, w_inf := norm_weights(area_ncu_ha_tot)]
log('4: multinomial scope: ', nrow(bc), ' NCUs, ',
    round(100 * sum(bc$area_ncu_ha_tot) / sum(bwc$area_ncu_ha_tot), 1),
    '% of agricultural area')

m_common <- multinom(as.formula(paste('man_code ~', paste(COV, collapse = ' + '))),
                     data = bc, weights = w_inf, maxit = 500, trace = FALSE)
m_null   <- multinom(man_code ~ 1, data = bc, weights = w_inf,
                     maxit = 500, trace = FALSE)
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
    '  (NCUs = ', nrow(bc), ', converged = ', m_common$convergence == 0, ')')

# Fertilisation x SOC interaction, on the same estimation sample and weights.
m_int   <- multinom(
  as.formula(paste('man_code ~',
                   paste(setdiff(COV, c('cov_fert','cov_soc')), collapse = ' + '),
                   '+ cov_fert * cov_soc')),
  data = bc, weights = w_inf, maxit = 500, trace = FALSE)
mcf_int <- 1 - deviance(m_int) / deviance(m_null)

# Likelihood-ratio test for the nested comparison. nnet::multinom has no anova()
# method, so deviances are compared directly; lmtest::lrtest() gives the same.
lr_chi <- deviance(m_common) - deviance(m_int)
lr_df  <- m_int$edf - m_common$edf
lr_p   <- pchisq(lr_chi, df = lr_df, lower.tail = FALSE)
log('4: fert x SOC interaction: McFadden R2 ', round(mcf, 3), ' -> ',
    round(mcf_int, 3), ' (delta = ', round(mcf_int - mcf, 3), ')')
log('4: fert x SOC interaction: LR chi2 = ', round(lr_chi, 1),
    ', df = ', lr_df, ', p = ', signif(lr_p, 3))


# 5. Indicator-target decomposition + stacked bar chart --------------------
# Selection margin: the driver of a measure's selection is the indicator on
# which the best measure most out-scores the runner-up, i.e. the second-ranked
# single measure at the same NCU. sY_combi, sSOC_combi and sNsu_combi share one
# meaning - the fraction of the remaining gap closed, zero where there is no gap
# - and scorefun handles each indicator's direction internally, so the three are
# directly comparable and no re-orientation is applied.
# At uw = c(1,1,1) the composite is the unweighted mean of the three scores, so
# margin = 3 * (bipmc_best - bipmc_ru): this decomposes the ranking itself.
# Four classes: yield / soc / nsurplus carry the largest positive margin; 'tied'
# has no positive margin on any indicator, meaning the two leading measures
# score identically and the winner came from the seeded tie-break in runDST.
# 'tied' is a result, not a residual. Rationale and alternatives considered are
# in report Section 3.5.2.
if (has_scores && !is.null(ru)) {

  bwc <- merge(bwc, ru[, .(ncu, man_code_ru, sY_ru, sSOC_ru, sNsu_ru)],
               by = 'ncu', all.x = TRUE)
  bwc[, `:=`(m_Y   = sY_combi   - sY_ru,
             m_SOC = sSOC_combi - sSOC_ru,
             m_N   = sNsu_combi - sNsu_ru)]
  bwc[, margin := m_Y + m_SOC + m_N]          # = 3 * (bipmc_best - bipmc_ru)
  bwc[, m_max  := pmax(m_Y, m_SOC, m_N)]
  bwc[, m_min  := pmin(m_Y, m_SOC, m_N)]
  # second-largest of three values = sum - max - min
  bwc[, m_2nd  := margin - m_max - m_min]
  bwc[, driver := fcase(
    is.na(m_max) | m_max <= 0,      'tied',
    m_Y   >= m_SOC & m_Y   >= m_N,  'yield',
    m_SOC >= m_Y   & m_SOC >= m_N,  'soc',
    default =                        'nsurplus')]

  # Two confidence flags answering different questions:
  #   small_margin  was the choice of MEASURE clear?  (whole composite advantage)
  #   ambiguous     was the choice of DRIVER clear?   (largest margin vs second)
  bwc[, drv_sep      := m_max - m_2nd]
  bwc[, small_margin := driver != 'tied' & margin < 0.01]
  bwc[, ambiguous    := driver != 'tied' & drv_sep < 0.01]

  fwrite(bwc[, .(ncu, man_code, man_code_ru, area_ncu_ha_tot,
                 sY_combi, sSOC_combi, sNsu_combi,
                 m_Y, m_SOC, m_N, margin, drv_sep,
                 small_margin, ambiguous, driver)],
         P$driver_ncu)

  # Shares are area-weighted, on the same basis as §3 and §4.
  DRV_LEV <- c('yield','soc','nsurplus','tied')
  bwc[, driver := factor(driver, levels = DRV_LEV)]

  dd <- bwc[, .(area = sum(area_ncu_ha_tot), N = .N), by = .(man_code, driver)]
  dd[, pct := round(100 * area / sum(area), 1), by = man_code]
  fwrite(dcast(dd, man_code ~ driver, value.var = 'pct', fill = 0),
         'results/srq2_driver_by_measure.csv')

  ov <- bwc[, .(pct_area = round(100 * sum(area_ncu_ha_tot) /
                                   sum(bwc$area_ncu_ha_tot), 1),
                pct_ncu  = round(100 * .N / nrow(bwc), 1)), by = driver][order(-pct_area)]
  log('5: selection-margin driver split:')
  log(paste(capture.output(print(ov)), collapse = '\n'))

  n_tied <- bwc[driver == 'tied', .N]
  log('5: tied: ', n_tied, ' NCUs (',
      round(100 * n_tied / nrow(bwc), 1), '% of NCUs, ',
      round(100 * bwc[driver == 'tied', sum(area_ncu_ha_tot)] /
              sum(bwc$area_ncu_ha_tot), 1), '% of area) - the two leading ',
      'measures score identically on all three indicators and the reported ',
      'best measure is set by the tie-break at dst_functions_fert.R l.264.')

  # Confidence flags among the classified (non-tied) NCUs.
  n_cls <- bwc[driver != 'tied', .N]
  n_sm  <- bwc[small_margin == TRUE, .N]
  n_amb <- bwc[ambiguous    == TRUE, .N]
  log('5: of ', n_cls, ' classified NCUs: small margin (<0.01) ',
      round(100 * n_sm / n_cls, 1), '% - the measure choice was near-arbitrary; ',
      'ambiguous driver (drv_sep <0.01) ', round(100 * n_amb / n_cls, 1),
      '% - two indicators contribute almost equally.')

  # --- Co-occurrence of an N gap under a yield or SOC margin -----------------
  # ti_N is the framework's measure-independent flag for "the N surplus is at or
  # below its critical limit"; sNsu_combi > 0 is not, since it also requires the
  # best measure to reduce the surplus (Section 3.5.2).
  # The reciprocal statistic for the N-surplus class is forced (driver ==
  # 'nsurplus' implies ti_N == 0) and so is not reported. Reported here: of the
  # area whose margin comes from yield or SOC, the share that nonetheless has an
  # N gap - separating "no gap" from "a gap outweighed by a larger one".
  if ('ti_N' %in% names(bwc)) {
    a_ysoc    <- bwc[driver %in% c('yield','soc'), sum(area_ncu_ha_tot)]
    ngap_ysoc <- round(100 * bwc[driver %in% c('yield','soc') & ti_N == 0,
                                 sum(area_ncu_ha_tot)] / a_ysoc, 1)
    log('5: of the yield- or SOC-driven area, ', ngap_ysoc,
        '% has an N surplus above its critical limit (ti_N == 0) - the margin ',
        'came from yield or SOC despite it; ', round(100 - ngap_ysoc, 1),
        '% is already at or below the limit.')
  } else {
    ngap_ysoc <- NA_real_
    log('5: ti_N absent from out_best_n1 - N-gap co-occurrence not computed.')
  }

  # --- Driver by country ----------------------------------------------------
  # Area-weighted per-country split, including the tied class, with both
  # confidence flags as shares of the classified area.
  ncu_ctry_diag <- unique(d1[, .(ncu, country)], by = 'ncu')
  dbc <- merge(bwc, ncu_ctry_diag, by = 'ncu')[
    , .(n_ncu        = .N,
        area_kha     = round(sum(area_ncu_ha_tot) / 1000, 1),
        pct_yield    = round(100 * sum(area_ncu_ha_tot[driver == 'yield'])    / sum(area_ncu_ha_tot), 1),
        pct_soc      = round(100 * sum(area_ncu_ha_tot[driver == 'soc'])      / sum(area_ncu_ha_tot), 1),
        pct_nsurplus = round(100 * sum(area_ncu_ha_tot[driver == 'nsurplus']) / sum(area_ncu_ha_tot), 1),
        pct_tied     = round(100 * sum(area_ncu_ha_tot[driver == 'tied'])     / sum(area_ncu_ha_tot), 1),
        pct_small_margin = round(100 * sum(small_margin) / sum(driver != 'tied'), 1),
        pct_ambiguous    = round(100 * sum(ambiguous)    / sum(driver != 'tied'), 1)),
    by = country][order(-pct_nsurplus)]
  fwrite(dbc, 'results/srq2_driver_by_country.csv')

  # stacked bar: area share by driver x measure — Fig_driver_by_measure
  # All four classes are plotted, so every bar sums to 100%. Dropping 'tied'
  # would leave the minor measures as near-empty bars with no explanation.
  dd_plot <- copy(dd)
  dd_plot[, driver_label := factor(driver, levels = DRV_LEV,
                                   labels = c('Yield','SOC','N surplus','Tied'))]

  driver_cols <- c('Yield'     = '#E69F00',   # orange
                   'SOC'       = '#009E73',   # green
                   'N surplus' = '#0072B2',   # blue
                   'Tied'      = '#BBBBBB')   # grey

  p_drv <- ggplot(dd_plot,
                  aes(x = man_code, y = pct, fill = driver_label)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = driver_cols, drop = FALSE) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) +
    xlab('') + ylab('Share of agricultural area (%)') +
    labs(fill = 'Selection margin') +
    ggtitle('Indicator giving each measure its selection margin') +
    theme(text            = element_text(size = 18),
          plot.title      = element_text(hjust = 0.5, size = 14),
          legend.position = 'right',
          panel.border    = element_blank(),
          axis.line.x     = element_line())
  ggsave('maps/Fig_driver_by_measure.png', p_drv,
         width = 22, height = 18, units = 'cm', dpi = 400)
  log('5: Fig_driver_by_measure.png saved')

} else {
  log('5: driver decomposition skipped (scores or runner-up file absent).')
}


# 6. Portfolios — frequency and country concentration ----------------------
if (have_portfolios) {
  od <- readRDS(P$out_duo)
  ot <- readRDS(P$out_trio)

  if (!'area_ncu_ha_tot' %in% names(od) || !'area_ncu_ha_tot' %in% names(ot))
    stop('out_duo / out_trio lack area_ncu_ha_tot; required as the area weight.')

  # Portfolio frequency as an area share, on the same basis as §3-§5. runDST
  # filters pout4/pout5 to bipmcs == 1, so each holds one row per NCU.
  freq <- function(x) {
    f <- x[, .(area = sum(area_ncu_ha_tot), N = .N),
           by = .(portfolio_id = man_code)][order(-area)]
    f[, rank := .I]
    f[, pct  := round(100 * area / sum(area), 1)]
    f[]
  }
  duo_freq  <- freq(od)
  trio_freq <- freq(ot)
  top10_duos  <- duo_freq[rank  <= 10]
  top10_trios <- trio_freq[rank <= 10]

  ncu_country <- unique(d1[, .(ncu, country)], by = 'ncu')
  # Country concentration of each top-10 portfolio, area-weighted for
  # consistency with the leaderboard above: of the area on which a portfolio is
  # the best choice, the share falling in each country. The three largest are
  # retained per portfolio.
  conc <- function(x, top) {
    m  <- merge(x[, .(ncu, portfolio_id = man_code, area_ncu_ha_tot)],
                ncu_country, by = 'ncu')
    m  <- m[portfolio_id %in% top$portfolio_id]
    cc <- m[, .(n_ncu = .N, area_kha = round(sum(area_ncu_ha_tot) / 1000, 1)),
            by = .(portfolio_id, country)][order(portfolio_id, -area_kha)]
    cc[, share := round(100 * area_kha / sum(area_kha), 1), by = portfolio_id]
    cc[, head(.SD, 3), by = portfolio_id]
  }
  country_conc <- rbind(
    conc(od, top10_duos)[,  type := 'duo'],
    conc(ot, top10_trios)[, type := 'trio'])
  fwrite(country_conc, 'results/portfolio_country_concentration.csv')

  # No driver class is carried to the portfolios: §5's classes come from the
  # nmax = 1 optimisation, the duos and trios from separate nmax = 2 / 3 runs.
  fwrite(top10_duos,  'results/top10_duos.csv')
  fwrite(top10_trios, 'results/top10_trios.csv')
  log('6: top duo  ', top10_duos[1,  portfolio_id], ' (',
      top10_duos[1,  pct], '% of area);  top trio ', top10_trios[1, portfolio_id],
      ' (', top10_trios[1, pct], '% of area).')

  # 6b. Combination leaderboard — report Figure 4.1 --------------------------
  # Horizontal bar charts of the most frequent best-measure selections, as a
  # share of agricultural area: (a) single measures, (b) top-15 two-measure
  # combinations, (c) top-15 three-measure combinations. Combinations beyond the
  # top 15 are grouped as 'Other'.
  if (requireNamespace('patchwork', quietly = TRUE)) {
    library(patchwork)

    single_freq <- ob[, .(area = sum(area_ncu_ha_tot), N = .N),
                      by = man_code][order(-area)]
    single_freq[, pct := round(100 * area / sum(area), 1)]
    single_freq[, label := factor(as.character(man_code),
                                  levels = as.character(man_code)[order(pct)])]

    top_with_other <- function(fr, k = 15) {
      fr  <- copy(fr)[order(-area)]
      top <- fr[seq_len(min(k, nrow(fr)))]
      other_a <- fr[-seq_len(min(k, nrow(fr))), sum(area)]
      if (other_a > 0)
        top <- rbind(top,
                     data.table(portfolio_id = 'Other', area = other_a,
                                N = NA_integer_, rank = NA_integer_),
                     fill = TRUE)
      top[, pct := round(100 * area / sum(fr$area), 1)]
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
        xlab('') + ylab('Share of agricultural area (%)') + ggtitle(ttl) +
        theme(text = element_text(size = 12),
              plot.title   = element_text(size = 12, hjust = 0.5),
              panel.border = element_blank(), axis.line.x = element_line())
    }

    fig41 <- (bar_panel(single_freq, '(a) Single measure') |
                bar_panel(duo_top,     '(b) Two-measure combinations') |
                bar_panel(trio_top,    '(c) Three-measure combinations')) +
      plot_annotation(
        title = 'Most frequent best single, two, and three measure selections',
        theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))
    ggsave('maps/Fig_combo_leaderboard.png', fig41,
           width = 34, height = 16, units = 'cm', dpi = 400)
    log('6b: Fig_combo_leaderboard.png saved (report Figure 4.1; area shares). ',
        'Single-measure area shares: ',
        paste(single_freq$man_code, ' ', single_freq$pct, '%',
              sep = '', collapse = ', '))
  }

} else {
  log('6: skipped (out_duo / out_trio not available).')
}


# 7. Run metadata ----------------------------------------------------------
# The §5 quantities exist only when the driver decomposition ran, so they are
# recorded conditionally; a run without the runner-up file still writes a valid
# provenance record with those fields set to NA.
run_meta <- list(
  script            = 'Script_4_Spatial_Driver_Analyses_14072026.R',
  started_ts        = run_ts,
  finished_ts       = Sys.time(),
  r_version         = R.version.string,
  session_info      = sessionInfo(),
  input_hashes      = list(
    d1           = d1_hash_now,
    out_best     = digest::digest(file = P$out_best, algo = 'md5'),
    out_runnerup = if (file.exists(P$out_runnerup))
      digest::digest(file = P$out_runnerup, algo = 'md5') else NA_character_),
  n_ncus            = nrow(bwc),
  covariates        = COV,
  has_scores        = has_scores,
  has_runnerup      = !is.null(ru),
  unit_of_analysis   = 'NCU',
  area_weight        = 'area_ncu_ha_tot',
  covariate_reduction = 'per-covariate area-dominant class within NCU',
  multinom_mcfadden     = mcf,
  multinom_mcfadden_int = mcf_int,
  multinom_ncus         = nrow(bc),
  interaction_LRchi2 = lr_chi,
  interaction_df     = lr_df,
  interaction_p      = lr_p,
  tied_n             = if (exists('n_tied')) n_tied else NA_integer_,
  tied_pct           = if (exists('n_tied'))
    round(100 * n_tied / nrow(bwc), 1) else NA_real_,
  small_margin_pct   = if (exists('n_sm') && exists('n_cls'))
    round(100 * n_sm / n_cls, 1) else NA_real_,
  ambiguous_pct      = if (exists('n_amb') && exists('n_cls'))
    round(100 * n_amb / n_cls, 1) else NA_real_,
  ngap_under_ysoc_pct = if (exists('ngap_ysoc')) ngap_ysoc else NA_real_,
  ngap_under_ysoc_pct = if (exists('ngap_ysoc')) ngap_ysoc else NA_real_,
  upstream_meta     = list(script2a_finished = script2a_meta$finished_ts)
)
saveRDS(run_meta, P$meta_s4)
log('7: Saved: ', P$meta_s4)


# 8. Map — indicator-target driver by NCU ----------------------------------
# Conventions follow Script 3 / dst_plotfigures.R (geom_tile, coord_sf crs 3035,
# scale_fill_manual, 25 x 25 cm / 400 dpi). Colours are Okabe-Ito, chosen to be
# colourblind-safe and distinct from the Tableau_10 measure palette of Fig 7.
# Tied NCUs are mapped, not blanked: they are a result, not missing data.
if (has_scores && !is.null(ru) && requireNamespace('terra', quietly = TRUE)) {

  library(terra)

  ncu_raster_path <- if (file.exists('data/gncu2010_ext.asc')) 'data/gncu2010_ext.asc'
  else if (file.exists('gncu2010_ext.asc')) 'gncu2010_ext.asc'
  else NA_character_

  if (!is.na(ncu_raster_path)) {

    driver_ncu <- fread(P$driver_ncu)[, .(ncu, driver)]
    driver_ncu[, driver_int := fcase(
      driver == 'yield',    1L,
      driver == 'soc',      2L,
      driver == 'nsurplus', 3L,
      driver == 'tied',     4L,
      default =             NA_integer_)]

    # 2-column classify matrix: [NCU value, becomes driver_int]
    rcl      <- as.matrix(driver_ncu[!is.na(driver_int), .(ncu, driver_int)])
    r        <- rast(ncu_raster_path)
    r_driver <- classify(r, rcl, others = NA)
    names(r_driver) <- 'driver'

    df <- as.data.frame(r_driver, xy = TRUE)
    colnames(df) <- c('x', 'y', 'driver_int')
    df <- df[!is.na(df$driver_int), ]
    df$driver_label <- factor(df$driver_int, levels = c(1,2,3,4),
                              labels = c('Yield','SOC','N surplus','Tied'))

    driver_cols <- c('Yield'     = '#E69F00',
                     'SOC'       = '#009E73',
                     'N surplus' = '#0072B2',
                     'Tied'      = '#BBBBBB')

    p_map <- ggplot() +
      geom_tile(data = df, aes(x = x, y = y, fill = driver_label)) +
      coord_sf(crs = 3035, lims_method = 'box') +
      scale_fill_manual(values      = driver_cols,
                        na.translate = FALSE,
                        na.value    = 'white',
                        drop        = FALSE) +
      xlab('') + ylab('') +
      labs(fill = 'Selection margin') +
      ggtitle('Best-measure selection driver') +
      theme(text              = element_text(size = 28),
            legend.text       = element_text(size = 20),
            legend.position   = c(0.28, 0.88),
            legend.background = element_rect(fill = 'white', color = 'white'),
            panel.border      = element_blank(),
            plot.title        = element_text(hjust = 0.5))
    ggsave('maps/Fig_driver_map.png', p_map,
           width = 25, height = 25, units = 'cm', dpi = 400)
    log('8: indicator-target driver map saved to maps/Fig_driver_map.png')

  } else {
    log('8: gncu2010_ext.asc not found in data/ or working directory — map skipped.')
  }

} else if (!has_scores || is.null(ru)) {
  log('8: skipped (per-indicator scores or runner-up absent — §5 must run first).')
} else {
  log('8: terra not installed — map skipped.')
}

# ============================================================================
# Fig 4.4 — combined two-panel figure: (a) driver bar + (b) driver map
# Requires patchwork; p_drv (§5) and p_map (§8) must exist.
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

  fig44 <- (a | b) +
    plot_layout(widths = c(0.65, 1.35), guides = 'collect') +
    plot_annotation(
      title = "Indicator behind each measure's selection margin",
      theme = theme(plot.title = element_text(size = TITLE + 2, hjust = 0.5))) &
    theme(legend.position = 'bottom')

  ggsave('maps/Fig_driver_panel.png', fig44,
         width = 32, height = 18, units = 'cm', dpi = 400)
  log('Fig 4.4: Fig_driver_panel.png saved (combined a+b)')
}

# =============================================================================
log('=== Script 4 complete: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ===')
log('Duration: ', sprintf('%.1f min',
                          as.numeric(difftime(Sys.time(), run_ts, units = 'mins'))))
close(report_con)
# End of Script 4 -------------------------------------------------------------
