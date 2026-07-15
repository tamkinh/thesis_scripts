# Functions to run DST of Maddy
#**** this is automatically sourced by dst_prepare_input_fert

#Parameters to define as input when calling runDST
#data tables are imported in dst_main

#' @param db (data.table) a data.table of INTEGRATOR/EUROSTAT data
#' @param dt.m (data.table) a data.table of meta-analytical estimates for Y, Nsu and SOC per management and ncu
#' @param uw (vector) (user weight) is an argument that gives the relative importance for yield, soc or n-surplus (a value above 0, likely a fraction between 0 and 2). Default is c(1,1,1).
#' @param simyear (num) value for the default year for simulation of impacts (default is 5 years)
#' @param output (string) optional arguments to select different types of outputs. Options include: 'total_impact','best_impact','score_single','score_due','score_best' or 'all'.
#' @param quiet (boolean) option to show progress bar to see progress of the function
#' @param nmax (integer) the max number of measure combinations evaluated
#' @param nopt (boolean) option to prefer lower number of measures above higher number when equal in score
#'
#' @details
#'  run the optimizer with the following options for output:
#'
#'  'all' gives the five following outputs:
#'
#'  'total_impact' lists all impacts for all combinations of measures
#'                 returns various numbers of options per NCU with changes listed in all indicators
#'                 maps of impacts for each and/or all measure options that can be applied
#'
#'
#'  'best_impact' gives the one best option per NCU for total change in SOC, Nsp and yield
#'                what are the changes in indicators when only the best measure or best combo of measures applied
#'
#'  'score_single' gives the order of all single measures by how they contribute to reaching desired targets for indicators
#'                 returns one ncu row with 6 columns ranking all measures individually
#'
#'  'score_duo' gives the order of all combinations of 2 measures by how they contribute to reaching desired targets
#'              returns one ncu row with 16 columns ranking each combo of measures
#'
#'  'score_best' gives the best combination to reach desired targets
#'
#'  the uw input argument is always in the order of (1) yield, (2) soc and (3) N surplus.
#'
#'  @export ??


# # # copy here inputs to run line-by-line, then run each line within the function rather than calling it
# db = d1[ncu==5]
# dt.m = dt.m
# output = 'best_impact'
# uw = c(2,1,1)
# simyear = 5
# quiet = FALSE
# nmax=1
# nopt=FALSE


# CHECKING results 1 NCU at a time - sim$total_impact[ncu==1830]

runDST <- function(db, dt.m, output = 'all',uw = c(1,1,1), simyear = 5, quiet = TRUE,nmax = NULL,nopt = TRUE,
                   yield_factor = 1, nsu_factor = 1){

  # make local copy of the INTEGRATOR/Eurostat database
  d2 <- copy(db)

  # add fraction for potential no-till area that overlaps with reduced till (both cannot be applied simultaneously)
  # assume that the potential corresponds to the proportion of NT/RT already applied (in NUTS2 zone)
  d2[,fr_rtnt := pmin(1,sum(parea.rtct,na.rm = T) / sum( parea.ntct,na.rm=T)), by = ncu]
  d2[is.na(fr_rtnt), fr_rtnt := 0]

  # subset only the columns needed to estimate the DISTANCE TO TARGET for Yield, SOC and N surplus
  d2 <- d2[,.(ncu,crop_name,NUTS2,area_ncu,area_ncu_ha,yield_ref,yield_target,density,
              soc_ref,soc_target,n_sp_ref,n_sp_sw_crit,n_sp_gw_crit,c_man_ncu,fr_rtnt)]
  d2[,area_ncu_ha_tot := sum(area_ncu_ha,na.rm = TRUE), by =.(ncu)]

  # re-arrange the table of MA models (object dt.m) to facilitate joining with integrator data
  dtm.mean <- dcast(dt.m,ncu + man_code ~ indicator, value.var = c('mmean','msd'))

  # merge the impact of measures from the MA models with the integrator db by NCU
  d3 <- merge(d2,dtm.mean,by=c('ncu'),allow.cartesian = TRUE)

  # correct effects for overlapping treatments via area correction on NCU level
  d3[man_code == "RT-CT" ,mmean_Nsu := mmean_Nsu * (1 - fr_rtnt)]
  d3[man_code == "RT-CT" ,mmean_SOC := mmean_SOC * (1 - fr_rtnt)]
  d3[man_code == "RT-CT" ,mmean_Y := mmean_Y * (1 - fr_rtnt)]
  d3[man_code == "NT-CT" ,mmean_Nsu := mmean_Nsu * fr_rtnt]
  d3[man_code == "NT-CT" ,mmean_SOC := mmean_SOC * fr_rtnt]
  d3[man_code == "NT-CT" ,mmean_Y := mmean_Y * fr_rtnt]

  # add estimated change in indicators for a period of X years (given as argument simyear)
  # converts from percentage change (from the ML models) to multiplied factor for yield, SOC and N surplus

  d3[, dY   := yield_factor * mmean_Y   * 0.01]
  d3[, dSOC := simyear      * mmean_SOC * 0.01]
  d3[, dNsu := nsu_factor   * mmean_Nsu * 0.01]

  # prevent that the change in SOC due to manure input exceeds available C (kg ha-1 * 1000/ kg soil ha-1 = g C / kg soil * 0.1 = % SOC * 0.01 = fraction change)
  d3[grepl('^CF|^OF', man_code) & dSOC > 0, dSOC := pmin(dSOC, c_man_ncu * 1000 *.1/(100 * 100 *.25 * density))]

  # adaptation for missing/ NA models
  d3[is.na(dY),dY := 0]
  d3[is.na(dSOC),dSOC := 0]
  d3[is.na(dNsu),dNsu := 0]

  # add metric for INITIAL distance to target for inspecting maps
  d3[, dist_Y := yield_ref / yield_target ]
  d3[, dist_C := soc_ref / soc_target ]
  d3[, n_sp_crit := pmin(n_sp_sw_crit,n_sp_gw_crit)]
  d3[, dist_N := n_sp_ref / pmin(n_sp_sw_crit,n_sp_gw_crit)] #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  d3[is.na(dist_N), dist_N := 1]

  # estimate overall impact per measure & NCU given the different area coverage based on crop types
  # the outcome is a weighted mean based on sub-ncu crop areas for all parameters below
  cols <- colnames(d3)[grepl('^dY|^sY|^sSOC|^dSOC|^sNsu|^dNsu|^dist_Y|^dist_C|^dist_N|^yield_ref|^soc_ref|^n_sp_ref|^density|^yield_target|^soc_target|^n_sp_crit',colnames(d3))]
  d3 <- d3[,lapply(.SD,function(x) weighted.mean(x,area_ncu,na.rm=T)),.SDcols = cols,by=c('ncu','man_code','area_ncu_ha_tot')]

  # create a set with all combinations of measures

  # what is the list of available measures (man_code)
  measures <- unique(d3$man_code) # 7 measures

  # what is the desired number of allowed combinations of measures
  # when nmax=1, no combo allowed, nmax=2 means 2 measures can be applied together, etc.
  if(is.null(nmax)){nmax = length(measures)}

  # make a data.table with all possible combinations of the measures
  # when combo is 2, meas^2 = 7x7 = 49 options, etc.
  meas.combi <- do.call(CJ,replicate(nmax,measures,FALSE))

  # add an unique ID per measure combi (including duplicates)
  meas.combi[,cid := .I]

  # melt the data.table to facilitate calculations per NCU per measure combi
  meas.combi <- melt(meas.combi,id.vars ='cid',variable.name = 'measure')

  # remove duplicates
  meas.combi <- unique(meas.combi, by = c('cid','value'))

  # create an exponential number to enable unique sum of individual measures
  nfc <- 10^(1:length(measures))

  # add an unique value for each measure
  meas.combi[,value2 := nfc[as.factor(value)]]

  # sum the total values of the measures (so that each measure combination has an unique code)
  # cid unique code for each combination
  meas.combi[,value3 := sum(value2),by=cid]

  # add an unique ID for each unique measure combination
  meas.combi[,cgid := .GRP,by=.(value,value3)]

  # remove duplicates
  meas.combi <- unique(meas.combi,by = 'cgid')

  # update the unique ID per measurement combi
  # gives unique combo of measures that are possible
  meas.combi[,cgid := .GRP,by=.(value3)]

  # select only the relevant columns
  meas.combi <- meas.combi[,.(cgid,man_code = value)]

  # make an unique table with each combination of measures per unique ID
  dt.meas.combi <- meas.combi[,list(man_code = paste(man_code,collapse = '-'),
                                    man_n = .N),by='cgid']

  # remove options that can not be applied together
  # dt.meas.combi <- dt.meas.combi[!grepl('RT-CT-NT-CT|NT-CT-RT-CT',man_code)]
  # dt.meas.combi <- dt.meas.combi[!grepl('OF-MF-CF-MF|CF-MF-OF-MF',man_code)]
  # dt.meas.combi <- dt.meas.combi[!grepl('EE-OF-MF|OF-MF-EE',man_code)]

  # update combined restriction for 3 combinations
  dt.meas.combi <- dt.meas.combi[!(grepl('RT-CT',man_code) & grepl('NT-CT',man_code))]
  dt.meas.combi <- dt.meas.combi[!(grepl('OF-MF',man_code) & grepl('CF-MF',man_code))]
  dt.meas.combi <- dt.meas.combi[!(grepl('EE',man_code) & grepl('OF-MF',man_code))]


  # combine all measurement combinations per NCU
  dt <- merge.data.table(d3,meas.combi,by='man_code', all= TRUE,allow.cartesian = TRUE)


  # CALCULATE SCORES AND MEASURE ORDER

  # define output list to store results
  dt.out <- list()

  # make progress bar to know how the run time - at different points in the function it is updated while running
  if(!quiet) {pb <- txtProgressBar(min = 0, max = 1199, style = 3);j=0}

  # make a sequence to split the database
  ncu_steps <- unique(round(seq(0,max(dt$ncu),length.out = 60)))
  #--------
  # for loop goes through NCUs in subsets to save memory
  for(i in 2:length(ncu_steps)){

    # i=6
    # select the row numbers to subset
    ncu_min <- ncu_steps[i-1]
    ncu_max <- ncu_steps[i]

    # can define min/max when running separate
    # subset the dataset from the first ncu after the previous step up to the last ncu of this step
    dt.ss <- dt[ncu > ncu_min & ncu <= ncu_max]

    # when score 0 replace value with very low number, so that output is not zero from all the multiplying
    dt.ss[dY == 0, dY := 1e-4]
    dt.ss[dNsu == 0, dNsu := 1e-4]
    dt.ss[dSOC == 0, dSOC := 1e-4]

    # add random noise to avoid same rank for measures with similar impacts
    set.seed(123)
    cols <- c('dY','dSOC','dNsu')
    dt.ss[, c(cols) := lapply(.SD, function(x)  x * (1 + rnorm(.N,0,0.0001))),.SDcols = cols]

    # add a relative score for the impact of each measure based on meta-models only
    # rank of absolute effects of measures; this shows order of magnitude for each impact Y/C/N
    dt.ss[, c('odY','odSOC','oNsu') := lapply(.SD,function(x) frankv(abs(x),order=-1)),.SDcols = c('dY','dSOC','dNsu'),by=.(ncu,cgid)]
    #----------
    # update progress bar
    if(!quiet) {j = j+1; setTxtProgressBar(pb, j)}

    # score function
    scorefun <- function(v_change,v_ref,v_target,type='yield'){

      if(type %in% c('yield','soc')){out <- pmin(1,pmax(0,ifelse(v_ref <= v_target,v_change * v_ref  / (v_target - v_ref),0)))}
      if(type =='nsp'){out <- pmin(1,pmax(0,ifelse(v_ref <= v_target,0,v_change * v_ref  / (v_target - v_ref))))}
      return(out)
    }

    # estimate total change and the total "change" in view of distance to target, accounting for user weights
    dt.ss2 <- dt.ss[, list(dY = sum(dY/odY),
                           dSOC = sum(dSOC/odSOC),
                           dNsu = sum(dNsu/oNsu),
                           sY_combi = scorefun(v_change = sum(dY/odY),v_ref = yield_ref[1],v_target = yield_ref[1]/dist_Y[1],type='yield'),
                           sSOC_combi = scorefun(v_change = sum(dSOC/odSOC),v_ref = soc_ref[1],v_target = soc_ref[1]/dist_C[1],type='soc'),
                           sNsu_combi = scorefun(v_change = sum(dNsu/oNsu),v_ref=n_sp_ref[1],v_target = n_sp_ref[1]/dist_N[1],type='nsp'),
                           dist_Y=dist_Y[1],
                           dist_C=dist_C[1],
                           dist_N=dist_N[1],
                           yield_ref=yield_ref[1],
                           soc_ref=soc_ref[1],
                           n_sp_ref=n_sp_ref[1],
                           area_ncu_ha_tot=area_ncu_ha_tot[1],
                           bd=density[1],
                           yield_target=yield_target[1], #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                           soc_target=soc_target[1],
                           n_sp_crit=n_sp_crit[1]),
                    by=.(ncu,cgid)]
    #-------------
    # add total score of individual and combined measures, evaluated as total distance to target
    dt.ss2[,iuw := sum(uw)]
    dt.ss2[,bipmc := (uw[1]*sY_combi + uw[2]*sSOC_combi+uw[3]*sNsu_combi)/iuw]

    # add scaling function (value - min) / (max - min)        ==> when Scombi scores all ZERO for entire ncu (min & max)
    mmsfun <- function(x){(x - min(x))/(max(x+1e-4)-min(x))}     # ==> denominator=0 and bipmc=NaN / added 0.0001 to max

    # add total score of individual and combined measures including user weighing
    dt.ss2[iuw !=3, bipmc := (mmsfun(sY_combi) * uw[1] + mmsfun(sSOC_combi)* uw[2] + mmsfun(sNsu_combi)* uw[3])/iuw,by=.(ncu)]
    #-------------
    # add the names of measures taken back in from cgid (saves some memory)
    dt.ss2 <- merge(dt.ss2,dt.meas.combi,by='cgid',all.x = TRUE)

    # remove combinations of measures that are not allowed
    dt.ss2 <- dt.ss2[!is.na(man_code)]

    # add a correction score for the number of measures
    if(nopt == TRUE){dt.ss2[,nmcf := 0.05/man_n]} else {dt.ss2[,nmcf := 0]}

    # add some uncertainty on the scoring value (in particular when zero)
    dt.ss2[, bipmc := bipmc + sample(1:.N,.N)/1e10,by=ncu]

    # add a ranking based on the integral score for each ncu
    dt.ss2[,bipmcs := as.integer(frankv(bipmc+nmcf,order = -1)),by=ncu]

    # save into a list, this output gives all unique management combinations (cgid) per ncu along with score (bipmcs), impacts, distances to target
    # dt.out[[i]] <- copy(dt.ss2[,.(ncu,area_ncu_ha_tot,cgid,man_code,man_n,dY,dSOC,dNsu,dist_Y,dist_C,dist_N,bipmc,bipmcs,yield_ref,soc_ref,n_sp_ref,bd,yield_target,soc_target,n_sp_crit)])

    ## sY_combi/sSOC_combi/sNsu_combi exported for thesis
    dt.out[[i]] <- copy(dt.ss2[,.(ncu,area_ncu_ha_tot,cgid,man_code,man_n,dY,dSOC,dNsu,sY_combi,sSOC_combi,sNsu_combi,dist_Y,dist_C,dist_N,bipmc,bipmcs,yield_ref,soc_ref,n_sp_ref,bd,yield_target,soc_target,n_sp_crit)])
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  }


  # converts the list into a table (rowbind)
  dt.out <- rbindlist(dt.out)

  # estimate whether the target has been met on NCU level
  # dist_Y = current yield / yield target
  # if target already met OR change in indicator will bring it above/below
  # dt.out[,tm_Y := fifelse(dist_Y >= 1 | dY >= (1 - dist_Y), 1,0)]
  # dt.out[,tm_C := fifelse(dist_C >= 1 | dSOC >= (1 - dist_C),1,0)]
  # dt.out[,tm_N := fifelse(dist_N <= 1 | dNsu <= (1 - dist_N),1,0)]
  dt.out[,tm_Y := fifelse(dY < (1 - dist_Y),0,1)]
  dt.out[,tm_C := fifelse(dSOC < (1 - dist_C),0,1)]
  dt.out[,tm_N := fifelse(dNsu > (1 - dist_N),0,1)]


  # compare to initial targets met
  dt.out[,ti_Y := fifelse(dist_Y >= 1,1,0)]
  dt.out[,ti_C := fifelse(dist_C >= 1,1,0)]
  dt.out[,ti_N := fifelse(dist_N <= 1,1,0)]

  # dt.out[,dist_Yt :=  yield_ref_w / yield_targ_w ]
  # d2.yRT$dist_Yt = d2.yRT$yield_targ_w / d2.yRT$yield_ref_w
  # d2.yRT$diff_Y = d2.yRT$yield_targ_w - d2.yRT$yield_ref_w
  # d2.yRT$perc_Yt = (d2.yRT$yield_targ_w - d2.yRT$yield_ref_w)/d2.yRT$yield_targ_w

  # OUTPUT DATA COLLECTION

  # DT.OUT is always used which CONTAINS ALL DATA AND FILTERING/AGGREGATING/ADAPTING TO WHAT OUTPUT WE WANT
  # get familiar with dt.out structure to adapt these lines myself

  # collect relevant output for case that all measures have been applied
  if(sum(grepl('total_impact|all',output))>0){

    # select relevant data and sort
    pout1 <- dt.out[man_n == max(man_n),.(ncu,man_code,dY,dSOC,dNsu,tm_Y,tm_C,tm_N)]
    # pout1 <- dt.out[bipmcs==1,.(ncu,man_code,dY,dSOC,dNsu)]
    setorder(pout1,ncu)
  } else {pout1 = NULL}

  # the case where only best measure (or best combination of measures) have been applied
  # the best ranking is chosen, whether the best is one measure or a combination
  # show the changes in indicators as well as distance to target status

  # ***** Question *****
  # Within an NCU, reference values change by crop types. Y and Nsu by all crop types. SOC and BD only different for grassland
  # Where should area-weighted means for reference values be calculated?
  # ********************
  if(sum(grepl('best_impact|all',output))>0){

    # select relevant data and sort
    # BIPMCS MUST BE CHANGED TO 1 for best measure and 11 for worst measure
    # **May 2024 weighted reference values and total ha per ncu
    ## carry gap-closure scores through to the best-single output (out_best_n1)
    pout2 <- dt.out[bipmcs==1,.(ncu,area_ncu_ha_tot,man_code,dY,dist_Y,dSOC,dist_C,dNsu,dist_N,sY_combi,sSOC_combi,sNsu_combi,tm_Y,tm_C,tm_N,
                                ti_Y,ti_C,ti_N,yield_ref,soc_ref,n_sp_ref,bd,yield_target,soc_target,n_sp_crit,bipmc,bipmcs)]
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    setorder(pout2,ncu)
    ## runner-up (2nd-best single) per-indicator scores, for the SRQ2 selection-margin decomposition
    poutRU <- dt.out[bipmcs==2,.(ncu,man_code_ru=man_code,sY_ru=sY_combi,sSOC_ru=sSOC_combi,sNsu_ru=sNsu_combi,bipmc_ru=bipmc)]
    setorder(poutRU,ncu)

    ## full single-measure score ladder per NCU
    ## (margin vs next *distinct* alternative). man_n == 1
    ## keeps this to single measures under any nmax; bipmcs is the within-NCU rank.
    poutALL <- dt.out[man_n == 1,.(ncu,man_code,bipmcs,bipmc,sY_combi,sSOC_combi,sNsu_combi)]
    setorder(poutALL,ncu,bipmcs)
    } else {pout2 = NULL; poutRU = NULL; poutALL = NULL}

  # collect the ORDER/RANKING of the single measures if applied separately
  # ranking of each individual measure in a column
  # filtering by single measures
  if(sum(grepl('score_single|all',output))>0){

    # select relevant data and sort
    pout3 <- dt.out[man_n == 1,.(ncu,man_code,bipmcs,dist_Y,dist_C,dist_N)][,bipmcs := frankv(bipmcs),by=ncu]
    #,tm_Y,tm_C,tm_N
    # change into table format (with the number varying from 1 (the best) to 7 (the lowest impact))
    #  TO ADD THOSE columns
    pout3 <- dcast(pout3,ncu+dist_Y+dist_C+dist_N~bipmcs,value.var = 'man_code')
    #+tm_Y+tm_C+tm_N
    # pout3 <- dcast(pout3,ncu+dist_Y+dist_C+dist_N+tm_Y+tm_C+tm_N~bipmcs,value.var = 'man_code')
  } else {pout3 = NULL}

  # collect the ORDER/RANKING of TWO combinations of measures
  # selects the output where nmax is 2, ranks by bipmcs metric, selects #1
  if(sum(grepl('score_duo|all',output))>0 & nmax >= 2){

    # select relevant data and sort
    #SELECT NO. OF MEASURES 2; SELECT COLUMNS; OVERWRITE COLUMN BIPMCS (ORDER OF ALL COMBO MEASURES)
    #SO MAKING A SUBSET OF THIS BASED ON A CRITERIA; THEN RECALCULATE THE ORDER BASED ON "MISSING RANK VALUES"
    # **May 2024 weighted reference values and total ha per ncu
    pout4 <- dt.out[man_n == 2,.(ncu,area_ncu_ha_tot,man_code,bipmc,bipmcs,dY,dSOC,dNsu,dist_Y,dist_C,dist_N,tm_Y,tm_C,tm_N,
                                 ti_Y,ti_C,ti_N,yield_ref,soc_ref,n_sp_ref,bd)][,bipmcs := frankv(bipmcs),by=ncu]
    pout4 = pout4[bipmcs==1]
    # change into table format (with the number varying from 1 (the best) to 7 (the lowest impact))
    #<<<<<<< HEAD
    # pout4 <- dcast(pout4,ncu+dist_Y+dist_C+dist_N+tm_Y+tm_C+tm_N~bipmcs,value.var = 'man_code')
    #=======
    pout4a <- dcast(pout4,ncu+dist_Y+dist_C+dist_N~bipmcs,value.var = 'man_code')

    pout4[,tm_all := tm_Y + tm_C + tm_N]
    pout4b <- dcast(pout4,ncu~man_code,value.var = 'tm_all')

    pout4 = pout4[bipmcs==1]
    pout4[,]
    #>>>>>>> 4b2de830595cff818eaf3e359f6d23e71b5119d2

    # pout4 <- dt.out[man_n == 2,.(ncu,man_code,bipmcs,dY,dSOC,dNsu,dist_Y,dist_C,dist_N)][,bipmcs := frankv(bipmcs),by=ncu]
    # pout4 <- dt.out[bipmcs==1,.(ncu,man_code,dY,dist_Y,dSOC,dist_C,dNsu,dist_N)]
    # pout4 <- dcast(pout4,ncu+dist_Y+dist_C+dist_N~bipmcs,value.var = 'man_code')

  } else {pout4 = NULL}

  # collect the ORDER/RANKING of THREE combinations of measures
  # selects the output where nmax is 3, ranks by bipmcs metric, selects #1
  if(sum(grepl('score_trio|all',output))>0 & nmax >= 3){

    # **May 2024 weighted reference values and total ha per ncu
    pout5 <- dt.out[man_n == 3,.(ncu,area_ncu_ha_tot,man_code,bipmc,bipmcs,dY,dSOC,dNsu,dist_Y,dist_C,dist_N,tm_Y,tm_C,tm_N,
                                 ti_Y,ti_C,ti_N,yield_ref,soc_ref,n_sp_ref,bd)][,bipmcs := frankv(bipmcs),by=ncu]
    pout5 = pout5[bipmcs==1]

    # change into table format (with the number varying from 1 (the best) to 7 (the lowest impact))
    # pout5 <- dcast(pout5,ncu+dist_Y+dist_C+dist_N+tm_Y+tm_C+tm_N~bipmcs,value.var = 'man_code')
  } else {pout5 = NULL}

  # THIS IS SAME AS POUT2 BUT FOR BEST COMBO
  # if(sum(grepl('score_best|all',output))>0){
  #
  #   # select relevant data and sort
  #   pout5 <- dt.out[bipmcs==1,.(ncu,man_code,dist_Y,dist_C,dist_N)]
  #   setorder(pout5,ncu)
  # } else {pout5 = NULL}


  # update progress bar
  if(!quiet) {j = j+1; setTxtProgressBar(pb, j)}

  # select relevant output to be returned
  out <- list(impact_total = pout1,
              impact_best = pout2,
              impact_runnerup = poutRU,
              impact_allscores = poutALL,
              score_single = pout3,
              score_duo = pout4,
              score_trio = pout5)

  # update progress bar
  if(!quiet) {j = j+1; setTxtProgressBar(pb, j);close(pb)}

  # return output
  return(out)
}

#' Function to run the DST including Monte Carlo simulations to estimate the mean and stand deviations of the DST output

#' @param db (data.table) a data.table with INTEGRATOR data used for the DST
#' @param measures (data.table) a vector with different management options to be estimated. Options include: c('CC','RES','ROT','NT-CT','RT-CT','CF-MF','OF-MF')
#' @param uw (vector) (user weight) is an argument that gives the relative importance for yield, soc or n-surplus (a value above 0, likely a fraction between 0 and 2). Default is c(1,1,1).
#' @param nsim (integer) the number of Monte Carlo simulations to be runned
#' @param mam (list) a list containing the generic meta-analytical model information
#' @param simyear (num) value for the default year for simulation of impacts (default is 5 years)
#' @param covar (boolean) is it desired to make use of the ma-models including covariates. Options: TRUE/FALSE.
#' @param output (string) optional arguments to select different types of outputs. Options include: 'total_impact','best_impact','score_single','score_due','score_best' or 'all'.
#' @param nmax (integer) the max number of measure combinations evaluated
#'
#' @export

# RUNNING DST A SET NUMBER OF TIMES AND GENERATING ERROR OF MEAN OUTPUT (BASED ON VARIANCE OF INPUTS)
# DEFAULTS SET HERE BUT CAN ADAPT WHEN RUNNING IN IN DST MAIN

# for troubleshooting
# db = d1
# mam = ma_models
# nsim = 4
# covar = FALSE
# simyear = 5
# uw = c(1,1,1)
# measures = c('CF-MF','OF-MF','EE','RFR','RFT','RFP')
# output = 'all'
# nmax=1

runMC_DST <- function(db, uw = c(1,1,1),simyear = 5,
                      measures = c('CF-MF','OF-MF','EE','RFR','RFT','RFP'),
                      mam = ma_models, nsim = 10, covar = FALSE,
                      output = 'all',nmax = 2){

  #we can adapt nsim later to be what needed (at least 100 but needs a lot of time)
  #output is the different options in runDST

  # make local copy of the input database integrator
  d1.rmc <- copy(db)

  # make list to store output of join MA- impact models per management measure as well ass Monte Carlo simulations
  dt.list = sim.list = sim.list1 = sim.list2 = sim.list3 = sim.list4 = sim.list5 = list()

  # function to calculate the model
  #table is counting frequency of occurance within a column
  #stores the most frequent option
  #FOR CATEGORICAL VARIABLES BECAUSE NO NUMERIC VALUE. HOW OFTEN ARE EACH RANKED HIGHEST, ETC.
  fmod <- function(x) as.integer(names(sort(table(x),decreasing = T)[1]))
  fmods <- function(x) names(sort(table(x),decreasing = T)[1])

  # make progress bar
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)


  # do Monte Carlo simulations
  for(i in 1:nsim){

    dt.list = list()

    # join MA-impact models per management measure (those defined above in measures)
    # cIMAm gives different output each time if mc set to true; here running various sims
    # list of separate data tables per measure
    for(j in measures){dt.list[[j]] <- cIMAm(management= j,db = d1.rmc, mam = mam, covar = covar, montecarlo = TRUE)}

    # combine all measures and their impacts into one data.table
    dt.m <- rbindlist(dt.list)

    sim <- runDST(db = d1.rmc, dt.m = dt.m, output = output,uw = uw,simyear = simyear,
                  quiet = FALSE, nmax = nmax)

    # ?? CHECKING 1 NCU sim$total_impact[ncu==1830]

    #adding column with sim number
    if(output !='all'){
      #impact_best is currently the results used
      sim.list[[i]] <- copy(sim$impact_total)[,sim:=i]
    } else {
      #if you want all, multiple outputs each in a list
      sim.list1[[i]] <- sim$impact_total[,sim :=i]
      sim.list2[[i]] <- sim$impact_best[,sim :=i]
      sim.list3[[i]] <- sim$score_single[,sim :=i]
      if(is.null(sim$score_duo)){
        sim.list4[[i]] <- NULL
      } else {
        sim.list4[[i]] <- sim$score_duo[,sim :=i]
      }

      sim.list5[[i]] <- sim$score_best[,sim :=i]
    }

    # update progress bar (just for simulation)
    setTxtProgressBar(pb, i)

  }

  # combine all in three objects
  if(output !='all'){
    out <- rbindlist(sim.list)
  } else{
    #summarizing the outputs
    out1 <- rbindlist(sim.list1)
    out1 <- out1[,list(dYmean = mean(dY), dSOCmean = mean(dSOC),dNsumean = mean(dNsu),
                       dYsd = sd(dY),dSOCsd = sd(dSOC), dNsusd = sd(dNsu)), by= .(ncu,man_code)]
    out2 <- rbindlist(sim.list2)
    out2 <- out2[,list(dYmean = mean(dY), dSOCmean = mean(dSOC),dNsumean = mean(dNsu),
                       dYsd = sd(dY),dSOCsd = sd(dSOC), dNsusd = sd(dNsu)), by= .(ncu,man_code)]
    out3 <- rbindlist(sim.list3)
    out3 <- out3[,lapply(.SD,fmod),by=ncu]
    out4 <- rbindlist(sim.list4)
    if(nrow(out4)>0){out4 <- out4[,lapply(.SD,fmod),by=ncu]}
    out5 <- rbindlist(sim.list5)
    out5 <- out5[, list(modal = fmods(man_code)),by=c('ncu')]

    # combine all output
    out <- list(impact_total = out1,
                impact_best = out2,
                score_single = out3,
                score_duo = out4,
                score_best = out5)
  }

  # close progressbar
  close(pb)

  # return output
  return(out)
}

#' Function to to link Integrator to MA models being available.

#' @param db (data.table) a data.table with INTEGRATOR data used for the DST
#' @param management (string) a management code with different management options to be estimated. Options include: c('CC','RES','ROT','NT-CT','RT-CT','CF-MF','OF-MF')
#' @param mam (list) a list containing the generic meta-analytical model information
#' @param covar (boolean) is it desired to make use of the ma-models including covariates. Options: TRUE/FALSE.
#' @param montecarlo (boolean) is it desired to make use of MC estimates of the ma-model estimates. Options: TRUE/FALSE.
#'
#' @export
#

# management = 'EE'
# db = d1
# mam = ma_models
# montecarlo = FALSE
# covar = TRUE

cIMAm <- function(management,db = d1, mam = ma_models,montecarlo = FALSE, covar = FALSE){
  # ADAPTED FROM MAIN CODE TO MATCH TOTAL NCU AREA TO THE PRACTICE (WE DO NOT NEED BASELINE MANAGEMENT /CROP AREAS)
  #STILL NEED EFFECT OF FERTILIZER VERSUS NO FERTILIZER SO THIS MIGHT AFFECT THE CODE

  # add nitrogen use calculations to update the potential applicable measures for fertilizers
  # add area for EE, RFT, RFP and RFR to the database
  db[,nup := n_fert+n_man + n_fix + n_dep - n_sp_ref]
  db[,nue := nup/(n_fert+n_man + n_fix + n_dep)]
  db[,parea.rfr := pmin(1,pmax(0,1-nue)) * area_ncu_ha]
  db[,parea.ee := pmin(1,yield_target / yield_ref - 1) * area_ncu_ha]
  db[,parea.rft := pmin(1,yield_target / yield_ref - 1) * area_ncu_ha]
  db[,parea.rfp := pmin(1,yield_target / yield_ref - 1) * area_ncu_ha]

  # select the right column for the area
  # one of the following 11 measures is selected when running the function
  if(management=='CC'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.cc)] # just selecting area from integrator
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = parea.cc)] # separate selection adding covariate property
  } else if(management=='RES'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.cres)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = parea.cres)]
  } else if(management=='ROT'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.cr)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = parea.cr)]
  } else if(management=='NT-CT'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.ntct)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = parea.ntct)]
  } else if(management=='RT-CT'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.rtct)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = parea.rtct)]
  } else if(management=='EE'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.ee)] # add total ncu area from integrator
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = area_ncu_ha)] # separate selection adding covariate property
  } else if(management=='RFP'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.rfp)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = area_ncu_ha)]
  } else if(management=='RFR'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.rfr)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = area_ncu_ha)]
  } else if(management=='RFT'){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.rft)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = area_ncu_ha)]
  } else if(management %in% c('OF-NF','CF-NF','MF-NF')){
    dt.m1 <- db[,.(ncu,ha_m1 = parea.nifnof)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = parea.nifnof)]
  } else if(management=='CF-MF'){
    # impact of partly replacing mineral fertilizers by organic ones
    # since impact is lower on highly fertilized soils, I halve the area (being identical to half the impact)
    # correction for highly fertilized soils- ha_m1 = parea.nifnof+parea.hifnof+0.75*parea.mifnof+0.75*parea.hifhof+0.25*parea.mifhof
    dt.m1 <- db[,.(ncu,ha_m1 = parea.nifnof+parea.hifnof+0.75*parea.mifnof+0.75*parea.hifhof+0.25*parea.mifhof)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = parea.nifnof+parea.hifnof+0.75*parea.mifnof+0.75*parea.hifhof+0.25*parea.mifhof)]
  } else if(management=='OF-MF'){
    # replacing inorganic with organic has only impact on soils with no organic input and limited inorganic
    # correction for highly fertilized soils- ha_m1 = parea.nifnof+parea.hifnof+0.75*parea.mifnof+0.75*parea.hifhof+0.25*parea.mifhof
    dt.m1 <- db[,.(ncu,ha_m1 = parea.nifnof+parea.hifnof+0.75*parea.mifnof+0.75*parea.hifhof+0.25*parea.mifhof)]
    dt.cov.m1 <- db[,.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc, ha_m1 = parea.nifnof+parea.hifnof+0.75*parea.mifnof+0.75*parea.hifhof+0.25*parea.mifhof)]
  }

  #OUTPUT of 2 columns and 7 columns specifying areas within each ncu matched to

  # unlist the meta-analytical models (extract from list of objects already loaded in workspace)
  ma_mean <- mam$ma_mean
  ma_sd <- mam$ma_sd
  ma_cov_mean <- mam$ma_cov_mean
  ma_cov_sd <- mam$ma_cov_sd

  # join the integrator database with the main meta-analytical models
  # adds MA models to the ncu areas
  dt.m1 <- dt.m1[,.(ha_m1 = sum(ha_m1),man_code=management),by=ncu] #selecting area and adding man code
  dt.m1 <- merge(dt.m1,ma_mean[man_code==management],by='man_code',all.x = TRUE) #merging mean estimate
  dt.m1 <- merge(dt.m1,ma_sd[man_code==management],by='man_code',all.x = TRUE) #again for SD

  #extract names for indicators (N, C, Y), by removing mean/sd text - selects those present in your db
  cols <- unique(gsub('^mean_','',colnames(dt.m1)[grepl('^mean',colnames(dt.m1))]))
  cols <- unique(gsub('^mean_|^sd_','',colnames(dt.m1)[grepl('^mean|^sd',colnames(dt.m1))]))

  #changing columns to row for each indicator (instead of wide it becomes long with a variable to specify the indicator)
  #also col for mean and sd
  # output is 3 rows per NCU with
  dt.m1 <- melt(dt.m1,id.vars = c('ncu','ha_m1'),
                measure=patterns("^mean_", "^sd_"),
                variable.factor = TRUE,
                value.name = c('mean','sd'),
                variable.name = 'indicator')
  dt.m1[,indicator := cols[as.integer(indicator)]] #adds name for each number assigned to the indicator
  dt.m1[ha_m1 == 0,  c('mean','sd') := list(mean = 0, sd = 0)]


  # estimate the mean impact given normal distribution
  if(montecarlo==TRUE){
    suppressWarnings(dt.m1[,mean := rnorm(.N,mean[1],sd = sd[1]),by='indicator'])
    #create normal distribution (N=number samples or rows for each indicator); select first row because all the same
    #using mean and sd of global models for each indicator-measure combon as basis for MC simulation
    dt.m1[ha_m1==0,mean := 0]
  }

  #some NCUs have no area for a measure;if so the impact should be zero because in end we sum impacts by area to get total impact
  #at this point connected: NCU + area of measure potentially applied + impacts of measures

  #similar operation for covariate models
  if(covar==TRUE){

    # join the integrator database with the meta-analytical models that uses covariables
    dt.cov.m1 <- dt.cov.m1[,.(ha_m1 = sum(ha_m1),unid = .GRP),by=.(ncu,cov_soil,cov_clim,cov_crop,cov_fert,cov_soc)]
    #unique id for all combos of covar with area
    dt.cov.m1[,man_code := management] #input management code you selected
    dt.cov.m1 <- melt(dt.cov.m1,id.vars = c('ncu','ha_m1','man_code','unid'),value.name = 'group') #make columns to rows
    dt.cov.m1[,variable := NULL] #removed column for covariable name to save memory because codes for each group is unique to variable
    dt.cov.m1 <- merge(dt.cov.m1,ma_cov_mean[man_code==management],by=c('group','man_code'),all.x = TRUE)
    #long list of models and filter the management code and join based on right groups (e.g. CC and cov_crop )
    dt.cov.m1 <- unique(dt.cov.m1)[,mods:=NULL] #remove some duplicated NCUs
    dt.cov.m1 <- merge(dt.cov.m1,ma_cov_sd[man_code==management],by=c('group','man_code'),all.x = TRUE)
    #same long list operation for SD
    dt.cov.m1 <- unique(dt.cov.m1)

    #now have predicted change for models including covariates

    #takes an average of the covariate models applicaple for that NCU
    cols <- colnames(dt.cov.m1)[grepl('mean_|sd_',colnames(dt.cov.m1))]
    #not sure why a weighted mean used? is it for different areas of crops in each NCU?
    #considering coarse groups of meta-models, there are different types of each group in the list per NCU
    dt.cov.m1 <- dt.cov.m1[,lapply(.SD,function(x) weighted.mean(x,w=ha_m1,na.rm=T)),.SDcols = c(cols),by=c('ncu','man_code')]
    # if area missing then set to zero change, and add some noise to replace the NA (being zero now) later
    dt.cov.m1[,c(cols):= lapply(.SD,function(x) fifelse(is.na(x),0,x+0.00001)),.SDcols = cols] #GERARD CHECKS THIS LINE (SEE BELOW)
    #reshaping columns to rows
    cols <- unique(gsub('^mean_|^sd_','',colnames(dt.cov.m1)[grepl('^mean|^sd',colnames(dt.cov.m1))]))
    dt.cov.m1 <- melt(dt.cov.m1,id.vars = c('ncu','man_code'),
                      measure=patterns("^mean_", "^sd_"),
                      variable.factor = TRUE,
                      value.name = c('cov_mean','cov_sd'),
                      variable.name = 'indicator')
    dt.cov.m1[,indicator := cols[as.integer(indicator)]] #assigning indicator name instead of integer

    #
    # estimate the mean impact given normal distribution
    if(montecarlo==TRUE){
      suppressWarnings(dt.cov.m1[,cov_mean := rnorm(.N,cov_mean[1],sd = cov_sd[1]),by=indicator])
      dt.cov.m1[cov_sd==0,cov_mean := 0]
    }



    # merge the models with and without covariates and select covariates when available
    # we select the available covar model over the global model; when covar missing take global one
    # covar model at this point is average of the site properties applicable for that NCU
    dt.m1.fin <- merge(dt.m1,dt.cov.m1,by=c('ncu','indicator'),all.x=TRUE)
    dt.m1.fin[,mmean := fifelse(cov_mean==0,mean,cov_mean)] #IF COVARIATE MODEL IS NA THEN TAKE GLOBAL; BUT ABOVE CHANGED NAs to 0
    dt.m1.fin[,msd := fifelse(cov_mean==0,sd,cov_sd)] #MADDY NOW CHANGED IT TO REPLACE ZERO WITH GLOBAL MODELS
    dt.m1.fin <- dt.m1.fin[,.(ncu,ha_m1,indicator,man_code,mmean,msd)] #define final output = #NCUs * #indicators = 88000

  } else {
    dt.m1.fin <- dt.m1[,.(ncu,ha_m1,indicator,man_code = management,mmean = mean,msd = sd)]
  }


  # return output
  return(dt.m1.fin)
}


# function to load the default meta-analytical models from a given filename
# and the different models are combined using inverse weighting on SD

#' @param fname (string) a filename where to find the excel with ma model results
#'
#' @export
lmam <- function(fname){

  # load meta-analytical models for main models without covariates
  m1.main <- as.data.table(readxl::read_xlsx(fname,sheet='main'))
  setnames(m1.main,gsub('-1','',colnames(m1.main)))

  # remove cases that are not correct
  m1.main<- m1.main[!(n==1 & dyr==1 & SEyr==0.1)] # "markers" in the template file

  # load meta-analytical models for main models with covariates, and select only relevant columns
  m1.covar <- as.data.table(readxl::read_xlsx(fname,sheet='covariates'))
  setnames(m1.covar,gsub('-1','',colnames(m1.covar)))
  m1.covar <- m1.covar[,.(ind_code, mods = `moderator/factor` ,man_code,group,SEyr,dyr,n)]

  # remove models without SE
  m1.covar <- m1.covar[!is.na(SEyr)]

  # add simplified grouping for moderators to simplify joining with integrator data later
  # groups are unique
  m1.covar[grepl('texture',mods), mods := 'cov_soil']
  m1.covar[grepl('crop',mods), mods := 'cov_crop']
  m1.covar[grepl('SOC',mods), mods := 'cov_soc']
  m1.covar[grepl('rate',mods), mods := 'cov_fert']
  m1.covar[grepl('pH',mods), mods := 'cov_ph']
  m1.covar[grepl('durat',mods), mods := 'cov_duration']
  m1.covar[grepl('climate',mods), mods := 'cov_clim']
  m1.covar <- m1.covar[!mods %in% c('cov_ph','cov_duration')]
  m1.covar[grepl('cov_soil',mods) & grepl('coarse|sand',group), group2 := 'coarse']
  m1.covar[grepl('cov_soil',mods) & grepl('medium|loam',group), group2 := 'medium']
  m1.covar[grepl('cov_soil',mods) & grepl('fine|clay|silt',group), group2 := 'fine']
  m1.covar[grepl('cov_crop',mods) & grepl('maiz',group), group2 := 'maize']
  m1.covar[grepl('cov_crop',mods) & grepl('cerea',group), group2 := 'cereals']
  m1.covar[grepl('cov_crop',mods) & grepl('root',group), group2 := 'rootcrops']
  m1.covar[grepl('cov_clim',mods) & grepl('subtropical',group), group2 := 'subtropical']
  m1.covar[grepl('cov_clim',mods) & grepl('temperate',group), group2 := 'temperate']
  m1.covar[grepl('cov_soc|cov_fert',mods) & grepl('low|no fert',group), group2 := 'low']
  m1.covar[grepl('cov_soc|cov_fert',mods) & grepl('medium',group), group2 := 'medium']
  m1.covar[grepl('cov_soc|cov_fert',mods) & grepl('high',group), group2 := 'high']
  m1.covar[grepl('cov_fert',mods),group2 := paste0('f',group2)]
  m1.covar[grepl('cov_soc',mods),group2 := paste0('c',group2)]
  m1.covar <- m1.covar[!is.na(group2)][,group := group2][,group2 := NULL]

  # calculate the weighed mean models per measure
  m1.main[,wm.sd := 1/sqrt(sum(1 / SEyr^2)) , by=c('ind_code','man_code')]
  m1.main[,wm.mean := sum((dyr / SEyr^2)) / sum(1 / SEyr^2) , by=c('ind_code','man_code')]
  m1.covar[,wm.sd := 1/sqrt(sum(1 / SEyr^2)) , by=c('ind_code','man_code','mods','group')]
  m1.covar[,wm.mean := sum((dyr / SEyr^2)) / sum(1 / SEyr^2) , by=c('ind_code','man_code','mods','group')]

  # remove duplicates
  m2.main <- unique(m1.main[,.(ind_code,man_code,wm.mean,wm.sd)])
  m2.covar <- unique(m1.covar[,.(ind_code,man_code,mods,group,wm.mean,wm.sd)])

  # mean model responses for main ma-models (without covariates)
  m2.mean <- dcast(m2.main,man_code ~ ind_code, value.var = 'wm.mean')
  setnames(m2.mean,colnames(m2.mean)[-1],paste0('mean_',colnames(m2.mean)[-1]))
  m2.sd <- dcast(m2.main,man_code ~ ind_code, value.var = 'wm.sd')
  setnames(m2.sd,colnames(m2.sd)[-1],paste0('sd_',colnames(m2.sd)[-1]))

  # mean model responses for main ma-models with covariates
  m2.cov.mean <- dcast(m2.covar, man_code + mods + group ~ ind_code, value.var = 'wm.mean')
  setnames(m2.cov.mean,colnames(m2.cov.mean)[-c(1:3)],paste0('mean_',colnames(m2.cov.mean)[-c(1:3)]))
  m2.cov.sd <- dcast(m2.covar, man_code + mods + group ~ ind_code, value.var = 'wm.sd')
  setnames(m2.cov.sd,colnames(m2.cov.sd)[-c(1:3)],paste0('sd_',colnames(m2.cov.sd)[-c(1:3)]))

  # Monte Carlo requires SD values, which are missing for covariates
  # Here I tried to remove them to see if that
  # if (covar==FALSE) {
  #   # add the meta-analytical models into one list
  #   out = list(ma_mean = m2.mean, ma_sd = m2.sd)
  #
  # } else {
  # add the meta-analytical models into one list
  out = list(ma_mean = m2.mean, ma_sd = m2.sd,ma_cov_mean = m2.cov.mean, ma_cov_sd = m2.cov.sd)
  # }


  # return
  return(out)
}

# *************************************************************************************************
# helper function to re-allocate cropping categories on NUTS level to NCU level given likelihoods
# *************************************************************************************************

dsRA <- function(area_ncu,area_nuts2,lkh, nuts2){

  # make local copy of the inputs
  dt <- data.table(area_ncu, area_nuts2,lkh, nuts2)

  # convert unit for integrator NCU (km2) to match eurostat (ha)
  dt[,area_ncu := area_ncu * 100]

  # add unique id
  dt[, id := .I]

  # estimate total likely area on NUTS2 level given NCU crop area and likelihood for measure X
  dt[,a_1 := sum(lkh * area_ncu, na.rm = TRUE), by = nuts2]

  # correct this proportial to the area that is likely used for management measure X
  # the area of some NCUs already maximized/filled up, so cannot go over the total crop area
  # if likely area is less than NCU area, it saves this (see below)
  # gives an NCU likely area which is proportional to total area at the NUTS2 level
  dt[,a_2 := pmin(area_ncu,lkh * area_ncu * area_nuts2 / a_1,na.rm = TRUE), by = nuts2]

  # how much NUTS2 area leftover should still be reallocated for management measure X
  dt[,a_3 := area_nuts2 - sum(a_2), by = nuts2]

  # how much NCU area is still left within crops that have a likelihood for measure X
  # calculates area where more may be applied, if lkh is 0 then assigns 0 value
  dt[,a_4 := fifelse(lkh > 0,(area_ncu - a_2),0), by = nuts2]

  # redistribute the area that need to be reallocated to the area that is left
  # calculates exact proportion of NUTS2 area to go to each; if lkh is 0 then a_5 becomes 0
  dt[,a_5 := a_3 * a_4/sum(a_4),by = nuts2]

  # estimate total area of measure X to add to each NCU; if no more room then 0 added
  dt[,a_6 := a_2 + pmax(0,a_5,na.rm=T)]  # ??? why is pmax function needed ???

  # ensure correct order
  setorder(dt,id)

  # extract value
  value <- dt[,a_6]

  # return value
  return(value)

}

# *************************************************************************************************
# updated helper function to re-allocate cropping categories on NUTS level to NCU level given likelihoods
# *************************************************************************************************

dsRA_new <- function(area_ncu,area_nuts2,lkh, nuts2,area_nuts2_all = NULL){

  # make local copy of the inputs
  dt <- data.table(area_ncu, area_nuts2,lkh, nuts2,area_nuts2_all)

  # assume that total area for the measure equals the measure to be reallocated
  if(is.null(area_nuts2_all)){
    dt[,area_nuts2_all := area_nuts2]
  }

  # convert unit for integrator NCU (km2) to match eurostat (ha)
  dt[,area_ncu := area_ncu * 100]

  # add unique id
  dt[, id := .I]

  # estimate the total area in ncu on nutslevel
  dt[,a_0 := sum(area_ncu,na.rm = T),by=nuts2]

  # estimate total likely area on NUTS2 level given NCU crop area and likelihood for measure X
  dt[,a_1 := sum(lkh * area_ncu, na.rm = TRUE), by = nuts2]

  # correct this proportial to the area that is likely used for management measure X
  # the area of some NCUs already maximized/filled up, so cannot go over the total crop area
  # if likely area is less than NCU area, it saves this (see below)
  # gives an NCU likely area which is proportional to total area at the NUTS2 level
  dt[,a_2 := pmin(area_ncu,lkh * area_ncu * area_nuts2_all / a_1,na.rm = TRUE), by = nuts2]

  # how much NUTS2 area leftover should still be reallocated for management measure X
  dt[,a_3 := area_nuts2_all - sum(a_2), by = nuts2]

  # how much NCU area is still left within crops that have a likelihood for measure X
  # calculates area where more may be applied, if lkh is 0 then assigns 0 value
  dt[,a_4 := fifelse(lkh > 0,(area_ncu - a_2),0), by = nuts2]

  # redistribute the area that need to be reallocated to the area that is left
  # calculates exact proportion of NUTS2 area to go to each; if lkh is 0 then a_5 becomes 0
  dt[,a_5 := a_3 * a_4/sum(a_4),by = nuts2]

  # estimate total area of measure X to add to each NCU; if no more room then 0 added
  dt[,a_6 := a_2 + pmax(0,a_5,na.rm=T)]

  # estimate the fraction of measure x
  dt[,afin := a_6 * area_nuts2 / area_nuts2_all]

  # ensure correct order
  setorder(dt,id)

  # extract value
  value <- dt[,afin]

  # return value
  return(value)

}
