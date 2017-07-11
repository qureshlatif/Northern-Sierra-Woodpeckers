setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")

#################################################
# Tabulate descriptive statistics of predictors #
#################################################

library(dplyr)

load("Data_compiled.RData")

## Summary for all species ##
tab.sum <- dat %>% 
  mutate(SPECIES = replace(SPECIES, is.na(SPECIES), "RAND")) %>%
  filter(SPECIES %in% c("BBWO", "HAWO", "WHWO", "NOFL", "RBSA", "MOBL", "RAND")) %>%
  group_by(SPECIES) %>%
  summarise(DBH = paste0(round(mean(DBH), digits = 1), " (", round(sd(DBH), digits = 1), ")"),
            SP_FIR = round(sum(SP_FIR)/n(), digits = 2),
            SP_PINE = round(sum(SP_PINE)/n(), digits = 2),
            SNAG = round(sum(SNAG)/n(), digits = 2),
            SNAG_INT = round(sum(SNAG_INT)/n(), digits = 2),
            SNAG_DEC = round(sum(SNAG_DEC)/n(), digits = 2),
            BRKN_before = round(sum(BRKN_before)/n(), digits = 2),
            BRKN_after = round(sum(BRKN_after)/n(), digits = 2),
            BRKN = round(sum(BRKN)/n(), digits = 2),
            SnagDens_23to38 = paste0(round(mean(SnagDens_23to38), digits = 1), " (",
                                     round(sd(SnagDens_23to38), digits = 1), ")"),
            SnagDens_38to50 = paste0(round(mean(SnagDens_38to50), digits = 1), " (",
                                     round(sd(SnagDens_38to50), digits = 1), ")"),
            SnagDens_gt50 = paste0(round(mean(SnagDens_gt50), digits = 1), " (",
                                     round(sd(SnagDens_gt50), digits = 1), ")"),
            SnagDens_pine = paste0(round(mean(SnagDens_pine), digits = 1), " (",
                                     round(sd(SnagDens_pine), digits = 1), ")"),
            SnagDens_fir = paste0(round(mean(SnagDens_fir), digits = 1), " (",
                                     round(sd(SnagDens_fir), digits = 1), ")"),
            SnagDens_all = paste0(round(mean(SnagDens_all), digits = 1), " (",
                                  round(sd(SnagDens_all), digits = 1), ")"),
            slope = paste0(round(mean(slope), digits = 1), " (", round(sd(slope), digits = 1), ")"),
            sinasp = paste0(round(mean(sinasp), digits = 1), " (", round(sd(sinasp), digits = 1), ")"),
            cosasp = paste0(round(mean(cosasp), digits = 1), " (", round(sd(cosasp), digits = 1), ")"),
            ccmort_loc = paste0(round(mean(ccmort_loc), digits = 1), " (",
                                round(sd(ccmort_loc), digits = 1), ")"),
            blk_lndcc = paste0(round(mean(blk_lndcc), digits = 1), " (",
                               round(sd(blk_lndcc), digits = 1), ")"),
            canlo_loc = paste0(round(mean(canlo_loc), digits = 1), " (",
                               round(sd(canlo_loc), digits = 1), ")"),
            canlo_lnd = paste0(round(mean(canlo_lnd), digits = 1), " (",
                               round(sd(canlo_lnd), digits = 1), ")"),
            canmd_loc = paste0(round(mean(canmd_loc), digits = 1), " (",
                               round(sd(canmd_loc), digits = 1), ")"),
            canmd_lnd = paste0(round(mean(canmd_lnd), digits = 1), " (",
                               round(sd(canmd_lnd), digits = 1), ")"),
            canhi_loc = paste0(round(mean(canhi_loc), digits = 1), " (",
                               round(sd(canhi_loc), digits = 1), ")"),
            canhi_lnd = paste0(round(mean(canhi_lnd), digits = 1), " (",
                               round(sd(canhi_lnd), digits = 1), ")"),
            sizsm_loc = paste0(round(mean(sizsm_loc), digits = 1), " (",
                                round(sd(sizsm_loc), digits = 1), ")"),
            sizlrg_loc = paste0(round(mean(sizlrg_loc), digits = 1), " (",
                                round(sd(sizlrg_loc), digits = 1), ")"),
            sizlrg_lnd = paste0(round(mean(sizlrg_lnd), digits = 1), " (",
                                round(sd(sizlrg_lnd), digits = 1), ")"),
            fir_1km = paste0(round(mean(fir_1km), digits = 1), " (",
                             round(sd(fir_1km), digits = 1), ")"))

ord <- c("RAND", "BBWO", "HAWO", "WHWO", "NOFL", "RBSA", "MOBL") # set desired order
tab.sum <- tab.sum %>%
  slice(match(ord, SPECIES))

write.csv(tab.sum, "Predictor_summary_stats.csv", row.names = F)

## Summaries by species ##
spp <- c("BBWO", "HAWO", "WHWO", "NOFL", "RBSA", "MOBL")
for(sp in spp) {
  dta <- eval(as.name(paste0("dat.",sp)))
  vars <- names(dta)[which(names(dta) %in% dimnames(scale.factors)[[1]])]
  for(v in vars) dta[, v] <- dta[, v]*scale.factors[v, "SD"] + scale.factors[v, "mean"]
  if(sp %in% c("BBWO", "HAWO", "WHWO", "NOFL")) {
    dta <- dta %>% 
      group_by(TYPE) %>%
      summarise(DBH = paste0(round(mean(DBH), digits = 1), " (", round(sd(DBH), digits = 1), ")"),
                SP_FIR = round(sum(SP_FIR)/n(), digits = 2),
                SP_PINE = round(sum(SP_PINE)/n(), digits = 2),
                SNAG_INT = round(sum(SNAG_INT)/n(), digits = 2),
                SNAG_DEC = round(sum(SNAG_DEC)/n(), digits = 2),
                BRKN_before = round(sum(BRKN_before)/n(), digits = 2),
                BRKN_after = round(sum(BRKN_after)/n(), digits = 2),
                BRKN = round(sum(BRKN)/n(), digits = 2),
                SnagDens_23to38 = paste0(round(mean(SnagDens_23to38), digits = 1), " (",
                                         round(sd(SnagDens_23to38), digits = 1), ")"),
                SnagDens_38to50 = paste0(round(mean(SnagDens_38to50), digits = 1), " (",
                                         round(sd(SnagDens_38to50), digits = 1), ")"),
                SnagDens_gt50 = paste0(round(mean(SnagDens_gt50), digits = 1), " (",
                                       round(sd(SnagDens_gt50), digits = 1), ")"),
                SnagDens_pine = paste0(round(mean(SnagDens_pine), digits = 1), " (",
                                       round(sd(SnagDens_pine), digits = 1), ")"),
                SnagDens_fir = paste0(round(mean(SnagDens_fir), digits = 1), " (",
                                      round(sd(SnagDens_fir), digits = 1), ")"),
                SnagDens_all = paste0(round(mean(SnagDens_all), digits = 1), " (",
                                      round(sd(SnagDens_all), digits = 1), ")"),
                slope = paste0(round(mean(slope), digits = 1), " (", round(sd(slope), digits = 1), ")"),
                sinasp = paste0(round(mean(sinasp), digits = 1), " (", round(sd(sinasp), digits = 1), ")"),
                cosasp = paste0(round(mean(cosasp), digits = 1), " (", round(sd(cosasp), digits = 1), ")"),
                ccmort_loc = paste0(round(mean(ccmort_loc), digits = 1), " (",
                                    round(sd(ccmort_loc), digits = 1), ")"),
                blk_lndcc = paste0(round(mean(blk_lndcc), digits = 1), " (",
                                   round(sd(blk_lndcc), digits = 1), ")"),
                canlo_loc = paste0(round(mean(canlo_loc), digits = 1), " (",
                                   round(sd(canlo_loc), digits = 1), ")"),
                canlo_lnd = paste0(round(mean(canlo_lnd), digits = 1), " (",
                                   round(sd(canlo_lnd), digits = 1), ")"),
                canmd_loc = paste0(round(mean(canmd_loc), digits = 1), " (",
                                   round(sd(canmd_loc), digits = 1), ")"),
                canmd_lnd = paste0(round(mean(canmd_lnd), digits = 1), " (",
                                   round(sd(canmd_lnd), digits = 1), ")"),
                canhi_loc = paste0(round(mean(canhi_loc), digits = 1), " (",
                                   round(sd(canhi_loc), digits = 1), ")"),
                canhi_lnd = paste0(round(mean(canhi_lnd), digits = 1), " (",
                                   round(sd(canhi_lnd), digits = 1), ")"),
                sizsm_loc = paste0(round(mean(sizsm_loc), digits = 1), " (",
                                   round(sd(sizsm_loc), digits = 1), ")"),
                sizlrg_loc = paste0(round(mean(sizlrg_loc), digits = 1), " (",
                                    round(sd(sizlrg_loc), digits = 1), ")"),
                sizlrg_lnd = paste0(round(mean(sizlrg_lnd), digits = 1), " (",
                                    round(sd(sizlrg_lnd), digits = 1), ")"),
                fir_1km = paste0(round(mean(fir_1km), digits = 1), " (",
                                 round(sd(fir_1km), digits = 1), ")"))
    
  }
  else{
    dta <- dta %>% 
      group_by(TYPE) %>%
      summarise(DBH = paste0(round(mean(DBH), digits = 1), " (", round(sd(DBH), digits = 1), ")"),
                SP_FIR = round(sum(SP_FIR)/n(), digits = 2),
                SP_PINE = round(sum(SP_PINE)/n(), digits = 2),
                SnagDens_23to38 = paste0(round(mean(SnagDens_23to38), digits = 1), " (",
                                         round(sd(SnagDens_23to38), digits = 1), ")"),
                SnagDens_38to50 = paste0(round(mean(SnagDens_38to50), digits = 1), " (",
                                         round(sd(SnagDens_38to50), digits = 1), ")"),
                SnagDens_gt50 = paste0(round(mean(SnagDens_gt50), digits = 1), " (",
                                       round(sd(SnagDens_gt50), digits = 1), ")"),
                SnagDens_pine = paste0(round(mean(SnagDens_pine), digits = 1), " (",
                                       round(sd(SnagDens_pine), digits = 1), ")"),
                SnagDens_fir = paste0(round(mean(SnagDens_fir), digits = 1), " (",
                                      round(sd(SnagDens_fir), digits = 1), ")"),
                SnagDens_all = paste0(round(mean(SnagDens_all), digits = 1), " (",
                                      round(sd(SnagDens_all), digits = 1), ")"),
                slope = paste0(round(mean(slope), digits = 1), " (", round(sd(slope), digits = 1), ")"),
                sinasp = paste0(round(mean(sinasp), digits = 1), " (", round(sd(sinasp), digits = 1), ")"),
                cosasp = paste0(round(mean(cosasp), digits = 1), " (", round(sd(cosasp), digits = 1), ")"),
                ccmort_loc = paste0(round(mean(ccmort_loc), digits = 1), " (",
                                    round(sd(ccmort_loc), digits = 1), ")"),
                blk_lndcc = paste0(round(mean(blk_lndcc), digits = 1), " (",
                                   round(sd(blk_lndcc), digits = 1), ")"),
                canlo_loc = paste0(round(mean(canlo_loc), digits = 1), " (",
                                   round(sd(canlo_loc), digits = 1), ")"),
                canlo_lnd = paste0(round(mean(canlo_lnd), digits = 1), " (",
                                   round(sd(canlo_lnd), digits = 1), ")"),
                canmd_loc = paste0(round(mean(canmd_loc), digits = 1), " (",
                                   round(sd(canmd_loc), digits = 1), ")"),
                canmd_lnd = paste0(round(mean(canmd_lnd), digits = 1), " (",
                                   round(sd(canmd_lnd), digits = 1), ")"),
                canhi_loc = paste0(round(mean(canhi_loc), digits = 1), " (",
                                   round(sd(canhi_loc), digits = 1), ")"),
                canhi_lnd = paste0(round(mean(canhi_lnd), digits = 1), " (",
                                   round(sd(canhi_lnd), digits = 1), ")"),
                sizsm_loc = paste0(round(mean(sizsm_loc), digits = 1), " (",
                                   round(sd(sizsm_loc), digits = 1), ")"),
                sizlrg_loc = paste0(round(mean(sizlrg_loc), digits = 1), " (",
                                    round(sd(sizlrg_loc), digits = 1), ")"),
                sizlrg_lnd = paste0(round(mean(sizlrg_lnd), digits = 1), " (",
                                    round(sd(sizlrg_lnd), digits = 1), ")"),
                fir_1km = paste0(round(mean(fir_1km), digits = 1), " (",
                                 round(sd(fir_1km), digits = 1), ")"))
    
  }
  dta <- dta %>% mutate(TYPE = paste0(sp, "_", TYPE))
  assign(paste0("tabsum.", sp), dta)
}

tab.sum <- bind_rows(tabsum.BBWO, tabsum.HAWO, tabsum.WHWO, tabsum.NOFL, tabsum.RBSA, tabsum.MOBL)
write.csv(tab.sum, "Predictor_summary_stats.csv", row.names = F)
