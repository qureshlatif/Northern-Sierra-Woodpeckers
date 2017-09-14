###################################################################################################
# Purpose: Initial data compilation for Northern Sierra post-fire woodpecker habitat modeling.    #

# Content:                                                                                        #
# 1. Imports and merges tree-level and patch-level field-collected data, and remotely sensed data #
# for nest and random sites.                                                                      #
# 2. Compiles variables considered for weighted logistic regression habitat models.               #
# 3. Saves workspace containing resulting data tables                                             #
###################################################################################################

library(dplyr)
library(stringr)
library(foreign)

setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")

## Nest tree data ##
dat <- read.csv("random_and_nest_trees_for_analysis_20170531.csv",
                header = T, stringsAsFactors = F) %>% tbl_df %>%
  select(SAMPLE_ID, EASTING, NORTHING, YEAR, SPECIES, TYPE, TREATMENT, TREE_SPECIES, DECAY, DBH, TOP) %>%
  filter(!is.na(DBH)) %>% # Get rid of records with missing DBH (16 randoms and 1 NOFL)
  mutate(SP_FIR = 1*(TREE_SPECIES %in% c("ABLO", "ABMA", "PSME", "ABIES", "FIRSP"))) %>%
  mutate(SP_FIR = SP_FIR %>% as.integer) %>%
  mutate(SP_PINE = 1*(TREE_SPECIES %in% c("PIPO", "PILA", "PICO", "PINUS", "PIYE", "PIJE"))) %>%
  mutate(SP_PINE = SP_PINE %>% as.integer) %>%
  select(-TREE_SPECIES) %>%
  mutate(SNAG = 1*(DECAY %in% c(3, 4, 5, 6, 7))) %>%
  mutate(SNAG = SNAG %>% as.integer) %>%
  mutate(SNAG_INT = 1*(DECAY %in% c(3, 4))) %>%
  mutate(SNAG_INT = SNAG_INT %>% as.integer) %>%
  mutate(SNAG_DEC = 1*(DECAY %in% c(5, 6, 7))) %>%
  mutate(SNAG_DEC = SNAG_DEC %>% as.integer)
dat <- dat %>%
  mutate(BRKN_before = 1*(TOP %in% c("BB", "F BB"))) %>%
  mutate(BRKN_before = replace(BRKN_before, which(dat$TOP %in% c("B", "BB") & dat$DECAY %in% c(4, 5)), 1)) %>%
  mutate(BRKN_before = BRKN_before %>% as.integer) %>%
  mutate(BRKN_after = 1*(TOP %in% c("BA", "F BA"))) %>%
  mutate(BRKN_after = replace(BRKN_after, which(dat$TOP %in% c("B", "BB") & dat$DECAY == 3), 1)) %>%
  mutate(BRKN_after = BRKN_after %>% as.integer) %>%
  mutate(BRKN = (BRKN_before == 1 | BRKN_after == 1)*1) %>%
  select(-DECAY) %>%
  select(-TOP)

dat <- dat[-which(duplicated(dat$SAMPLE_ID)),] # Remove duplicates

## Patch (snag density) variables ##
dat.plot <- read.csv("snag_plots_for_analysis_20170531.csv",
                     header = T, stringsAsFactors = F) %>% tbl_df %>%
  select(SAMPLE_ID, TREE_SPECIES, DECAY, DBH) %>%
  filter(!DECAY == 4) %>% #remove live trees
  group_by(SAMPLE_ID) %>%
  summarise(SnagDens_23to38 = sum(DBH < 38),
            SnagDens_38to50 = sum(DBH >= 38 & DBH < 50),
            SnagDens_gt50 = sum(DBH >= 50),
            SnagDens_23to50 = sum(DBH < 50),
            SnagDens_gt38 = sum(DBH >= 38),
            SnagDens_all = sum(DBH>=23),
            SnagDens_pine = sum(TREE_SPECIES %in% c("PIPO", "PILA", "PICO", "PINUS", "PIYE", "PIJE")),
            SnagDens_fir = sum(TREE_SPECIES %in% c("ABLO", "ABMA", "PSME", "ABIES", "FIRSP")))

dat <- dat %>% left_join(dat.plot, by = "SAMPLE_ID") %>%
  mutate(SnagDens_23to38 = replace(SnagDens_23to38, which(is.na(SnagDens_23to38)), 0)) %>%
  mutate(SnagDens_38to50 = replace(SnagDens_38to50, which(is.na(SnagDens_38to50)), 0)) %>%
  mutate(SnagDens_gt50 = replace(SnagDens_gt50, which(is.na(SnagDens_gt50)), 0)) %>%
  mutate(SnagDens_23to50 = replace(SnagDens_23to50, which(is.na(SnagDens_23to50)), 0)) %>%
  mutate(SnagDens_gt38 = replace(SnagDens_gt38, which(is.na(SnagDens_gt38)), 0)) %>%
  mutate(SnagDens_all = replace(SnagDens_all, which(is.na(SnagDens_all)), 0)) %>%
  mutate(SnagDens_pine = replace(SnagDens_pine, which(is.na(SnagDens_pine)), 0)) %>%
  mutate(SnagDens_fir = replace(SnagDens_fir, which(is.na(SnagDens_fir)), 0))
rm(dat.plot)

## Remotely sensed variables ##
dat.remote <- read.dbf("E:/GISData/PtBlue_Sierra/NR_points.dbf", as.is = T) %>% tbl_df %>%
  select(SAMPLE_ID, slope:pine_1km) %>%
  mutate(blk_lndcc = blk_lndcc*100) %>%
  mutate(blk_lnddnb = blk_lnddnb*100) %>%
  mutate(grn_lndcc = grn_lndcc*100) %>%
  mutate(grn_lnddnb = grn_lnddnb*100) %>%
  mutate(canlo_loc = 100*canlo_loc) %>%
  mutate(canmd_loc = 100*canmd_loc) %>%
  mutate(canhi_loc = 100*canhi_loc) %>%
  mutate(canlo_lnd = 100*canlo_lnd) %>%
  mutate(canmd_lnd = 100*canmd_lnd) %>%
  mutate(canhi_lnd = 100*canhi_lnd) %>%
  mutate(sizpl_loc = 100*sizpl_loc) %>%
  mutate(sizsm_loc = 100*sizsm_loc) %>%
  mutate(sizlrg_loc = 100*sizlrg_loc) %>%
  mutate(sizpl_lnd = 100*sizpl_lnd) %>%
  rename(sizsm_lnd = sizesm_lnd) %>%
  mutate(sizsm_lnd = 100*sizsm_lnd) %>%
  mutate(sizlrg_lnd = 100*sizlrg_lnd) %>%
  mutate(fir_1km = 100*fir_1km) %>%
  mutate(pine_1km = 100*pine_1km)
dat <- dat %>% left_join(dat.remote, by = "SAMPLE_ID")
rm(dat.remote)

# Check sum of remotely sensed tree size variables
#hist(dat$sizpl_loc + dat$sizsm_loc + dat$sizlrg_loc)
#hist(dat$sizpl_lnd + dat$sizsm_lnd + dat$sizlrg_lnd)

# Correlation matrix (all species)
# write.csv(cor(dat[,-c(1:5)], use = "complete.obs"), "CorMat_predictors.csv")

# Add time since fire #
dat <- dat %>%
  mutate(Site = substr(SAMPLE_ID, 1, 2)) %>%
  mutate(fire_year = 2008)
dat <- dat %>%
  mutate(fire_year = replace(fire_year, which(dat$Site == "CH"), 2012)) %>%
  mutate(fire_year = replace(fire_year, which(dat$Site == "ML"), 2007)) %>%
  mutate(TimSincFire = YEAR - fire_year) %>%
  select(-fire_year)

save.image("Data_compiled0.RData")
