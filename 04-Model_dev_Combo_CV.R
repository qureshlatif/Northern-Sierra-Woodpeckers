setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")
#################################################################################################
# Purpose: Develop models with field-collected and remotely sensed variables for 4 species with #
# sufficient sample sizes for 3-fold cross-fire validation.                                     #                                                         #
#################################################################################################

load("Data_compiled.RData")

library(R.utils)
library(gtools)
library(dplyr)
library(PresenceAbsence)
library(stringr)
sites <- unique(dat.BBWO$Site)

### Load functions ###
source("Northern-Sierra-Woodpeckers/Functions.R")

# Species with sufficient sample sizes for cross-fire validation #
sites <- c("CB", "ML", "CH")
spp <- c("BBWO", "HAWO", "WHWO", "NOFL")

## Set bin parameters for calculating RPIs ##
# Review sample sizes by site for each species #
n.val <- getNbySpp(spp, sites, datList(spp))

# Try different unit sizes and bin widths to get n = 10 bins for each species
#n.val$NOFL # Review sample sizes for evaluation datasets and set bin parameters accordingly:
#s <- 3
#n <- n.val$NOFL
#unit.size <- c(10, 30, 7) # size of components of bins (HAWO)
#bin.width <- c(3, 3, 3) # in number of units
#unitID <- rep(1:ceiling(n[s]/unit.size[s]), each = unit.size[s])[1:n[s]]
#bin.st <- 1:(max(unitID) - bin.width[s] + 1)
#length(bin.st) # Targe value = 10
#rm(n)

unit.size <- list(BBWO = c(10, 32, 6), HAWO = c(9, 28, 5), WHWO = c(11, 30, 7),
                  NOFL = c(10, 30, 7)) # size of components of bins (BBWO)
bin.width <- list(BBWO = c(3, 3, 3), HAWO = c(3, 3, 4), WHWO = c(3, 3, 4),
                  NOFL = c(3, 3, 3)) # in number of units
sprmn.crit <- 0.564 # Critical value at alpha = 0.05 for Spearman correlation at n = 10

# List variables #
vars.rs <- list(BBWO = c("ccmort_loc", "blk_lndcc", "canhi_loc"),
                HAWO = c("ccmort_loc", "canhi_lnd", "sizlrg_loc"),
                WHWO = c("ccmort_loc", "blk_lndcc"),
                NOFL = c("slope", "ccmort_loc", "sizlrg_loc"))
vars.fc <- list(BBWO = c("DBH", "DBH + I(DBH^2)", "SnagDens_23to50"),
                HAWO = c("DBH", "DBH + I(DBH^2)", "SP_PINE", "BRKN",
                         "SnagDens_all", "SP_PINE + TimSincFire + I(SP_PINE*TimSincFire)"),
                WHWO = c("DBH", "DBH + I(DBH^2)", "BRKN", "SnagDens_gt50"),
                NOFL = c("DBH", "DBH + I(DBH^2)", "SP_FIR", "BRKN", "SnagDens_gt50"))
vars <- list()
for(sp in 1:length(spp)) vars[[sp]] <- c(vars.rs[[sp]], vars.fc[[sp]])
names(vars) <- spp
rm(sp)

# Fit models and calculate evaluation metrics
CMBmods <- AllFit(spp, sites, vars = vars)

## Review model lists and save selected model for each species ##

## Save model output tables for each species ##
write.csv(CMBmods$BBWO, "Models_for_review_BBWO.csv", row.names = F)
write.csv(CMBmods.HAWO, "Models_for_review_HAWO.csv", row.names = F)
write.csv(CMBmods.HAWO, "Models_for_review_WHWO.csv", row.names = F)
write.csv(CMBmods.HAWO, "Models_for_review_NOFL.csv", row.names = F)

# For each species, make sure selected model's coefficients are consistent in direction across folds.
#dat.fold <- dat.BBWO[which(dat.BBWO$Site!=sites[3]),]
#mod <- WLR_fit(dat.fold, formula = Nest ~ ccmort_loc + blk_lndcc + canhi_loc)
#summary(mod)
#rm(dat.fold)

# Save selected models #
mod <- WLR_fit(dat.BBWO, formula = Nest ~ ccmort_loc + blk_lndcc + canhi_loc + DBH + I(DBH^2) + SnagDens_23to50)
saveObject(mod, "Model_CMB_BBWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_CMB_BBWO.csv")

mod <- WLR_fit(dat.HAWO, formula = Nest ~ ccmort_loc + canhi_lnd + sizlrg_loc + BRKN +
                 SnagDens_all + SP_PINE + TimSincFire + I(SP_PINE*TimSincFire))
saveObject(mod, "Model_CMB_HAWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_CMB_HAWO.csv")

mod <- WLR_fit(dat.WHWO, formula = Nest ~ ccmort_loc + blk_lndcc + DBH + I(DBH^2) + BRKN)
saveObject(mod, "Model_CMB_WHWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_CMB_WHWO.csv")

mod <- WLR_fit(dat.NOFL, formula = Nest ~ slope + DBH + I(DBH^2) + BRKN)
saveObject(mod, "Model_CMB_NOFL") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_CMB_NOFL.csv")

rm(mod)

save.image("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/Wrkspc_CMBmods_CV.RData")
