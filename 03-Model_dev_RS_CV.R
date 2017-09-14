setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")
#############################################################################################
# Develops models with remotely sensed data only for 4 species with sufficient sample sizes #
# for 3-fold cross-fire validation                                                          #
#############################################################################################

load("Data_compiled.RData")

library(R.utils)
library(gtools)
library(dplyr)
library(PresenceAbsence)
library(stringr)

### Load functions ###
source("scripts/Functions.R")

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
#rm(n, s)

unit.size <- list(BBWO = c(10, 32, 6), HAWO = c(9, 28, 5), WHWO = c(11, 30, 7),
                  NOFL = c(10, 30, 7)) # size of components of bins (BBWO)
bin.width <- list(BBWO = c(3, 3, 3), HAWO = c(3, 3, 4), WHWO = c(3, 3, 4),
                  NOFL = c(3, 3, 3)) # in number of units
sprmn.crit <- 0.564 # Critical value at alpha = 0.05 for Spearman correlation at n = 10

# Remote-sensed variables only #
vars.rs <- names(dat.BBWO)[c(25:28, 30, 35:41)] # For canopy mortality severity metrics

# Fit models and calculate evaluation metrics
RSmods <- AllFit(spp, sites, vars = vars.rs)

## Review model lists and save selected model for each species ##

## Save model output tables for each species ##
write.csv(RSmods$BBWO, "Models_for_review_BBWO.csv", row.names = F)
write.csv(RSmods.HAWO, "Models_for_review_HAWO.csv", row.names = F)
write.csv(RSmods.HAWO, "Models_for_review_WHWO.csv", row.names = F)
write.csv(RSmods.HAWO, "Models_for_review_NOFL.csv", row.names = F)

# For each species, make sure selected model's coefficients are consistent in direction across folds.
#dat.fold <- dat.BBWO[which(dat.BBWO$Site!=sites[3]),]
#mod <- WLR_fit(dat.fold, formula = Nest ~ ccmort_loc + blk_lndcc + canhi_loc)
#summary(mod)
#rm(dat.fold)

  # Save selected models #
mod <- WLR_fit(dat.BBWO, formula = Nest ~ ccmort_loc + blk_lndcc + canhi_loc)
saveObject(mod, "Model_RS_BBWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_RS_BBWO.csv")

mod <- WLR_fit(dat.HAWO, formula = Nest ~ ccmort_loc + canhi_lnd + sizlrg_loc)
saveObject(mod, "Model_RS_HAWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_RS_HAWO.csv")

mod <- WLR_fit(dat.WHWO, formula = Nest ~ ccmort_loc + blk_lndcc)
saveObject(mod, "Model_RS_WHWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_RS_WHWO.csv")

mod <- WLR_fit(dat.NOFL, formula = Nest ~ slope + ccmort_loc + sizlrg_loc)
saveObject(mod, "Model_RS_NOFL") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_RS_NOFL.csv")

rm(mod)

save.image("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/Wrkspc_RSmods_CV.RData")
