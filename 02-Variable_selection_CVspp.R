setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")

####################################################################################################################
# Purpose: fit models within variable groups for each species to inform screening prior to full model construction #
####################################################################################################################
# Content:                                                                                                         #
#                                                                                                                  #
# Field-collected variable groups:                                                                                 #
#  Nest tree variables: Tree species, decay class (intact or decayed snag), DBH (linear or quadratic),             #
#      top condition (intact, broken before, or broken after)                                                      #
#  Nest patch variables: Snag density by size class (23-38cm, 38-50cm, >50cm, 23-50cm, >38cm),                     #
#      Snag density by species (pine, fir)                                                                         #
#                                                                                                                  #
# Remotely sensed variable groups:                                                                                 #
#  Topography: slope, aspect (cosine + sine)                                                                       #
#  Burn severity (MTBS): local (3x3 cell) median RdNBR-derived % canopy mortality,                                 #
#     landscape (1km radius) proportion high severity (canopy mort > 64%)                                          #
#  Canopy cover (CWHR): local- or landscape scale % low (<25%), % moderate (25-40%), and/or % high (>40%)          #
#  Tree size dominance (CWHR): local- or landscape-scale % small (11-24") or % mid-large (>24") dominated forest   #
#  Forest type (CWHR): landscape proportion fir-dominated forest (dropped pine-dominated forest due to limited     #
#     variation)                                                                                                   #
#                                                                                                                  #
# Fits models representing all combinations within groups and generates AIC tables for review                      #
####################################################################################################################

load("Data_compiled.RData")

library(R.utils)
library(gtools)
library(dplyr)
library(PresenceAbsence)
library(stringr)
sites <- unique(dat.BBWO$Site)

### Load functions ###
source("scripts/Functions.R")

spp <- c("BBWO", "HAWO", "WHWO", "NOFL", "RBSA", "MOBL")

## Remote-sensed variables ##

# Topography #
vars <- c("slope", "sinasp", "cosasp") # For canopy mortality severity metrics
Group_label <- "Topo"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

# Burn severity #
vars <- c("ccmort_loc", "blk_lndcc") # For canopy mortality severity metrics
Group_label <- "BSev"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

# Canopy cover #
vars <- c("canlo_loc", "canlo_lnd", "canmd_loc", "canmd_lnd", "canhi_loc", "canhi_lnd")
Group_label <- "CCov"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

# Tree size #
vars <- c("sizsm_loc", "sizlrg_loc", "sizlrg_lnd")
Group_label <- "TrSize"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

# Forest type #
vars <- c("fir_1km")
Group_label <- "Fortyp"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

## Field-collected variables for excavators ##
spp <- c("BBWO", "HAWO", "WHWO", "NOFL")

# Nest tree #
vars <- c("DBH", "DBH + I(DBH^2)", "SP_FIR", "SP_PINE", "SNAG_DEC",
          "BRKN_before", "BRKN_after", "BRKN",
          "BRKN_after + TimSincFire + I(BRKN_after*TimSincFire)",
          "BRKN_before + TimSincFire + I(BRKN_before*TimSincFire)",
          "SP_FIR + TimSincFire + I(SP_FIR*TimSincFire)",
          "SP_PINE + TimSincFire + I(SP_PINE*TimSincFire)")
Group_label <- "Tree"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

# Patch structure #
vars <- c("SnagDens_23to38", "SnagDens_38to50", "SnagDens_gt50", "SnagDens_23to50",
          "SnagDens_gt38", "SnagDens_all", "SnagDens_pine", "SnagDens_fir")
Group_label <- "Patch"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

## Field-collected variables for excavators ##
spp <- c("RBSA", "MOBL")

# Nest tree #
vars <- c("DBH", "DBH + I(DBH^2)", "SP_FIR", "SP_PINE", 
          "SP_FIR + TimSincFire + I(SP_FIR*TimSincFire)",
          "SP_PINE + TimSincFire + I(SP_PINE*TimSincFire)")
Group_label <- "Tree"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

# Patch structure #
vars <- c("SnagDens_23to38", "SnagDens_38to50", "SnagDens_gt50", "SnagDens_23to50",
          "SnagDens_gt38", "SnagDens_all", "SnagDens_pine", "SnagDens_fir")
Group_label <- "Patch"

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  maxK <- min(sum(dat.spp$Nest)*2 - 1, length(vars))
  out <- AIC_table(dat.spp, vars, maxK)
  assign(paste0(Group_label, ".", spp[sp]), out)
}

## Combine tables by species and save ##
out.BBWO <- bind_rows(Topo.BBWO, BSev.BBWO, CCov.BBWO, TrSize.BBWO, Fortyp.BBWO,
                      Tree.BBWO, Patch.BBWO)
write.csv(out.BBWO, "Variable_selection_models_BBWO.csv", row.names = F)

out.HAWO <- bind_rows(Topo.HAWO, BSev.HAWO, CCov.HAWO, TrSize.HAWO, Fortyp.HAWO,
                      Tree.HAWO, Patch.HAWO)
write.csv(out.HAWO, "Variable_selection_models_HAWO.csv", row.names = F)

out.WHWO <- bind_rows(Topo.WHWO, BSev.WHWO, CCov.WHWO, TrSize.WHWO, Fortyp.WHWO,
                      Tree.WHWO, Patch.WHWO)
write.csv(out.WHWO, "Variable_selection_models_WHWO.csv", row.names = F)

out.NOFL <- bind_rows(Topo.NOFL, BSev.NOFL, CCov.NOFL, TrSize.NOFL, Fortyp.NOFL,
                      Tree.NOFL, Patch.NOFL)
write.csv(out.NOFL, "Variable_selection_models_NOFL.csv", row.names = F)

out.RBSA <- bind_rows(Topo.RBSA, BSev.RBSA, CCov.RBSA, TrSize.RBSA, Fortyp.RBSA,
                      Tree.RBSA, Patch.RBSA)
write.csv(out.RBSA, "Variable_selection_models_RBSA.csv", row.names = F)

out.MOBL <- bind_rows(Topo.MOBL, BSev.MOBL, CCov.MOBL, TrSize.MOBL, Fortyp.MOBL,
                      Tree.MOBL, Patch.MOBL)
write.csv(out.MOBL, "Variable_selection_models_MOBL.csv", row.names = F)
