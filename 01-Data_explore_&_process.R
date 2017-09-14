######################################################################################################
# Purpose: Some data exploration and compilation of species-specific datasets                        #
######################################################################################################
# Content:                                                                                           #
# 1. Scatter plots and correlation coefficients used to screen variables to avoid multi-collinearity #
# 2. Screen non-nest points for min distance to nest of 130m,                                        #
# 3. Screen non-nest points centered on live trees to match snag:live ratios for nests               #
# 4. Scale (z-score) variables to mean = 0, SD = 1                                                   #
# 5. Compile species-specific datasets and save workspace containing these                           #
######################################################################################################

setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")
load("Data_compiled0.RData")

library(dplyr)
library(gcookbook)

# Trim variables based on correlations
dat <- dat %>%
  select(-grn_lndcc) %>% select(-grn_lnddnb) %>%
  select(-blk_lnddnb) %>%
  select(-sizpl_loc) %>% select(-sizpl_lnd) %>%
  select(-pine_1km)

## All species scatter plots of remaining variables ##
#RSdat <- dat %>% select(slope:fir_1km)
#pairs(RSdat)

#FCdat <- dat %>% select(DBH:SnagDens_fir)
#pairs(FCdat %>% select(DBH, SnagDens_23to38:SnagDens_fir))

#rm(FCdat, RSdat)

# Correlation matrix (all species)
#write.csv(cor(dat[,-c(1:7, 37:38)], use = "complete.obs"), "CorMat_predictors.csv")

## Function for finding nearest nest distance for random points ##
library(flexclust)
nearNestDist <- function(dat) {
  nests <- dat[which(dat$TYPE == "NEST"), c("EASTING", "NORTHING")]
  d <- dist2(dat[, c("EASTING", "NORTHING")], nests)
  d <- apply(d, 1, min)
  return(d)
}
minDist <- 130 # Minimum distance (m) between nest and random points for a given species

## z-score all continuous predictors and store scaling factors ##
vars.cont <- names(dat)[c(8, 17:41, 43)]
scale.factors <- matrix(NA, nrow = length(vars.cont), ncol = 2,
                        dimnames = list(vars.cont, c("mean", "SD")))
datz <- dat
for(i in 1:length(vars.cont)) {
  scale.factors[i, "mean"] <- dat[which(dat$TYPE == "RAND"), vars.cont[i]] %>% as.matrix %>% mean
  scale.factors[i, "SD"] <- dat[which(dat$TYPE == "RAND"), vars.cont[i]] %>% as.matrix %>% sd
  datz[, vars.cont[i]] <- (dat[, vars.cont[i]] - scale.factors[i, "mean"]) / scale.factors[i, "SD"]
}

# Species-specific datasets
# BBWO
dat.BBWO <- datz %>% filter(TYPE == "RAND" | SPECIES == "BBWO") %>%
  select(SAMPLE_ID, EASTING, NORTHING, YEAR, TYPE, Site, TimSincFire, DBH:fir_1km) %>%
  mutate(Nest = (TYPE == "NEST")*1)
dat.BBWO <- dat.BBWO %>% mutate(dist = nearNestDist(dat.BBWO)) %>%
  filter(!(TYPE == "RAND" & dist < 130)) %>%
  select(-dist)
  # filter out non-snag randoms to equalize with nests
SnagNest <- sum(dat.BBWO$TYPE == "NEST" & dat.BBWO$SNAG == 1)
LiveNest <- sum(dat.BBWO$TYPE == "NEST" & dat.BBWO$SNAG == 0)
SnagRand <- sum(dat.BBWO$TYPE == "RAND" & dat.BBWO$SNAG == 1)
LiveRand <- round(SnagRand * (LiveNest / SnagNest))
keep <- c(which(dat.BBWO$TYPE == "NEST"),
    which(dat.BBWO$TYPE == "RAND" & dat.BBWO$SNAG == 1),
    sample(which(dat.BBWO$TYPE == "RAND" & dat.BBWO$SNAG == 0), LiveRand))
dat.BBWO <- dat.BBWO %>% slice(keep)

#RSdat <- dat.BBWO %>% select(Nest, slope:fir_1km)
#pairs(RSdat) # Look for separation (none found)

#FCdat <- dat.BBWO %>% select(Nest, DBH:SnagDens_fir)
#pairs(FCdat) # Look for separation (none found)

# HAWO
dat.HAWO <- datz %>% filter(TYPE == "RAND" | SPECIES == "HAWO") %>%
  select(SAMPLE_ID, EASTING, NORTHING, YEAR, TYPE, Site, TimSincFire, DBH:fir_1km) %>%
  mutate(Nest = (TYPE == "NEST")*1)
dat.HAWO <- dat.HAWO %>% mutate(dist = nearNestDist(dat.HAWO)) %>%
  filter(!(TYPE == "RAND" & dist < 130)) %>%
  select(-dist)
  # filter out non-snag randoms to equalize with nests
SnagNest <- sum(dat.HAWO$TYPE == "NEST" & dat.HAWO$SNAG == 1)
LiveNest <- sum(dat.HAWO$TYPE == "NEST" & dat.HAWO$SNAG == 0)
SnagRand <- sum(dat.HAWO$TYPE == "RAND" & dat.HAWO$SNAG == 1)
LiveRand <- round(SnagRand * (LiveNest / SnagNest))
keep <- c(which(dat.HAWO$TYPE == "NEST"),
          which(dat.HAWO$TYPE == "RAND" & dat.HAWO$SNAG == 1),
          sample(which(dat.HAWO$TYPE == "RAND" & dat.HAWO$SNAG == 0), LiveRand))
dat.HAWO <- dat.HAWO %>% slice(keep)

#RSdat <- dat.HAWO %>% select(Nest, slope:fir_1km)
#pairs(RSdat) # Look for separation (none found)

#FCdat <- dat.HAWO %>% select(Nest, DBH:SnagDens_fir)
#pairs(FCdat) # Look for separation (none found)

# WHWO
dat.WHWO <- datz %>% filter(TYPE == "RAND" | SPECIES == "WHWO") %>%
  select(SAMPLE_ID, EASTING, NORTHING, YEAR, TYPE, Site, TimSincFire, DBH:fir_1km) %>%
  mutate(Nest = (TYPE == "NEST")*1)
dat.WHWO <- dat.WHWO %>% mutate(dist = nearNestDist(dat.WHWO)) %>%
  filter(!(TYPE == "RAND" & dist < 130)) %>%
  select(-dist)
  # filter out non-snag randoms to equalize with nests
SnagNest <- sum(dat.WHWO$TYPE == "NEST" & dat.WHWO$SNAG == 1)
LiveNest <- sum(dat.WHWO$TYPE == "NEST" & dat.WHWO$SNAG == 0)
SnagRand <- sum(dat.WHWO$TYPE == "RAND" & dat.WHWO$SNAG == 1)
LiveRand <- round(SnagRand * (LiveNest / SnagNest))
keep <- c(which(dat.WHWO$TYPE == "NEST"),
          which(dat.WHWO$TYPE == "RAND" & dat.WHWO$SNAG == 1),
          sample(which(dat.WHWO$TYPE == "RAND" & dat.WHWO$SNAG == 0), LiveRand))
dat.WHWO <- dat.WHWO %>% slice(keep)

#RSdat <- dat.WHWO %>% select(Nest, slope:fir_1km)
#pairs(RSdat) # Look for separation (none found)

#FCdat <- dat.WHWO %>% select(Nest, DBH:SnagDens_fir)
#pairs(FCdat) # Look for separation (none found)

# NOFL
dat.NOFL <- datz %>% filter(TYPE == "RAND" | SPECIES == "NOFL") %>%
  select(SAMPLE_ID, EASTING, NORTHING, YEAR, TYPE, Site, TimSincFire, DBH:fir_1km) %>%
  mutate(Nest = (TYPE == "NEST")*1)
dat.NOFL <- dat.NOFL %>% mutate(dist = nearNestDist(dat.NOFL)) %>%
  filter(!(TYPE == "RAND" & dist < 130)) %>%
  select(-dist)
  # filter out non-snag randoms to equalize with nests
SnagNest <- sum(dat.NOFL$TYPE == "NEST" & dat.NOFL$SNAG == 1)
LiveNest <- sum(dat.NOFL$TYPE == "NEST" & dat.NOFL$SNAG == 0)
SnagRand <- sum(dat.NOFL$TYPE == "RAND" & dat.NOFL$SNAG == 1)
LiveRand <- round(SnagRand * (LiveNest / SnagNest))
keep <- c(which(dat.NOFL$TYPE == "NEST"),
          which(dat.NOFL$TYPE == "RAND" & dat.NOFL$SNAG == 1),
          sample(which(dat.NOFL$TYPE == "RAND" & dat.NOFL$SNAG == 0), LiveRand))
dat.NOFL <- dat.NOFL %>% slice(keep)

#RSdat <- dat.NOFL %>% select(Nest, slope:fir_1km)
#pairs(RSdat) # Look for separation (none found)

#FCdat <- dat.NOFL %>% select(Nest, DBH:SnagDens_fir)
#pairs(FCdat) # Look for separation (none found)

# RBSA
dat.RBSA <- datz %>% filter(TYPE == "RAND" | SPECIES == "RBSA") %>%
  select(SAMPLE_ID, EASTING, NORTHING, YEAR, TYPE, Site, TimSincFire, DBH:SNAG, SnagDens_23to38:fir_1km) %>%
  mutate(Nest = (TYPE == "NEST")*1)
dat.RBSA <- dat.RBSA %>% mutate(dist = nearNestDist(dat.RBSA)) %>%
  filter(!(TYPE == "RAND" & dist < 130)) %>%
  select(-dist)
  # filter out non-snag randoms to equalize with nests
SnagNest <- sum(dat.RBSA$TYPE == "NEST" & dat.RBSA$SNAG == 1)
LiveNest <- sum(dat.RBSA$TYPE == "NEST" & dat.RBSA$SNAG == 0)
SnagRand <- sum(dat.RBSA$TYPE == "RAND" & dat.RBSA$SNAG == 1)
LiveRand <- round(SnagRand * (LiveNest / SnagNest))
keep <- c(which(dat.RBSA$TYPE == "NEST"),
          which(dat.RBSA$TYPE == "RAND" & dat.RBSA$SNAG == 1),
          sample(which(dat.RBSA$TYPE == "RAND" & dat.RBSA$SNAG == 0), LiveRand))
dat.RBSA <- dat.RBSA %>% slice(keep)

#RSdat <- dat.RBSA %>% select(Nest, slope:fir_1km)
#pairs(RSdat) # Look for separation (none found)

#FCdat <- dat.RBSA %>% select(Nest, DBH:SnagDens_fir)
#pairs(FCdat) # Look for separation (none found)

# MOBL
dat.MOBL <- datz %>% filter(TYPE == "RAND" | SPECIES == "MOBL") %>%
  select(SAMPLE_ID, EASTING, NORTHING, YEAR, TYPE, Site, TimSincFire, DBH:SNAG, SnagDens_23to38:fir_1km) %>%
  mutate(Nest = (TYPE == "NEST")*1)
dat.MOBL <- dat.MOBL %>% mutate(dist = nearNestDist(dat.MOBL)) %>%
  filter(!(TYPE == "RAND" & dist < 130)) %>%
  select(-dist)
  # filter out non-snag randoms to equalize with nests
dat.MOBL <- dat.MOBL %>% filter(!(TYPE == "RAND" & SNAG != 1))

#RSdat <- dat.MOBL %>% select(Nest, slope:fir_1km)
#pairs(RSdat) # Look for separation (none found)

#FCdat <- dat.MOBL %>% select(Nest, DBH:SnagDens_fir)
#pairs(FCdat) # Look for separation (none found)

rm(nearNestDist, i, vars.cont, keep, LiveNest, LiveRand, SnagNest, SnagRand)
save.image("Data_compiled.RData")
