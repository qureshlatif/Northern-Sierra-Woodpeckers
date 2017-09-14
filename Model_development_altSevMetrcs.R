setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")
load("Data_compiled.RData")

library(R.utils)
library(gtools)
library(dplyr)
library(PresenceAbsence)
library(foreign)
sites <- unique(dat.BBWO$Site)

dat.add <- read.dbf("E:/GISData/PtBlue_Sierra/NR_points_altSevMets.dbf", as.is = T) %>% tbl_df %>%
  select(SAMPLE_ID, ccmort_1km:cmgt64_500) %>%
  mutate(cmgt64_500 = 100*cmgt64_500)
dat.BBWO <- dat.BBWO %>% left_join(dat.add, by = "SAMPLE_ID")
rm(dat.add)


w <- rep(1, nrow(dat.BBWO))
w[which(dat.BBWO$Nest == 0)] <- sum(dat.BBWO$Nest == 1) / sum(dat.BBWO$Nest == 0)

## Correlations ##
cor(dat.BBWO[,c("ccmort_loc", "blk_lndcc", "ccmort_1km", "ccmort_500", "cmgt64_500")])

## mean ccmort 1km ##
mod <- glm(Nest ~ ccmort_loc + ccmort_1km + canhi_loc, family = binomial, weights = w, data = dat.BBWO)
summary(mod)

## mean ccmort 500 m ##
mod <- glm(Nest ~ ccmort_loc + ccmort_500 + canhi_loc, family = binomial, weights = w, data = dat.BBWO)
summary(mod)

## mean cmgt64 500 m ##
mod <- glm(Nest ~ ccmort_loc + cmgt64_500 + canhi_loc, family = binomial, weights = w, data = dat.BBWO)
summary(mod)

## quadratic cmort local ##
dat.BBWO <- dat.BBWO %>% mutate(ccmort_loc.z = (ccmort_loc - mean(ccmort_loc))/sd(ccmort_loc))
mod <- glm(Nest ~ ccmort_loc.z + I(ccmort_loc.z^2) + canhi_loc, family = binomial, weights = w, data = dat.BBWO)
summary(mod)

## quadratic cmort 1km ##
dat.BBWO <- dat.BBWO %>% mutate(ccmort_lnd.z = (ccmort_1km - mean(ccmort_1km))/sd(ccmort_1km))
mod <- glm(Nest ~ ccmort_loc.z + ccmort_lnd.z + I(ccmort_lnd.z^2) + canhi_loc, family = binomial, weights = w, data = dat.BBWO)
summary(mod)

## cmgt64 1km by site ##
dat.site <- dat.BBWO %>% filter(Site != "CB")
w <- rep(1, nrow(dat.site))
w[which(dat.site$Nest == 0)] <- sum(dat.site$Nest == 1) / sum(dat.site$Nest == 0)
mod <- glm(Nest ~ ccmort_loc + ccmort_1km + canhi_loc, family = binomial, weights = w, data = dat.site)
summary(mod)

dat.site <- dat.BBWO %>% filter(Site != "ML")
w <- rep(1, nrow(dat.site))
w[which(dat.site$Nest == 0)] <- sum(dat.site$Nest == 1) / sum(dat.site$Nest == 0)
mod <- glm(Nest ~ ccmort_loc + ccmort_1km + canhi_loc, family = binomial, weights = w, data = dat.site)
summary(mod)

dat.site <- dat.BBWO %>% filter(Site != "CH")
w <- rep(1, nrow(dat.site))
w[which(dat.site$Nest == 0)] <- sum(dat.site$Nest == 1) / sum(dat.site$Nest == 0)
mod <- glm(Nest ~ ccmort_loc + ccmort_1km + canhi_loc, family = binomial, weights = w, data = dat.site)
summary(mod)
