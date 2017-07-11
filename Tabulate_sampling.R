setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")
load("Data_compiled.RData")

library(R.utils)
library(gtools)
library(dplyr)

dat <- dat %>%
  mutate(SPECIES = replace(SPECIES, which(is.na(dat$SPECIES)), "RAND")) %>%
  group_by(SPECIES) %>%
  summarise(CB = sum(Site == "CB"), ML = sum(Site == "ML"), CH = sum(Site == "CH"),
            Y1 = sum(TimSincFire == 1), Y2 = sum(TimSincFire == 2), Y3 = sum(TimSincFire == 3),
            Y4 = sum(TimSincFire == 4), Y5 = sum(TimSincFire == 5))

write.csv(n.summary, "Sampling_allSpp_X_site,year.csv")

sum.BBWO <- dat.BBWO %>%
  mutate(TYPE = paste0("BBWO.",TYPE)) %>%
  group_by(TYPE) %>%
  summarise(CB = sum(Site == "CB"), ML = sum(Site == "ML"), CH = sum(Site == "CH"),
            Y1 = sum(TimSincFire == 1), Y2 = sum(TimSincFire == 2), Y3 = sum(TimSincFire == 3),
            Y4 = sum(TimSincFire == 4), Y5 = sum(TimSincFire == 5))

sum.HAWO <- dat.HAWO %>%
  mutate(TYPE = paste0("HAWO.",TYPE)) %>%
  group_by(TYPE) %>%
  summarise(CB = sum(Site == "CB"), ML = sum(Site == "ML"), CH = sum(Site == "CH"),
            Y1 = sum(TimSincFire == 1), Y2 = sum(TimSincFire == 2), Y3 = sum(TimSincFire == 3),
            Y4 = sum(TimSincFire == 4), Y5 = sum(TimSincFire == 5))

sum.WHWO <- dat.WHWO %>%
  mutate(TYPE = paste0("WHWO.",TYPE)) %>%
  group_by(TYPE) %>%
  summarise(CB = sum(Site == "CB"), ML = sum(Site == "ML"), CH = sum(Site == "CH"),
            Y1 = sum(TimSincFire == 1), Y2 = sum(TimSincFire == 2), Y3 = sum(TimSincFire == 3),
            Y4 = sum(TimSincFire == 4), Y5 = sum(TimSincFire == 5))

sum.NOFL <- dat.NOFL %>%
  mutate(TYPE = paste0("NOFL.",TYPE)) %>%
  group_by(TYPE) %>%
  summarise(CB = sum(Site == "CB"), ML = sum(Site == "ML"), CH = sum(Site == "CH"),
            Y1 = sum(TimSincFire == 1), Y2 = sum(TimSincFire == 2), Y3 = sum(TimSincFire == 3),
            Y4 = sum(TimSincFire == 4), Y5 = sum(TimSincFire == 5))

sum.MOBL <- dat.MOBL %>%
  mutate(TYPE = paste0("MOBL.",TYPE)) %>%
  group_by(TYPE) %>%
  summarise(CB = sum(Site == "CB"), ML = sum(Site == "ML"), CH = sum(Site == "CH"),
            Y1 = sum(TimSincFire == 1), Y2 = sum(TimSincFire == 2), Y3 = sum(TimSincFire == 3),
            Y4 = sum(TimSincFire == 4), Y5 = sum(TimSincFire == 5))

sum.RBSA <- dat.RBSA %>%
  mutate(TYPE = paste0("RBSA.",TYPE)) %>%
  group_by(TYPE) %>%
  summarise(CB = sum(Site == "CB"), ML = sum(Site == "ML"), CH = sum(Site == "CH"),
            Y1 = sum(TimSincFire == 1), Y2 = sum(TimSincFire == 2), Y3 = sum(TimSincFire == 3),
            Y4 = sum(TimSincFire == 4), Y5 = sum(TimSincFire == 5))

sum.all <- rbind(sum.BBWO, sum.HAWO, sum.WHWO, sum.NOFL, sum.MOBL, sum.RBSA)

write.csv(sum.all, "Sampling_spp_X_site,year.csv")

Y.summary <- dat %>% group_by(Site) %>% filter(SPECIES == "RAND") %>%
  summarise(Y1 = sum(TimSincFire == 1), Y2 = sum(TimSincFire == 2), Y3 = sum(TimSincFire == 3),
            Y4 = sum(TimSincFire == 4), Y5 = sum(TimSincFire == 5))

write.csv(Y.summary, "Sampling_site_X_year.csv")
