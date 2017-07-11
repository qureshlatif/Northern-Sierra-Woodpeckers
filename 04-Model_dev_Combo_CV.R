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
source("scripts/Functions.R")

# Species with sufficient sample sizes for cross-fire validation #
spp <- c("BBWO", "HAWO", "WHWO", "NOFL")

## Set bin parameters for calculating RPIs ##
# Review sample sizes by site for each species #
n.val <- list()
for(sp in 1:length(spp)) {
  n <- integer(length = length(sites)) # sizes of evaluation datasets
  for(s in 1:length(sites)) n[s] <- eval(as.name(paste0("dat.",spp[sp]))) %>%
      filter(Site == sites[s]) %>% nrow
  n.val[[length(n.val)+1]] <- n
}
names(n.val) <- spp

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

# Field-collected variable combos #
vars.rs <- list(BBWO = c("ccmort_loc", "blk_lndcc", "canhi_loc"),
                HAWO = c("ccmort_loc", "canhi_lnd", "sizlrg_loc"),
                WHWO = c("ccmort_loc", "blk_lndcc"),
                NOFL = c("slope", "ccmort_loc", "sizlrg_loc", "sizlrg_lnd"))
vars.fc <- list(BBWO = c("DBH", "DBH + I(DBH^2)", "SnagDens_all"),
                HAWO = c("DBH", "DBH + I(DBH^2)", "SP_PINE", "BRKN",
                         "SnagDens_all", "SP_PINE + TimSincFire + I(SP_PINE*TimSincFire)"),
                WHWO = c("DBH", "DBH + I(DBH^2)", "SP_FIR", "BRKN", "SnagDens_gt50", "SnagDens_all"),
                NOFL = c("DBH", "DBH + I(DBH^2)", "SP_FIR", "BRKN", "SnagDens_all"))

for(sp in 1:length(spp)) {
  dat.spp <- eval(as.name(paste0("dat.", spp[sp])))
  vars <- c(vars.rs[[sp]], vars.fc[[sp]])
  maxK <- min(floor(sum(dat.spp$Nest)/10), length(vars))
  mods <- model.construct(vars, maxK)
  
  out <- data.frame(model = mods, stringsAsFactors = F)
  out$RPI_CH <-  out$RPI_ML <- out$RPI_CB <- out$AUC_CH <-  out$AUC_ML <- out$AUC_CB <- out$K <- 0
  for(s in 1:length(sites)) {
    dat.cal <- dat.spp %>% filter(Site != sites[s])
    dat.val <- dat.spp %>% filter(Site == sites[s])
    w <- rep(1, nrow(dat.cal))
    w[which(dat.cal$Nest == 0)] <- sum(dat.cal$Nest == 1) / sum(dat.cal$Nest == 0)
    
    # Parameters for equal-area moving window bins to calculate RPI
    unitID <- rep(1:ceiling(n.val[[sp]][s]/unit.size[[sp]][s]), each = unit.size[[sp]][s])[1:n.val[[sp]][s]]
    bin.st <- 1:(max(unitID) - bin.width[[sp]][s] + 1)
    bin.end <- bin.width[[sp]][s]:max(unitID)
    
    for(m in 1:nrow(out)) {
      mod <- glm(as.formula(paste0("Nest ~ ", out$model[m])) , family = binomial, weights = w, data = dat.cal)
      out$K[m] <- mod$rank
      obs <- dat.val$Nest
      p <- as.vector(plogis(predict(mod, dat.val)))
      out[m, s + 2] <- auc(cbind(rep(1, nrow(dat.val)), dat.val$Nest, p))$AUC
      p.norm <- (p - min(p))/(max(p) - min(p)) # normalized predictions
      obs <- obs[order(p.norm)]
      p.norm <- sort(p.norm)
      obs.bin <- pred.bin <- numeric(length = length(bin.st))
      for(b in 1:length(bin.st)) {
        st <- min(which(unitID == bin.st[b]))
        end <- max(which(unitID == bin.end[b]))
        obs.bin[b] <- sum(obs[st:end])/length(st:end)
        pred.bin[b] <- mean(p.norm[st:end])
      }
      rpi <- cor(obs.bin, pred.bin, method = "spearman")
      ifelse(rpi >= sprmn.crit, out[m, s + 5] <- rpi, out[m, s + 5] <- NA)
    }
  }
  out$AUC <- apply(as.matrix(out[,3:5]), 1, mean)
  out$RPI <- apply(as.matrix(out[,6:8]), 1, mean)
  out$RPI_sd <- apply(as.matrix(out[,6:8]), 1, sd)
  
  out <- out[which(!is.na(out$RPI)),]
  out$SSS_CH <- out$SSS_ML <- out$SSS_CB <- out$Sens_CH <- out$Sens_ML <- out$Sens_CB <-
    out$HSIthrsh <- out$maxSSS <- numeric(length = nrow(out))
  out$mod_coefs_sign <- ""
  out$AIC <- NA
  out$deviance <- NA
  
  for(m in 1:nrow(out)) {
    w <- rep(1, nrow(dat.spp))
    w[which(dat.spp$Nest == 0)] <- sum(dat.spp$Nest == 1) / sum(dat.spp$Nest == 0)
    mod <- glm(as.formula(paste0("Nest ~ ", out$model[m])) , family = binomial, weights = w, data = dat.spp)
    out$AIC[m] <- mod$aic
    out$deviance[m] <- mod$deviance
    coef.names <- str_split(out$model[m], fixed("+"))[[1]]
    coef.est <- round(mod$coefficients[-1], digits = 2)
    coef.p <- summary(mod)$coefficients[-1,"Pr(>|z|)"]
    coef.sign <- character(length = length(coef.p))
    coef.sign[which(coef.p < 0.05)] <- "*"
    coef.sign[which(coef.p >= 0.05 & coef.p < 0.1)] <- "."
    out$mod_coefs_sign[m] <- paste(paste0(coef.names, " (", coef.est, coef.sign, ")"), collapse = " + ")
    p <- as.vector(plogis(predict(mod, dat.spp)))
    obs <- dat.spp$Nest
    maxSSS <- SnsPlsSpcMax(obs, p)
    if(length(maxSSS$HSI_thrshld) == 1){
      out$maxSSS[m] <- maxSSS$sensitivity_at_SPS + maxSSS$specificity_at_SPS
      out$HSIthrsh[m] <- maxSSS$HSI_thrshld
    }
    else {
      ind <- which(maxSSS$HSI_thrshld == min(maxSSS$HSI_thrshld)) # Favor sensitivity by picking lower HSI when more than one maxSSS threshold
      out$maxSSS[m] <- maxSSS$sensitivity_at_SPS[ind] + maxSSS$specificity_at_SPS[ind]
      out$HSIthrsh[m] <- maxSSS$HSI_thrshld[ind]
    }
    for(s in 1:length(sites)) {
      dat.val <- dat.spp %>% filter(Site == sites[s])
      p <- as.vector(plogis(predict(mod, dat.val)))
      obs <- dat.val$Nest
      out[m, paste0("Sens_", sites[s])] <- sum(p[which(obs == 1)] >= out$HSIthrsh[m]) / sum(obs == 1)
      out[m, paste0("SSS_", sites[s])] <- sum(p[which(obs == 1)] >= out$HSIthrsh[m]) / sum(obs == 1) +
        sum(p[which(obs == 0)] < out$HSIthrsh[m]) / sum(obs == 0)
    }
  }
  out <- out[order(out$RPI, decreasing = T),] # Sort by RPI

  check <- numeric() # To indicate whether all coefficients have p < 0.05
  coefs <- unique(str_trim(c(unlist(vars.rs), unlist(str_split(unlist(vars.fc), fixed("+"))))))
  for(i in 1:nrow(out)) {
    v.ind <- str_detect(fixed(out$model[i]), fixed(coefs))
    names(v.ind) <- coefs
    if(v.ind["I(DBH^2)"] == T) v.ind["DBH"] <- F
    if(v.ind["I(SP_PINE*TimSincFire)"] == T)
      v.ind["SP_PINE"] <- v.ind["TimSincFire"] <- F

    # Trim check of statistical significance to required parameters
        #(e.g., highest order terms in the presence of quadratics or interactions)
    chk <- out$mod_coefs_sign[i] %>%
      str_split(fixed("+"), simplify = T) %>% str_trim %>%
      str_sub(start = -2, end = -1) == "*)"
    names(chk) <- out$model[i] %>%
      str_split(fixed("+"), simplify = T) %>% str_trim
    v.ind <- v.ind[which(v.ind & names(v.ind) %in% names(chk))]
    chk <- chk[which(names(chk) %in% names(v.ind))]
    
    check <- c(check, (sum(chk) == length(chk))*1)
  }

  out <- out %>% tbl_df %>%
    mutate(Sens_check = (Sens_CB >= 0.5 & Sens_ML >= 0.5 & Sens_CH >= 0.5)) %>%
    mutate(MaxGain_check = (SSS_CB > 1 & SSS_ML > 1 & SSS_CH > 1)) %>%
    mutate(All_coefs_sign = check) %>%
    mutate(AIC_top = AIC < (min(AIC) + 6)) %>%
    select(model, mod_coefs_sign, AIC_top, Sens_check, All_coefs_sign, MaxGain_check, K, RPI, RPI_sd,
           AUC, AIC, deviance, AUC_CB:RPI_CH, maxSSS:SSS_CH)
  assign(paste0("CMBmods.", spp[sp]), out)
}

#Cleanup
rm(out, b, bin.end, bin.st, check, coef.est, coef.names, coef.p, coef.sign, end, i, ind,
   m, maxK, maxSSS, mod, mods, n, obs, obs.bin, p, p.norm, pred.bin, rpi, s, sp,
   st, unitID, w, coefs, chk, v.ind, vars, vars.fc, vars.rs)

## Manually review remaining models, and select models that:
  ## BBWO ##
write.csv(CMBmods.BBWO, "Models_for_review_BBWO.csv", row.names = F)
# Make sure selected model's coefficients are consistent in direction when any one fire is left out.
dat.fold <- dat.BBWO[which(dat.BBWO$Site!=sites[3]),]
w <- rep(1, nrow(dat.fold))
w[which(dat.fold$Nest == 0)] <- sum(dat.fold$Nest == 1) / sum(dat.fold$Nest == 0)
mod <- glm(Nest ~ ccmort_loc+blk_lndcc+DBH + I(DBH^2) + SnagDens_pine , family = binomial, weights = w,
           data = dat.fold)
summary(mod)
rm(dat.fold, w)

  # Save selected model #
w <- rep(1, nrow(dat.BBWO))
w[which(dat.BBWO$Nest == 0)] <- sum(dat.BBWO$Nest == 1) / sum(dat.BBWO$Nest == 0)
mod <- glm(Nest ~ ccmort_loc+blk_lndcc+DBH + I(DBH^2) + SnagDens_pine , family = binomial, weights = w, data = dat.BBWO)
summary(mod)
saveObject(mod, "Model_CMB_BBWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_CMB_BBWO.csv")


## HAWO ##
write.csv(CMBmods.HAWO, "Models_for_review_HAWO.csv", row.names = F)
  # Make sure selected model's coefficients are consistent in direction when any one fire is left out.
dat.fold <- dat.HAWO[which(dat.HAWO$Site!=sites[3]),]
w <- rep(1, nrow(dat.fold))
w[which(dat.fold$Nest == 0)] <- sum(dat.fold$Nest == 1) / sum(dat.fold$Nest == 0)
mod <- glm(Nest ~ ccmort_loc+DBH+BRKN_before+SnagDens_23to38 , family = binomial, weights = w,
           data = dat.fold)
summary(mod)
rm(dat.fold, w)

  # Save selected model #
w <- rep(1, nrow(dat.HAWO))
w[which(dat.HAWO$Nest == 0)] <- sum(dat.HAWO$Nest == 1) / sum(dat.HAWO$Nest == 0)
mod <- glm(Nest ~ ccmort_loc+DBH+BRKN_before+SnagDens_23to38 , family = binomial, weights = w, data = dat.HAWO)
summary(mod)
saveObject(mod, "Model_CMB_HAWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_CMB_HAWO.csv")


## WHWO ##
write.csv(CMBmods.WHWO, "Models_for_review_WHWO.csv", row.names = F)

  # Make sure selected model's coefficients are consistent in direction when any one fire is left out.
dat.fold <- dat.WHWO[which(dat.WHWO$Site!=sites[3]),]
w <- rep(1, nrow(dat.fold))
w[which(dat.fold$Nest == 0)] <- sum(dat.fold$Nest == 1) / sum(dat.fold$Nest == 0)
mod <- glm(Nest ~ ccmort_loc + blk_lndcc + DBH + SP_PINE + BRKN_before , family = binomial, weights = w,
           data = dat.fold)
summary(mod)
rm(dat.fold, w)

  # Save selected model #
w <- rep(1, nrow(dat.WHWO))
w[which(dat.WHWO$Nest == 0)] <- sum(dat.WHWO$Nest == 1) / sum(dat.WHWO$Nest == 0)
mod <- glm(Nest ~ ccmort_loc + blk_lndcc + DBH + SP_PINE + BRKN_before , family = binomial, weights = w, data = dat.WHWO)
summary(mod)
saveObject(mod, "Model_CMB_WHWO") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_CMB_WHWO.csv")


## NOFL ##
write.csv(CMBmods.NOFL, "Models_for_review_NOFL.csv", row.names = F)

# Tabulate top condition by time since fire
#dat.NOFL %>% filter(Nest == 1) %>% select(BRKN_before, BRKN_after, TimSincFire) %>%
#  mutate(TimSincFire = TimSincFire*scale.factors["TimSincFire", "SD"] +
#           scale.factors["TimSincFire", "mean"]) %>% group_by(TimSincFire) %>%
#  summarise(BRKN_before = sum(BRKN_before) / n(), BRKN_after = sum(BRKN_after) / n(), n = n())

  # Make sure selected model's coefficients are consistent in direction when any one fire is left out.
dat.fold <- dat.NOFL[which(dat.NOFL$Site!=sites[3]),]
w <- rep(1, nrow(dat.fold))
w[which(dat.fold$Nest == 0)] <- sum(dat.fold$Nest == 1) / sum(dat.fold$Nest == 0)
mod <- glm(Nest ~ slope + DBH + I(DBH^2) + SNAG_DEC + BRKN_after , family = binomial, weights = w,
           data = dat.fold)
summary(mod)
rm(dat.fold, w)

  # Save selected model #
w <- rep(1, nrow(dat.NOFL))
w[which(dat.NOFL$Nest == 0)] <- sum(dat.NOFL$Nest == 1) / sum(dat.NOFL$Nest == 0)
mod <- glm(Nest ~ slope + DBH + I(DBH^2) + SNAG_DEC + BRKN_after , family = binomial, weights = w, data = dat.NOFL)
summary(mod)
saveObject(mod, "Model_CMB_NOFL") # Save selected model #
write.csv(summary(mod)$coefficients, "Model_CMB_NOFL.csv")
