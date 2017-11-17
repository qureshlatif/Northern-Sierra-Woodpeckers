##################################################################
# Contains functions used in model-fitting and plotting scripts. #
# Store this file in workspace designated in primary scripts     #
##################################################################

require(dplyr)
require(stringr)
require(gtools)
require(stats)

######### Functions for Sierra woodpecker nest habitat modeling project ########################
datList <- function(spp) {
  dat <- list()
  for(sp in 1:length(spp)) dat[[sp]] <- eval(as.name(paste0("dat.", spp[sp])))
  names(dat) <- spp
  return(dat)
}

getNbySpp <- function(spp, sites, dat = datList(spp)) { # dat = list of length(spp)
  n.val <- list()
  for(sp in 1:length(spp)) {
    n <- integer(length = length(sites)) # sizes of evaluation datasets
    for(s in 1:length(sites)) n[s] <- dat[[sp]] %>%
        filter(Site == sites[s]) %>% nrow
    n.val[[length(n.val)+1]] <- n
  }
  names(n.val) <- spp
  return(n.val)
}

model.construct <- function(vars, maxK) { # vars = vector of variable names
  mods <- c()
  for(j in 1:maxK) {
    index <- combinations(length(vars),j)
    for (a in 1:nrow(index)) {
      v <- vars[index[a,]]
      k <- v %>% paste0(collapse = "+") %>% str_split(fixed("+"))
      k <- k[[1]] %>% length
      if(k <= maxK)
      if(!(any(v == "DBH") & any(v == "DBH + I(DBH^2)")))
      if(!(any(v == "SP_PINE") & any(v == "SP_FIR")))
      if(!(any(v == "SP_FIR") & any(v == "SP_FIR + TimSincFire + I(SP_FIR*TimSincFire)")))
      if(!(any(v == "SP_PINE") & any(v == "SP_PINE + TimSincFire + I(SP_PINE*TimSincFire)")))
      if(!((any(v == "BRKN_after") | any(v == "BRKN_before")) & any(v == "BRKN")))
      if(!((any(v == "BRKN_after") | any(v == "BRKN_before") | any(v == "BRKN")) & any(v == "BRKN_after + TimSincFire + I(BRKN_after*TimSincFire)")))
      if(!((any(v == "BRKN_after") | any(v == "BRKN_before") | any(v == "BRKN")) & any(v == "BRKN_before + TimSincFire + I(BRKN_before*TimSincFire)")))
      if(!(any(v == "BRKN_after") & any(v == "BRKN_before")))
      if(!(any(v == "BRKN_after + TimSincFire + I(BRKN_after*TimSincFire)") & any(v == "BRKN_before + TimSincFire + I(BRKN_before*TimSincFire)")))
      if(!(any(v == "SnagDens_fir") & any(v == "SnagDens_23to38")))
      if(!(any(v == "SnagDens_fir") & any(v == "SnagDens_38to50")))
      if(!(any(v == "SnagDens_fir") & any(v == "SnagDens_pine")))
      if(!((any(v == "SnagDens_23to38") | any(v == "SnagDens_38to50") | any(v == "SnagDens_gt38")) & any(v == "SnagDens_23to50")))
      if(!((any(v == "SnagDens_23to38") | any(v == "SnagDens_38to50") | any(v == "SnagDens_23to50")) & any(v == "SnagDens_gt38")))
      if(!((any(v == "SnagDens_23to38") | any(v == "SnagDens_38to50") | any(v == "SnagDens_gt50") |
            any(v == "SnagDens_23to50") | any(v == "SnagDens_gt38") | any(v == "SnagDens_fir") |
            any(v == "SnagDens_pine")) &
           any(v == "SnagDens_all")))
      if(!(any(v == "canhi_loc") & (any(v == "canmd_loc") | any(v == "canlo_loc"))))
      if(!(any(v == "canhi_lnd") & (any(v == "canmd_lnd") | any(v == "canlo_lnd"))))
      if(!(any(v == "canmd_lnd") & (any(v == "can60_loc") | any(v == "can60_lnd"))))
      if(!(any(v == "canhi_loc") & (any(v == "can60_loc"))))
      if(!(any(v == "canhi_lnd") & (any(v == "can60_lnd"))))
      if(!(any(v == "sizlrg_loc") & any(v == "sizsm_loc"))) 
      if(!(any(v == "sizlrg_lnd") & any(v == "sizsm_lnd"))) 
      if(!(any(v == "cosasp") & (!any(v == "sinasp"))) &
        !(!any(v == "cosasp") & (any(v == "sinasp")))) {
          m <- paste(v,collapse="+")
            mods <- c(mods, m)
        }
      }
    }
  return(mods)
  }

SnsPlsSpcMax <- function(Obs, HSI) {        # Function for finding threshold that maximizes sensitivity + specificity
  HSI.nst <- HSI[which(Obs==1)]
  HSI.lnd <- HSI[which(Obs==0)]
  HSI.cnds <- seq(0.01,0.99,by=0.01)
  sens <- spcf <- Sum <- numeric(length=length(HSI.cnds))
  for (ii in 1:length(HSI.cnds)) {
    sens[ii] <- sum(HSI.nst>=HSI.cnds[ii])/length(HSI.nst)
    spcf[ii] <- sum(HSI.lnd<HSI.cnds[ii])/length(HSI.lnd)
    Sum[ii] <- sens[ii]+spcf[ii]
  } 
  SPS_th <- HSI.cnds[which(Sum==max(Sum))]
  SPS_sns <- sens[which(Sum==max(Sum))]
  SPS_spc <- spcf[which(Sum==max(Sum))]
  MSPS<-list(SPS_th,SPS_sns,SPS_spc)
  names(MSPS)<-c("HSI_thrshld","sensitivity_at_SPS","specificity_at_SPS")
  return(MSPS)
}

WLR_fit <- function(dat, formula, Obs = "Nest") { # Fit weighted logistic regression
  w <- rep(1, nrow(dat))
  w[which(dat[, Obs] == 0)] <- sum(dat[, Obs] == 1) / sum(dat[, Obs] == 0)
  dat$w <- w
  mod <- glm(formula, family = binomial, weights = w, data = dat)
  return(mod)
}

calcBins <- function(unit.size, n, bin.width) { # bin.width = number of units
  unitID <- rep(1:ceiling(n/unit.size), each = unit.size)
  unitID <- unitID[1:n]
  bin.st <- 1:(max(unitID) - bin.width + 1)
  bin.end <- bin.width:max(unitID)
  return(list(unitID = unitID, st = bin.st, end = bin.end))
}

AIC_table <- function(dat.spp, vars, maxK) {
  mods <- model.construct(vars, maxK)
  out <- data.frame(model = c("Intercept-only", mods), stringsAsFactors = F)
  out$K <- 0
  out$AIC <- NA

  # Intercept-only model
  mod <- WLR_fit(dat.spp, formula = Nest ~ 1)
  out$K[1] <- mod$rank
  out$AIC[1] <- mod$aic
  
  # Other models
  for(m in 2:nrow(out)) {
    mod <- WLR_fit(dat.spp, formula = as.formula(paste0("Nest ~ ", out$model[m])))
    out$K[m] <- mod$rank
    out$AIC[m] <- mod$aic
  }
  out <- out %>% tbl_df %>% arrange(AIC)
  return(out)
}


AllFit <- function(spp, sites, vars , dat = datList(spp)) { # dat = list of datasets with length = length(spp)
  modList <- list()
  for(sp in 1:length(spp)) {
    dat.spp <- dat[[sp]]
    if(is.list(vars)) {
      maxK <- min(floor(sum(dat.spp$Nest)/10), length(vars[[sp]]))
      mods <- model.construct(vars[[sp]], maxK)
    } else {
      maxK <- min(floor(sum(dat.spp$Nest)/10), length(vars))
      mods <- model.construct(vars, maxK)
    }
    
    out <- data.frame(model = mods, stringsAsFactors = F)
    out$RPI_CH <-  out$RPI_ML <- out$RPI_CB <- out$AUC_CH <-  out$AUC_ML <- out$AUC_CB <- out$K <- 0
    
    # Leave out each wildfire (site) for cross-validation #
    for(s in 1:length(sites)) {
      dat.cal <- dat.spp %>% filter(Site != sites[s]) # Calibration dataset
      dat.val <- dat.spp %>% filter(Site == sites[s]) # Validation dataset (held out fire)
      
      # Parameters for equal-area moving window bins to calculate RPI
      bins <- calcBins(unit.size[[sp]][s], n.val[[sp]][s], bin.width[[sp]][s])

      for(m in 1:nrow(out)) {
        mod <- WLR_fit(dat = dat.cal, formula = as.formula(paste0("Nest ~ ", out$model[m])))
        out$K[m] <- mod$rank
        obs <- dat.val$Nest # observed validation
        p <- as.vector(plogis(predict(mod, dat.val))) # predicted validation
        
        # AUC #
        out[m, s + 2] <- auc(cbind(rep(1, nrow(dat.val)), dat.val$Nest, p))$AUC
        
        # Calculate RPI #
        p.norm <- (p - min(p))/(max(p) - min(p)) # normalized predictions
        obs <- obs[order(p.norm)]
        p.norm <- sort(p.norm)
        obs.bin <- pred.bin <- numeric(length = length(bins$st))
        for(b in 1:length(bins$st)) {
          st <- min(which(bins$unitID == bins$st[b]))
          end <- max(which(bins$unitID == bins$end[b]))
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
    
    out <- out[which(!is.na(out$RPI)),] # Drop models with RPI < critical value at any one wildfire
    
    # Selection criteria calculated for remaining models fitted to all data #
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
    coefs <- unique(str_trim(unlist(vars, fixed("+"))))
    for(i in 1:nrow(out)) {
      v.ind <- str_detect(fixed(out$model[i]), fixed(coefs))
      names(v.ind) <- coefs
      if(any(names(v.ind) == "I(DBH^2)")) if(v.ind["I(DBH^2)"] == T) v.ind["DBH"] <- F
      if(any(names(v.ind) == "I(SP_PINE*TimSincFire)")) if(v.ind["I(SP_PINE*TimSincFire)"] == T)
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
      select(mod_coefs_sign, Sens_check, MaxGain_check, All_coefs_sign, AIC_top, K, RPI, RPI_sd,
             AUC, AIC, deviance, AUC_CB:RPI_CH, maxSSS:SSS_CH) %>%
      rename(model = mod_coefs_sign)
    modList[[length(modList) + 1]] <- out
  }
  names(modList) <- spp
  return(modList)
}

################################################################################################
