require(dplyr)
require(stringr)
require(gtools)

######### Functions for Sierra woodpecker nest habitat modeling project ########################
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
      if(!((any(v == "BRKN_after") | any(v == "BRKN_before")) & any(v == "BRKN_after + TimSincFire + I(BRKN_after*TimSincFire)")))
      if(!((any(v == "BRKN_after") | any(v == "BRKN_before")) & any(v == "BRKN_before + TimSincFire + I(BRKN_before*TimSincFire)")))
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
      if(!(any(v == "sizlrg_loc") & any(v == "sizsm_loc"))) 
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

AIC_table <- function(dat.spp, vars, maxK) {
  mods <- model.construct(vars, maxK)
  out <- data.frame(model = c("Intercept-only", mods), stringsAsFactors = F)
  out$K <- 0
  out$AIC <- NA

  # Intercept-only model
  w <- rep(1, nrow(dat.spp))
  w[which(dat.spp$Nest == 0)] <- sum(dat.spp$Nest == 1) / sum(dat.spp$Nest == 0)
  mod <- glm(Nest ~ 1 , family = binomial, weights = w, data = dat.spp)
  out$K[1] <- mod$rank
  out$AIC[1] <- mod$aic
  
  # Other models
  for(m in 2:nrow(out)) {
    w <- rep(1, nrow(dat.spp))
    w[which(dat.spp$Nest == 0)] <- sum(dat.spp$Nest == 1) / sum(dat.spp$Nest == 0)
    mod <- glm(as.formula(paste0("Nest ~ ", out$model[m])) , family = binomial, weights = w, data = dat.spp)
    out$K[m] <- mod$rank
    out$AIC[m] <- mod$aic
  }
  out <- out %>% tbl_df %>% arrange(AIC)
  return(out)
}
################################################################################################
