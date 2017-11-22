setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")

#################################################
# Tabulate descriptive statistics of predictors #
#################################################

library(dplyr)
library(R.utils)
require(stringr)

load("Data_compiled.RData")

RS.BBWO <- loadObject("Model_RS_BBWO")
RS.HAWO <- loadObject("Model_RS_HAWO")
RS.WHWO <- loadObject("Model_RS_WHWO")
RS.NOFL <- loadObject("Model_RS_NOFL")
CB.BBWO <- loadObject("Model_CMB_BBWO")
CB.HAWO <- loadObject("Model_CMB_HAWO")
CB.WHWO <- loadObject("Model_CMB_WHWO")
CB.NOFL <- loadObject("Model_CMB_NOFL")

pars <- c(coef(RS.BBWO) %>% names,
          coef(RS.HAWO) %>% names,
          coef(RS.WHWO) %>% names,
          coef(RS.NOFL) %>% names,
          coef(CB.BBWO) %>% names,
          coef(CB.HAWO) %>% names,
          coef(CB.WHWO) %>% names,
          coef(CB.NOFL) %>% names) %>% unique

spp <- c("BBWO", "HAWO", "WHWO", "NOFL")
cols <- c(paste0("RS.", spp), paste0("CB.", spp))
out <- matrix("", nrow = length(pars), ncol = length(cols))
dimnames(out) <- list(pars, cols)

for(c in 1:length(cols)) {
  m <- eval(as.name(cols[c]))
  par <- coef(m) %>% names
  p <- summary(m)$coefficients[, "Pr(>|z|)"]
  signif <- character(length = length(p))
  signif[which(p <= 0.1)] <- "*"
  signif[which(p <= 0.05)] <- "**"
  out[par, cols[c]] <- summary(m)$coefficients[, c("Estimate")] %>%
    round(digits = 2) %>%
    paste0(" (",
           summary(m)$coefficients[, c("Std. Error")] %>%
             round(digits = 2),
           ")",
           signif)
}

write.csv(out, "Model_param_table.csv")
