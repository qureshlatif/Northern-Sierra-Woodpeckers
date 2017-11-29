#####################################################################
# Purpose: Find maxSSS thresholds for classifying suitable habitat. #
#####################################################################

require(R.utils)
require(dplyr)

setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")
load("Data_compiled.RData") # Workspace containing data
source("Northern-Sierra-Woodpeckers/Functions.R") ## Script file with functions for model-fitting and plotting ##

rows <- c("RS.BBWO", "RS.HAWO", "RS.WHWO", "RS.NOFL",
          "CB.BBWO", "CB.HAWO", "CB.WHWO", "CB.NOFL")
cols <- c("Threshold", "sensitivity", "specificity")
out <- matrix(NA, nrow = length(rows), ncol = length(cols))
dimnames(out) <- list(rows, cols)

## Remotely sensed models ##
# Models #
mod.BBWO <- loadObject("Model_RS_BBWO")
mod.HAWO <- loadObject("Model_RS_HAWO")
mod.WHWO <- loadObject("Model_RS_WHWO")
mod.NOFL <- loadObject("Model_RS_NOFL")

# Apply models #
dat.BBWO <- dat.BBWO %>%
  mutate(HSI = predict(mod.BBWO, ., type = "response"))
dat.HAWO <- dat.HAWO %>%
  mutate(HSI = predict(mod.HAWO, ., type = "response"))
dat.WHWO <- dat.WHWO %>%
  mutate(HSI = predict(mod.WHWO, ., type = "response"))
dat.NOFL <- dat.NOFL %>%
  mutate(HSI = predict(mod.NOFL, ., type = "response"))

# maxSSS thresholds #
calc <- SnsPlsSpcMax(dat.BBWO$Nest, dat.BBWO$HSI)
out["RS.BBWO", ] <- unlist(calc)
calc <- SnsPlsSpcMax(dat.HAWO$Nest, dat.HAWO$HSI)
out["RS.HAWO", ] <- unlist(calc)
calc <- SnsPlsSpcMax(dat.WHWO$Nest, dat.WHWO$HSI)
out["RS.WHWO", ] <- unlist(calc)
calc <- SnsPlsSpcMax(dat.NOFL$Nest, dat.NOFL$HSI)
out["RS.NOFL", ] <- unlist(calc)

## Remotely sensed + field collected models ##
# Models #
mod.BBWO <- loadObject("Model_CMB_BBWO")
mod.HAWO <- loadObject("Model_CMB_HAWO")
mod.WHWO <- loadObject("Model_CMB_WHWO")
mod.NOFL <- loadObject("Model_CMB_NOFL")

# Apply models #
dat.BBWO <- dat.BBWO %>%
  mutate(HSI = predict(mod.BBWO, ., type = "response"))
dat.HAWO <- dat.HAWO %>%
  mutate(HSI = predict(mod.HAWO, ., type = "response"))
dat.WHWO <- dat.WHWO %>%
  mutate(HSI = predict(mod.WHWO, ., type = "response"))
dat.NOFL <- dat.NOFL %>%
  mutate(HSI = predict(mod.NOFL, ., type = "response"))

# maxSSS thresholds #
calc <- SnsPlsSpcMax(dat.BBWO$Nest, dat.BBWO$HSI)
out["CB.BBWO", ] <- unlist(calc)
calc <- SnsPlsSpcMax(dat.HAWO$Nest, dat.HAWO$HSI)
out["CB.HAWO", ] <- unlist(calc)
calc <- SnsPlsSpcMax(dat.WHWO$Nest, dat.WHWO$HSI)
out["CB.WHWO", ] <- unlist(calc)
calc <- SnsPlsSpcMax(dat.NOFL$Nest, dat.NOFL$HSI)
out["CB.NOFL", ] <- unlist(calc)

write.csv(out, "MaxSSS_thresholds.csv")
