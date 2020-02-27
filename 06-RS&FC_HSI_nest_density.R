############################################################################################
# Purpose: Relate observed nest densities for 4 species (BBWO, HAWO, WHWO, NOFL) with HSIs #
# from remote-sensed & field-collected models.                                             #
############################################################################################
##***Assumes non-nest points represent proportion area

library(foreign)
library(R.utils)
#devtools::install_github("qureshlatif/WoodpeckerHSI") # Run this upon first use.
  # Note: You might get some disconcerting warnings upon first install.
    #To avoid them, restart R after installing package.
library(WoodpeckerHSI)
#library(dplyr)
#library(ggplot2)
#library(cowplot)

#____________________________________ Inputs _______________________________________#
setwd("C:/Users/Quresh.Latif/files/projects/prior/PtBlue_Sierra/")
load("Data_compiled.RData") # Workspace containing data
NR_points <- "C:/Users/Quresh.Latif/files/GIS/prior/PtBlue_Sierra/NR_points.dbf" # Needs to contain nest points with HSIs (randoms not used from here)

spp <- c("BBWO", "HAWO", "WHWO", "NOFL")
area <- 9 * 4 #Area surveyed in 100 ha units X study duration.

transects <- unique(substr(dat$SAMPLE_ID, 1, 4)) # Needed for bootstrapping at transect level
R <- 5000 # Desired number of bootstrapped samples for calculating Density uncertainty 
thresholds <- list(BBWO = c(0.2, 0.55), HAWO = c(0.3, 0.55),
                   WHWO = c(0.3, 0.63), NOFL = c(0.68))
#___________________________________________________________________________________#

dat.grid <- datz %>% filter(TYPE == "RAND")

# Models #
mod.BBWO <- loadObject("Model_CMB_BBWO")
mod.HAWO <- loadObject("Model_CMB_HAWO")
mod.WHWO <- loadObject("Model_CMB_WHWO")
mod.NOFL <- loadObject("Model_CMB_NOFL")

# Apply models #
dat.grid <- dat.grid %>%
  mutate(HSI_BBWO = predict(mod.BBWO, ., type = "response")) %>%
  mutate(HSI_HAWO = predict(mod.HAWO, ., type = "response")) %>%
  mutate(HSI_WHWO = predict(mod.WHWO, ., type = "response")) %>%
  mutate(HSI_NOFL = predict(mod.NOFL, ., type = "response")) %>%
  mutate(Transect = str_sub(SAMPLE_ID, 1, 4)) %>%
  select(SAMPLE_ID, Transect, HSI_BBWO:HSI_NOFL)
nest.BBWO <- dat.BBWO %>%
  filter(Nest == 1) %>%
  mutate(HSI = predict(mod.BBWO, ., type = "response")) %>%
  mutate(Transect = str_sub(SAMPLE_ID, 1, 4)) %>%
  select(SAMPLE_ID, Transect, HSI)
nest.HAWO <- dat.HAWO %>%
  filter(Nest == 1) %>%
  mutate(HSI = predict(mod.HAWO, ., type = "response")) %>%
  mutate(Transect = str_sub(SAMPLE_ID, 1, 4)) %>%
  select(SAMPLE_ID, Transect, HSI)
nest.WHWO <- dat.WHWO %>%
  filter(Nest == 1) %>%
  mutate(HSI = predict(mod.WHWO, ., type = "response")) %>%
  mutate(Transect = str_sub(SAMPLE_ID, 1, 4)) %>%
  select(SAMPLE_ID, Transect, HSI)
nest.NOFL <- dat.NOFL %>%
  filter(Nest == 1) %>%
  mutate(HSI = predict(mod.NOFL, ., type = "response")) %>%
  mutate(Transect = str_sub(SAMPLE_ID, 1, 4)) %>%
  select(SAMPLE_ID, Transect, HSI)

#_____________ Get functions (delete after replacing with WoodpeckerHSI package) _____________#
#source("Northern-Sierra-Woodpeckers/Functions.R") ## Script file with functions for model-fitting and plotting ##
#source("Northern-Sierra-Woodpeckers/Functions_plotting.R") ## Additional functions used in this script only ##
#_________________________________________#

# Bins for densities #
bins <- calcBins(50, nrow(dat.grid), 3)

# Tabulate values for plotting observed densities for moving window HSI bins #
for(s in spp) {
  nests <- eval(as.name(paste0("nest.", s)))$HSI
  grid <- dat.grid %>%
    select(matches(paste0("HSI_", s)))
  names(grid)[1] <- "HSI"
  grid <- grid$HSI
  assign(paste0("tab.", s),
         calcBinDensities(nests, grid, bins, area))
}
rm(s, nests, grid)

# for(s in spp) {
#   thrs <- thresholds[[s]]
#   nests <- eval(as.name(paste0("nest.", s)))$HSI
#   datg <- dat.grid %>% mutate(HSI = eval(as.name(paste0("HSI_", s))))
#   grid <- datg$HSI
#   dat.class <- calcClassDensities(nests, grid, thrs, area)
# 
#   # Add bootstrapped CIs #
#   dat.class <- dat.class %>% HSIClassDensityBS(dat.sample = eval(as.name(paste0("nest.", s))),
#                                          datg, transects, thrs, area, R, UnitID = "Transect")
#   assign(paste0("dat.class.", s), dat.class)
# }
# rm(thrs, dat.class, nests, grid, datg)

# Cache and retrieve #
#write.csv(dat.class.BBWO, "Plot_cache_CMB_BBWO_class.csv", row.names = F)
#write.csv(dat.class.HAWO, "Plot_cache_CMB_HAWO_class.csv", row.names = F)
#write.csv(dat.class.WHWO, "Plot_cache_CMB_WHWO_class.csv", row.names = F)
#write.csv(dat.class.NOFL, "Plot_cache_CMB_NOFL_class.csv", row.names = F)
dat.class.BBWO <- read.csv("Plot_cache_CMB_BBWO_class.csv", header = T)
dat.class.HAWO <- read.csv("Plot_cache_CMB_HAWO_class.csv", header = T)
dat.class.WHWO <- read.csv("Plot_cache_CMB_WHWO_class.csv", header = T)
dat.class.NOFL <- read.csv("Plot_cache_CMB_NOFL_class.csv", header = T)

#__________________________________ PLOTTING (manuscript version) ________________________________________#
##### BBWO #####
#___Inputs___#
s <- "BBWO"
thresholds <- c(0.2, 0.55) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 20

# Parameters for additional plot labels
classLabSize <- 10
labxy <- rbind(Low = c(x = 0.1, y = 11), 
               Moderate = c(x = 0.375, y = 11),
               High = c(x = 0.85, y = 11),
               spp = c(x = 0.1, y = 12.8))
sppLabSize <- 11
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.class <- eval(as.name(paste0("dat.class.", s)))
nests <- eval(as.name(paste0("nest.", s)))$HSI

plt <- plotDens(dat.plot, sampleHSIs = nests, dat.class, thresholds, binPntSize,
                classPntSize, tickLabSize = tickLabSize,
                BS = T, ylabel = NULL, xlabel = NULL)
plt <- plt +
  annotate("text", x = labxy[1, 1], y = labxy[1, 2], label = "Low", size = classLabSize) +
  annotate("text", x = labxy[2, 1], y = labxy[2, 2], label = "Moderate", size = classLabSize) +
  annotate("text", x = labxy[3, 1], y = labxy[3, 2], label = "High", size = classLabSize) +
  annotate("text", x = labxy[4, 1], y = labxy[4, 2], label = s, size = sppLabSize)

assign(paste0("plt.", s), plt)
rm(plt)

##### HAWO #####
#___Inputs___#
s <- "HAWO"
thresholds <- c(0.3, 0.55) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 20

# Parameters for additional plot labels
classLabSize <- 10
labxy <- rbind(Low = c(x = 0.18, y = 13.5), 
               Moderate = c(x = 0.425, y = 13.5),
               High = c(x = 0.775, y = 13.5),
               spp = c(x = 0.15, y = 15))
sppLabSize <- 11
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.class <- eval(as.name(paste0("dat.class.", s)))
nests <- eval(as.name(paste0("nest.", s)))$HSI

plt <- plotDens(dat.plot, sampleHSIs = nests, dat.class, thresholds, binPntSize,
                classPntSize, tickLabSize = tickLabSize,
                BS = T, ylabel = NULL, xlabel = NULL)
plt <- plt +
  annotate("text", x = labxy[1, 1], y = labxy[1, 2], label = "Low", size = classLabSize) +
  annotate("text", x = labxy[2, 1], y = labxy[2, 2], label = "Moderate", size = classLabSize) +
  annotate("text", x = labxy[3, 1], y = labxy[3, 2], label = "High", size = classLabSize) +
  annotate("text", x = labxy[4, 1], y = labxy[4, 2], label = s, size = sppLabSize)

assign(paste0("plt.", s), plt)
rm(plt)

##### WHWO #####
#___Inputs___#
s <- "WHWO"
thresholds <- c(0.3, 0.63) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 20

# Parameters for additional plot labels
classLabSize <- 10
labxy <- rbind(Low = c(x = 0.18, y = 25), 
               Moderate = c(x = 0.465, y = 25),
               High = c(x = 0.9, y = 25),
               spp = c(x = 0.1, y = 28))
sppLabSize <- 11
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.class <- eval(as.name(paste0("dat.class.", s)))
nests <- eval(as.name(paste0("nest.", s)))$HSI

plt <- plotDens(dat.plot, sampleHSIs = nests, dat.class, thresholds, binPntSize,
                classPntSize, tickLabSize = tickLabSize,
                BS = T, ylabel = NULL, xlabel = NULL)
plt <- plt +
  annotate("text", x = labxy[1, 1], y = labxy[1, 2], label = "Low", size = classLabSize) +
  annotate("text", x = labxy[2, 1], y = labxy[2, 2], label = "Moderate", size = classLabSize) +
  annotate("text", x = labxy[3, 1], y = labxy[3, 2], label = "High", size = classLabSize) +
  annotate("text", x = labxy[4, 1], y = labxy[4, 2], label = s, size = sppLabSize)

assign(paste0("plt.", s), plt)
rm(plt)

##### NOFL #####
#___Inputs___#
s <- "NOFL"
thresholds <- c(0.68) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 20

# Parameters for additional plot labels
classLabSize <- 10
labxy <- rbind(Low = c(x = 0.34, y = 23), 
               High = c(x = 0.75, y = 23),
               spp = c(x = 0.1, y = 26))
sppLabSize <-11
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.class <- eval(as.name(paste0("dat.class.", s)))
nests <- eval(as.name(paste0("nest.", s)))$HSI

plt <- plotDens(dat.plot, sampleHSIs = nests, dat.class, thresholds, binPntSize,
                classPntSize, tickLabSize = tickLabSize,
                BS = T, ylabel = NULL, xlabel = NULL)
plt <- plt +
  annotate("text", x = labxy[1, 1], y = labxy[1, 2], label = "Low", size = classLabSize) +
  annotate("text", x = labxy[2, 1], y = labxy[2, 2], label = "High", size = classLabSize) +
  annotate("text", x = labxy[3, 1], y = labxy[3, 2], label = s, size = sppLabSize)

assign(paste0("plt.", s), plt)
rm(plt)

#### Final Plot ####
p <- ggdraw() + 
  draw_plot(plt.BBWO, x = 0.05, y = 0.525, width = 0.475, height = 0.475) +
  draw_plot(plt.HAWO, x = 0.525, y = 0.525, width = 0.475, height = 0.475) +
  draw_plot(plt.WHWO, x = 0.05, y = 0.05, width = 0.475, height = 0.475) +
  draw_plot(plt.NOFL, x = 0.525, y = 0.05, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("Density (nests per 100 ha)", "Habitat Suitability Index (HSI)"),
                  size=c(30, 30), x=c(0.01, 0.15), y=c(0.2, 0.05), angle = c(90, 0))
#p

save_plot("manuscript/Figure_CMBmod_HSI_NestDens.tiff", p, ncol = 3.3, nrow = 3.3, dpi = 1200)

