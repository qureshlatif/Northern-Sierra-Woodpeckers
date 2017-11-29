#############################################################################################
# Purpose: Relate observed nest densities for 4 species (BBWO, HAWO, WHWO, NOFL) with HSIs. #
# This version is for GIS application tool user's manual.                                   #
#############################################################################################

library(foreign)
library(R.utils)
#devtools::install_github("qureshlatif/WoodpeckerHSI") # Run this upon first use.
# Note: You might get some disconcerting warnings upon first install.
#To avoid them, restart R after installing package.
library(WoodpeckerHSI)

#____________________________________ Inputs _______________________________________#
setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")
load("Data_compiled.RData") # Workspace containing data
NR_points <- "E:/GISData/PtBlue_Sierra/NR_points.dbf" # Needs to contain nest points with HSIs (randoms not used from here)
Grid_points <- "E:/GISData/PtBlue_Sierra/Trns50m_30m_grid.dbf" # Contains all 30-m grid points and associated HSIs within surveyed belt transects

spp <- c("BBWO", "HAWO", "WHWO", "NOFL")
sites <- c("ML", "CB", "CH")
area <- 9 * 4 #Area surveyed in 100 ha units X study duration.

transects <- unique(substr(dat$SAMPLE_ID, 1, 4)) # Needed for bootstrapping at transect level
R <- 5000 # Desired number of bootstrapped samples for calculating Density uncertainty 
thresholds <- list(BBWO = c(0.35, 0.65), HAWO = c(0.25, 0.51),
                   WHWO = c(0.4, 0.7), NOFL = c(0.3, 0.55))
#___________________________________________________________________________________#

#___________ Retrieve data_____________#
dat.hsi <- read.dbf(NR_points, as.is = T) %>%
  tbl_df() %>%
  select(SAMPLE_ID, inPlot:HSI_NOFL)
for(sp in spp) {
  dat.spp <- eval(as.name(paste0("dat.", sp))) %>%
    left_join(dat.hsi, by = "SAMPLE_ID") %>%
    filter(Nest == 1) %>%
    select(SAMPLE_ID:inPlot50, matches(paste0("HSI_", sp)))
  names(dat.spp)[length(dat.spp)] <- "HSI"
  dat.spp <- dat.spp %>% filter(inPlot50 == 1)
  assign(paste0("dat.", sp), dat.spp)
}
rm(dat.hsi, dat.spp, sp)

dat.grid <- read.dbf(Grid_points, as.is = T) %>%
  tbl_df() # Within 50 m of transects
#______________________________________#

# Bins for densities #
bins <- calcBins(1000, nrow(dat.grid), 3)

# Tabulate values for plotting observed densities for moving window HSI bins #
for(s in spp) {
  nests <- eval(as.name(paste0("dat.", s)))$HSI
  grid <- dat.grid %>%
    select(matches(paste0("HSI_", s)))
  names(grid)[1] <- "HSI"
  grid <- grid$HSI
  assign(paste0("tab.", s),
         calcBinDensities(nests, grid, bins, area))
}
rm(s, nests, grid)

# Retrieve cached tables with class densities #
dat.class.BBWO <- read.csv("Plot_cache_RS_BBWO_class.csv", header = T)
dat.class.HAWO <- read.csv("Plot_cache_RS_HAWO_class.csv", header = T)
dat.class.WHWO <- read.csv("Plot_cache_RS_WHWO_class.csv", header = T)
dat.class.NOFL <- read.csv("Plot_cache_RS_NOFL_class.csv", header = T)

#______________________________ PLOTTING (GIS tool manual version) _______________________________#
##### BBWO #####
#____ Inputs ____#
s <- "BBWO"
thresholds <- c(0.35, 0.65) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15

# Parameters for additional plot labels
classLabSize <- 5
labxy <- rbind(Low = c(x = 0.2, y = 10.7 * (1000 / 247.105)), 
               Moderate = c(x = 0.53, y = 10.7 * (1000 / 247.105)),
               High = c(x = 0.75, y = 10.7 * (1000 / 247.105)),
               spp = c(x = 0.1, y = 11 * (1000 / 247.105)))
sppLabSize <- 6
#________________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

dat.class <- eval(as.name(paste0("dat.class.", s)))
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95lo <- dat.class$Dens95lo * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95hi <- dat.class$Dens95hi * (1000 / 247.105) # Rescale to nests per 1000 acres

nests <- eval(as.name(paste0("dat.", s)))$HSI

plt <- plotDens(dat.plot, sampleHSIs = nests, dat.class, thresholds, binPntSize,
                classPntSize, tickLabSize,
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
thresholds <- c(0.25, 0.51) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15

# Parameters for additional plot labels
classLabSize <- 5
labxy <- rbind(Low = c(x = 0.15, y = 7 * (1000 / 247.105)), 
               Moderate = c(x = 0.35, y = 7 * (1000 / 247.105)),
               High = c(x = 0.7, y = 7 * (1000 / 247.105)),
               spp = c(x = 0.15, y = 7.3 * (1000 / 247.105)))
sppLabSize <- 6
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

dat.class <- eval(as.name(paste0("dat.class.", s)))
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95lo <- dat.class$Dens95lo * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95hi <- dat.class$Dens95hi * (1000 / 247.105) # Rescale to nests per 1000 acres

nests <- eval(as.name(paste0("dat.", s)))$HSI

plt <- plotDens(dat.plot, sampleHSIs = nests, dat.class, thresholds, binPntSize,
                classPntSize, tickLabSize,
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
thresholds <- c(0.4, 0.7) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15

# Parameters for additional plot labels
classLabSize <- 5
labxy <- rbind(Low = c(x = 0.25, y = 11.8 * (1000 / 247.105)), 
               Moderate = c(x = 0.55, y = 11.8 * (1000 / 247.105)),
               High = c(x = 0.85, y = 11.8 * (1000 / 247.105)),
               spp = c(x = 0.1, y = 12.1 * (1000 / 247.105)))
sppLabSize <- 6
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

dat.class <- eval(as.name(paste0("dat.class.", s)))
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95lo <- dat.class$Dens95lo * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95hi <- dat.class$Dens95hi * (1000 / 247.105) # Rescale to nests per 1000 acres

nests <- eval(as.name(paste0("dat.", s)))$HSI

plt <- plotDens(dat.plot, sampleHSIs = nests, dat.class, thresholds, binPntSize,
                classPntSize, tickLabSize,
                BS = T, ylabel = NULL, xlabel = NULL)
plt <- plt +
  annotate("text", x = labxy[1, 1], y = labxy[1, 2], label = "Low", size = classLabSize) +
  annotate("text", x = labxy[2, 1], y = labxy[2, 2], label = "Moderate", size = classLabSize) +
  annotate("text", x = labxy[3, 1], y = labxy[3, 2], label = "High", size = classLabSize) +
  annotate("text", x = labxy[4, 1], y = labxy[4, 2], label = s, size = sppLabSize)

assign(paste0("plt.", s), plt)
rm(plt)

#### Final Plot ####
p <- ggdraw() + 
  draw_plot(plt.BBWO, x = 0.05, y = 0.525, width = 0.475, height = 0.475) +
  draw_plot(plt.HAWO, x = 0.525, y = 0.525, width = 0.475, height = 0.475) +
  draw_plot(plt.WHWO, x = 0.05, y = 0.05, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("Density (nests per 1000 ac)", "Habitat Suitability Index (HSI)"),
                  size=c(30, 30), x=c(0.01, 0.15), y=c(0.2, 0.05), angle = c(90, 0))
#p

save_plot("F:/research stuff/FS_PostDoc/Model_application_tool/Figure_NSierra_HSIDens.tiff", p, ncol = 3.3, nrow = 3.3, dpi = 100)

#___________________ Tabulate class densities for manual _________________________#
class.densities <- array(NA, dim = c(3, 8, 3)) # Recepticle for class densities for table in GIS application tool manual
dimnames(class.densities) <-
  list(NULL, c("no_nests", "area", "Density", "Dens95lo", "Dens95hi", "perc", "perc95lo",
               "perc95hi"), spp[-which(spp == "NOFL")])
for(s in spp[1:3]) {
  class.densities[, c("no_nests", "Density", "Dens95lo", "Dens95hi", "perc", "perc95lo",
                      "perc95hi"), s] <-
    eval(as.name(paste0("dat.class.", s)))[, c("no_nests", "Density", "Dens95lo", "Dens95hi",
                                               "Perc", "Perc95lo", "Perc95hi")] %>%
    as.matrix
}

class.densities[, c("Density", "Dens95lo", "Dens95hi"), ] <-
  round(class.densities[, c("Density", "Dens95lo", "Dens95hi"), ] * (1000/247.105), digits = 1)
class.densities[, c("perc", "perc95lo", "perc95hi"), ] <- 
  round(class.densities[, c("perc", "perc95lo", "perc95hi"), ])
class.densities[, "area",] <- (class.densities[, "no_nests",] / class.densities[, "Density",]) * 1000
class.densities[, "area",] <- round(class.densities[, "area",], digits = 1)

cols <- c("Species", "Quantity", "low", "moderate", "high")
tab.dens <- matrix(NA, nrow = 24, ncol = length(cols), dimnames = list(NULL, cols)) %>%
  data.frame() %>% tbl_df %>%
  mutate(Species = Species %>% as.character) %>%
  mutate(Quantity = Species %>% as.character) %>%
  mutate_all(funs(as.numeric(.)), names(.[, sapply(., is.logical)])) %>%
  mutate(Species = rep(c("BBWO", "HAWO", "WHWO"), each = 8)) %>%
  mutate(Quantity = rep(c("No. nests", "Area surveyed (acres)",
                          "Density (per 1000 acres)", "Dens95lo", "Dens95hi", "PercNest",
                          "Perc95lo", "Perc95hi"), 3)) %>%
  mutate(low = class.densities[1,,] %>% as.numeric) %>%
  mutate(moderate = class.densities[2,,] %>% as.numeric) %>%
  mutate(high = class.densities[3,,] %>% as.numeric) %>%
  mutate(low = low %>% as.character) %>%
  mutate(moderate = moderate %>% as.character) %>%
  mutate(high = high %>% as.character)

for(s in spp) {
  tab.dens[which(tab.dens$Quantity == "Density (per 1000 acres)" & tab.dens$Species == s),
           c("low", "moderate", "high")] <-
    paste0(tab.dens[which(tab.dens$Quantity == "Density (per 1000 acres)" & tab.dens$Species == s),
                    c("low", "moderate", "high")],
           " (",
           tab.dens[which(tab.dens$Quantity == "Dens95lo" & tab.dens$Species == s),
                    c("low", "moderate", "high")],
           ",",
           tab.dens[which(tab.dens$Quantity == "Dens95hi" & tab.dens$Species == s),
                    c("low", "moderate", "high")],
           ")")
  tab.dens[which(tab.dens$Quantity == "PercNest" & tab.dens$Species == s),
           c("low", "moderate", "high")] <-
    paste0(tab.dens[which(tab.dens$Quantity == "PercNest" & tab.dens$Species == s),
                    c("low", "moderate", "high")],
           " (",
           tab.dens[which(tab.dens$Quantity == "Perc95lo" & tab.dens$Species == s),
                    c("low", "moderate", "high")],
           ",",
           tab.dens[which(tab.dens$Quantity == "Perc95hi" & tab.dens$Species == s),
                    c("low", "moderate", "high")],
           ")")
}
tab.dens <- tab.dens[-which(tab.dens$Quantity %in% c("Dens95lo", "Dens95hi", "Perc95lo", "Perc95hi")), ]

write.csv(tab.dens, "Class_densities_manual.csv", row.names = F)
