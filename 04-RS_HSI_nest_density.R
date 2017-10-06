############################################################################################
# Purpose: Relate observed nest densities for 4 species (BBWO, HAWO, WHWO, NOFL) with HSIs #
# from remote-sensed models.                                                               #
############################################################################################

library(foreign)
library(dplyr)
library(ggplot2)
library(cowplot)

#____________________________________ Inputs _______________________________________#
setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/PtBlue_Sierra/")
load("Data_compiled.RData") # Workspace containing data
NR_points <- "E:/GISData/PtBlue_Sierra/NR_points.dbf" # Needs to contain nest points with HSIs (randoms not used from here)
Grid_points <- "E:/GISData/PtBlue_Sierra/Trns50m_30m_grid.dbf" # Contains all 30-m grid points and associated HSIs within surveyed belt transects
within_50m <- TRUE # Indicates if densities are to be calculated within 50 m of transects (TRUE) or strictly within transects (FALSE)

spp <- c("BBWO", "HAWO", "WHWO", "NOFL")
sites <- c("ML", "CB", "CH")
no.yrs <- 4
area <- 9 # For manuscript: set area surveyed in units of 100 hectares

transects <- unique(substr(dat$SAMPLE_ID, 1, 4)) # Needed for bootstrapping at transect level
R <- 5000 # Desired number of bootstrapped samples for calculating Density uncertainty 
#___________________________________________________________________________________#

#___________ Retrieve data_____________#
dat.hsi <- read.dbf(NR_points, as.is = T) %>%
  tbl_df() %>%
  select(SAMPLE_ID, inPlot:HSI_NOFL)
for(sp in spp) {
  dat.spp <- eval(as.name(paste0("dat.", sp))) %>%
    left_join(dat.hsi, by = "SAMPLE_ID") %>%
    filter(Nest == 1) %>%
    select(SAMPLE_ID:inPlot50, get(paste0("HSI_", sp)))
  names(dat.spp)[length(dat.spp)] <- "HSI"
  ifelse(within_50m,
         dat.spp <- dat.spp %>% filter(inPlot50 == 1),
         dat.spp <- dat.spp %>% filter(inPlot == 1))
  assign(paste0("dat.", sp), dat.spp)
}
rm(dat.hsi, dat.spp, sp)

dat.grid <- read.dbf(Grid_points, as.is = T) %>%
  tbl_df() # Within 50 m of transects
#______________________________________#

#_____________ Get functions _____________#
source("Northern-Sierra-Woodpeckers/Functions.R") ## Script file with functions for model-fitting and plotting ##
source("Northern-Sierra-Woodpeckers/Functions_plotting.R") ## Additional functions used in this script only ##
#_________________________________________#

# Bins for densities #
if(within_50m) {
  bins <- calcBins(600, nrow(dat.grid), 4)
} else {
  bins <- calcBins(250, nrow(dat.grid), 4)
  }

# Tabulate values for plotting observed densities for moving window HSI bins #
for(s in spp) {
  nests <- eval(as.name(paste0("dat.", s)))$HSI
  grid <- getGridHSIs(s, dat.grid)
  assign(paste0("tab.", s),
         calcBinDensities(nests, grid, bins, area, no.yrs))
}
rm(s, nests, grid)

# Tabulate HSI class-specific values $
thresholds <- list(BBWO = c(0.41, 0.65), HAWO = c(0.51, 0.7),
                   WHWO = c(0.49, 0.75), NOFL = c(0.44, 0.6))
for(s in spp) {
  thrs <- thresholds[[s]]
  nests <- eval(as.name(paste0("dat.", s)))$HSI
  grid <- getGridHSIs(s, dat.grid)
  dat.class <- HSIClassDensities(nests, grid, thrs, area)

  # Add bootstrapped CIs #
  dat.class$Perc <- (((dat.class$Density) / sum(dat.class$Density)) *
    100) %>% round
  dat.class <- dat.class %>% HSIClassAddBS(dat.nest = eval(as.name(paste0("dat.", s))) %>%
                                           filter(Nest == 1),
                                         dat.grid, transects, thrs, area, R)
  assign(paste0("dat.class.", s), dat.class)
}
rm(thrs, dat.class, nests, grid)

#__________________________________ PLOTTING (manuscript version) ________________________________________#
##### BBWO #####
#___Inputs___#
s <- "BBWO"
thresholds <- c(0.41, 0.65) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15
classLabSize <- 4
sppLabSize <- 6
# Coordinates for plot labels
labxy <- rbind(maxSSS = c(x = 0.39, y = 9.7),
               Low = c(x = 0.2, y = 10.5), 
               Moderate = c(x = 0.53, y = 10.5),
               High = c(x = 0.75, y = 10.5),
               spp = c(x = 0.1, y = 10.8))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.class <- eval(as.name(paste0("dat.class.", s)))

plt <- plotNestDens(dat.plot, nests = eval(as.name(paste0("dat.", s)))$HSI,
                    dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize, BS = T)
assign(paste0("plt.", s), plt)
rm(plt)

##### HAWO #####
#___Inputs___#
s <- "HAWO"
thresholds <- c(0.51, 0.7) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15
classLabSize <- 4
sppLabSize <- 6
labxy <- rbind(maxSSS = c(x = 0.49, y = 12.7),
               Low = c(x = 0.3, y = 13.3), 
               Moderate = c(x = 0.6, y = 13.3),
               High = c(x = 0.77, y = 13.3),
               spp = c(x = 0.17, y = 13.6))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.class <- eval(as.name(paste0("dat.class.", s)))

plt <- plotNestDens(dat.plot, nests = eval(as.name(paste0("dat.", s)))$HSI,
                    dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize, BS = T)
assign(paste0("plt.", s), plt)
rm(plt)

##### WHWO #####
#___Inputs___#
s <- "WHWO"
thresholds <- c(0.49, 0.75) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15
classLabSize <- 4
sppLabSize <- 6
labxy <- rbind(maxSSS = c(x = 0.47, y = 12.5),
               Low = c(x = 0.25, y = 13.7), 
               Moderate = c(x = 0.63, y = 13.7),
               High = c(x = 0.87, y = 13.7),
               spp = c(x = 0.13, y = 14))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.class <- eval(as.name(paste0("dat.class.", s)))

plt <- plotNestDens(dat.plot, nests = eval(as.name(paste0("dat.", s)))$HSI,
                    dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize, BS = T)
assign(paste0("plt.", s), plt)
rm(plt)

##### NOFL #####
#___Inputs___#
s <- "NOFL"
thresholds <- c(0.44, 0.6) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15
classLabSize <- 4
sppLabSize <- 6
labxy <- rbind(maxSSS = c(x = 0.42, y = 16.5),
               Low = c(x = 0.25, y = 17.5), 
               Moderate = c(x = 0.52, y = 17.5),
               High = c(x = 0.75, y = 17.5),
               spp = c(x = 0.22, y = 18))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.class <- eval(as.name(paste0("dat.class.", s)))

plt <- plotNestDens(dat.plot, nests = eval(as.name(paste0("dat.", s)))$HSI,
                    dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize, BS = T)
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

save_plot("manuscript/Figure_HSI_NestDens.tiff", p, ncol = 3.3, nrow = 3.3, dpi = 100)




#______________________________ PLOTTING (GIS tool manual version) _______________________________#
##### BBWO #####
#___Inputs___#
s <- "BBWO"
thresholds <- c(0.41, 0.65) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15
classLabSize <- 4
sppLabSize <- 6
# Coordinates for plot labels
labxy <- rbind(maxSSS = c(x = 0.39, y = 9.7 * (1000 / 247.105)),
               Low = c(x = 0.2, y = 10.5 * (1000 / 247.105)), 
               Moderate = c(x = 0.53, y = 10.5 * (1000 / 247.105)),
               High = c(x = 0.75, y = 10.5 * (1000 / 247.105)),
               spp = c(x = 0.1, y = 10.8 * (1000 / 247.105)))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

dat.class <- eval(as.name(paste0("dat.class.", s)))
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95lo <- dat.class$Dens95lo * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95hi <- dat.class$Dens95hi * (1000 / 247.105) # Rescale to nests per 1000 acres

plt <- plotNestDens(dat.plot, nests = eval(as.name(paste0("dat.", s)))$HSI,
                    dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize, BS = T)

assign(paste0("plt.", s), plt)
rm(plt)

##### HAWO #####
#___Inputs___#
s <- "HAWO"
thresholds <- c(0.51, 0.7) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15
classLabSize <- 4
sppLabSize <- 6
labxy <- rbind(maxSSS = c(x = 0.49, y = 12.7 * (1000 / 247.105)),
               Low = c(x = 0.3, y = 13.3 * (1000 / 247.105)), 
               Moderate = c(x = 0.6, y = 13.3 * (1000 / 247.105)),
               High = c(x = 0.77, y = 13.3 * (1000 / 247.105)),
               spp = c(x = 0.17, y = 13.6 * (1000 / 247.105)))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

dat.class <- eval(as.name(paste0("dat.class.", s)))
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95lo <- dat.class$Dens95lo * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95hi <- dat.class$Dens95hi * (1000 / 247.105) # Rescale to nests per 1000 acres

plt <- plotNestDens(dat.plot, nests = eval(as.name(paste0("dat.", s)))$HSI,
                    dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize, BS = T)

assign(paste0("plt.", s), plt)
rm(plt)

##### WHWO #####
#___Inputs___#
s <- "WHWO"
thresholds <- c(0.49, 0.75) # Thresholds for low, moderate, and high suitability classes
binPntSize <- 2
classPntSize <- 5
tickLabSize <- 15
classLabSize <- 4
sppLabSize <- 6
labxy <- rbind(maxSSS = c(x = 0.47, y = 12.5 * (1000 / 247.105)),
               Low = c(x = 0.25, y = 13.7 * (1000 / 247.105)), 
               Moderate = c(x = 0.63, y = 13.7 * (1000 / 247.105)),
               High = c(x = 0.87, y = 13.7 * (1000 / 247.105)),
               spp = c(x = 0.13, y = 14 * (1000 / 247.105)))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

dat.class <- eval(as.name(paste0("dat.class.", s)))
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95lo <- dat.class$Dens95lo * (1000 / 247.105) # Rescale to nests per 1000 acres
dat.class$Dens95hi <- dat.class$Dens95hi * (1000 / 247.105) # Rescale to nests per 1000 acres

plt <- plotNestDens(dat.plot, nests = eval(as.name(paste0("dat.", s)))$HSI,
                    dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize, BS = T)

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


#ylab("Density (nests / 100 ha)") + xlab("Habitat suitability index")

#___________________ Tabulate class densities for manual _________________________#
class.densities <- array(NA, dim = c(3, 8, 3)) # Recepticle for class densities for table in GIS application tool manual
dimnames(class.densities) <-
  list(NULL, c("no_nests", "area", "Density", "Dens95lo", "Dens95hi", "perc", "perc95lo",
               "perc95hi"), spp[-which(spp == "NOFL")])
for(s in spp[1:3]) {
  class.densities[, c("no_nests", "Density", "Dens95lo", "Dens95hi", "perc", "perc95lo",
                      "perc95hi"), s] <-
    eval(as.name(paste0("dat.class.", s)))[, c("no_nests", "Density", "Dens95lo", "Dens95hi",
                                               "Perc", "perc95lo", "perc95hi")] %>%
    as.matrix
}

class.densities[, c("Density", "Dens95lo", "Dens95hi"), ] <-
  round(class.densities[, c("Density", "Dens95lo", "Dens95hi"), ] * (1000/247.105), digits = 1)
class.densities[, "area",] <- (class.densities[, "no_nests",] / class.densities[, "Density",]) * 1000
class.densities[, "area",] <- round(class.densities[, "area",], digits = 1)

cols <- c("Species", "Quantity", "low", "moderate", "high")
tab.dens <- matrix(NA, nrow = 24, ncol = length(cols), dimnames = list(NULL, cols)) %>%
  data.frame() %>% tbl_df %>%
  mutate(Species = Species %>% as.character) %>%
  mutate(Quantity = Species %>% as.character) %>%
  mutate_each_(funs(as.numeric(.)), names(.[, sapply(., is.logical)])) %>%
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
