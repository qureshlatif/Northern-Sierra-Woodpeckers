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
#___________________________________________________________________________________#

#___________ Retrieve data_____________#
dat.hsi <- read.dbf(NR_points, as.is = T) %>%
  tbl_df() %>%
  select(SAMPLE_ID, inPlot:HSI_NOFL)

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
  nests <- getNestHSIs(s, dat.hsi, within_50m)
  grid <- getGridHSIs(s, dat.grid)
  assign(paste0("tab.", s),
         calcBinDensities(nests, grid, bins, area, no.yrs))
}
rm(s)

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
labxy <- rbind(maxSSS = c(x = 0.39, y = 6.7),
               Low = c(x = 0.2, y = 7.4), 
               Moderate = c(x = 0.53, y = 7.4),
               High = c(x = 0.75, y = 7.4),
               spp = c(x = 0.1, y = 7.9))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
nests <- getNestHSIs(s, dat.hsi, within_50m)
grid <- getGridHSIs(s, dat.grid)
dat.class <- HSIClassDensities(nests, grid, thresholds, area)

plt <- plotNestDens(dat.plot, dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize)
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
labxy <- rbind(maxSSS = c(x = 0.49, y = 7.7),
               Low = c(x = 0.3, y = 8.4), 
               Moderate = c(x = 0.6, y = 8.4),
               High = c(x = 0.77, y = 8.4),
               spp = c(x = 0.17, y = 8.9))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
nests <- getNestHSIs(s, dat.hsi, within_50m)
grid <- getGridHSIs(s, dat.grid)
dat.class <- HSIClassDensities(nests, grid, thresholds, area)

plt <- plotNestDens(dat.plot, dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize)
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
labxy <- rbind(maxSSS = c(x = 0.47, y = 8),
               Low = c(x = 0.25, y = 9.15), 
               Moderate = c(x = 0.63, y = 9.15),
               High = c(x = 0.87, y = 9.15),
               spp = c(x = 0.13, y = 9.65))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
nests <- getNestHSIs(s, dat.hsi, within_50m)
grid <- getGridHSIs(s, dat.grid)
dat.class <- HSIClassDensities(nests, grid, thresholds, area)

plt <- plotNestDens(dat.plot, dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize)
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
labxy <- rbind(maxSSS = c(x = 0.42, y = 10.5),
               Low = c(x = 0.25, y = 11.3), 
               Moderate = c(x = 0.52, y = 11.3),
               High = c(x = 0.75, y = 11.3),
               spp = c(x = 0.22, y = 11.8))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
nests <- getNestHSIs(s, dat.hsi, within_50m)
grid <- getGridHSIs(s, dat.grid)
dat.class <- HSIClassDensities(nests, grid, thresholds, area)

plt <- plotNestDens(dat.plot, dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize)
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
class.densities <- array(NA, dim = c(3, 4, 3)) # Recepticle for class densities for table in GIS application tool manual
dimnames(class.densities) <- list(NULL, c("no_nests", "area", "Density", "perc"),
                                  spp[-which(spp == "NOFL")])

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
labxy <- rbind(maxSSS = c(x = 0.39, y = 28),
               Low = c(x = 0.2, y = 29), 
               Moderate = c(x = 0.53, y = 29),
               High = c(x = 0.75, y = 29),
               spp = c(x = 0.1, y = 30))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

nests <- getNestHSIs(s, dat.hsi, within_50m)
grid <- getGridHSIs(s, dat.grid)
dat.class <- HSIClassDensities(nests, grid, thresholds, area)
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

class.densities[, c("no_nests", "Density"), s] <- dat.class[, c("no_nests", "Density")] %>% as.matrix

plt <- plotNestDens(dat.plot, dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize)
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
labxy <- rbind(maxSSS = c(x = 0.49, y = 33),
               Low = c(x = 0.3, y = 34), 
               Moderate = c(x = 0.6, y = 34),
               High = c(x = 0.77, y = 34),
               spp = c(x = 0.17, y = 35))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

nests <- getNestHSIs(s, dat.hsi, within_50m)
grid <- getGridHSIs(s, dat.grid)
dat.class <- HSIClassDensities(nests, grid, thresholds, area)
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

class.densities[, c("no_nests", "Density"), s] <- dat.class[, c("no_nests", "Density")] %>% as.matrix

plt <- plotNestDens(dat.plot, dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize)
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
labxy <- rbind(maxSSS = c(x = 0.47, y = 36),
               Low = c(x = 0.25, y = 37), 
               Moderate = c(x = 0.63, y = 37),
               High = c(x = 0.87, y = 37),
               spp = c(x = 0.13, y = 38))
#____________#

dat.plot <- eval(as.name(paste0("tab.", s)))
dat.plot$Density <- dat.plot$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

nests <- getNestHSIs(s, dat.hsi, within_50m)
grid <- getGridHSIs(s, dat.grid)
dat.class <- HSIClassDensities(nests, grid, thresholds, area)
dat.class$Density <- dat.class$Density * (1000 / 247.105) # Rescale to nests per 1000 acres

class.densities[, c("no_nests", "Density"), s] <- dat.class[, c("no_nests", "Density")] %>% as.matrix

plt <- plotNestDens(dat.plot, dat.class, thresholds, binPntSize, classPntSize, tickLabSize,
                    classLabSize, labxy, spp = s, sppLabSize)
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
class.densities[, "Density",] <- round(class.densities[, "Density",], digits = 1)
class.densities[, "area",] <- (class.densities[, "no_nests",] / class.densities[, "Density",]) * 1000
class.densities[, "area",] <- round(class.densities[, "area",], digits = 1)
class.densities[, "perc",] <-
  round(t(t(class.densities[, "no_nests",]) / colSums(class.densities[, "no_nests",]))*100)

cols <- c("Species", "Quantity", "low", "moderate", "high")
tab.dens <- matrix(NA, nrow = 12, ncol = length(cols), dimnames = list(NULL, cols)) %>%
  data.frame() %>% tbl_df %>%
  mutate(Species = Species %>% as.character) %>%
  mutate(Quantity = Species %>% as.character) %>%
  mutate_each_(funs(as.numeric(.)), names(.[, sapply(., is.logical)])) %>%
  mutate(Species = rep(c("BBWO", "HAWO", "WHWO"), each = 4)) %>%
  mutate(Quantity = rep(c("No. nests", "Area surveyed (acres)",
                          "Density (per 1000 acres)", "Percent of total nests"), 3)) %>%
  mutate(low = class.densities[1,,] %>% as.numeric) %>%
  mutate(moderate = class.densities[2,,] %>% as.numeric) %>%
  mutate(high = class.densities[3,,] %>% as.numeric)

write.csv(tab.dens, "Class_densities_manual.csv", row.names = F)
