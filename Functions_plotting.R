##############################################################
# Contains additional functions used in plotting scripts.    #
# Store this file in workspace designated in primary scripts #
##############################################################

require(ggplot2)
require(cowplot)
require(dplyr)

getNestHSIs <- function(spp, dat.hsi, within_50m = TRUE) {
  nests <- eval(as.name(paste0("dat.", spp))) %>%
    left_join(dat.hsi, by = "SAMPLE_ID") %>%
    filter(Nest == 1)
  ifelse(within_50m,
         nests <- nests %>% filter(inPlot50 == 1),
         nests <- nests %>% filter(inPlot50 == 1))
  nests <- nests %>% select(eval(as.name(paste0("HSI_", s))))
  names(nests)[1] <- "HSI"
  nests <- nests$HSI
  return(nests)
}

getGridHSIs <- function(spp, dat.grid) { # dat.grid must contain species-specific HSI fields labeled as "HSI_[spp]"
  grid <- dat.grid %>%
    select(eval(as.name(paste0("HSI_", spp))))
  names(grid)[1] <- "HSI"
  grid <- grid$HSI
  return(grid)
}

calcBinDensities <- function(sampleHSI, bgroundHSI, bins, area, no.yrs) { # bins = output from calcBins
  cell.area <- area/length(bgroundHSI)
  HSI.sort <- sort(bgroundHSI)
  
  tab <- data.frame(HSI.md = numeric(length = length(bins$st)))
  tab$HSI.end <- tab$HSI.st <- 0
  tab$Density <- NA
  
  for(b in 1:length(bins$st)) {
    tab$HSI.md[b] <-
      mean(HSI.sort[which(bins$unitID %in% bins$st[b]:bins$end[b])])
    tab$HSI.st[b] <-
      min(HSI.sort[which(bins$unitID %in% bins$st[b]:bins$end[b])])
    tab$HSI.end[b] <-
      max(HSI.sort[which(bins$unitID %in% bins$st[b]:bins$end[b])])
    if(b == 1) tab$HSI.st[b] <- 0 # To make sure all nests with lowest HSIs are included
    if(b == nrow(tab)) tab$HSI.end[b] <- 1 # To make sure all nests with highest HSIs are included
    no.nests <- sum(sampleHSI >= tab$HSI.st[b] & sampleHSI <= tab$HSI.end[b])
    bin.area <- sum(bgroundHSI >= tab$HSI.st[b] & bgroundHSI <= tab$HSI.end[b])*cell.area
    tab$Density[b] <- (no.nests / bin.area) / no.yrs
  }
  return(tab)
}

HSIClassDensities <- function(nestHSIs, gridHSIs, thresholds, area) {
  dat.class <- data.frame(HSI.md = c(mean(gridHSIs[which(gridHSIs < thresholds[1])]),
                                     mean(gridHSIs[which(gridHSIs >= thresholds[1] & gridHSIs < thresholds[2])]),
                                     mean(gridHSIs[which(gridHSIs >= thresholds[2])])),
                          HSI.st = c(0, thresholds[1], thresholds[2]),
                          HSI.end = c(thresholds[1], thresholds[2], 1))
  dat.class$no_nests <- NA
  dat.class$Density <- NA
  for(r in 1:nrow(dat.class)) {
    dat.class$no_nests[r] <- sum(nestHSIs >= dat.class$HSI.st[r] & nestHSIs < dat.class$HSI.end[r])
    dat.class$Density[r] <-
      sum(nestHSIs >= dat.class$HSI.st[r] & nestHSIs < dat.class$HSI.end[r]) /
      (sum(gridHSIs >= dat.class$HSI.st[r] & gridHSIs < dat.class$HSI.end[r]) * (area / length(gridHSIs)) * no.yrs)
  }
  return(dat.class)
}

plotNestDens <- function(dat.plot, dat.class, thresholds, binPntSize = 2, classPntSize = 5,
                         tickLabSize = 15, classLabSize = 5, lab.params, spp, sppLabSize = 8) {
  theme_set(theme_bw())
  plt <- ggplot(dat.plot, aes(HSI.md, Density)) +
    geom_point(size = binPntSize, shape = 1) +
    geom_point(data = dat.class, aes(x = HSI.md, y = Density), size = classPntSize, shape = 16) +
    geom_rug(data = data.frame(nests), aes(x = nests, y = NULL), colour = "black", alpha = 0.3, size = 1) +
    geom_vline(xintercept = thresholds[1], linetype = "dashed") +
    geom_vline(xintercept = thresholds[2], linetype = "dashed") +
    ylab(NULL) + xlab(NULL) +
    theme(axis.text.x=element_text(size=tickLabSize)) +
    theme(axis.text.y=element_text(size=tickLabSize)) +
    annotate("text", x = lab.params[1, 1], y = lab.params[1, 2], label = "maxSSS", angle = 90) +
    annotate("text", x = lab.params[2, 1], y = lab.params[2, 2], label = "Low", size = classLabSize) +
    annotate("text", x = lab.params[3, 1], y = lab.params[3, 2], label = "Moderate", size = classLabSize) +
    annotate("text", x = lab.params[4, 1], y = lab.params[4, 2], label = "High", size = classLabSize) +
    annotate("text", x = lab.params[5, 1], y = lab.params[5, 2], label = spp, size = sppLabSize)
  return(plt)
}