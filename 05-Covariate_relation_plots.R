#############################################################################################
# Purpose: Relate observed nest densities for 4 species (BBWO, HAWO, WHWO, NOFL) with HSIs. #
# This version is for manuscript for peer-reviewed publication.                             #
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
#___________________________________________________________________________________#

##########################
# Remotely sensed models #
##########################

###_________BBWO, remotely sensed model___________###
mod <- loadObject("Model_RS_BBWO")
covs <- c("ccmort_loc", "blk_lndcc", "canhi_loc")
cov.names <- c("LocBurn", "LandBurn", "LocCC")

for(i in 1:length(covs)) {
  dat.plot <- data.frame(x = seq(min(dat[, covs[i]]), max(dat[, covs[i]]), length.out = 20))
  dat.plotz <- matrix(0, nrow = 20, ncol = length(covs)) %>% data.frame
  names(dat.plotz) <- covs
  dat.plotz[, covs[i]] <- (dat.plot$x - scale.factors[covs[i], "mean"]) /
    scale.factors[covs[i], "SD"]
  dat.plot$HSI <- predict(mod, dat.plotz, type = "response")
  
  plt <- ggplot(dat.plot, aes(x, HSI)) +
    geom_line(size = 3) +
    ylim(0,1) +
    ylab(NULL) + xlab(cov.names[i]) +
    theme(axis.text.x=element_text(size=15)) +
    theme(axis.text.y=element_text(size=15)) +
    theme(axis.title.x=element_text(size=20)) +
    theme(axis.title.y=element_text(size=20))
  assign(paste0("plt.", cov.names[i]), plt)
}

# Put plots together #

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(plt.LocBurn, x = 0.05, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LandBurn, x = 0.525, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LocCC, x = 0.05, y = 0, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("Habitat Suitability Index (HSI)", "BBWO, remotely sensed model"),
                  size=c(27, 30), x=c(0.01, 0.05), y=c(0.06, 1), angle = c(90, 0))

#p

save_plot("manuscript/Figure_BBWO_RS_relations.tiff", p, ncol = 3, nrow = 3, dpi = 200)
###_______________________________________________###

###_________HAWO, remotely sensed model___________###
mod <- loadObject("Model_RS_HAWO")
covs <- c("ccmort_loc", "canhi_lnd", "sizlrg_loc")
cov.names <- c("LocBurn", "LandCC", "LocSizeLrg")

for(i in 1:length(covs)) {
  dat.plot <- data.frame(x = seq(min(dat[, covs[i]]), max(dat[, covs[i]]), length.out = 20))
  dat.plotz <- matrix(0, nrow = 20, ncol = length(covs)) %>% data.frame
  names(dat.plotz) <- covs
  dat.plotz[, covs[i]] <- (dat.plot$x - scale.factors[covs[i], "mean"]) /
    scale.factors[covs[i], "SD"]
  dat.plot$HSI <- predict(mod, dat.plotz, type = "response")
  
  plt <- ggplot(dat.plot, aes(x, HSI)) +
    geom_line(size = 3) +
    ylim(0,1) +
    ylab(NULL) + xlab(cov.names[i]) +
    theme(axis.text.x=element_text(size=15)) +
    theme(axis.text.y=element_text(size=15)) +
    theme(axis.title.x=element_text(size=20)) +
    theme(axis.title.y=element_text(size=20))
  assign(paste0("plt.", cov.names[i]), plt)
}

# Put plots together #

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(plt.LocBurn, x = 0.05, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LandCC, x = 0.525, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LocSizeLrg, x = 0.05, y = 0, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("Habitat Suitability Index (HSI)", "HAWO, remotely sensed model"),
                  size=c(27, 30), x=c(0.01, 0.05), y=c(0.06, 1), angle = c(90, 0))

#p

save_plot("manuscript/Figure_HAWO_RS_relations.tiff", p, ncol = 3, nrow = 3, dpi = 200)
###_______________________________________________###

###_________WHWO, remotely sensed model___________###
mod <- loadObject("Model_RS_WHWO")
covs <- c("ccmort_loc", "blk_lndcc")
cov.names <- c("LocBurn", "LandBurn")

for(i in 1:length(covs)) {
  dat.plot <- data.frame(x = seq(min(dat[, covs[i]]), max(dat[, covs[i]]), length.out = 20))
  dat.plotz <- matrix(0, nrow = 20, ncol = length(covs)) %>% data.frame
  names(dat.plotz) <- covs
  dat.plotz[, covs[i]] <- (dat.plot$x - scale.factors[covs[i], "mean"]) /
    scale.factors[covs[i], "SD"]
  dat.plot$HSI <- predict(mod, dat.plotz, type = "response")
  
  plt <- ggplot(dat.plot, aes(x, HSI)) +
    geom_line(size = 3) +
    ylim(0,1) +
    ylab(NULL) + xlab(cov.names[i]) +
    theme(axis.text.x=element_text(size=15)) +
    theme(axis.text.y=element_text(size=15)) +
    theme(axis.title.x=element_text(size=20)) +
    theme(axis.title.y=element_text(size=20))
  assign(paste0("plt.", cov.names[i]), plt)
}

# Put plots together #

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(plt.LocBurn, x = 0.05, y = 0, width = 0.475, height = 0.9) +
  draw_plot(plt.LandBurn, x = 0.525, y = 0, width = 0.475, height = 0.9) +
  draw_plot_label(label = c("Habitat Suitability Index (HSI)", "WHWO, remotely sensed model"),
                  size=c(27, 30), x=c(0.01, 0.3), y=c(0.08, 1), angle = c(90, 0), hjust = c(0, 0))

#p

save_plot("manuscript/Figure_WHWO_RS_relations.tiff", p, ncol = 3, nrow = 1.5, dpi = 200)
###_______________________________________________###

###_________NOFL, remotely sensed model___________###
mod <- loadObject("Model_RS_NOFL")
covs <- c("slope", "ccmort_loc", "sizlrg_loc", "sizlrg_lnd")
cov.names <- c("SLOPE", "LocBurn", "LocSizeLrg", "LandSizeLrg")

for(i in 1:length(covs)) {
  dat.plot <- data.frame(x = seq(min(dat[, covs[i]]), max(dat[, covs[i]]), length.out = 20))
  dat.plotz <- matrix(0, nrow = 20, ncol = length(covs)) %>% data.frame
  names(dat.plotz) <- covs
  dat.plotz[, covs[i]] <- (dat.plot$x - scale.factors[covs[i], "mean"]) /
    scale.factors[covs[i], "SD"]
  dat.plot$HSI <- predict(mod, dat.plotz, type = "response")
  
  plt <- ggplot(dat.plot, aes(x, HSI)) +
    geom_line(size = 3) +
    ylim(0,1) +
    ylab(NULL) + xlab(cov.names[i]) +
    theme(axis.text.x=element_text(size=15)) +
    theme(axis.text.y=element_text(size=15)) +
    theme(axis.title.x=element_text(size=20)) +
    theme(axis.title.y=element_text(size=20))
  assign(paste0("plt.", cov.names[i]), plt)
}

# Put plots together #

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(plt.SLOPE, x = 0.05, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LocBurn, x = 0.525, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LocSizeLrg, x = 0.05, y = 0, width = 0.475, height = 0.475) +
  draw_plot(plt.LandSizeLrg, x = 0.525, y = 0, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("Habitat Suitability Index (HSI)", "NOFL, remotely sensed model"),
                  size=c(27, 30), x=c(0.01, 0.05), y=c(0.06, 1), angle = c(90, 0))

#p

save_plot("manuscript/Figure_NOFL_RS_relations.tiff", p, ncol = 3, nrow = 3, dpi = 200)
###_______________________________________________###

######################
# Cobmination models #
######################

###_________BBWO, combination model___________###
mod <- loadObject("Model_CMB_BBWO")
covs <- c("ccmort_loc", "blk_lndcc", "DBH", "SnagDens_23to50")
cov.names <- c("LocBurn", "LandBurn", "DBH", "SnagDens23to50")

for(i in 1:length(covs)) {
  dat.plot <- data.frame(x = seq(min(dat[, covs[i]]), max(dat[, covs[i]]), length.out = 20))
  dat.plotz <- matrix(0, nrow = 20, ncol = length(covs)) %>% data.frame
  names(dat.plotz) <- covs
  dat.plotz[, covs[i]] <- (dat.plot$x - scale.factors[covs[i], "mean"]) /
    scale.factors[covs[i], "SD"]
  dat.plot$HSI <- predict(mod, dat.plotz, type = "response")
  
  plt <- ggplot(dat.plot, aes(x, HSI)) +
    geom_line(size = 3) +
    ylim(0,1) +
    ylab(NULL) + xlab(cov.names[i]) +
    theme(axis.text.x=element_text(size=15)) +
    theme(axis.text.y=element_text(size=15)) +
    theme(axis.title.x=element_text(size=20)) +
    theme(axis.title.y=element_text(size=20))
  assign(paste0("plt.", cov.names[i]), plt)
}

# Put plots together #

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(plt.LocBurn, x = 0.05, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LandBurn, x = 0.525, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.DBH, x = 0.05, y = 0, width = 0.475, height = 0.475) +
  draw_plot(plt.SnagDens23to50, x = 0.525, y = 0, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("Habitat Suitability Index (HSI)", "BBWO, combination model"),
                  size=c(27, 30), x=c(0.01, 0.15), y=c(0.06, 1), angle = c(90, 0))

#p

save_plot("manuscript/Figure_BBWO_CMB_relations.tiff", p, ncol = 3, nrow = 3, dpi = 200)
###_______________________________________________###

###_________HAWO, combination model___________###
mod <- loadObject("Model_CMB_HAWO")
covs <- c("ccmort_loc", "sizlrg_loc", "DBH", "BRKN")
cov.names <- c("LocBurn", "LocSizeLrg", "DBH", "BRKN")
categorical <- which(!covs %in% row.names(scale.factors))

for(i in 1:length(covs)) {
  if(i %in% categorical) {
    dat.plot <- data.frame(x = c(0, 1))
    dat.plotz <- matrix(0, nrow = 2, ncol = length(covs)) %>% data.frame
    for(cv in categorical)
      dat.plotz[, cv] <- mean(dat[, covs[cv]] %>% as.matrix %>% as.numeric %>% mean)
    names(dat.plotz) <- covs
    dat.plotz[, covs[i]] <- dat.plot[, "x"]
  } else {
    dat.plot <- data.frame(x = seq(min(dat[, covs[i]]), max(dat[, covs[i]]), length.out = 20))
    dat.plotz <- matrix(0, nrow = 20, ncol = length(covs)) %>% data.frame
    for(cv in categorical)
      dat.plotz[, cv] <- mean(dat[, covs[cv]] %>% as.matrix %>% as.numeric %>% mean)
    names(dat.plotz) <- covs
    dat.plotz[, covs[i]] <- (dat.plot$x - scale.factors[covs[i], "mean"]) /
      scale.factors[covs[i], "SD"]
  } 
  dat.plot$HSI <- predict(mod, dat.plotz, type = "response")
  
  if(!i %in% categorical) {
    plt <- ggplot(dat.plot, aes(x, HSI)) +
      geom_line(size = 3) +
      ylim(0,1) +
      ylab(NULL) + xlab(cov.names[i]) +
      theme(axis.text.x=element_text(size=15)) +
      theme(axis.text.y=element_text(size=15)) +
      theme(axis.title.x=element_text(size=20)) +
      theme(axis.title.y=element_text(size=20))
  } else {
    plt <- ggplot(dat.plot, aes(x, HSI)) +
      geom_point(size = 10, shape = 16) +
      ylim(0,1) +
      scale_x_continuous(limits = c(-0.5, 1.5), breaks = c(0, 1)) +
      ylab(NULL) + xlab(cov.names[i]) +
      theme(axis.text.x=element_text(size=15)) +
      theme(axis.text.y=element_text(size=15)) +
      theme(axis.title.x=element_text(size=20)) +
      theme(axis.title.y=element_text(size=20))
  }
  assign(paste0("plt.", cov.names[i]), plt)
}

# Put plots together #

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(plt.LocBurn, x = 0.05, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LocSizeLrg, x = 0.525, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.DBH, x = 0.05, y = 0, width = 0.475, height = 0.475) +
  draw_plot(plt.BRKN, x = 0.525, y = 0, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("Habitat Suitability Index (HSI)", "HAWO, combination model"),
                  size=c(27, 30), x=c(0.01, 0.15), y=c(0.06, 1), angle = c(90, 0))

#p

save_plot("manuscript/Figure_HAWO_CMB_relations.tiff", p, ncol = 3, nrow = 3, dpi = 200)
###_______________________________________________###

###_________WHWO, combination model___________###
mod <- loadObject("Model_CMB_WHWO")
covs <- c("ccmort_loc", "blk_lndcc", "DBH", "BRKN")
cov.names <- c("LocBurn", "LandBurn", "DBH", "BRKN")
categorical <- which(!covs %in% row.names(scale.factors))

for(i in 1:length(covs)) {
  if(i %in% categorical) {
    dat.plot <- data.frame(x = c(0, 1))
    dat.plotz <- matrix(0, nrow = 2, ncol = length(covs)) %>% data.frame
    for(cv in categorical)
      dat.plotz[, cv] <- mean(dat[, covs[cv]] %>% as.matrix %>% as.numeric %>% mean)
    names(dat.plotz) <- covs
    dat.plotz[, covs[i]] <- dat.plot[, "x"]
  } else {
    dat.plot <- data.frame(x = seq(min(dat[, covs[i]]), max(dat[, covs[i]]), length.out = 20))
    dat.plotz <- matrix(0, nrow = 20, ncol = length(covs)) %>% data.frame
    for(cv in categorical)
      dat.plotz[, cv] <- mean(dat[, covs[cv]] %>% as.matrix %>% as.numeric %>% mean)
    names(dat.plotz) <- covs
    dat.plotz[, covs[i]] <- (dat.plot$x - scale.factors[covs[i], "mean"]) /
      scale.factors[covs[i], "SD"]
  } 
  dat.plot$HSI <- predict(mod, dat.plotz, type = "response")
  
  if(!i %in% categorical) {
    plt <- ggplot(dat.plot, aes(x, HSI)) +
      geom_line(size = 3) +
      ylim(0,1) +
      ylab(NULL) + xlab(cov.names[i]) +
      theme(axis.text.x=element_text(size=15)) +
      theme(axis.text.y=element_text(size=15)) +
      theme(axis.title.x=element_text(size=20)) +
      theme(axis.title.y=element_text(size=20))
  } else {
    plt <- ggplot(dat.plot, aes(x, HSI)) +
      geom_point(size = 10, shape = 16) +
      ylim(0,1) +
      scale_x_continuous(limits = c(-0.5, 1.5), breaks = c(0, 1)) +
      ylab(NULL) + xlab(cov.names[i]) +
      theme(axis.text.x=element_text(size=15)) +
      theme(axis.text.y=element_text(size=15)) +
      theme(axis.title.x=element_text(size=20)) +
      theme(axis.title.y=element_text(size=20))
  }
  assign(paste0("plt.", cov.names[i]), plt)
}

# Put plots together #

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(plt.LocBurn, x = 0.05, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.LocSizeLrg, x = 0.525, y = 0.475, width = 0.475, height = 0.475) +
  draw_plot(plt.DBH, x = 0.05, y = 0, width = 0.475, height = 0.475) +
  draw_plot(plt.BRKN, x = 0.525, y = 0, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("Habitat Suitability Index (HSI)", "WHWO, combination model"),
                  size=c(27, 30), x=c(0.01, 0.15), y=c(0.06, 1), angle = c(90, 0))

#p

save_plot("manuscript/Figure_WHWO_CMB_relations.tiff", p, ncol = 3, nrow = 3, dpi = 200)
###_______________________________________________###

###_________NOFL, combination model___________###
mod <- loadObject("Model_CMB_NOFL")
covs <- c("DBH", "BRKN")
cov.names <- c("DBH", "BRKN")
categorical <- which(!covs %in% row.names(scale.factors))

for(i in 1:length(covs)) {
  if(i %in% categorical) {
    dat.plot <- data.frame(x = c(0, 1))
    dat.plotz <- matrix(0, nrow = 2, ncol = length(covs)) %>% data.frame
    for(cv in categorical)
      dat.plotz[, cv] <- mean(dat[, covs[cv]] %>% as.matrix %>% as.numeric %>% mean)
    names(dat.plotz) <- covs
    dat.plotz[, covs[i]] <- dat.plot[, "x"]
  } else {
    dat.plot <- data.frame(x = seq(min(dat[, covs[i]]), max(dat[, covs[i]]), length.out = 500))
    dat.plotz <- matrix(0, nrow = 500, ncol = length(covs)) %>% data.frame
    for(cv in categorical)
      dat.plotz[, cv] <- mean(dat[, covs[cv]] %>% as.matrix %>% as.numeric %>% mean)
    names(dat.plotz) <- covs
    dat.plotz[, covs[i]] <- (dat.plot$x - scale.factors[covs[i], "mean"]) /
      scale.factors[covs[i], "SD"]
  } 
  dat.plot$HSI <- predict(mod, dat.plotz, type = "response")
  
  if(!i %in% categorical) {
    plt <- ggplot(dat.plot, aes(x, HSI)) +
      geom_line(size = 3) +
      ylim(0,1) +
      ylab(NULL) + xlab(cov.names[i]) +
      theme(axis.text.x=element_text(size=15)) +
      theme(axis.text.y=element_text(size=15)) +
      theme(axis.title.x=element_text(size=20)) +
      theme(axis.title.y=element_text(size=20))
  } else {
    plt <- ggplot(dat.plot, aes(x, HSI)) +
      geom_point(size = 10, shape = 16) +
      ylim(0,1) +
      scale_x_continuous(limits = c(-0.5, 1.5), breaks = c(0, 1)) +
      ylab(NULL) + xlab(cov.names[i]) +
      theme(axis.text.x=element_text(size=15)) +
      theme(axis.text.y=element_text(size=15)) +
      theme(axis.title.x=element_text(size=20)) +
      theme(axis.title.y=element_text(size=20))
  }
  assign(paste0("plt.", cov.names[i]), plt)
}

# Put plots together #

theme_set(theme_bw())
p <- ggdraw() + 
  draw_plot(plt.DBH, x = 0.05, y = 0, width = 0.475, height = 0.9) +
  draw_plot(plt.BRKN, x = 0.525, y = 0, width = 0.475, height = 0.9) +
  draw_plot_label(label = c("Habitat Suitability Index (HSI)", "WHWO, combination model"),
                  size=c(27, 30), x=c(0.01, 0.3), y=c(0.08, 1), angle = c(90, 0),
                  hjust = c(0, 0))

#p

save_plot("manuscript/Figure_NOFL_CMB_relations.tiff", p, ncol = 3, nrow = 1.5, dpi = 200)
###_______________________________________________###
