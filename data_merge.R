## Preliminary example file ##
# Author: Brent Campos
# ***Not in use***

##=====================================================================================================================
## LOAD NECESSARY PACKAGES

require(foreign)
require(plyr)


##=====================================================================================================================
## SUMMARIZE SNAG PLOTS AND COMBINE DATA INTO FLAT FILE WITH EXPORTED DATA

snags <- read.csv(file="//prbo.org/Data/Home/Petaluma/bcampos/Documents/PRBO/Analysis/PLAS/Habitat Suitability/For RMRS/snag_plots_for_analysis_20170425.csv")
pivots <- read.csv(file="//prbo.org/Data/Home/Petaluma/bcampos/Documents/PRBO/Analysis/PLAS/Habitat Suitability/For RMRS/random_and_nest_trees_for_analysis_20170425.csv")

## Summarize the snag plot data into plot-level metrics.
mean.na <- function(x) mean(x, na.rm=TRUE)
snagsum <- ddply(snags,.(SAMPLE_ID), function(x) apply(x[,names(snags) %in% c("DBH","TREE_HEIGHT","DECAY","SCORCH")],2,mean.na))

## Asks whether Plot was duplicated in this summarized dataframe. If so, that means a single plot was assigned to two plot types.
duplicated(snagsum$SAMPLE_ID)

## Create some variables for analyses
MEAN_SNAG_DBH <- ddply(snags,.(SAMPLE_ID), summarise, MEAN_SNAG_DBH = mean(DBH, na.rm=TRUE))
SD_SNAG_DBH <- ddply(snags,.(SAMPLE_ID), summarise, SD_SNAG_DBH = sd(DBH, na.rm=TRUE))
MED_SNAG_DECAY <- ddply(snags,.(SAMPLE_ID), summarise, MED_SNAG_DECAY = median(DECAY, na.rm=TRUE))
SNAG23_ <- ddply(snags,.(SAMPLE_ID), function(x) nrow(x[!(x$TREE_SPECIES=="NONE" | x$TREE_SPECIES=="" | is.na(x$TREE_SPECIES)),]))
names(SNAG23_) <- c("SAMPLE_ID","SNAG23_")
SNAG23_49 <- ddply(snags,.(SAMPLE_ID), function(x) nrow(x[( x$DBH >= 23 & x$DBH < 50 & !is.na(x$TREE_SPECIES) & !(x$TREE_SPECIES=="NONE") & !(x$TREE_SPECIES=="")),]))
names(SNAG23_49) <- c("SAMPLE_ID","SNAG23_49")
SNAG50_74 <- ddply(snags,.(SAMPLE_ID), function(x) nrow(x[( x$DBH >= 50 & x$DBH < 75 & !is.na(x$TREE_SPECIES) & !(x$TREE_SPECIES=="NONE") & !(x$TREE_SPECIES=="")),]))
names(SNAG50_74) <- c("SAMPLE_ID","SNAG50_74")
SNAG75_99 <- ddply(snags,.(SAMPLE_ID), function(x) nrow(x[( x$DBH >= 75 & x$DBH < 100 & !is.na(x$TREE_SPECIES) & !(x$TREE_SPECIES=="NONE") & !(x$TREE_SPECIES=="")),]))
names(SNAG75_99) <- c("SAMPLE_ID","SNAG75_99")
SNAG100_ <- ddply(snags,.(SAMPLE_ID), function(x) nrow(x[( x$DBH >= 100 & !is.na(x$TREE_SPECIES) & !(x$TREE_SPECIES=="NONE") & !(x$TREE_SPECIES=="")),]))
names(SNAG100_) <- c("SAMPLE_ID","SNAG100_")
snagss <- join_all(list(MEAN_SNAG_DBH, SD_SNAG_DBH, MED_SNAG_DECAY, SNAG23_, SNAG23_49, SNAG50_74, SNAG75_99, SNAG100_), c("SAMPLE_ID"), match="all")
View(snagss)

## Merge the summarized snag plot data with the random snag and nest cavity data.
data <- merge(snagss, pivots, by=c("SAMPLE_ID"), all=T)
data <- data[order(data$YEAR, as.character(data$SAMPLE_ID)),]

write.csv(data, file="//prbo.org/Data/Home/Petaluma/bcampos/Documents/PRBO/Analysis/PLAS/Habitat Suitability/For RMRS/summarized_data_for_analysis_20170425.csv", row.names=FALSE)
