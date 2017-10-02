#####################################################################################################
# Title: Remotely sensed predictors                                                                 #
# Content: Generates all spatial predictors used to model post-fire habitat for nesting woodpeckers #
#       at 3 wildfires (Cub, Chips, and Moonlight) in the Northern Sierras.                         #
# Author: Quresh S. Latif                                                                           #
#####################################################################################################

import arcpy
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
import os, sys, shutil

def checkPath(path):
    try: os.makedirs(path)
    except: pass

################### Inputs ###########################################
workspace = "E:/GISData/PtBlue_Sierra/" # base directory for all work
env.workspace = workspace
fire_perim = "fire_perims.shp" # Polygon file with fire perimeter
coord_reference = "NR_points.shp" # Set coordinate system to match this file.
# Need to add if/then statement here that sets coord_reference to fire_perim if not provided by user. 
dem = workspace + "/topo/elev" # Raw DEM that covers entire fire perimeter
burnsev = "RAVG/ccmort_mosaic.tif" # Raw burn severity file. Should be % canopy mortality derived from RdNBR
CWHR_polys = ["E:/GISData/PtBlue_Sierra/CWHR/Chips/ExistVegNorInterior1999_2009_clip.shp", # List of CWHR polygons - source for canopy and size class variables
              "E:/GISData/PtBlue_Sierra/CWHR/Chips/ExistVegNorSierra2000_2009_clip.shp",
              "E:/GISData/PtBlue_Sierra/CWHR/Cub_&_Moonlight/EvegTile16B.shp",
              "E:/GISData/PtBlue_Sierra/CWHR/Cub_&_Moonlight/EvegTile17A.shp",
              "E:/GISData/PtBlue_Sierra/CWHR/Cub_&_Moonlight/EvegTile9A.shp",
              "E:/GISData/PtBlue_Sierra/CWHR/Cub_&_Moonlight/EvegTile9B.shp"]
######################################################################

arcpy.env.overwriteOutput = True
scratch_ws = workspace + "_scratch"
checkPath(scratch_ws)
arcpy.env.scratchWorkspace = scratch_ws

#Set environmental settings for projection
arcpy.env.outputCoordinateSystem = arcpy.Describe(coord_reference).spatialReference

# Buffer fire perimeters by 1500m
SA = arcpy.Buffer_analysis(fire_perim, scratch_ws + "/SA.shp", 1000)

# Set output folder for all predictor layers
pred_out = workspace + "predictors/"
checkPath(pred_out)

###### Burn severity ######
# Generate median dNBR for 3x3 neighborhood
ccmort = arcpy.gp.FocalStatistics_sa(burnsev, scratch_ws + "/ccmort_loc0.tif", "Rectangle 3 3 CELL", "MEDIAN")
ccmort = arcpy.gp.ExtractByMask_sa(burnsev, SA, pred_out + "/ccmort_loc.tif")
arcpy.env.snapRaster = pred_out + "/ccmort_loc.tif" # Set everything else to snap to this.

# Proportion and high severity within 1km radius neighborhood
ccmort_mask = arcpy.gp.ExtractByMask_sa(burnsev, SA)
ccmort64 = arcpy.gp.Reclassify_sa(ccmort_mask, "Value", "-1 64 0;64 100 1")
ccmort64_1km_mean = arcpy.gp.FocalStatistics_sa(ccmort64, scratch_ws + "/hisev_1km.tif", "Circle 1000 MAP", "MEAN")
ccmort64_1km = Raster(ccmort64_1km_mean) * 100
ccmort64_1km = arcpy.gp.ExtractByMask_sa(ccmort64_1km, ccmort)
Raster(ccmort64_1km).save(pred_out + "/ccmort_gt64.tif")

###### Topography #######
# Calculate slope from DEM
slp = arcpy.gp.Slope_sa(dem, scratch_ws + "/slope.tif", "PERCENT_RISE", "8.99994364370288E-06")
slp = arcpy.gp.ExtractByMask_sa(slp, ccmort)
Raster(slp).save(pred_out + "/slope.tif")

###### Veg (CWHR data) ######
# Convert polygon tiles to rasters for Type, Density, and Size (Done as batch but doesn't seem batch structure was retained in snippet, i.e., only one example shown here)
HighDens_rasters = []; LrgSize_rasters = []
for p in range(0, len(CWHR_polys)):
    # Reclass density to high #
    dens_class_rast = arcpy.PolygonToRaster_conversion(in_features = CWHR_polys[p], value_field="WHRDENSITY", cellsize="30")
    dens_class_rast = Raster(dens_class_rast)
    rows = arcpy.SearchCursor(dens_class_rast, "", "", "VALUE; WHRDENSITY")
    rcls_table = []
    for r in rows:
        v = r.getValue("VALUE")
        c = r.getValue("WHRDENSITY")
        if c == 'M' or c == 'D':
            rcls_table = rcls_table + [str(v) + ' 1']
        else:
            rcls_table = rcls_table + [str(v) + ' 0']
    rcls_table = ";".join(rcls_table)
    highDensR = arcpy.gp.Reclassify_sa(dens_class_rast, "Value", rcls_table)
    HighDens_rasters = HighDens_rasters + [highDensR]
    # Reclass size to large #
    size_class_rast = arcpy.PolygonToRaster_conversion(in_features = CWHR_polys[p], value_field="WHRSIZE", cellsize="30")
    size_class_rast = Raster(size_class_rast)
    rows = arcpy.SearchCursor(size_class_rast, "", "", "VALUE; WHRSIZE")
    rcls_table = []
    for r in rows:
        v = r.getValue("VALUE")
        c = r.getValue("WHRSIZE")
        if c == '5':
            rcls_table = rcls_table + [str(v) + ' 1']
        else:
            rcls_table = rcls_table + [str(v) + ' 0']
    rcls_table = ";".join(rcls_table)
    lrgSizeR = arcpy.gp.Reclassify_sa(size_class_rast, "Value", rcls_table)
    LrgSize_rasters = LrgSize_rasters + [lrgSizeR]

HighDens = arcpy.MosaicToNewRaster_management(HighDens_rasters, scratch_ws, "HighDens.tif", number_of_bands = 1)
canhi_loc = arcpy.gp.FocalStatistics_sa(HighDens, scratch_ws + "/canhi_loc.tif", "Rectangle 3 3 CELL", "MEAN")
canhi_loc = Raster(canhi_loc)*100
canhi_loc = arcpy.gp.ExtractByMask_sa(canhi_loc, ccmort)
Raster(canhi_loc).save(pred_out + "canhi_loc.tif")
canhi_lnd = arcpy.gp.FocalStatistics_sa(HighDens, scratch_ws + "/canhi_lnd.tif", "Circle 1000 MAP", "MEAN")
canhi_lnd = Raster(canhi_lnd)*100
canhi_lnd = arcpy.gp.ExtractByMask_sa(canhi_lnd, ccmort)
Raster(canhi_lnd).save(pred_out + "canhi_lnd.tif")

LrgSize = arcpy.MosaicToNewRaster_management(LrgSize_rasters, scratch_ws, "LrgSize.tif", number_of_bands = 1)
sizlrg_loc = arcpy.gp.FocalStatistics_sa(LrgSize, scratch_ws + "/sizlrg_loc.tif", "Rectangle 3 3 CELL", "MEAN")
sizlrg_loc = Raster(sizlrg_loc)*100
sizlrg_loc = arcpy.gp.ExtractByMask_sa(sizlrg_loc, ccmort)
Raster(sizlrg_loc).save(pred_out + "sizlrg_loc.tif")
sizlrg_lnd = arcpy.gp.FocalStatistics_sa(LrgSize, scratch_ws + "/sizlrg_lnd.tif", "Circle 1000 MAP", "MEAN")
sizlrg_lnd = Raster(sizlrg_lnd)*100
sizlrg_lnd = arcpy.gp.ExtractByMask_sa(sizlrg_lnd, ccmort)
Raster(sizlrg_lnd).save(pred_out + "sizlrg_lnd.tif")

## Cleanup ##
del(sizlrg_lnd, sizlrg_loc, canhi_lnd, canhi_loc, LrgSize, HighDens, highDensR, lrgSizeR, v, c, p, r, rows, rcls_table, ccmort, slp, SA, ccmort64_1km,
    ccmort_mask, size_class_rast, dens_class_rast, HighDens_rasters, LrgSize_rasters, fire_perim, burnsev, ccmort64, ccmort64_1km_mean,
    coord_reference, dem, pred_out, Slope, arcpy.env.scratchWorkspace, CWHR_polys, env)
arcpy.Delete_management(scratch_ws + "/SA.shp")
arcpy.env.workspace = scratch_ws
filesToRemove = arcpy.ListRasters()
for f in filesToRemove:
    arcpy.Delete_management(f)
#shutil.rmtree(scratch_ws) # Some file is staying locked, so this command is not working.
    #Not seeing what's left in the environment that is locking the file.
