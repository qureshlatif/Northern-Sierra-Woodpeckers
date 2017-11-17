#############################################################################
# Title: Remotely sensed predictors                                         #
# Content: Generates alternative percent canopy cover >60% predictors       #
#       at 3 wildfires (Cub, Chips, and Moonlight) in the Northern Sierras. #
# Author: Quresh S. Latif                                                   #
#############################################################################

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

# Set snap raster to local burn severity layer
arcpy.env.snapRaster = pred_out + "/ccmort_loc.tif" # Set everything else to snap to this.
ccmort = Raster(pred_out + "/ccmort_loc.tif")

###### Veg (CWHR data) ######
# Convert polygon tiles to rasters for Type, Density, and Size (Done as batch but doesn't seem batch structure was retained in snippet, i.e., only one example shown here)
HighDens_rasters = []
for p in range(0, len(CWHR_polys)):
    # Reclass density to high #
    dens_class_rast = arcpy.PolygonToRaster_conversion(in_features = CWHR_polys[p], value_field="WHRDENSITY", cellsize="30")
    dens_class_rast = Raster(dens_class_rast)
    rows = arcpy.SearchCursor(dens_class_rast, "", "", "VALUE; WHRDENSITY")
    rcls_table = []
    for r in rows:
        v = r.getValue("VALUE")
        c = r.getValue("WHRDENSITY")
        if c == 'D':
            rcls_table = rcls_table + [str(v) + ' 1']
        else:
            rcls_table = rcls_table + [str(v) + ' 0']
    rcls_table = ";".join(rcls_table)
    highDensR = arcpy.gp.Reclassify_sa(dens_class_rast, "Value", rcls_table)
    HighDens_rasters = HighDens_rasters + [highDensR]

HighDens = arcpy.MosaicToNewRaster_management(HighDens_rasters, scratch_ws, "HighDens.tif", number_of_bands = 1)
canhi_loc = arcpy.gp.FocalStatistics_sa(HighDens, scratch_ws + "/canhi_loc.tif", "Rectangle 3 3 CELL", "MEAN")
canhi_loc = Raster(canhi_loc)*100
canhi_loc = arcpy.gp.ExtractByMask_sa(canhi_loc, ccmort)
Raster(canhi_loc).save(pred_out + "canhi60_loc.tif")
canhi_lnd = arcpy.gp.FocalStatistics_sa(HighDens, scratch_ws + "/canhi_lnd.tif", "Circle 1000 MAP", "MEAN")
canhi_lnd = Raster(canhi_lnd)*100
canhi_lnd = arcpy.gp.ExtractByMask_sa(canhi_lnd, ccmort)
Raster(canhi_lnd).save(pred_out + "canhi60_lnd.tif")

## Cleanup ##
del(canhi_lnd, canhi_loc, HighDens, highDensR, v, c, p, r, rows, rcls_table, SA,
    dens_class_rast, HighDens_rasters, fire_perim,
    coord_reference, pred_out, Slope, arcpy.env.scratchWorkspace, CWHR_polys, env)
arcpy.Delete_management(scratch_ws + "/SA.shp")
arcpy.env.workspace = scratch_ws
filesToRemove = arcpy.ListRasters()
for f in filesToRemove:
    arcpy.Delete_management(f)
#shutil.rmtree(scratch_ws) # Some file is staying locked, so this command is not working.
    #Not seeing what's left in the environment that is locking the file.
