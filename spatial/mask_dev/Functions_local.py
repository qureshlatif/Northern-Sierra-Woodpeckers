import os, arcpy, math
from arcpy.sa import *
from math import ceil
from math import pi

def checkPath(path):
    try: os.makedirs(path)
    except: pass

def makeWHRTYPE_mask(whrtype, vals): #(CWHR raster, list of values to accept in mask)
    rows = arcpy.SearchCursor(whrtype, "", "", "VALUE; WHRTYPE")
    rcls_table = []
    for r in rows:
        v = r.getValue("VALUE")
        c = r.getValue("WHRTYPE")
        if c in vals:
            rcls_table = rcls_table + [str(v) + ' 1']
        else:
            rcls_table = rcls_table + [str(v) + ' 0']
    rcls_table = ";".join(rcls_table)
    classR = arcpy.gp.Reclassify_sa(whrtype, "Value", rcls_table)
    return(classR)

# Clump and eliminate #
def ClumpAndEliminate(rast, minPixels):
    regionout = RegionGroup(rast, "FOUR", "WITHIN", "NO_LINK", "")
    reg_select = Con(Lookup(regionout, "Count") >= minPixels, regionout) 
    rast_CE = Nibble(rast, reg_select, "DATA_ONLY")
    return(rast_CE)

# Extend binary raster by specified 
def extendMask(rast, distance):
    rast_1km_sum = FocalStatistics(rast, NbrCircle(distance, "MAP"), "SUM", "DATA")
    pixel_width = int(str(arcpy.GetRasterProperties_management(rast, "CELLSIZEX")))
    maxVal = ceil((pi*(distance**2))/(pixel_width**2))
    mask = Reclassify(rast_1km_sum, "VALUE", RemapRange([[0, "NODATA"], [1, maxVal, 1]]))
    return(mask)
