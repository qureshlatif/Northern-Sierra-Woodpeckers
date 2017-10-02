import os, arcpy, math
from arcpy.sa import *
from Functions_local import *

######## Inputs ########
base = "E:/GISData/Mdl_appl_tlbx/mask_dev/CA/" # set workspace #
arcpy.env.workspace = base
out_loc = "E:/GISData/Mdl_appl_tlbx/TOOLBOX/masks/NSierra/" # Where final masks are to be stored #

# CWHR geodatabases and polygons within these containing raw data
CWHRgdbs = ["ExistingVegGreatBasin1999_2009_v1.gdb/", "ExistingVegNorInterior1999_2009_v1.gdb/",
            "ExistingVegNorSierra2000_2014_v1.gdb/", "ExistingVegSouthSierra2000_2008_v1.gdb/"]

CWHRpolys = ["ExistingVegGreatBasin1999_2009_v1", "ExistingVegNorInterior1999_2009_v1",
             "ExistingVegNorSierra2000_2014_v1", "ExistingVegSouthSierra2000_2008_v1"]
########################

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = base + "/_scratch"
checkPath(arcpy.env.scratchWorkspace)

Reclass_list = []

for i in range(0, len(CWHRgdbs)):
    whrtype = arcpy.PolygonToRaster_conversion(in_features = base + CWHRgdbs[i] + CWHRpolys[i], value_field = "WHRTYPE",
                                               out_rasterdataset = base + CWHRgdbs[i] + "whrtype", cellsize = 30)
    whrtype = Raster(whrtype)
    # Reclassify chunk #
    vals = ['SMC', 'WFR', 'RFR']
    Reclass = makeWHRTYPE_mask(whrtype, vals)
    Reclass = Raster(Reclass)
    Reclass.save(base + CWHRgdbs[i] + "Reclass")
    Reclass_list = Reclass_list + [base + CWHRgdbs[i] + "Reclass"]
    del(Reclass, vals)

Reclass = arcpy.MosaicToNewRaster_management(input_rasters = Reclass_list,
                                   output_location = base,
                                   raster_dataset_name_with_extension = "Reclass.tif",
                                   number_of_bands = 1)

Reclass_CE = ClumpAndEliminate(Reclass, 3491) #Threshold currently set at 314 ha (3491 pixels).

mask = extendMask(Reclass_CE, 1000)
mask.save(out_loc + "mask")
