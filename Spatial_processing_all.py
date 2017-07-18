#####################################################################################################
# Title: Remotely sensed predictors                                                                 #
# Content: Steps taken to compile predictors for post-fire woodpecerk habitat suitability models    #
# Conducted steps in ArcGIS. Then copied Python snippets from geoprocessing results window to here. #
#####################################################################################################

## Generate nest and random point file from spatial coordinates file ##
#*** Done in ArcMap. Snippet not copied. ***#

## Fire perimeters ##
# Combine fire perimeters into one shapefile
arcpy.Union_analysis(in_features="E:\GISData\PtBlue_Sierra\MTBS\Chips_ca4001012127720120729\ca4001012127720120729_20110801_20130721_burn_bndy.shp #;E:\GISData\PtBlue_Sierra\MTBS\Cub_ca4019412149120080621\ca4019412149120080621_20070721_20090710_burn_bndy.shp #;E:\GISData\PtBlue_Sierra\MTBS\Moonlight_ca4022012073620070903\ca4022012073620070903_20070705_20080707_burn_bndy.shp #", out_feature_class="E:/GISData/PtBlue_Sierra/fire_perim.shp", join_attributes="ONLY_FID", cluster_tolerance="", gaps="GAPS")

# Buffer fire perimeters by 1500m
#*** Done in ArcMap. Snippet not copied. ***#

## Topography ##
# Calculate slope from DEM
arcpy.gp.Slope_sa("E:/GISData/PtBlue_Sierra/topo/elev", "E:/GISData/PtBlue_Sierra/topo/slope.tif", "PERCENT_RISE", "8.99994364370288E-06")

# Re-project topographic slope & aspect
arcpy.ProjectRaster_management(in_raster="E:/GISData/PtBlue_Sierra/topo/aspect", out_raster="E:/GISData/PtBlue_Sierra/topo/aspect_n83.tif", out_coor_system="PROJCS['WGS_1984_UTM_Zone_10N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',-123.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]", resampling_type="NEAREST", cell_size="27.6663645096706 27.6663645096705", geographic_transform="WGS_1984_(ITRF00)_To_NAD_1983", Registration_Point="", in_coor_system="GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],VERTCS['Unknown VCS',VDATUM['Unknown'],PARAMETER['Vertical_Shift',0.0],PARAMETER['Direction',1.0],UNIT['Meter',1.0]]")
arcpy.ProjectRaster_management(in_raster="E:/GISData/PtBlue_Sierra/topo/slope.tif", out_raster="E:/GISData/PtBlue_Sierra/topo/slope_wgs84z10.tif", out_coor_system="PROJCS['WGS_1984_UTM_Zone_10N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',-123.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]", resampling_type="NEAREST", cell_size="27.6663645096706 27.6663645096706", geographic_transform="WGS_1984_(ITRF00)_To_NAD_1983", Registration_Point="", in_coor_system="GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],VERTCS['Unknown VCS',VDATUM['Unknown'],PARAMETER['Vertical_Shift',0.0],PARAMETER['Direction',1.0],UNIT['Meter',1.0]]")


# Calculate sine and cosine aspect (Not sure that these snippets will run. May need to break down into steps compiled in model builder)
# Tools for these can be found in the "Misc_Tools.tbx"
arcpy.Model(Aspect_Raster="E:/GISData/PtBlue_Sierra/topo/aspect_wgs84z10.tif", Output_Raster__Cosine_Degrees_="E:/GISData/PtBlue_Sierra/topo/cosasp.tif")
arcpy.Model1(Aspect_Raster="E:/GISData/PtBlue_Sierra/topo/aspect_wgs84z10.tif", Output_Raster__Sine_Degrees_="E:/GISData/PtBlue_Sierra/topo/sinasp.tif")

# Mask topography to study area
arcpy.gp.ExtractByMask_sa("E:/GISData/PtBlue_Sierra/topo/slope_wgs84z10.tif", "E:/GISData/PtBlue_Sierra/fire_perims_1500m.shp", "E:/GISData/PtBlue_Sierra/predictors/slope.tif")
arcpy.gp.ExtractByMask_sa("E:/GISData/PtBlue_Sierra/topo/sinasp.tif", "E:/GISData/PtBlue_Sierra/fire_perims_1500m.shp", "E:/GISData/PtBlue_Sierra/predictors/sinasp.tif")
arcpy.gp.ExtractByMask_sa("E:/GISData/PtBlue_Sierra/topo/cosasp.tif", "E:/GISData/PtBlue_Sierra/fire_perims_1500m.shp", "E:/GISData/PtBlue_Sierra/predictors/cosasp.tif")

## Burn severity ##
# Mosaic dNBR layers
arcpy.MosaicToNewRaster_management(input_rasters="E:\GISData\PtBlue_Sierra\MTBS\Chips_ca4001012127720120729\ca4001012127720120729_20110801_20130721_dnbr.tif;E:\GISData\PtBlue_Sierra\MTBS\Cub_ca4019412149120080621\ca4019412149120080621_20070721_20090710_dnbr.tif;E:\GISData\PtBlue_Sierra\MTBS\Moonlight_ca4022012073620070903\ca4022012073620070903_20070705_20080707_dnbr.tif", output_location="E:/GISData/PtBlue_Sierra/MTBS", raster_dataset_name_with_extension="dnbr_mosaic.tif", coordinate_system_for_the_raster="", pixel_type="16_BIT_SIGNED", cellsize="30", number_of_bands="1", mosaic_method="LAST", mosaic_colormap_mode="FIRST")

# Re-project dNBR to WGS84 z10 to match points
arcpy.ProjectRaster_management(in_raster="E:/GISData/PtBlue_Sierra/MTBS/dnbr_mosaic.tif", out_raster="E:/GISData/PtBlue_Sierra/MTBS/dnbr_mosaic_wgs84z10.tif", out_coor_system="PROJCS['WGS_1984_UTM_Zone_10N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',-123.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]", resampling_type="NEAREST", cell_size="30 30", geographic_transform="WGS_1984_(ITRF00)_To_NAD_1983", Registration_Point="", in_coor_system="PROJCS['USA_Contiguous_Albers_Equal_Area_Conic_USGS_version',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Albers'],PARAMETER['false_easting',0.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',-96.0],PARAMETER['standard_parallel_1',29.5],PARAMETER['standard_parallel_2',45.5],PARAMETER['latitude_of_origin',23.0],UNIT['Meter',1.0]]")

# Mask dNBR to 1500m around fire perimeters
arcpy.gp.ExtractByMask_sa("E:/GISData/PtBlue_Sierra/MTBS/dnbr_mosaic_wgs84z10.tif", "E:/GISData/PtBlue_Sierra/fire_perims_1500m.shp", "E:/GISData/PtBlue_Sierra/MTBS/dnbr_mosaic_wgs84z10_mask.tif")

# Generate median dNBR for 3x3 neighborhood
arcpy.gp.FocalStatistics_sa("E:/GISData/PtBlue_Sierra/MTBS/dnbr_mosaic_wgs84z10_mask.tif", "E:/GISData/PtBlue_Sierra/predictors/dnbr_3x3md.tif", "Rectangle 3 3 CELL", "MEDIAN", "DATA")

# Proportion and high and low severity within 1km radius neighborhood
arcpy.gp.Reclassify_sa("E:/GISData/PtBlue_Sierra/MTBS/dnbr_mosaic_wgs84z10_mask.tif", "Value", "-32768 366 0;366 2000 1", "E:/GISData/PtBlue_Sierra/MTBS/high_severity.tif", "DATA")
arcpy.gp.Reclassify_sa("E:/GISData/PtBlue_Sierra/MTBS/dnbr_mosaic_wgs84z10_mask.tif", "Value", "-32768 176 1;176 2000 0", "E:/GISData/PtBlue_Sierra/MTBS/low_severity.tif", "DATA")
arcpy.gp.FocalStatistics_sa("E:/GISData/PtBlue_Sierra/MTBS/dnbr_gt366.tif", "E:/GISData/PtBlue_Sierra/predictors/sum_hisev.tif", "Circle 1000 MAP", "SUM", "DATA")
arcpy.gp.FocalStatistics_sa("E:/GISData/PtBlue_Sierra/MTBS/dnbr_lt176.tif", "E:/GISData/PtBlue_Sierra/predictors/sum_losev.tif", "Circle 1000 MAP", "SUM", "DATA")

## Veg (CWHR data) ##
# Convert polygon tiles to rasters for Type, Density, and Size (Done as batch but doesn't seem batch structure was retained in snippet, i.e., only one example shown here)
arcpy.PolygonToRaster_conversion(in_features="E:/GISData/PtBlue_Sierra/CWHR/Cub_&_Moonlight/EvegTile17A.shp", value_field="WHRTYPE", out_rasterdataset="E:/GISData/PtBlue_Sierra/CWHR/Cub_&_Moonlight/WHRType_Tile17A.tif", cell_assignment="CELL_CENTER", priority_field="NONE", cellsize="30")

# Reclassify to 0/1 categoricals
arcpy.gp.Reclassify_sa("E:/GISData/PtBlue_Sierra/CWHR/CanDensity_CWHRclass.tif", "Value", "1 4 0;5 1;6 0", "E:/GISData/PtBlue_Sierra/CWHR/CanDens_low.tif", "DATA")

# Generate predictors with focal statistics (batch structure not retained; 2 examples shown here) #
arcpy.gp.FocalStatistics_sa("E:/GISData/PtBlue_Sierra/CWHR/SizeClass_pole.tif", "E:/GISData/PtBlue_Sierra/predictors/size_pole_3x3.tif", "Rectangle 3 3 CELL", "SUM", "DATA")
arcpy.gp.FocalStatistics_sa("E:/GISData/PtBlue_Sierra/CWHR/CanDens_low.tif", "E:/GISData/PtBlue_Sierra/predictors/can_low_3x3.tif", "Rectangle 3 3 CELL", "SUM", "DATA")

# Generate Density, Size, and Type mosaics (batch structure not retained; only 1 example shown here) #
#***Note need to extract density, size, and forest type classification layers for each polygon file first and then mosaic.***
arcpy.MosaicToNewRaster_management(input_rasters="E:\GISData\PtBlue_Sierra\CWHR\Chips\WHRDensity_EVInterior.tif;E:\GISData\PtBlue_Sierra\CWHR\Chips\WHRDensity_EVSierra.tif;E:\GISData\PtBlue_Sierra\CWHR\Cub_&_Moonlight\WHRDensity_Tile17A.tif;E:\GISData\PtBlue_Sierra\CWHR\Cub_&_Moonlight\WHRDensity_Tile9A.tif;E:\GISData\PtBlue_Sierra\CWHR\Cub_&_Moonlight\WHRDensity_Tile9B.tif", output_location="E:/GISData/PtBlue_Sierra/CWHR", raster_dataset_name_with_extension="CanDensity_CWHRclass.tif", coordinate_system_for_the_raster="PROJCS['WGS_1984_UTM_Zone_10N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',-123.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]", pixel_type="8_BIT_UNSIGNED", cellsize="", number_of_bands="1", mosaic_method="LAST", mosaic_colormap_mode="FIRST")


# Extract predictor values to points
arcpy.gp.ExtractMultiValuesToPoints_sa("E:/GISData/PtBlue_Sierra/NR_points.shp", "E:\GISData\PtBlue_Sierra\predictors\fir_1km.tif fir_1km;E:\GISData\PtBlue_Sierra\predictors\pine_1km.tif pine_1km", "NONE")
