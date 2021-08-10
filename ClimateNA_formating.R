#####################################################
##                                                 ##
##  Redundancy Analysis (RDA): a Swiss-army knife  ##
##             for landscape genomics              ##
##                                                 ##
##              Capblancq & Forester               ##
##                                                 ##
##  Formatting ClimateNA rasters for WGS84         ##
##                                                 ##
#####################################################

# Required libraries
library(pegas)
library(raster)
library(rgdal)

# Set working directory
setwd("./Data/")


##############################
####   Period 1961-1990   ####

# Load 33 Climate NA bioclimatic variables downloaded here: https://adaptwest.databasin.org/pages/adaptwest-climatena/
ras_past <- stack(lapply(lapply(list.files("./NA_NORM_6190_Bioclim_ASCII/", pattern = ".asc$", full.names = T), read.asciigrid), raster))
names(ras_past) <- unlist(strsplit(unlist(lapply(strsplit(list.files("./NA_NORM_6190_Bioclim_ASCII", pattern = ".asc$", full.names = T), split = "./NA_NORM_6190_Bioclim_ASCII/"), function(x) x[2])), split = ".asc"))

# Set the projection based on the reference .asc file
crs(ras_past) <- crs(raster("./ClimateNA_Reference/ClimateNA_ID.asc"))

# Reproject the rasters into WGS84 coordinate system
ras_past_proj <- list()
for(i in 1:length(names(ras_past))){
  ras_past_proj[[i]] <- projectRaster(ras_past[[i]], crs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +no_defs")
}
ras_past_proj <- stack(ras_past_proj)

# Crop the rasters to the wanted extent
ras_past_proj <- crop(ras_past_proj, extent(-150, -103, 30, 65))

# Exporting the final rasters for later use
lapply(1:dim(ras_past_proj)[3], function(x) writeRaster(ras_past_proj[[x]], filename = as.character(paste(names(ras_past_proj)[x], "6190.img", sep = "_")), format = "HFA", overwrite = T))

# Clean the R environment
rm(list = ls())

##############################


##############################
####   Predictions 2050   ####

# Load Climate NA bioclimatic variables predicted for 2050 using RCP8.5 climate scenario and the ensemble projection of 15 CMIP 5 AOGCMs
# The rasters can be downloaded here: https://adaptwest.databasin.org/pages/adaptwest-climatena/
ras_2050 <- stack(lapply(lapply(list.files("./NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII/", pattern = ".asc$", full.names = T), read.asciigrid), raster))
names(ras_2050) <- unlist(strsplit(unlist(lapply(strsplit(list.files("./NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII", pattern = ".asc$", full.names = T), split = "./NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII/"), function(x) x[2])), split = ".asc"))

# Set the projection based on the reference .asc file
crs(ras_2050) <- crs(raster("./ClimateNA_Reference/ClimateNA_ID.asc"))

# Reproject the rasters into WGS84 coordinate system
ras_2050_proj <- list()
for(i in 1:length(names(ras_2050))){
  ras_2050_proj[[i]] <- projectRaster(ras_2050[[i]], crs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +no_defs")
}
ras_2050_proj <- stack(ras_2050_proj)

# Crop the rasters to the wanted extent
ras_2050_proj <- crop(ras_2050_proj, extent(-150, -103, 30, 65))

# Exporting the final rasters for later use
lapply(1:dim(ras_2050_proj)[3], function(x) writeRaster(ras_2050_proj[[x]], filename = as.character(paste(names(ras_2050_proj)[x], "2050_85.img", sep = "_")), format = "HFA", overwrite = T))

# Clean the R environment
rm(list = ls())

##############################


##############################
####   Predictions 2080   ####

# Load Climate NA bioclimatic variables predicted for 2080 using RCP8.5 climate scenario and the ensemble projection of 15 CMIP 5 AOGCMs
# The rasters can be downloaded here: https://adaptwest.databasin.org/pages/adaptwest-climatena/
ras_2080 <- stack(lapply(lapply(list.files("./NA_ENSEMBLE_rcp85_2080s_Bioclim_ASCII/", pattern = ".asc$", full.names = T), read.asciigrid), raster))
names(ras_2080) <- unlist(strsplit(unlist(lapply(strsplit(list.files("./NA_ENSEMBLE_rcp85_2080s_Bioclim_ASCII", pattern = ".asc$", full.names = T), split = "./NA_ENSEMBLE_rcp85_2080s_Bioclim_ASCII/"), function(x) x[2])), split = ".asc"))

# Set the projection based on the reference .asc file
crs(ras_2080) <- crs(raster("./ClimateNA_Reference/ClimateNA_ID.asc"))

# Reproject the rasters into WGS84 coordinate system
ras_2080_proj <- list()
for(i in 1:length(names(ras_2080))){
  ras_2080_proj[[i]] <- projectRaster(ras_2080[[i]], crs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +no_defs")
}
ras_2080_proj <- stack(ras_2080_proj)

# Crop the rasters to the wanted extent
ras_2080_proj <- crop(ras_2080_proj, extent(-150, -103, 30, 65))

# Export the final rasters for downstream analyses
lapply(1:dim(ras_2080_proj)[3], function(x) writeRaster(ras_2080_proj[[x]], filename = as.character(paste(names(ras_2080_proj)[x], "2080_85.img", sep = "_")), format = "HFA", overwrite = T))

# Clean the R environment
rm(list = ls())

##############################
