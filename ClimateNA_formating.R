#####################################################
##                                                 ##
##  Redundancy Analysis (RDA): a Swiss-army knife  ##
##             for landscape genomics              ##
##                                                 ##
##              Capblancq & Forester               ##
##                                                 ##
##     Formatting ClimateNA rasters for WGS84      ##
##                                                 ##
#####################################################

# Begin by downloading the three data sets from ClimateNA.
# Here we use a climate normal for the period 1961-1990
# as well as two ensemble projections for 2050 and 2080.

# Go to this site: https://adaptwest.databasin.org/pages/adaptwest-climatena-cmip5/

# Scroll down to the following data sets and download them:
# *** In all cases, download the ASCII format ***
  # (1) Reference files, Elevation, ID
  # (2) Climate Normals, 1961-1990 period, 27 Bioclimatic variables
  # (3) Projection, Ensemble of 15 CMIP5 AOGCMs, RCP 8.5, 2050s, 27 Bioclimatic variables
  # (4) Projection, Ensemble of 15 CMIP5 AOGCMs, RCP 8.5, 2080s, 27 Bioclimatic variables

# Unzip using the default folder names & put these into a folder called "ClimateNA" located in the folder "Data"
# and you should be ready to run the below code (check your path if you run into trouble).

# After you've run the code below, all of the layers will be output to your "Data/ClimateNA" folder.
# You should have 162 total files for the three time periods, two files for each layer (XYZ.img and XYZ.img.aux.xml)
# You can now proceed to the RDA_Pinus_contorta.Rmd file.


# Required libraries
library(pegas)
library(raster)
library(rgdal)

# Set working directory
setwd("./Data/")


##############################
####   Period 1961-1990   ####

# Load 27 ClimateNA bioclimatic variables for 1961-1990:
ras_past <- stack(lapply(lapply(list.files("./NA_NORM_6190_Bioclim_ASCII/", pattern = ".asc$", full.names = T), read.asciigrid), raster))
names(ras_past) <- unlist(strsplit(unlist(lapply(strsplit(list.files("./NA_NORM_6190_Bioclim_ASCII", pattern = ".asc$", full.names = T), split = "./NA_NORM_6190_Bioclim_ASCII/"), function(x) x[2])), split = ".asc"))

# Set the projection based on the reference .asc file
crs(ras_past) <- crs(raster("./NA_Reference_files_ASCII/ClimateNA_ID.asc"))

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

# Load 27 ClimateNA bioclimatic variables predicted for 2050 using RCP8.5 climate scenario and the ensemble projection of 15 CMIP5 AOGCMs:
ras_2050 <- stack(lapply(lapply(list.files("./NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII/", pattern = ".asc$", full.names = T), read.asciigrid), raster))
names(ras_2050) <- unlist(strsplit(unlist(lapply(strsplit(list.files("./NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII", pattern = ".asc$", full.names = T), split = "./NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII/"), function(x) x[2])), split = ".asc"))

# Set the projection based on the reference .asc file
crs(ras_2050) <- crs(raster("./NA_Reference_files_ASCII/ClimateNA_ID.asc"))

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

# Load 27 ClimateNA bioclimatic variables predicted for 2080 using RCP8.5 climate scenario and the ensemble projection of 15 CMIP5 AOGCMs:
ras_2080 <- stack(lapply(lapply(list.files("./NA_ENSEMBLE_rcp85_2080s_Bioclim_ASCII/", pattern = ".asc$", full.names = T), read.asciigrid), raster))
names(ras_2080) <- unlist(strsplit(unlist(lapply(strsplit(list.files("./NA_ENSEMBLE_rcp85_2080s_Bioclim_ASCII", pattern = ".asc$", full.names = T), split = "./NA_ENSEMBLE_rcp85_2080s_Bioclim_ASCII/"), function(x) x[2])), split = ".asc"))

# Set the projection based on the reference .asc file
crs(ras_2080) <- crs(raster("./NA_Reference_files_ASCII/ClimateNA_ID.asc"))

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
