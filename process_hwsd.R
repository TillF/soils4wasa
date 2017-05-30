# process HWSD raster file
# download from http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HTML/HWSD_Data.html?sb=4
# data is cropped to the target region
# specified soil property raster files and raster with soil IDs and FAO90 names are loaded into the respective location in GRASS
#
# grass location needs to exist already
# 
# HWSD database 2 options:
# 1) HWSD database in SQLite format (converted from MS Access)
# 2) database as csv table with all relevant information already extracted from database
#
# ATTENTION: depending on the size of your target area the script might need lots of computer ressources (time + memory)
# If an error occurs during executing of GRASS commands, manual garbage collection
# by applying 'gc()' multiple times in a row sometimes helps.
#
# ERROR in HWSD.mdb
# in HWSD.mdb I downloaded there were errors in the relations of the database:
# For tables D_Symbol_* the column CODE was related to SU_SYM* in HWSD_DATA but it
# should be related to SU_CODE*. This has to be fixed, otherwise errors in queries
# will occur.
#
# Information on BULK DENSITY
# 2 methods:
# REFERENCE BULK DENSITY ('*__REF_BULK_DENSITY'):
# equations developed by Saxton et al. (1986) that relate to the texture of the
# soil only -> overestimation especially for soils with high OM content or high porosity.
# SOTWIS BULK DENSITY ('*_BULK_DENSITY'):
# based on information in SOTWIS database of soil texture, organic matter content
# and porosity
#
# Information on Reference soil depth:
# Reference depth of the soil unit. Reference soil depth of all soil units are set 
# at 100 cm, except for Rendzinas and Rankers of FAO-74 and Leptosols of FAO-90,
# where the reference soil depth is set at 30 cm, and for Lithosols of FAO-74 and
# Lithic Leptosols of FAO-90, where it is set at 10 cm6. An approximation of actual
# soil depth can be derived through accounting for relevant depth limiting soil
# phases, obstacles to roots and occurrence of impermeable layers (the latter two
# refer to ESDB only).


library(raster)
library(spgrass6)
library(rgdal)
#library(RSQLite) # not needed when you use a csv table
library(maptools)

### SETTINGS ###
setwd("/home/tobias/Promotion/Modellierung/Bengue/Preprocessing/soil_parameters/HWSD/extract_data/")

# start GRASS session
initGRASS( gisBase="/opt/grass",
           home=tempdir(),
           location="Bengue",
           mapset="PERMANENT",
           gisDbase="/home/tobias/Promotion/Modellierung/Bengue/grassdata",
           override=T)

# overwrite when writing output into GRASS location?
ov_flag=T

# raster data file (the .hwd and .hdr files have to be in the directory as well!)
hwsd_file <- "/home/tobias/Promotion/Modellierung/Jaguaribe/Data/HWSD/data_raw/hwsd.bil"

# projection string of hwsd file (shouldn't need to be changed)
hwsd_proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# name of the HWSD database (has to be in SQLite format (*.db) OR a csv table (*.csv))
db_file <- "/home/tobias/Promotion/Modellierung/Jaguaribe/Data/HWSD/HWSD_DATA.csv"

# specify soil properties that should be extracted from the database
# a raster file for each soil property is created and written into grass as "HWSD_[name]"
# additionally "MU_GLOBAL", "SHARE" and the FAO90 naming "VALUE_90" and "SYMBOL_90" are needed
soil_tables <- c("REF_DEPTH",
                 "T_GRAVEL", "T_SILT", "T_CLAY", "T_BULK_DENSITY", "T_OC", "T_SAND",
                 "S_GRAVEL", "S_SILT", "S_CLAY", "S_BULK_DENSITY", "S_OC", "S_SAND")

# mask to define GRASS region
mask_reg <- "MASK_corr"


### CALCULATIONS ###
# remove existing mask from location
execGRASS("r.mask", flags=c("r"))
# use standard region setting
execGRASS("g.region", flags=c("d"))

# read HWSD raster data
hwsd <- raster(hwsd_file)

# assign projection (by default longlat with datum WGS84)
proj4string(hwsd) <- hwsd_proj

# grass region extent in longlat coordinates
loc_info <- attr(execGRASS("g.region", flags=c("b", "g")), "resOut")
loc_info <- as.numeric(gsub("[a-z_=]*", "", loc_info))

# crop hwsd to extent of grass location (slightly larger in case of small deviations from utm projection)
hwsd_crop <- crop(hwsd, extent(loc_info[3]-.2, loc_info[4]+.2, loc_info[2]-.2, loc_info[1]+.2))

# re-project to projection of grass location
proj_raster <- projectExtent(hwsd_crop, CRS(getLocationProj()))
res(proj_raster) <- c(gmeta6()$nsres, gmeta6()$ewres)
gc() # manual garbage collection
gc() # manual garbage collection
hwsd_crop_proj <- projectRaster(hwsd_crop, proj_raster, method="ngb")
gc() # manual garbage collection
gc() # manual garbage collection

# set values=zero (oceans) to NA
values(hwsd_crop_proj)[which(getValues(hwsd_crop_proj) == 0)] <- NA

# the following workaround is cumbersome but I had problems directly converting raster to SGDF
# save as raster file
writeRaster(hwsd_crop_proj, "temp/soils", format="EHdr", overwrite=ov_flag)
# load raster file
hwsd_crop_proj_sgdf <- readGDAL("temp/soils")
# convert to integer
hwsd_crop_proj_sgdf@data[[1]] <- as.integer(hwsd_crop_proj_sgdf@data[[1]])
# load map with soil IDs into grass
gc() # manual garbage collection
gc() # manual garbage collection
gc() # manual garbage collection
gc() # manual garbage collection
writeRAST6(hwsd_crop_proj_sgdf, "HWSD_soils", overwrite=ov_flag, flags=c("o"))

# specify as mask
#execGRASS("r.mask", parameters=list(input="HWSD_soils"))

# access database
if (grepl("*.db", db_file)) {
  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = db_file)
} else if (grepl("*.csv", db_file)) {
  dat_all <- read.table(db_file, header=T, sep=",")
} else {
  stop("Couldn't detect database file. Is the format + file name ending correct (*.db or *.csv)?")
}

# FAO90 names as category labels in raster file
# loop over soil types
write(x=NULL,file="temp/soil_cat",append=F)
for (s in unique(hwsd_crop_proj)) {
  # main soil type (first occurence of certain soil type (MU_GLOBAL))
  cat_string <- paste(as.vector(as.matrix(dat_all[which(dat_all$MU_GLOBAL == s),c("SYMBOL_90","VALUE_90")][1,])), collapse=" - ")
  cat_string <- paste(s,cat_string,sep=":")
  
  # create reclassification file
  write(x=cat_string,file="temp/soil_cat",append=T)
}

# classify raster file
execGRASS("r.category", parameters=list(map="HWSD_soils", rules="temp/soil_cat"))


# create necessary tables
for (t in soil_tables) {
  rast <- hwsd_crop_proj
  rast_vals <- values(rast)
  
  # loop over soil types
  for (s in unique(rast)) {
  
    # get data
    if (grepl("*.db", db_file)) {
      dat <- dbGetQuery(con, paste0("select * from HWSD_DATA where MU_GLOBAL = ", s))
    } else if (grepl("*.csv", db_file)) {
      dat <- dat_all[which(dat_all[,"MU_GLOBAL"] == s),]
    }
    
    # calculate area-weighted mean of soil property
    dat <- weighted.mean(dat[,t], dat[,"SHARE"], na.rm=T)
    
    # replace soil type values in raster file by area-weighted mean
    rast_vals[which(rast_vals == s)] <- dat
  }
  
  # replace values in raster file
  values(rast) <- rast_vals
  
  # the following workaround is cumbersome but I had problems directly converting raster to SGDF
  # save as raster file
  writeRaster(rast, paste0("temp/",t), format="EHdr", overwrite=ov_flag)
  
  # load raster file
  rast <- readGDAL(paste0("temp/",t))
  
  # write into grass location
  writeRAST6(rast, paste0("HWSD_", t), overwrite=ov_flag, flags=c("o"))
}

closeAllConnections()
