
# TODO:
# - what value of field capacity to use? 316 hPa, 63 hPa or 330 hPa (Wikipedia)?
# - nfk should rather be fk - theta_pwp insteadt of fk - theta_r (Wikipedia)?
# - in doc of pft.rawls suction: problem, maybe of units? All calculated values larger than ranges in Maidment Table 5.5.5)
#   - in pft.rawls: ret_val[, "S_f"] = ret_val[, "S_f"] * 10 -> values in mm? But nevertheless values are still high
# - ksat sometimes zero (very low clay and/or silt contents)?!
# - fk sometimes larger than theta_s?!
# - theta_r sometimes larger than pwp
# - water content values sometimes unrealistic (compared to ranges given in Maidment)

#install.packages("soilwaterfun", repos="http://R-Forge.R-project.org")
library(soilwaterfun)
#install.packages("soilwaterptf", repos="http://R-Forge.R-project.org")
library(soilwaterptf)
# or subversion (Linux commands):
# svn checkout svn://scm.r-forge.r-project.org/svnroot/soilwater/
# if your R version is too old and you don't want to update R:
# get an older package version with 'svn update PATH -r REV' PATH=path to download dir, REV=revision number you want
# install that version: 'R CMD INSTALL FILE' FILE=tar.gz file with package source data

#install.packages("soiltexture")
library(soiltexture)


### INPUT ###

setwd("/home/tobias/Promotion/Modellierung/Bengue/Preprocessing/soil_parameters/HWSD/calc_pars/")


# output directory
dir_out <- "out"

# particle size [mm] classes definition of your soil data source
part_def <- data.frame(class_id=c(1,2,3),
                       description=c("clay", "silt", "sand"),
                       upper_limit=c(0.002, 0.05, 2))



library(spgrass6)
# start GRASS session
initGRASS( gisBase="/opt/grass",
           home=tempdir(),
           location="Bengue",
           mapset="PERMANENT",
           gisDbase="/home/tobias/Promotion/Modellierung/Bengue/grassdata",
           override=T)

# remove existing mask from location
execGRASS("r.mask", flags=c("r"))
# use standard region setting
execGRASS("g.region", flags=c("d"))


#get non-horizon-specific data (soil IDs + further data as available)
layer2get=c("HWSD_soils", "HWSD_REF_DEPTH")
rasts=readRAST6(vname=layer2get, return_SGDF=FALSE)

# soil IDs
soils <- rasts$dataList[[1]]

# soil descriptions (raster categories of soil types)
soil_desc <- attr(execGRASS("r.category", parameters=list(map=layer2get[1],fs=":")),"resOut")
soil_desc <- grep("^[0]:", soil_desc, invert=T, value=T)
soil_desc <- as.data.frame(matrix(unlist(strsplit(soil_desc, ":")), ncol=2, byrow=T))
names(soil_desc) <- c("ID", "descr")
soil_desc$ID <- as.integer(as.character(soil_desc$ID))

# further data
total_depth <- rasts$dataList[[2]]

# agregate non-horizon specific properties to soil types
soil_sum <- aggregate(x=total_depth*10, by=list(ID=soils), FUN=mean, na.rm=TRUE) #aggregate according to soil_id
soil_sum <- merge(soil_desc,soil_sum, by="ID")
names(soil_sum) <- c("ID", "descr", "depth")

rm(total_depth, rasts, soil_desc)
gc()
gc()



### CALCULATIONS ###

# prepare horizon specific output object
hor_par_out <- matrix(data=NA, nrow=1, ncol=18,
                       dimnames=list(NULL, c("pid","description","soil_id","position","theta_r",
                                             "theta_pwp","fk","fk63","nfk","theta_s","thickness",
                                             "ks","suction","pore_size_i","bubb_pres","coarse_frag",
                                             "shrinks", "soil_dens")))

# prepare soil profile specific output object
soil_par_out <- matrix(data=NA, nrow=nrow(soil_sum), ncol=13,
                       dimnames=list(NULL, c("pid", "description", "bedrock_flag", "alluvial_flag", "b_om",
                                             "a_musle_k", "a_clay", "a_silt", "a_sand", "a_f_csand", "a_f_cl_si",
                                             "a_f_orgc", "a_f_hisand")))

soil_par_out[,"pid"] <- soil_sum[, "ID"]
soil_par_out[,"description"] <- paste0(soil_sum[, "descr"])
soil_par_out[,"bedrock_flag"] <- 1
soil_par_out[,"alluvial_flag"] <- ifelse(grepl("fluvisol", paste0(soil_sum[, "descr"]), ignore.case=T), 1, 0)

# prepare particle distribution output object
part_distr_out <- matrix(data=NA, nrow=nrow(soil_sum)*3, ncol=3,
                         dimnames=list(NULL,c("soil_id","class_id","fraction")))

prefix="HWSD_"
layer2get=c("CLAY","SILT","BULK_DENSITY","OC","GRAVEL","SAND")

i <- 0

# number of layers for this soil
nlay <- 2

# loop over layers
for (l in 1:nlay) {
  
  ll <- c("T_","S_")[l]
  
  # load rasters from GRASS location
  rasts=readRAST6(vname=paste0(prefix,ll,layer2get), return_SGDF=FALSE)
  
  #reshape to dataframe as required by ptf.wosten()
  soil_attributes=data.frame(clay   =rasts$dataList[[1]],
                             silt   =rasts$dataList[[2]],
                             bulkD  =rasts$dataList[[3]],
                             om     =1.74 * rasts$dataList[[4]], #convert OC to OM [%]
                             coarse_frag =rasts$dataList[[5]],
                             topSoil=as.numeric(ll=="T_"),
                             sand   =rasts$dataList[[6]])                             
  
  rm(rasts)
  gc()
  gc()
  
  # remove NA values
  is_na=
    soil_attributes$silt   ==0 | !is.finite(soil_attributes$silt) |
    soil_attributes$clay   ==0 | !is.finite(soil_attributes$clay) |
    soil_attributes$om     ==0 | !is.finite(soil_attributes$om)   |
    soil_attributes$bulkD  ==0 | !is.finite(soil_attributes$bulkD |
                                              soil_attributes$sand   ==0 | !is.finite(soil_attributes$sand)) 
  
  soil_dat <- soil_attributes[!is_na,]
  
  #aggregate properties from basic horizon input data
  soil_sum2 = aggregate(x=soil_dat, by=list(ID=soils[!is_na]), FUN=mean, na.rm=TRUE) #aggregate according to soil_id
  
  rm(soil_dat, soil_attributes)
  gc()
  gc()
  
  # merge with other soil information
  soil_sum2 <- merge(soil_sum2, soil_sum[grep("depth", names(soil_sum), invert=T)], by="ID")
  
  # organic mater content (soil profile specific; value only for top layer)
  if(l == 1)
    soil_par_out[,"b_om"] <- merge(soil_par_out, soil_sum2[,c("ID", "om")], by.x="pid", by.y="ID")$om/100
  
  # collection information about particle distribution (topsoil only)
  if(l == 1) {
    part_distr_out <- data.frame(soil_id=soil_sum2$ID, class_id=rep(1:3, each=nrow(soil_sum2)),
                 fraction=c(soil_sum2$clay/100, soil_sum2$silt/100, soil_sum2$sand/100))
    if(any(round(tapply(part_distr_out$fraction, part_distr_out$soil_id, sum), 3) != 1))
      stop("Fractions of particle classes do not sum up to one for every soil type! Check your input!")
  }

  # depth (thickness) of horizon [mm]
  if (l == 1) {
    depth <- 300
  } else {
    depth <- soil_sum[,"depth"] - 300
  }

  
  # residual volumetric water content - theta_r
  theta_r <- pft.rawls(soilprop = soil_sum2[,c("clay", "bulkD", "silt", "om")], parameters="theta_r", h=0)[,"theta_r"]
  
  # vol. water content at permanent wilting point (1500 kPa/15000 cm) - theta_pwp
  theta_pwp <- pft.rawls(soilprop = soil_sum2[,c("clay", "bulkD", "silt", "om")], parameters="theta", h=15000)[,"theta"]
  
  # vol. water content at field capacity (316 hPa / pF=2.6) - fk
  fk <- pft.rawls(soilprop = soil_sum2[,c("clay", "bulkD", "silt", "om")], parameters="theta", h=316)[,"theta"]
  
  # vol. water content at field capacity (63 hPa / pF=1.8) - fk63
  fk63 <- pft.rawls(soilprop = soil_sum2[,c("clay", "bulkD", "silt", "om")], parameters="theta", h=63)[,"theta"]
  
  # usable field capacity - nfk
  nfk <- fk - theta_r
  
  # Saturated volumetric soil water content - theta_s
  theta_s <- ptf.wosten.theta.s(soil_sum2$clay, soil_sum2$bulkD, soil_sum2$silt,
                                soil_sum2$om, ifelse(l==1,1,0))

  # make sure that values are plausible (use values from wikipedia, https://en.wikipedia.org/wiki/Water_content, and Maidment: Handbook of Hydrology)
  if(any(theta_r > theta_pwp))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[theta_r > theta_pwp],collapse=", "), ", layer ", l, ": residual water content is greater than permanent wilting point!"))
  
  if(any(theta_pwp > fk))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[theta_pwp > fk],collapse=", "), ", layer ", l, ": permanent wilting point is greater than field capacity!")) 
  
  if(any(fk > theta_s))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[fk > theta_s],collapse=", "), ", layer ", l, ": field capacity is greater than saturated water content!"))
  
  if(any(theta_s > 0.5) | any(theta_s < 0.2))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[theta_s > 0.5 | theta_s < 0.2], collapse=", "), ", layer ", l, ": value of saturated water content is ", paste(theta_s[theta_s > 0.5 | theta_s < 0.2],collapse=", "), " wich seems rather unlikely!"))
    
  if(any(fk > 0.35) | any(fk < 0.1))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[fk > 0.35 | fk < 0.1], collapse=", "), ", layer ", l, ": value of field capacity is ", paste(fk[fk > 0.35 | fk < 0.1],collapse=", "), " wich seems rather unlikely!"))
  
  if(any(theta_pwp > 0.25) | any(theta_pwp < 0.01))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[theta_pwp > 0.25 | theta_pwp < 0.01], collapse=", "), ", layer ", l, ": value of permanent wilting point is ", paste(theta_pwp[theta_pwp > 0.25 | theta_pwp < 0.01],collapse=", "), " wich seems rather unlikely!"))
  
  if(any(theta_r > 0.1) | any(theta_r < 0.001))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[theta_r > 0.1 | theta_r < 0.001], collapse=", "), ", layer ", l, ": value of residual water content is ", paste(theta_r[theta_r > 0.1 | theta_r < 0.001],collapse=", "), " wich seems rather unlikely!"))
  
  
  # saturated hydraulic conductivity - ks [mm/day]
    ks <- ptf.wosten.ksat(soil_sum2$clay, soil_sum2$bulkD, soil_sum2$silt,
                          soil_sum2$om, ifelse(l==1,1,0), c("mm", "day"))

  if(any(ks > 6000 | ks < 10))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[ks > 6000 | ks < 10],collapse=", "), ", layer ", l, ": value of saturated hydraulic conductivity is ", paste(ks[ks > 6000 | ks < 10],collapse=", "), " mm/day wich seems rather unlikely!"))
  
  
  # suction suction at the wetting front - suction [cm]
  suction <- pft.rawls(soilprop = soil_sum2[,c("clay", "bulkD", "silt", "om")], parameters="S_f", h=0)[,"S_f"]
  
  if(any(suction > 150 | suction < 1))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[suction > 150 | suction < 1],collapse=", "), ", layer ", l, ": value of suction at wetting front is ", paste(suction[suction > 150 | suction < 1],collapse=", "), " cm wich seems rather unlikely!"))
  
  
  # pore-size-index - pore_size_i
  lambda <- pft.rawls(soilprop = soil_sum2[,c("clay", "bulkD", "silt", "om")], parameters="lambda", h=0)[,"lambda"]
  
  if(any(lambda > 1 | lambda < 0.04))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[lambda > 1 | lambda < 0.04],collapse=", "), ", layer ", l, ": value of pore-size-index is ", paste(lambda[lambda > 1 | lambda < 0.04],collapse=", "), " wich seems rather unlikely!"))
  
  
  # bubbling pressure - bubb_pres [cm]
  bubb_pres <- pft.rawls(soilprop = soil_sum2[,c("clay", "bulkD", "silt", "om")], parameters="h_b", h=0)[,"h_b"]
  
  if(any(bubb_pres > 190 | bubb_pres < 1))
    warning(paste0("For soil no. ", paste(soil_sum2$ID[bubb_pres > 190 | bubb_pres < 1],collapse=", "), ", layer ", l, ": value of bubbling pressure is ", paste(bubb_pres[bubb_pres > 190 | bubb_pres < 1],collapse=", "), " cm wich seems rather unlikely!"))
  
  
  # write into output object
  out_t <- data.frame(pid=NA, description=soil_sum2[,"descr"], soil_id=soil_sum2[, "ID"], position=l,
            theta_r, theta_pwp, fk, fk63, nfk, theta_s, thickness=depth, ks, suction, pore_size_i=lambda, 
            bubb_pres, coarse_frag=soil_sum2[,"coarse_frag"]/100, shrinks=0, soil_dens=soil_sum2$bulkD*1000)
  hor_par_out <- rbind(hor_par_out, out_t)
}


hor_par_out <- hor_par_out[-1,]
hor_par_out$pid <- 1:nrow(hor_par_out)

# write output
write.table(hor_par_out, file=paste0(dir_out, "/horizons.dat"), quote=F, sep="\t",
            row.names=F)

write.table(soil_par_out, file=paste0(dir_out, "/soil.dat"), quote=F, sep="\t",
            row.names=F)

write.table(part_distr_out, file=paste0(dir_out, "/r_soil_contains_particles.dat"), quote=F, sep="\t",
            row.names=F)

write.table(part_def, file=paste0(dir_out, "/particle_classes.dat"), quote=F, sep="\t",
            row.names=F)
