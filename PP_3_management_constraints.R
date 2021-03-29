####################################################################################################
#### Management constraints
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 02/17/2019
#### Last Modified: 07/09/2020
####################################################################################################
# REQUIRES GDAL installation (developed and tested with GDAl 3.0.2)
# GDAL is free and open source. Get it here: http://www.gisinternals.com/release.php
#
# This script is designed to create seven products related to treatment feasibility and cost. 
# The three fuel treatment options currently being considered are thin only, thin and 
# prescribed (Rx) fire, and Rx fire only. Each treatment needs a cost and feasibility raster
# (binary). For now, we're assuming any area thinned will be suitable for Rx fire, so the 
# thin and Rx fire feasibility is the same as the thin only. We also assume the cost of thin
# and Rx fire is the sum of the thin only and Rx fire costs. Also, since there is overlap in
# the feasibility rasters, we need another raster to represent combined feasibility so we 
# don't double count the same acres.
#
# Thin only
# 1) Cost as function of distance from roads and slope
# 2) Feasibility (restricted from wilderness and Upper Tier roadless and parts of Rocky Mountain
#    National Park)
#
# Rx fire only
# 3) Cost assumed to be uniform...for now, based on communication with local fuels and fire
#    planners (Bryan Karchut (ARP) and James White (CLRD)). Current costs are close to 
#    $1000/ac but they're hoping to get down to $500/ac by increasing the size of burns.
# 4) Feasibility (restricted from within 250 m of structures, wet forest types, and areas with
#    predicted > 30% crown fraction burned [CFB] under 70th percentile weather conditions)
#    NOTE: we are currently ignoring the CFB constraint because recent Rx fires have shown this
#    wasn't being adhered to.
# 
# Thin and Rx fire
# 5) Cost assumed to be the sum of thin only and Rx fire only costs
# 6) Feasibility assumed to be the same as thin only
# 
# Any treatment
# 7) Combined (total) feasibility (sum of thin only and Rx fire only feasibility)
#
####################################################################################################
#-> Get working and packages directory paths from command call
initial.options <- commandArgs(trailingOnly=F)
setwd(dirname(sub('--file=','',initial.options[grep('--file=',initial.options)])))
#setwd('C:/Users/bgannon/Desktop/tRADS/scripts')
wd <- getwd()
pd <- paste(c(unlist(strsplit(wd,'/'))[1:(length(unlist(strsplit(wd,'/')))-1)],
            'Portable_R/Packages'),collapse='/')
####################################################################################################

###########################################START MESSAGE############################################
cat('Management constraints\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages
.libPaths(pd)
packages <- c('raster','rgdal','rgeos','plyr')
for(package in packages){
	if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
		install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,lib.loc=pd,character.only=T))
	}
}

#-> Maximize raster processing speed
rasterOptions(maxmemory = 10^9)

#-> Load in spatial data
setwd(paste0(wd,'/INPUT/SPATIAL/RAW_LANDFIRE')) #===> Change directory
slope <- raster('us_slp2010')
cc <- raster('cc')
evt <- raster('evt')
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
nothin <- readOGR(dsn='VECTOR_INPUT.gdb',layer='exclude_thinning',verbose=F)
roads <- readOGR(dsn='VECTOR_INPUT.gdb',layer='Trans_RoadSegment',verbose=F)
WUI_structures <- readOGR(dsn='VECTOR_INPUT.gdb',layer='MS_Bing_structures',verbose=F)

#-> Load in tabular
setwd(paste0(wd,'/INPUT')) #===> Change directory
Rx_evt_feas <- read.csv('Rx_evt_feasibility.csv',header=T)

#-> Define projection
proj <- CRS('+proj=utm +zone=13 +datum=NAD83')

#-> Define Con function to match ArcGIS syntax
Con <- function(condition,trueValue,falseValue){
	return(condition*trueValue + (!condition)*falseValue)
}

#-> Rasterize vector in GDAL function
# inSPDF = spatial data frame you want to rasterize (if data in R)
# inPath = alternatively, path to shapefile on drive
# tr = template raster
# outPath = temporary drive folder to interface with GDAL
gdal_rasterize <- function(inSPDF=NULL,inPath=NULL,tr,outPath,returnRaster=T){
	vp <- Sys.getenv('PROJ_LIB') # Get original variable path
	Sys.setenv(PROJ_LIB = 'C:/Program Files (x86)/GDAL/projlib') # Temporarilly reset
	path_2_gdal_function <- 'C:/Program Files (x86)/GDAL/gdal_rasterize'
	if(!is.null(inSPDF)){
		inSPDF <- spTransform(inSPDF,crs(tr)) # Match spatial references
		suppressWarnings(writeOGR(inSPDF,dsn=outPath,layer='inVector',driver='ESRI Shapefile',
						 overwrite=T))
		inVector <- paste0(outPath,'/inVector.shp')				 
	}else{
		inVector <- inPath
	}
	outRaster <- paste0(outPath,'/rasterized.tif')
	#-> Set extent (-te <xmin> <ymin> <xmax> <ymax>)
	etext <- paste('-te',extent(tr)@xmin,extent(tr)@ymin,extent(tr)@xmax,extent(tr)@ymax)			
	#-> Set resolution (-tr <xres> <yres>)
	rtext <- paste('-tr',res(tr)[1],res(tr)[2])														
	args <- paste('-burn 1 -a_nodata 0',etext,rtext,inVector,outRaster)
	system2(path_2_gdal_function,args) # Call GDAL function (will return stdout)
	#-> Read in resulting raster if requested
	if(returnRaster==T){
		return(raster(paste0(outPath,'/rasterized.tif')))
	}
	Sys.setenv(PROJ_LIB = vp) # Reset original variable path
}

#-> Calculate Euclidean distance in GDAL function
# inRaster = raster representing the from extent (if in R)
# inPath = alternatively, path to raster on drive
# outPath = temporary drive folder to interface with GDAL
gdal_proximity <- function(inRaster=NULL,inPath=NULL,outPath,returnRaster=T){
	path_2_gdal_function <- '"C:/Program Files (x86)/GDAL/gdal_proximity.py"' # Quotes b/c space
	vp <- Sys.getenv('PROJ_LIB') # Get original variable path
	Sys.setenv(PROJ_LIB = 'C:/Program Files (x86)/GDAL/projlib') # Temporarilly reset
	if(!is.null(inRaster)){
		writeRaster(inRaster,paste0(outPath,'/inRaster.tif'),format='GTiff',overwrite=T)
		inRaster <- paste0(outPath,'/inRaster.tif')				 
	}else{
		inRaster <- inPath
	}
	outRaster <- paste0(outPath,'/distance_raster.tif')
	args <- paste(path_2_gdal_function,'-distunits GEO',inRaster,outRaster)
	system2('Python',args)
	#-> Read in resulting raster if requested
	if(returnRaster==T){
		return(raster(paste0(outPath,'/distance_raster.tif')))
	}
	Sys.setenv(PROJ_LIB = vp) # Reset original variable path
}

#############################################END SET UP#############################################

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/INPUT/SPATIAL/Constraints'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Empty output folder
	d_list <- list.files()
	if(length(d_list) > 0){
		unlink(d_list,recursive=T)
	}
}

##########################################END CLEAR OUTPUT##########################################

###########################################START ANALYSIS###########################################

#-> Crop data 
nothin <- crop(nothin,cc)
roads <- crop(roads,cc)
WUI_structures <- crop(WUI_structures,cc)

###---> Thining treatment costs

#-> Model variables 
MaxCost = 10000 # $/ac                                   
BaseCost = 2500 # $/ac
slopeTH = 40 # percent
slopeMAX = 200 # percent  
distTH = 800 # meters            
distMAX = 6400 # meters

#-> Create TEMP processing folder
setwd(paste0(wd,'/INTERMEDIATE')) #===> Change directory
dir.create('TEMP')
outPath <- paste0(wd,'/INTERMEDIATE/TEMP')

#-> Rasterize roads
gdal_rasterize(inSPDF=roads,inPath=NULL,tr=cc,outPath=outPath,returnRaster=F)

#-> Distance from roads
rdist <- gdal_proximity(inRaster=NULL,inPath=paste0(outPath,'/rasterized.tif'),outPath=outPath,
                        returnRaster=T)

#-> Calculate road-distance additional costs
m1 <- (MaxCost-BaseCost)/(distMAX-distTH)
rcost <- Con(rdist <= distTH,0,(rdist-distTH)*m1)

#-> Convert slope to percent (in degrees)
PerSlope <- tan(slope*(pi/180))*100
m2 <- (MaxCost-BaseCost)/(slopeMAX-slopeTH)
scost <- Con(PerSlope <= slopeTH,0,(PerSlope-slopeTH)*m2)

#-> Combine base, road, and slope costs
mocost <- rcost + scost + BaseCost
mocost <- Con(mocost > MaxCost,MaxCost,mocost)

#-> Save
setwd(paste0(wd,'/INPUT/SPATIAL/Constraints')) #===> Change directory
crs(mocost) <- proj
writeRaster(mocost,'mocost.tif',format='GTiff',overwrite=T)

###---> Thin and Rx fire treatment costs

setwd(paste0(wd,'/INPUT/SPATIAL/Constraints')) #===> Change directory
writeRaster((mocost+1000),'mRxcost.tif',format='GTiff',overwrite=T)

###---> Rx fire treatment costs

setwd(paste0(wd,'/INPUT/SPATIAL/Constraints')) #===> Change directory
Rxcost <- cc*0 + 1000
writeRaster(Rxcost,'Rxcost.tif',format='GTiff',overwrite=T)

###---> Thin only treatment feasibility

#-> Rasterize
mofeas <- rasterize(nothin,cc,field=rep(0,nrow(nothin)))
mofeas[is.na(mofeas)] <- 1

#-> Filter to forest
forest <- Con(cc >= 10,1,0)
mofeas <- mofeas*forest

#-> Save
setwd(paste0(wd,'/INPUT/SPATIAL/Constraints')) #===> Change directory
crs(mofeas) <- proj																	
writeRaster(mofeas,'mofeas.tif',format='GTiff',overwrite=T)

###---> Thin and Rx fire feasibility

writeRaster(mofeas,'mRxfeas.tif',format='GTiff',overwrite=T)

###---> Rx fire only feasibility
# Note that previous versions used a fire effects filter based on crown fraction burned.
# This was removed based on observations that the USFS was burning in many areas we mapped
# as infeasible. Similar criteria may be appropirate to include in other projects.

#-> Rasterize structures
gdal_rasterize(inSPDF=WUI_structures,inPath=NULL,tr=cc,outPath=outPath,returnRaster=F)

#-> Distance from structures
sdist <- gdal_proximity(inRaster=NULL,inPath=paste0(outPath,'/rasterized.tif'),outPath=outPath,
                        returnRaster=T)

#-> Create structure distance filter
sdfilt <- Con(sdist <= 250,0,1)

#-> Restrict to ecologically-appropriate forest types
rcl <- Rx_evt_feas[,c('VALUE','FEASIBLE')]
eafilt <- reclassify(evt,rcl)

#-> Combine
Rxofeas <- sdfilt*eafilt

#-> Save
setwd(paste0(wd,'/INPUT/SPATIAL/Constraints')) #===> Change directory
crs(Rxofeas) <- proj	
writeRaster(Rxofeas,'Rxfeas.tif',format='GTiff',overwrite=T)

###---> Combined feasibility

#-> Add together thin and Rx only rasters
cfeas <- mofeas + Rxofeas
cfeas <- Con(cfeas > 0,1,0)

#-> Save
setwd(paste0(wd,'/INPUT/SPATIAL/Constraints')) #===> Change directory
crs(cfeas) <- proj	
writeRaster(cfeas,'cfeas.tif',format='GTiff',overwrite=T)

############################################END ANALYSIS############################################

######################################START CLEAR INTERMEDIATE######################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/INTERMEDIATE'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	unlink('TEMP',recursive=T)
}

#######################################END CLEAR INTERMEDIATE#######################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
#############################################END MESSAGE############################################
