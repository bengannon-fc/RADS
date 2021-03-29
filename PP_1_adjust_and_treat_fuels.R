####################################################################################################
#### Adjust and treat fuels
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 07/27/2019
#### Last Modified: 08/20/2020
####################################################################################################
# Summary: this script adjusts the LANDFIRE 2016 data to represent current (baseline) conditions, 
# and hypothetical fuel treatments applied to current conditions. The baseline conditions are 
# adjusted in lodgepole pine to intensify fire behavior. The output includes individual raster 
# files for each layer and a multi-band GeoTIFF fuelscape file.
#
# Data sources:
# CC, CH, CBH, CBD, and FBFM40 from LANDFIRE 2016 refresh
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
cat('Adjust and treat fuels\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages
.libPaths(pd)
packages <- c('raster','rgdal','rgeos')
for(package in packages){
	if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
		install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,lib.loc=pd,character.only=T))
	}
}

#-> Maximize raster processing speed
rasterOptions(maxmemory = 10^9)

#-> Load in LANDFIRE data
setwd(paste0(wd,'/INPUT/Spatial/RAW_LANDFIRE')) #===> Change directory
cc <- raster('cc')
ch <- raster('ch')
cbh <- raster('cbh')
cbd <- raster('cbd')
fbfm <- raster('fbfm40')
dem <- raster('us_dem2010')
asp <- raster('us_asp2010')
slp <- raster('us_slp2010')
evt <- raster('evt')

#-> Load in fuel treatment data data
setwd(paste0(wd,'/INPUT/Spatial/Fuel_treatments/Compiled')) #===> Change directory
CanTrts <- shapefile('CanEff_flattened.shp')
SurfTrts <- shapefile('SurfEff_flattened.shp')

#-> Relabel projection information so it plays well with Arc
proj <- '+proj=utm +zone=13 +ellps=GRS80 +datum=NAD83'
crs(cc) <- CRS(proj)
crs(ch) <- CRS(proj)
crs(cbh) <- CRS(proj)
crs(cbd) <- CRS(proj)
crs(fbfm) <- CRS(proj)
crs(dem) <- CRS(proj)
crs(asp) <- CRS(proj)
crs(slp) <- CRS(proj)
crs(evt) <- CRS(proj)
crs(CanTrts) <- CRS(proj)
crs(SurfTrts) <- CRS(proj)

#-> Load in treatment tables
setwd(paste0(wd,'/INPUT')) #===> Change directory
canopy_effs <- read.csv('canopy_effects.csv',header=T)
fbfm_rcl <- read.csv('surface_effects.csv',header=T)

#-> Define Con function to match ArcGIS syntax
Con <- function(condition,trueValue,falseValue){
	return(condition*trueValue + (!condition)*falseValue)
}

#############################################END SET UP#############################################	

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios'))) #===> Change directory

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

#-> Add GIS subdirectory
dir.create('Baseline')
dir.create('Thin')
dir.create('RxFire')
dir.create('Complete')

##########################################END CLEAR OUTPUT##########################################

###########################################START ANALYSIS###########################################

###---> Baseline (for prospective analysis)
setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios/Baseline')) #===> Change directory

#-> Adjust lodgepole cbh
m_cbh <- Con(evt==7050,cbh*0.7,cbh) # Lowered canopy base height 30%

#-> Adjust lodgepole fbfm
m_fbfm <- Con(evt==7050 & (fbfm==181 | fbfm==183),185,fbfm) # Changed timber litter fuel models

#-> Merge canopy effects to spatial data
CanTrts <- merge(CanTrts,by.x='CanEff',canopy_effs,by.y='Treatment',all.x=T)

#-> Canopy bulk density
AF <- rasterize(CanTrts,cc,field='cbd_AF',background=1)
cbd <- cbd*AF
writeRaster(cbd,'cbd.tif',format='GTiff',overwrite=T)

#-> Canopy base height
AF <- rasterize(CanTrts,cc,field='cbh_AF',background=1)
m_cbh <- m_cbh*AF
writeRaster(m_cbh,'cbh.tif',format='GTiff',overwrite=T)

#-> Canopy cover
AF <- rasterize(CanTrts,cc,field='cc_AF',background=1)
cc <- cc*AF
writeRaster(cc,'cc.tif',format='GTiff',overwrite=T)

#-> Canopy height
AF <- rasterize(CanTrts,cc,field='ch_AF',background=1)
ch <- ch*AF
writeRaster(ch,'ch.tif',format='GTiff',overwrite=T)

#-> Update surface fuels
strts <- unique(SurfTrts$SurfEff)
for(i in 1:length(strts)){
	trast <- rasterize(SurfTrts[SurfTrts$SurfEff==strts[i],],m_fbfm)
	rcl <- fbfm_rcl[,c('FBFM40',strts[i])]
	tfbfm <- reclassify(m_fbfm,rcl)
	m_fbfm <- Con(is.na(trast),m_fbfm,tfbfm)
}
writeRaster(m_fbfm,'fbfm.tif',format='GTiff',overwrite=T)

#-> Save as multiband GeoTIFF
setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios')) #===> Change directory
writeRaster(brick(dem,slp,asp,m_fbfm,cc,ch,m_cbh,cbd),'baseline.tif',format='GTiff',
            overwrite=T)

###---> Simulate fuel treatments

#-> Define treat fire behavior fuel model funtion
# Treat_rcl = fbfm reclassification table
# tfield = treatment name in reclassification table
# tname = appended to output raster
treatSurface <- function(fbfm,Treat_rcl,tfield,flist){
	rcl <- Treat_rcl[,c('FBFM40',tfield)]
	flist[[4]] <- reclassify(fbfm,rcl)
	writeRaster(flist[[4]],paste0('fbfm.tif'),format='GTiff',overwrite=T)
	return(flist)
}

#-> Define treat canopy fuels function 
# CFAFs = vector of cc_AF, ch_AF, cbh_AF, cbd_AF
# baseline = list of cc, ch, cbh, cbd rasters
# tname = appended to output raster
treatCanopy <- function(baseline,CFAFs,flist){
	isForest <- Con(baseline[[1]]*CFAFs[1] >= 10,1,0) # Check if post-trt is non-forest
	bnames <- c('cc','ch','cbh','cbd')
	for(i in 1:length(baseline)){
		flist[[4+i]] <- Con(isForest==1,baseline[[i]]*CFAFs[i],0)
		writeRaster(flist[[4+i]],paste0(bnames[i],'.tif'),format='GTiff',overwrite=T)
	}
	return(flist)
}	

#-> Create baseline list of canopy attributes
baseline <- list(cc,ch,m_cbh,cbd)

#-> Base list
blist <- list(dem,slp,asp)

#-> Canopy adjustment factor order
corder <- c('cc_AF','ch_AF','cbh_AF','cbd_AF')

#-> Thin treatment
# Fule et al. 2012, Ziegler et al. 2017
setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios/Thin')) #===> Change directory
flist <- treatSurface(m_fbfm,Treat_rcl=fbfm_rcl,tfield='Manage',flist=blist)
CFAFs <- as.numeric(canopy_effs[canopy_effs$Treatment=='Thin',corder]) 
flist <- treatCanopy(baseline,CFAFs,flist)
setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios')) #===> Change directory
writeRaster(brick(flist),'thin.tif',format='GTiff',overwrite=T)

#-> Rx fire treatment
# Stephens and Moghaddas 2005
setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios/RxFire')) #===> Change directory
flist <- treatSurface(m_fbfm,Treat_rcl=fbfm_rcl,tfield='RxFire',flist=blist)
CFAFs <- as.numeric(canopy_effs[canopy_effs$Treatment=='RxFire',corder])
flist <- treatCanopy(baseline,CFAFs,flist)
setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios')) #===> Change directory
writeRaster(brick(flist),'RxFire.tif',format='GTiff',overwrite=T)

#-> "Complete" treatment (mechanical followed by Rx fire)
# Fule et al. 2012, Ziegler et al. 2017
setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios/Complete')) #===> Change directory
flist <- treatSurface(m_fbfm,Treat_rcl=fbfm_rcl,tfield='RxFire',flist=blist)
CFAFs <- as.numeric(canopy_effs[canopy_effs$Treatment=='Complete',corder])
flist <- treatCanopy(baseline,CFAFs,flist)
setwd(paste0(wd,'/INPUT/Spatial/Fuel_scenarios')) #===> Change directory
writeRaster(brick(flist),'complete.tif',format='GTiff',overwrite=T)

############################################END ANALYSIS############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
