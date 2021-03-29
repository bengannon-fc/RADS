####################################################################################################
#### Process fuel treatment data
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 07/02/2020
#### Last Modified: 07/02/2020
####################################################################################################
# Summary: this script assembles historical fuel treatment data to update LANDFIRE 2016 to current 
# conditions.
#
# Data sources:
# USFS hazardous fuel treatment polygons (accessed 03/17/2020)
# DOI NFPORS hazardous fuel treatment polygons (current to end of 2019) - none in planning extent!
# CSFS fuel treatment polygons (current to end of 2017)
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
cat('Process fuel treatment data\n',sep='')
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
rasterOptions(maxmemory=10^9)

#-> Load in spatial data
setwd(paste0(wd,'/INPUT/Spatial/Fuel_treatments/RAW_agency_data')) #===> Change directory
USFS <- shapefile('USFS_Activity_HazFuelTrt_PL.shp')
#DOI <- shapefile('DOI_AllTreatmentPoly.shp')
CSFS <- shapefile('CSFS_Stands_Master_13_17.shp')

#-> Relabel projection information so it plays well with Arc
proj <- '+proj=utm +zone=13 +ellps=GRS80 +datum=NAD83'
crs(USFS) <- CRS(proj)
#crs(DOI) <- CRS(proj)
crs(CSFS) <- CRS(proj)

#-> Load in classification tables
setwd(paste0(wd,'/INPUT/Spatial/Fuel_treatments')) #===> Change directory
tcw <- read.csv('Treatment_crosswalk_table.csv',header=T)

#-> Define Con function to match ArcGIS syntax
Con <- function(condition,trueValue,falseValue){
	return(condition*trueValue + (!condition)*falseValue)
}

#############################################END SET UP#############################################	

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/INPUT/Spatial/Fuel_treatments/Compiled'))) #===> Change directory

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

###---> USFS

#-> Subset to date range
# FISCAL_Y_1 = FISCAL_YEAR_COMPLETED
# FY_AWARDED = FISCAL_YEAR_AWARDED
# NOTE: being inclusive due to ambiguity in stewardship contracting
USFS$FY_COMP <- USFS$FISCAL_Y_1 # Rename for clarity
hc <- USFS[USFS$FY_COMP >= 2016,]
mc <- USFS[USFS$FY_COMP == 0 & USFS$FY_AWARDED >= 2016,]
USFS <- rbind(hc,mc)

#-> Clean up attributes
USFS$Agency <- 'USFS'
USFS$Year <- USFS$FY_COMP
USFS$Year <- ifelse(USFS$Year == 0,USFS$FY_AWARDED,USFS$Year)
USFS$Project <- USFS$NEPA_DOC_N
USFS$Unit <- NA
USFS$Treatment <- USFS$TREATMENT_
USFS <- USFS[,c('Agency','Year','Project','Unit','Treatment')]

#-> Subset to valid treatments
USFS <- USFS[!is.na(USFS$Treatment),]

#-> Cross walk to canopy and surface effects
tcws <- tcw[tcw$Agency=='USFS',c('Treatment','CanEff','SurfEff')]
USFS <- merge(USFS,tcws,by='Treatment',all.x=T)

# ###---> DOI

# #-> Subset to date range
# DOI <- DOI[DOI$YEAR >= 2016,]

# #-> Clean up attributes
# DOI$Agency <- DOI$BUREAUNAME
# DOI$Year <- DOI$YEAR
# DOI$Project <- DOI$PROJECTNAM 
# DOI$Unit <- NA
# DOI$Treatment <- DOI$TYPE_NAME
# DOI <- DOI[,c('Agency','Year','Project','Unit','Treatment')]

# #-> Subset to valid treatments
# DOI <- DOI[DOI$Treatment!='Chemical',]

# #-> Cross walk to canopy and surface effects
# tcws <- tcw[tcw$Agency=='DOI',c('Treatment','CanEff','SurfEff')]
# DOI <- merge(DOI,tcws,by='Treatment',all.x=T)

###---> CSFS
# NOTE: due to the inconsistent data, we used only reported silviculture and fire actions
# as indicators of thinning and prescribed fire, respectively.

#-> Subset to date range
CSFS$Fire <- suppressWarnings(as.numeric(CSFS$Fire))
CSFS$Fire[is.na(CSFS$Fire)] <- 0
CSFS$Silv <- suppressWarnings(as.numeric(CSFS$Silv))
CSFS$Silv[is.na(CSFS$Silv)] <- 0
CSFS <- CSFS[CSFS$Fire >= 2016 | CSFS$Silv >= 2016,]

#-> Clean up attributes
CSFS$Agency <- 'CSFS'
CSFS$Project <- NA
CSFS$Unit <- NA
CSFS$Year <- ifelse(CSFS$Fire >= CSFS$Silv,CSFS$Fire,CSFS$Silv)
CSFS$Treatment <- NA

#-> Assign effects
CSFS$CanEff <- ifelse(CSFS$Fire==0 & CSFS$Silv>0,'Thin',NA)
CSFS$CanEff <- ifelse(CSFS$Fire>0 & CSFS$Silv>0,'Thin and Rx fire',CSFS$CanEff)
CSFS$CanEff <- ifelse(CSFS$Fire>0 & CSFS$Silv==0,'Rx fire',CSFS$CanEff)
CSFS$SurfEff <- ifelse(CSFS$Fire==0,'None',NA)
CSFS$SurfEff <- ifelse(CSFS$Fire>0,'Rx fire',CSFS$SurfEff)

#-> Subset to fields
CSFS <- CSFS[,c('Agency','Year','Project','Unit','Treatment','CanEff','SurfEff')]

###---> Combine

#-> Merge and save before flattening
ftrts <- rbind(USFS,CSFS) # Removed DOI since empty
writeOGR(ftrts,dsn='.',layer='Compiled_wattrs_notflat',driver='ESRI Shapefile',overwrite=T)

#-> Dissolve by canopy treatment type (with safety checks)
if(nrow(ftrts[ftrts$CanEff=='Thin',]) > 0){
	ct_t <- gUnaryUnion(ftrts[ftrts$CanEff=='Thin',])
}else{
	ct_t <- NULL
}
if(nrow(ftrts[ftrts$CanEff=='Thin and Rx fire',]) > 0){
	ct_tRx <- gUnaryUnion(ftrts[ftrts$CanEff=='Thin and Rx fire',])
}else{
	ct_tRx <- NULL
}
if(nrow(ftrts[ftrts$CanEff=='Rx fire',]) > 0){
	ct_Rx <- gUnaryUnion(ftrts[ftrts$CanEff=='Rx fire',])
}else{
	ct_Rx <- NULL
}

#-> Enforce canopy treatment priorities
if(!(is.null(ct_t) | is.null(ct_tRx))){ # Complete treatment trumps all
	ct_t <- gDifference(ct_t,ct_tRx)
}
if(!(is.null(ct_Rx) | is.null(ct_tRx))){
	ct_Rx <- gDifference(ct_Rx,ct_tRx)
}
if(!(is.null(ct_t) | is.null(ct_Rx))){
	make_tRx <- intersect(ct_t,ct_Rx) # Combine thin and Rx fire
	if(length(make_tRx)>0){
		ct_t <- gDifference(ct_t,make_tRx)
		ct_Rx <- gDifference(ct_Rx,make_tRx)
		if(!is.null(ct_tRx)){
			ct_tRx <- gUnaryUnion(rbind(ct_tRx,make_tRx))
		}else{
			ct_tRx <- gUnaryUnion(make_tRx)
		}
	}
}

#-> Merge into flattened file (mainly for critique)
ct.l <- list()
if(!is.null(ct_t)){
	ct.l[[1]] <- SpatialPolygonsDataFrame(ct_t,data.frame(CanEff='Thin'))
}
if(!is.null(ct_tRx)){
	ct.l[[2]] <- SpatialPolygonsDataFrame(ct_tRx,data.frame(CanEff='Complete'))
}
if(!is.null(ct_Rx)){
	ct.l[[3]] <- SpatialPolygonsDataFrame(ct_Rx,data.frame(CanEff='RxFire'))
}
ct.l <- ct.l[lengths(ct.l) != 0]
ct <- do.call('rbind',ct.l)
writeOGR(ct,dsn='.',layer='CanEff_flattened',driver='ESRI Shapefile',overwrite=T)

#-> Dissolve by surface treatment type (with safety checks)
if(nrow(ftrts[ftrts$SurfEff=='Manage',]) > 0){
	st_man <- gUnaryUnion(ftrts[ftrts$SurfEff=='Manage',])
}else{
	st_man <- NULL
}
if(nrow(ftrts[ftrts$SurfEff=='Lop and scatter'|ftrts$SurfEff=='Masticate',]) > 0){
	st_rea <- gUnaryUnion(ftrts[ftrts$SurfEff=='Lop and scatter'|ftrts$SurfEff=='Masticate',])
}else{
	st_rea <- NULL
}
if(nrow(ftrts[ftrts$SurfEff=='Rx fire',]) > 0){
	st_Rx <- gUnaryUnion(ftrts[ftrts$SurfEff=='Rx fire',])
}else{
	st_Rx <- NULL
}

#-> Enforce surface treatment priorities
if(!(is.null(st_man) | is.null(st_Rx))){ # Rx fire trumps all
	st_man <- gDifference(st_man,st_Rx) 
}
if(!(is.null(st_rea) | is.null(st_Rx))){
	st_rea <- gDifference(st_rea,st_Rx)
}
if(!(is.null(st_rea) | is.null(st_man))){ # Manage trumps rearrange
	st_rea <- gDifference(st_rea,st_man) 
}

#-> Merge into flattened file (mainly for critique)
st.l <- list()
if(!is.null(st_man)){
	st.l[[1]] <- SpatialPolygonsDataFrame(st_man,data.frame(SurfEff='Manage'))
}
if(!is.null(st_rea)){
	st.l[[2]] <- SpatialPolygonsDataFrame(st_rea,data.frame(SurfEff='Rearrange'))
}
if(!is.null(st_Rx)){
	st.l[[3]] <- SpatialPolygonsDataFrame(st_Rx,data.frame(SurfEff='RxFire'))
}
st.l <- st.l[lengths(st.l) != 0]
st <- do.call('rbind',st.l)
writeOGR(st,dsn='.',layer='SurfEff_flattened',driver='ESRI Shapefile',overwrite=T)

############################################END ANALYSIS############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
