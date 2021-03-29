####################################################################################################
#### Exposure analysis
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 12/31/2019
#### Last Modified: 07/24/2020
####################################################################################################
# This is a generic workflow designed to calculate Expected Area Burned for a Highly Valued Resource
# or Asset (HVRA) with the option to bin by cNVC.
####################################################################################################
#-> Get working directory path from command call
initial.options <- commandArgs(trailingOnly=F)
setwd(dirname(sub('--file=','',initial.options[grep('--file=',initial.options)])))
#setwd('C:/Users/bgannon/Desktop/tRADS/scripts')
wd <- getwd()
pd <- paste(c(unlist(strsplit(wd,'/'))[1:(length(unlist(strsplit(wd,'/')))-1)],
            'Portable_R/Packages'),collapse='/')
####################################################################################################

###########################################START MESSAGE############################################
cat('Exposure analysis\n',sep='')
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

#-> Load in HVRA and relative importance tables
setwd(paste0(wd,'/INPUT')) #===> Change directory
rs <- read.csv('HVRA_settings.csv',header=T)
rs <- rs[rs$Include==1,]

#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
rae <- raster('ra_extent.tif')

#-> Load in burn probability
setwd(paste0(wd,'/INPUT/SPATIAL/Fire_simulation')) #===> Change directory
bp <- raster('BURNPROBABILITY.tif')

#-> Set path to cNVC rasters
inrasters <- paste0(wd,'/OUTPUT/cNVC/GIS')

#############################################END SET UP#############################################

#########################################START PROCESSING###########################################

#-> Clean up burn probability
bp[is.na(bp)] <- 0

#-> Correct to planning period burn probability
ppBP <- 1 - (1-bp)^25

#-> Limit to analysis extent
ppBP <- ppBP*rae

#-> Generic exposure function
# Default EAB unit is hectare
# HVRA is a binary raster (either 0/1 or NA/1)
# BP is a raster surface of burn probability
# c_raster is a co-variate raster to summarize EAB by
# c_bins are the breakpoint thresholds for the c_raster summary categories
exposure <- function(HVRA,BP,c_raster=NULL,c_bins=NULL){
	if(is.null(c_raster)){
		EAB <- cellStats(HVRA*BP*0.09,'sum',na.rm=T)
		return(EAB)
	}else{
		rcl <- data.frame(Low=c_bins[1:(length(c_bins)-1)],High=c_bins[-1],
		                  Class=1:(length(c_bins)-1))
		c_class <- reclassify(c_raster,rcl,include.lowest=T)
		EABt <- data.frame(zonal(BP*0.09,c_class,'sum',na.rm=T))
		colnames(EABt) <- c('Class','EAB_ha')
		EABt <- merge(EABt,rcl,by='Class',all.x=T)
		total <- sum(EABt$EAB_ha)
		EABt$Class <- NULL
		EABt$Class <- paste(EABt$Low,'to',EABt$High)
		EABt[nrow(EABt)+1,] <- NA
		EABt[nrow(EABt),'Class'] <- 'Total'
		EABt[nrow(EABt),'EAB_ha'] <- total
		return(EABt[,c('Class','EAB_ha')])
	}

}

#-> Create empty list to store EAB tables
EABts <- list()

#-> Process eNVC by category
for(i in 1:nrow(rs)){
	
	#-> Read in raster
	c_raster <- raster(paste0(inrasters,'/cNVC_',i,'.tif'))
	c_raster <- c_raster*rae
	
	#-> Convert to HVRA extent
	HVRA <- c_raster
	HVRA[!is.na(HVRA)] <- 1
	
	#-> Get exposure stats
	EABt <- exposure(HVRA,BP=ppBP,c_raster,c_bins=seq(-100,100,20))
	
	#-> Add HVRA name and category
	EABt$HVRA <- as.character(rs$HVRA[i])
	EABt$Category <- as.character(rs$Category[i])
	
	#-> Save to category list
	EABts[[i]] <- EABt[,c('HVRA','Category','Class','EAB_ha')]

}

#-> Combine list
EABt <- do.call('rbind',EABts)

#-> Save output
setwd(paste0(wd,'/OUTPUT')) #===> Change directory
write.csv(EABt,'Exposure_analysis.csv',row.names=F)

##########################################END PROCESSING############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
