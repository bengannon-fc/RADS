####################################################################################################
#### Process spatial data
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 05/06/2019
#### Last Modified: 02/24/2021
####################################################################################################
# This is a generic spatial analysis workflow designed to read in vector or raster spatial
# data on HVRA locations and convert them into a standardized raster format. It will also 
# calculate the area of each resource needed for the relative importance weighting.
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
log <- file(paste0(wd,'/OUTPUT/HVRA_extents/RA1.log'))
sink(file=log,append=T,type=c('output','message'),split=T)
cat('Process spatial data\n',sep='')
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

#-> Load in HVRA table
setwd(paste0(wd,'/INPUT')) #===> Change directory
rs <- read.csv('HVRA_settings.csv',header=T)
rs <- rs[rs$Include==1,]
ms <- read.csv('map_settings.csv',header=T)

#-> Load in feature classes
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
pextent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='risk_analysis_extent',verbose=F)
extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='mapping_extent',verbose=F)
ra_extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='risk_analysis_extent',verbose=F)
roads_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_roads',verbose=F)
rivers_fc <- suppressWarnings(readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_streams',verbose=F))
contours_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_contours',verbose=F)

#-> Set paths to geodatabase and raster folder
ingdb <- paste0(wd,'/INPUT/SPATIAL/VECTOR_INPUT.gdb')
inrasters <- paste0(wd,'/INPUT/SPATIAL/HVRA_rasters')

#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
ma_extent <- raster('ma_extent.tif') # Extent for mapping w/ buffer around analysis area
ra_extent <- raster('ra_extent.tif') # Risk analysis area for relative extent calcs

#-> Load in extent feature classes
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory														
extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='mapping_extent',verbose=F)

#############################################END SET UP#############################################

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/OUTPUT/HVRA_extents'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Empty output folder
	d_list <- list.files()
	if(length(d_list) > 0){
		unlink(d_list[grep('RA1.log',d_list,invert=T)],recursive=T)
	}
}

#-> Add GIS subdirectory
dir.create('GIS')

##########################################END CLEAR OUTPUT##########################################

##########################################START PROCESSING##########################################

#-> Define plot reference function
# This function will add a scale bar and north arrow to map
# bb = bounding box polygon (does not have to be box)
# pos = position for reference information (allowable options are bottomright and bottomleft
# blength = scale bar length in miles
plotRef <- function(bb,pos,blength){
	#-> Scale bar
	if(pos=='bottomleft'){
		sbx <- xmin(bb) + 0.01*(xmax(bb)-xmin(bb))
	}
	if(pos=='bottomright'){
		xprop <- ((xmax(bb)-xmin(bb)) - (blength*1600))/(xmax(bb)-xmin(bb)) - 0.01
		sbx <- xmin(bb) + xprop*(xmax(bb)-xmin(bb))
	}
	sby <- ymin(bb) + 0.01*(ymax(bb)-ymin(bb))
	scalebar((blength*1600),xy=c(sbx,sby),type='line',divs=2,lwd=2,label=paste(blength,'miles'))
	#-> North arrow
	if(pos=='bottomright'){
		sbx <- xmin(bb) + 0.99*(xmax(bb)-xmin(bb))
	}
	y1 <- ymin(bb) + 0.025*(ymax(bb)-ymin(bb))
	y2 <- ymin(bb) + 0.075*(ymax(bb)-ymin(bb))
	arrows(sbx,y1,sbx,y2,length=0.1,lwd=2)	
}

#-> Create table to store "relative" area for relative importance weighting
ra <- rs[,c('Layer','HVRA')]
ra$Area_ac <- NA

#-> Crop reference data
roads_fc <- gIntersection(roads_fc,extent_fc)
rivers_fc <- gIntersection(rivers_fc,extent_fc)
contours_fc <- gIntersection(contours_fc,extent_fc)

#-> Start map book
setwd(paste0(wd,'/OUTPUT/HVRA_extents')) #===> Change directory
pdf('HVRA_map_book.pdf',height=(ms$Height_pixels/200),width=(ms$Width_pixels/200))

#-> Processing loop
for(i in 1:nrow(rs)){
	
	#-> Print message
	cat(paste0('Processing ',rs$HVRA[i],'\n'))
	
	###---> Process vector data (if needed)
	if(rs$Type[i] %in% c('Polygon','Polyline','Point')){
		
		#-> Read in feature class [suppressWarnings() used for z dimensions]
		fc <- suppressWarnings(readOGR(dsn=ingdb,layer=paste0(rs$FeatureClass[i]),verbose=F))							   
		
		#-> Buffer (if needed)
		if(rs$Buffer_m[i]>0){
			fc <- gBuffer(fc,width=paste(rs$Buffer_m[i]),byid=T)
		}
				
		#-> Add zone field
		fc$zone <- 1
		
		#-> Convert to raster
		zone <- rasterize(fc,ma_extent,'zone')
		
		#-> Clip to extent
		zone <- zone*ma_extent
		
	}	
	
	###---> Process raster data (if needed)
	if(rs$Type[i]=='Raster'){
		
		if(rs$Represents[i]=='Location'){
			
			#-> Read in raster
			zone <- raster(paste0(inrasters,'/',rs$FeatureClass[i]))
			
			#-> Standardize zero to null
			zone[zone==0] <- NA
			
		}
		
		if(rs$Represents[i]=='cNVC'){
			
			#-> Read in raster
			zone <- raster(paste0(inrasters,'/',rs$FeatureClass[i]))
			
			#-> Convert to zone format
			zone[!is.na(zone)] <- 1
			
		}
		
		#-> Crop to match extent
		zone <- crop(zone,ma_extent)
			
		#-> Mask with extent
		zone <- zone*ma_extent
	}

	#-> Save raster
	setwd(paste0(wd,'/OUTPUT/HVRA_extents/GIS')) #===> Change directory
	writeRaster(zone,paste0('HVRA_',rs$Layer[i],'.tif'),format='GTiff',overwrite=T)
	
	#-> Store area WITHIN Risk Analysis Area
	ra$Area_ac[i] <- cellStats(zone*ra_extent,'sum',na.rm=T)*
	                           (((res(zone)[1]^2)/10000)*2.47105)
	
	#-> Make simple map
	plot(zone,col='green',axes=F,legend=F,box=F,main=paste(rs$HVRA[i]))
	plot(contours_fc,col='grey80',lwd=0.5,add=T)
	plot(rivers_fc,col='blue',add=T)
	plot(roads_fc,col='red',add=T)
	plot(pextent_fc,border='black',lwd=1.5,add=T)
	plot(extent_fc,border='black',add=T)
	if(rs$Type[i]=='Polyline' | rs$Type[i]=='Point'){
		plot(fc,col='green',border='green',lwd=1,add=T)
	}
	plotRef(bb=extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)

}	

g <- dev.off()

#-> Save area
setwd(paste0(wd,'/OUTPUT/HVRA_extents')) #===> Change directory
ra$RE <- ra$Area_ac/sum(ra$Area_ac)
write.csv(ra,'relative_area.csv',row.names=F)

###########################################END PROCESSING###########################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
sink()
cat('Contents saved to log file.\n',sep='')
cat('Close command window to proceed!\n',sep='')
#############################################END MESSAGE############################################
