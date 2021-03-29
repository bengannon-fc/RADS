####################################################################################################
#### Calculate conditional net value change (cNVC) like Technosylva
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 05/07/2019
#### Last Modified: 08/11/2020
####################################################################################################
# This is a generic workflow designed to calculate cNVC by HVRA from a table of response
# functions. Instead of using conditional flame length probability rasters, this uses flame
# lengths modeled for four fire weather scenarios (25th, 50th, 90th, and 97th percentiles). 
# cNVC is calculated for each scenario and then a weighted mean is calculated 
# (weights = 0.01, 0.09. 0.2, 0.7). Alternatively, it copies a cNVC raster to the proper
# folder if provided one.
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
cat('Calculate conditional net value change (cNVC)\n',sep='')
log <- file(paste0(wd,'/OUTPUT/cNVC/RA2.log'))
sink(file=log,append=T,type=c('output','message'),split=T)
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

#-> Load in input tables
setwd(paste0(wd,'/INPUT')) #===> Change directory
rs <- read.csv('HVRA_settings.csv',header=T)
rs <- rs[rs$Include==1,]
ms <- read.csv('map_settings.csv',header=T)

#-> Load in feature classes
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='mapping_extent',verbose=F)
ra_extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='risk_analysis_extent',verbose=F)
roads_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_roads',verbose=F)
rivers_fc <- suppressWarnings(readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_streams',verbose=F))

#-> Set path to input rasters folder for any user-provided cNVCs
incNVC <- paste0(wd,'/INPUT/SPATIAL/HVRA_rasters')

#-> Set path to HVRA binary rasters
inrasters <- paste0(wd,'/OUTPUT/HVRA_extents/GIS')

#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
ma_extent <- raster('ma_extent.tif')

#-> Load in flame length rasters
setwd(paste0(wd,'/INPUT/SPATIAL/Fire_simulation')) #===> Change directory
FL25 <- raster('baseline_25_FLAMELENGTH.tif')*3.28084 # Convert to feet
FL50 <- raster('baseline_50_FLAMELENGTH.tif')*3.28084
FL90 <- raster('baseline_90_FLAMELENGTH.tif')*3.28084
FL97 <- raster('baseline_97_FLAMELENGTH.tif')*3.28084

#-> Load in hillshade
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
hillshade <- raster('hillcl.tif')
hillshade <- crop(hillshade,extent_fc)
SoG <- colorRampPalette(c('grey0','grey100'))(50)

#############################################END SET UP#############################################

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/OUTPUT/cNVC'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Empty output folder
	d_list <- list.files()
	if(length(d_list) > 0){
		unlink(d_list[grep('RA2.log',d_list,invert=T)],recursive=T)
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
	scalebar((blength*1600),xy=c(sbx,sby),type='line',divs=2,lwd=8,label=paste(blength,'miles'))
	#-> North arrow
	if(pos=='bottomright'){
		sbx <- xmin(bb) + 0.99*(xmax(bb)-xmin(bb))
	}
	y1 <- ymin(bb) + 0.025*(ymax(bb)-ymin(bb))
	y2 <- ymin(bb) + 0.075*(ymax(bb)-ymin(bb))
	arrows(sbx,y1,sbx,y2,length=0.1,lwd=8)	
}

#-> Crop reference data
hillshade <- crop(hillshade,extent_fc)
roads_fc <- gIntersection(roads_fc,extent_fc)
rivers_fc <- gIntersection(rivers_fc,extent_fc)

#-> Process each HVRA
for(i in 1:nrow(rs)){
	
	#-> Print message
	cat(paste0('Calculating cNVC for ',rs$HVRA[i],'\n'))
	
	#-> Read in raster
	zone <- raster(paste0(inrasters,'/HVRA_',rs$Layer[i],'.tif'))
	
	if(rs$Represents[i]=='Location'){
	
		#-> Prep response function
		FILbins <- c(0,2,4,6,8,12,1000)
		rf <- as.numeric(rs[i,c('FIL1','FIL2','FIL3','FIL4','FIL5','FIL6')])
		rcl <- data.frame(from=FILbins[1:6],to=FILbins[2:7],rf)
		
		#-> Calculate cNVC by scenario
		cNVC25 <- reclassify(FL25,rcl)
		cNVC50 <- reclassify(FL50,rcl)
		cNVC90 <- reclassify(FL90,rcl)
		cNVC97 <- reclassify(FL97,rcl)
		
		#-> Weight scenarios and limit to extent
		cNVC <- zone*(0.01*cNVC25 + 0.09*cNVC50 + 0.2*cNVC90 + 0.7*cNVC97)
		
	}
	
	if(rs$Represents[i]=='cNVC'){
		
		#-> Read in user-provided cNVC raster
		cNVC <- raster(paste0(incNVC,'/',rs$FeatureClass[i]))
		
		#-> Crop to match extent
		cNVC <- crop(cNVC,zone)
		
		#-> Add function to rescale non-compliant cNVC measures
		# Linear stretch with percentile clip
		if(cellStats(cNVC,'min') < -100 | cellStats(cNVC,'max') > 100){
			nzvals <- abs(values(cNVC))
			nzvals <- nzvals[!is.na(nzvals) & nzvals!=0]
			tval <- quantile(nzvals,0.975)
			cNVC[cNVC < tval*(-1)] <- tval*(-1) 
			cNVC[cNVC > tval] <- tval
			cNVC <- cNVC*(100/tval)
		}
		
	}
	
	#-> Save raster
	setwd(paste0(wd,'/OUTPUT/cNVC/GIS')) #===> Change directory
	writeRaster(cNVC,paste0('cNVC_',rs$Layer[i],'.tif'),format='GTiff',overwrite=T)
	
	#-> Mask to project area
	cNVC <- mask(cNVC,extent_fc)
	
	#-> Make map
	setwd(paste0(wd,'/OUTPUT/cNVC')) #===> Change directory
	tiff(paste0('cNVC_',rs$Layer[i],'.tif'),width=ms$Width_pixels,height=ms$Height_pixels,
	     pointsize=ms$Pointsize,compression='lzw',type='windows')
	bins <- seq(-100,100,10)
	llabels <- c('-100 to -90','-90 to -80','-80 to -70','-70 to -60','-60 to -50','-50 to -40',
	             '-40 to -30','-30 to -20','-20 to -10','-10 to 0','0 to 10','10 to 20','20 to 30',
	             '30 to 40','40 to 50','50 to 60','60 to 70','70 to 80','80 to 90','90 to 100')
	cols <- colorRampPalette(unlist(strsplit(as.character(ms$NVC_cols),',')))(length(bins)-1)
	par(mar=c(0.6,0.6,2.5,1.6)) # Larger right margin for response function barplot
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	plot(cNVC,breaks=bins,col=cols,alpha=0.6,axes=F,maxpixels=5000000,legend=F,box=F,add=T)
	title(main=paste(rs$HVRA[i]),cex.main=1.75)
	plot(rivers_fc,col='blue',add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(ra_extent_fc,border='black',lwd=2,add=T)
	plot(extent_fc,border='grey20',add=T)
	plotRef(bb=extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	legend(as.character(ms$Legend_pos),llabels,fill=cols,title=expression(bold('cNVC')),bty='n',
	       cex=0.75)
	
	#-> Response function
	names <- c('0-2','2-4','4-6','6-8','8-12','>12')
	rf <- as.numeric(rs[i,c('FIL1','FIL2','FIL3','FIL4','FIL5','FIL6')])
	if(!is.na(max(rf))){
		par(fig=c(0.70,1,0.80,0.95),mar=c(1,2.1,1,0.6),mgp=c(0.8,0.2,0),tck=-0.05,new=T,
		    cex=0.8,cex.axis=0.6,cex.lab=0.6,cex.main=0.8)
		barplot(rf,ylim=c(-100,100),col='grey',border=NA,ylab='Value Change',
		        main='Response Function')
		axis(1,at=(0.7+1.2*seq(0,5)),labels=names,tick=F,line=-0.5)
	}
	
	g <- dev.off()

}	

###########################################END PROCESSING###########################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
sink()
cat('Contents saved to log file.\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
