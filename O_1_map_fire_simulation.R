####################################################################################################
#### Map fire simulation output
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 11/20/2019
#### Last Modified: 07/13/2020
####################################################################################################
# Makes simple maps of fire model output by fuel moisture/weather scenario.
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
cat('Map fire simulation output\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages
.libPaths(pd)
packages <- c('raster','rgdal','rgeos','plyr','magick')
for(package in packages){
	if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
		install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,lib.loc=pd,character.only=T))
	}
}

#-> Maximize raster processing speed
rasterOptions(maxmemory = 10^9)

#-> Load in mapping settings
setwd(paste0(wd,'/INPUT')) #===> Change directory
ms <- read.csv('map_settings.csv',header=T)

#-> Load in feature classes
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='risk_analysis_extent',verbose=F)
roads_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_roads',verbose=F)
rivers_fc <- suppressWarnings(readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_streams',verbose=F))
							  
#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
rextent <- raster('ra_extent')

#-> Load in conditional fire intensity rasters
setwd(paste0(wd,'/INPUT/SPATIAL/Fire_simulation')) #===> Change directory
BP <- raster('BP_Short')
FL <- list()
FL[[1]] <- raster('baseline_25_FLAMELENGTH.tif')*3.28084 # Convert to feet
FL[[2]] <- raster('baseline_50_FLAMELENGTH.tif')*3.28084
FL[[3]] <- raster('baseline_90_FLAMELENGTH.tif')*3.28084
FL[[4]] <- raster('baseline_97_FLAMELENGTH.tif')*3.28084
CFA <- list()
CFA[[1]] <- raster('baseline_25_CROWNSTATE.tif')
CFA[[2]] <- raster('baseline_50_CROWNSTATE.tif')
CFA[[3]] <- raster('baseline_90_CROWNSTATE.tif')
CFA[[4]] <- raster('baseline_97_CROWNSTATE.tif')

#-> Load in hillshade
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
hillshade <- raster('hillcl.tif')
hillshade <- crop(hillshade,extent_fc)
hillshade <- mask(hillshade,extent_fc)
SoG <- colorRampPalette(c('grey0','grey100'))(50)

#############################################END SET UP#############################################

###########################################START GRAPHICS###########################################
setwd(paste0(wd,'/OUTPUT/Fire_simulation'))	#===> Change directory

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

###---> Burn probability

#-> Load in and crop data
BP <- BP*rextent

#-> Make map
rmax <- ceiling(cellStats(BP,'max',na.rm=T)*1000)/1000
bins <- rev(c(rmax,rmax/c(2,4,8,16,32,64),0)) # Even seq(0,rmax,rmax/7)
llabels <- c('Low: 0','','','Moderate','','',paste('High:',round(rmax,3)))
cols <- colorRampPalette(c('blue','cornsilk','red'))(length(bins)-1)
tiff('Burn_probability.tif',width=ms$Width_pixels,height=ms$Height_pixels,pointsize=ms$Pointsize,
     compression='lzw',type='windows')
par(mar=c(0.6,0.6,2.5,0.6))
plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
plot(BP,breaks=bins,col=cols,alpha=0.6,axes=F,maxpixels=5000000,legend=F,box=F,add=T)
title(main='Burn Probability',cex.main=1.75)
plot(rivers_fc,col='blue',add=T)
plot(roads_fc,col='grey20',lwd=1.5,add=T)
plot(extent_fc,border='black',add=T)
plotRef(bb=extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
legend(as.character(ms$Legend_pos),llabels,fill=cols,title=expression(bold('Mean Annual BP')),
       bty='n',cex=0.8)	      
g <- dev.off()

###---> Flame lengths

titles <- c('Flame Length - Low Scenario','Flame Length - Moderate Scenario',
            'Flame Length - High Scenario','Flame Length - Extreme Scenario')

for(i in 1:length(titles)){

	#-> Crop
	sFL <- FL[[i]]*rextent
	
	#-> Bins and labels
	bins <- c(0,2,4,6,8,12,1000)
	llabels <- c('0-2','2-4','4-6','6-8','8-12','> 12')
	cols <- colorRampPalette(c('forestgreen','yellow','red'))(length(bins)-1)
	
	tiff(paste0('FL_',i,'.tif'),width=ms$Width_pixels,height=ms$Height_pixels,
	     pointsize=ms$Pointsize,compression='lzw',type='windows')
	
	#-> Map
	par(mar=c(0.6,0.6,2.5,1.6)) # Larger right margin for histograms
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	plot(sFL,breaks=bins,col=cols,alpha=0.6,axes=F,maxpixels=5000000,legend=F,box=F,add=T)
	title(main=titles[i],cex.main=1.75,adj=0.05)
	plot(rivers_fc,col='blue',add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(extent_fc,border='black',add=T)
	plotRef(bb=extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	
	#-> Histogram legend
	counts <- hist(sFL,breaks=bins,plot=F)$count
	rcounts <- 100*(counts/sum(counts))
	par(fig=c(0.70,1,0.80,0.95),mar=c(1,2.1,1,0.6),mgp=c(0.8,0.2,0),tck=-0.05,new=T,
		    cex=0.8,cex.axis=0.6,cex.lab=0.6,cex.main=0.8)
	barplot(rcounts,ylim=c(0,100),ylab='Frequency (%)',main='Flame Length (ft)',names=llabels,
	        col=cols)
			
	g <- dev.off()
	
}

###---> Crown fire activity

titles <- c('Crown Fire Activity - Low Scenario','Crown Fire Activity - Moderate Scenario',
            'Crown Fire Activity - High Scenario','Crown Fire Activity - Extreme Scenario')

for(i in 1:length(titles)){

	#-> Crop
	sCFA <- CFA[[i]]*rextent
	
	#-> Bins and labels
	bins <- seq(-0.5,3.5,1)
	llabels <- c('Unburned','Surface','Passive','Active')
	cols <- colorRampPalette(c('forestgreen','yellow','orange','red'))(length(bins)-1)
	
	tiff(paste0('CFA_',i,'.tif'),width=ms$Width_pixels,height=ms$Height_pixels,
	     pointsize=ms$Pointsize,compression='lzw',type='windows')
	
	#-> Map
	par(mar=c(0.6,0.6,2.5,1.6)) # Larger right margin for histograms
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	plot(sCFA,breaks=bins,col=cols,alpha=0.6,axes=F,maxpixels=5000000,legend=F,box=F,add=T)
	title(main=titles[i],cex.main=1.75,adj=0.05)
	plot(rivers_fc,col='blue',add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(extent_fc,border='black',add=T)
	plotRef(bb=extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	
	#-> Histogram legend
	counts <- hist(sCFA,breaks=bins,plot=F)$count
	rcounts <- 100*(counts/sum(counts))
	par(fig=c(0.70,1,0.80,0.95),mar=c(1,2.1,1,0.6),mgp=c(0.8,0.2,0),tck=-0.05,new=T,
		    cex=0.8,cex.axis=0.5,cex.lab=0.6,cex.main=0.8)
	barplot(rcounts,ylim=c(0,100),ylab='Frequency (%)',main='Crown Fire Activity',names=llabels,
	        col=cols)
			
	g <- dev.off()
	
}

############################################END GRAPHICS############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
