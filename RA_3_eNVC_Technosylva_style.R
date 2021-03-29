####################################################################################################
#### Calculate expected net value change (eNVC) like Technosylva
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 05/07/2019
#### Last Modified: 02/24/2021
####################################################################################################
# This is a generic workflow designed to calculate eNVC from cNVC rasters, relative 
# importance weights, and relative extent. For presentation purposes, eNVC will be calculated
# first by HVRA "category" and then summed to calculate the total eNVC. It will also output
# a composite cNVC raster.
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
cat('Calculate expected net value change (eNVC)\n',sep='')
log <- file(paste0(wd,'/OUTPUT/eNVC/RA3.log'))
sink(file=log,append=T,type=c('output','message'),split=T)
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

#-> Load in input tables
setwd(paste0(wd,'/INPUT')) #===> Change directory
rs <- read.csv('HVRA_settings.csv',header=T)
rs <- rs[rs$Include==1,]
ms <- read.csv('map_settings.csv',header=T)
ri <- read.csv('relative_importance.csv',header=T)

#-> Load in feature classes
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='mapping_extent',verbose=F)
ra_extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='risk_analysis_extent',verbose=F)
roads_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_roads',verbose=F)
rivers_fc <- suppressWarnings(readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_streams',verbose=F))

#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
ma_extent <- raster('ma_extent.tif')
ra_extent <- raster('ra_extent.tif')

#-> Set path to cNVC rasters
inrasters <- paste0(wd,'/OUTPUT/cNVC/GIS')

#-> Load in relative extent raster
setwd(paste0(wd,'/OUTPUT/HVRA_extents')) #===> Change directory
re <- read.csv('relative_area.csv',header=T)

#-> Load in burn probability
setwd(paste0(wd,'/INPUT/SPATIAL/Fire_simulation')) #===> Change directory
bp <- raster('BURNPROBABILITY.tif')

#-> Load in hillshade
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
hillshade <- raster('hillcl.tif')
hillshade <- crop(hillshade,extent_fc)
SoG <- colorRampPalette(c('grey0','grey100'))(50)

#############################################END SET UP#############################################

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/OUTPUT/eNVC'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Empty output folder
	d_list <- list.files()
	if(length(d_list) > 0){
		unlink(d_list[grep('RA3.log',d_list,invert=T)],recursive=T)
	}
}

#-> Add GIS subdirectory
dir.create('GIS')

##########################################END CLEAR OUTPUT##########################################

#####################################START RELATIVE IMPORTANCE######################################
# Note that no assumptions are made that the weights add up to 100%. Instead, the category
# RIs are coverted to relative measures that sum to one. Then, the HVRA RIs are relativized
# by category and they are combined with multiplication.

#-> Check that all relative importance categories are represented
ri_ri <- unique(ri$Category)
ri_rs <- unique(rs$Category)
extra_ri <- ri_ri[!(ri_ri %in% ri_rs)]
extra_rs <- ri_rs[!(ri_rs %in% ri_ri)]
if(length(extra_ri) > 0){
	if(length(extra_ri)==1){
		cat('Warning:',paste0(extra_ri),
		    'is assigned relative importance but not included.\n')
	}
	if(length(extra_ri)>1){
		cat('Warning:',paste(extra_ri,collapse=', '),
		    'are assigned relative importance but not included.\n')
	}
}
if(length(extra_rs) > 0){
	if(length(extra_rs)==1){
		cat('Warning:',paste0(extra_rs),
		    'is included but not assigned relative importance.\n')
	}
	if(length(extra_rs)>1){
		cat('Warning:',paste(extra_rs,collapse=', '),
		    'are included but not assigned relative importance.\n')
	}
}

#-> Calculate Relative Importance/Relative Extent weights
ri$RI_cat <- ri$RI/sum(ri$RI) # Calculate total relative importance for category
w <- rs[,c('Layer','Category','RI_HVRA')] # Subset important information from rs table
ct <- ddply(w,.(Category),summarize,
            CatTot = sum(RI_HVRA))
w <- merge(w,ct,by='Category',all.x=T)	 
w$RI_HVRA <- w$RI_HVRA/w$CatTot
w$CatTot <- NULL
w <- merge(w,ri[,c('Category','RI_cat')],by='Category',all.x=T)
w$RI <- w$RI_HVRA*w$RI_cat
w <- merge(w,re[,c('Layer','RE')],by='Layer',all.x=T)
w$RIoRE <- w$RI/w$RE
rs <- merge(rs,w[,c('Layer','RIoRE')],by='Layer',all.x=T)

#-> Save outputs for review
kfields <- c('Layer','HVRA','Category','RI_HVRA','RIoRE')
write.csv(rs[,kfields],'RelImp_over_RelExt_check.csv',row.names=F)

######################################END RELATIVE IMPORTANCE#######################################

#########################################START PROCESSING###########################################

#-> Clean up burn probability
bp[is.na(bp)] <- 0

#-> Get vector of HVRA categories
cats <- unique(rs$Category)

#-> Create empty list to store category-level cNVC
cNVCs <- list()

#-> Create empty list to store risk totals
rTots <- list()

#-> Process eNVC by category
setwd(paste0(wd,'/OUTPUT/eNVC/GIS')) #===> Change directory
for(i in 1:length(cats)){
	
	cat('Calculating eNVC for ',paste(cats[i]),'\n',sep='')
	
	#-> Subset table
	catrs <- rs[rs$Category==paste(cats[i]),]
	
	#-> Create empty list to store HVRA-level eNVC
	hvra.l <- list()
	
	#-> Create list to store coverage rasters
	cov.l <- list()
	
	#-> Calculate weighted eNVC for each HVRA
	for(j in 1:nrow(catrs)){
		
		#-> Read in raster
		cNVC <- raster(paste0(inrasters,'/cNVC_',catrs$Layer[j],'.tif'))
		
		#-> Save coverage
		coverage <- cNVC; coverage[!is.na(coverage)] <- 1; coverage[is.na(coverage)] <- 0
		cov.l[[j]] <- coverage
		
		#-> Fill null with zero
		cNVC[is.na(cNVC)] <- 0
		
		#-> Calculate cNVC with relative importance adjustment
		# RIoRE = relative importance/relative extent
		hvra.l[[j]] <- cNVC*catrs$RIoRE[j]
		
		#-> Get total risk associated with HVRA
		rTot <- catrs[j,c('Layer','HVRA','Category')]
		rTot$Tot_cNVC <- cellStats(cNVC*catrs$RIoRE[j]*ra_extent,'sum')
		rTot$Tot_eNVC <- cellStats(bp*cNVC*catrs$RIoRE[j]*ra_extent,'sum')
		rTots[[length(rTots)+1]] <- rTot
	
	}
	
	if(length(hvra.l) > 1){
	
		#-> Calculate cNVC for category
		cNVC <- sum(stack(hvra.l))
		
		#-> Get coverage
		coverage <- sum(stack(cov.l))
		
	}else{
		
		cNVC <- hvra.l[[1]]
		coverage <- cov.l[[1]]
		
	}
	
	#-> Set cNVC raster to null where there is no coverage
	cNVC[coverage==0] <- NA
	
	#-> Save category raster to file
	writeRaster(suppressWarnings(cNVC*bp*ma_extent),paste0('eNVC_category_',i,'.tif'),
	            format='GTiff',overwrite=T)
	
	#-> Save to category list
	cNVCs[[i]] <- cNVC

}

#-> Calculate integrated cNVC across categories
if(length(cNVCs) > 1){
	ts <- stack(cNVCs) # Stack rasters
	ts[is.na(ts)] <- 0	# Fill nulls with zeros
	cNVC <- sum(ts)
}else{
	ts <- cNVCs[[1]]
	ts[is.na(ts)] <- 0	# Fill nulls with zeros
}
	
#-> Save composite cNVC raster
writeRaster(cNVC*ma_extent,'Total_cNVC.tif',format='GTiff',overwrite=T)

#-> Save composite eNVC raster
writeRaster(suppressWarnings(cNVC*bp*ma_extent),'Total_eNVC.tif',format='GTiff',overwrite=T)

#-> Save category key
Number <- seq(1,length(cats),1)
ckey <- data.frame(Number,Category=cats)
write.csv(ckey,'category_key.csv',row.names=F)

###---> Compile and save risk report
setwd(paste0(wd,'/OUTPUT/eNVC')) #===> Change directory

#-> Compile list to data frame
rdf <- do.call('rbind',rTots)

#-> Save
write.csv(rdf,'Total_risk_report.csv',row.names=F)

##########################################END PROCESSING############################################

##########################################START GRAPHICS############################################
setwd(paste0(wd,'/OUTPUT/eNVC')) #===> Change directory

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

#-> Trim dead legend space with image magick
# allows use of raster plot instead of image to get better resolution
trimLegend <- function(fname,baseX,baseY,adjFact){
	ri <- image_read(fname)
	ic <- image_crop(ri,paste0(ceiling(baseX*adjFact),'x',baseY))  
	image_write(ic,fname)
}

#-> Organize rasters
cNVCs[[length(cNVCs)+1]] <- cNVC

#-> Crop reference data
hillshade <- crop(hillshade,extent_fc)
roads_fc <- gIntersection(roads_fc,extent_fc)
rivers_fc <- gIntersection(rivers_fc,extent_fc)

#-> Make pretty maps
titles <- c(as.character(ckey$Category),'Composite Wildfire Risk')
bins <- c(-1000,-1,-0.1,-0.01,-0.001,0.001,0.01,0.1,1,1000)
llabels <- c('Negative','','','','Neutral','','','','Positive')
cols <- colorRampPalette(unlist(strsplit(as.character(ms$NVC_cols),',')))(length(bins)-1)

for(i in 1:length(cNVCs)){
	
	#-> Load in raster, clip to extent, and apply BP
	eNVC <- suppressWarnings(cNVCs[[i]]*ma_extent*bp)
	
	#-> Mask to extent
	eNVC <- mask(eNVC,extent_fc)
	
	#-> Map
	tiff(paste0(gsub(' ','_',titles[i]),'_eNVC.tif'),width=ms$Width_pixels,height=ms$Height_pixels,
	     pointsize=ms$Pointsize,compression='lzw',type='windows')
	par(mar=c(0.6,0.6,2.5,0.6))
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	plot(eNVC,breaks=bins,col=cols,alpha=0.6,axes=F,maxpixels=5000000,legend=F,box=F,add=T)
	title(main=paste(titles[i]),cex.main=1.5)
	plot(rivers_fc,col='blue',lwd=2,add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(ra_extent_fc,border='black',lwd=2,add=T)
	plot(extent_fc,border='grey20',add=T)
	plotRef(bb=extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	legend(as.character(ms$Legend_pos),llabels,fill=cols,title=expression(bold('eNVC')),bty='n',
	       cex=0.75)
	g <- dev.off()
	trimLegend(paste0(gsub(' ','_',titles[i]),'_eNVC.tif'),ms$Width_pixels,ms$Height_pixels,0.925)

}	

#-> Save composite risk to kml
# peNVC <-  projectRaster(eNVC,crs=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'),method='bilinear')
# KML(peNVC,'eNVC.kml',breaks=bins,col=cols,maxpixels=5000000,blur=1,overwrite=T) 

###---> Composite cNVC map

#-> Define bins
#bins <- c(-1000000,-10000,-1000,-100,-10,10,100,1000,10000,1000000)
bins <- c(-1000000,-1000,-100,-10,-1,1,10,100,1000,1000000)

#-> Mask to extent
cNVC <- mask(cNVC,extent_fc)

#-> Map
tiff('Composite_Wildfire_cNVC.tif',width=ms$Width_pixels,height=ms$Height_pixels,
	     pointsize=ms$Pointsize,compression='lzw',type='windows')
par(mar=c(0.6,0.6,2.5,0.6))
plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
plot(cNVC,breaks=bins,col=cols,alpha=0.6,axes=F,maxpixels=5000000,legend=F,box=F,add=T)
title(main='Composite Conditional Net Value Change',cex.main=1.5,adj=0.25)
plot(rivers_fc,col='blue',lwd=2,add=T)
plot(roads_fc,col='grey20',lwd=1.5,add=T)
plot(ra_extent_fc,border='black',lwd=2,add=T)
plot(extent_fc,border='grey20',add=T)
plotRef(bb=extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
legend(as.character(ms$Legend_pos),llabels,fill=cols,title=expression(bold('cNVC')),bty='n',
       cex=0.75)
g <- dev.off()
trimLegend('Composite_Wildfire_cNVC.tif',ms$Width_pixels,ms$Height_pixels,0.925)

###########################################END GRAPHICS#############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
sink()
cat('Contents saved to log file.\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
