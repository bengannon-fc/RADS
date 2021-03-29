####################################################################################################
#### Spatial optimization maps
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 10/21/2019
#### Last Modified: 08/11/2020
####################################################################################################
# Summary: this completes spatial analyses and mapping to complement the optimization script. Some
# of it is repetitive to avoid unecessary file creation. This workflow is separated from the
# optimization so as not to slow it down.
#
# Data sources:
# Raster inputs are controlled by treatment in the so_treatment_specs.csv table
# Treatment units are provided in raster and shapefile form to avoid unnecessary processing
# Budgets are read in from the so_budgets.csv file
# Min and max project size for each treatment unit is controlled with the so_project_size.csv file
# Many other data sources are pulled in from elsewhere using relative path names. They are not all 
# listed in the set up section because some are meant to be dynamic (flexible to provided treatment 
# type).
#
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
cat('Spatial optimization maps\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages
.libPaths(pd)
packages <- c('raster','rgdal','rgeos','magick')
for(package in packages){
	if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
		install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,lib.loc=pd,character.only=T))
	}
}

#-> Maximize raster processing speed
rasterOptions(maxmemory = 10^9)

#-> Load in feature classes
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
ra_extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='risk_analysis_extent',verbose=F)
roads_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_roads',verbose=F)
rivers_fc <- suppressWarnings(readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_streams',verbose=F))
							  
#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
ra_extent <- raster('ra_extent.tif')

#-> Load in hillshade
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
hillshade <- raster('hillcl.tif')
hillshade <- crop(hillshade,ra_extent_fc)
hillshade <- mask(hillshade,ra_extent_fc)
SoG <- colorRampPalette(c('grey0','grey100'))(50)

#-> Load in run specs
setwd(paste0(wd,'/INPUT')) #===> Change directory
tspecs <- read.csv('so_treatment_specs.csv',header=T)
budgets <- read.csv('so_budgets.csv',header=T)
ms <- read.csv('map_settings.csv',header=T)

#############################################END SET UP#############################################

########################################START CLEAR OUTPUT##########################################

e <- try(setwd(paste0(wd,'/OUTPUT/Treatment_plan'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Add GIS subdirectory
	dir.create('Maps')
}

#########################################END CLEAR OUTPUT###########################################

##########################################START GRAPHICS############################################

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

#-> Log bins
logBins <- function(xmin=0.001,xmax){
	ll <- log10(xmin)
	ul <- ceiling(log10(xmax))
	bins <- 0
	for(i in (ll:ul)){
		bins <- c(bins,10^i)
	}
	return(bins)
}		

#-> Define transparency function
# This function was created by 
# http://menugget.blogspot.com/2012/04/adding-transparent-image-layer-to-plot.html
add.alpha <- function(COLORS,ALPHA){
   if(missing(ALPHA)) stop('provide a value for alpha between 0 and 1')
   RGB <- col2rgb(COLORS, alpha=TRUE)
   RGB[4,] <- round(RGB[4,]*ALPHA)
   NEW.COLORS <- rgb(RGB[1,],RGB[2,],RGB[3,],RGB[4,],maxColorValue=255)
   return(NEW.COLORS)
}

#-> Clip reference layers to extent
roads_fc <- gIntersection(roads_fc,ra_extent_fc)
rivers_fc <- gIntersection(rivers_fc,ra_extent_fc)

###---> Feasibility

#-> Read in data
setwd(paste0(wd,'/INPUT/Spatial/Constraints')) #===> Change directory
Fs <- list()
for(i in 1:nrow(tspecs)){
	feas <- raster(paste0(tspecs$T_Feas[i]))
	feas[feas==0] <- NA
	Fs[[i]] <- feas*ra_extent
}

#-> Define bins and colors
llabels <- 'Feasible'
cols <- add.alpha('green',0.6) # Not sure why, but necessary

#-> Maps
setwd(paste0(wd,'/OUTPUT/Treatment_plan/Maps')) #===> Change directory
for(i in 1:nrow(tspecs)){
	tiff(paste0(tspecs$FB_Code[i],'_feasibility.tif'),width=ms$Width_pixels,height=ms$Height_pixels,
	     pointsize=ms$Pointsize,compression='lzw',type='windows')
	par(mar=c(0.6,0.6,2.5,0.6))
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	plot(Fs[[i]],col=cols,axes=F,legend=F,box=F,maxpixels=5000000,add=T)
	title(main=paste0(tspecs$Treatment[i],' Feasibility'),cex.main=1.75)
	plot(rivers_fc,col='blue',lwd=1.5,add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(ra_extent_fc,border='black',lwd=2,add=T)
	legend(as.character(ms$Legend_pos),'Feasibile',fill=cols,bty='n')
	plotRef(bb=ra_extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	g <- dev.off()
	trimLegend(paste0(tspecs$FB_Code[i],'_feasibility.tif'),ms$Width_pixels,ms$Height_pixels,0.925)

}

###---> Risk reduction

#-> Read in data
setwd(paste0(wd,'/OUTPUT/Supplementary')) #===> Change directory
RR <- list()
vmax <- vector()
for(i in 1:nrow(tspecs)){
	RR[[i]] <- raster(paste0(tspecs$FB_Code[i],'_RR.tif'))*ra_extent
	RR[[i]] <- RR[[i]]*Fs[[i]] # Limit to feasible
	vmax[i] <- quantile(values(RR[[i]]),0.99,na.rm=T)
}

#-> Define bins and colors
bins <- logBins(xmin=0.0001,xmax=max(vmax))
llabels <- format(bins[-1],scientific=F,big.L=3,big.mark=',',trim=T,drop0trailing=T)
l2 <- rep('',length(llabels)) # Add high and low labels for clarity
l2[1] <- ' - low'
l2[length(l2)] <- ' - high'
llabels <-paste0(llabels,l2)
cols <- colorRampPalette(c('forestgreen','yellow','red'))(length(bins)-1)

#-> Maps
setwd(paste0(wd,'/OUTPUT/Treatment_plan/Maps')) #===> Change directory
for(i in 1:nrow(tspecs)){
	tiff(paste0(tspecs$FB_Code[i],'_risk_reduction.tif'),width=ms$Width_pixels,
	     height=ms$Height_pixels,pointsize=ms$Pointsize,compression='lzw',type='windows')
	par(mar=c(0.6,0.6,2.5,0.6))
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	RR[[i]][RR[[i]]>max(vmax)] <- max(vmax[i]) # Drop top 1% for display purposes
	plot(RR[[i]],breaks=bins,col=cols,alpha=0.6,axes=F,legend=F,box=F,maxpixels=5000000,add=T)
	title(main=paste0(tspecs$Treatment[i],' Risk Reduction'),cex.main=1.75)
	plot(rivers_fc,col='blue',lwd=1.5,add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(ra_extent_fc,border='black',lwd=2,add=T)
	legend(as.character(ms$Legend_pos),llabels,fill=cols,title=expression(bold('Risk Red. (eNVC)')),
		   bty='n')
	plotRef(bb=ra_extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	g <- dev.off()
	trimLegend(paste0(tspecs$FB_Code[i],'_risk_reduction.tif'),ms$Width_pixels,ms$Height_pixels,
	           0.925)
}

###---> Cost

#-> Read in data
setwd(paste0(wd,'/INPUT/Spatial/Constraints')) #===> Change directory
Cs <- list()
vmax <- vector()
for(i in 1:nrow(tspecs)){
	if(suppressWarnings(!is.na(as.numeric(as.character(tspecs$T_Cost[i]))))){
		Cs[[i]] <- ra_extent*as.numeric(as.character(tspecs$T_Cost[i]))
	}else{
		Cs[[i]] <- raster(as.character(tspecs$T_Cost[i]))
	}
	Cs[[i]] <- Cs[[i]]*Fs[[i]]*ra_extent # Limit to feasible and extent
	vmax[i] <- max(values(Cs[[i]]),na.rm=T)
}

#-> Define bins and colors
bins <- c(0,1000,2500,3000,3500,4000,4500,round(max(vmax),-3))
llabels <- format(bins[-1],scientific=F,big.L=3,big.mark=',',drop0trailing=T)
cols <- colorRampPalette(c('forestgreen','yellow','red'))(length(bins)-1)

#-> Maps
setwd(paste0(wd,'/OUTPUT/Treatment_plan/Maps')) #===> Change directory
for(i in 1:nrow(tspecs)){
	tiff(paste0(tspecs$FB_Code[i],'_cost.tif'),width=ms$Width_pixels,height=ms$Height_pixels,
	     pointsize=ms$Pointsize,compression='lzw',type='windows')
	par(mar=c(0.6,0.6,2.5,0.6))
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	plot(Cs[[i]],breaks=bins,col=cols,alpha=0.6,axes=F,legend=F,box=F,maxpixels=5000000,add=T)
	title(main=paste0(tspecs$Treatment[i],' Cost'),cex.main=1.75)
	plot(rivers_fc,col='blue',lwd=1.5,add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(ra_extent_fc,border='black',lwd=2,add=T)
	legend(as.character(ms$Legend_pos),llabels,fill=cols,title=expression(bold('Cost (USD/ac)')),
		   bty='n')
	plotRef(bb=ra_extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	g <- dev.off()
	trimLegend(paste0(tspecs$FB_Code[i],'_cost.tif'),ms$Width_pixels,ms$Height_pixels,0.925)
}

###---> Cost effectiveness

#-> Combine data
CE <- list()
vmax <- vector()
for(i in 1:nrow(tspecs)){
	CE[[i]] <- RR[[i]]/Cs[[i]]
	vmax[i] <- quantile(values(CE[[i]]),0.99,na.rm=T)
}

#-> Define bins and colors
bins <- logBins(xmin=0.00001,xmax=max(vmax))
llabels <- format(bins[-1],scientific=F,big.L=3,big.mark=',',drop0trailing=T)
l2 <- rep('',length(llabels)) # Add high and low labels for clarity
l2[1] <- ' - low'
l2[length(l2)] <- ' - high'
llabels <-paste0(llabels,l2)
cols <- colorRampPalette(c('forestgreen','yellow','red'))(length(bins)-1)

#-> Maps
setwd(paste0(wd,'/OUTPUT/Treatment_plan/Maps')) #===> Change directory
for(i in 1:nrow(tspecs)){
	tiff(paste0(tspecs$FB_Code[i],'_cost_effectiveness.tif'),width=ms$Width_pixels,
	     height=ms$Height_pixels,pointsize=ms$Pointsize,compression='lzw',type='windows')
	par(mar=c(0.6,0.6,2.5,0.6))
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	CE[[i]][CE[[i]]>max(vmax)] <- max(vmax[i]) # Drop top 1% for display purposes
	plot(CE[[i]],breaks=bins,col=cols,alpha=0.6,axes=F,legend=F,box=F,maxpixels=5000000,add=T)
	title(main=paste0(tspecs$Treatment[i],' Cost Effectiveness'),cex.main=1.75)
	plot(rivers_fc,col='blue',lwd=1.5,add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(ra_extent_fc,border='black',lwd=2,add=T)
	legend(as.character(ms$Legend_pos),llabels,fill=cols,
	       title=expression(bold('Cost Effectiveness')),bty='n',cex=0.9)
	plotRef(bb=ra_extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	g <- dev.off()
	trimLegend(paste0(tspecs$FB_Code[i],'_cost_effectiveness.tif'),ms$Width_pixels,
	           ms$Height_pixels,0.925)
}

###---> Treatment maps by budget

#-> Read in total feasibility
setwd(paste0(wd,'/OUTPUT/Supplementary')) #===> Change directory
tfeas <- raster('tfeas.tif')*ra_extent
tfeas[tfeas==0] <- NA

#-> Read in treatment plans
setwd(paste0(wd,'/OUTPUT/Treatment_plan')) #===> Change directory
TPs <- list()
for(i in 1:nrow(budgets)){
	TPs[[i]] <- shapefile(paste0('tplan_',budgets$Budget[i],'.shp'))
}

#-> Maps
setwd(paste0(wd,'/OUTPUT/Treatment_plan/Maps')) #===> Change directory
for(i in 1:nrow(budgets)){
	tiff(paste0('tplan_',budgets$Budget[i],'.tif'),width=ms$Width_pixels,
	     height=ms$Height_pixels,pointsize=ms$Pointsize,compression='lzw',type='windows')
	par(mar=c(0.6,0.6,2.5,0.6))
	plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
	pfeas <- mask(tfeas,TPs[[i]])
	plot(pfeas,col=add.alpha('green',0.6),axes=F,legend=F,box=F,maxpixels=5000000,add=T)
	plot(TPs[[i]],border='black',lwd=2,add=T)
	title(main=paste0('$',format(budgets$Budget[i],big.L=2,big.mark=','),' Treatment Priorities'),
	      cex.main=1.75)
	plot(rivers_fc,col='blue',lwd=1.5,add=T)
	plot(roads_fc,col='grey20',lwd=1.5,add=T)
	plot(ra_extent_fc,border='black',lwd=2,add=T)
	plotRef(bb=ra_extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
	g <- dev.off()
	trimLegend(paste0('tplan_',budgets$Budget[i],'.tif'),ms$Width_pixels,ms$Height_pixels,0.925)
}

###---> Flattened treatment plan

#-> Read in treatment plans
setwd(paste0(wd,'/OUTPUT/Treatment_plan')) #===> Change directory
tpriority <- shapefile('treatment_priority.shp')

#-> Rasterize
rftp <- rasterize(tpriority,tfeas,'PCode')
rftp <- rftp*tfeas # Limit to inclusive definition of feasible
rftp[rftp==0] <- NA

#-> Define colors
cols <- colorRampPalette(c('red','yellow'))(nrow(budgets))

#-> Maps
setwd(paste0(wd,'/OUTPUT/Treatment_plan/Maps')) #===> Change directory
tiff('treatment_priority.tif',width=ms$Width_pixels,height=ms$Height_pixels,pointsize=ms$Pointsize,
     compression='lzw',type='windows')
par(mar=c(0.6,0.6,2.5,0.6))
plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
plot(rftp,col=cols,alpha=0.6,axes=F,legend=F,box=F,maxpixels=5000000,add=T)
title(main='Fuel Treatment Priorities',cex.main=1.75)
plot(rivers_fc,col='blue',lwd=1.5,add=T)
plot(roads_fc,col='grey20',lwd=1.5,add=T)
plot(ra_extent_fc,border='black',lwd=2,add=T)
legend(as.character(ms$Legend_pos),as.character(budgets$Priority),fill=cols,bty='n')
plotRef(bb=ra_extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
g <- dev.off()
trimLegend('treatment_priority.tif',ms$Width_pixels,ms$Height_pixels,0.925)

#-> Save flattened treatment plan to kml
prftp <-  projectRaster(rftp,crs=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'),method='bilinear')
KML(prftp,'ftp.kml',col=cols,maxpixels=5000000,blur=1,overwrite=T) 

###---> Treatment type map
# It would be nice to do something fancier than this, but it is difficult to precisely allocate 
# multiple treatment types within a catchment.

#-> Select highest budget level
TP <- TPs[[length(TPs)]]

#-> Identify dominant treatment type in catchment
cnames <- paste0('T',seq(1,nrow(tspecs),1),'_ac')
ta <- data.frame(TP)[,cnames]
ta$Dtrt <- NA
for(i in 1:nrow(ta)){
	ta$DT_name[i] <- as.character(tspecs$Treatment)[which(ta[i,cnames]==max(ta[i,cnames]))]
}
TP$DT_name <- ta$DT_name

#-> Map
cols <- c('forestgreen','red','orange','sienna','purple','cyan')
crm <- data.frame(Treatment=tspecs$Treatment,Color=cols[1:nrow(tspecs)])
TP <- merge(TP,by.x='DT_name',crm,by.y='Treatment',all.x=T)
tiff('dominant_treatment.tif',width=ms$Width_pixels,height=ms$Height_pixels,pointsize=ms$Pointsize,
     compression='lzw',type='windows')
par(mar=c(0.6,0.6,2.5,0.6))
plot(hillshade,col=SoG,axes=F,box=F,legend=F,asp=1,maxpixels=5000000)
plot(TP,col=add.alpha(TP$Color,0.6),border='grey40',add=T)
title(main='Dominant Treatment Type',cex.main=1.75)
plot(rivers_fc,col='blue',lwd=1.5,add=T)
plot(roads_fc,col='grey20',lwd=1.5,add=T)
plot(ra_extent_fc,border='black',lwd=2,add=T)
plotRef(bb=ra_extent_fc,pos=as.character(ms$ScaleArrow_pos),blength=ms$ScaleBar_mi)
legend(as.character(ms$Legend_pos),as.character(crm$Treatment),fill=as.character(crm$Color),
       title=expression(bold('Treatment')),bty='n')
g <- dev.off()
trimLegend('dominant_treatment.tif',ms$Width_pixels,ms$Height_pixels,0.925)

###########################################END GRAPHICS#############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
