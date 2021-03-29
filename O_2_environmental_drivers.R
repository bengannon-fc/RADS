####################################################################################################
#### Environmental drivers of fire behavior and risk
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 07/30/2019
#### Last Modified: 07/24/2020
####################################################################################################
# This will map the fire simulation products including burn probability, flame lengths by
# scenario, and weighted flame length across scenarios.
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
cat('Visualize fire simulation products\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

###########################################START SET UP#############################################

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

#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
rextent <- raster('ma_extent.tif')

#-> Load in conditional fire intensity rasters
setwd(paste0(wd,'/INPUT/SPATIAL/Fire_simulation')) #===> Change directory
BP <- raster('BURNPROBABILITY.tif')
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

#-> Load in composite risk map
setwd(paste0(wd,'/OUTPUT/eNVC/GIS')) #===> Change directory
eNVC <- raster('Total_eNVC.tif')

#-> Load in environmental layers
setwd(paste0(wd,'/INPUT/Spatial/RAW_LANDFIRE')) #==> Change directory
evt <- raster('evt',RAT=T)*rextent
evt_rat <- read.csv('evt_rat.csv',header=T)
dem <- raster('us_dem2010')*rextent
slp <- raster('us_slp2010')*rextent

############################################END SET UP##############################################

##########################################START ANALYSIS############################################
setwd(paste0(wd,'/OUTPUT')) #===> Change directory

###---> Fire behavior by existing vegetation type

#-> Identify top vegetation types
nevt <- freq(evt)
nevt <- data.frame(EVT=nevt[,'value'],Count=nevt[,'count'])
nevt <- nevt[order(-nevt$Count),]
nevt <- nevt[!is.na(nevt$EVT),]
nevt$CCount <- cumsum(nevt$Count)
nevt$CPer <- 100*(nevt$CCount/max(nevt$CCount))
nevt$Per <- 100*(nevt$Count/max(nevt$CCount))
rat <- nevt[1:20,]
rat <- merge(rat,by.x='EVT',evt_rat[,c('VALUE','EVT_Name')],by.y='VALUE',all.x=T)

#-> Summarize flame lengths and CFA by scenario
tiff('FB_by_EVT.tif',width=1200,height=1600,compression='lzw',type='windows')
plot(NULL,xlim=c(0,1),ylim=c(0,1),ylab='',xlab='',axes=F)
scenarios <- c('Low - 25th Percentile','Moderate = 50th Percentile','High = 90th Percentile',
               'Extreme = 97th Percentile')
for(i in 1:4){
	
	#-> Flame lengths
	zt <- data.frame(zonal(FL[[i]],evt))
	srat <- merge(rat,by.x='EVT',zt,by.y='zone',all.x=T)
	
	#-> CFA
	ct <- crosstab(evt,CFA[[i]])
	ct <- data.frame(EVT=as.numeric(rownames(ct)),Unburn=ct[,'0'],Surface=ct[,'1'],
					 Passive=ct[,'2'],Active=ct[,'3'])
	ct$RTot <- rowSums(ct[,-1])
	ct$P0 <- 100*(ct$Unburn/ct$RTot)
	ct$P1 <- 100*(ct$Surface/ct$RTot)
	ct$P2 <- 100*(ct$Passive/ct$RTot)
	ct$P3 <- 100*(ct$Active/ct$RTot)
	ct <- ct[,c('EVT','P0','P1','P2','P3')]
	srat <- merge(srat,by.x='EVT',ct,by.y='EVT',all.x=T)
	
	#-> Reorder
	srat <- srat[order(srat$Per),]
	
	#-> Plot mean flame lengths
	par(mar=c(5.1,32.1,4.1,2.1),fig=c(0,0.6,(1-i*0.25),(1-(i-1)*0.25)),new=T)
	barplot(srat$mean,horiz=T,names=srat$EVT_Name,main='Mean Flame Length',
	        xlab='Flame Length (ft)',las=1,xlim=c(0,70))
	
	#-> Plot CFA barplots
	par(mar=c(5.1,2.1,4.1,2.1),fig=c(0.6,1.0,(1-i*0.25),(1-(i-1)*0.25)),new=T)
	m <- t(as.matrix(srat[,c('P0','P1','P2','P3')]))
	barplot(m,col=c('green','yellow','orange','red'),horiz=T,main='Crown Fire Activity',
	        xlab='Percent',axisnames=F)
	
	#-> Add scenario title
	par(mar=c(5.1,4.1,4.1,2.1),fig=c(0,0.3,(1-i*0.25),(1-(i-1)*0.25)),new=T)
	plot(NULL,xlim=c(0,1),ylim=c(0,1),ylab='',xlab='',axes=F,main=scenarios[i],cex.main=2)
	
}
g <- dev.off()

###---> Expected area burned by existing vegetation type

#-> Convert BP to expected area burned in acres
EAB <- BP*(prod(res(BP))/10000)*2.47105

#-> Summarize by evt
zt <- data.frame(zonal(EAB,evt,'sum'))
srat <- merge(rat,by.x='EVT',zt,by.y='zone',all.x=T)

#-> Reorder
srat <- srat[order(srat$Per),]
	
#-> Plot EAB
tiff('BP_by_EVT.tif',width=960,height=480,compression='lzw',type='windows',pointsize=14)
par(mar=c(5.1,32.1,4.1,2.1))
barplot(srat$sum,horiz=T,names=srat$EVT_Name,main='Annual Expected Area Burned',
	        xlab='Annual Expected Area Burned (ac)',las=1)	
g <- dev.off()

###---> Risk by existing vegetation type

#-> Summarize by evt
zt <- data.frame(zonal(eNVC,evt,'sum'))
srat <- merge(rat,by.x='EVT',zt,by.y='zone',all.x=T)

#-> Reorder
srat <- srat[order(srat$Per),]
	
#-> Plot eNVC
tiff('Risk_by_EVT.tif',width=960,height=480,compression='lzw',type='windows',pointsize=14)
par(mar=c(5.1,32.1,4.1,2.1))
barplot(srat$sum,horiz=T,names=srat$EVT_Name,main='Expected Net Value Change (eNVC)',
	        xlab='Expected Net Value Change (eNVC)',las=1)	
g <- dev.off()

###---> Fire behavior by elevation

#-> Classify elevation
dem <- dem*3.28084
ewidth <- 500
lelev <- floor(cellStats(dem,'min')/ewidth)*ewidth
uelev <- ceiling(cellStats(dem,'max')/ewidth)*ewidth
bins <- seq(lelev,uelev,ewidth)
rcl <- data.frame(From=bins[1:(length(bins)-1)],To=bins[-1],Code=seq(1,(length(bins)-1),1))
ebins <- reclassify(dem,rcl)
edf <- rcl
edf$Label <- paste(format(edf$From,big.L=3,big.mark=','),'to',
                   format(edf$To,big.L=3,big.mark=','),'ft')

#-> Summarize flame lengths and CFA by scenario
tiff('FB_by_elev.tif',width=1200,height=1600,compression='lzw',type='windows')
plot(NULL,xlim=c(0,1),ylim=c(0,1),ylab='',xlab='',axes=F)
scenarios <- c('Low - 25th Percentile','Moderate = 50th Percentile','High = 90th Percentile',
               'Extreme = 97th Percentile')
for(i in 1:4){
	
	#-> Flame lengths
	zt <- data.frame(zonal(FL[[i]],ebins))
	sedf <- merge(edf,by.x='Code',zt,by.y='zone',all.x=T)
	
	#-> CFA
	ct <- crosstab(ebins,CFA[[i]])
	ct <- data.frame(ebin=as.numeric(rownames(ct)),Unburn=ct[,'0'],Surface=ct[,'1'],
					 Passive=ct[,'2'],Active=ct[,'3'])
	ct$RTot <- rowSums(ct[,-1])
	ct$P0 <- 100*(ct$Unburn/ct$RTot)
	ct$P1 <- 100*(ct$Surface/ct$RTot)
	ct$P2 <- 100*(ct$Passive/ct$RTot)
	ct$P3 <- 100*(ct$Active/ct$RTot)
	ct <- ct[,c('ebin','P0','P1','P2','P3')]
	sedf <- merge(sedf,by.x='Code',ct,by.y='ebin',all.x=T)
	
	#-> Reorder
	sedf <- sedf[order(sedf$From),]
	
	#-> Plot mean flame lengths
	par(mar=c(5.1,12.1,4.1,2.1),fig=c(0,0.6,(1-i*0.25),(1-(i-1)*0.25)),new=T)
	barplot(sedf$mean,horiz=T,names=sedf$Label,main='Mean Flame Length',
	        xlab='Flame Length (ft)',las=1,xlim=c(0,70))
	
	#-> Plot CFA barplots
	par(mar=c(5.1,2.1,4.1,2.1),fig=c(0.6,1.0,(1-i*0.25),(1-(i-1)*0.25)),new=T)
	m <- t(as.matrix(sedf[,c('P0','P1','P2','P3')]))
	barplot(m,col=c('green','yellow','orange','red'),horiz=T,main='Crown Fire Activity',
	        xlab='Percent',axisnames=F)
	
	#-> Add scenario title
	par(mar=c(5.1,4.1,4.1,2.1),fig=c(0,0.3,(1-i*0.25),(1-(i-1)*0.25)),new=T)
	plot(NULL,xlim=c(0,1),ylim=c(0,1),ylab='',xlab='',axes=F,main=scenarios[i],cex.main=2)
	
}
g <- dev.off()

###---> Expected area burned by elevation

#-> Summarize by evt
zt <- data.frame(zonal(EAB,ebins,'sum'))
sedf <- merge(edf,by.x='Code',zt,by.y='zone',all.x=T)

#-> Reorder
sedf <- sedf[order(sedf$From),]
	
#-> Plot EAB
tiff('BP_by_elev.tif',width=960,height=480,compression='lzw',type='windows',pointsize=14)
par(mar=c(5.1,12.1,4.1,2.1))
barplot(sedf$sum,horiz=T,names=sedf$Label,main='Annual Expected Area Burned',
	        xlab='Annual Expected Area Burned (ac)',las=1)	
g <- dev.off()

###---> Risk by elevation

#-> Summarize by evt
zt <- data.frame(zonal(eNVC,ebins,'sum'))
sedf <- merge(edf,by.x='Code',zt,by.y='zone',all.x=T)

#-> Reorder
sedf <- sedf[order(sedf$From),]
	
#-> Plot eNVC
tiff('Risk_by_elev.tif',width=960,height=480,compression='lzw',type='windows',pointsize=14)
par(mar=c(5.1,12.1,4.1,2.1))
barplot(sedf$sum,horiz=T,names=sedf$Label,main='Expected Net Value Change (eNVC)',
	        xlab='Expected Net Value Change (eNVC)',las=1)
g <- dev.off()

########################################END ANALYSIS#########################################

#############################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
#######################################END LOGGING###########################################
