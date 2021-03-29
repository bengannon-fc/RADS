####################################################################################################
#### Calculate post-treatment expected net value change (eNVC) like Technosylva
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 07/27/2019
#### Last Modified: 03/15/2021
####################################################################################################
# This is a generic workflow designed to calculate post-treatment eNVC for a suite of HVRA 
# from a table of response functions. The workflow is simplified here to only produce the
# final composite eNVC raster for each treatment instead of creating numerous intermediate
# products.
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
cat('Calculate treated expected net value change (eNVC)\n',sep='')
log <- file(paste0(wd,'/OUTPUT/PT_eNVC/TB1.log'))
sink(file=log,append=T,type=c('output','message'),split=T)
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
rasterOptions(maxmemory=10^9)

#-> Load in HVRA and relative importance tables
setwd(paste0(wd,'/INPUT')) #===> Change directory
rs <- read.csv('HVRA_settings.csv',header=T)
rs <- rs[rs$Include==1,]
ri <- read.csv('relative_importance.csv',header=T)

#-> Set path to input rasters folder for any user-provided cNVCs
incNVC <- paste0(wd,'/INPUT/SPATIAL/HVRA_rasters')

#-> Set path to HVRA binary rasters
inrasters <- paste0(wd,'/OUTPUT/HVRA_extents/GIS')

#-> Load in relative area table
setwd(paste0(wd,'/OUTPUT/HVRA_extents')) #===> Change directory
re <- read.csv('relative_area.csv',header=T)

#-> Load in burn probability
setwd(paste0(wd,'/INPUT/SPATIAL/Fire_simulation')) #===> Change directory
bp <- raster('BURNPROBABILITY.tif')

#-> Load in treatment specifications
setwd(paste0(wd,'/INPUT')) #===> Change directory
tspecs <- read.csv('so_treatment_specs.csv',header=T)

#############################################END SET UP#############################################

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/OUTPUT/PT_eNVC'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Empty output folder
	d_list <- list.files()
	if(length(d_list) > 0){
		unlink(d_list[grep('TB1.log',d_list,invert=T)],recursive=T)
	}
}

##########################################END CLEAR OUTPUT##########################################

######################################START RELATIVE IMPORTANCE#####################################
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

#######################################END RELATIVE IMPORTANCE######################################

##########################################START PROCESSING##########################################

for(i in 1:nrow(tspecs)){
	
	#-> Print message
	cat(paste0('#-> Calculating eNVC for ',tspecs$Treatment[i],'\n'))
	
	#-> Load in flame length rasters
	setwd(paste0(wd,'/INPUT/SPATIAL/Fire_simulation')) #===> Change directory
	FL25 <- raster(paste0(tspecs$FB_Code[i],'_25_FLAMELENGTH.tif'))*3.28084 # Convert to feet
	FL50 <- raster(paste0(tspecs$FB_Code[i],'_50_FLAMELENGTH.tif'))*3.28084
	FL90 <- raster(paste0(tspecs$FB_Code[i],'_90_FLAMELENGTH.tif'))*3.28084
	FL97 <- raster(paste0(tspecs$FB_Code[i],'_97_FLAMELENGTH.tif'))*3.28084
	
	###---> Loop thru HVRAs to calculate cNVC
	cNVC.l <- list()
	for(j in 1:nrow(rs)){
	
		#-> Print message
		cat(paste0('Calculating cNVC for ',rs$HVRA[j],'\n'))
		
		#-> Read in raster
		zone <- raster(paste0(inrasters,'/HVRA_',rs$Layer[j],'.tif'))
		
		if(rs$Represents[j]=='Location'){
		
			#-> Prep response function
			FILbins <- c(0,2,4,6,8,12,1000)
			rf <- as.numeric(rs[j,c('FIL1','FIL2','FIL3','FIL4','FIL5','FIL6')])
			rcl <- data.frame(from=FILbins[1:6],to=FILbins[2:7],rf)
			
			#-> Calculate cNVC by scenario
			cNVC25 <- reclassify(FL25,rcl)
			cNVC50 <- reclassify(FL50,rcl)
			cNVC90 <- reclassify(FL90,rcl)
			cNVC97 <- reclassify(FL97,rcl)
			
			#-> Weight scenarios and limit to extent
			cNVC <- zone*(0.01*cNVC25 + 0.09*cNVC50 + 0.2*cNVC90 + 0.7*cNVC97)
			
		}
		
		if(rs$Represents[j]=='cNVC'){
			
			#-> Read in user-provided cNVC raster
			cNVC <- raster(paste0(incNVC,'/',tspecs$FB_Code[i],'_',rs$FeatureClass[j]))
			
			#-> Crop to match extent
			cNVC <- crop(cNVC,zone)
			
		}
		
		#-> Save raster to list
		cNVC.l[[j]] <- cNVC
		
	}
	
	###---> Calculate eNVC

	#-> Clean up burn probability
	bp[is.na(bp)] <- 0

	#-> Get vector of HVRA categories
	cats <- unique(rs$Category)

	#-> Create empty list to store category-level eNVC
	eNVCs <- list()

	#-> Process eNVC by category
	for(j in 1:length(cats)){
		
		cat('Calculating eNVC for ',paste(cats[j]),'\n',sep='')
		
		#-> Subset table
		catrs <- rs[rs$Category==paste(cats[j]),]
		
		#-> Create empty list to store HVRA-level eNVC
		hvra.l <- list()
		
		#-> Calculate weighted eNVC for each HVRA
		for(k in 1:nrow(catrs)){
			
			#-> Read in raster
			cNVC <- cNVC.l[[which(rs$Layer==catrs$Layer[k])]]
			
			#-> Fill null with zero
			cNVC[is.na(cNVC)] <- 0
			
			#-> Calculate eNVC with relative importance adjustment
			# RIoRE = relative importance/relative extent
			hvra.l[[k]] <- suppressWarnings(bp*cNVC)*catrs$RIoRE[k]
		
		}
		
		if(length(hvra.l) > 1){
		
			#-> Stack rasters
			cs <- stack(hvra.l)
			
			#-> Calculate eNVC for category
			eNVC <- sum(cs)
			
		}else{
			
			eNVC <- hvra.l[[1]]
			
		}
		
		#-> Save to category list
		eNVCs[[j]] <- eNVC

	}

	#-> Stack rasters and Calculate integrated eNVC across categories
	if(length(cats) > 1){
		ts <- stack(eNVCs)
		eNVC <- sum(ts)
	}else{
		eNVC <- eNVCs[[1]]
	}
	#-> Save raster
	setwd(paste0(wd,'/OUTPUT/PT_eNVC')) #===> Change directory
	writeRaster(eNVC,paste0(tspecs$FB_Code[i],'_eNVC.tif'),format='GTiff',overwrite=T)

	cat(paste0('Finished ',tspecs$Treatment[i],' eNVC','\n\n'))
	
}

###########################################END PROCESSING###########################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
sink()
cat('Contents saved to log file.\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
