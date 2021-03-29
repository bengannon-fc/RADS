####################################################################################################
#### Batch FlamMap basic fire behavior
#### Date Created: 07/01/2019
#### Last Modified: 03/22/2021
####################################################################################################
# Summary: this script automates FlamMap basic fire behavior analyses for multiple scenarios.
# NOTE: it uses the RMRS fire behavior library maintained by Stuart Brittain. This modeling can 
# otherwise be completed with the public distribution of FlamMap. 
####################################################################################################
#-> Get working and packages directory paths from command call
initial.options <- commandArgs(trailingOnly=F)
setwd(dirname(sub('--file=','',initial.options[grep('--file=',initial.options)])))
#setwd('C:/Users/bgannon/Desktop/tRADS/scripts')
wd <- getwd()
pd <- paste(c(unlist(strsplit(wd,'/'))[1:(length(unlist(strsplit(wd,'/')))-1)],
            'Portable_R/Packages'),collapse='/')
fmp <- paste(c(unlist(strsplit(wd,'/'))[1:(length(unlist(strsplit(wd,'/')))-1)],
            'FB_Library/bin/TestFlamMap'),collapse='/')
####################################################################################################

###########################################START MESSAGE############################################
cat('Batch FlamMap basic fire behavior\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load in scenarios table
fms <- read.csv('./INPUT/fire_scenarios.csv',header=T)

#-> Specify fuelscape path
fsp <- paste0(wd,'/INPUT/Spatial/Fuel_scenarios')

#############################################END SET UP#############################################

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/INPUT/SPATIAL/Fire_simulation'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Empty output folder
	d_list <- list.files()
	d_list <- d_list[!(d_list %in% c('BURNPROBABILITY.tif','BURNPROBABILITY.tfw'))]
	if(length(d_list) > 0){
		unlink(d_list,recursive=T)
	}
}

##########################################END CLEAR OUTPUT##########################################

###########################################START FlamMap############################################

for(i in 1:nrow(fms)){
	
	cat('###---> Running FlamMap for scenario',as.character(fms$Scenario[i]),'\n')
	
	###---> FlamMap input file
	
	#-> Subset
	fm <- fms[i,]
	
	#-> Start input file
	fmi <- file(paste0(getwd(),'/FlamMap.input'))
	sink(file=fmi,append=T,type=c('output','message'),split=T)
	
	#-> Header
	cat('FlamMap-Inputs-File-Version-1\n')
	cat('\n')
	
	#-> Specify fuel moistures
	cat('FUEL_MOISTURES_DATA: 1\n')
	cat('0',fm$FM_1hr,fm$FM_10hr,fm$FM_100hr,fm$FM_herb,fm$FM_woody,sep=' ')
	cat('\n')
	
	#-> Foliar moisture
	cat('FOLIAR_MOISTURE_CONTENT: 100\n')
	
	#-> Crown fire method
	cat('CROWN_FIRE_METHOD:',as.character(fm$CROWN_FIRE_METHOD))
	cat('\n')
	
	#-> Number of processors
	cat('NUMBER_PROCESSORS: 1\n')
	
	#-> Wind speed
	cat('WIND_SPEED:',fm$WIND_SPEED)
	cat('\n')
	
	#-> Wind direction
	# Note that -1 is uphill and -2 is downhill
	cat('WIND_DIRECTION:',fm$WIND_DIRECTION)
	cat('\n')
	
	#-> WindNinja settings
	if(fm$GRIDDED_WINDS_GENERATE=='Yes'){
		#-> Turn on WindNinja
		cat('GRIDDED_WINDS_GENERATE: Yes\n')
		#-> Set WindNinja resolution
		cat('GRIDDED_WINDS_RESOLUTION:',fm$GRIDDED_WINDS_RESOLUTION)
		cat('\n')
	}
	
	#-> Outputs
	# Options: FLAMELENGTH, CROWNSTATE, SPREADRATE,INTENSITY,HEATAREA,MIDFLAME,HORIZRATE,
	# MAXSPREADDIR,ELLIPSEDIM_A, ELLIPSEDIM_B, ELLIPSEDIM_C, MAXSPOT_DIR, MAXSPOT_DX, 
	# SOLARRADIATION, FUELMOISTURE1, FUELMOISTURE10, FUELMOISTURE100, FUELMOISTURE1000, 
	# WINDDIRGRID, WINDSPEEDGRID, CROWNFRACTIONBURNED
	outs <- unlist(strsplit(as.character(fm$Outputs),', '))
	for(j in 1:length(outs)){
		cat(outs[j],': 1\n',sep='')
	}

	#-> Close file connection
	sink()
	close(fmi)
	
	###---> FlamMap command file
	
	#-> Start command file
	fmc <- file(paste0(getwd(),'/FMcommand.txt'))
	sink(file=fmc,append=T,type=c('output','message'),split=T)
	
	#-> Specify landscape file, FlamMap input file, output path, output format 
	# (0 = both, 1 = ASCII grid, 2 = FlamMap binary grid [GeoTIFF])
	cat(paste0(fsp,'/',fm$LCP),
	    paste0(getwd(),'/FlamMap.input'),
	    paste0(getwd(),'/',fm$Scenario),2)
	
	#-> Close file connection
	sink()
	close(fmc)
	
	###---> Run scenario in FlamMap
	system2(fmp,args=paste0(getwd(),'/FMcommand.txt'))
	
	#-> Clean up temp files
	unlink(c('FlamMap.input','FMcommand.txt'),recursive=T)
		
	cat('\n\n')

}

############################################END FlamMap#############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
