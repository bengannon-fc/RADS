####################################################################################################
#### Spatially optimize treatments to minimize risk
#### Date Created: 07/20/2019
#### Last Modified: 03/15/2021
####################################################################################################
# Summary: this completes both the spatial analyses to parameterize the optimization model and the
# optimization. Optimization is performed with lpSolve using the lpSolveAPI package.
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
log <- file(paste0(wd,'/OUTPUT/Treatment_plan/so_treatments.log'))
sink(file=log,append=T,type=c('output','message'),split=T)
cat('Spatially optimize treatments\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages
.libPaths(pd)
packages <- c('raster','rgdal','rgeos','plyr','lpSolveAPI')
for(package in packages){
	if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
		install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,lib.loc=pd,character.only=T))
	}
}

#-> Load in run specifications files
setwd(paste0(wd,'/INPUT')) #===> Change directory
budgets <- read.csv('so_budgets.csv',header=T)
tspecs <- read.csv('so_treatment_specs.csv',header=T)
MinArea <- read.csv('so_project_size.csv',header=T)[1,'MinArea']
MaxArea <- read.csv('so_project_size.csv',header=T)[1,'MaxArea']

#-> Load in treatment units
setwd(paste0(wd,'/INPUT/Spatial/Treatment_units')) #===> Change directory
TUs <- shapefile('treatment_units.shp')
rTUs <- raster('raster_treatment_units.tif')

#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
ma_extent <- raster('ma_extent.tif')

#-> Load in baseline 
setwd(paste0(wd,'/OUTPUT/eNVC/GIS')) #===> Change directory
eNVC <- raster('Total_eNVC.tif')

#-> Calculate area variables
p2acres <- (prod(res(rTUs))/10000)*2.47105

#-> Define projection
proj <- CRS('+proj=utm +zone=13 +datum=NAD83')

#-> Define Con function to match ArcGIS syntax
Con <- function(condition, trueValue, falseValue){
	return(condition*trueValue + (!condition)*falseValue)
}

#############################################END SET UP#############################################

#########################################START CLEAR OUTPUT#########################################

###---> Empty output folder
e <- try(setwd(paste0(wd,'/OUTPUT/Treatment_plan'))) #===> CD

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Empty output folder
	d_list <- list.files()
	d_list <- d_list[!(d_list %in% 'so_treatments.log')]
	if(length(d_list) > 0){
		unlink(d_list,recursive=T)
	}
}

##########################################END CLEAR OUTPUT##########################################

#######################################START SPATIAL ANALYSIS#######################################

#-> Start data.frame to store treatment unit attributes
tua <- freq(rTUs,useNA='no')
tua[,2] <- tua[,2]*p2acres # Convert to acres
tua <- data.frame(UID=tua[,1],Area_ac=tua[,2])

###---> Feasibility

#-> Load in rasters
setwd(paste0(wd,'/INPUT/Spatial/Constraints')) #==> Change directory
Fs <- list()
for(i in 1:nrow(tspecs)){
	if(tspecs$T_Feas[i]=='all' | tspecs$T_Feas[i]=='All'){
		Fs[[i]] <- Con(!is.na(rTUs),1,0)
	}else{
		Fs[[i]] <- raster(as.character(tspecs$T_Feas[i]))
		crs(Fs[[i]]) <- crs(rTUs)
	}
}

#-> Calculate feasible area (ac) by treatment
for(i in 1:nrow(tspecs)){
	ct <- crosstab(rTUs,Fs[[i]])
	tua$X <- tua$Area_ac*(ct[,'1']/rowSums(ct))
	colnames(tua)[ncol(tua)] <- paste0('T',i,'_feas')
}

#->  Calculate total feasible area (ac) for all treatments
tfeas <- sum(stack(Fs))
tfeas[tfeas>0] <- 1
setwd(paste0(wd,'/OUTPUT/Supplementary')) #===> Change directory
writeRaster(tfeas,'tfeas.tif',format='GTiff',overwrite=T)
ct <- crosstab(rTUs,tfeas)
tua$Tot_feas <- tua$Area_ac*(ct[,'1']/rowSums(ct))

###---> Risk reduction
# Note that this math assumes the input raster is risk per pixel. It is more appropriate to take
# the mean of non-NA pixels if the input raster is a rate (risk per acre) AND the risk reduction
# raster has no NA values other than what we mask out due to feasibility.

#-> Load in rasters
setwd(paste0(wd,'/OUTPUT/PT_eNVC')) #===> Change directory
PTRs <- list()
for(i in 1:nrow(tspecs)){
	PTRs[[i]] <- raster(paste0(tspecs$FB_Code[i],'_eNVC.tif'))
	crs(PTRs[[i]]) <- crs(rTUs)
}

#-> Calculate risk reduction per acre for each treatment
for(i in 1:nrow(tspecs)){
	nfeas <- Fs[[i]] # Create null version of feasibility raster
	nfeas[nfeas==0] <- NA
	rr <- PTRs[[i]] - eNVC
	rr[rr<0] <- 0 # Do not allow risk to increase
	rr <- rr*nfeas
	zt <- zonal(rr,rTUs,'sum',na.rm=T) # Taking sum and dividing by feasible is safer
	tua$X <- zt[,'sum']/tua[,paste0('T',i,'_feas')]
	tua$X[is.na(tua$X)] <- 0 # Set NAs to zero value
	colnames(tua)[ncol(tua)] <- paste0('T',i,'_RR')
}

#-> Calculate total risk while you're at it
zt <- zonal(eNVC,rTUs,'sum',na.rm=T)
tua$eNVC_Tot <- zt[,'sum']
tua$eNVC_mean <- tua$eNVC_Tot/tua$Area_ac

###---> Cost

#-> Load in rasters
setwd(paste0(wd,'/INPUT/Spatial/Constraints')) #===> Change directory
Cs <- list()
for(i in 1:nrow(tspecs)){
	if(suppressWarnings(!is.na(as.numeric(as.character(tspecs$T_Cost[i]))))){
		Cs[[i]] <- Con(!is.na(rTUs),as.numeric(as.character(tspecs$T_Cost[i])),0)
	}else{
		Cs[[i]] <- raster(as.character(tspecs$T_Cost[i]))
		crs(Cs[[i]]) <- crs(rTUs)
	}
}

#-> Calculate cost for each treatment
for(i in 1:nrow(tspecs)){
	nfeas <- Fs[[i]] # Create null version of feasibility raster
	nfeas[nfeas==0] <- NA
	tc <- Cs[[i]]*nfeas
	zt <- zonal(tc,rTUs,'mean',na.rm=T)
	tua$X <- zt[,'mean']
	tua$X[is.na(tua$X)] <- 999999 # Set NAs to extreme cost
	colnames(tua)[ncol(tua)] <- paste0('T',i,'_Cost')
}

#-> Assemble long form data.frame of treatment benefits and constraints for treatment units
# with feasible area to treat
dv.list <- list()
for(i in 1:nrow(tspecs)){
	dv <- tua[,c('UID','Tot_feas',paste0('T',i,'_feas'),paste0('T',i,'_RR'),
	             paste0('T',i,'_Cost'))]
	dv$TrtType <- i
	colnames(dv) <- c('UID','TotFeasAcre','FeasAcre','RedPerAcre','CostPerAcre','TrtType')
	dv.list[[i]] <- dv
}
lfdv <- do.call(rbind,dv.list)
lfdv <- lfdv[lfdv$TotFeasAcre > 0,] # Subset to feasible

#-> Enforce project size requirements
lfdv$FeasAcre <- ifelse(lfdv$FeasAcre < MinArea,0,lfdv$FeasAcre)
lfdv$TotFeasAcre <- ifelse(lfdv$TotFeasAcre > MaxArea,MaxArea,lfdv$TotFeasAcre)

###---> Calculate Risk Reduction and BCR by treatment and save output
setwd(paste0(wd,'/OUTPUT/Supplementary')) #===> Change directory
for(i in 1:nrow(tspecs)){
	rr <- PTRs[[i]] - eNVC
	rr[rr<0] <- 0 # Do not allow risk to increase
	rrsave <- rr*ma_extent # Make version to save
	crs(rrsave) <- proj
	writeRaster(rrsave,paste0(tspecs$FB_Code[i],'_RR.tif'),format='GTiff',overwrite=T)
	bcr <- rr/Cs[[i]]
	bcr[Fs[[i]]==0] <- NA
	bcr <- bcr*ma_extent
	crs(bcr) <- proj
	writeRaster(bcr,paste0(tspecs$FB_Code[i],'_BCR.tif'),format='GTiff',overwrite=T)
}

###---> Calculate BCR and BCR percentiles for treatment units

#-> Benefit-cost ratio
for(i in 1:nrow(tspecs)){
	tua$X <- tua[,paste0('T',i,'_RR')]/tua[,paste0('T',i,'_Cost')]
	colnames(tua)[ncol(tua)] <- paste0('T',i,'_BCR')
}

#-> Percentiles of BCR
eCDF <- function(X){
	X <- X[!is.na(X)]
	X <- X[order(X)]
	rank <- seq(1,length(X),1)
	for(i in 2:length(X)){ # Break ties with lower value
		if(X[i]==X[i-1]){
			rank[i] <- rank[i-1]
		}
	}
	pp <- rank/(length(X)+1)
	return(data.frame(X=X,PP=pp))
}
for(i in 1:nrow(tspecs)){
	edf <- eCDF(tua[,paste0('T',i,'_BCR')])
	colnames(edf) <- c(paste0('T',i,'_BCR'),paste0('T',i,'_Per'))
	edf <- unique(edf)
	tua <- merge(tua,edf,by=paste0('T',i,'_BCR'),all.x=T)
}

#-> Merge to spatial and save
setwd(paste0(wd,'/OUTPUT/Supplementary')) #===> Change directory
TU_stats <- merge(TUs,tua,by='UID',all.x=T)
writeOGR(obj=TU_stats,dsn='.',layer='TrtUnit_metics',driver='ESRI Shapefile',overwrite=T)
write.csv(tua,'TrtUnit_metics.csv',row.names=F)

########################################END SPATIAL ANALYSIS########################################

#########################################START OPTIMIZATION#########################################

###---> Formulate model to maximize risk reduction reduction
# The model constraints include:
# treatment-specific feasibilities
# combined feasibility
# treatment budget proportion constraint 
# total budget constraint

#-> Make an lp object
# arguments are nrow (# of constraints) and ncol (# of decision variables)
rows <- nrow(lfdv)+nrow(lfdv[lfdv$TrtType==1,])+nrow(tspecs)+1
cols <- nrow(lfdv)
lp <- make.lp(nrow=rows,ncol=cols)

#-> Set to maximize
g <- lp.control(lp,sense='max')

#-> Fill in objective function coefficients
set.objfn(lp,lfdv$RedPerAcre)

#-> Fill in constraints for treatment-specific feasibilities
for(i in 1:nrow(lfdv)){
	set.row(lp,i,1,i)
}

#-> Fill in constraints for combined feasibility
br <- nrow(lfdv) # Base row to start index from
tn <- nrow(lfdv[lfdv$TrtType==1,]) # Number of columns per treatment
xt <- rep(1,nrow(tspecs)) # Vector of constraint matrix values
ci <- rep(tn,nrow(tspecs))*seq(0,(nrow(tspecs)-1),1) # Base vector constraint column indices
for(i in 1:tn){
	set.row(lp,br+i,xt,ci+i)
}

#-> Fill in treatment budget proportion constraints
for(i in 1:nrow(tspecs)){
	ti.cpa <- ifelse(lfdv$TrtType==i,lfdv$CostPerAcre,0)
	set.row(lp,(nrow(lfdv)+nrow(lfdv[lfdv$TrtType==1,])+i),ti.cpa)
}

#-> Fill in total budget constraint
set.row(lp,(nrow(lfdv)+nrow(lfdv[lfdv$TrtType==1,])+nrow(tspecs)+1),lfdv$CostPerAcre)

#-> Fill in constraint type
set.constr.type(lp, rep('<=',rows))

#-> Fill in rhs values of constraints
# Final zero is placeholder for program budget
rhsv <- c(lfdv$FeasAcre,lfdv[lfdv$TrtType==1,'TotFeasAcre'],rep(0,nrow(tspecs)),0) 
set.rhs(lp,rhsv)

###---> Solve for provided budgets
setwd(paste0(wd,'/OUTPUT/Treatment_plan')) #===> Change directory

#-> Set up budget summary table
bsum <- data.frame(Budget=budgets$Budget)
bsum$ObjVal <- NA
for(i in 1:nrow(tspecs)){
	bsum$X <- NA
	colnames(bsum)[ncol(bsum)] <- paste0('T',i,'_acres')
	bsum$X <- NA
	colnames(bsum)[ncol(bsum)] <- paste0('T',i,'_USD')
}

#-> loop thru budgets
for(i in 1:nrow(bsum)){
	#-> Solve model
	for(j in 1:nrow(tspecs)){ # Set budget proportions
		rhsv[(nrow(lfdv)+nrow(lfdv[lfdv$TrtType==1,])+j)] <- bsum$Budget[i]*tspecs$MaxBudgetProp[j]
	}
	rhsv[length(rhsv)] <- bsum$Budget[i] # Set total budget
	set.rhs(lp,rhsv)
	solve(lp) # Solve model
	stable <- cbind(lfdv,Acres=get.variables(lp)) # Extract solution
	stable$Acres <- ifelse(stable$Acres > 0 & stable$Acres < MinArea,0,stable$Acres) # Enforce min
	#-> Save output as table
	write.csv(stable,paste0('tplan_',bsum$Budget[i],'.csv'),row.names=F)
	#-> Save output as shapefile
	stt <- data.frame(UID=tua$UID)
	for(j in 1:nrow(tspecs)){
		stt <- merge(stt,stable[stable$TrtType==j,c('UID','Acres')],by='UID',all.x=T)
		colnames(stt)[ncol(stt)] <- paste0('T',j,'_ac')
	}
	stt$Tot_ac <- rowSums(stt[,paste0('T',seq(1,nrow(tspecs)),'_ac')],na.rm=T)
	stplan <- merge(TUs,stt,by='UID',all.x=T) 
	writeOGR(obj=stplan[!is.na(stplan$Tot_ac) & stplan$Tot_ac>0,],dsn='.',
	         layer=paste0('tplan_',bsum$Budget[i]),driver='ESRI Shapefile',overwrite=T)
	#-> Fill in budget summary table
	bsum$ObjVal[i] <- get.objective(lp) # Extract objective value
	for(j in 1:nrow(tspecs)){
		bsum[i,paste0('T',j,'_acres')] <- sum(stable[stable$TrtType==j,'Acres'])
		bsum[i,paste0('T',j,'_USD')] <- sum(stable[stable$TrtType==j,'Acres']*
			stable[stable$TrtType==j,'CostPerAcre'])
	}
}	

#-> Calculate percent of total risk
TRisk <- -sum(tua$eNVC_Tot)
bsum$PerRR <- 100*(bsum$ObjVal/TRisk)

#-> Save budget summary
write.csv(bsum,'budget_summary.csv',row.names=F)

###---> Flatten treatment plans

budgets <- budgets[order(budgets$Budget),]
budgets$PCode <- seq(1,nrow(budgets),1)
for(i in nrow(budgets):1){ # Note reversal of order
	tp <- read.csv(paste0('tplan_',bsum$Budget[i],'.csv'),header=T)
	tp <- ddply(tp,.(UID),summarize,
				Acres = sum(Acres))
	if(exists('ftp')){
		ftp$PCode <- ifelse(tp$Acres>0,i,ftp$PCode)
	}else{
		ftp <- tp
		ftp$PCode <- ifelse(ftp$Acres>0,i,0)
		ftp$Acres <- NULL
	}
}
ftp <- merge(TUs,ftp,by='UID')
ftp <- merge(ftp,budgets,by='PCode')
ftp <- ftp[!is.na(ftp$Priority),]
ftp$Shape_Length <- NULL; ftp$Shape_Area <- NULL
writeOGR(obj=ftp,dsn='.',layer='treatment_priority',driver='ESRI Shapefile',overwrite=T)

##########################################END OPTIMIZATION##########################################

########################################START AVOIDED IMPACT########################################

###---> Avoided Impact Analysis
#-> Estimate maximum budget for avoided impact curve
mtab <- ddply(lfdv,.(TrtType),summarize,
              TotCost = sum(FeasAcre*CostPerAcre,na.rm=T))
maxBudget <- max(mtab$TotCost,na.rm=T)

#-> Determine budget intervals
int <- maxBudget/100
aibudgets <- seq(0,maxBudget,int)

#-> Solve model and save results for each budget level
aia <- data.frame(Budget=aibudgets)
aia$ObjVal <- NA
for(i in 1:nrow(tspecs)){
	aia$X <- NA
	colnames(aia)[ncol(aia)] <- paste0('T',i)
}

for(i in 1:nrow(aia)){
	for(j in 1:nrow(tspecs)){ # Set budget proportions
		rhsv[(nrow(lfdv)+nrow(lfdv[lfdv$TrtType==1,])+j)] <- 
			aia$Budget[i]*tspecs$MaxBudgetProp[j]
	}
	rhsv[length(rhsv)] <- aia$Budget[i] # Set total budget
	set.rhs(lp,rhsv)
	solve(lp) # Solve model
	stable <- cbind(lfdv,Acres=get.variables(lp)) # Extract solution
	aia$ObjVal[i] <- get.objective(lp) # Extract objective value
	for(j in 1:nrow(tspecs)){
		aia[i,paste0('T',j)] <- sum(stable[stable$TrtType==j,'Acres'])
	}
}	
	
trts <- paste0('T',seq(1:nrow(tspecs)))
aia$Total <- rowSums(aia[,trts],na.rm=T)	

#-> Summary figure
tiff('Avoided_impact_curve.tif',width=1000,height=1100,compression='lzw',pointsize=22,
     type='windows')
par(mfrow=c(2,1))

cols <- c('forestgreen','red','orange','sienna','purple','cyan')
xtics <- seq(0,floor(maxBudget/100000000),1)*100000000

#-> Avoided impact curve
ylim <- c(0,max(aia$ObjVal)*1.05)
plot(ObjVal~Budget,data=aia,type='l',xaxs='i',yaxs='i',ylim=ylim,lwd=2,axes=F,
     ylab='Risk Reduction',xlab='Fuel treatment budget ($)',main='Risk Reduction by Budget')
axis(1,at=xtics,labels=paste0('$',xtics/1000000,'M'))
axis(2,at=pretty(aia$ObjVal),labels=format(pretty(aia$ObjVal),big.mark=',',big.interval=3L,
     trim=T,scientific=F))
box()
 
#-> Treatment allocation curves
ylim <- c(0,max(aia$Total,na.rm=T)*1.05)
plot(Total~Budget,data=aia,xaxs='i',yaxs='i',ylim=ylim,lwd=2,col=NA,axes=F,
     ylab='Treated Acres',xlab='Fuel treatment budget ($)',main='Treatment Allocation')
for(i in 1:nrow(tspecs)){
	lines(aia$Budget,aia[,paste0('T',i)],lwd=2,col=cols[i])
}
lines(Total~Budget,data=aia,lwd=2,col='black')	 
axis(1,at=xtics,labels=paste0('$',xtics/1000000,'M'))
axis(2,at=pretty(c(0,max(ylim))),labels=format(pretty(c(0,max(ylim))),big.mark=',',
     big.interval=3L,trim=T,scientific=F))
box()

legend('topleft',c(as.character(tspecs$Treatment),'Total'),
       col=c(cols[1:nrow(tspecs)],'black'),lwd=2,bty='n')
g <- dev.off()
		
#########################################END AVOIDED IMPACT#########################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
sink()
cat('Contents saved to log file.\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################
