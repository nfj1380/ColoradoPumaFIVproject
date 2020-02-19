###################################################
#                                                 #
#   Landscape genomics analysis for mule deer     #
#     Adaptation of Zeller 2017 code by ML        #
#                                                 #
###################################################




#### 2. Transform raw environmental layers to resistance surfaces ####
# Zeller approach: take the raw data and transform it to resistance using a suite of transformations
# (monomolecular, ricker, linear inverse, etc.) to provide a number of different resistance hypotheses

library(raster)
library(rgdal)
library(scales)

### CONTINUOUS VARIABLES ###

## Create resistance surfaces from rasters


#create list of rasters to test from directory
file.list <- list.files("Environ_rasters", pattern=".asc",full.names = T)
file.nms <- list.files("Environ_rasters", pattern=".asc",full.names = F)
file.nms <- gsub("\\..*","",file.nms )

# set output directory
dir.create('TransformedRasters')
out.dir<-"TransformedRasters/"

#Set resistance scale
a <- 1
b <- 100

#Decide what pixel size you want to test
# for MD: smallest that is computationally feasible is 90m (this is default for elevation also), so thinking
# to do 90m, 180m, 450m, 900m (1x, 2x, 5x, 10x)
cell <- c(90,180,450,900,4500) #need this to rename transformed rasters
grain <- c(1,2,5,10,50) #need this to tell the aggregate function how much to increase cell size

for (i in 1:length(file.list)){ #to loop through multiple variables
  for (j in grain) { #run transformations on each pixel size
    rast_gauss<-raster(file.list[i])
    cell <- j*90
    if (j==1) { #way to deal with not changing grain for 90m
      ext<-rast_gauss
      # rescale from 1-100
      rast_gauss<-(((b-a)*(rast_gauss-min(values(rast_gauss), na.rm=T)))/((max(values(rast_gauss), na.rm=T))-(min(values(rast_gauss),na.rm = T))))+a
      
      ### POSITIVE RELATIONSHIPS
      
      ## apply transformation
      ## linear
      resis<-rast_gauss^1
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_positive_linear",sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ## positive concave monomolecular
      resis<--100*(1-exp(0.07*rast_gauss))
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_positive_monomolec_concave",sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ##  positive convex monomolecular
      resis<-100*(1-exp(-0.07*rast_gauss))
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_positive_monomolec_convex",sep=""),format="GTiff",overwrite=TRUE)
      
      ### NEGATIVE RELATIONSHIPS
      
      ## apply transformation
      ## linear
      resis<-(-1*rast_gauss)^1
      resis<-100+resis
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_negative_linear",sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ## negative monomolecular concave
      resis<--100*(1-exp(-0.07*rast_gauss))
      #resis<--1*resis
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_negative_monomolec_concave",sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ## negative monomolecular convex
      resis<-100*(1-exp(0.07*rast_gauss))
      #resis<--1*resis
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_negative_monomolec_convex", sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ## inverse ricker
      resis<--1*rast_gauss*(exp(-0.0275*rast_gauss))
      #resis<--1*resis
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_inverse_ricker", sep=""),format="GTiff",overwrite=TRUE)
    } 
    else {
      # aggregate to alternative grain size
      rast_gauss <- aggregate(rast_gauss,fact=j,fun=mean) 
      ext<-rast_gauss #reset extent since it's different now
      # rescale from 1-100
      rast_gauss<-(((b-a)*(rast_gauss-min(values(rast_gauss), na.rm=T)))/((max(values(rast_gauss), na.rm=T))-(min(values(rast_gauss),na.rm = T))))+a
      
      ### POSITIVE RELATIONSHIPS
      
      ## apply transformation
      ## linear
      resis<-rast_gauss^1
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_positive_linear",sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ## positive concave monomolecular
      resis<--100*(1-exp(0.07*rast_gauss))
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_positive_monomolec_concave",sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ##  positive convex monomolecular
      resis<-100*(1-exp(-0.07*rast_gauss))
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_positive_monomolec_convex",sep=""),format="GTiff",overwrite=TRUE)
      
      ### NEGATIVE RELATIONSHIPS
      
      ## apply transformation
      ## linear
      resis<-(-1*rast_gauss)^1
      resis<-100+resis
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_negative_linear",sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ## negative monomolecular concave
      resis<--100*(1-exp(-0.07*rast_gauss))
      #resis<--1*resis
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_negative_monomolec_concave",sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ## negative monomolecular convex
      resis<-100*(1-exp(0.07*rast_gauss))
      #resis<--1*resis
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_negative_monomolec_convex", sep=""),format="GTiff",overwrite=TRUE)
      
      ## apply transformation
      ## inverse ricker
      resis<--1*rast_gauss*(exp(-0.0275*rast_gauss))
      #resis<--1*resis
      ## rescale from 1 to 100
      resis<-(((b-a)*(resis-min(values(resis), na.rm=T)))/((max(values(resis), na.rm=T))-(min(values(resis),na.rm = T))))+a
      resis<-mask(resis,ext)
      ## write out
      writeRaster(resis, filename=paste(out.dir, file.nms[i], "_", cell, "m", "_inverse_ricker", sep=""),format="GTiff",overwrite=TRUE)
    }
  }
}

#Check what they look like if you want
x <- raster("TransformedRasters/b_ap_r_180m_negative_linear.tif")
plot(x)



### CATEGORICAL VARIABLES ###
#Need a different approach for categorical variables - above re-scaling from 1-100 doesn't work for categorical
# Transform polygons to raster with 0/1 binary coding (1 for feature, 0 for absence of feature) - for
# cateogrical data with more than 1 type (e.g., landcover), need to create separate raster for each 
# category (so binary for forests, sage, developed, etc. separately)



## Highways ##

## Turn polygons/lines into raster
# require(rgeos)
# # Try using national elevation raster and clipping to buffered area around WY (prob need to do this anyway)
# #import state shape
# WD <- "/Users/melanielacava/Desktop/WorkingFolder"
# setwd(WD)
# state <- readOGR(dsn=paste(getwd(),"/LandscapeFeatures",sep=""), layer="state")
# proj4string(state) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# #import national highways data
# nat.highways <- readOGR(dsn=paste(getwd(),"/LandscapeFeatures",sep=""), layer="tl_2017_us_primaryroads", stringsAsFactors=F)
# plot(nat.highways)
# proj4string(nat.highways)
# nat.highways <- spTransform(nat.highways,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# #add buffer around state to crop to (state + 50km buffer)
# temp <- spTransform(state,CRS("+proj=utm +zone=12 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
# buf <- gBuffer(temp,width=50000) #buffer state by 50km
# plot(buf,col="gray")
# plot(temp,add=T)
# buf <- spTransform(buf,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# proj4string(buf)==proj4string(nat.highways)
# #crop to buffered area
# wy.highway <- crop(nat.highways,buf)
# plot(elev)
# plot(wy.highway,add=T)
# #rasterize
# road <- rasterize(wy.highway,elev,field=1,fun="last")
# names(road) <- "interstate"
# plot(road,col="black")
# #export as raster
# writeRaster(road,filename="/Users/melanielacava/Data/MD/LGanalysis/RastersToTest/interstate",format="GTiff",overwrite=TRUE)

#Import rasterized highway layer
WD <- "/Users/melanielacava/Data/MD/LGanalysis/RastersToTest"
setwd(WD)
highways <- raster("interstate_buf50k_res90m.tif")

# Apply Gaussian smooth to surface and transformations to resistance 1-100
library(raster)
library(scales)
#install.packages("smoothie")
library(smoothie)

# read in raw GeoTiff
rast<-highways

# set bandwidth (smoothie uses number of cells for bandwith (esssentially the SD of a Gaussian kernel) 
# so you'll have to calculate distance in number of cells for whatever bandwidth you want to smooth at)
# NOTE: left at what Zeller used - 8
cells<-8

# gaussian smooth with smoothie package
# set up matrix
cellsize<-res(rast)[1]
zmat <- as.matrix(rast)

# run smoothing - takes ~5 min
f <- kernel2dsmooth(zmat, kernel.type="gauss", nx=nrow(rast), ny=ncol(rast),
                    sigma=cells, cellsize=cellsize)

# assign matrix values back to raster
rast_gauss <- rast
values(rast_gauss) <- f

#plot
plot(rast_gauss)

#export
writeRaster(rast_gauss,filename="/Users/melanielacava/Data/MD/LGanalysis/RastersToTest/interstate_buf50k_res90m_smooth",format="GTiff",overwrite=TRUE)

#now you can put the smoothed raster for the categorical variable into the transformation for loop





#STEP 3 
# Move all transformed rasters from laptop to Teton
# Use ResistDist.R and RD_sbatch to run this code on Teton

#### 3. Calculate pairwise resistance distances between all sample pairs ####
# Zeller approach: You can use cost distance (least cost path), commute distance, or Circuitscape 
# to calculate resistance distances. Zeller used cost distance to save time

library(gdistance)


# get raster file paths  
in.dir <- "TransformedRasters"
res.files <- list.files(in.dir, pattern=".tif", full.names = T)

# specify out dir
dir.create('ResistanceDistances')
out.dir <- "ResistanceDistances/"

# read in sample locations from GDM code

duplicated(envWS$Long)
duplicated(envWS$Lat)

dspWS <- SpatialPoints(envWS[,c("Long","Lat")])

sites <-dspWS #from the GDM code
#check that sites overlap raster
plot(x)
plot(sites,add=T)

# loop through each surface -least cost distances
for (i in 1:length(res.files)){
  
  start.time <- Sys.time() # For fun let's time stamp
  
  rast <- raster(res.files[i])
  
  ## Get name
  namegrid <- basename(res.files[i])
  namegrid<-gsub("\\..*","",namegrid)
  
  ## make transition layer
  
  #when you get Error: vector memory exhausted (limit reached?)
  #In terminal ~ create file called .Renviron with code "R_MAX_VSIZE=32Gb" - re-open R studio
  # Check Sys.getenv('R_MAX_VSIZE')
  trCost <- transition(x=1/rast,transitionFunction=mean,directions=8) #Transition object from the raster
  trCost<- geoCorrection(trCost, type="c") #scl=T so that neighbor probabilities normalize to sum to 1 (higher probabilities in transition matrix than without scl=T)
  #plot(raster(trCost)) #to check what transition raster looks like
  
  ## calculate cost distances
  cost_dist <- costDistance(x=trCost,fromCoords=coordinates(sites), toCoords=coordinates(sites))
  #cost_dist[1:9,1:9]
  #curr_dist<- commuteDistance(x=trCost, coords=coordinates(sites)) #crashed my laptop and on Teton
  #try Circuitscape instead of costDistance
  
  ## write out
  write.csv(cost_dist,file=paste(out.dir,"CDMAT_",namegrid,".csv",sep=""))
  
  print("costDistance() complete in ")
  print(Sys.time() - start.time)
  
}

#took 10 mins to run



#### 4. Run MPLE on each resistance surface to determine which transformation is best ####

#GET CODE FROM ERICK TO RUN MLPE CODE

#NOTE: need to compare best model with IBD null model - Zeller 2017 paper says "To determine 
# whether the variables explained the genetic distance among individuals more than Euclidean 
# distance alone, we also ran a regression model with a simple Euclidean distance matrix among 
# sample locations. This resulted in a ΔAICc of278, which was much higher than any other model, 
# indicating the environmental variables ex- plained the genetic distance among individuals 
# better than Euclidean distance alone"
# -- need to figure out how you are able to incorporate this model performance into model
# comparison so we can compare delta AICc values

library(lme4)
library(plyr)
library(MuMIn)
library(ecodist)

#genetic and geographic (Euclidean) distance matrices
genetic <- lower(hostGenWSraw)
euc <- distance(as.matrix(envWS[,2:3]), method = "euclidean")
dim(euc)
### Functions from ResistanceGA ###
# Make to-from population list
To.From.ID<-function(POPS){
  tmp <- matrix(nrow=POPS,ncol=POPS)
  dimnames(tmp) <- list( 1:POPS, 1:POPS)  
  tmp2 <- as.data.frame( which( row(tmp) < col(tmp), arr.ind=TRUE))  
  tmp2[[2]] <-dimnames(tmp)[[2]][tmp2$col]
  tmp2[[1]] <-dimnames(tmp)[[2]][tmp2$row]
  colnames(tmp2)<-c("pop1","pop2")
  as.numeric(tmp2$pop1);as.numeric(tmp2$pop2)
  ID<-arrange(tmp2,as.numeric(pop1),as.numeric(pop2))
  #   ID<-tmp2[with(tmp2, order(pop1, pop2)), ]
  p1<-ID[POPS-1,1]; p2<-ID[POPS-1,2]
  ID[POPS-1,1]<-p2; ID[POPS-1,2]<-p1
  ID$pop1 <- factor(ID$pop1)
  ID$pop2 <- factor(ID$pop2)
  return(ID)
}

# Create ZZ matrix for mixed effects model
ZZ.mat <- function(ID) {
  Zl <- lapply(c("pop1","pop2"), function(nm) Matrix::fac2sparse(ID[[nm]],"d", drop=FALSE))
  ZZ <- Reduce("+", Zl[-1], Zl[[1]])
  return(ZZ)
}


### Read in each cost distance matrix (from geographic distance code) and run the MLPE models
## get file list

in.dir<-"ResistanceDistances/"
file.list<-list.files(in.dir)

## set up df
df<-as.data.frame(matrix(data=NA, nrow=length(file.list), ncol=2))
df[,1]<-file.list

## do MLPE for each in loop
for (i in 1:length(file.list)){
  mm<-read.csv(paste0(in.dir,file.list[i]))
  #mm <- mm[c(-10,-11,-15),c(-1,-11, -12, -16)]
  mm<-mm[,c(-1)]
  m<-nrow(mm)
  #cost<-as.matrix(cost)
  mm<-lower(mm)
  #mm <- mm[which(mm!=-1)]
  
  ID<-To.From.ID(POPS=m)
  
  ZZ<-ZZ.mat(ID=ID)
  
  #cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
  
  #dat<-data.frame(ID,resistance=cs.matrix,response=genetic)
  dat<-data.frame(ID,resistance=mm,response=genetic)
  colnames(dat)<-c("pop1","pop2","resistance","response")
  
  # Fit model
  mod <- lFormula(response ~ resistance + (1|pop1), data=dat, REML=F) #(1|pop1) is using individuals as a random effect
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  MOD <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))   
  df[i,2]<-AICc(MOD)
}

## write out
write.csv(df, "WS_Results.csv",row.names=F)
#df <- read.csv("/Users/melanielacava/Data/MD/LGanalysis/Univariate_MLPE_AICc_Results_Elev180m.csv")

### Figure out best transformation for each variable
df #check for lowest AICc value
df<- df[order(df$V2),] #likelihood


#Need to figure out where to incorprate how many variables are in each model once I combine multiple
# resistance surfances to single surface (AICc accounts for this, but I don't see where it is indicated here)
