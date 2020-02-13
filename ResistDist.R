library(gdistance)
library(raster)
library(rgdal)

# get raster file paths
in.dir <- "/project/wildgen/mlacava/analysis/muledeer/lg/TransformedRasters"
res.files <- list.files(in.dir, pattern=".tif", full.names = T)

# specify out dir
out.dir <- "/project/wildgen/mlacava/analysis/muledeer/lg/ResistanceDistances/"

# read in sample locations
WD <- "/project/wildgen/mlacava/analysis/muledeer/lg"
setwd(WD)
locs<-read.csv("LatLong_406ind.csv", header=F)
locs<-locs[,2:3]
sites <- as.matrix(locs)
sites <- SpatialPoints(sites)

# loop through each surface
for (i in 1:length(res.files)){
  
  start.time <- Sys.time() # For fun let's time stamp
  
  rast <- raster(res.files[i])
  
  ## Get name
  namegrid <- basename(res.files[i])
  namegrid<-gsub("\\..*","",namegrid)
  
  ## make transition layer
  trCost <- transition(x=1/rast,transitionFunction=mean,directions=8) #Transition object from the raster
  trCost<- geoCorrection(trCost, type="c") #scl=T so that neighbor probabilities normalize to sum to 1 (higher probabilities in transition matrix than without scl=T)
  
  ## calculate cost distances
  cost_dist <- costDistance(x=trCost,fromCoords=coordinates(sites), toCoords=coordinates(sites))
  #curr_dist<- commuteDistance(x=trCost, coords=coordinates(sites)) crashed computer -- don't do!!! for now...
  
  ## write out
  write.csv(cost_dist,file=paste(out.dir,"CDMAT_",namegrid,".csv",sep=""))
  
  print("costDistance() complete in ")
  print(Sys.time() - start.time)
  
}
