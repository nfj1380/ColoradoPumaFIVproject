#-------------------------------------------------------------------#
############GDM modelling#######################
#-------------------------------------------------------------------#


#Authors: Nick Fountain-Jones (nfj@umn.edu) & Simon Dellicour 
#Date May 2019

#prelims
rm(list = ls())


#-------------------------------------------------------------------#
####################Making resistance matrices#######################
#-------------------------------------------------------------------#

library(fields)
library(diagram)
library(vioplot)
library(seraphim)
library(lubridate)
library(tidyr)

#----------------------------------------------------------------------------
#make host genomic data a resistance surface layer for the Western Slope
#----------------------------------------------------------------------------
library(igraph)
#load full host matrix
dpsWS <- read.csv('ws_dps.csv', head=T)

# Make undirected so that graph matrix will be symmetric
gWS <- graph.data.frame(dpsWS , directed=FALSE)

# add value as a weight attribute
gmatWS <- get.adjacency(g, attr="prpshrd", sparse=FALSE)
dim(gmat)
#1

#load coords
hostSpatialWS <- read.csv('ws6.csv')
#check for duplicated values
duplicated(hostSpatialWS$longitude)
duplicated(hostSpatialWS$latitude)

#these individuals were resampled so therefore remove
hostSpatialWS <- hostSpatialWS[-c(30),] 
hostSpatialWS <- hostSpatialWS[-c(77),] 
#add a small ammount of noise to a individual sampled at the same longitude
hostSpatial$latitude[28] <- hostSpatial$longitude[28] +0.0001

dspWS <- SpatialPoints(hostSpatial[,c("longitude","latitude")])

dim(dpsWSmat)

#mkrig to relate spatial data to genomic data
tpsWS <- mKrig(coordinates(dspWS),gmat , find.trA =T,  na.rm=TRUE)

#this is one of the rasters from our data set - should I resent the raster vaules Simon?
tpsInterp <- interpolate(envVariables[[2]], tpsWS)

plot(tpsInterp)
str(hostGenWSraw)

writeRaster(tpsInterp, "TPS_interpolationWS.asc") 

#make host genomic data a resistance surface layer for the FR

#load full host matrix
dpsFR <- read.csv('fr_dps.csv', head=T)

# Make undirected so that graph matrix will be symmetric
gFR <- graph.data.frame(dpsFR , directed=FALSE)

# add value as a weight attribute
gmatFR <- get.adjacency(gFR, attr="prpshrd", sparse=FALSE)
dim(gmatFR)

#load coords
hostSpatialFR <- read.csv('FRpumas5.csv')

#these individuals were resampled so therefore remove
hostSpatialFR <- hostSpatialFR[-c(24),] 
hostSpatialFR <- hostSpatialFR[-c(25),] 

#check for duplicated values
duplicated(hostSpatialFR$longitude)
duplicated(hostSpatialFR$latitude)
#add a small ammount of noise to a individual sampled at the same longitude
hostSpatialFR$longitude[52] <- hostSpatialFR$longitude[52] +0.0001
hostSpatialFR$longitude[47] <- hostSpatialFR$longitude[47] +0.0001
hostSpatialFR$latitude[52] <- hostSpatialFR$latitude[52] +0.0001
hostSpatialFR$latitude[47] <- hostSpatialFR$latitude[47] +0.0001

#collect spatial data
dspFR <- SpatialPoints(hostSpatialFR[,c("longitude","latitude")])

#mkrig to relate spatial data to genomic data
tpsFR <- mKrig(coordinates(dspFR),gmatFR , find.trA =T,  na.rm=TRUE)

#this is one of the rasters from our data set - should I resent the raster vaules Simon?
envVariables[[2]] <- clearValues(envVariables[[2]])
plot(envVariables[[2]])
tpsInterpFR <- interpolate(envVariables[[2]], tpsFR)

plot(tpsInterpFR)
str(hostGenWSraw)

#need to add path to EnvironRasters
writeRaster(tpsInterpFR, "TPS_interpolationFR.asc", overwrite=T)

#different k (these matrices are already)

for (i in 1:2)
{
  rasters = list.files("Environ_rasters")
  rasters = rasters[which(grepl(paste0("b_"),rasters))]
  rasters = gsub(".asc","",rasters)
  envVariables = list()
  resistances = list()
  avgResistances = list()
  c = 0
  for (k in c(10,100,1000))
  { #resistances
    for (j in 1:length(rasters))
    {
      c = c+1
      rast = raster(paste("Environ_rasters/",rasters[j],".asc",sep=""))
      rast[rast[]<0] = 0
      M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
      names(rast) = paste(rasters[j], "_k", k, sep="")
      envVariables[[c]] = rast; names(envVariables[[c]]) = paste(rasters[j],"_k",k,sep="")
      resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
    } #conductances
    for (j in 1:length(rasters))
    {
      c = c+1
      rast = raster(paste("Environ_rasters/", rasters[j],".asc",sep=""))
      rast[rast[]<0] = 0
      M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
      names(rast) = paste(rasters[j], "_k", k, sep="")
      envVariables[[c]] = rast; names(envVariables[[c]]) = paste(rasters[j],"_k",k,sep="")
      resistances[[c]] = FALSE; avgResistances[[c]] = FALSE
    }
  }}
#checked worked
plot(envVariables[[33]]) # 8(r), 24(c), 40(r), 56(c), 72(r), 88(c) FR host relatedness
#write files

envWSnoSite <-envWSnoSite[,c(2,1)] #change long/lat around
envFRnoSite <-envFRnoSite[,c(2,1)]
colnames(envFRnoSite) <- NULL
envWSnoSite <- tibble::rownames_to_column(envWSnoSite, "ID")
colnames(envWSnoSite)<- NULL
write.table(envFR, "FRlongLat_test.txt")
write.table(envWSnoSite, "FRlongLat_test.txt")

writeRaster(envVariables[[33]], filename="ap_100.asc", datatype='ascii', overwrite=TRUE)

#everything from here to 189 is potentially useful but not incorporated
#-------------------------------------------------------
test <- 14
resistances_cs <- list()
for (i, in test)#length(envVariables)))
{circuitScape(i,
              names(paste[i], "cs", sep="_"),
              resistance = TRUE,
              avgResistance = TRUE,
              fourCells = FALSE,
              fromCoor,
              toCoor,
              OS = "Windows",
              prefix = "",
              ID = "",
              nberOfCores_CS=1)
  resistances_cs[[i]]}

test <- circuitScape(envVariables[[6]],
                     names(test1_cs),
                     resistance = TRUE,
                     avgResistance = TRUE,
                     fourCells = FALSE,
                     fromCoor,
                     toCoor,
                     OS = "Windows",
                     prefix = "",
                     ID = "",
                     nberOfCores_CS=1)


#-------------------------------------------------------------------#
####################Constructing GDM #######################
#-------------------------------------------------------------------#

####---------------------FR--------------------------##########
library(gdm)

#-------------------------------------------------------------------#
############Load data #######################
#-------------------------------------------------------------------#

SitesFR <- read.csv("SiteFR.csv", header = TRUE) #in this case they are animal codes

#patristic disatnce

commFRdist <-  read.csv("PatristicDistanceFRMay2019noLabels.csv", header = FALSE)
commFR <- sqrt(commFRdist) #sqrt transform comm matrix improve matrix properties (Vienee et al 2011)
commFRDissim <- (commFR/max(commFR)) #make into a dissimilarity

#have to cloumn bind to sites

commFR <- cbind(SitesFR, commFRDissim)

#predictor matrices (with model 4)

EVIFR <- read.csv("EVIFR.csv", header = FALSE)
EVIFR <- (EVIFR/max(EVIFR))
EVIFR <- cbind(SitesFR, EVIFR)

ccFR <- read.csv("ccFR.csv", header = FALSE)
ccFR <- (ccFR/max(ccFR))
ccFR <- cbind(SitesFR, ccFR)

impFR <- read.csv("impFR.csv", header = FALSE)
impFR <- (impFR/max(impFR))
impFR <- cbind(SitesFR, impFR)

lcFR <- read.csv("lcFR1.csv", header = FALSE)
lcFR <- (lcFR/max(lcFR))
lcFR <- cbind(SitesFR, lcFR)

tempFR <- read.csv("tempFR.csv", header = FALSE)
tempFRdissim <- (tempFR/max(tempFR))
tempFR <- cbind(SitesFR, tempFRdissim)

apFRraw10 <- read.csv("FR_ap10.csv", header = FALSE)
apFRdissim10 <- (apFRraw10/max(apFRraw10))
apFR10 <- cbind(SitesFR, tempFRdissim)

apFRraw100 <- read.csv("FR_ap100.csv", header = FALSE)
apFRdissim100 <- (apFRraw100/max(apFRraw100))
apFR10 <- cbind(SitesFR, tempFRdissim100)

apFRraw10con <- read.csv("FR_ap10con.csv", header = FALSE)
apFRdissim10Con <- (apFRraw10con/max(apFRraw10con))
apFR10con <- cbind(SitesFR, tempFRdissim10con)

roadsFR <- read.csv("roadsFR.csv", header = FALSE)
roadsFR <- (roadsFR/max(roadsFR))
roadsFR <- cbind(SitesFR, roadsFR)

strmsFR <- read.csv("strmsFR.csv", header = FALSE)
strmsFR <- (strmsFR/max(strmsFR))
strmsFR <- cbind(SitesFR, strmsFR)

tpFR <- read.csv("tpFR.csv", header = FALSE)
tpFR <- (tpFR/max(tpFR))
tpFR <- cbind(SitesFR, tpFR)

hostGenFRraw <- read.csv("FRHostGenUpdatedMaywithMiss1.csv", header = FALSE)
hostGenFRDisim <- (hostGenFRraw/max(hostGenFRraw))
hostGenFR <- cbind(SitesFR, hostGenFRDisim)

FRhostGenCSraw <- read.csv('FRhostCS.csv', header= FALSE)
FRhostGenCSdissim <-  (FRhostGenCSraw/max(FRhostGenCSraw))                      
FRhostGen <- cbind(SitesFR, FRhostGenCSdissim)

FRhostGenCSrawCond <- read.csv('FRhostCScond.csv', header= FALSE)
FRhostGenCondSdissim <-  (FRhostGenCSrawCond/max(FRhostGenCSrawCond))                      
FRhostGenCond <- cbind(SitesFR, FRhostGenCondSdissim )
#there are some individuals super connected (really low restance values 0.00001)
FRhostGenCSrawCond1000 <- read.csv('FRhostCScond1000.csv', header= FALSE)
FRhostGenCondSdissim1000 <-  (FRhostGenCSrawCond1000/max(FRhostGenCSrawCond1000))                      
FRhostGenCond1000 <- cbind(SitesFR, FRhostGenCondSdissim1000 )

#Landscape data (predictor vectors)
envFRnoSite <- read.csv("envFR.csv", header = TRUE) 
envFR<- cbind(SitesFR, envFRnoSite)

#-------------------------------------------------------------------#
############Weighting variables by space #######################
#-------------------------------------------------------------------#

#weighted host genetic distance
library(enmSdm)
#remove non-spatial data
envFRnoSite$Year <- NULL
str(envFRnoSite)

#turn coordinates into distances
SpatialDistFR <- as.data.frame(pointDist(envFRnoSite, distFunct = NULL, longLat = c( 'Long', 'Lat')))

#turn into a dissim matrix for Mantel tests
SpatialDistFRdissim <- (SpatialDistFR/max(SpatialDistFR))

#elementwise weighting of genetic data
WeightedHostGenFR <- hostGenFRraw   *  SpatialDistFR

#turn into a dissim
WeightedHostGenFRdissim <- (WeightedHostGenFR/max(WeightedHostGenFR))
WeightedHostGenFR <- cbind(SitesFR, WeightedHostGenFRdissim )


#-------------------------------------------------------------------#
############ Mantel tests to explore correlations#######################
#-------------------------------------------------------------------#
#check for correlations using Mantel without the cbind step above 

library(vegan)
#need to do this without adding site names

mantel(tempFRdissim , SpatialDistFRdissim,  method="pearson", permutations=999) 

#very high correlation in the FR (0.99) no not include in final model
mantel(WeightedHostGenFRdissim , commFRDissim,  method="pearson", permutations=999) 

mantel(apFRdissim100 , apFRdissim10, method="pearson", permutations=999)
mantel(hostGenFRDisim , FRhostGenCondSdissim1000 , method="pearson", permutations=999)
mantel(tempFRdissim , commFRDissim,  method="pearson", permutations=999) 

#-------------------------------------------------------------------#
#-----------Prepare GDM site-pair tables#######################
#-------------------------------------------------------------------#

#bioFormat3 is for already prepared dissim matrices (two step process)

#Step 1 create base GDM object
gdmTabFR <- formatsitepair(commFR, bioFormat=3, XColumn="Lat", YColumn="Long",
                          siteColumn="Site", predData=envFR)

# Step 2 now add the matrices
gdMTabwmatrixFR<- formatsitepair(gdmTabFR, bioFormat=4, siteColumn="Site", predData=envFR, distPreds=list(strmsFR,tpFR, roadsFR, tempFR, lcFR, impFR, ccFR, EVIFR, FRhostGenCond1000 , hostGenFR))

#reduced formulation as per Trumbo et al 2019. Removed temp as it was interaction with geographic distance
gdMTabwmatrixFR<- formatsitepair(gdmTabFR, bioFormat=4, siteColumn="Site", predData=envFR, distPreds=list(strmsFR, roadsFR, impFR, ccFR, FRhostGenCond1000, hostGenFR))

#run the analysis. Splines = null etc means that you accept the default three splines

gdmPumaFR <- gdm(gdMTabwmatrixFR, geo=T) 

summary.gdm(gdmPumaFR)

#basic plot
plot(gdmPumaFR, plot.layout=c(3,4))

#this a more useful plot

plotUncertainty(gdMTabwmatrixFR  ,bsIters=99, sampleSites= 0.9, geo=TRUE, plot.layout=c(3,4))

#backward elimination variable selection

#make sure nperm is atleast 100 for the final model - but try with lower numbers first (say 50). This can take some time. 200 permutations didn't do change the results. Reducing proportion sampled didn't change much - should I keep at 1?

mod.testFR <- gdm.varImp(gdMTabwmatrixFR, geo=T, splines = NULL, knots = NULL, 
           fullModelOnly = FALSE, nPerm = 100, parallel = TRUE, cores = 3, sampleSites = 0.95, sampleSitePairs =1,
           outFile = NULL)

data1a <- as.data.frame(mod.testFR[1]) #model summary
data2a <- as.data.frame(mod.testFR[3]) #model p values
data3a <- as.data.frame(mod.testFR[2]) #variable deviance explained

#save file
save(mod.testFR, file='mod.testFR')
load('mod.testFR')

#plot
par(mfrow = c(1,1)) #resets par
barplot(sort(mod.testFR[[2]][,1], decreasing=T))


####---------------------WS--------------------------##########

#load data sets

SitesWS <- read.csv("SiteWS.csv", header = TRUE)

# Patristic distance 
commWSdist <-  read.csv("PatristicDistanceWSMay2019nolabels.csv", header = FALSE)
#convert to dissim
commWS <- sqrt(commWSdist) #sqrt transform comm matrix improve matrix properties (Vienee et al 2011)
commWSdissim <- (commWS/max(commWS))

#add sites to the com matrix
commWS <- cbind(SitesWS, commWSdissim)

#predictor matrices (with model 4)

EVIWS<- read.csv("EVIWS.csv", header = FALSE)
EVIWS <- (EVIWS/max(EVIWS))
EVIWS<- cbind(SitesWS, EVIWS)

ccWS <- read.csv("ccWS.csv", header = FALSE)
ccWS <- (ccWS/max(ccWS))
ccWS <- cbind(SitesWS, ccWS)

impWS <- read.csv("impWS.csv", header = FALSE)
impWS <- (impWS/max(impWS))
impWS <- cbind(SitesWS, impWS)

lcWS <- read.csv("lcWS.csv", header = FALSE)
lcWS <- (lcWS/max(lcWS))
lcWS <- cbind(SitesWS, lcWS)

tempWS<- read.csv("tempWS.csv", header = FALSE)
tempWS <- (tempWS/max(tempWS))
tempWS<- cbind(SitesWS, tempWS)

roadsWS <- read.csv("roadsWS.csv", header = FALSE)
roadsWS <- (roadsWS/max(roadsWS))
roadsWS <- cbind(SitesWS, roadsWS)

strmsWS <- read.csv("strmsWS.csv", header = FALSE)
strmsWS <- (strmsWS/max(strmsWS))
strmsWS <- cbind(SitesWS, strmsWS)

tpWS<- read.csv("tpWS.csv", header = FALSE)
tpWS<- (tpWS/max(tpWS))
tpWS<- cbind(SitesWS, tpWS)

hostGenWSraw <- read.csv("WSHostGen1.csv", header = FALSE)
hostGenWSdissim <- (hostGenWSraw/max(hostGenWSraw))
hostGenWS <- cbind(SitesWS, hostGenWSdissim)

WShostGenCSraw <- read.csv('WShostCS.csv', header= FALSE)
WShostGenCSdissim <-  (WShostGenCSraw/max(WShostGenCSraw))                      
WShostGen <- cbind(SitesWS, WShostGenCSdissim)
  
#Landscape data (predictor vectors)
envWSnoSite<- read.csv("WSenv.csv", header = TRUE) 
envWS<- cbind(SitesWS, envWSnoSite)

#-------------------------------------------------------------------#
############Weighting variables by space #######################
#-------------------------------------------------------------------#

#weighted host genetic distance
library(enmSdm)

#remove non-spatial data
envWSnoSite$Year <- NULL
str(envWSnoSite)

#turn coordinates into distances
SpatialDistWS <- as.data.frame(pointDist(envWSnoSite, distFunct = NULL, longLat = c( 'Long', 'Lat')))

#turn into a dissim matrix for Mantel tests
SpatialDistWSdissim <- (SpatialDistWS/max(SpatialDistWS))

#elementwise weighting of host genetic data (Haddamard product)
WeightedHostGenWS <- hostGenWSraw   *  SpatialDistWS
 
#turn into a dissim
WeightedHostGenWSdissim <- (WeightedHostGenWS/max(WeightedHostGenWS))
WeightedHostGenWS <- cbind(SitesWS, WeightedHostGenWSdissim )

# high correlation with spatial data (Mantel p = 0.95) - so not appropriate

#-------------------------------------------------------------------#
############ Mantel tests #######################
#-------------------------------------------------------------------#
#check for correlations using Mantel without the cbind step above 

library(vegan)
mantel(WeightedHostGenWSdissim, SpatialDistWSdissim,  method="pearson", permutations=999) 
#also very high correlation coefficent

mantel(SpatialDistWSdissim, WShostGenCSraw,  method="pearson", permutations=999) 

# compare to Mantel

#need to do this without adding site names
mantel(WeightedHostGenWSdissim , commWSdissim, method="pearson", permutations=999)


#-------------------------------------------------------------------#
#-----------Prepare GDM site-pair tables#######################
#-------------------------------------------------------------------#

#prepare site-pair table - bioFormat3 is for already prepared dissim matrices (two step process)
#Step 1 create base GDM object
gdmTabWS <- formatsitepair(commWS, bioFormat=3, XColumn="Lat", YColumn="Long",
                          siteColumn="Site", predData=envWS)
#now add the matrices
gdMTabwmatrixWS <- formatsitepair(gdmTabWS, bioFormat=4, siteColumn="Site", predData=envWS, distPreds=list(tpWS, strmsWS, roadsWS, tempWS, lcWS, impWS, ccWS, EVIWS, hostGenWS, WShostGen))

#reduced formulation as per Trumbo et al 2019
gdMTabwmatrixWS <- formatsitepair(gdmTabWS, bioFormat=4, siteColumn="Site", predData=envWS, distPreds=list(strmsWS, roadsWS, tempWS, ccWS, hostGenWS, WShostGen))


#run the analysis. Splines = null etc means that you accet the default three
gdmPumaWS <- gdm(gdMTabwmatrixWS, geo=T) 
summary.gdm(gdmPumaWS)

#plot
plot(gdmPumaWS, plot.layout=c(3,3))

#plots with uncertainty
plotUncertainty(gdMTabwmatrixWS  ,bsIters=99, sampleSites= 0.9, geo=TRUE, plot.layout=c(3,4))

#variable importance #make sure bs is atleast 99 for the final model - but try with lower numbers first (say 50). This can take some time.

mod.testWS <- gdm.varImp(gdMTabwmatrixWS, geo=T, splines = NULL, knots = NULL, 
                        fullModelOnly = FALSE, nPerm = 100, parallel = TRUE, cores = 2, sampleSites = 0.95, sampleSitePairs =0.95,
                        outFile = NULL)
data4 <- as.data.frame(mod.testWS[1])# model summary
data5 <- as.data.frame(mod.testWS[3]) #model p values
data6 <- as.data.frame(mod.testWS[2]) # variable deviance explained

#save file
save(mod.testWS, file='mod.testWS')
load('mod.testWS')

par(mfrow = c(1,1)) #resets par
barplot(sort(mod.testWS[[2]][,1], decreasing=T))
str(mod.testWS)

#-------------------------------------------------------------------#
#-----------#######################Heat map#######################
#-------------------------------------------------------------------#

#heatmap needs ggplot and reshape2. Heatmap.csv is a file  of the model results (full model) - this is clumsy and culd be improved!
library(ggplot2)
library(reshape2)

dev<-  read.csv("Heatmap.csv", header = T)
dev.m <- melt(dev)
str(dev.m)
p <- ggplot(dev.m, aes(x=variable,y=X, fill=value) 
            
p <- p+geom_tile(aes(colour = "white", size = 0.1)
p <- ggplot(dev.m, aes(x=variable,y=X) ) + geom_tile(aes(fill = value, width = .1), colour = "grey") + scale_fill_gradient(low = "white",high = "red")
base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "", y = "")+ scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))+theme( axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.8,  hjust = 0, colour = "grey50")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank()))

