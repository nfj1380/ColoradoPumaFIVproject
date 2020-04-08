#-------------------------------------------------------------------#
############ Landscape genetics/phylodynamics analysis pipeline modelling#######################
#-------------------------------------------------------------------#


#Authors: Nick Fountain-Jones (nfj@umn.edu), Erick Gagne & Simon Dellicour 
#Date April 2020

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
#Make resistance matrices with different levels of K
#----------------------------------------------------------------------------

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
plot(envVariables[[33]]) # 8(r), 24(c), 40(r), 56(c), 72(r), 88(c) FR host relatedne

#----------------------------------------------------------------------------

#ResiastanceGA routine to optimize resistance/conductance surfaces absed on host genetic distance

#----------------------------------------------------------------------------

library(devtools)
# devtools::install_github("wpeterman/ResistanceGA", build_vignettes=T)

library(doMC)
library(raster)
library(ResistanceGA)

# note: you cannot have any space " " in getwd() !!
# to do: sudo chmod 755 /usr/local/bin/csrun.py

files = list.files("Original_rasters_1")
for (i in 1:length(files))
{
  if (i == 1)
  {
    template = raster(paste0("Original_rasters_1/",files[i]))
    writeRaster(template,paste0("Original_rasters_2/",files[i]), overwrite=T)
  }	else		{
    rast = raster(paste0("Original_rasters_1/",files[i]))
    if ((dim(rast)[1]==dim(template)[1])&(dim(rast)[2]==dim(template)[2]))
    {
      extent(rast) = extent(template)
    }	else	{
      rast = resample(rast, template)
    }
    writeRaster(rast,paste0("Original_rasters_2/",files[i]), overwrite=T)
  }	
}

datasets = c("FR","WS")
rastersToDiscard = c("b_HostWS.asc","b_HostFR.asc")

# 1. Univariate analyses

for (i in 1:length(datasets))
{
  samples = read.csv(paste0(datasets[i],"_coordinates.csv"), header=T)[,2:1]
  write.table(samples, paste0(datasets[i],"_coordinates.txt"), sep="\t", col.names=F, quote=F)
  samplingPoints = SpatialPoints(samples)
  geneticDistances = read.csv(paste0(datasets[i],"_host_dgen.csv"), header=F)
  geneticDistances = geneticDistances[lower.tri(geneticDistances)]
  CS_Point.File = paste0(datasets[i],"_coordinates.txt")
  CS.program = "'/usr/local/bin/csrun.py'"
  files = list.files("Original_rasters_2")
  files = files[which(files!=rastersToDiscard[i])]
  registerDoMC(cores=8); buffer = list()
  buffer = foreach(j = 1:length(files)) %dopar% {
    # for (j in 1:length(files)) {
    rast = raster(paste0("Original_rasters_2/",files[j]))
    resultsDir = paste0(getwd(),"/ResistanceGA_",datasets[i],"/",names(rast),"/"); dir.create(resultsDir, showWarnings=F)
    GA.inputs = GA.prep(ASCII.dir=rast, Results.dir=resultsDir, method="AIC", select.trans="A")
    CS.inputs = CS.prep(n.Pops=length(samplingPoints), response=geneticDistances, CS_Point.File=CS_Point.File, CS.program=CS.program, platform="other")
    SS_result = SS_optim(CS.inputs=CS.inputs, GA.inputs=GA.inputs); # plot(SS_result)
    j
  }
}

# 2. Multivariate analyses

for (i in 1:length(datasets))
{
  samples = read.csv(paste0(datasets[i],"_coordinates.csv"), header=T)[,2:1]
  write.table(samples, paste0(datasets[i],"_coordinates.txt"), sep="\t", col.names=F, quote=F)
  samplingPoints = SpatialPoints(samples)
  geneticDistances = read.csv(paste0(datasets[i],"_host_dgen.csv"), header=F)
  geneticDistances = geneticDistances[lower.tri(geneticDistances)]
  CS_Point.File = paste0(datasets[i],"_coordinates.txt")
  CS.program = "'/usr/local/bin/csrun.py'"
  files = list.files("Original_rasters_2"); stack = list()
  files = files[which(files!=rastersToDiscard[i])]
  for (j in 1:length(files)) stack[[j]] = raster(paste0("Original_rasters_2/",files[j]))
  stack = stack(stack); resultsDir = paste0(getwd(),"/ResistanceGA_",datasets[i],"/")
  GA.inputs = GA.prep(ASCII.dir=stack, Results.dir=resultsDir, method="AIC")
  CS.inputs = CS.prep(n.Pops=length(samplingPoints), response=geneticDistances, CS_Point.File=CS_Point.File, CS.program=CS.program, platform="other")
  MM_result = SS_optim(CS.inputs=CS.inputs, GA.inputs=GA.inputs); plot(MM_result)
}



#----------------------------------------------------------------------------
#HOST GENETIC RESISTANCE
#make host genomic data a resistance surface layer 
#----------------------------------------------------------------------------
library(igraph)
#load full host matrix
dpsWS <- read.csv('ws_dps.csv', head=T) 

# Make undirected so that graph matrix will be symmetric
gWS <- graph.data.frame(dpsWS , directed=FALSE)

# add value as a weight attribute
gmatWS <- get.adjacency(gWS, attr="prpshrd", sparse=FALSE)
dim(gmatWS)
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
hostSpatialWS$latitude[28] <- hostSpatialWS$latitude[28] +0.0001
hostSpatialWS$longitude[28] <- hostSpatialWS$longitude[28] +0.0001
dspWS <- SpatialPoints(hostSpatialWS[,c("longitude","latitude")])

dim(dpsWS)

#mkrig to relate spatial data to genomic data
tpsWS <- mKrig(coordinates(dspWS),gmatWS , find.trA =T,  na.rm=TRUE)

#this is one of the rasters from our data set - should I resend the raster vaules Simon?
tpsInterp <- interpolate(envVariables[[2]], tpsWS) #envVariables in this case sets the ratser spatial coordinates

plot(tpsInterp)
str(hostGenWSraw)

writeRaster(tpsInterp, "TPS_interpolationWS.asc") # this raster was then taken 
#and run through the circuitscape process externally to R

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

#write files

envWSnoSite <-envWSnoSite[,c(2,1)] #change long/lat around
envFRnoSite <-envFRnoSite[,c(2,1)]
colnames(envFRnoSite) <- NULL
envWSnoSite <- tibble::rownames_to_column(envWSnoSite, "ID")
colnames(envWSnoSite)<- NULL
write.table(envFR, "FRlongLat_test.txt")
write.table(envWSnoSite, "FRlongLat_test.txt")

writeRaster(envVariables[[33]], filename="ap_100.asc", datatype='ascii', overwrite=TRUE)


#-------------------------------------------------------------------#
####################Constructing GDM #######################
#-------------------------------------------------------------------#

####---------------------WUI--------------------------##########
library(gdm)

#-------------------------------------------------------------------#
############Load data #######################
#-------------------------------------------------------------------#

SitesFR <- read.csv("SiteFR.csv", header = TRUE) #in this case they are animal codes

#patristic disatnce

commFRdist <-  read.csv("PatristicDistanceFRMay2019noLabels.csv", header = FALSE)
commFR <- sqrt(commFRdist) #sqrt transform comm matrix improve matrix properties (Vienee et al 2011)
commFRDissim <- (commFR/max(commFR)) #make into a dissimilarity

#write.csv(commFRDissim, "commFRDissim.csv" ) #make this a csv file for MLPE analysis
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
write.csv(FRhostGenCSrawCond, "FRhostGenCSrawCond.csv")
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
############ Mantel tests to explore correlations#######################
#-------------------------------------------------------------------#
#check for correlations using Mantel without the cbind step above 

library(vegan)
#need to do this without adding site names

mantel(tempFRdissim , SpatialDistFRdissim,  method="pearson", permutations=999) 

#very high correlation in the FR (0.99) no not include in final model
mantel(WShostGenCSraw, roadsWS  ,  method="pearson", permutations=999) 

mantel(apFRdissim100 , apFRdissim10, method="pearson", permutations=999)
mantel(WShostGenC , WShostGen , method="pearson", permutations=999)
mantel(WShostGenC, commWSdist,  method="pearson", permutations=999) 

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


####---------------------UB--------------------------##########

#load data sets

SitesWS <- read.csv("SiteWS.csv", header = TRUE)

# Patristic distance 
commWSdist <-  read.csv("PatristicDistanceWSMay2019nolabels.csv", header = FALSE)
#convert to dissim
commWS <- sqrt(commWSdist) #sqrt transform comm matrix improve matrix properties (Vienee et al 2011)
commWSdissim <- (commWS/max(commWS))
write.csv(commWSdissim , "commWSdissim.csv ")
#add sites to the com matrix
commWS <- cbind(SitesWS, commWSdissim)

#predictor matrices (with model 4)

EVIWSraw<- read.csv("EVIWS.csv", header = FALSE)
EVIWS <- (EVIWSraw/max(EVIWSraw))
EVIWS<- cbind(SitesWS, EVIWS)

ccWSraw <- read.csv("ccWS.csv", header = FALSE)
ccWS <- (ccWSraw/max(ccWSraw))
ccWS <- cbind(SitesWS, ccWS)

impWSraw <- read.csv("impWS.csv", header = FALSE)
impWS <- (impWSraw/max(impWSraw))
impWS <- cbind(SitesWS, impWS)

lcWSraw <- read.csv("lcWS.csv", header = FALSE)
lcWS <- (lcWSraw/max(lcWSraw))
lcWS <- cbind(SitesWS, lcWS)

tempWSraw<- read.csv("tempWS.csv", header = FALSE)
tempWS <- (tempWSraw/max(tempWSraw))
tempWS<- cbind(SitesWS, tempWS)

roadsWS <- read.csv("roadsWS.csv", header = FALSE)
roadsWS <- (roadsWS/max(roadsWS))
roadsWS <- cbind(SitesWS, roadsWS)

strmsWSraw <- read.csv("strmsWS.csv", header = FALSE)
strmsWS <- (strmsWSraw/max(strmsWSraw))
strmsWS <- cbind(SitesWS, strmsWS)

tpWSraw<- read.csv("tpWS.csv", header = FALSE)
tpWS<- (tpWSraw/max(tpWSraw))
tpWS<- cbind(SitesWS, tpWS)

hostGenWSraw <- read.csv("WSHostGen1.csv", header = FALSE)
hostGenWSdissim <- (hostGenWSraw/max(hostGenWSraw))
hostGenWS <- cbind(SitesWS, hostGenWSdissim)

WShostGenCSraw <- read.csv('WShostCS.csv', header= FALSE)
WShostGenCSdissim <-  (WShostGenCSraw/max(WShostGenCSraw))                      
WShostGenC <- cbind(SitesWS, WShostGenCSdissim)
  
#Landscape data (predictor vectors)
envWSnoSite<- read.csv("WSenv.csv", header = TRUE) 
envWS<- cbind(SitesWS, envWSnoSite)

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

gdMTabwmatrixWS <- formatsitepair(gdmTabWS, bioFormat=4, siteColumn="Site", predData=envWS, distPreds=list(strmsWS))


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
p <- ggplot(dev.m, aes(x=variable,y=X, fill=value)) 
            
p <- p+geom_tile(aes(colour = "white", size = 0.1))
p <- ggplot(dev.m, aes(x=variable,y=X) ) + geom_tile(aes(fill = value, width = .1), colour = "grey") + scale_fill_gradient(low = "white",high = "red")
base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "", y = "")+ scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))+theme( axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.8,  hjust = 0, colour = "grey50")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank()))

#-------------------------------------------------------------------#
#-----------Is FIV prevalence different in each study area? #########
#-------------------------------------------------------------------#

res <- prop.test(x = c(43, 40), n = c(74, 97))
# Printing the results
res 

#-------------------------------------------------------------------#
#-----------MLPE----------------#
# Thanks to Robert Peterman for helping with this code
#-------------------------------------------------------------------#
#to add resistanceGA

## Code to fit MLPE models

library(ResistanceGA)
library(AICcmodavg)
library(MuMIn)
library(lme4)


#---------------------------------
#For the UB
#---------------------------------
csv_files <- list.files("WS_orig_resistance", ".csv", full.names = T)
csv_files #point the analysis to the folder where the resistance surfaces are
commDiss <- read.csv(csv_files[3])[,-1] 

csv_list <- lapply(csv_files[-3], read.csv, header = F) %>%
  lapply(., lower)

names(csv_list) <- basename(csv_files[-3]) %>% sub('.csv', "", .)

id <- To.From.ID(nrow(commDiss))

df <- scale(do.call(cbind, csv_list)) %>% as.data.frame()
df$commDiss <- lower(commDiss)
df$pop <- id$pop1

mlpe_list <- vector('list', length(csv_list))
names(mlpe_list) <- names(csv_list)
for(i in 1:length(mlpe_list)){
  mlpe_list[[i]] <- mlpe_rga(commDiss ~ scale(df[,i]) + (1|pop),
                             data = df,
                             REML = F)
}

mod_sel <- AICcmodavg::bictab(mlpe_list,
                   modnames = names(mlpe_list),
                   nobs = nrow(commDiss))
mod_sel

m1 <- mlpe_rga(commDiss ~ WShostCS  +roadsWS +  (1|pop),
                         data = df,
                         REML = F, na.action = na.fail)

m2 <- mlpe_rga(commDiss ~ WShostCS +  (1|pop),
               data = df,
               REML = F, na.action = na.fail)

m3 <- mlpe_rga(commDiss ~ roadsWS +  (1|pop),
               data = df,
               REML = F, na.action = na.fail)

multi_listWS <- as.list(m1,m2,m3)
dredge(m1,m2,m3)

#multivariate


#---------------------------------
#For the WUI
#---------------------------------
csv_files_FR <- list.files("FR_orig_resistance", ".csv", full.names = T) 

csv_files_FR#point the analysis to the folder where the resistance surfaces are

commDiss_FR <- read.csv(csv_files_FR[2])[,-1] 

csv_list_FR <- lapply(csv_files_FR[-2], read.csv, header = F) %>%
  lapply(., lower)

names(csv_list_FR) <- basename(csv_files_FR[-2]) %>% sub('.csv', "", .)

id_FR <- To.From.ID(nrow(commDiss_FR))

df_FR <- scale(do.call(cbind, csv_list_FR)) %>% as.data.frame()
df_FR$commDiss_FR <- lower(commDiss_FR)
df_FR$pop <- id_FR$pop1

mlpe_list_FR <- vector('list', length(csv_list_FR))
names(mlpe_list_FR) <- names(csv_list_FR)
for(i in 1:length(mlpe_list_FR)){
  mlpe_list_FR[[i]] <- mlpe_rga(commDiss_FR ~ scale(df_FR[,i]) + (1|pop),
                             data = df_FR,
                             REML = F)
}

mod_sel_FR <- AICcmodavg::bictab(mlpe_list_FR,
                              modnames = names(mlpe_list_FR),
                              nobs = nrow(commDiss_FR))
mod_sel_FR 

#all variable from univariate analysis within 2 AIC
multivariate_FR <- mlpe_rga(commDiss_FR ~ roadsFR  +FRHostGen+ tpFR  + SpatialDistFR+ strmsFR + (1|pop),
                         data = df_FR,
                         REML = F, na.action = na.fail)

#multivariate
library(cAIC4)
multi_step_FR <- stepcAIC(multivariate_FR, direction = "both", trace = TRUE, data = df, groupCandidates=c("pop"))

#-------------------------------------------------------------------#
############Weighting variables by space (not used in manuscript) #######################
#-------------------------------------------------------------------#

#weighted host genetic distance
library(enmSdm) #not working right - some conflict with raster

#remove non-spatial data
envWSnoSite$Year <- NULL
str(envWSnoSite)

#turn coordinates into distances
SpatialDistWS <- as.data.frame(pointDist(envWSnoSite, distFunct = NULL, longLat = c( 'Long', 'Lat')))

#write.csv(SpatialDistWS, "SpatialDistWS.csv", row.names=FALSE) #for MLPE analysis
#turn into a dissim matrix for Mantel tests
SpatialDistWSdissim <- (SpatialDistWS/max(SpatialDistWS))

#elementwise weighting of host genetic data (Haddamard product)
WeightedHostGenWS <- hostGenWSraw   *  SpatialDistWS

#turn into a dissim
WeightedHostGenWSdissim <- (WeightedHostGenWS/max(WeightedHostGenWS))
WeightedHostGenWS <- cbind(SitesWS, WeightedHostGenWSdissim )

# high correlation with spatial data (Mantel p = 0.95) - so not appropriate
