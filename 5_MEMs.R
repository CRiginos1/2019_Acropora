#Undertaking MEM based analsyes
#Note that outputs not summarized - record as you go - not really a script
#Largely modified from Borcard, Gillet, & Legendre. 2011. Numerical Ecology with R. Springer.

library(adespatial)
library(ade4)
library(spdep)
library(SoDA)
library(vegan)

#Load in the appropiate set of files
load(file="CLEANDATA/tenuis/acroten.genind.RData")
load(file="CLEANDATA/tenuis/locations.aten.RData")
aten=TRUE


load(file="CLEANDATA/millepora/acromill.genind.RData")
load(file="CLEANDATA/millepora/locations.mill.RData")
aten=FALSE

#PREPARE GEOGRAPHIC AND GENETIC DATA FILES 
data.xy<-as.data.frame(geoXY(locations$Lat, locations$Long, unit=1000))  #unit is KM
data.xy.c <- scale(data.xy, center=TRUE, scale=FALSE) #centered xy values

#convert to chi-square distances 
acro.pop<-genind2genpop(acro.genind)
acro.chi<-decostand(acro.pop$tab, "chi.square") #Tranform to chi-square genetic distances 
anova(rda(acro.chi, data.xy), step=1000)  # Result: significant trend, 5% and 10% variance constrained
response.matrix.det<- resid(lm(as.matrix(acro.chi) ~ ., data=data.xy))  #detrended data - use for MEMs but not AEM



##MEMS
##GOAL/STRATEGY: take various distance matrices and undertake MEM analysis (see Borcard Pg 265)
#1- for each distance model, select MEM variables to retain based on a) lowest AICc, b) forward selection
#2- compare across distance models based on a) lowest AICc, or b) highest adjusted R2


#Distance models
# 1 - Delaunay triangulation connections (neighbors connected), binary:not scaled by distance
# 2 - Delaunay triangulation connections (neighbors connected), scaled by distance raised to power of 1-3
# 3 - PCNM: connections within minimum spanning tree distance and others weighted by 4*Dmin
# 4 - All connections by distance


# MEM analysis of the detrended data - much of the scripting copied from Borcard et al. , adespatial
# ************************************************

if (aten==TRUE) {
  sink("CLEANDATA/tenuis/MEMresults.txt")
} else {
  sink("CLEANDATA/millepora/MEMresults.txt")
}


# 1. Search based on Delaunay triangulation.
#    No weighting matrix (binary weights); function test.W selects among the MEM variables constructed on the basis of the Delaunay triangulation.
cat("*Binary Delaunay*", "\r")

data.del <- tri2nb(data.xy)   #make connections between nearest neighbors, nb object
#plot(data.del, data.xy, col="red", pch=20, cex=1)
dist.del<-nb2listw(data.del, style="B") #creates weights based on neighboring links
scores.del <- scores.listw(dist.del)  #computes MEMs from spatial weighting matrix
del.del.res <- test.W(response.matrix.det, data.del) 
R2.del <- del.del.res$best$AIC$R2[which.min(del.del.res$best$AIC$AICc)] # This line returns the R^2 of the model with the smallest AICc value
cat("Adj R2 best spatial: ", RsquareAdj(R2.del, nrow(data.xy), del.del.res$all$NbVar), "\r\r\r")  #numb models examined, 1 variable in fitted model -data.del.res reports number of variables; numb models is #row of data.del.res$best$values or #AIC tests+1


# 2. Delaunay triangulation weighted by a function of distance.
#    Distances are ranged to maximum 1, and raised to power alpha
cat("*Weighted Delaunay*", "\r\r")

f2 <- function(x, dmax, y) { 1 - (x ^ y) / (dmax) ^ y}  #standard approach
f3 <-function(D, dmax) {1-log(D)/log(dmax)}   #custom log model
max.d1 <- max(unlist(nbdists(data.del, as.matrix(data.xy))))  # Largest Euclidean distance on links belonging to the Delaunay triangulation

cat("a)Power function weightings:", "\r")
del.f2 <- test.W(response.matrix.det, data.del, f=f2, y=1:3, dmax=max.d1, xy=as.matrix(data.xy))  ## Power is set from 1:3
R2.f2.del <- del.f2$best$AIC$R2[which.min(del.f2$best$AIC$AICc)] # This line returns the R^2 of the model with the smallest AICc value
cat("Adj R2 best spatial: ",RsquareAdj(R2.f2.del, nrow(data.xy), del.del.res$all$NbVar), "\r\r")  #numb models examined, 1 variable in fitted model -data.del.res reports number of variables; numb models is #row of data.del.res$best$values or #AIC tests+1

cat("b)Custom log weighting:", "\r")
del.f3 <- test.W(response.matrix.det, data.del, f=f3, dmax=max.d1, xy=as.matrix(data.xy))  
R2.f3.del <- del.f3$best$AIC$R2[which.min(del.f3$best$AIC$AICc)] # This line returns the R^2 of the model with the smallest AICc value
cat("Adj R2 best spatial: ", RsquareAdj(R2.f3.del, nrow(data.xy), del.del.res$all$NbVar),  "\r\r\r")  #numb models examined, 1 variable in fitted model -data.del.res reports number of variables; numb models is #row of data.del.res$best$values or #AIC tests+1


# 3 PCNM:connectivity matrix based on a distance (radius around points)
# Orig Borcard instructions have multiple thresholds created within the larger distance; here I just include one as with original PCNM
cat("*PCNM*", "\r\r")

thresh<-give.thresh(dist(data.xy))  #404.8 km for AT, 287 for AM
list10nb <- lapply(thresh, dnearneigh, x=as.matrix(data.xy), d1=0)
#print(listw2mat(nb2listw(list10nb[[1]], style="B"))[1:10,1:10], digits=1) # Display an excerpt of the first neighbourhood matrix

cat("a)Power function weightings:", "\r")
data.thresh.f2 <- lapply(list10nb, function(x) test.W(x, Y=response.matrix.det, f=f2, 
                                                      y=1:3, dmax=max(unlist(nbdists(x, as.matrix(data.xy)))), 
                                                      xy=as.matrix(data.xy)))
R2.f2.del <- data.thresh.f2[[1]]$best$AIC$R2[which.min(data.thresh.f2[[1]]$best$AIC$AICc)] # This line returns the R^2 of the model with the smallest AICc value
cat("Adj R2 best spatial: ", RsquareAdj(R2.f2.del, nrow(data.xy), del.del.res$all$NbVar),  "\r\r")  

cat("b)Custom log weighting:", "\r")
data.thresh.f3 <- lapply(list10nb, function(x) test.W(x, Y=response.matrix.det, f=f3, 
                                                    dmax=max(unlist(nbdists(x, as.matrix(data.xy)))), 
                                                      xy=as.matrix(data.xy)))
R2.f3.del <- data.thresh.f3[[1]]$best$AIC$R2[which.min(data.thresh.f3[[1]]$best$AIC$AICc)] # This line returns the R^2 of the model with the smallest AICc value
cat("Adj R2 best spatial: ", RsquareAdj(R2.f3.del, nrow(data.xy), del.del.res$all$NbVar),  "\r\r\r")  


## 4 Saturated distance matrix (a: y=1, b: y=2, c: y=3, d: 1-log(d)/max(log(d))
#done long hand with custom matrices
cat("*Custom saturated*", "\r")
cat("(d = 4 is the log function)", "\r\r")

data.threshAll <- dnearneigh(as.matrix(data.xy), 0, 2000)
# plot(data.threshAll, data.xy, col="red", pch=20, cex=0.5, lwd=0.1)

# Conversion of a "nb" object into a "listw" object; "B" is for "binary" - no weighting for distance
data.threshAll.lw <- nb2listw(data.threshAll, style="B")
#print(listw2mat(data.threshAll.lw)[1:10,1:10], digits=1)

# Creation of a spatial weighting matrix W = Hadamard product of B and A
# Replace "1" by Euclidean distances in the connectivity matrix
max.d2 <- max(unlist(nbdists(data.threshAll, as.matrix(data.xy)))) #1621 = longest distance in saturated network
data.threshAll.d2 <- nbdists(data.threshAll, as.matrix(data.xy))

#Distance matrix versions:
data.inv.dist.a <- lapply(data.threshAll.d2, function(x) 1-x/max.d2 )  #f2, y=1
data.inv.dist.b <- lapply(data.threshAll.d2, function(x) 1-(x/max.d2)^2 )  #f2, y=2
data.inv.dist.c <- lapply(data.threshAll.d2, function(x) 1-(x/max.d2)^3 )  #f2, y=3
data.inv.dist.d <- lapply(data.threshAll.d2, function(x) 1-(log(x)/log(max.d2)) )  #f3

dist.weightings<-list(data.inv.dist.a, data.inv.dist.b,data.inv.dist.c, data.inv.dist.d)

for (d in 1:4){
  
  data.inv.dist<-dist.weightings[[d]]
  
  # Creation of spatial weighting matrix W. Argument "B" stands for 
  # "binary" but concerns the links themselves, not their weights
  data.invdist.lw <- nb2listw(data.threshAll, glist=data.inv.dist, style="B")  #adds spatial weights
  #print(listw2mat(data.invdist.lw)[1:10,1:10], digits=2)
  
  # Computation of MEM variables (from an object of class listw)
  data.invdist.MEM <- scores.listw(data.invdist.lw)  #creates Moran's eigenvectors
  #pass listw object to ortho.AIC()
  data.invdist.AIC<-ortho.AIC(response.matrix.det, data.invdist.MEM, ord.var = TRUE)
  
  min(data.invdist.AIC[[1]], na.rm=TRUE)  
  var.numb<-which.min(data.invdist.AIC[[1]])  #ortho.AIC orders the MEMs by their contribution; this returns the numb at which AIC ceases to support adding extra variables
  R2.thresh.cust <- data.invdist.AIC[[4]][which.min(data.invdist.AIC[[1]])]
  cat("\r")
  cat("d = ", paste(d), "\r")
  cat("numb variables: ", paste(var.numb), "\r")
  cat("AICc = ", paste(min(data.invdist.AIC[[1]], na.rm=TRUE)), "\r")
  cat("Adj R2 best spatial: ", RsquareAdj(R2.thresh.cust, nrow(data.xy), var.numb), "\r\r")
}

cat("AICc for the null model: ",data.invdist.AIC$AICc0 , "\r")

sink()


#Tidy up
rm(list=ls())
