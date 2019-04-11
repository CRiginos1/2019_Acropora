#AEM analyses
#Follows Blanchet, F. G., Legendre, P., Maranger, R., Monti, D., & Pepin, P. (2011). Modelling the effect of directional spatial ecological processes at different scales. Oecologia, 166(2), 357–368. http://doi.org/10.1007/s00442-010-1867-y
#Created by Cynthia Riginos 2016-2018


library(AEM)  
library(vegan)
library(adegenet)


#### LOAD IN THE APPROPRIATE EDGE AND WEIGHT MATRICES, THEN PROCEED WITH ANALYSES####
#semimanually created sites by edges matrix and edge weights, see A_SiteByEdgesConstruction_AEM.R, B_EdgeWeightConstruction_AEM.R

filepath<-"CLEANDATA/AEM_matrices/"
filenames<-c("SitesByEdges_NWtoSE_downstream_AT_rel50", "SitesByEdges_SEtoNW_upstream_AT_rel50","SitesByEdges_downstream_AM_rel50", "SitesByEdges_upstream_AM_rel50")
species<-c("tenuis", "tenuis", "millepora", "millepora")
direction<-c("NtoS", "both", "NtoS","both")

for(a in 1:4){
 
#read in genetic data files 
   if (species[a]=="tenuis"){
     load(file="CLEANDATA/tenuis/acroten.genind.RData")
    
     #make species specific edge and weight files
     edge1<-as.matrix(read.csv(paste(filepath,filenames[1],".csv",sep=""))[,-1], header=FALSE)
     weight1<-read.csv(paste(filepath,filenames[1],"_weights.csv",sep=""),header=TRUE)
     
     edge2<-as.matrix(read.csv(paste(filepath,filenames[2],".csv",sep=""))[,-1], header=FALSE)
     weight2<-read.csv(paste(filepath,filenames[2],"_weights.csv",sep=""),header=TRUE)
 
    } else {
    load(file="CLEANDATA/millepora/acromill.genind.RData")
      
     #make species specific edge and weight files
    edge1<-as.matrix(read.csv(paste(filepath,filenames[3],".csv",sep=""))[,-1], header=FALSE)
    weight1<-read.csv(paste(filepath,filenames[3],"_weights.csv",sep=""),header=TRUE)
      
    edge2<-as.matrix(read.csv(paste(filepath,filenames[4],".csv",sep=""))[,-1], header=FALSE)
    weight2<-read.csv(paste(filepath,filenames[4],"_weights.csv",sep=""),header=TRUE)
      
  }

#Create chi-square distances 
acro.pop<-genind2genpop(acro.genind)
acro.chi<-decostand(acro.pop$tab, "chi.square") #Tranform to chi-square genetic distances 

if(direction[a]=="NtoS"){
  edge<-edge1
  weight<-weight1
} else {
  edge<-cbind(edge1, edge2)
  weight<-rbind(weight1, weight2)
}

sink(paste0("CLEANDATA/", species[a], "/AEMresults", direction[a], ".txt"))
cat("**Note that because adj R2 is based on permuations, there will be slight differences in P values among runs!", "\r\r")

#### AEM ANALYSES######

#1- BINARY (all connections weighted equally)
cat("*Binary Delaunay*", "\r\r")

pos.vecs<-floor(nrow(edge)/2)
edge.aem <- aem(binary.mat=edge)   #aem constructs AEM maps, n-1 eigenvectors
edge.aem.vec <- edge.aem$vectors

acro.aem <- rda(acro.chi ~ ., as.data.frame(edge.aem.vec[,1:pos.vecs]))
null<-rda(acro.chi ~ 1, data=as.data.frame(edge.aem.vec[,1:pos.vecs]))
forward.AEM<-ordiR2step(null, scope=formula(acro.aem), directon="forward", psteps=1000, trace=FALSE)  #ordiR2step implements Blanchets stopping criterion

cat(paste(forward.AEM$call[[2]]), "\r")
cat("Adj R2: ", RsquareAdj(forward.AEM)[[2]], "\r\r")


#2 - WEIGHTED (e.g. stepping stone)
cat("*Weighted Delaunay*", "\r\r")
#weightings: SS, max flow, reliable path
edge_weight.SS<-weight[,"SSDist"]
edge_weight.MF<-1-log(weight[,"max_flow"])
edge_weight.R<-1-log(weight[,"reliable_path"])
edge_weights<-list(edge_weight.SS,edge_weight.MF,edge_weight.R)
weightings<-c("Stepping Stones", "Maximum Flow", "Reliable Path")

edgesize<-ncol(edge)

#sink(paste0("CLEANDATA/", species[a], "/AEMresults", direction[a], ".txt"))

for(w in 1:3) {
  cat(weightings[w], "\r\r")
  edge_weight.v<-unlist(edge_weights[w])  
  #transform by weighting scheme
  edge_weight.vec.a<-1-edge_weight.v/max(edge_weight.v) #function a
  edge_weight.vec.b<-1-edge_weight.v/max(edge_weight.v)^2 #function b
  edge_weight.vec.c<-1-edge_weight.v/max(edge_weight.v)^3 #function c
  #edge_weight.vec.d<-1-log(edge_weight.v+1)/max(log(edge_weight.v+1))#function d
  
  edge.weightings<-list(edge_weight.vec.a, edge_weight.vec.b, edge_weight.vec.c)
  for (e in 1:3){
    cat("y = ", e, "\r")
    edge_weight.vec<-edge.weightings[[e]]
    
    #multiply edge matrix by weighting vector
    geog_edge<-data.frame(matrix(ncol = edgesize, nrow = nrow(edge)))
    for (c in 1:edgesize) geog_edge[,c]<-edge_weight.vec[c]*edge[,c]
    geog_edge<-as.matrix(geog_edge)
  
    #svd
    geog_edge.c<-scale(geog_edge, center = TRUE, scale = FALSE)
    geog_edge.svd<-svd(geog_edge.c)
    
    #model selection
    acro.aem <- rda(acro.chi ~ ., as.data.frame(geog_edge.svd$u[,1:pos.vecs]))
    null<-rda(acro.chi ~ 1, data=as.data.frame(geog_edge.svd$u[,1:pos.vecs]))
    forward.AEM<-ordiR2step(null, scope=formula(acro.aem), directon="forward", psteps=1000, trace=FALSE)  #ordiR2step implements Blanchets stopping criterion
  
    cat(paste(forward.AEM$call[[2]]), "\r")
    cat("Adj R2: ", RsquareAdj(forward.AEM)[[2]], "\r\r")
  
  }
}
sink()

}
