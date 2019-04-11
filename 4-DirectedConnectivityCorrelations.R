#Testing for correlations between biophysical model predictors and migration as estimated by DivMigrate.
#Written by Cynthia Riginos 2017-2019


library(vegan)



#READIN BIOPHYSICAL MIGRATION MATRICES AND DIVMIGRATE OUTPUTS
# tenuis
load(file="CLEANDATA/tenuis/locations.aten.RData")
bpdist1<-read.csv("ORIGDATA/tenuis_bp_connections/steppingstones_directed.csv", header=FALSE)
bpdist2<-read.csv("ORIGDATA/tenuis_bp_connections/max_flow_directed.csv", header=FALSE)
bpdist3<-read.csv("ORIGDATA/tenuis_bp_connections/reliable_path_directed.csv", header=FALSE)
eucdist<-read.csv("ORIGDATA/tenuis_bp_connections/geog_dist.csv", header=FALSE)
divmig<-read.table("CLEANDATA/PairwiseMigration-Aten/divMigrate_aten1928.txt")
aten=TRUE

# millepora
load(file="CLEANDATA/millepora/locations.mill.RData")
bpdist1<-read.csv("ORIGDATA/millepora_bp_connections/steppingstones_directed.csv", header=FALSE)
bpdist2<-read.csv("ORIGDATA/millepora_bp_connections/max_flow_directed.csv", header=FALSE)
bpdist3<-read.csv("ORIGDATA/millepora_bp_connections/reliable_path_directed.csv", header=FALSE)
eucdist<-read.csv("ORIGDATA/millepora_bp_connections/geog_dist.csv", header=FALSE)
divmig<-read.table("CLEANDATA/PairwiseMigration-Aten-Amill/divMigrate_amil918.txt")
aten=FALSE

#TIDY UP AND UNFOLD DIRECTED MATRICES TO MAKE A DATAFRAME 
bplist<-list(bpdist1, bpdist2, bpdist3, divmig)
for (d in 1:4){
  bpdist<-bplist[[d]]
  diag(bpdist) <- NA
  bpdist<-as.data.frame(bpdist)
  bpdist[bpdist=="Inf"] <-NA #equal to max
  if (aten==TRUE) {  #this is due to clustering individuals from proximate reefs for Aten only
  colnames(bpdist)<- c(1:23, 25, 27:38, 40, 42, 43)
  rownames(bpdist)<- c(1:23, 25, 27:38, 40, 42, 43)
  } else {
    colnames(bpdist)<- seq(1,19)
    rownames(bpdist)<- seq(1,19)
  }
  bplist[[d]]<-bpdist
}

bpdist1<-bplist[[1]]
bpdist2<-bplist[[2]]
bpdist3<-bplist[[3]]
divmig<-bplist[[4]]

# #for adding EucDist
eucdist<-as.data.frame(eucdist)
if (aten==TRUE) {  #this is due to clustering individuals from proximate reefs for Aten only
  colnames(eucdist)<- c(1:23, 25, 27:38, 40, 42, 43) 
  rownames(eucdist)<- c(1:23, 25, 27:38, 40, 42, 43)
} else {
  colnames(eucdist)<- seq(1,19)
  rownames(eucdist)<- seq(1,19)
}

#make a dataframe
bpdist.unfolded<-data.frame((matrix(ncol = 11, nrow = (nrow(bpdist))^2)))
names(bpdist.unfolded)<-c("origin", "destination","O_DistToCoastPC","D_DistToCoastPC", "O_ParallelToCoastPC", "D_ParallelToCoastPC", "SSdist", "max_flow", "reliable_path", "EucDist", "DivMig")
#rows are origins, columns are destinations
u=0  #unfolded matrix row number counter
for (c in 1:ncol(bpdist)) {
  for (r in 1:nrow(bpdist)) {
    u=u+1  #unfolded matrix row number
    bpdist.unfolded$origin[u]<-row.names(bpdist[r,])
    bpdist.unfolded$destination[u]<-colnames(bpdist[c])
    bpdist.unfolded$O_DistToCoastPC[u]<-locations[r,"DistToCoastPC"]
    bpdist.unfolded$D_DistToCoastPC[u]<-locations[c,"DistToCoastPC"]
    bpdist.unfolded$O_ParallelToCoastPC[u]<-locations[r,"ParallelToCoastPC"]
    bpdist.unfolded$D_ParallelToCoastPC[u]<-locations[c,"ParallelToCoastPC"]
    bpdist.unfolded$SSdist[u]<-bpdist1[r,c]
    bpdist.unfolded$max_flow[u]<-bpdist2[r,c]
    bpdist.unfolded$reliable_path[u]<-bpdist3[r,c]
    bpdist.unfolded$EucDist[u]<-eucdist[r,c]
    bpdist.unfolded$DivMig[u]<-divmig[r,c]
  }
} 

# Parse distances relative to coastline
bpdist.unfolded$upcoast<-bpdist.unfolded$D_ParallelToCoastPC -bpdist.unfolded$O_ParallelToCoastPC  #pos = NW, neg = SE
bpdist.unfolded$offshore<-bpdist.unfolded$D_DistToCoastPC -bpdist.unfolded$O_DistToCoastPC  #pos = offshore, neg = inshore

# Transformations for reliable path and max flow
bpdist.unfolded$max_flow_logdist<-1-log(bpdist.unfolded$max_flow)
bpdist.unfolded$reliable_path_logdist<-1-log(bpdist.unfolded$reliable_path)


#TEST FOR CORRELATIONS BETWEEN BP PREDICTORS AND DIVMIGRATE
directed<-c("SSdist", "max_flow_logdist", "reliable_path_logdist")

#bp and divmigrate
if (aten==TRUE) {
  sink("CLEANDATA/tenuis/directed_dist_divmigr_corr_tenuis.txt")
} else {
  sink("CLEANDATA/millepora/directed_dist_divmigr_corr_mill.txt")
}
for (d in 1:length(directed)) {
  distances<-as.data.frame(cbind(bpdist.unfolded[,paste(directed[[d]])],bpdist.unfolded$DivMig))
  distances[distances=="-Inf"] <-NA
  distances[distances=="Inf"] <-NA
  row.has.na <- apply(distances, 1, function(x){any(is.na(x))})
  noNA.dist<- distances[!row.has.na,]
  
  plot(noNA.dist[,2], noNA.dist[,1], xlab=paste(directed[[d]]), ylab="Relative migration", pch=20, col="gray", cex=1.5)
  cat(paste(directed[[d]]), "\r")
  cat(cor(noNA.dist[,2], noNA.dist[,1]), "\r")
  
  #Permutation test 
  obs<-cor(noNA.dist[,2], noNA.dist[,1])
  permuted.dist<-vector()
  for (i in 1:999)  {
    p<-sample(noNA.dist[,1], length(noNA.dist[,1]), replace = FALSE)
    permuted.dist[i]<-cor(p, noNA.dist[,2])
  }
  cat(quantile(permuted.dist, probs = c(0.001, 0.01, 0.05, 0.95, 0.99, 0.999)))
  cat("\r")
}
sink()


#Save files
if(isTRUE(aten)){
  write.csv(bpdist.unfolded, "CLEANDATA/tenuis/bpdist.unfolded.tenuis39.csv")
} else {
  write.csv(bpdist.unfolded, "CLEANDATA/millepora/bpdist.unfolded.mill19.csv")
}


#optional - corelations between metrics
# distances<-as.data.frame(cbind(bpdist.unfolded$SSdist, bpdist.unfolded$max_flow_logdist))
# distances<-as.data.frame(cbind(bpdist.unfolded$SSdist, bpdist.unfolded$reliable_path_logdist))
# distances<-as.data.frame(cbind(bpdist.unfolded$max_flow_logdist, bpdist.unfolded$reliable_path_logdist))
# distances[distances=="Inf"] <-NA
# distances[distances==0] <-NA
# row.has.na <- apply(distances, 1, function(x){any(is.na(x))})
# noNA.dist<- distances[!row.has.na,]
# cor(noNA.dist[,2], noNA.dist[,1])
# permuted.dist<-vector()
# #for (i in 1:999)  {
# p<-sample(noNA.dist[,1], length(noNA.dist[,1]), replace = FALSE)
# permuted.dist[i]<-cor(p, noNA.dist[,2])
# }
# #quantile(permuted.dist, probs = c(0.001, 0.01, 0.05, 0.95, 0.99, 0.999))
# #plot(noNA.dist[,2], noNA.dist[,1])

rm(list=ls())

