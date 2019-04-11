#Making AEM files - second step is to calculate the weights. Here we use stepping stone distances and geographic distances
#Using reefID numbers corresponding to biophysical files
#Takes connections from Site by Edges matrices and looks up weights from unfolded distance matrices (#4)
#Edge weight vectors follow Blanchet, F. G., Legendre, P., Maranger, R., Monti, D., & Pepin, P. (2011). Modelling the effect of directional spatial ecological processes at different scales. Oecologia, 166(2), 357–368. http://doi.org/10.1007/s00442-010-1867-y
#Created by Cynthia Riginos 2016-2018


filepath<-"CLEANDATA/AEM_matrices/" 



#Site by Edges matrices, created by #6a&b
filenames<-c("SitesByEdges_NWtoSE_downstream_AT_rel50", "SitesByEdges_SEtoNW_upstream_AT_rel50","SitesByEdges_downstream_AM_rel50", "SitesByEdges_upstream_AM_rel50")
species<-c("AT", "AT", "AM", "AM")

for(f in 1:4){

site_edgeM<-read.csv(paste(filepath,filenames[f],".csv",sep=""))
comparisons<-colnames(site_edgeM)[-1]

#set up weight matrix
site_weight_matrix<-data.frame(matrix(NA, nrow = length(comparisons), ncol = 7))
site_weight_matrix[,1]<-comparisons
colnames(site_weight_matrix)<-c("edge_name", "From", "To", "GeogDist", "SSDist", "reliable_path", "max_flow")

#read in stepping stone data created in #7
if (species[f]=="AT"){
  bpdist.unfolded<-read.csv("CLEANDATA/tenuis/bpdist.unfolded.tenuis39.csv")
} else {
  bpdist.unfolded<-read.csv("CLEANDATA/millepora/bpdist.unfolded.mill.csv")
}

#pull out node info - origin and destination
for (r in 1:nrow(site_weight_matrix)){
  path<-site_weight_matrix[r,1]
  path<-gsub("X","", path)
  path<-unlist(strsplit(path, "\\."))
  site_weight_matrix[r,2]<-path[1]
  site_weight_matrix[r,3]<-path[2]
}
#limit dataframe to nonzero locations
real<-site_weight_matrix[as.numeric(site_weight_matrix[,"From"]) !=0 ,]
for (r in 1:nrow(real)){
  real[r,"SSDist"]<-bpdist.unfolded[real[r,"From"]==bpdist.unfolded[,"origin"] & real[r,"To"]==bpdist.unfolded[,"destination"], "SSdist"]
  real[r,"max_flow"]<-bpdist.unfolded[real[r,"From"]==bpdist.unfolded[,"origin"] & real[r,"To"]==bpdist.unfolded[,"destination"], "max_flow"]  
  real[r,"reliable_path"]<-bpdist.unfolded[real[r,"From"]==bpdist.unfolded[,"origin"] & real[r,"To"]==bpdist.unfolded[,"destination"], "reliable_path"]  
  real[r,"GeogDist"]<- bpdist.unfolded[real[r,"From"]==bpdist.unfolded[,"origin"] & real[r,"To"]==bpdist.unfolded[,"destination"], "EucDist"] 
}

site_weight_matrix<-real
write.csv(site_weight_matrix, file=(paste(filepath,filenames[f],"_weights.csv",sep="")))

}
