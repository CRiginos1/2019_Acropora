#Making AEM files - first construct the edges semi manuallly. In this script, I start with the most downstream locations and iteratively work up the network
#Using reefID numbers corresponding to biophysical files
#Connections are limited to the 50% most reliable connections
#Site by edges matrix construction follows Blanchet, F. G., Legendre, P., Maranger, R., Monti, D., & Pepin, P. (2011). Modelling the effect of directional spatial ecological processes at different scales. Oecologia, 166(2), 357–368. http://doi.org/10.1007/s00442-010-1867-y
#Created by Cynthia Riginos 2016-2018

#walk upstream function
sitebyedges_expansion <- function(dataframe, upstream_number, downstream_vector){
  aggregate_downstream<-data.frame(matrix(NA, nrow = nrow(dataframe), ncol = length(downstream_vector)))
  
  for (c in 1:length(downstream_vector)){
    pattern<-paste("^", downstream_vector[c], ">\\w+", sep="")  #text for later grepl search
    colnames(aggregate_downstream)[c]<-paste(upstream_number,">",downstream_vector[c], sep="")
    downstream_res<-data.frame(matrix(NA, nrow = nrow(dataframe), ncol = 0))
    for (s in 1:ncol(dataframe)) { #find previous columns further downstream from downstream vector element
      if (grepl(pattern, colnames(dataframe)[s])) {downstream_res[,ncol(downstream_res)+1]<-dataframe[,s]}  
    }
    aggregate_downstream[,c]<-rowSums(downstream_res) #combine all the downstream effects
    aggregate_downstream[downstream_vector[c],c]<-1
  }
  dataframe<-cbind(dataframe,aggregate_downstream )
  
  return(dataframe)
} #creates connections from upstream_number to all items in downstream_vector and propogates through dataframe

#######
# 1 NWtoSE all connections for Acropora millepora, reliable>50% for 19 locations (downstream)

#make data frame
SitesByEdges_NWtoSE_AM<-data.frame(matrix(NA, nrow = 19, ncol = 0))

#define columns for fully downstream (end) locations
upstream<-c(5,9,15, 17,18, 19)
SitesByEdges_NWtoSE_AM[,1:length(upstream)]<-0
downstream<-"x"
for (c in 1:length(upstream)) {
  colnames(SitesByEdges_NWtoSE_AM)[c]<-paste(upstream[c], ">", downstream, sep="")
  SitesByEdges_NWtoSE_AM[upstream[c],c]<-1
}

SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 14, c(17))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 16, c(14,17,19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 11, c(16,14,17,18,19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 12, c(11,14,15,16,17,18,19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 13, c(11,12,14,15,16,17,18,19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 6, c(11,12,13,14,15,16,17,18,19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 10, c(6,9,11:17, 19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 8, c(6,7,11:17, 19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 7, c(6,9:17, 19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 4, (6:19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 3, c(4,6:19))
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 2, c(3:17))  
SitesByEdges_NWtoSE_AM<-sitebyedges_expansion(SitesByEdges_NWtoSE_AM, 1, c(2:16))


#Tidy up
temp<-SitesByEdges_NWtoSE_AM
temp[temp>1]<-1
temp<-temp[,(length(upstream)+1):ncol(temp)]
temp->SitesByEdges_NWtoSE_AM
rm(temp)
write.csv(SitesByEdges_NWtoSE_AM, file="CLEANDATA/AEM_matrices/SitesByEdges_downstream_AM_rel50.csv")


############

#2 - SE to NW all connections for Acropora millepora, reliable>50% for 19 locations (upstream)

#make data frame
SitesByEdges_SEtoNW_AM<-data.frame(matrix(NA, nrow = 19, ncol = 0))

#define columns for fully downstream (end) locations
upstream<-c(1,2,5)
SitesByEdges_SEtoNW_AM[,1:length(upstream)]<-0
downstream<-"x"
for (c in 1:length(upstream)) {
  colnames(SitesByEdges_SEtoNW_AM)[c]<-paste(upstream[c], ">", downstream, sep="")
  SitesByEdges_SEtoNW_AM[upstream[c],c]<-1
}



SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 3, c(2,5))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 4, c(3,2,5))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 7, c(3,4,5))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 8, c(3,4,5,7))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 10, c(4,7,8))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 6, c(3,4,5))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 9, c(7,8,10))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 13, c(3:10))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 12, c(3:10,13))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 11, c(3:10,12,13))
SitesByEdges_SEtoNW_AM<-sitebyedges_expansion(SitesByEdges_SEtoNW_AM, 19, 17)


#Tidy up
temp<-SitesByEdges_SEtoNW_AM
temp[temp>1]<-1
temp<-temp[,(length(upstream)+1):ncol(temp)]
temp->SitesByEdges_SEtoNW_AM
rm(temp)
write.csv(SitesByEdges_SEtoNW_AM, file="CLEANDATA/AEM_matrices/SitesByEdges_upstream_AM_rel50.csv")

