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

#because some locations lumped together we drop some numbers for AT
excluded<-c(24,26,39,41)

#######
#1 - NWtoSE all connections for Acropora tenuis in the top 50% of reliability for 39 locations (downstream)

#make data frame
SitesByEdges_NWtoSE_AT<-data.frame(matrix(NA, nrow = 43, ncol = 0))

#define columns for fully downstream (end) locations
upstream<-c(43, 40)
downstream<-"x"
SitesByEdges_NWtoSE_AT[,1:length(upstream)]<-0
for (c in 1:length(upstream)) {
  colnames(SitesByEdges_NWtoSE_AT)[c]<-paste(upstream[c], ">", downstream, sep="")
  SitesByEdges_NWtoSE_AT[upstream[c],c]<-1
}


#Swains
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 42, c(40, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 36, c(40,42,43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 37, c(33,34,36,40,42,43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 34, c(36,40,42,43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 35, c(33,34,36,37,40,42,43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 33, c(34,35,36,37,38,40,42,43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 38, c(33,34,35,36,37,40,42,43))

#Southern outer, adding in Heron & Keppel as downstream
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 31, 43)
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 32, c(31, 33:38, 40, 42, 43))

SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 19, c(31:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 27, c(19,31:38, 40, 42, 43))

#Central - outer and midshelf
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 18, c(19, 27, 31:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 17, c(18, 19, 27, 31:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 30, c(19, 27, 31:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 29, c(17:19,27, 31:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 23, c(19,27, 30:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 28, c(17:19,27, 29,31:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 22, c(19,27, 30:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 20, c(17:19, 23,27:35,38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 7, c(17:19, 27:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 25, c(17, 18, 27:30))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 21, c(17:19, 22, 23, 27:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 6, c(7,17:20, 23,27:35,37, 38, 40, 42, 43))

#Northern
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 12, c(6,7, 17:23, 25, 27:32, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 9, c(6,7, 17:23, 25, 27:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 8, c(6,7,9, 17:23, 25, 27:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 10, c(6,7, 12, 17:23, 25, 27:30, 32,  43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 14,  c(6:9, 17:23, 25, 27:38, 40, 42, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 16, c(6,7, 10, 12, 17:23, 25, 27:30,32,40, 42, 43))

#Far northern
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 11, c(6:10, 12, 14, 16:23, 25, 27:32, 43))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 13, c(6:12, 14, 16:23, 25, 27:33, 43))

SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 5, c(6, 10, 11, 12, 15:18, 21, 22, 25, 27:30))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 3, c(5, 11,15))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 4, c(5, 11,15))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 2, c(3:5, 8,9, 11, 13:16))
SitesByEdges_NWtoSE_AT<-sitebyedges_expansion(SitesByEdges_NWtoSE_AT, 1, c(2:5, 8:16, 18, 22, 30))


#Tidy up
temp<-SitesByEdges_NWtoSE_AT
temp<-temp[! c(1:43) %in% excluded,]
temp[temp>1]<-1
temp<-temp[,(length(upstream)+1):ncol(temp)]
temp->SitesByEdges_NWtoSE_AT
rm(temp)
write.csv(SitesByEdges_NWtoSE_AT, file="CLEANDATA/AEM_matrices/SitesByEdges_NWtoSE_downstream_AT_rel50.csv")


############
#2 - SE to NW all connections for Acropora tenuis (upstream) in the top 50% of reliability

#make data frame
SitesByEdges_SEtoNW_AT<-data.frame(matrix(NA, nrow = 43, ncol = 0))

#define columns for fully downstream (end) locations
upstream<-c(1,2)  #these are actually the fully downstream locations - ends of the line
downstream<-"x"  #dummy downstream population
SitesByEdges_SEtoNW_AT[,1:length(upstream)]<-0

for (c in 1:length(upstream)) {
  colnames(SitesByEdges_SEtoNW_AT)[c]<-paste(upstream[c], ">", downstream, sep="")
  SitesByEdges_SEtoNW_AT[upstream[c],c]<-1
}


#Far northern
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 3, 1) 
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 4, 1) 
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 5, 4) 
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 13, 5) 
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 11, 5) 

#northern
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 14, c(5,11,15)) 
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 16, c(5,11,15)) 
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 8, c(5,11,14,15)) 
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 10, c(5,11, 14:16)) 
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 9, c(5,8,11, 14, 15))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 12, c(8,10,11, 14:16))

#central
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 21, c(8:12, 14:16))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 6, c(8, 9, 11, 14:16))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 25, 14)
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 20, c(6, 8,9, 11, 14:16))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 7, c(8:11, 14:16, 21))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 22, c(8, 9, 11, 14:16, 21))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 23, c(8, 9, 14:16, 21, 22))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 28, c(6:12, 14:16, 20, 21))

SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 30, c(6:12, 14:16, 20:23, 25))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 29, c(6:9, 14:16, 20, 21, 28))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 17, c(6:9, 14:16, 20, 21, 28,29))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 18, c(6:9, 14:17, 20, 21, 28, 29))

SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 27, c(6:12, 14:18, 20:23, 25, 28:30))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 19, c(6:12, 14:18,20:23, 25, 27:30))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 32, c(6,7, 14, 17:23, 25, 27:30))
#31 is a sink

#Swains, keppels, etc
#43 and 40 are sinks
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 34, c(35, 37))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 36, c(35, 37))
SitesByEdges_SEtoNW_AT<-sitebyedges_expansion(SitesByEdges_SEtoNW_AT, 42, c(35, 37))


#Tidy up
temp<-SitesByEdges_SEtoNW_AT
temp<-temp[! c(1:43) %in% excluded,]
temp[temp>1]<-1
temp<-temp[,(length(upstream)+1):ncol(temp)]
temp->SitesByEdges_SEtoNW_AT
rm(temp)
write.csv(SitesByEdges_SEtoNW_AT, file="CLEANDATA/AEM_matrices/SitesByEdges_SEtoNW_upstream_AT_rel50.csv")


rm(list=ls())
