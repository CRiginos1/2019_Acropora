#Initial reading and tidying of data - exported to genind  objects
#Written by Cynthia Riginos 2017-2019
#CR is running this on R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"


library(adegenet)

NA_num <- function(v) {sum(is.na(v))}    #counts the number of NA values

#READ IN FILES FOR ADEGENET
#tenuis
acro.data<-read.table("ORIGDATA/Atenuis_genind_1939_Inds_with_ReefID.txt", header =TRUE)
locations<-read.csv("ORIGDATA/ReefID_centroids_tenuis.csv", header =TRUE)  #these are lat and long for centers of each reef location
#millepora
#acro.data<-read.table("ORIGDATA/Amille_genind_918_Inds_with_ReefID.txt", header =TRUE)
#locations<-read.csv("ORIGDATA/ReefID_centroids_millepora.csv", header =TRUE)  #these are lat and long for centers of each reef location

#read in tables, convert to genind objects, and replace lat and long with values for reef centroids
acro.genind<-df2genind(acro.data[,6:ncol(acro.data)], sep=":", ind = acro.data[,1], pop= acro.data[,3] , type="codom")
for (r in 1:nrow(acro.data))
{
  acro.data[r,"centLat"]<-locations[acro.data[r,"ReefID"],"Lat"]
  acro.data[r,"centLong"]<-locations[acro.data[r,"ReefID"],"Long"]
}
acro.genind$other$latlong <- acro.data[,c("centLat", "centLong")]   #adds lat and long genind object


#HOW MUCH MISSING DATA BY LOCUS AND BY INDIVIDUAL?
acro.df<-genind2df(acro.genind, usepop = FALSE)
rownames(acro.df)<-indNames(acro.genind)
missingdata<-apply(acro.df, 1, NA_num)  #ID individuals with lots of missing data
hist(missingdata, main ="Missing data by individuals", xlab="Amount of data missing")
# tenuis: one individual with 4 loci missing; 10 with 3; 68 with 2; 220 with 1
# mille: 40 individuals with 1 missing locus: ignore

ind_cutoff<-2   #remove individuals with more than two loci missing
length(missingdata[missingdata>ind_cutoff]) 
highNAs<-indNames(acro.genind[missingdata>ind_cutoff])
acro.genind<-acro.genind[!indNames(acro.genind) %in% highNAs]  #AT: 1928 individuals; AM: 918 individuals


#look for low polymorphism loci (singletons, doubletons)
lowPM.loci<-which(!isPoly(acro.genind,by="allele", thres=20/(2*length(indNames(acro.genind)))))  #present in three copies or fewer
#both species: zero low PM loci

#SAVE and TIDY
save(acro.genind, file="CLEANDATA/tenuis/acroten.genind.RData")
#save(acro.genind, file="CLEANDATA/millepora/acromill.genind.RData")


rm(list=ls())
