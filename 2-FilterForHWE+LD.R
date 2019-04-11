#Testing for LD and HWE in populations
#See Waples 2014 Heredity for extensive discussion on reasonable approaches 
#Written by Cynthia Riginos 2017-2019

library(adegenet)
library(pegas)

ToNumeric<-function(x) {as.numeric(as.character(x))}

#Note - using pegas functions because they can handle multiallelic data

load(file="CLEANDATA/tenuis/acroten.genind.RData")
locations<-read.csv("ORIGDATA/ReefID_centroids_tenuis.csv", header =TRUE)

#load(file="CLEANDATA/millepora/acromill.genind.RData")
#locations<-read.csv("ORIGDATA/ReefID_centroids_millepora.csv", header =TRUE)


#HARDY WEINBERG
#Loop through populations and record HWE probs for each locus-population combination
pop.list<-seppop(acro.genind)
HWres<-data.frame()

for (p in 1:length(pop.list)) {
  pop<-pop.list[[p]]
  pop.loci<-genind2loci(pop)
  temp.hw.res<-as.data.frame(cbind(hw.test(pop.loci, B=0), names(pop.list)[p]))
  HWres<-rbind(HWres, temp.hw.res)
}

#Rejig and make sure numerics are numerics, etc
temp<-as.data.frame(sapply(HWres[,1:3], ToNumeric))
temp$V4<-HWres$V4
names(temp)<-c("chi^2", "df", "prob", "pop")
HWres<-temp
rm(temp) ; rm(temp.hw.res)
#Proportion of significant tests for all locus by population combinations
length(HWres[HWres$prob<0.05,])/nrow(HWres)  #AT: 0.010 tests are below 0.05 probability. AM: 0.021 below 0.05.  No need to worry about HW deviations!

#LINKAGE DISEQUILIBRIA
#LD: because we are using microsatellites we need to use a multi-allelic method
#LD2 method from pegas is appropriate
#However, NA's are not handled so the following script drops out individuals with NA's per pairwise LD calc

#pop.list<-seppop(acro)
LDres<-data.frame(stringsAsFactors=FALSE)
Popres<-data.frame()

#loop by populations and loci
for (r in 1:length(popNames(acro.genind))) {
  popx<-pop.list[[r]]
  
  #loop by loci
  for (loc1 in 1:(length(locNames(acro.genind))-1)){  #loop over locus 1
    for (loc2 in (loc1+1):length(locNames(acro.genind))) {
      #subset to focal loci and remove individuals with NA's
      temp.genind<-popx[,loc=c(loc1,loc2)]
      row.has.na <- apply(temp.genind$tab, 1, function(x){any(is.na(x))})
      temp.genind<-temp.genind[!row.has.na,]
      temp.loci<-genind2loci(temp.genind)
      
      #coerce results into a data frame, LDres
      #There is probably a better way of doing this but the mixture of data classes kept messing things up
      temp.ld.res<-as.list(c(as.numeric(as.character(LD2(temp.loci, details = FALSE)))))
      
      temp.ld.res<-as.data.frame(t(temp.ld.res))
      LDres<-rbind(LDres, temp.ld.res) 
      res<-as.data.frame(t(c(popNames(acro.genind)[r], locNames(acro.genind)[loc1], locNames(acro.genind)[loc2])))
      Popres<-rbind(Popres,res)
      
    }
  }
}

temp<-as.data.frame(sapply(LDres, ToNumeric))
temp<-cbind(Popres, temp)
LDres<-temp
names(LDres)<-c("pop", "loc1", "loc2", "T2", "df", "P")
#Proportion of significant tests for all locus by population combinations
length(LDres[LDres$P<0.05,])/nrow(LDres)  #AT: 0.003 tests are below 0.05 probability; AM: 0.007 No need to worry about LD either!
rm(temp); rm(temp.ld.res); rm(temp.loci); rm(res); rm(Popres)

#save files
write.csv(HWres, file="CLEANDATA/tenuis/HW_tests.csv")
write.csv(LDres, file="CLEANDATA/tenuis/LD_tests.csv")

# write.csv(HWres, file="CLEANDATA/millepora/HW_tests.csv")
# write.csv(LDres, file="CLEANDATA/millepora/LD_tests.csv")

rm(list=ls())


