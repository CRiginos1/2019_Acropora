# 2019_Acropora

##Asymmetric dispersal is a critical element of concordance between biophysical dispersal models and spatial genetic structure in Great Barrier Reef corals

## Genetic data filtering, migration directions, spatial eigenvector mapping

Scripts and data to support analyses from:

Riginos, Hock, Matias, Mumby, van Oppen & Lukoschek. 2019.  *Diversity & Distributions*



*Electronic notebook author:* Cynthia Riginos
***


##Organisation and workflow

###Original microsatellite data sources
*Acropora tenuis*: Lukoschek, V., Riginos, C., & van Oppen, M. J. H. (2016). Congruent patterns of connectivity can inform management for broadcast spawning corals on the Great Barrier Reef. Molecular Ecology, 25(13), 3065–3080. 

*Acropora millepora*: van Oppen, M. J. H., Peplow, L. M., Kinnimonth, S., & Berkelmans, R. (2011). Historical and contemporary factors shape the population genetic structure of the broadcast spawning coral, *Acropora millepora*, on the Great Barrier Reef. Molecular Ecology, 20(23), 4899–4914. 

###Manual preprocessing
To align the genetic data against the larval dispersal model:
For *A. tenuis*: we excluded three locations from the original study (22_SN_PIT, 23_HIGH_RK, 68_MAG_CAY) and merged some sites within reefs or islands as follows: Yonge (sites 10+11); No Name (sites 12+13); Palm Islands (sites 38-41,44); Heron and One Tree Island (sites 60-65). 
For *A. millepora*, we merged Sudbury 1 with Sudbury 2 and removed one locus (Am2_002) due to high missing data.

Manual preprocessed data files are in ORIGDATA folder.

### DivMigrate
To estimate migration following the divMigrate approach (Sundqvist et al. 2016. *Ecology and Evolution*, 6(11), 3461–3475), genind formatted R files were exported (adegenet::genind2df) and then transformed into genpop format using PDGSpider ver. 2.1.0.0 (Lischer & Excoffier. 2012 *Bioinformatics* 28: 298-299). Resultant genpop formatted files were imported into the online version of divMigrate (https://popgen.shinyapps.io/divMigrate-online/) and the Gst option was used to create pairwise migration estimates, with the resultant file exported.


##R Files
Most can be run as scripts. Use hash tag to toggle between *A. tenuis* and *A. millepora* analsyes to read in appropriate files (top of script) and write results (bottom of script).

###1 - ReadInFiles-Convert2Genind.R
Converts .txt genotype files to genind objects; replaces field georeferencing with reef ID centroids to align to larval dispersal model distances; filters out individuals with high missing data. Resultant genind files are placed in CLEANDATA folder (acroten.genind.RData and acromill.genind.RData).

###2 - FilterForHWE+LD.R
Genetic data are checked for systematic or large deviations from Hardy Weinberg expectation and for linkage disequilibria. Since no large or systemic patterns is detected, no further genetic data filtering is undertaken.

###3 - GeographicInformation.R
The geographic information for sites with genetic data are characterised. This includes assigning color values by latitudinal position and creating ordinations of distance to coastline and position relative to coastline (roughly NW to SE, using 50m coastline shape file from Natural Earth: www.naturalearthdata.com). Resultant dataframe files are placed in CLEANDATA folder (locations.aten.RData and locations.mill.RData).

####4 - DirectedConnectivityCorrelations.R
Distance matrices of connectivity and distance between sites are imported, turned into vectors, and organised in a data frame. Correlations between biophysical measures and divMigrate estimates are evaluated. Resultant dataframe files are placed in CLEANDATA folder (bpdist.unfolded.tenuis39.csv and bpdist.unfolded.mill19.csv).

####5 - MEMs.R
Ability of Moran Eigenvector Models to describe spatial genetic structure is evaluated largely based on approaches described in detail by Borcard, Gillet, & Legendre. 2011. Numerical Ecology with R. Springer. Distance models include binary Delauney, distance scaled (weighted) Delauney, PCNM (principal coordinate analysis of neighbour matrices), and a custom saturated model. Result summary files are placed in the CLEANDATA folder (tenuis/MEMresults.txt and millepora/MEMresults.txt).

####6a-SiteByEdgesConstruction_AEM_tenuis>50percentreliable
####6b-SiteByEdgesConstruction_AEM_millepora>50percentreliable
Site by edges matrix construction in preparation for AEMs follows Blanchet, F. G., Legendre, P., Maranger, R., Monti, D., & Pepin, P. (2011). Modelling the effect of directional spatial ecological processes at different scales. Oecologia, 166(2), 357–368. http://doi.org/10.1007/s00442-010-1867-y. Only connections in the top 50% of reliable connections are retained. Site by edges matrices are made both for species and for N->S and S->N connections. Resultant matrix files are placed in the CLEANDATA/AEM\_matrices folder (SitesByEdges\_NWtoSE\_downstream\_AT\_rel50.csv, SitesByEdges\_NWtoSE_downstream\_AT\_rel50.csv, SitesByEdges\_downstream\_AM\_rel50.csv, SitesByEdges\_upstream\_AM\_rel50.csv).

####7-EdgeWeightConstructionAEM.R
Creates edge weightings in preparation for AEMs. Links are extracted from site by edges matrices (#6a & b)and distance weightings are extracted from unfolded distance dataframes (#4). Resultant files are placed in the CLEANDATA/AEM\_matrices folder with the same naming convention as #6 but with "_weights" appended to the name.

####8-AEManalyses.R
Conducts the asymmetric eigenvector analyses based on a variety of distance weightings and using connections N->S and in both directions (N->S and S->N). Summary result files are placed in CLEANDATA/tenuis and CLEANDATA/millepora folders.


***

