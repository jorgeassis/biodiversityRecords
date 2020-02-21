# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
#
# ------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
source("Dependencies.R")

# ------------------------------------------------------------------------------------

phylum <- "Cnidaria"
taxa <- "Alcyonacea" # 
taxa <- "Pennatulacea" # 
taxa <- "Helioporacea" # 
taxa <- "Gorgonacea" # 

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Get Taxonomy by higher ranked level

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Future vulnerability of rhodolith beds under multiple environmental stressors/Data/Occurrence records/"
taxaAll <- c("Lithothamnion" , "Lithophyllum" , "Mesophyllum" )
rank <- "higherTaxonomicRank"

# -------------------------------

Cluster <- makeCluster( 16 ) 
registerDoParallel(Cluster)

species <- foreach(taxa=taxaAll, .verbose=FALSE,  .combine=rbind, .packages=c("worms","worrms")) %dopar% {  return(getTaxonomyNewSpecies(taxa,rank)) }

stopCluster(Cluster) ; rm(Cluster)

species <- species[which( sapply(1:nrow(species), function(x) { TRUE %in% (taxaAll %in% species[x,]) } ) ),]
species <- species[species$status == "accepted",]
nrow(species)
colnames(species)

write.table(species,file=paste0(directory,"speciesDatabase.csv"), row.names = FALSE, quote=FALSE,sep = ";",col.names = TRUE)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Get species taxonomy

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Future vulnerability of rhodolith beds under multiple environmental stressors/Data/Occurrence records/"
fileSpecies <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Future vulnerability of rhodolith beds under multiple environmental stressors/Data/Occurrence records/speciesDatabase.csv"

species <- read.table(file=fileSpecies, sep = ";", stringsAsFactors = FALSE)[,1]
taxonomy <- do.call("rbind",sapply(species,function(x) { wormsbynames(taxon_names=x) } ))
speciesDB <- cbind(taxa ,  )

write.table(speciesDB,file=paste0(directory,"speciesDatabase.csv"), row.names = FALSE, quote=FALSE,sep = ";",col.names = TRUE)
