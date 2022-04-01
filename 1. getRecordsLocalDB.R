# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
#
# ------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
source("mainFunctions.R")
for(file in list.files("../sourceFunctions", full.names = TRUE)) { source(file) } 

dumpDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/"
dumpfile <- "biodiversityDataBaseLocal.csv"

## ---------------------
## ---------------------

## 1. Get species records

sourceDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/_ Sources/"

database <- sf::st_read(paste0(sourceDirectory,"/UNEP-WCMC/WCMC001_ColdCorals2017_Pt_v5_1.shp"))
database <- data.frame(st_coordinates(database)[which(!duplicated(st_coordinates(database)[,3])),1:2],data.frame(database)[,-which(colnames(data.frame(database))=="geometry")])
head(database)

database <- database[,c("X","Y","SPECIES","START_DATE")]
colnames(database) <- c("decimalLongitude","decimalLatitude","name","date")
database$bibliographicCitation <- "Freiwald A, Rogers A, Hall-Spencer J, Guinotte JM, Davies AJ, Yesson C, Martin CS, Weatherdon LV (2021). Global distribution of cold-water corals (version 5.1). Fifth update to the dataset in Freiwald et al. (2004) by UNEP-WCMC, in collaboration with Andre Freiwald and John Guinotte. Cambridge (UK): UN Environment Programme World Conservation Monitoring Centre. DOI: https://doi.org/10.34892/72×9-rt61"
database$datasetName <- "Global Distribution of Cold-water Corals"
database$basisOfRecord <- "Biodiversity Information Facility"
database$year <- as.numeric(substr(database$date,7,10))
database$month <- as.numeric(substr(database$date,4,5))
database$day <- as.numeric(substr(database$date,1,2))
database <- database[,-which(colnames(database) %in% c("date"))]

## Remove no data
database <- database[database$name != "Not Reported",]

## Add sp. to one word species name
speciesCountWords <- sapply(database$name,function(x) { countWords(x) } )
database$name[which(speciesCountWords == 1)] <- paste0(database$name[which(speciesCountWords == 1)]," sp.")

## Change multiple assignments
database$name <- str_replace_all(database$name, " sp. cf.", " sp.")
database$name <- str_replace_all(database$name, " cf. sp.", " sp.")
database$name <- str_replace_all(database$name, " spp.", " sp.")
database$name <- str_replace_all(database$name, " spp", " sp.")
database$name <- str_replace_all(database$name, " n. sp. [123456789]", " sp.")
database$name <- str_replace_all(database$name, " n. sp. [ABCDEFGHIJ]", " sp.")
database$name <- str_replace_all(database$name, " sp. [ABCDEFGHIJ]", " sp.")
database$name <- str_replace_all(database$name, " sp\\. [ABCDEFGHIJ]", " sp.")
database$name <- str_replace_all(database$name, "\\?", "")
database$name <- str_replace_all(database$name, "sp./6", "sp.")
database$name <- str_replace_all(database$name, "[123456789]", "")
database$name <- str_replace_all(database$name, " ,", "")
database$name <- gsub("\\s*\\([^\\)]+\\)","",database$name)
database$name <- gsub("sp\\..*","sp.",database$name)
database$name <- gsub("c\\.f\\..*","sp.",database$name)
database$name <- gsub("cf\\..*","sp.",database$name)
database$name <- gsub(",","",database$name)
database$name <- str_trim(database$name, "right")
database$name <- str_trim(database$name, "left")
database$name[which(sapply(database$name,function(x) { str_locate_all(x, "sp")[[1]][2] - str_locate_all(x, "sp")[[1]][1] == 1 & str_locate_all(x, "sp")[[1]][2] == nchar(x) } ))] <- paste0(database$name[which(sapply(database$name,function(x) { str_locate_all(x, "sp")[[1]][2] - str_locate_all(x, "sp")[[1]][1] == 1 & str_locate_all(x, "sp")[[1]][2] == nchar(x) } ))],".")
database$name <- word(database$name , 1 , 2)
colnames(database)

## ---------------------
## ---------------------

# get data from additional databases

sourceDirectory <- "../../../Data/Genetic Data/Genetic Data [Literature Review]/Coldwater corals/Data/"
sourceFiles <- list.files(sourceDirectory, recursive = TRUE, full.names = TRUE, pattern = "_data")

database.literature <- data.frame()

for(f in sourceFiles) {
  
  data.i <- read.xls(f)
  database.literature <- rbind(database.literature,data.i[,c("Species","Lon","Lat","depth","Reference")])

}

colnames(database.literature) <- c("name","decimalLongitude","decimalLatitude","verbatimDepth","bibliographicCitation")
database.literature$basisOfRecord <- "Peer-reviewed literature"
colnames(database.literature)

## ---------------------
## ---------------------

sourceDirectory <- "../../../Data/Biodiversity Data/Occurrence records/Corals/"
sourceFiles <- list.files(sourceDirectory, recursive = TRUE, full.names = TRUE, pattern = ".csv")

database.literature2 <- data.frame()

for(f in sourceFiles) {
  
  data.i <- read.csv(f, sep=";")
  database.literature2 <- rbind(database.literature2,data.i[,c("speciesName","Lon","Lat","depth","minimumDepthInMeters","maximumDepthInMeters","country","locality","year","month","day","originalSourceType","originalSource","originalSourceDOI")])
  
}

colnames(database.literature2)[which(colnames(database.literature2) %in% c("Lon"))] <- c("decimalLongitude")
colnames(database.literature2)[which(colnames(database.literature2) %in% c("Lat"))] <- c("decimalLatitude")
colnames(database.literature2)[which(colnames(database.literature2) %in% c("originalSource"))] <- c("bibliographicCitation")
colnames(database.literature2)[which(colnames(database.literature2) %in% c("speciesName"))] <- c("name")
colnames(database.literature2)[which(colnames(database.literature2) %in% c("depth"))] <- c("verbatimDepth")

database.literature2 <- database.literature2[database.literature2$originalSourceType != "ExternalDB",]
database.literature2$basisOfRecord <- "Peer-reviewed literature"
database.literature2$bibliographicCitation <- paste0(database.literature2$bibliographicCitation," DOI: ",database.literature2$originalSourceDOI)
database.literature2 <- database.literature2[,-which(colnames(database.literature2) %in% c("originalSourceType","originalSourceDOI")),]
database.literature <- rbind.fill(database.literature2,database.literature2)
colnames(database.literature)

## ---------------------
## ---------------------

sourceDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/_ Sources/"
database.noaa <- read.csv(paste0(sourceDirectory,"NOAA/deep_sea_corals_6cb7_510e_0a83.csv"))
database.noaa <- database.noaa[-1,]
database.noaa$year <- as.numeric(substr(database.noaa$ObservationDate,1,4))
database.noaa$month <- as.numeric(substr(database.noaa$ObservationDate,6,7))
database.noaa$day <- as.numeric(substr(database.noaa$ObservationDate,9,10))

database.noaa$bibliographicCitation <- paste0("Hourigan T (2020). NOAA Deep Sea Corals Research and Technology Program. Version 1.6. United States Geological Survey. Data provider: ",database.noaa$DataProvider)
database.noaa$datasetName <- "NOAA Deep Sea Corals Research and Technology Program"
database.noaa$basisOfRecord <- "Biodiversity Information Facility"

database.noaa <- database.noaa[,-which(colnames(database.noaa) %in% c("DataProvider","ObservationDate","CatalogNumber","TaxonRank","Station","VernacularNameCategory","Repository","SurveyID","DepthMethod","SamplingEquipment","RecordType","ShallowFlag","EventID","SampleID","DatasetID"))]
colnames(database.noaa) <- c("name","decimalLatitude","decimalLongitude","verbatimDepth","locality","coordinatePrecision","identifiedBy","year","month","day","bibliographicCitation","datasetName","basisOfRecord")
colnames(database.noaa)

## ---------------------
## ---------------------

sourceDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/_ Sources/"
database.emodnet <- read.csv(paste0(sourceDirectory,"EMODNET/Emod_net_all.csv"), sep=";")
database.emodnet <- database.emodnet[,1:15]
colnames(database.emodnet) <- c("name","decimalLongitude","decimalLatitude","georeferenceProtocol","coordinateUncertaintyInMeters","country","locality","verbatimDepth","minimumDepthInMeters" ,"maximumDepthInMeters","year","month","day","Remove1","bibliographicCitation")
database.emodnet <- database.emodnet[,-which(colnames(database.emodnet) %in% c("Remove1"))]

database.emodnet$bibliographicCitation <- "European Marine Observation Data Network (EMODnet) Biology project (www.emodnet-biology.eu), funded by the European Commission’s Directorate - General for Maritime Affairs and Fisheries (DG MARE)."
database.emodnet$datasetName <- "EMODnet Biology"
database.emodnet$basisOfRecord <- "Biodiversity Information Facility"
colnames(database.emodnet)

## ------------

availableSpecies <- unique(c(database$name,database.literature$name,database.literature2$name,database.noaa$name))
database.emodnet <- database.emodnet[which(database.emodnet$name %in% availableSpecies) , ]

## ---------------------------------------------------
## ---------------------------------------------------
## Merge databases

database <- rbind.fill(database,database.literature,database.noaa,database.emodnet)
colnames(database)

## ------------------------------------------------------------
## ------------------------------------------------------------

# Add Aquamaps records

filename <- "../../../Data/Biodiversity Data/Occurrence records/Biodiversity database Aquamaps.db"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver, dbname = filename)
dbListTables(db)
taxa <- dbReadTable(db,"taxa")
taxaSp <- paste0(taxa$Genus," ",taxa$Species)
nativemaps <- dbReadTable(db,"nativemaps")
cells <- dbReadTable(db,"hcaf")
occurrence <- dbReadTable(db,"occ")

speciesToAddRecords <- sort(unique(database$name))
database.aquamaps <- data.frame()

for( sp in speciesToAddRecords) {
  
  cat(sp,"\n")
  sp.id <- taxa[which(taxaSp == sp),"SPECIESID"]
  
  if( length(sp.id) == 0) { next }

  sp.records <- occurrence[occurrence$SpeciesID == sp.id & occurrence$GoodCell == 1,c("CenterLong","CenterLat")]
  database.aquamaps <- rbind(database.aquamaps,data.frame( decimalLongitude=sp.records[,1],decimalLatitude=sp.records[,2] , name=sp))
  
}

database.aquamaps$bibliographicCitation <- "Kaschner, K., J. Rius-Barile, K. Kesner-Reyes, C. Garilao, S.O. Kullander, T. Rees and R. Froese. 2010. AquaMaps: Predicted range maps for aquatic species. World wide web electronic publication, www.aquamaps.org, Version 08/2010."
database.aquamaps$datasetName <- "AquaMaps"
database.aquamaps$basisOfRecord <- "Biodiversity Information Facility"
head(database.aquamaps)
nrow(database.aquamaps)

## ---------------------------------------------------
## ---------------------------------------------------
## Merge databases

database <- rbind.fill(database,database.aquamaps)
colnames(database)

## ------------------------------------------------------------
## ------------------------------------------------------------

database$decimalLongitude <- as.numeric(as.character(database$decimalLongitude) )
database$decimalLatitude <- as.numeric(as.character(database$decimalLatitude ) )
database$year <- as.numeric(as.character(database$year ) )
database$month <- as.numeric(as.character(database$month ) )
database$day <- as.numeric(as.character(database$day ) )
database$coordinateUncertaintyInMeters <- as.numeric(as.character(database$coordinateUncertaintyInMeters ) )

# Remove records at the genus level

database <- database[-which(grepl("sp\\.",database$name)),]

## ------------------------------------------------------------
## ------------------------------------------------------------

## Clean by worms; checking taxonomy

speciesList <- unique(database$name)
taxonomyDataBase <- data.frame()

for( sp in speciesList ) { taxonomy <- getTaxonomyWorms(sp, verbose=FALSE);  if( ! is.null(taxonomy) ) { taxonomyDataBase <- rbind(taxonomyDataBase,taxonomy) } }

taxonomyDataBase <- rbind(taxonomyDataBase[1,],taxonomyDataBase)
taxonomyDataBase[1,] <- NA
colnames(taxonomyDataBase) <- c("aphiaID","name","scientificNameAuthorship","taxonomicStatus","kingdom","phylum","class","order","family","genus","revisionByWormsDate","acceptedAphiaID","acceptedName")
head(taxonomyDataBase)

# Inject taxonomy

database <- cbind(database,
                  data.frame(taxonomyDataBase[ sapply(database$name, function(x) { row <- which( taxonomyDataBase$name == x); ifelse(length(row) == 0, 1, row)}),-which(colnames(taxonomyDataBase) %in% c("name","revisionByWormsDate","acceptedAphiaID"))]))

head(database)
colnames(database)

## ------------------------------------
## ------------------------------------

# Remove unresolved or wrong taxa

taxClass <- c("Anthozoa","Hydrozoa")
database <- database[database$class %in% taxClass,]

# Species level only

database <- database[which( !is.na(database$acceptedName) ),]
database <- database[which( countWords(database$acceptedName) > 1 ),]

## ------------------------------------
## ------------------------------------

# data.frame to DarwinCore

source("https://raw.githubusercontent.com/jorgeassis/marineforestsDB/master/sourceMe.R")
darwinColumns <- extractDataset("brownAlgae",pruned=TRUE)
darwinColumns <- darwinColumns[1,]; darwinColumns[1,] <- NA
darwinColumns <- colnames(darwinColumns)

colnames(database)[which( ! colnames(database) %in% darwinColumns )]
unique(database$basisOfRecord)

## -----------
## -----------

for(c in colnames(database) ) { database[,c] <- gsub(";",",",database[,c]) }

head(database)
length(unique(database$acceptedName))

# Export data.frame

write.csv(database,file=paste0(dumpDirectory,dumpfile),row.names = FALSE)

## ---------------------------------------------------
## ---------------------------------------------------