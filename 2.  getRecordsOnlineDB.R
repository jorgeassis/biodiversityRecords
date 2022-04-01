# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
source("mainFunctions.R")
for(file in list.files("../sourceFunctions", full.names = TRUE)) { source(file) } 

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Open Marine Forests Database

source("https://raw.githubusercontent.com/jorgeassis/marineforestsDB/master/sourceMe.R")

dataset <- extractDataset("brownAlgae",pruned=TRUE)

unique(dataset$name)
dataset <- dataset[dataset$name != "Sargassum fluitans",]
dataset <- dataset[dataset$name != "Sargassum natans",]
dataset <- dataset[dataset$name != "Sargassum pusillum",c("decimalLongitude","decimalLatitude")]
dataset <- cbind(as.numeric(as.character(dataset[,1])),as.numeric(as.character(dataset[,2])))

plot(dataset)
write.table(dataset,file=paste0(directory,"/occurrenceRecordsMF.csv"), row.names = FALSE, quote=FALSE,sep = ";",col.names = TRUE)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Get Obis and GBIF records

dumpDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/"
dumpfile <- "biodiversityDataBaseOnline.csv"

myQueryFile <- read.csv(paste0(dumpDirectory,"biodiversityDataBaseLocal.csv"), sep=";")
myQuerySpecies <- myQueryFile[, "acceptedName"]
myQuerySpecies <- unique(myQuerySpecies)
myQuerySpecies <- myQuerySpecies[!is.na(myQuerySpecies)]
sum(countWords(myQuerySpecies) == 1)

tempDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/_ Sources/GBIF OBIS INAT/"

# ---------

rm(myQueryFile)
gc(reset=TRUE, full = TRUE)

# ---------

getFrom <- c("Gbif","Obis")

Cluster <- makeCluster( detectCores() ) 
registerDoParallel(Cluster)

parallelGetRecords <- foreach(sp=myQuerySpecies, .verbose=FALSE, .packages=c("stringi","plyr","RJSONIO","robis","raster","stringr","rgbif","rinat","dismo","data.table","dplyr")) %dopar% { 
  
  records.1 <- NULL
  records.2 <- NULL
  records.2 <- NULL
  
  if( ! file.exists(paste0(tempDirectory,"/",sp,".csv")) ) {
    
    if( "INaturalist" %in% getFrom ) { records.1 <- getExternalDataINaturalist(sp) }
    if(!is.null(records.1) ) { 
      records.1$bibliographicCitation <- ""
      records.1$datasetName <- "INaturalist"
    }
    
    # ---------
    
    if( "Obis" %in% getFrom ) { records.2 <- getExternalDataObis(sp) }
    if(!is.null(records.2)) { 
      if(nrow(records.2) > 0) { 
        records.2$datasetName <- "Ocean Biogeographic Information System"
      }
    }
    
    # ---------
    
    if( "Gbif" %in% getFrom ) { records.3 <-  getExternalDataGbif(sp) }
    if(!is.null(records.3) ) { 
      if(nrow(records.3) > 0) { 
        records.3$datasetName <- "Global Biodiversity Information Facility"  
      }
    }
    
    # ---------

    records <- data.frame()
    
    if( "INaturalist" %in% getFrom & !is.null(records.1)) { records <- rbind.fill(records,records.1) }
    if( "Obis" %in% getFrom & !is.null(records.2)) { records <- rbind.fill(records,records.2) }
    if( "Gbif" %in% getFrom & !is.null(records.3)) { records <- rbind.fill(records,records.3) }
    
    # ---------
    
    if( nrow(records) > 0 ) {
      
      records$name <- records$scientificName
      
      if(!is.null(records$date_year)) {
        records$year[is.na(records$year)] <- records$date_year[is.na(records$year)]
      }
      
      records <- records[which(!is.na(records$decimalLongitude) & !is.na(records$decimalLatitude) ) , ]
      records <- records[,-which(colnames(records) %in% c("typeStatus","taxonRemarks","disposition","higherGeography","nomenclaturalCode","collectionKey","typifiedName","suborderid","suborder","aphiaID","dynamicProperties","waterBody","otherCatalogNumbers","municipality","continent","islandGroup","island","taxonomicStatus","taxonConceptID","vernacularName","taxonID","higherClassification","specificEpithet","taxonRank","scientificNameAuthorship","species","class","originalScientificName","brackish","order","scientificName","dropped","date_year","brackish","genusid","species","scientificNameID","institutionKey","gbifID","inCluster","isInCluster","lastInterpreted","iucnRedListCategory","genericName","acceptedTaxonKey","speciesKey","acceptedScientificName","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey","taxonKey","hostingOrganizationKey","crawlId","lastParsed","lastCrawled","protocol","publishingCountry","installationKey","publishingOrgKey","datasetKey","key","publishingOrgKey","forma","formaid","footprintWKT","flags","sst","shoredistance","sss","node_id","bathymetry","fieldNumber","terrestrial","modified","familyid","language","subclassid","phylumid"   ,"basisOfRecord","id","speciesid","genus","phylum","subclass","family","kingdomid" ,"kingdom","classid","orderid","date_mid","marine"))]
      
      records <- unique(records)
      rownames(records) <- NULL

      for(c in colnames(records)) { records[,c] <- gsub(";",",",records[,c]) }
            
      write.csv(records,file=paste0(tempDirectory,"/",sp,".csv"),row.names = FALSE)

    }
    
  }
  
  records.1 <- NULL
  records.2 <- NULL
  records.2 <- NULL
  gc(reset=TRUE)
  
  return(NULL)
  
}

stopCluster(Cluster) ; rm(Cluster)
closeAllConnections()

## ---------------------------------------------------
## ---------------------------------------------------
## Merge databases

rm(list=(ls()[ls()!="v"]))
source("mainFunctions.R")
for(file in list.files("../sourceFunctions", full.names = TRUE)) { source(file) } 

dumpDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/"
dumpfile <- "biodiversityDataBaseOnline.csv"

# flaged species : with data for known depth and distributional range
flagedAcceptedName <- read.csv(file="../../../Data/Biodiversity Data/Cold Water Corals/acceptedNameToKeep.csv")[,1]

tempDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/_ Sources/GBIF OBIS INAT/"
databaseExternalFiles <- list.files(tempDirectory, full.names = T, pattern=".csv")

myDBStructure <- c("kingdom","phylum","class","order","family","genus","name","acceptedName","scientificNameAuthorship","taxonomicStatus","aphiaID","decimalLongitude","decimalLatitude","coordinatePrecision","coordinateUncertaintyInMeters","georeferenceProtocol","verbatimDepth","minimumDepthInMeters","maximumDepthInMeters","depthAccuracy","country","locality","year","month","day","absence","individualCount","identifiedBy","basisOfRecord","datasetName","license","bibliographicCitation")

for(f in databaseExternalFiles) { 
  
  cat(which(databaseExternalFiles == f),"out of",length(databaseExternalFiles),"\n")
  
  database.i <- read.csv(f)
  colnames(database.i)[which(colnames(database.i) == "depth") ] <- "verbatimDepth"
  database.i <- database.i[, unlist(sapply( myDBStructure , function(x) which(colnames(database.i) == x) )) ]
  
  if( TRUE %in% (unique(database.i$name) %in% speciesInDB) ) { next }
  
  taxonomy <- getTaxonomyWorms(unique(database.i$name), verbose=FALSE)$acceptedName
  taxonomy <- c(unique(database.i$name),taxonomy[!is.na(taxonomy)])
  
  if( ! TRUE %in% (taxonomy %in% flagedAcceptedName) ) { next }
  
  if( file.exists(paste0(dumpDirectory,dumpfile)) ) { 
    
    database <- read.csv(file=paste0(dumpDirectory,dumpfile)) 
    database <- rbind.fill(database,database.i)
    
    }
  
  if( ! file.exists(paste0(dumpDirectory,dumpfile)) ) { 
    
    database <- database.i
  
  }
  
  speciesInDB <- unique(database$acceptedName)
  write.csv(database,file=paste0(dumpDirectory,dumpfile),row.names = FALSE)   
  rm(database); rm(database.i)
  gc(reset=TRUE)
    
}

## ------------------

database <- read.csv(file=paste0(dumpDirectory,dumpfile))
nrow(database)
colnames(database)

## ------------------

database$verbatimDepth[which(is.na(database$verbatimDepth))] <- database$verbatimDepth.1[which(is.na(database$verbatimDepth))] 
database <- database[,-which(colnames(database) == c("verbatimDepth.1")),]
database$basisOfRecord <- "Biodiversity Information Facility"
colnames(database)

## ------------------
## accessRights

unique(database$license)
inproperLicences <- c("http://data.gc.ca/eng/open-government-licence-canada & http://www.canadensys.net/norms")
database <- database[-which(database$license %in% inproperLicences),]
unique(database$license)

## ------------------------------------------------------------
## ------------------------------------------------------------

database$decimalLongitude <- as.numeric(as.character(database$decimalLongitude) )
database$decimalLatitude <- as.numeric(as.character(database$decimalLatitude ) )
database$year <- as.numeric(as.character(database$year ) )
database$month <- as.numeric(as.character(database$month ) )
database$day <- as.numeric(as.character(database$day ) )
database$coordinateUncertaintyInMeters <- as.numeric(as.character(database$coordinateUncertaintyInMeters ) )

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

database <- database[which( !is.na(database$acceptedName) ),]
database <- database[which( countWords(database$acceptedName) > 1 ),]
nrow(database)

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