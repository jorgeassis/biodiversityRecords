# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
#
# ------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
source("mainFunctions.R")
for(file in list.files("../sourceFunctions", full.names = TRUE)) { source(file) } 

# ---------------------------
# ---------------------------

dumpDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/"

database <- rbind.fill(read.csv(file=paste0(dumpDirectory,"biodiversityDataBaseLocal.csv")),
                       read.csv(file=paste0(dumpDirectory,"biodiversityDataBaseOnline.csv")))

myDBStructure <- c("kingdom","phylum","class","order","family","genus","name","acceptedName","scientificNameAuthorship","taxonomicStatus","aphiaID","absence","decimalLongitude","decimalLatitude","coordinatePrecision","coordinateUncertaintyInMeters","georeferenceProtocol","verbatimDepth","minimumDepthInMeters","maximumDepthInMeters","depthAccuracy","country","locality","year","month","day","individualCount","identifiedBy","basisOfRecord","datasetName","bibliographicCitation")
database <- database[,unlist(sapply( myDBStructure , function(x) which(colnames(database) ==x) ))]
colnames(database)
nrow(database)

database <- database[which(!is.na(database$decimalLongitude)),]
database <- database[which(!is.na(database$decimalLatitude)),]

## ---------------------
## ---------------------

database$MeasurementOrFact <- '{"flagLand": "0", "flagVerticalRange": "0", "flagDistributionalRange": "0"}'

## ---------------------
## ---------------------
## Remove duplicated records

database <- database[-which(duplicated(database[,c("decimalLongitude","decimalLatitude","verbatimDepth","acceptedName","year","month","day")])),]
nrow(database)

## ---------------------
## ---------------------
## Flag records 

# Depth range and over land

speciesToFlag <- sort(unique(database$acceptedName))

bathymetry <- raster("../../../Data/Spatial information/Rasters/GEBCO Bathymetry Global.tif")
knownDepthRange <- read.csv(file="../../../Data/Biodiversity Data/Cold Water Corals/VerticalDistributions/coldWaterCorals.csv")
knownDistributionalRange <- read.csv(file="../../../Data/Biodiversity Data/Cold Water Corals/expertReview/expertReview.csv", sep=";", header = TRUE)

faoAreas <- shapefile("dataSources/faoAreas/World_Fao_Areas.shp")

filename <- "../../../Data/Biodiversity Data/Occurrence records/Biodiversity database Aquamaps.db"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver, dbname = filename)
dbListTables(db)
taxa <- dbReadTable(db,"taxa")
taxaSp <- paste0(taxa$Genus," ",taxa$Species)
nativemaps <- dbReadTable(db,"nativemaps")
cells <- dbReadTable(db,"hcaf")

# https://www.iucnredlist.org/resources/spatial-data-download
db1 <- sf::st_read("../../../Data/Biodiversity Data/Cold Water Corals/_ Sources/IUCN Corals/REEF_FORMING_CORALS_PART1.shp")
db2 <- sf::st_read("../../../Data/Biodiversity Data/Cold Water Corals/_ Sources/IUCN Corals/REEF_FORMING_CORALS_PART2.shp")
db3 <- sf::st_read("../../../Data/Biodiversity Data/Cold Water Corals/_ Sources/IUCN Corals/REEF_FORMING_CORALS_PART3.shp")

for( sp.i in 1:length(speciesToFlag)  ) {
  
  sp <- speciesToFlag[sp.i]
  
  cat(sp, "::" , which(speciesToFlag == sp),"out of",length(speciesToFlag),"\n")

  # grepl(gsub("([A-Za-z]+).*", "\\1", sp),speciesToFlag)
  
  database.i <- database[which(database$acceptedName == sp),]
  bathymetry.i <- as.numeric(extract(bathymetry,database.i[,c("decimalLongitude","decimalLatitude")]))

  if(! NA %in% knownDepthRange[knownDepthRange$Species == sp,c("minDepth","maxDepth")] ) {
    
    depthRange <- c(as.numeric(database.i$verbatimDepth), knownDepthRange[knownDepthRange$Species == sp,"minDepth"],knownDepthRange[knownDepthRange$Species == sp,"maxDepth"] )

    minDepthRange <- min(depthRange, na.rm=T)
    maxDepthRange <- max(depthRange, na.rm=T)
    
  }

  if( NA %in% knownDepthRange[knownDepthRange$Species == sp,c("minDepth","maxDepth")] ) {
    
    minDepthRange <- NA
    maxDepthRange <- NA
    
  }

  flags <- database.i[,"MeasurementOrFact"]
  
  # -------
  # -------

  # Land
  
  for( f in 1:nrow(database.i) ) {
    
    flagsFlatten <- fromJSON(flags[f])
    flagsFlatten[1] <- ifelse( bathymetry.i[f] > 0 , "-1" , "1" )
    flagsFlatten <- gsub("[\r\n]","", toJSON(flagsFlatten) )
    flags[f] <- flagsFlatten
    
  }
  
  # -------
  # -------
  
  # Outside Depth Range
  
  if( ! is.na(minDepthRange) | ! is.na(maxDepthRange)) {
    
      for( f in 1:nrow(database.i) ) {
        
        flagsFlatten <- fromJSON(flags[f])
        flagsFlatten[2] <- ifelse( abs(bathymetry.i[f]) > maxDepthRange | abs(bathymetry.i[f]) < minDepthRange , "-1" , "1" )
        flagsFlatten <- gsub("[\r\n]","", toJSON(flagsFlatten) )
        flags[f] <- flagsFlatten
        
      }
      
  }
  
  # -------
  # -------
  
  # Outside Distributional Range
  
  inFaoZones <- extract(faoAreas,database.i[,c("decimalLongitude","decimalLatitude")])$zone
  
  # SeaLifeBase & Literature review
  
  doesOccursInSL <- webScrapSeaLifeFAORegions(sp)
  if( nrow(doesOccursInSL) > 0 ) { doesOccursInSL <- unique(as.numeric(strsplit(doesOccursInSL$FAOAreas,",")[[1]])) }
  
  doesOccursInLt <- knownDistributionalRange[knownDistributionalRange$Species == sp,"FAORegionsOccurs"]
  if( length(doesOccursInLt) > 0 ) { doesOccursInLt <- unique(as.numeric(strsplit(doesOccursInLt,",")[[1]])) }
  
  doesOccursIn <- unique(c(unlist(doesOccursInSL),unlist(doesOccursInLt)))

  if( length(doesOccursIn) > 0 ) {
    
    for( f in 1:nrow(database.i) ) {
      
      flagsFlatten <- fromJSON(flags[f])
      flagsFlatten[3] <- ifelse(inFaoZones[f] %in% doesOccursIn , "1" , "-1" )
      flagsFlatten <- gsub("[\r\n]","", toJSON(flagsFlatten) )
      flags[f] <- flagsFlatten
      
    }
  }
  
  # Aquamaps review

  sp.id <- taxa[which(taxaSp == sp),"SPECIESID"]
  
  if( length(sp.id) != 0) { 
    
    sp.code <- nativemaps[nativemaps$SpeciesID == sp.id ,"CsquareCode"]
    sp.cells <- cells[cells$CsquareCode %in% sp.code, c("CenterLong","CenterLat")]

    distance <- sapply(1:nrow(database.i), function(x) { min(spDistsN1( as.matrix(sp.cells) , as.matrix(database.i[x,c("decimalLongitude","decimalLatitude")],ncol=2) , longlat = TRUE )) })
    distanceThreshold <- 500
    
    if( length(distance) > 0 ) {
      for( f in 1:nrow(database.i) ) {
        
        flagsFlatten <- fromJSON(flags[f])
        flagsFlatten[3] <- ifelse(distance[f] > distanceThreshold , "-1" , "1" )
        flagsFlatten <- gsub("[\r\n]","", toJSON(flagsFlatten) )
        flags[f] <- flagsFlatten
        
      }
    }
  }
  
  # IUCN review
  
  sp.id1 <- which(db1$binomial == sp)
  sp.id2 <- which(db2$binomial == sp)
  sp.id3 <- which(db3$binomial == sp)
  
  if( length(sp.id1) != 0 & length(sp.id2) != 0 & length(sp.id3) != 0 ) { 
    
    if( length(sp.id1) != 0 ) { dataToUse <- db1[sp.id1,] }
    if( length(sp.id2) != 0 ) { dataToUse <- db2[sp.id2,] }
    if( length(sp.id3) != 0 ) { dataToUse <- db3[sp.id3,] }
    
    database.i.sp <- database.i[,c("decimalLongitude","decimalLatitude")]
    coordinates(database.i.sp) <- ~Lon+Lat
    crs(database.i.sp) <- crs(dataToUse)
    
    distance <- gDistance(database.i.sp, as(dataToUse, "Spatial"),byid=TRUE)
    distanceThreshold <- 500
    
    if( length(distance) > 0 ) {
      for( f in 1:nrow(database.i) ) {
        
        flagsFlatten <- fromJSON(flags[f])
        flagsFlatten[3] <- ifelse(distance[f] > distanceThreshold , "-1" , "1" )
        flagsFlatten <- gsub("[\r\n]","", toJSON(flagsFlatten) )
        flags[f] <- flagsFlatten
        
      }
    }
  }
  
  # Literature records
  
  literatureRecords <- which(database.i$basisOfRecord == "Peer-reviewed literature")
  
  if( length(literatureRecords) > 0 ) {
    
    for( f in literatureRecords ) {
      
      flagsFlatten <- fromJSON(flags[f])
      flagsFlatten[3] <- "1"
      flagsFlatten <- gsub("[\r\n]","", toJSON(flagsFlatten) )
      flags[f] <- flagsFlatten
      
    }
  }
  
  # -------
  # -------
  
  database[which(database$acceptedName == sp),"MeasurementOrFact"] <- flags
  
}

database.Bk <- database

# -----------------
# -----------------
# Read flags

flags <- database$MeasurementOrFact
flagLand <- numeric(length(flags))
flagVerticalRange <- numeric(length(flags))
flagDistributionalRange <- numeric(length(flags))

for(i in 1:length(flags)) {
  
  flagsFlatten <- fromJSON(flags[i])
  
  if(class(flagsFlatten) != "list") {
    flagLand[i] <- as.numeric(flagsFlatten[1])
    flagVerticalRange[i] <- as.numeric(flagsFlatten[2])
    flagDistributionalRange[i] <- as.numeric(flagsFlatten[3])
  }
  if(class(flagsFlatten) == "list") {
    flagLand[i] <- ifelse(length(as.numeric(flagsFlatten$flagLand)) > 0 , as.numeric(flagsFlatten$flagLand) , 0)
    flagVerticalRange[i] <- ifelse(length(as.numeric(flagsFlatten$flagVerticalRange)) > 0 , as.numeric(flagsFlatten$flagVerticalRange) , 0)
    flagDistributionalRange[i] <- ifelse(length(as.numeric(flagsFlatten$flagDistributionalRange)) > 0 , as.numeric(flagsFlatten$flagDistributionalRange) , 0)
    
  }
}

flagedAcceptedName <- data.frame( acceptedName=unique(database$acceptedName),
                      flagDistributionalRange=sapply(unique(database$acceptedName),function(x) { sum(flagDistributionalRange[which(database$acceptedName == x)] != 0 )  }),
                      flagVerticalRange=sapply(unique(database$acceptedName),function(x) { sum(flagVerticalRange[which(database$acceptedName == x)] != 0) }))

## ---------------

acceptedNameToKeep <- flagedAcceptedName[which(flagedAcceptedName$flagDistributionalRange != 0 & flagedAcceptedName$flagVerticalRange != 0 ),"acceptedName"]
database <- database[which(database$acceptedName %in% acceptedNameToKeep),]

acceptedNameToKeep <- data.frame(acceptedName=acceptedNameToKeep)
write.csv(acceptedNameToKeep,file="../../../Data/Biodiversity Data/Cold Water Corals/acceptedNameToKeep.csv", row.names = FALSE )

## ---------------

database[which(database$year == -999),"year"] <- NA
database <- data.frame(id=1:nrow(database),database)

write.csv(database,file=paste0(dumpDirectory,"biodiversityDataBaseMergedFlagged.csv"), row.names = FALSE)

## ---------------

database <- rbind.fill(read.csv(file=paste0(dumpDirectory,"biodiversityDataBaseLocal.csv")),
                       read.csv(file=paste0(dumpDirectory,"biodiversityDataBaseOnline.csv")))

database[which( database$year == -999),"year"] <- NA
database <- data.frame(id=1:nrow( database), database)

write.csv( database,file=paste0(dumpDirectory,"biodiversityDataBaseMergedUnflagged.csv"), row.names = FALSE)

## -----------------------------------------------------
## -----------------------------------------------------
