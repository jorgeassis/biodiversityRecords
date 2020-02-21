# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
#
# ------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
source("Dependencies.R")

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Marine Forests Database

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
# Obis and GBIF

library(gtools)
rm(ExternalDataObis)
rm(ExternalDataGbif)

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Future vulnerability of rhodolith beds under multiple environmental stressors/Data/Occurrence records/"
taxa <- read.table(file=paste0(directory,"speciesDatabase.csv"), sep = ";", stringsAsFactors = FALSE)[,1]

# ------------------------------

for( i in 1:length(taxa) ) {
  
  taxa.i <- taxa[i]
  cat( "\n" , which(taxa == taxa.i) ,"out of" , length(taxa)  )
  
  error <- TRUE
  while(error) {
    error <- FALSE
    tryCatch(  ExternalDataObis.i <- getExternalDataObis(taxa.i) , error=function(e) { error <<- TRUE } )
  }
  if( exists("ExternalDataObis") & !is.null(ExternalDataObis.i$decimalLongitude)) { ExternalDataObis <- mySmartBind(ExternalDataObis,ExternalDataObis.i) }
  if( ! exists("ExternalDataObis") & !is.null(ExternalDataObis.i$decimalLongitude) ) { ExternalDataObis <- ExternalDataObis.i }
  print(nrow(ExternalDataObis))
  
  # ---------
  
  error <- TRUE
  while(error) {
    error <- FALSE
    tryCatch(  ExternalDataGbif.i <- getExternalDataGbif(taxa.i) , error=function(e) { error <<- TRUE })
  }
  if( exists("ExternalDataGbif") & !is.null(ExternalDataGbif.i$decimalLongitude) ) { ExternalDataGbif <- mySmartBind(ExternalDataGbif,ExternalDataGbif.i) }
  if( ! exists("ExternalDataGbif") & !is.null(ExternalDataGbif.i$decimalLongitude) ) { ExternalDataGbif <- ExternalDataGbif.i }
  
}

nrow(ExternalDataObis)
nrow(ExternalDataGbif)

# Merge datasets
colnames(ExternalDataObis) ; colnames(ExternalDataGbif)
ExternalData <- mySmartBind(ExternalDataObis,ExternalDataGbif)
colnames(ExternalData)

# Lon Lat Columns 

# Save dataset
write.table(ExternalData,file=paste0(directory,"/occurrenceRecordsDB.csv"), row.names = FALSE, quote=FALSE,sep = ";",col.names = TRUE)

## ----------------------------------
