# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
#
# ------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
source("Dependencies.R")

# ------------------------------------------------------------------------------------

library(devtools)
install_github("raquamaps/aquamapsdata")

library(aquamapsdata)
library(purrr)

download_db(force = TRUE)

library(RSQLite)

filename <- "/Volumes/Jellyfish/Dropbox/Data/aquamaps.db"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver, dbname = filename)

taxa <- dbReadTable(db,"taxa")
taxa[taxa$Genus == "Octopus",]

occurrences <- dbReadTable(db,"occ")
occurrences.i <- occurrences[occurrences$SpeciesID %in% taxa[taxa$Genus == "Octopus","SPECIESID"] , ]
plot(occurrences.i[,c("CenterLong","CenterLat")])

hspen <- dbReadTable(db,"hspen")
hspen.i <- hspen[hspen$SpeciesID %in% taxa[taxa$Genus == "Octopus","SPECIESID"] , ]
head(hspen)

dbListTables(db)
temp <- dbReadTable(db,"nativemaps")
head(temp)

# Lon Lat Columns 
