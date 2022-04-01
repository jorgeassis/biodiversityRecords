# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
#
# ------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))

library(ggplot2)
library(jsonlite)

dumpDirectory <- "../../../Data/Biodiversity Data/Cold Water Corals/"
database <- read.csv(file=paste0(dumpDirectory,"biodiversityDataBaseMergedFlagged.csv"))

# -----------------
# Summary data

colnames(database)

nrow(database)

length(unique(database$class))
length(unique(database$genus))
length(unique(database$name))
length(unique(database$acceptedName))

range(database$year,na.rm=T)
unique(database$year)

unique(database$basisOfRecord)
sum(database$basisOfRecord == "Biodiversity Information Facility") / nrow(database)
sum(database$basisOfRecord == "Peer-reviewed literature") / nrow(database)

# -----------------
# Records per year

perYear <- aggregate(rep(1,nrow(database)), by=list(year=database$year), FUN=sum)
perYear$cumulative <- cumsum(perYear$x)

plot(perYear$year,perYear$cumulative)

range(database$year,na.rm=T)
2021-1758
hist(database$year, breaks= 263)

# ------------------------
# ------------------------
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

# -----

unique(database$class)

sum(flagLand == 1, na.rm=T) / nrow(database)
sum(flagLand == -1, na.rm=T) / nrow(database)

sum( flagLand[which(database$class == "Anthozoa")] == 1, na.rm=T) / sum(database$class == "Anthozoa")
sum( flagLand[which(database$class == "Anthozoa")] != 1, na.rm=T) / sum(database$class == "Anthozoa")

sum(flagVerticalRange == 1, na.rm=T) / length(flagVerticalRange)
sum(flagVerticalRange == -1, na.rm=T) / length(flagVerticalRange)

sum(flagDistributionalRange == 1, na.rm=T) / length(flagDistributionalRange)
sum(flagDistributionalRange == -1, na.rm=T) / length(flagDistributionalRange)

# -----------------
# Mapping

# https://github.com/exaexa/scattermore
# equal area projection , maybe "albers" in ggplot() + coord_map()

# https://rud.is/b/2015/07/24/a-path-towards-easier-map-projection-machinations-with-ggplot2/

database[ which( flagLand == -1 | flagVerticalRange == -1 | flagDistributionalRange == -1 )    ,c("decimalLongitude","decimalLatitude"    )]

database[ which( flagLand == 1 & flagVerticalRange == 1 & flagDistributionalRange == 1 )    ,c("decimalLongitude","decimalLatitude"    )]


