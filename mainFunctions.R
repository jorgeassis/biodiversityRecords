
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------

library(stringr)
library(stringi)
library(plyr)
library(maptools)
library(raster)
library(sf)
library(ggplot2)
library(rgbif)
library(dismo)
library(rinat)
require(robis)
library(RJSONIO)
library(parallel)
library(doParallel)
library(RSQLite)
library(rgeos)
library(gdata)
library(readxl)
# library(rfishbase)
library(rnaturalearth)

# ---------------------------------------------------------------------------------------------------------------

countWords <- function(x) { str_count(x, '\\w+') }

# ---------------------------------------------------------------------------------------------------------------
