# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
# One Pipeline for Modelling the distribution of marine Species
#
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# devtools::install_github('adamlilith/enmSdm')

packages.to.use <- c(  "ecospat",
                       "spatstat",
                       "enmSdm",
                       "mboost",
                       "ecodist"
                       , "raster"
                       , "verification"
                       , "gdata"
                       , "leaflet"
                       , "dismo"
                       , "gbm"
                       , "sp"
                       , "sm"
                       , "biomod2"
                       , "adehabitatHS"
                       , "SDMTools"
                       , "parallel"
                       , "doParallel"
                       , "biganalytics"
                       , "nicheROVER"
                       , "vegan"
                       , "parallel"
                       , "rgeos"
                       , "rgdal"
                       , "sdmpredictors"
                       , "usdm"
                       , "doSNOW"
                       , "ENMeval"
                       , "maptools"
                       , "rgl"
                       , "monmlp"
                       # , "rJava"
                       , "gstat")

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! "enmSdm" %in% rownames(installed.packages()) ) { devtools::install_github('adamlilith/enmSdm') }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

## -----------------------

listAllFiles <- function(directory,fileType) {
  
  files <- character()
  
  repeat {
    
    directoryTemp <- character()
    
    for( d in directory) {
      
      directoryTemp <- c(list.files(d,full.names = TRUE),directoryTemp)
      
    }
    
    files <- c(files,directoryTemp)
    directory <- directoryTemp
    
    if(length(list.files(directoryTemp[length(directoryTemp)],full.names = TRUE)) == 0) { break }
    
  }
  
  files <- files[grepl(dataLayersFileType,files)]
  return(files)
  
}

## -----------------------

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

## -----------------------

plotMap <- function(coordinates,radius,color) {
  
  set.seed(42)
  
  m <- leaflet()
  m <- addTiles(m)
  # m <- addMarkers(m, lng=coordinates[,1], lat=coordinates[,2], popup=paste0( "Species record ") , icon = greenLeafIcon)
  
  m <- addCircleMarkers(m, lng=coordinates[,1], lat=coordinates[,2], 
                        popup=paste0( "Species record ") , 
                        radius = radius, color = color , 
                        stroke = FALSE, fillOpacity = 0.5 )
  
  return(m)
  
}

## -----------------------

readRecords <- function(dataRecordsFile,separator,excludeLowRes,extraColumns,maskRecordsShapefile) {
  
  if( grepl("xls",dataRecordsFile) ) {  presences <- read.xls( dataRecordsFile , sheet=1) }
  if( grepl("csv",dataRecordsFile) ) {  presences <- read.csv( dataRecordsFile , sep=";", header = T,row.names = NULL ) }
  
  extraColumns <- c("Lon","Lat",extraColumns)
  columnsUse <- which( colnames(presences) %in% extraColumns)
  presences <- presences[,columnsUse]
  
  cat("\n")
  cat("Initial: ",nrow(presences), "occurrence records")
  cat("\n")
  
  records.to.keep <- which( ! is.na( presences$Lon ) & ! is.na( presences$Lat ) )
  presences <- presences[records.to.keep,]
  
  presences <- presences[  which( ! duplicated(presences[,c(which(colnames(presences) == "Lon"),which(colnames(presences) == "Lat"))] , fromLast = TRUE ) ) , ]
  
  if(excludeLowRes) {
    presences <- presences[which( sapply( presences$Lon , function(x) { decimalplaces(x) } ) >= 2),]
  }
  
  presences <- presences[which(presences$Lon <= 180 & presences$Lon >= -180 & presences$Lat < 90 & presences$Lat > -90 ),]
  
  cat("Final: ",nrow(presences), "occurrence records")
  cat("\n")
  cat("\n")
  
  if( !is.null(maskRecordsShapefile)) {
    
    presence.records.mask.shp <- shapefile(maskRecordsShapefile)
    crs(presence.records.mask.shp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    
    presences.pts <- presences[,1:2]
    coordinates(presences.pts) <- c("Lon", "Lat")
    crs(presences.pts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    points.not.over <- which(is.na( sp::over(presences.pts,presence.records.mask.shp) ) )
    presences <- presences[-points.not.over,]
  }
  
  return(presences)
  
}

## -----------------------

list.sdm.predictors <- function() {
  
  raster.predictors <- list_layers()
  raster.predictors <- raster.predictors[raster.predictors$dataset_code == "Bio-ORACLE",2:3]
  return(raster.predictors)
}

## -----------------------

processLayers <- function(rasters,occurrenceRecords,regionBuffer,ecoregionsLevel,minDepth,maxDepth,intertidalRegion) {
  
  set.seed(42)
  
  if( ! is.null(regionBuffer) ) {
    
    final.extent.rasters <- c( min(occurrenceRecords[,1] , na.rm=T) - regionBuffer[1],
                               max(occurrenceRecords[,1] , na.rm=T) + regionBuffer[2],
                               min(occurrenceRecords[,2] , na.rm=T) - regionBuffer[3],
                               max(occurrenceRecords[,2] , na.rm=T) + regionBuffer[4] )
    
    final.extent.rasters <- c(  ifelse(final.extent.rasters[1] < -180 , -180 , final.extent.rasters[1]),
                                ifelse(final.extent.rasters[2] > 180 , 180 , final.extent.rasters[2]),
                                ifelse(final.extent.rasters[3] < -90 , -90 , final.extent.rasters[3]),
                                ifelse(final.extent.rasters[4] > 90 , 90 , final.extent.rasters[4]) )
    
    rasters <- crop(rasters, extent(final.extent.rasters))
    
  }
  
  ## -----------------------
  
  if( !is.null(ecoregionsLevel) ) {  
    
    ecoregions <- shapefile("Dependencies/Data/Shapefiles/marine_ecoregions.shp")
    
    if( ecoregionsLevel == 1 ) {
      
      ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$RLM_CODE)
      eco.names <- names(ecoregions.shp)
      
    }
    
    if( ecoregionsLevel == 2 ) {
      
      ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$PROV_CODE)
      eco.names <- names(ecoregions.shp)
      
    }
    
    if( ecoregionsLevel == 3 ) {
      
      ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$ECO_CODE)
      eco.names <- names(ecoregions.shp)
      
    }
    
    results <- cbind(raster::extract(ecoregions.shp,occurrenceRecords[,1:2]),data.frame(occurrenceRecords))
    eco.regions.count <- aggregate(rep(1,length(results$poly.ID)), by=list(results$poly.ID), FUN=sum)
    
    adjacent.pairs <- gTouches(ecoregions.shp, byid=TRUE)
    polys <- eco.regions.count[which(eco.regions.count[,2] > 0),1]
    final.ecoregions <- numeric(0)
    
    for( poly in polys) {
      
      final.ecoregions <- c(final.ecoregions,as.numeric(eco.names[adjacent.pairs[ poly ,]]))
      
    }
    
    final.ecoregions <- unique(c(as.numeric(eco.names[ polys ]),final.ecoregions))
    
    if( ecoregionsLevel == 1 ) {
      
      ecoregions <- ecoregions[ecoregions$RLM_CODE %in% final.ecoregions,]
      
    }
    
    if( ecoregionsLevel == 2 ) {
      
      ecoregions <- ecoregions[ecoregions$PROV_CODE %in% final.ecoregions,]
      
    }
    
    if( ecoregionsLevel == 3 ) {
      
      ecoregions <- ecoregions[ecoregions$ECO_CODE %in% final.ecoregions,]
      
    }
    
    rasters <- mask(rasters, ecoregions)
    
    if(sum(is.na(raster::extract(rasters,occurrenceRecords[,1:2]))) > 0) { stop("Some records outside marine ecoregion areas")}
    
  }
  
  ## -----------------------
  
  if( intertidalRegion ) {
    
    if( ! exists("predictName")) { predictName <- "" }
    if( grepl("LGM", predictName) & ! grepl("LGM", coastLineDataLayer) ) { stop("Error: 1999") }
    
    intertidal <- raster(coastLineDataLayer)
    intertidal <- crop(intertidal, final.extent.rasters)
    
    rastersIntertidal <- mask(rasters,intertidal)
    
  }
  
  ## -----------------------
  
  if( ! is.null(minDepth) | ! is.null(maxDepth) ) {
    
    if( is.null(minDepth) ) { minDepth <- 0 }
    if( is.null(maxDepth) ) { maxDepth <- 99998 }
    
    minDepth <- minDepth * (-1)
    maxDepth <- maxDepth * (-1)
    
    rclmat <- data.frame(from=c(-99999,maxDepth,minDepth),to=c(maxDepth,minDepth,99999),value=c(NA,1,NA))
    
    bathymetry <- raster(bathymetryDataLayer)
    bathymetry <- stack(crop(bathymetry, extent(rasters)))
    
    bathymetry <- reclassify(bathymetry, rclmat)
    rastersBathymetry <- mask(rasters,bathymetry)
    
  }
  
  if( intertidalRegion & exists("rastersBathymetry") ) {
    
    rasters <- stack( sapply(1:length(names(rasters)),function(x) {  calc(stack(subset(rastersBathymetry,x),subset(rastersIntertidal,x)) , mean,na.rm=TRUE ) }) )
    names(rasters) <- names(rastersBathymetry)
    
  }
  
  if( intertidalRegion & ! exists("rastersBathymetry") ) { rasters <- rastersIntertidal }
  
  if( ! intertidalRegion & exists("rastersBathymetry") ) { rasters <- rastersBathymetry }
  
  plot(subset(rasters,1))
  
  rasters.range <- data.frame()
  
  for(i in 1:length(names(rasters))) {
    
    rasters.range.i <- data.frame( Name= names(subset(rasters,i)) ,
                                   Min=min(getValues(subset(rasters,i)),na.rm=TRUE) ,
                                   Max=max(getValues(subset(rasters,i)),na.rm=TRUE)  )
    
    rasters.range <- rbind( rasters.range , rasters.range.i )
    
  }
  
  cat( paste0("\n"))
  cat( paste0("......................................................................"))
  cat( paste0("\n"))
  cat( paste0("Range of raster processed:"))
  cat( paste0("\n"))
  cat( paste0("\n"))
  print(rasters.range)
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("......................................................................"))
  
  diff.calc <- apply(data.frame(minValue(rasters),maxValue(rasters)),1,diff)		
  if( 0 %in% diff.calc ) { stop( names(rasters)[which(diff.calc == 0)] , " has no variation")   }		
  
  return(stack(rasters))
  
}

## -----------------------

correlatedPairs <- function(rasters,threhold) { 
  
  list.of.cor <- data.frame()
  
  for( i in 1:(length(names(rasters))-1)) {
    
    for( j in (i+1):length(names(rasters))) {
      
      corr.val <- abs(cor(getValues(subset(rasters,i)),getValues(subset(rasters,j)),use = "pairwise.complete.obs"))
      
      if( corr.val >= threhold) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasters,i)),
                                                                               Var.2=names(subset(rasters,j)),
                                                                               Cor=corr.val,stringsAsFactors = FALSE)) }
      
    }
    
  }
  
  return(list.of.cor)
  
  
}

## -----------------------

relocateNACoords <- function(occurrenceRecords,relocateType,relocateSpeciesDistance,relocateSpeciesDepth,nCores) {
  
  set.seed(42)
  
  if(relocateType == "nonDistance") {
    
    vals <- raster::extract(subset(rasterLayers,1),occurrenceRecords)
    occurrenceRecords <- occurrenceRecords[which(!is.na(vals)),]
    
  }
  
  
  if(relocateType == "distance") {
    
    
    if(   relocateSpeciesDepth ) { shape <- subset(rastersRelocation,1) } 
    if( ! relocateSpeciesDepth ) { shape <- subset(rasterLayers,1) } 
    
    to.relocate <- unique(which(is.na(raster::extract(shape,occurrenceRecords))))
    coordinates.to.relocate <- occurrenceRecords[to.relocate,]
    
    old.presences <- occurrenceRecords[ (1:nrow(occurrenceRecords))[! 1:nrow(occurrenceRecords) %in% to.relocate] ,]
    
    seqListing <- round(seq(1,ncell(shape),by =1000000))
    parallelChunks <- data.frame(from = ifelse(length(seqListing[-length(seqListing)]) > 0 , seqListing[-length(seqListing)] , 1), to = c(seqListing[-c(1,length(seqListing))] , ncell(shape) ) )
    
    cl.2 <- makeCluster(nCores)
    registerDoParallel(cl.2)
    
    result <- foreach(i=1:nrow(parallelChunks), .combine=c, .verbose=FALSE, .packages=c("ncdf4","gstat","gdata","raster","data.table","bigmemory","FNN")) %dopar% { 
      
      values.temp <- shape[parallelChunks[i,1]:parallelChunks[i,2]]
      values.temp <- (parallelChunks[i,1]:parallelChunks[i,2])[which(!is.na(values.temp))]
      
      return(values.temp)
      
    }
    
    stopCluster(cl.2) ; rm(cl.2)
    
    correct.points <- xyFromCell(shape, result)
    
    if( nrow(coordinates.to.relocate) > 0 ) { 
      
      cat( paste0("Relocating ",length(to.relocate)," Points that were falling out of range"))
      cat( paste0("\n"))
      
      near.cells <- numeric(nrow(coordinates.to.relocate))
      
      for(p in 1:nrow(coordinates.to.relocate)) {
        
        near.cell.p <- spDistsN1( as.matrix(correct.points), as.matrix(coordinates.to.relocate[p,]),longlat=TRUE)
        
        if( near.cell.p[which.min(near.cell.p)] <= sqrt(sum(relocateSpeciesDistance^2,relocateSpeciesDistance^2)) ) {
          
          near.cell.p <- which.min(near.cell.p)
          
        } else {   near.cell.p <- NA }
        
        near.cells[p] <- near.cell.p
        
      }
      
      relocated <- which(!is.na(near.cells))
      
      if( length(relocated) > 0) {
        
        near.cells <- correct.points[near.cells[relocated],]
        
        colnames(near.cells) <- c("Lon","Lat")
        
        occurrenceRecords <- rbind(old.presences,near.cells)
        
      }
      
    }
    
    if( nrow(coordinates.to.relocate) == 0) { 
      
      cat( paste0("None to Relocate"))
      cat( paste0("\n"))
      
    }
    
    ## -----------------------
  }
  
  return( occurrenceRecords )
  
}

## -----------------------

correctLayer <- function(rasters,name,type,value.is,value.will) {		
  
  name.i <- which(names(rasters) == name)
  temp.r <- rasters[[name.i]]		
  
  if(type =="he") { temp.r[temp.r >= value.is] <- value.will }		
  if(type =="le") { temp.r[temp.r <= value.is] <- value.will }		
  if(type =="h") { temp.r[temp.r > value.is] <- value.will }		
  if(type =="l") { temp.r[temp.r < value.is] <- value.will }		
  
  names(temp.r) <- name
  
  rasters <- stack(subset(rasters,(1:(name.i-1))) , temp.r , subset(rasters,((name.i+1):dim(rasters)[3])) )
  
  return(rasters)		
  
}		

## -----------------------		

correlatedPairs <- function(rasters,threhold) { 
  
  list.of.cor <- data.frame()
  
  for( i in 1:(length(names(rasters))-1)) {
    
    for( j in (i+1):length(names(rasters))) {
      
      corr.val <- abs(cor(getValues(subset(rasters,i)),getValues(subset(rasters,j)),use = "pairwise.complete.obs"))
      
      if( corr.val >= threhold) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasters,i)),
                                                                               Var.2=names(subset(rasters,j)),
                                                                               Cor=corr.val,stringsAsFactors = FALSE)) }
      
    }
    
  }
  
  return(list.of.cor)
  
}

## -----------------------		

spatialAutocorrelation <- function(occurrenceRecords,rasterLayers,autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationSignif) {
  
  set.seed(42)
  
  # --------
  
  presences.t <- occurrenceRecords
  
  if(nrow(presences.t) > 1000) { presences.t <- presences.t[sample(1:nrow(presences.t),1000,replace=FALSE),]}
  
  # --------
  
  corrPairs <- correlatedPairs(rasterLayers,0.95)
  
  if( nrow(corrPairs) > 0) { rasters.i <- subset(rasterLayers,(1:length(names(rasterLayers)))[-which(names(rasterLayers) %in% corrPairs[,1])]) }
  if( nrow(corrPairs) == 0) { rasters.i <- rasterLayers }
  
  # --------
  
  presences.environment <- raster::extract(rasters.i,presences.t)
  to.remove <- unique(which(is.na(presences.environment),arr.ind = T)[,1])
  
  if(length(to.remove) > 0) {
    
    presences.t <- presences.t[-to.remove,]
    presences.environment <- raster::extract(rasterLayers,presences.t)
    
  }
  
  space <- spDists(as.matrix(presences.t),as.matrix(presences.t),longlat=TRUE)
  data <- ecodist::distance( presences.environment , method = "mahalanobis")
  
  if( is.null(autocorrelationMaxDistance) ) { autocorrelationMaxDistance <- max(space) }
  
  n.class <- round(autocorrelationMaxDistance / autocorrelationClassDistance)
  
  correlogram <- mgram( as.dist(data), as.dist(space), breaks = seq(1,autocorrelationMaxDistance,length.out = n.class) , nperm = 99 , alternative= "one.sided" )
  correlogram <- as.data.frame(correlogram$mgram) 
  
  vect.signif <- as.numeric(correlogram[,4] >= autocorrelationSignif)
  if(vect.signif[1] == 1) { vect.signif[1] <- 0}
  vect.signif.pch <- c(19,1)
  
  par(mar = c(4.5, 5.5, 4.5, 4.5) , bg = "#FFFFFF")
  plot(correlogram[,1],correlogram[,3],axes=FALSE,xlab="Distance (Km)",ylab="")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#F6F6F6")
  lines(correlogram[,1],correlogram[,3],lty=2,col="#000000",type="l",xlab="Distance (Km)",ylab="")
  
  points(correlogram[,1],correlogram[,3],pch=vect.signif.pch[vect.signif+1],col="#5B5B5B")
  axis(2,las=2,col="White",col.ticks="Black")
  axis(1,las=0,col="White",col.ticks="Black")
  box()
  title(ylab="Correlation (Mantel R)",mgp=c(4,1,0)) 
  abline(h = 0,lty=1)
  
  distance <- round( correlogram[ which(correlogram[,4] >= autocorrelationSignif)  , 1 ][1] )
  if( is.na(distance)) { distance <- autocorrelationMaxDistance }
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("First non-correlated distance: ",distance," km"))
  
  return(distance)
  
}

## -----------------------

spatialThinning <- function(occurrenceRecords,minDistance) {
  
  coordinates.t <- occurrenceRecords
  
  space <- spDists(as.matrix(coordinates.t),as.matrix(coordinates.t),longlat=TRUE)
  diag(space) <- NA
  
  reclass <- space <= minDistance
  reclass[lower.tri(reclass, diag=TRUE)] <- NA
  
  v <- colSums(reclass, na.rm=TRUE) == 0
  coordinates.t <- coordinates.t[v,]
  
  # Number of All occurrences and number to keep
  
  cat( paste0("\n"))
  cat( paste0("\n"))
  
  cat( paste0("Input Records: ",nrow(occurrenceRecords)))
  cat( paste0("\n"))
  cat( paste0("Final Records: ",nrow(coordinates.t)))
  
  # Remove from main dataset of occurrences
  
  return(coordinates.t)
  
}

## -----------------------

pseudoAbsences <- function(occurrenceRecords,rasterLayers,polygonPA,polygonPAType,paMindistance,paType,paBiasKernelSurface,paBiasKernelSurfaceProb,paBiasKernelSurfaceExport,paRatio,paEnvironmentStratification,paSpatialBalance,plotFigure) {
  
  options(warn=-1)
  
  shape <- subset(rasterLayers,1)
  
  if( ! is.null(polygonPA) ) {
    
    if( polygonPAType == "include" ) {
      shape <- mask(shape,polygonPA)
    }
    
    if( polygonPAType == "exclude" ) {
      cells <- cellFromPolygon(shape, polygonPA, weights=FALSE)
      shape[unique(unlist(cells))] <- NA
    }
    
  }
  
  if( paSpatialBalance ) {
    
    roi <- extent(shape) 
    
    if( cvType=="blocksLongitudinal" ) {
      
      # Same pa amount in extreme blocks
      
      centralPoint <- extent(shape)[1] + (max(extent(shape)[2],extent(shape)[1]) - min(extent(shape)[2],extent(shape)[1])) / 2
      
      roi.i <- roi
      roi.i[2] <- centralPoint
      cellsExtreme.1 <- as.data.frame(crop(shape,roi.i),xy=T)
      cellsExtreme.1 <- cellsExtreme.1[!is.na(cellsExtreme.1[3]),1:2]
      
      roi.i <- roi
      roi.i[1] <- centralPoint
      cellsExtreme.2 <- as.data.frame(crop(shape,roi.i),xy=T)
      cellsExtreme.2 <- cellsExtreme.2[!is.na(cellsExtreme.2[3]),1:2]
      
      if( nrow(cellsExtreme.1) < nrow(cellsExtreme.2) ) {
        resampleCellsExtreme <- (1:nrow(cellsExtreme.2))[-sample( 1:nrow(cellsExtreme.2), nrow(cellsExtreme.1),replace=FALSE )]
        shape[cellFromXY(shape,cellsExtreme.2[resampleCellsExtreme,])] <- NA
      }
      if( nrow(cellsExtreme.2) < nrow(cellsExtreme.1) ) {
        resampleCellsExtreme <- (1:nrow(cellsExtreme.1))[-sample( 1:nrow(cellsExtreme.1), nrow(cellsExtreme.2),replace=FALSE )]
        shape[cellFromXY(shape,cellsExtreme.1[ resampleCellsExtreme,])] <- NA
      }
    }
    
    if( cvType=="blocksLatitudinal" ) {
      
      # Same pa amount in extreme blocks
      
      centralPoint <- extent(shape)[3] + (max(extent(shape)[4],extent(shape)[3]) - min(extent(shape)[4],extent(shape)[3])) / 2
      roi.i <- roi
      roi.i[4] <- centralPoint
      
      cellsExtreme.1 <- as.data.frame(crop(shape,roi.i),xy=T)
      cellsExtreme.1 <- cellsExtreme.1[!is.na(cellsExtreme.1[3]),1:2]
      
      roi.i <- roi
      roi.i[3] <- centralPoint
      cellsExtreme.2 <- as.data.frame(crop(shape,roi.i),xy=T)
      cellsExtreme.2 <- cellsExtreme.2[!is.na(cellsExtreme.2[3]),1:2]
      
      if( nrow(cellsExtreme.1) < nrow(cellsExtreme.2) ) {
        resampleCellsExtreme <- (1:nrow(cellsExtreme.2))[-sample( 1:nrow(cellsExtreme.2), nrow(cellsExtreme.1),replace=FALSE )]
        shape[cellFromXY(shape,cellsExtreme.2[resampleCellsExtreme,])] <- NA
      }
      if( nrow(cellsExtreme.2) < nrow(cellsExtreme.1) ) {
        resampleCellsExtreme <- (1:nrow(cellsExtreme.1))[-sample( 1:nrow(cellsExtreme.1), nrow(cellsExtreme.2),replace=FALSE )]
        shape[cellFromXY(shape,cellsExtreme.1[ resampleCellsExtreme,])] <- NA
      }
      
    }
    
  }
  
  ## ------------
  
  nonNACells <- Which(!is.na(shape), cells=TRUE) 
  sink.points <- xyFromCell(shape, nonNACells)
  
  ## ------------
  
  if( paRatio <= 1 ) { final.paRatio <- round(nrow(occurrenceRecords) * (1/paRatio)) }
  if( paRatio > 1 ) { final.paRatio <- paRatio }
  
  ## ------------
  
  if( paBiasKernelSurface ) {
    
    # Generates bias surface kernel
    
    bias <- cellFromXY(shape, as.matrix(occurrenceRecords))
    cells <- unique(sort(bias))
    kernelXY <- xyFromCell(shape, cells)
    samps <- as.numeric(table(bias))
    KDEsur <- sm.density(kernelXY, weights=samps, display="none", hmult=2, ylim=c(extent(shape)[3]-2.5,extent(shape)[4]), xlim=c(extent(shape)[1]-2.5,extent(shape)[2]-2.5), nbins=NA)
    KDErast <- SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
    KDErast <- SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, length(KDEsur$estimate))))
    KDErast <- raster(KDErast)
    KDErast <- raster::resample(KDErast, shape, method="bilinear")
    Bias.surface <- mask(KDErast,shape) ; rm(KDErast)
    projection(Bias.surface) <- CRS("+proj=longlat +datum=WGS84")
    Bias.surface <- Bias.surface / max(getValues(Bias.surface),na.rm=TRUE)
    
    sink.points.p <- raster::extract(Bias.surface,sink.points)
    sink.points <- sink.points[which(sink.points.p > quantile(getValues(Bias.surface),probs = paBiasKernelSurfaceProb,na.rm=T)),]
    
    print(plot(Bias.surface))
    
    # paBiasKernelSurfaceExport
    
  }
  
  ## ------------
  
  if( ! is.null(paMindistance) ) {
    
    # Removes those closer than paDist
    
    sink.points.poly <- as.data.frame(occurrenceRecords)
    coordinates( sink.points.poly ) <- c( "Lon", "Lat" )
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    
    sink.points.poly <- gBuffer( sink.points.poly, width=paMindistance / 111.699, byid=TRUE )
    # plot(sink.points.poly)
    
    sink.points.pts <- as.data.frame(sink.points)
    colnames( sink.points.pts ) <- c( "Lon", "Lat" )
    coordinates( sink.points.pts ) <- c( "Lon", "Lat" )
    proj4string( sink.points.pts ) <- CRS( "+proj=longlat +datum=WGS84" )
    
    to.remove.id <- sp::over(sink.points.pts,sink.points.poly)
    to.keep <- which(is.na(to.remove.id))
    sink.points <- sink.points[to.keep,]
    
  }
  
  ## ------------
  
  if( paEnvironmentStratification ) {
    
    # Generates environmental surface for stratification
    
    to.generate <- min( c( final.paRatio , nrow(sink.points) ) )
    
    if( final.paRatio < nrow(sink.points) - (nrow(sink.points) * 0.1) ) {
      
      pa.environment <- raster::extract( rasterLayers,sink.points[,1:2] )
      pa.environment <- kmeans(pa.environment,centers=to.generate,iter.max = 10, nstart = 1)$cluster
      pa.clusters <- unique( pa.environment )
      
    }
    if( final.paRatio >= nrow(sink.points) - (nrow(sink.points) * 0.1) ) {
      
      pa.environment <- 1:nrow(sink.points)
      
    }
    
    finalKMClusters <- sapply(1:to.generate,function(x){ sample(which(pa.environment == x),1,replace=F) } )
    sink.points <- sink.points[finalKMClusters,]
    
  }
  
  if( ! paEnvironmentStratification ) {
    
    to.generate <- min( c( final.paRatio , nrow(sink.points) ) )
    sink.points <- sink.points[sample(1:nrow(sink.points),to.generate,replace=F),]
    
  }
  
  absences <- sink.points
  
  if(plotFigure) {
    
    plot(absences[,1:2],col="grey", pch=20 , cex=1 , main="Dataset : Presences and Pseudo-absences")
    points(occurrenceRecords,col="black", pch=20 , cex=1)
    
  }
  
  absences <- absences[,1:2]
  colnames(absences) <- c("Lon","Lat")
  return(absences)
  
  options(warn=0)
  
}

## -----------------------

pseudoAbsencesSimple <- function(occurrenceRecords,rasterLayers,paMindistance,paRatio,plotFigure) {
  
  shape <- subset(rasterLayers,1)
  shape[!is.na(shape)] <- 1
  shape[is.na(shape)] <- 1
  
  if( paRatio <= 1 ) { final.paRatio <- round(nrow(occurrenceRecords) * (1/paRatio)) ; pa.number <- (1/paRatio) * nrow(occurrenceRecords) }
  if( paRatio > 1 ) { final.paRatio <- paRatio ; pa.number <- final.paRatio }
  
  pa.object <- biomod2::BIOMOD_FormatingData(  resp.var = rep(1,nrow(occurrenceRecords)),
                                               expl.var = stack(rasterLayers),
                                               resp.xy = occurrenceRecords,
                                               resp.name = paste0("Species"),
                                               PA.nb.rep = 1,
                                               PA.nb.absences = pa.number,
                                               PA.strategy = "disk",
                                               PA.dist.min = paMindistance * 1000)  # plot(pa.object)
  
  sink.points <- pa.object@coord[ -(1:nrow(occurrenceRecords)) ,]
  return(sink.points)
}

## -----------------------

pseudoAbsencesSRE <- function(occurrenceRecords,rasterLayers,paMindistance,paRatio,plotFigure) {
  
  shape <- subset(rasterLayers,1)
  shape[!is.na(shape)] <- 1
  shape[is.na(shape)] <- 1
  
  if( paRatio <= 1 ) { final.paRatio <- round(nrow(occurrenceRecords) * (1/paRatio)) }
  if( paRatio > 1 ) { final.paRatio <- paRatio }
  
  
  pa.number <- (1/paRatio) * nrow(occurrenceRecords)
  
  pa.object <- biomod2::BIOMOD_FormatingData(  resp.var = rep(1,nrow(occurrenceRecords)),
                                               expl.var = stack(rasters),
                                               resp.xy = occurrenceRecords[,c("lon","lat")],
                                               resp.name = paste0("Species"),
                                               PA.nb.rep = 1,
                                               PA.nb.absences = pa.number,
                                               PA.strategy = "sre",
                                               PA.dist.min = NULL,
                                               PA.sre.quant = 0.025 )  # plot(pa.object)
  
  sink.points <- pa.object@coord[ -(1:nrow(occurrenceRecords)) ,]
  return(sink.points)
}

## -----------------------

dataPartitioning <- function(presences,absences,type,k) {
  
  if( type=="blocksLongitudinal" ) {
    
    min.Lon <- min( presences$Lon )
    max.Lon <- max( presences$Lon )
    
    bands.vect.p <- numeric(nrow(presences))
    bands.vect.a <- numeric(nrow(absences))
    
    bands <- seq(from=min.Lon,to=max.Lon,length.out = k+1)
    bands <- data.frame(lon.from=bands[-length(bands)],lat.to=bands[-1])
    
  }
  
  if( type=="blocksLatitudinal" ) {
    
    min.lat <- min( presences$Lat )
    max.lat <- max( presences$Lat )
    
    bands.vect.p <- numeric(nrow(presences))
    bands.vect.a <- numeric(nrow(absences))
    
    bands <- seq(from=min.lat,to=max.lat,length.out = k+1)
    bands <- data.frame(lat.from=bands[-length(bands)],lat.to=bands[-1])
    
  }
  
  return( bands )
  
}

## -----------------------

drawPolygon <- function(rasterLayers,occurrenceRecords,name) {
  
  x <- 0
  dataRecordsDirectory.i <- as.numeric(unlist(gregexpr("/",dataRecordsFile)))
  dataRecordsDirectory <- paste0(substr(dataRecordsFile,1,dataRecordsDirectory.i[length(dataRecordsDirectory.i)-1]),"Spatial/")
  
  if( ! dir.exists(dataRecordsDirectory) ) { dir.create(dataRecordsDirectory,recursive = TRUE) }
  
  if( paste0("Polygon.",name,".shp") %in% list.files(dataRecordsDirectory) ) {   
    
    cat("\n")
    cat("Polygon for region of interest alreay exists.")
    cat("\n")
    cat("\n")
    cat("1. Use available polygon")
    cat("\n")
    cat("2. Generate new polygon")
    cat("\n")
    while( x != 1 & x != 2) {
      
      x <- readline(":")  
      
    }
    
  }
  
  if( x == 2 | x == 0 ) {   
    
    r.1 <- subset(rasterLayers,1)
    r.1[!is.na(r.1)] <- 1
    
    plot(r.1,col="#000000")
    points(occurrenceRecords,col="red",pch=16)
    
    poly <- spatstat::clickpoly(add=TRUE)
    
    p = Polygon(cbind(poly$bdry[[1]]$x,poly$bdry[[1]]$y))
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    plot(sps,add=TRUE,col="#727272")
    
    crs(sps) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
    sps <- as(sps, "SpatialPolygons")
    rgdal::writeOGR(as(sps, "SpatialPolygonsDataFrame"), dataRecordsDirectory ,paste0("Polygon.",name), driver="ESRI Shapefile",overwrite_layer=TRUE)
    
  }
  
  if( x == 1 ) {   
    
    r.1 <- subset(rasterLayers,1)
    r.1[!is.na(r.1)] <- 1
    
    plot(r.1,col="#000000")
    
    sps = shapefile(paste0(dataRecordsDirectory,"/","Polygon.",name,".shp"))
    plot(sps,add=TRUE,col="#727272")
    points(occurrenceRecords,col="red",pch=20)
    
    sps <- as(sps, "SpatialPolygons")
    
  }
  
  return( sps )
  
}

## -----------------------

nicheModel <- function(modelType,simplifyModels,occurrenceRecords,absences,rasters,crossValidation,cvIndex,cvRemoveEdges,dataLayersMonotonocity,printResults) {
  
  ## -----------------------
  
  if( exists("model") ) { rm(model) }
  
  ## -----------------------
  
  predicted.accuracy <- NULL
  
  sorensenIndex <- function(obs,pred) {
    
    tp <- sum( obs == 1 &  pred == 1 ) / sum( obs == 1 )
    fp <- sum( obs == 0 &  pred == 1 ) / sum( obs == 0 )
    fn <- sum( obs == 1 &  pred == 0 ) / sum( obs == 1 )
    
    sorensen.i <- (2 * tp) / ( fn + (2 * tp) + fp )
    return(sorensen.i)
    
  }
  
  accuracy.area.c <- function(model,observed,predicted,predicted.map) {
    
    options(warn=-1)
    
    predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )
    predicted.accuracy <- data.frame(sedi=sediWeighted(predicted[which(observed == 1)] ,predicted[which(observed == 0)] , thresholds = seq(0, 1-0.01, by= 0.01) ),predicted.accuracy)
    
    predicted.distribution.area <- NA
    model.deviance <- NA
    aicc <- NA
    
    if( class(model)[1] == "MaxEnt") { n.param <- get.params(model) ; aicc <- calc.aicc(n.param, presences.test, predicted.map)$AICc }
    
    if( class(model)[1] == "gbm") { model.deviance <- 1 - model$cv.statistics$deviance.mean }
    
    if(cvIndex == "sedi") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$sedi),] }
    
    if(cvIndex == "tss" ) { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
    
    if(cvIndex == "deviance") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$deviance),] }
    
    if(cvIndex == "aicc") { predicted.accuracy <- predicted.accuracy[which.min(predicted.accuracy$aicc),] }
    
    if(cvIndex == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
    
    if(cvIndex == "area") { 
      
      predicted.accuracy <- predicted.accuracy[predicted.accuracy$sensitivity >= 0.95,] 
      predicted.accuracy <- predicted.accuracy[nrow(predicted.accuracy),] 
      
      if( nrow(predicted.accuracy) == 0 ) { predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )[1,] }
      
      predicted.distribution <- crop(predicted.map,extent( min( occurrenceRecords[,"Lon"] ) , max( occurrenceRecords[,"Lon"] ) , min( occurrenceRecords[,"Lat"] ) , max( occurrenceRecords[,"Lat"] )))
      
      full.area <- predicted.distribution
      full.area[!is.na(full.area)] <- 1
      full.area <- raster::resample(r.area,full.area)
      full.area <- raster::mask(full.area,predicted.distribution)
      full.area <- sum(getValues(full.area),na.rm=T)
      
      predicted.distribution[ predicted.distribution > predicted.accuracy$threshold ] <- 1
      predicted.distribution[ predicted.distribution <= predicted.accuracy$threshold ] <- NA
      predicted.distribution.area <- raster::resample(r.area,predicted.distribution)
      predicted.distribution.area <- raster::mask(predicted.distribution.area,predicted.distribution)
      predicted.distribution.area <- 1 - (sum(getValues(predicted.distribution.area),na.rm=T) / full.area)
      
    }
    
    threshold <- predicted.accuracy$threshold 
    auc <- predicted.accuracy$AUC 
    specificity <- predicted.accuracy$specificity 
    sensitivity <- predicted.accuracy$sensitivity 
    tss <- predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1 
    sedi <- predicted.accuracy$sedi
    sorensen <- sorensenIndex( observed, ifelse(predicted >= threshold , 1 , 0) )
    
    if( length(threshold) == 0 ) { threshold <- NA }
    if( length(auc) == 0 ) { auc <- NA }
    if( length(specificity) == 0 ) { specificity <- NA }
    if( length(sensitivity) == 0 ) { sensitivity <- NA }
    if( length(tss) == 0 ) { tss <- NA }
    if( length(sorensen) == 0 ) { sorensen <- NA }
    
    predicted.accuracy <- data.frame( sedi=sedi,
                                      sorensen=sorensen,
                                      threshold = threshold ,
                                      auc = auc ,
                                      specificity = specificity ,
                                      sensitivity = sensitivity ,
                                      tss = tss ,
                                      area = predicted.distribution.area ,
                                      aicc = aicc,
                                      deviance = model.deviance )
    
    options(warn=0)
    
    return(predicted.accuracy)
    
  }
  
  ## ------------------------------------------------------------------------------
  ## Clip rasters for performance
  
  extent.t <- extent( min(c(absences[,1],occurrenceRecords[,1])) - 1 , max(c(absences[,1],occurrenceRecords[,1])) + 1 , min(c(absences[,2],occurrenceRecords[,2])) - 1 , max(c(absences[,2],occurrenceRecords[,2])) + 1 )
  rasters <- crop(rasters,extent.t)
  predictors <- names(rasters)
  
  ## -----------------------
  
  cv.k.span <- 1:cvKFolds
  if( cvRemoveEdges ) { cv.k.span <- cv.k.span[c(-1,-length(cv.k.span))] }
  
  ## -----------------------
  
  if( cvIndex == "area" ) { r.area <- raster::area(subset(rasters,1)) }
  
  ## -----------------------------------------------------------------------------------
  
  cross.validation.rounds <- list()
  
  ## ----------------
  
  if( modelType == "MaxEnt" ) {
    
    stop("Error 02: Check Section")
    
    comb = expand.grid(cv.k = cv.k.span, feature.class=sapply(1:length(maxent.feature.class.span),function(x) { paste0(t(maxent.feature.class.span[x]),collapse=",")} ),betamultiplier=maxent.betamultiplier.span )
    
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat('\n')
    cat('Model: MaxEnt \n')
    cat('Cross-validation rounds:', nrow(comb))
    cat('\n')
    cat('Number of parameters:', length(maxent.feature.class.span), '+' , length(maxent.betamultiplier.span))
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    
    cv.accuracy <- foreach(c = 1:nrow(comb), .packages = c("dismo", "rJava","SDMTools","ENMeval")) %dopar% {
      
      cv <- comb[c,1]
      feature <- as.character(unlist(strsplit(as.character(comb[c,2]),split=",")))
      beta <- comb[c,3]
      args <- c(paste0("betamultiplier=",beta),feature)
      
      # ----------------------
      
      if( TRUE %in% grepl("lat", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lat <= crossValidation[cv,1] | presences$Lat >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lat <= crossValidation[cv,1] | absences$Lat >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lat >= crossValidation[cv,1] & presences$Lat <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lat >= crossValidation[cv,1] & absences$Lat <= crossValidation[cv,2] , ]
        
      }
      
      if( TRUE %in% grepl("lon", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lon <= crossValidation[cv,1] | presences$Lon >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lon <= crossValidation[cv,1] | absences$Lon >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lon >= crossValidation[cv,1] & presences$Lon <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lon >= crossValidation[cv,1] & absences$Lon <= crossValidation[cv,2] , ]
        
      }  
      
      presences.train <- presences.train[!is.na(presences.train[,1]),]
      absences.train <- absences.train[!is.na(absences.train[,1]),]
      
      presences.test <- presences.test[!is.na(presences.test[,1]),]
      absences.test <- absences.test[!is.na(absences.test[,1]),]
      
      if( nrow(presences.train) > 0 &nrow(absences.train) > 0 & nrow(presences.test) > 0 & nrow(absences.test) > 0 ) { 
        
        train.dataset <- data.frame( PA = c( rep(1,nrow(presences.train)) , rep(0,nrow(absences.train)) ) , raster::extract( rasters , rbind( presences.train, absences.train) ) )
        test.dataset <- data.frame( PA = c( rep(1,nrow(presences.test)) , rep(0,nrow(absences.test)) ) , raster::extract( rasters , rbind( presences.test, absences.test) ) )
        
        train.dataset <- train.dataset[complete.cases(train.dataset),]
        test.dataset <- test.dataset[complete.cases(test.dataset),]
        
        model <- dismo::maxent(x = rasters , p = presences.train, a=absences.train , args=c("-P","autorun=false",args) )
        
        observed <- test.dataset$PA
        predicted <- predict( model, test.dataset[,-1])
        
        if( cvIndex == "area" ) { 
          predicted.map <- predict( rasters , model )
        } else { predicted.map <- NULL}
        
        predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
        
        predicted.accuracy <- data.frame( cv.round=cv,
                                          feature=paste(feature,collapse="+"),
                                          beta=beta,
                                          predicted.accuracy)
        
        return(predicted.accuracy)
        
      }
      
    }
    
    stopCluster(cl); rm(cl) ; gc(reset=TRUE)
    
    # ------------------
    
    cv.accuracy <- do.call(rbind,cv.accuracy)
    cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold),]
    
    # ------------------
    
    if(cvIndex == "aicc") { 
      best.model <- aggregate(list(aicc=cv.accuracy$aicc), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
      best.model.i <- which.min(best.model[,"aicc"])
    }
    if(cvIndex == "tss") { 
      best.model <- aggregate(list(tss=cv.accuracy$tss), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
      best.model.i <- which.max(best.model[,"tss"])
    }
    if(cvIndex == "auc") { 
      best.model <- aggregate(list(auc=cv.accuracy$auc), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
      best.model.i <- which.max(best.model[,"auc"])
    }
    if(cvIndex == "area") { 
      best.model <- aggregate(list(area=cv.accuracy$area), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
      best.model.i <- which.max(best.model[,"area"])
    }
    
    # ------------------
    
    best.model <- best.model[best.model.i,]
    
    best.model.feature <- as.character(best.model$feature)
    best.model.beta <- as.numeric(best.model$beta)
    
    best.model.metrics <- cv.accuracy[cv.accuracy$feature == best.model.feature & cv.accuracy$beta == best.model.beta ,  ]
    best.cross.validation <- best.model.metrics[,cvIndex]
    
    model <- dismo::maxent(x = rasters , p = presences, a=absences , args=c("-P","autorun=false",c(paste0("betamultiplier=",best.model.beta),unlist(strsplit(best.model.feature, "[+]")))) )
    
    if( simplifyModels ) {
      
      model.importance <- summary.model(model,print.data=FALSE)
      model.simplift.vars <- as.character(model.importance$Variable)[model.importance$Permutation <= 1]
      
      model <- dismo::maxent(x = dropLayer(rasters,which(names(rasters) %in% model.simplift.vars)) , p = presences, a=absences , args=c("-P","autorun=false",c(paste0("betamultiplier=",best.model.beta),unlist(strsplit(best.model.feature, "[+]")))) )
      model.predicted <- predict( rasters , model )
      
    }
    
  }
  
  ## ----------------------- 
  ## -----------------------
  
  if( modelType == "BRT" ) {
    
    dataLayersMonotonocity.i <- dataLayersMonotonocity
    
    if( length(names(rasters)) == 1 ) {  
      
      dummy.raster <- rasters
      dummy.raster[1:length(rep(1:2,length(dummy.raster) / 2))] <- rep(1:2,length(dummy.raster) / 2)
      names(dummy.raster) <- "layer"
      rasters <- stack(rasters,dummy.raster)
      dataLayersMonotonocity.i <- c(dataLayersMonotonocity.i,+1)
      predictors <- c(predictors,"layer")
      
    }
    
    comb = expand.grid(cv.k = cv.k.span, learning.complex=brtLearning,tree.depth=brtTreeDepth , bag=brtBagFraction )
    
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat('\n')
    cat('Model: BRT \n')
    cat('Cross-validation rounds:', nrow(comb))
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind, .export=c("brtMaxTrees"), .packages = c("dismo","SDMTools","ENMeval","enmSdm")) %dopar% {
      
      cv <- comb[c,1]
      l.rate <- comb[c,2]
      tree.c <- comb[c,3]
      bag <- comb[c,4]
      
      presences <- as.data.frame(occurrenceRecords)
      absences <- as.data.frame(absences)
      
      if( TRUE %in% grepl("lat", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lat <= crossValidation[cv,1] | presences$Lat >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lat <= crossValidation[cv,1] | absences$Lat >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lat >= crossValidation[cv,1] & presences$Lat <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lat >= crossValidation[cv,1] & absences$Lat <= crossValidation[cv,2] , ]
        
      }
      
      if( TRUE %in% grepl("lon", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lon <= crossValidation[cv,1] | presences$Lon >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lon <= crossValidation[cv,1] | absences$Lon >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lon >= crossValidation[cv,1] & presences$Lon <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lon >= crossValidation[cv,1] & absences$Lon <= crossValidation[cv,2] , ]
        
      }  
      
      presences.train <- presences.train[!is.na(presences.train[,1]),]
      absences.train <- absences.train[!is.na(absences.train[,1]),]
      
      presences.test <- presences.test[!is.na(presences.test[,1]),]
      absences.test <- absences.test[!is.na(absences.test[,1]),]
      
      if( nrow(presences.train) > 0 & nrow(absences.train) > 0 & nrow(presences.test) > 0 & nrow(absences.test) > 0 ) { 
        
        train.dataset <- data.frame( PA = c( rep(1,nrow(presences.train)) , rep(0,nrow(absences.train)) ) , raster::extract( rasters , rbind( presences.train, absences.train) ) )
        train.dataset[train.dataset == "NaN"] <- NA
        colnames(train.dataset) <- c("PA",names(rasters))
        test.dataset <- data.frame( PA = c( rep(1,nrow(presences.test)) , rep(0,nrow(absences.test)) ) , raster::extract( rasters , rbind( presences.test, absences.test) ) )
        test.dataset[test.dataset == "NaN"] <- NA
        colnames(test.dataset) <- c("PA",names(rasters))
        
        model <- gbm.step( data=train.dataset, 
                           gbm.x = which(colnames(train.dataset) %in% predictors),
                           gbm.y = 1, 
                           family = "bernoulli", 
                           plot.main = FALSE,
                           tree.complexity = tree.c, 
                           learning.rate = l.rate, 
                           bag.fraction = bag, 
                           n.folds=10,
                           step.size=50,
                           max.trees=brtMaxTrees,
                           silent=TRUE,
                           var.monotone = dataLayersMonotonocity.i ,
                           verbose=FALSE)
        
        if( ! is.null(model) ) {
          
          num.tress <- model$gbm.call$best.trees
          observed <- test.dataset$PA
          predicted <- predict( model , test.dataset[,-1] , n.trees=num.tress,type="response")
          
          if( cvIndex == "area" ) { 
            predicted.map <- predict( rasters , model , n.trees=num.tress,type="response")
          } else { predicted.map <- NULL}
          
          predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
          
        }
        
        if( is.null(model)) {
          
          predicted.accuracy <- data.frame( sedi=0,sorensen=0,threshold=0,auc=0,specificity=0,sensitivity=0,tss=0,area=0,aicc=+Inf,deviance=0)
          
        }
        
        predicted.accuracy <- data.frame( cv.round=cv,tree.c=tree.c,l.rate=l.rate,bag=bag,predicted.accuracy)
        return(predicted.accuracy)
        
      } else { return(NULL) }
      
    }
    
    stopCluster(cl); rm(cl) ; gc(reset=TRUE)
    
    # ------------------
    
    cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
    
    if( nrow(cv.accuracy) != 0 ) { 
      
      best.model <- cv.accuracy
      best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(bag=best.model$bag,tree.c=best.model$tree.c,l.rate=best.model$l.rate), mean)
      
      # ------------------
      
      if(cvIndex == "sedi" | cvIndex == "deviance" | cvIndex == "tss" | cvIndex == "auc" | cvIndex == "area") { 
        best.model.i <- which.max(best.model[,"cvIndex"])
      }
      
      # ------------------
      
      best.model <- best.model[best.model.i,]
      
      best.model.tc <- best.model$tree.c
      best.model.lr <- best.model$l.rate
      best.model.bag <- best.model$bag
      
      best.model.metrics <- cv.accuracy[cv.accuracy$tree.c == best.model.tc & cv.accuracy$l.rate == best.model.lr & cv.accuracy$bag == best.model.bag ,  ]
      best.cross.validation <- best.model.metrics[,cvIndex]
      
      train.dataset <- data.frame( PA = c( rep(1,nrow(presences)) , rep(0,nrow(absences)) ) , raster::extract( rasters , rbind( presences, absences) ) )
      train.dataset[train.dataset == "NaN"] <- NA
      
      model <- gbm.step( data=train.dataset, 
                         gbm.x = which(colnames(train.dataset) %in% predictors),
                         gbm.y = 1, 
                         family = "bernoulli", 
                         plot.main = FALSE,
                         tree.complexity = best.model.tc, 
                         learning.rate = best.model.lr, 
                         bag.fraction = best.model.bag, 
                         n.folds=10,
                         step.size=50,
                         max.trees=brtMaxTrees,
                         silent=TRUE,
                         var.monotone = dataLayersMonotonocity.i ,
                         verbose=FALSE)
      
      if( simplifyModels ) {
        
        stop("Error 01: Check Section")
        
        model.simplify <- gbm.simplify(model, n.drops = length(predictors)-2) 
        model.simplify.dev <- which(model.simplify$deviance.summary$mean < 0) 
        
        if(length(model.simplify.dev) > 0) {
          
          how.many.to.drop <- max(model.simplify.dev)
          predictors <- names(rasters)[model.simplify$pred.list[[how.many.to.drop]] - 1]
          
          model <- gbm.step( data=train.dataset, 
                             gbm.x = which(colnames(train.dataset) %in% predictors),
                             gbm.y = 1, 
                             family = "bernoulli", 
                             plot.main = FALSE,
                             tree.complexity = best.model.tc, 
                             learning.rate = best.model.lr, 
                             bag.fraction = 0.5, 
                             n.folds=10,
                             step.size=50,
                             max.trees=1000,
                             silent=TRUE,
                             var.monotone = dataLayersMonotonocity.i ,
                             verbose=FALSE)
          
        }
      }
      
    }
    
    if( nrow(cv.accuracy) == 0 ) {
      
      model <- NULL 
      best.cross.validation <- rep(0,10)
      best.model.metrics <- data.frame(threshold = 1 ,
                                       cv <- rep(0,10) ,
                                       sedi = 0 ,
                                       auc= 0 ,
                                       specificity = 0 ,
                                       sensitivity = 0 ,
                                       tss = 0 ,
                                       aicc = 0 ,
                                       area = 1 ,
                                       deviance = 0 )
      
    }
    
  }
  
  ## -----------------------
  ## -----------------------
  
  if( modelType == "MBOOST" ) {
    
    dataLayersMonotonocity.i <- dataLayersMonotonocity
    
    comb = expand.grid(cv.k = cv.k.span, shrinkage=mboostShrinkage,df=mboostDF,mstop=mboostIterations )
    
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat('\n')
    cat('Model: MBOOST \n')
    cat('Cross-validation rounds:', nrow(comb))
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    cv.accuracy <- foreach(c = 1:nrow(comb), .combine = rbind, .packages = c("dismo","SDMTools","ENMeval","mboost","enmSdm")) %dopar% {
      
      cv <- comb[c,1]
      shrinkage <- comb[c,2]
      df <- comb[c,3]
      mstop <- comb[c,4]
      
      presences <- as.data.frame(occurrenceRecords)
      absences <- data.frame(absences)
      
      if( TRUE %in% grepl("lat", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lat <= crossValidation[cv,1] | presences$Lat >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lat <= crossValidation[cv,1] | absences$Lat >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lat >= crossValidation[cv,1] & presences$Lat <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lat >= crossValidation[cv,1] & absences$Lat <= crossValidation[cv,2] , ]
        
      }
      
      if( TRUE %in% grepl("lon", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lon <= crossValidation[cv,1] | presences$Lon >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lon <= crossValidation[cv,1] | absences$Lon >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lon >= crossValidation[cv,1] & presences$Lon <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lon >= crossValidation[cv,1] & absences$Lon <= crossValidation[cv,2] , ]
        
      }  
      
      presences.train <- presences.train[!is.na(presences.train[,1]),]
      absences.train <- absences.train[!is.na(absences.train[,1]),]
      
      presences.test <- presences.test[!is.na(presences.test[,1]),]
      absences.test <- absences.test[!is.na(absences.test[,1]),]
      
      if( nrow(presences.train) > 0 &nrow(absences.train) > 0 & nrow(presences.test) > 0 & nrow(absences.test) > 0 ) { 
        
        train.dataset <- data.frame( PA = c( rep(1,nrow(presences.train)) , rep(0,nrow(absences.train)) ) , raster::extract( rasters , rbind( presences.train, absences.train) ) )
        train.dataset[train.dataset == "NaN"] <- NA
        train.dataset <- train.dataset[complete.cases(train.dataset),]
        colnames(train.dataset) <- c("PA",names(rasters))
        
        train.dataset[,1] <- as.factor(train.dataset[,1])
        
        test.dataset <- data.frame( PA = c( rep(1,nrow(presences.test)) , rep(0,nrow(absences.test)) ) , raster::extract( rasters , rbind( presences.test, absences.test) ) )
        test.dataset[test.dataset == "NaN"] <- NA
        colnames(test.dataset) <- c("PA",names(rasters))
        
        predictors <- colnames(train.dataset)[-1]
        response <- colnames(train.dataset)[1]
        
        constr_mono <- as.character(dataLayersMonotonocity.i)
        constr_mono[constr_mono == "-1"] <- "decreasing"
        constr_mono[constr_mono == "1"] <- "increasing"
        
        rhs <- paste(c(paste("bmono(", predictors, ", constraint = \"", constr_mono,"\", df = ",df,")", sep = "")),collapse = " + ")
        fm_mono <- as.formula(paste(response, " ~ ", rhs, collapse = ""))
        
        ctrl <- boost_control(mstop = mstop, trace = FALSE, nu=shrinkage,stopintern=TRUE)
        
        tryCatch( model <- mboost(fm_mono, data = train.dataset, control = ctrl, family = Binomial(type = "adaboost" , link = "probit" )) , error=function(e) { Error <<- TRUE })
        
        if( exists("model") ) { 
          
          observed <- test.dataset$PA
          predicted <- predict(model,test.dataset)
          
          if( cvIndex == "area" ) { 
            predicted.map <- predict(rasters,model)
            predicted.map <- predicted.map + ( min(getValues(predicted.map),na.rm=T) * (-1))
            predicted.map <- predicted.map / max(getValues(predicted.map),na.rm=T)
          } else { predicted.map <- NULL}
          
          predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
          
        }
        
        if( ! exists("model") ) { 
          
          predicted.accuracy <- data.frame( threshold=NA, auc=NA, specificity=NA, sensitivity=NA, tss=NA, area=NA, aicc=NA, deviance=NA )
          
        }
        
        
        predicted.accuracy <- data.frame( cv.round=cv,
                                          shrinkage=shrinkage,
                                          df=df,
                                          n.iterations=mstop,
                                          predicted.accuracy)
        
        return(predicted.accuracy)
        
      }
      
    }
    
    stopCluster(cl); rm(cl) ; gc(reset=TRUE)
    
    # ------------------
    
    best.model <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
    best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(shrinkage=best.model$shrinkage,df=best.model$df,n.iterations=best.model$n.iterations), mean)
    
    # ------------------
    
    if(cvIndex == "sedi" | cvIndex == "deviance" | cvIndex == "tss" | cvIndex == "auc" | cvIndex == "area") { 
      best.model.i <- best.model[sort(best.model$cvIndex, decreasing = TRUE, index.return = TRUE)$ix,]
    }
    
    # ------------------
    
    model <- NULL
    m.i <- 0
    
    while( is.null(model) ){
      
      m.i <- m.i + 1
      
      best.model.shrinkage <- best.model.i[m.i,]$shrinkage
      best.model.df <- best.model.i[m.i,]$df
      best.model.iterations <- best.model.i[m.i,]$n.iterations
      
      best.model.metrics <- cv.accuracy[cv.accuracy$shrinkage == best.model.shrinkage & cv.accuracy$df == best.model.df & cv.accuracy$n.iterations == best.model.iterations ,  ]
      best.cross.validation <- best.model.metrics[,cvIndex]
      
      train.dataset <- data.frame( PA = c( rep(1,nrow(presences)) , rep(0,nrow(absences)) ) , raster::extract( rasters , rbind( presences, absences) ) )
      colnames(train.dataset) <- c("PA",names(rasters))
      train.dataset <- train.dataset[complete.cases(train.dataset),]
      
      train.dataset[,1] <- as.factor(train.dataset[,1])
      
      predictors <- colnames(train.dataset)[-1]
      response <- colnames(train.dataset)[1]
      
      constr_mono <- as.character(dataLayersMonotonocity.i)
      constr_mono[constr_mono == "-1"] <- "decreasing"
      constr_mono[constr_mono == "1"] <- "increasing"
      
      rhs <- paste(c(paste("bmono(", predictors, ", constraint = \"", constr_mono,"\", df = ",best.model.df,")", sep = "")),collapse = " + ")
      fm_mono <- as.formula(paste(response, " ~ ", rhs, collapse = ""))
      
      ctrl <- boost_control(mstop = best.model.iterations, trace = FALSE, nu=best.model.shrinkage,stopintern=TRUE)
      
      tryCatch( model <- mboost(fm_mono, data = train.dataset, control = ctrl, family = Binomial(type = "adaboost" , link = "probit" )) , error=function(e) { Error <- TRUE })
      
    }
    
    if( simplifyModels ) {
      
      stop("Needs programming")
      
    }
    
  }
  
  ## -----------------------------------------------
  ## -----------------------------------------------
  
  
  if( modelType == "MPNN" ) {
    
    dataLayersMonotonocity.i <- dataLayersMonotonocity
    
    comb = expand.grid(cv.k = cv.k.span, hidden=mpnnHidden,itereractions=mpnnItereractions )
    
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat('\n')
    cat('Model: MPNN \n')
    cat('Cross-validation rounds:', nrow(comb))
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind, .packages = c("monmlp","dismo","SDMTools","ENMeval","enmSdm")) %dopar% {
      
      cv <- comb[c,1]
      hidden <- comb[c,2]
      itereractions <- comb[c,3]
      
      presences <- as.data.frame(occurrenceRecords)
      absences <- as.data.frame(absences)
      
      if( TRUE %in% grepl("lat", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lat <= crossValidation[cv,1] | presences$Lat >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lat <= crossValidation[cv,1] | absences$Lat >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lat >= crossValidation[cv,1] & presences$Lat <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lat >= crossValidation[cv,1] & absences$Lat <= crossValidation[cv,2] , ]
        
      }
      
      if( TRUE %in% grepl("lon", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lon <= crossValidation[cv,1] | presences$Lon >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lon <= crossValidation[cv,1] | absences$Lon >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lon >= crossValidation[cv,1] & presences$Lon <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lon >= crossValidation[cv,1] & absences$Lon <= crossValidation[cv,2] , ]
        
      }  
      
      presences.train <- presences.train[!is.na(presences.train[,1]),]
      absences.train <- absences.train[!is.na(absences.train[,1]),]
      
      presences.test <- presences.test[!is.na(presences.test[,1]),]
      absences.test <- absences.test[!is.na(absences.test[,1]),]
      
      if( nrow(presences.train) > 0 & nrow(absences.train) > 0 & nrow(presences.test) > 0 & nrow(absences.test) > 0 ) { 
        
        train.dataset <- data.frame( PA = c( rep(1,nrow(presences.train)) , rep(0,nrow(absences.train)) ) , raster::extract( rasters , rbind( presences.train, absences.train) ) )
        train.dataset[train.dataset == "NaN"] <- NA
        colnames(train.dataset) <- c("PA",names(rasters))
        test.dataset <- data.frame( PA = c( rep(1,nrow(presences.test)) , rep(0,nrow(absences.test)) ) , raster::extract( rasters , rbind( presences.test, absences.test) ) )
        test.dataset[test.dataset == "NaN"] <- NA
        colnames(test.dataset) <- c("PA",names(rasters))
        
        correctLayers <- which(dataLayersMonotonocity.i == -1)
        
        for(correctLayers.i in correctLayers) {
          train.dataset[,correctLayers.i] <- train.dataset[,correctLayers.i] * (-1)
          test.dataset[,correctLayers.i] <- test.dataset[,correctLayers.i] * (-1)
        }
        
        model <- monmlp.fit(x = as.matrix(train.dataset[,-1]), y = as.matrix(train.dataset[,1]),  monotone = 1:(ncol(train.dataset)-1),  hidden1 = hidden, bag = TRUE, iter.max = itereractions, iter.stopped = 10)
        
        if( ! is.null(model) ) {
          
          observed <- test.dataset$PA
          predicted <- as.numeric(monmlp.predict(as.matrix(test.dataset[,-1]), model))
          
          if( cvIndex == "area" ) { 
            predicted.map <- predictDistribution(rasters,model, reclassToOne=FALSE)
          } else { predicted.map <- NULL}
          
          predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
          
        }
        
        if( is.null(model)) {
          
          predicted.accuracy <- data.frame( sedi=0,sorensen=0,threshold=0,auc=0,specificity=0,sensitivity=0,tss=0,area=0,aicc=+Inf,deviance=0)
          
        }
        
        predicted.accuracy <- data.frame( cv.round=cv,hidden=hidden,itereractions=itereractions,predicted.accuracy)
        return(predicted.accuracy)
        
      } else { return(NULL) }
      
    }
    
    stopCluster(cl); rm(cl) ; gc(reset=TRUE)
    
    # ------------------
    
    cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
    
    if( nrow(cv.accuracy) != 0 ) { 
      
      best.model <- cv.accuracy
      best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(hidden=best.model$hidden,itereractions=best.model$itereractions), mean)
      
      # ------------------
      
      if(cvIndex == "sedi" | cvIndex == "deviance" | cvIndex == "tss" | cvIndex == "auc" | cvIndex == "area") { 
        best.model.i <- which.max(best.model[,"cvIndex"])
      }
      
      # ------------------
      
      best.model <- best.model[best.model.i,]
      
      best.model.hidden <- best.model$hidden
      best.model.itereractions <- best.model$itereractions
      
      best.model.metrics <- cv.accuracy[cv.accuracy$hidden == best.model.hidden & cv.accuracy$itereractions == best.model.itereractions  ,  ]
      best.cross.validation <- best.model.metrics[,cvIndex]
      
      train.dataset <- data.frame( PA = c( rep(1,nrow(presences)) , rep(0,nrow(absences)) ) , raster::extract( rasters , rbind( presences, absences) ) )
      train.dataset[train.dataset == "NaN"] <- NA
      
      correctLayers <- which(dataLayersMonotonocity.i == -1)
      
      for(correctLayers.i in correctLayers) {
        train.dataset[,correctLayers.i] <- train.dataset[,correctLayers.i] * (-1)
      }
      
      model <- monmlp.fit(x = as.matrix(train.dataset[,-1]), y = as.matrix(train.dataset[,1]),  monotone = 1:(ncol(train.dataset)-1),  hidden1 = best.model.hidden, bag = TRUE, iter.max = best.model.itereractions, iter.stopped = 10)
      
      if( simplifyModels ) {
        
        stop("Error 01: Check Section")
        
      }
      
    }
    
    if( nrow(cv.accuracy) == 0 ) {
      
      model <- NULL 
      best.cross.validation <- rep(0,10)
      best.model.metrics <- data.frame(threshold = 1 ,
                                       cv <- rep(0,10) ,
                                       sedi = 0 ,
                                       auc= 0 ,
                                       specificity = 0 ,
                                       sensitivity = 0 ,
                                       tss = 0 ,
                                       aicc = 0 ,
                                       area = 1 ,
                                       deviance = 0 )
      
    }
    
  }
  
  ## -----------------------
  ## -----------------------
  
  if( printResults ) { 
    
    cat( paste0("\n"))
    cat( paste0("\n --------------------------------------------------------"))
    cat( paste0("\n"))
    cat( paste0("\n"))
    cat( paste0("Accuracy of model:"))
    cat( paste0("\n"))
    cat( paste0("SEDI: ",round( mean(best.model.metrics[,"sedi"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"sedi"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("AUC: ",round( mean(best.model.metrics[,"auc"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"auc"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("specificity: ",round( mean(best.model.metrics[,"specificity"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"specificity"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("sensitivity: ",round( mean(best.model.metrics[,"sensitivity"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"sensitivity"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("tss: ",round( mean(best.model.metrics[,"tss"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"tss"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("aicc: ",round( mean(best.model.metrics[,"aicc"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"aicc"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("deviance: ",round( mean(best.model.metrics[,"deviance"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"deviance"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("area: ",round( mean(best.model.metrics[,"area"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"area"]) , digits = 4)) )
    cat( paste0("\n"))
    
  }
  
  ## -----------------------------------------------
  
  best.model.metrics <- best.model.metrics[best.model.metrics$threshold != 0 ,]
  
  return( list(model=model, 
               cv=best.cross.validation,
               sedi=mean(best.model.metrics[,"sedi"]),
               auc=mean(best.model.metrics[,"auc"]),
               Specificity=mean(best.model.metrics[,"specificity"]),
               Sensitivity=mean(best.model.metrics[,"sensitivity"]),
               tss=mean(best.model.metrics[,"tss"]),
               area=mean(best.model.metrics[,"area"]),
               deviance=mean(best.model.metrics[,"deviance"]))
  )
  
}

## -----------------------

summaryModel <- function(model,print.data) {
  
  if( class(model) == "MaxEnt" ) { 
    
    data.to.plot <- ENMeval::var.importance(model)
    colnames(data.to.plot) <- c("Variable","Percentage","Permutation")
    rownames(data.to.plot) <- NULL
    
    data.to.plot <- data.to.plot[sort(data.to.plot$Percentage,index.return=T)$ix,]
    
    barplot(  t(data.to.plot[,2:3]), 
              main="Variable contribution", 
              horiz=TRUE,
              legend = TRUE,
              beside = TRUE, 
              names.arg=as.character(data.to.plot$Variable))
    
    if(print.data) { print(data.to.plot) }
    
    
  }
  
  if( class(model) == "gbm" ) { 
    
    data.to.plot <- summary(model)
    colnames(data.to.plot) <- c("Variable","Percentage")
    rownames(data.to.plot) <- NULL
    
    if(print.data) { print(data.to.plot) }
    
  }
  
  if( class(model) == "mboost" ) { 
    
    data.to.plot <- data.frame(mboost::varimp(model))
    data.to.plot <- data.frame(Variable=data.to.plot$variable,Percentage=(data.to.plot$reduction / sum(data.to.plot$reduction) ) * 100 )
    
    if(print.data) { print(data.to.plot) }
    
  }
  
  return(data.to.plot)
  
}

## -----------------------

modelPlot <- function(model,rasters,distribution,variable.to.plot,print.limiting,auto.limiting,distribution.threshold,distribution.predicted,val.limiting,export.data.to.plot,plot.x.lim,plot.y.lim) {
  
  if( missing(plot.x.lim) ) { plot.x.lim <- NULL }
  if( missing(plot.y.lim) ) { plot.y.lim <- NULL }
  if( missing(distribution.threshold) ) { distribution.threshold <- NULL }
  if( missing(distribution.predicted) ) { distribution.predicted <- NULL }
  
  options(warn=-1)
  
  ## -----------------------
  
  if( class(model) == "MaxEnt" ) { 
    
    names.on.model <- names(rasters)
    model.predictor <- names.on.model[variable.to.plot]
    
    if(variable.to.plot > length(names.on.model) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(names.on.model) ) { 
      
      data.to.plot <- data.frame(Variable=data.frame(response(model,variable.to.plot,expand=0))[,1] , Effect= data.frame(response(model,variable.to.plot,expand=0))[,2] )
      
      data.to.plot <- data.to.plot[ data.to.plot[,1] >= min(getValues(subset(rasters,variable.to.plot)),na.rm=T) &
                                      data.to.plot[,1] <= max(getValues(subset(rasters,variable.to.plot)),na.rm=T) , ]
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
      variable.name <- names.on.model[variable.to.plot]
      xlim <- plot.x.lim
      ylim <- plot.y.lim
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "gbm" ) { 
    
    model.predictor <- model$var.names[variable.to.plot]
    
    logit2prob <- function(logit){
      odds <- exp(logit)
      prob <- odds / (1 + odds)
      return(prob)
    }
    
    if(variable.to.plot > length(model$gbm.call$predictor.names) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(model$gbm.call$predictor.names) ) { 
      
      data.to.plot <- data.frame(Variable=plot(model,variable.to.plot, return.grid=TRUE)[,1] , Effect=scale(plot(model,variable.to.plot, return.grid=TRUE)[,2],scale=FALSE))
      
      if( distribution != "gaussian" ) {  data.to.plot[,2] <- logit2prob(data.to.plot[,2]) }
      
      variable.name <- model$var.names[variable.to.plot]
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "mboost" ) { 
    
    names.on.model <- names(rasters)
    model.predictor <- names.on.model[variable.to.plot]
    
    if(variable.to.plot > length(names.on.model) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(names.on.model) ) { 
      
      data.to.plot <- getValues(rasters)
      data.to.plot <- data.to.plot[!is.na(data.to.plot[,1]),]
      data.to.plot <- data.to.plot[sample(1:nrow(data.to.plot), ifelse( nrow(data.to.plot) > 100 , 100 , nrow(data.to.plot)  ) ,replace=FALSE),]
      
      data.to.plot.x <- data.to.plot[,variable.to.plot]
      
      means <- colMeans(matrix(data.to.plot[,-variable.to.plot],ncol=length(names(rasters) )-1))
      
      for(m in 1:length(means) ) {
        
        data.to.plot[,(1:length(names(rasters)))[-variable.to.plot][m]] <- means[m]
        
      }
      
      data.to.plot <- data.frame(Variable=data.to.plot.x , Effect= predict(model,as.data.frame(data.to.plot)) )
      data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
      variable.name <- names.on.model[variable.to.plot] 
      xlim <- plot.x.lim
      ylim <- plot.y.lim
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "list" ) { 
    
    model.predictor <- names(rasterLayers)[variable.to.plot]
    
    logit2prob <- function(logit){
      odds <- exp(logit)
      prob <- odds / (1 + odds)
      return(prob)
    }
    
    data.to.plot <- as.data.frame(rasterLayers)
    data.to.plot <- data.to.plot[complete.cases(data.to.plot),]
    data.to.plot <- data.to.plot[sample(1:nrow(data.to.plot),2500),]
    
    correctLayers <- which(dataLayersMonotonocity.i == -1)
    
    for(correctLayers.i in correctLayers) {
      data.to.plot[,correctLayers.i] <- data.to.plot[,correctLayers.i] * (-1)
    }
    
    for( var in (1:(length(names(rasterLayers))))[-variable.to.plot] ) {
      data.to.plot[,var] <- mean(data.to.plot[,var])
    }
    
    predicted <- monmlp.predict( as.matrix(data.to.plot), model)
    data.to.plot <- data.frame(Variable = as.numeric(data.to.plot[,variable.to.plot]) , Effect=as.numeric(predicted) )
    
    if(variable.to.plot == correctLayers) {
      data.to.plot[,"Variable"] <- data.to.plot[,"Variable"] * (-1)
    }
    
    data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
    
    if( distribution != "gaussian" ) {  data.to.plot[,2] <- logit2prob(data.to.plot[,2]) }
    
    variable.name <- model.predictor
    
    if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
      plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
    }
    
    if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
      plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
    }
    
  }
  
  ## -----------------------
  
  if( length(unique(plot.y.lim)) < 2 ) { plot.y.lim <- c(plot.y.lim[1]-0.1,plot.y.lim[1]+0.1) }
  
  ## -----------------------
  
  lab <- ifelse( exists("contribution"),paste0(variable.name," (",round(contribution,digits=2),"%)"),variable.name)
  
  par(mar = c(4.5, 5.5, 4.5, 4.5))
  plot(data.to.plot,ylim=plot.y.lim,xlim=plot.x.lim,lty=1,col="#8E8E8E",type="l",ylab="",xlab=lab,axes=FALSE)
  axis(2,las=2,col="White",col.ticks="Black")
  axis(1,las=0,col="White",col.ticks="Black")
  box()
  title(ylab="Effect on Response",mgp=c(4,1,0)) 
  
  span.vals <- data.to.plot[,2]
  
  if( length(unique(span.vals)) >= 2 ) { 
    
    span.vals.var <- numeric(length(span.vals)-1)
    for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- (abs(span.vals[i+1]) - abs(span.vals[i])) / (max(span.vals)-min(span.vals)) } 
    span.vals.var <- abs(span.vals.var)
    
    if( data.to.plot[1,2] > data.to.plot[nrow(data.to.plot),2]  ) { tail <- "Max" } else { tail <- "Min"}
    
    if( tail == "Min") {  
      
      tipping.point <- data.to.plot[which(span.vals.var >= 0.05)[1],1]
      
    }
    
    if( tail == "Max" ) {   
      
      tipping.point <- data.to.plot[which(span.vals.var >= 0.05)[length(which(span.vals.var >= 0.05))],1] 
      
    }
  }
  
  if( length(unique(span.vals)) < 2 ) { print.limiting <- FALSE ; tipping.point <- 0 }
  
  if( ! auto.limiting ) { tipping.point <- val.limiting }
  
  if( ! is.null(distribution.threshold)) {
    
    distribution.predicted[distribution.predicted >= distribution.threshold] <- 1
    distribution.predicted[distribution.predicted < distribution.threshold] <- NA
    distribution.predicted <- distribution.predicted * subset(rasters,variable.to.plot)
    
    tipping.point <- range(getValues(distribution.predicted),na.rm=T)[which.min(abs(range(getValues(distribution.predicted),na.rm=T) - tipping.point))]
    
  }
  
  
  cat( paste0("\n"))
  cat( paste0("Predictor: ",model.predictor ))
  cat( paste0("\n"))
  
  if( print.limiting ) {
    
    cat( paste0("Limiting point: ",tipping.point))
    cat( paste0("\n"))
    abline(v = tipping.point , lty=2 , col="#AC4F4F")
    
  }
  
  options(warn=0)
  
  if(export.data.to.plot) { return(data.to.plot)}
  
}

## -----------------------

predictDistribution <- function(rasters,model,reclassToOne) {
  
  if( class(model)[1] == "MaxEnt" ) {
    
    predicted.distribution <- predict( model , rasters )
    predicted.distribution <- predicted.distribution / max(getValues(predicted.distribution),na.rm=TRUE)
    
  }   
  
  if( class(model)[1] == "gbm" ) {
    
    num.tress <- model$gbm.call$best.trees
    if(is.null(num.tress)) { num.tress <- length(model$trees) }
    predicted.distribution <- predict( rasters , model , n.trees=num.tress,type="response")
    
  }
  
  if( class(model)[1] == "mboost" ) {
    
    options(warn=-1)
    predicted.distribution <- predict( rasters , model)
    options(warn=0)
    
  }
  
  if( class(model)[1] == "list" ) {
    
    shape <- subset(rasters,1)
    shape[!is.na(shape)] <- 1
    cellsPredict <- Which(!is.na(shape),cells=TRUE)
    cellsPredictValues <- as.matrix(rasters[cellsPredict])
    
    correctLayers <- which(dataLayersMonotonocity == -1)
    
    for(correctLayers.i in correctLayers) {
      cellsPredictValues[,correctLayers.i] <- cellsPredictValues[,correctLayers.i] * (-1)
    }
    
    predicted.distribution <- monmlp.predict( cellsPredictValues , model)
    shape[cellsPredict] <- as.numeric(predicted.distribution)
    predicted.distribution <- shape
    
  }
  
  if(reclassToOne) {
    predicted.distribution <- predicted.distribution + ( min(getValues(predicted.distribution),na.rm=T) * (-1))
    predicted.distribution <- predicted.distribution / max(getValues(predicted.distribution),na.rm=T)
  }
  
  return(predicted.distribution)
  
}

## -----------------------

accuracyPredicted <- function(predicted.distribution,presences,absences,type) {
  
  observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
  predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
  
  predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )
  
  predicted.accuracy <- data.frame(predicted.accuracy,sedi = sediWeighted(predicted[which(observed == 1)] ,predicted[which(observed == 0)] , thresholds = seq(0, 1-0.01, by = 0.01)))
  
  if(type == "sedi") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$sedi),] }
  if(type == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
  if(type == "tss") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
  if(type == "Kappa") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$Kappa),] }
  
  predicted.accuracy <- data.frame( threshold =  predicted.accuracy$threshold ,
                                    sedi = predicted.accuracy$sedi ,
                                    auc = predicted.accuracy$AUC ,
                                    specificity = predicted.accuracy$specificity ,
                                    sensitivity = predicted.accuracy$sensitivity ,
                                    tss = predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1 )
  
  return(predicted.accuracy)
  
}

## -----------------------

reclassifyPredicted <- function(predicted.distribution,presences,absences,method,reclassThreshold) {
  
  if( missing(presences) ) { presences <- NULL  }
  if( missing(absences) ) { absences <- NULL  }
  
  if(method == "directReclass") {
    
    
    predicted.distribution[predicted.distribution > reclassThreshold] <- 1
    predicted.distribution[predicted.distribution <= reclassThreshold] <- 0
    
  }
  
  if(method == "maxTSS") {
    
    observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
    predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
    
    predicted.accuracy <- accuracy( observed , predicted , threshold =100 )
    predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),]
    
    predicted.distribution[predicted.distribution > predicted.accuracy$threshold] <- 1
    predicted.distribution[predicted.distribution <= predicted.accuracy$threshold] <- 0
    
  }
  
  if(method == "minAREA") {
    
    observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
    predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
    
    predicted.accuracy <- accuracy( observed , predicted , threshold = 1000 )
    predicted.accuracy <- predicted.accuracy[predicted.accuracy$sensitivity > reclass.threshold,]
    predicted.accuracy <- predicted.accuracy[nrow(predicted.accuracy),]
    
    predicted.distribution[predicted.distribution > predicted.accuracy$threshold] <- 1
    predicted.distribution[predicted.distribution <= predicted.accuracy$threshold] <- 0
    
  }
  return(predicted.distribution)
}

## -----------------------

generate.points.from.prediction <- function(presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,threshold) {
  
  set.seed(42)
  
  # Best reclassification
  
  predictive.model.sp.1.v <- unique(getValues(predictive.model.sp.1))
  predictive.model.sp.2.v <- unique(getValues(predictive.model.sp.2))
  
  if( sum(!is.na(predictive.model.sp.1.v)) > 2 ) { stop("Wrong reclass file type") }
  if( sum(!is.na(predictive.model.sp.2.v)) > 2 ) { stop("Wrong reclass file type") }
  
  if( class(predictive.model.sp.1)[1] != "RasterLayer") { raster.predict.1 <- raster(predictive.model.sp.1) }
  if( class(predictive.model.sp.1)[1] == "RasterLayer") { raster.predict.1 <- predictive.model.sp.1 }
  
  predict.1 <- as.data.frame(raster.predict.1,xy=TRUE)
  presences.1 <- predict.1[which(predict.1[,3] == 1),1:2]
  
  if( class(predictive.model.sp.2)[1] != "RasterLayer") { raster.predict.2 <- raster(predictive.model.sp.2) }
  if( class(predictive.model.sp.2)[1] == "RasterLayer") { raster.predict.2 <- predictive.model.sp.2 }
  
  predict.2 <- as.data.frame(raster.predict.2,xy=TRUE)
  presences.2 <- predict.2[which(predict.2[,3] == 1),1:2]
  colnames(presences.2) <- colnames(presences.1)
  
  max.points <- min(nrow(presences.1),nrow(presences.2))
  
  set.seed(42)
  presences.1 <- cbind(presences.1[ sample(1:nrow(presences.1),max.points,replace=F) ,],rep(1, max.points )) ; colnames(presences.1) <- c("Lon","Lat","Sp")
  
  set.seed(42)
  presences.2 <- cbind(presences.2[ sample(1:nrow(presences.2),max.points,replace=F) ,],rep(2, max.points )) ; colnames(presences.2) <- c("Lon","Lat","Sp")
  
  m <- leaflet()
  m <- addTiles(m)
  m <- addCircleMarkers(m, lng=c(presences.1[,1],presences.2[,1]), lat=c(presences.1[,2],presences.2[,2]), popup=paste0( "Species record ") , radius = 3, color = c( rep("Red",nrow(presences.1)) , rep("Blue",nrow(presences.2))) )
  print(m)
  
  return(rbind(presences.1,presences.2))
  
}

## -----------------------

pca.ordination <- function(rasters.1,rasters.2,rnd.points.sp) {
  
  rnd.points.env.sp.1 <- as.data.frame(rasters.1,na.rm=TRUE)
  rnd.points.env.sp.1 <- rnd.points.env.sp.1[sample(1:nrow(rnd.points.env.sp.1),1000,replace=FALSE),]
  
  rnd.points.env.sp.2 <- as.data.frame(rasters.2,na.rm=TRUE)
  rnd.points.env.sp.2 <- rnd.points.env.sp.2[sample(1:nrow(rnd.points.env.sp.2),1000,replace=FALSE),]
  
  rnd.points.sp.1.env <- raster::extract(rasters.1,rnd.points.sp[rnd.points.sp$Sp==1,1:2])
  rnd.points.sp.1.env <- rnd.points.sp.1.env[complete.cases(rnd.points.sp.1.env),]
  rnd.points.sp.2.env <- raster::extract(rasters.2,rnd.points.sp[rnd.points.sp$Sp==2,1:2])
  rnd.points.sp.2.env <- rnd.points.sp.2.env[complete.cases(rnd.points.sp.2.env),]
  
  env.used <- rbind(rnd.points.env.sp.1,rnd.points.env.sp.2,rnd.points.sp.1.env,rnd.points.sp.2.env)
  # env.used <- rbind(rnd.points.sp.1.env,rnd.points.sp.2.env)
  pca.env.used <- dudi.pca(env.used,scannf=F,nf=2)
  ecospat.plot.contrib(contrib=pca.env.used$co, eigen=pca.env.used$eig)
  
  sp.vect <- c( rep("env.sp1",nrow(rnd.points.env.sp.1)) , 
                rep("env.sp2",nrow(rnd.points.env.sp.2)) , 
                rep("sp1",nrow(rnd.points.sp.1.env)) , 
                rep("sp2",nrow(rnd.points.sp.2.env)) )
  
  
  # PCA scores for the whole study area
  scores.globclim <- pca.env.used$li
  # PCA scores for the species native distribution
  scores.sp.1 <- suprow(pca.env.used,rnd.points.sp.1.env)$li
  scores.sp.2 <- suprow(pca.env.used,rnd.points.sp.2.env)$li
  
  scores.clim.sp.1 <- suprow(pca.env.used,rnd.points.env.sp.1)$li
  scores.clim.sp.2 <- suprow(pca.env.used,rnd.points.env.sp.2)$li
  
  grid.clim.sp.1 <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.sp.1,
                                          sp=scores.sp.1, R=100,
                                          th.sp=0)
  grid.clim.sp.2 <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.sp.2,
                                          sp=scores.sp.2, R=100,
                                          th.sp=0)
  
  newList <- list("grid.clim.sp.1" = grid.clim.sp.1, "grid.clim.sp.2" = grid.clim.sp.2)
  
  return(newList)
  
}

## -----------------------

prepare.data.niche.overlap <- function(rasters.1,rasters.2,presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,type) {
  
  if(type == 1) {
    
    rnd.points.sp <- data.frame( rbind(presences.sp.1,presences.sp.2) , Sp = c(rep(1,nrow(presences.sp.1)) , rep(2,nrow(presences.sp.2))))
    
  }
  
  if(type == 2) {
    
    rnd.points.sp <- generate.points.from.prediction(presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,0.95)
    
  }
  
  environment.used.by.species.1 <- raster::extract(rasters.1,rnd.points.sp[rnd.points.sp$Sp==1,1:2])
  environment.used.by.species.1 <- environment.used.by.species.1[complete.cases(environment.used.by.species.1),]
  environment.used.by.species.2 <- raster::extract(rasters.2,rnd.points.sp[rnd.points.sp$Sp==2,1:2])
  environment.used.by.species.2 <- environment.used.by.species.2[complete.cases(environment.used.by.species.2),]
  
  data <- data.frame(  species = c(rep(1,nrow(environment.used.by.species.1)),rep(2,nrow(environment.used.by.species.2))) ,
                       rbind(environment.used.by.species.1,environment.used.by.species.2) )
  
  print(aggregate(data[2:ncol(data)], data[1], mean))
  return(data)
  
}

## -----------------------

niche.overlap.plot <- function(data,nsamples) {
  
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  
  clrs <- c("black", "red")  # colors for each species
  nsamples <- 10
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  species.par.data <- tapply(1:nrow(data), data$species, function(ii) X = data[ii, 2:ncol(data)])
  niche.plot(niche.par = species.par, ndens=1000, niche.data = species.par.data, pfrac = 0, col = clrs, xlab = expression("Niche overlap"))
  
}

## -----------------------

mean.niche.overlap.plot <- function(data,nsamples,colors.sp) {
  
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  
  # Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher
  # accuracy.  the variable over.stat can be supplied directly to the
  # overlap.plot function
  
  over.stat <- nicheROVER::overlap(species.par, nreps = nsamples, nprob = 10000, alpha = c(0.95, 0.99))
  
  cat( paste0("\n"))
  cat( paste0("............................................................."))
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("The mean overlap metrics calculated across iteratations for both niche region sizes"))
  cat( paste0("\n"))
  
  over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
  print(round(over.mean[, ,  1], 2))
  
  cat( paste0("\n"))
  cat( paste0("............................................................."))
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("Credible interval"))
  cat( paste0("\n"))
  
  over.cred <- apply(over.stat * 100, c(1:2, 4), quantile, prob = c(0.025, 0.975), na.rm = TRUE)
  print(round(over.cred[, , , 1]),2)  # display alpha = .95 niche region
  
  # Overlap plot.Before you run this, make sure that you have chosen your
  # alpha ecoregionsLevel.
  
  clrs <- colors.sp  # colors for each species
  over.stat <- nicheROVER::overlap(species.par, nreps = nsamples, nprob = 1000, alpha = 0.95)
  overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE, xlab = "Overlap Probability (%) -- Niche Region Size: 95%")
  
}

## -----------------------

prepare.data.niche.overlap.2 <- function(rasters.1,rasters.2,presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,type=1) {
  
  
  if(type == 1) {
    
    rnd.points.sp <- data.frame( rbind(presences.sp.1,presences.sp.2) , Sp = c(rep(1,nrow(presences.sp.1)) , rep(2,nrow(presences.sp.2))))
    
  }
  
  if(type == 2) {
    
    rnd.points.sp <- generate.points.from.prediction(presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,0.95)
    
  }
  
  grid.clim <- pca.ordination(rasters.1,rasters.2,rnd.points.sp)
  
  return(grid.clim)
  
}

## -----------------------

Dsquared <- function(model = NULL, 
                     obs = NULL, 
                     pred = NULL, 
                     family = NULL, # needed only when 'model' not provided
                     adjust = FALSE, 
                     npar = NULL) { # needed only when 'model' not provided
  # version 1.4 (31 Aug 2015)
  
  model.provided <- ifelse(is.null(model), FALSE, TRUE)
  
  if (model.provided) {
    if (!("glm" %in% class(model))) stop ("'model' must be of class 'glm'.")
    if (!is.null(pred)) message("Argument 'pred' ignored in favour of 'model'.")
    if (!is.null(obs)) message("Argument 'obs' ignored in favour of 'model'.")
    obs <- model$y
    pred <- model$fitted.values
    
  } else { # if model not provided
    if (is.null(obs) | is.null(pred)) stop ("You must provide either 'obs' and 'pred', or a 'model' object of class 'glm'.")
    if (length(obs) != length(pred)) stop ("'obs' and 'pred' must be of the same length (and in the same order).")
    if (is.null(family)) stop ("With 'obs' and 'pred' arguments (rather than a model object), you must also specify one of two model family options: 'binomial' or 'poisson' (in quotes).")
    else if (!is.character(family)) stop ("Argument 'family' must be provided as character (i.e. in quotes: 'binomial' or 'poisson').")
    else if (length(family) != 1 | !(family %in% c("binomial", "poisson"))) stop ("'family' must be either 'binomial' or 'poisson' (in quotes).")
    
    if (family == "binomial") {
      if (any(!(obs %in% c(0, 1)) | pred < 0 | pred > 1)) stop ("'binomial' family implies that 'obs' data should be binary (with values 0 or 1) and 'pred' data should be bounded between 0 and 1.")
      link <- log(pred / (1 - pred))  # logit
    }  # end if binomial
    
    else if (family == "poisson") {
      if (any(obs %%1 != 0)) stop ("'poisson' family implies that 'obs' data should consist of whole numbers.")
      link <- log(pred)
    }  # end if poisson
    
    model <- glm(obs ~ link, family = family)
  }  # end if model not provided
  
  D2 <- (model$null.deviance - model$deviance) / model$null.deviance
  
  if (adjust) {
    if (model.provided) {
      n <- length(model$y)
      #p <- length(model$coefficients)
      p <- attributes(logLik(model))$df
    } else {
      if (is.null(npar)) stop ("Adjusted D-squared from 'obs' and 'pred' values (rather than a model object) requires specifying the number of parameters in the underlying model ('npar').")
      n <- length(na.omit(obs))
      p <- npar
    }  # end if model.provided else
    
    D2 <- 1 - ((n - 1) / (n - p)) * (1 - D2)
  }  # end if adjust
  
  return (D2)
}



dataTranformation <- function(vector,new.range) {
  
  exp.lm <- data.frame( input = sort(unique(vector)),
                        outup = seq(0,1,length.out = length(unique(vector))))
  
  fit <- lm(outup ~ input, data=exp.lm)
  
  new.data <- data.frame(input=vector)
  new.data <- unlist(predict.lm(fit,new.data ))
  new.data <- as.vector(new.data)
  new.data[which(new.data == min(new.data))] <- 0
  new.data[which(new.data == max(new.data))] <- 1
  
  return(new.data)
  
}

# ---------------------------------------------------------------------------------------------------------------

getTaxonomyNewSpecies <- function(taxa,rank) {
  
  library(worms)
  library(worrms)
  
  if( rank != "species" ) {
    
    error <- FALSE
    
    tryCatch({ list.of.species <- wormsbynames(taxon_names=taxa)$AphiaID }, error = function(e) { error <<- TRUE })
    
    if( ! error ) {     
      
      final.list <- data.frame()
      higher.taxa <- TRUE
      
      while(higher.taxa){
        
        list.of.childen <- data.frame()
        
        for(sp in 1:length(list.of.species)) {
          
          for(offset in seq(1,1000000,50)) {
            
            if( exists("list.of.childen.t") ) {  rm(list.of.childen.t)  }
            
            tryCatch( list.of.childen.t <- wm_children(id=list.of.species[sp],offset=offset,marine_only=FALSE) , error=function(e) { error <- TRUE })
            
            if( ! exists("list.of.childen.t") ) {  break  }
            
            if( exists("list.of.childen.t") ) { 
              
              list.of.childen <- rbind(list.of.childen,list.of.childen.t)
              rm(list.of.childen.t) 
              
            }
            
          }
          
        }
        
        list.of.species <- list.of.childen$AphiaID
        list.of.rank <- list.of.childen$rank
        
        to.include <- which(list.of.childen$rank == "Species")
        
        if( "Species" %in% list.of.rank ) { 
          
          final.list <- rbind(final.list,list.of.childen[to.include,] )
          
        }
        
        if( "Species" %in% list.of.rank & length(unique(list.of.rank)) == 1  ) { higher.taxa <- FALSE }
        if( "Species" %in% list.of.rank & length(unique(list.of.rank)) > 1  ) { higher.taxa <- TRUE ; list.of.species <- list.of.species[-to.include] }
        
        # wormsbyid(list.of.species)$scientificname
        
      }
      
      taxa <- final.list$AphiaID
      
    }
    
    if( error ) { taxa <- NULL }
    
  }
  
  # --------------------------
  
  if( rank == "species" ) {
    
    error <- FALSE
    
    tryCatch({ taxa <- wormsbynames(taxon_names=taxa)$AphiaID }, error = function(e) { error <<- TRUE })
    
    if( error ) { taxa <- NULL }
    
  }
  
  # ------------------ 
  
  if( !is.null(taxa)) {
    
    taxa <- unique(taxa)
    
    new.species <- data.frame()
    
    for( t in 1:length(taxa) ) {
      
      taxa.to.search <- taxa[t]
      new.species <- rbind(new.species,wormsbyid(taxa.to.search))
      
    }
    
    return(new.species)
    
  }
  
  # ------------------ 
  
  if( is.null(taxa)) {
    
    return(NULL)
    
  }
}

# ---------------------------------------------------------------------------------------------------------------

getExternalDataObis <- function(taxa) {
  
  if( missing(taxa)) { errormessage("no taxa (Worms name) introduced.") }
  
  library(robis)
  
  for( i in 1:length(taxa) ) { 
    
    print(i)
    
    taxa.i <- taxa[i]
    error <- FALSE
    
    tryCatch( my_occs_obis <- occurrence(scientificname = taxa.i ) , error=function(e) { error <<- TRUE })
    
    if( ! error ) { if( nrow(my_occs_obis) == 0 ) { my_occs_obis <- data.frame() } }
    
    if( error ) { my_occs_obis <- data.frame() }
    
    if( nrow(my_occs_obis) > 0) {
      
      my_occs_obis <- subset(my_occs_obis, my_occs_obis$decimalLongitude !=0 & my_occs_obis$decimalLatitude !=0)
      
    }
    
    if( nrow(my_occs_obis) > 0) {
      
      my_occs_obisInfo <- my_occs_obis$dataset_id
      my_occs_obisInfo <- unique(my_occs_obis$dataset_id)
      
      for(z in 1:length(my_occs_obisInfo) ) {
        
        error <- TRUE
        errortrials <- 0
        
        while(error & errortrials < 10) {
          error <- FALSE
          errortrials <- errortrials + 1
          tryCatch(  z.Res <- RJSONIO::fromJSON(paste0("https://api.obis.org/v3/dataset/",my_occs_obisInfo[z])) , error=function(e) { error <- TRUE })
        }
        
        if(!error) {  
          
          institutionCode <- my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"institutionCode"]
          collectionCode <- my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"collectionCode"]
          
          my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"accessRights"] <- z.Res$results[[1]]$intellectualrights
          
          z.Res <- paste0( z.Res$results[[1]]$citation,
                           " ",
                           ifelse(!is.na(institutionCode) | !is.null(institutionCode) , institutionCode , ""),
                           " ",
                           ifelse(!is.na(collectionCode) | !is.null(collectionCode) , collectionCode , ""),
                           " (Available: Ocean Biogeographic Information System. Intergovernmental Oceanographic Commission of UNESCO. www.iobis.org. Accessed: ", Sys.time())
          
          my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"bibliographicCitation"] <- z.Res
          
        }
        
      }
      
    }
    
  }
  
  return(my_occs_obis)
  
}

# ---------------------------------------------------------------------------------------------------------------

getExternalDataGbif <- function(taxa) {
  
  if( missing(taxa)) { errormessage("no taxa (Worms name) introduced.") }
  
  library(rgbif)
  
  nRecords <- gbif(strsplit(as.character(taxa.i), " ")[[1]][1], strsplit(as.character(taxa.i), " ")[[1]][2], geo=T, removeZeros=T , download=FALSE )
  
  if( nRecords > 300 ) {
    
    seqListing <- seq(0,nRecords,by =300)
    if(max(seqListing) < nRecords) { seqListing <- c(seqListing,nRecords) }
    parallelChunks <- data.frame(from = seqListing[-length(seqListing)], to = c(seqListing[-c(1,length(seqListing))] -1 , nRecords ) )
    
    my_occs_gbif <- data.frame()
    
    for(ch in 1:nrow(parallelChunks)) { 
      
      tmpfile <- paste(tempfile(), ".json", sep = "")
      download.file(paste0("http://api.gbif.org/v1/occurrence/search?scientificname=",strsplit(as.character(taxa.i), " ")[[1]][1],"+",strsplit(as.character(taxa.i), " ")[[1]][2],"&offset=",parallelChunks[ch,1],"&limit=300"), tmpfile, quiet = TRUE, method="curl") 
      json <- scan(tmpfile, what = "character", quiet = TRUE, sep = "\n", encoding = "UTF-8")
      json <- chartr("\a\v", "  ", json)
      x <- jsonlite::fromJSON(json)
      r <- x$results
      r <- r[, !sapply(r, class) %in% c("data.frame", "list")]
      rownames(r) <- NULL
      my_occs_gbif <- smartbind(my_occs_gbif,r)
      
    }
    
    if (length(my_occs_gbif) == 0) {
      my_occs_gbif <- data.frame()
    }
    if (length(my_occs_gbif) == 1) {
      my_occs_gbif <- my_occs_gbif[[1]]
    }
    
  }
  
  if( nRecords <= 300 ) {
    
    tmpfile <- paste(tempfile(), ".json", sep = "")
    test <- try(download.file(paste0("http://api.gbif.org/v1/occurrence/search?scientificname=",strsplit(as.character(taxa.i), " ")[[1]][1],"+",strsplit(as.character(taxa.i), " ")[[1]][2],"&offset=",0,"&limit=300"), tmpfile, quiet = TRUE))
    
    json <- scan(tmpfile, what = "character", quiet = TRUE, sep = "\n", encoding = "UTF-8")
    json <- chartr("\a\v", "  ", json)
    x <- jsonlite::fromJSON(json)
    r <- x$results
    
    if( length(r) > 0) {
      
      r <- r[, !sapply(r, class) %in% c("data.frame", "list")]
      rownames(r) <- NULL
      my_occs_gbif <- r 
      
    }
    
  }
  
  if( exists("my_occs_gbif") ) { if( is.null(my_occs_gbif) ) { my_occs_gbif <- data.frame() } }
  
  if( ! exists("my_occs_gbif") ) { my_occs_gbif <- data.frame() }
  
  if( ! is.null(my_occs_gbif$decimalLatitude) & nrow(my_occs_gbif) > 0 ) {
    
    my_occs_gbif <- subset(my_occs_gbif, decimalLatitude !=0 & decimalLongitude !=0)
    
  }
  
  if( nrow(my_occs_gbif) > 0 ) {
    
    my_occs_gbif_all <- unique(my_occs_gbif$datasetKey)
    
    for(z in 1:length(my_occs_gbif_all) ) {
      
      z.Res <- gbif_citation(x=my_occs_gbif_all[z])
      
      my_occs_gbif[my_occs_gbif$datasetKey == my_occs_gbif_all[z] ,"accessRights"] <- ifelse(!is.null(z.Res$rights),z.Res$rights,"")
      
      z.Res <- z.Res$citation$citation
      
      my_occs_gbif[my_occs_gbif$datasetKey == my_occs_gbif_all[z],"bibliographicCitation"] <- z.Res
      
    }
    
  }
  
  return(my_occs_gbif)
  
}

# ---------------------------------------------------------------------------------------------------------------

mySmartBind <- function(dfr1,dfr2) {
  
  common_cols <- intersect(colnames(dfr1), colnames(dfr2))
  common_cols <- rbind(
    subset(dfr1, select = common_cols), 
    subset(dfr2, select = common_cols)
  )
  return(common_cols)
}

