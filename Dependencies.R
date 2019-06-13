
## ------------------------------------------------------------------------------------------------------------------

library(readr)


packages.to.use <- c(
  "raster",
  "rgdal",
  "readr",
  "parallel",
  "doParallel",
  "devtools",
  "gdata",
  "leaflet",
  "data.table",
  "worms",
  "RCurl",
  "RJSONIO",
  "RSQLite",
  "taxize",
  "mapview",
  "htmlwidgets",
  "dismo",
  "maptools",
  "worrms",
  "robis",
  "fulltext",
  "devtools",
  "dplyr",
  "rcrossref",
  "DT",
  "shiny",
  "geosphere",
  "rotl",
  "worrms",
  "spatstat")

for(package in packages.to.use) {
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  # if( ! package %in% rownames(installed.packages()) & package == "rJava" ) { install.packages("rJava","http://rforge.net/",type="source") }
  
  if( package == "robis" & ! "robis" %in% rownames(installed.packages()) ) { devtools::install_github("iobis/robis") ; library(robis) }
  if( package == "worrms" & ! "worrms" %in% rownames(installed.packages()) ) { devtools::install_github("ropensci/worrms") ; library(worrms) }
  
  if( ! package %in% rownames(installed.packages()) ) { sink() ; stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

## ------------------------------------------------------------------------------------------------------------------

login.control <- function() {
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  dbExecute(sql,paste0("UPDATE RDB_Users SET LastLogDate='",Sys.time(),"' WHERE fid = ", UserIDSession)) 
  dbDisconnect(sql)
  
}

## ------------------------------------------------------------------------------------------------------------------

pop.taxa.information <- function(taxa) {
  
  main.taxa.table <- get.data.from.table("Taxa","fid, SpeciesName, Description, DescriptionSource, DescriptionSourceLink, Distribution, DistributionSource, CommonName","SpeciesName",taxa)
  
  for(t in 1:length(main.taxa.table$SpeciesName)) {
    
    cat('\014')
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------------------')
    cat('\n')
    cat('Getting information for species', t , ' out of ',length(main.taxa.table$SpeciesName))
    cat('\n')
    
    try( results <- get.taxa.description(main.taxa.table$SpeciesName[t]) , silent = TRUE)
    
    if(exists("results")) {
      
      fid.to.change <- main.taxa.table$fid[t]
      cat("\n Species: ",main.taxa.table$SpeciesName[t]," (fid:",fid.to.change,")")
      Description=gsub("'", "", results$text) 
      DescriptionSource=gsub("'", "", results$owner) 
      DescriptionSourceLink=gsub("'", "", results$owner.link) 
      Distribution=gsub("'", "", results$distribution) 
      DistributionSource=gsub("'", "", results$sources) 
      CommonName=gsub("'", "", results$en.common.names) 
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbExecute(sql,paste0("UPDATE RDB_Taxa SET Description='",Description,"' WHERE fid = ", fid.to.change))
      dbExecute(sql,paste0("UPDATE RDB_Taxa SET DescriptionSource='",DescriptionSource,"' WHERE fid = ", fid.to.change))
      dbExecute(sql,paste0("UPDATE RDB_Taxa SET DescriptionSourceLink='",DescriptionSourceLink,"' WHERE fid = ", fid.to.change))
      dbExecute(sql,paste0("UPDATE RDB_Taxa SET Distribution='",Distribution,"' WHERE fid = ", fid.to.change))
      dbExecute(sql,paste0("UPDATE RDB_Taxa SET DistributionSource='",DistributionSource,"' WHERE fid = ", fid.to.change))
      dbExecute(sql,paste0("UPDATE RDB_Taxa SET CommonName='",CommonName,"' WHERE fid = ", fid.to.change))
      dbDisconnect(sql)
      
      rm(results)
      
    }
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

get.taxa.description <- function(taxa) {
  
  try( id.sp <- wm_name2id(name=taxa) , silent = TRUE)
  
  if( exists("id.sp") ) {
    
    try( distribution <- wm_distribution(id=id.sp) , silent = TRUE)
    if( exists("distribution") ) { distribution <- paste(unique(distribution$higherGeography),collapse = "; ") } else { distribution <- "" }
    
    try( en.common.names <- wm_common_id(id=id.sp) , silent = TRUE)
    if( exists("en.common.names") ) { en.common.names <- paste(en.common.names$vernacular[en.common.names$language_code == "eng"],collapse = "; ") } else { en.common.names <- "" }
    
    try( sources <- wm_sources(id=id.sp) , silent = TRUE)
    if( exists("sources") ) { sources <- sources$reference[sources$use == "basis of record"] } else { sources <- "" }
    
  }
  
  link <- eol_search(terms=taxa)$link[1]
  
  if( ! is.null(link) ) {
    
    if( ! "InternalProc" %in% list.files("./Data/") ) { dir.create(file.path("Data/InternalProc")) }
    
    try(download.file(link,destfile=paste0("Data/InternalProc/main.search.txt"),method="auto",quiet=TRUE), silent = TRUE)
    main.parsed.file <- readLines(paste0( "Data/InternalProc/main.search.txt"))
    file.remove(list.files("./Data/InternalProc" , pattern=paste0( "main.search.txt") , full.names=TRUE))
    
    line.in.parse.start <- which(sapply(1:length(main.parsed.file),function(x) { gregexpr("<h1 class='scientific_name'>",main.parsed.file[x])[[1]][1] }) != -1)[1]
    line.in.parse.end <- (line.in.parse.start:length(main.parsed.file))[which(sapply(line.in.parse.start:length(main.parsed.file),function(x) { gregexpr("</i>",main.parsed.file[x])[[1]][1] }) != -1)[1]]
    species.name <- character()
    for(f in line.in.parse.start:line.in.parse.end) { species.name <- paste(species.name,main.parsed.file[f],sep="") }
    species.name <- gsub("<h1 class='scientific_name'>", "", species.name) 
    species.name <- gsub("<i>", "", species.name) 
    species.name <- gsub("</i>", "", species.name) 
    
    line.in.parse.start <- which(sapply(1:length(main.parsed.file),function(x) { gregexpr("<h4>Description</h4>",main.parsed.file[x])[[1]][1] }) != -1)[1]
    
    if(is.na(line.in.parse.start)) { link <- NULL } 
    
    if(!is.na(line.in.parse.start)) {
      
      line.in.parse.end <- (line.in.parse.start:length(main.parsed.file))[which(sapply(line.in.parse.start:length(main.parsed.file),function(x) { gregexpr("</div>",main.parsed.file[x])[[1]][1] }) != -1)[1]]
      text.string <- character()
      for(f in line.in.parse.start:line.in.parse.end) { text.string <- paste(text.string,main.parsed.file[f],sep="") }
      text.string <- gsub("<h4>Description</h4>", "", text.string) 
      text.string <- gsub("</div>", "", text.string) 
      text.string <- gsub("&nbsp;", "", text.string) 
      
      line.in.parse.start <- (line.in.parse.start:length(main.parsed.file))[which(sapply(line.in.parse.start:length(main.parsed.file),function(x) { gregexpr("<p class='owner'>",main.parsed.file[x])[[1]][1] }) != -1)[1]]
      line.in.parse.end <- (line.in.parse.start:length(main.parsed.file))[which(sapply(line.in.parse.start:length(main.parsed.file),function(x) { gregexpr("</p>",main.parsed.file[x])[[1]][1] }) != -1)[1]]
      text.string.owner <- character()
      for(f in line.in.parse.start:line.in.parse.end) { text.string.owner <- paste(text.string.owner,main.parsed.file[f],sep="") }
      text.string.owner <- gsub("<p class='owner'>", "", text.string.owner) 
      text.string.owner <- gsub("</p>", "", text.string.owner) 
      text.string.owner <- gsub("&copy; &nbsp;", "", text.string.owner) 
      
      line.in.parse.start <- (line.in.parse.start:length(main.parsed.file))[which(sapply(line.in.parse.start:length(main.parsed.file),function(x) { gregexpr("Source:",main.parsed.file[x])[[1]][1] }) != -1)[1]]
      line.in.parse.end <- (line.in.parse.start:length(main.parsed.file))[which(sapply(line.in.parse.start:length(main.parsed.file),function(x) { gregexpr("</p>",main.parsed.file[x])[[1]][1] }) != -1)[1]]
      text.string.owner.link <- character()
      for(f in line.in.parse.start:line.in.parse.end) { text.string.owner.link <- paste(text.string.owner.link,main.parsed.file[f],sep="") }
      text.string.owner.link <- gsub("Source:", "", text.string.owner.link) 
      text.string.owner.link <- gsub("<br>", "", text.string.owner.link) 
      text.string.owner.link <- gsub("</p>", "", text.string.owner.link) 
      text.string.owner.link <- substr(text.string.owner.link, unlist(unlist(gregexpr('"',text.string.owner.link)))[1]+1, unlist(gregexpr('">',text.string.owner.link))-1)
      
    }
    
  }
  
  if( is.null(link) ) { species.name <- "" }
  
  if( species.name != taxa ) {  text.string <- "" ; text.string.owner <- "" ; text.string.owner.link <- ""  }
  
  return(list(text=text.string,owner=text.string.owner,owner.link=text.string.owner.link,distribution=distribution,sources=sources,en.common.names=en.common.names))
  
}

## ------------------------------------------------------------------------------------------------------------------

get.fids.by.geopolygon <- function(fid) {
  
  options(warn=-1)
  
  queryied.data <- get.data.from.table("Record","fid, SpeciesID, Lon, Lat, CoordinateType, Depth, DateYear, SourceRefType , SourceRefID","fid",fid)
  queryied.data <- queryied.data[!is.na(queryied.data$Lon),]
  
  min.x <- min(queryied.data$Lon,na.rm=TRUE) - 5
  max.x <- max(queryied.data$Lon,na.rm=TRUE) + 5
  min.y <- min(queryied.data$Lat,na.rm=TRUE) - 5
  max.y <- max(queryied.data$Lat,na.rm=TRUE) + 5
  
  world.t <- crop(coast.line.lr,extent(min.x,max.x,min.y,max.y))
  
  coords <- cbind(queryied.data$Lon,queryied.data$Lat)
  plot(world.t,col="#C0C0C0",border=NA)
  points(coords,col="black",pch=20)
  
  poly <- clickpoly(add=TRUE)
  
  p = Polygon(cbind(poly$bdry[[1]]$x,poly$bdry[[1]]$y))
  ps = Polygons(list(p),1)
  sps = SpatialPolygons(list(ps))
  plot(sps,add=TRUE,col="#727272")
  
  coords <- as.data.frame(coords)
  colnames(coords) <- c("x", "y")
  coordinates(coords) <- c("x", "y")
  proj4string(coords) <- proj4string(sps)
  
  over.poly <- !is.na(over(coords, as(sps, "SpatialPolygons")))
  points(coords[over.poly],col="red",pch=20)
  
  return(queryied.data$fid[over.poly])
  
}

## ------------------------------------------------------------------------------------------------------------------

updated.version.control <- function() {
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  get.version <<- dbReadTable(sql, "Version")
  dbDisconnect(sql)
  
  updated.version.control <- data.frame(  Created = as.character(get.version$Created) ,
                                          Created.by.user = as.character(get.version$Created.by.user) , 
                                          Updated = as.character(Sys.time()) , 
                                          Last.update.by.user = UserIDSession , stringsAsFactors = FALSE )
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  dbWriteTable(sql, "Version", updated.version.control , append=FALSE, overwrite=TRUE )
  dbDisconnect(sql)
  
}

## ------------------------------------------------------------------------------------------------------------------

internet.available <- function() {
  
  sink("/dev/null") 
  tryCatch( i.test <- is.character(getURL("www.google.com"))
            , error=function(e) { } )
  sink()
  
  if(! exists("i.test") ) { errormessage("A internet connection is required") }
  
}

## ------------------------------------------------------------------------------------------------------------------

errormessage <- function(message) {
  
  if(missing(message)) { message <- "Error type 0" }
  
  cat("\f") ; cat("\n") ; cat("Error:",message) ; exit()
  
}

## ------------------------------------------------------------------------------------------------------------------

get.data.from.table <- function(type,field,query.field,query) {
  
  options(warn=-1)
  
  if(missing(type)) { errormessage("Table name requiered") }
  
  if(type == "Taxa") { table <- "RDB_Taxa" }
  if(type == "Literature") { table <- "RDB_Ref_Literature" }
  if(type == "User") { table <- "RDB_Users" }
  if(type == "Repository") { table <- "RDB_Ref_Repository" }
  if(type == "Record") { table <- "MDB_GeoData" }
  if(type == "Version") { table <- "Version" }
  if(type == "Abundance") { table <- "RDB_Abundance" }
  if(type == "Collection") { table <- "RDB_Sample_Collection" }
  
  if(missing(field)) { field <- "*" }
  if(missing(query.field)) { query.field <- NULL }
  if(missing(query)) { query <- NULL }
  
  if( class(query) == "data.frame" ) { query <- as.numeric(unlist(query))  }
  
  if( is.null(query.field) & ! is.null(query) ) { errormessage("query.field name requiered") }
  if( ! is.null(query.field) & is.null(query) ) { errormessage("query name requiered") }
  
  if( is.null(query.field) | is.null(query) ) {
    
    sql.parse.query <- paste0("SELECT ",field," FROM ",table)
    
  }
  
  if( ! is.null(query.field) & ! is.null(query) ) {
    
    query <- paste0("'", query, "'", collapse=", ")
    query <- paste0("(", query, ")")
    
    sql.parse.query <- paste0("SELECT ",field," FROM " , table , " WHERE " , query.field , " IN " , query,"")
    
  }
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  temp.data <- dbGetQuery(sql, sql.parse.query )
  dbDisconnect(sql)
  
  return(temp.data)
  
}

## ------------------------------------------------------------------------------------------------------------------

get.last.fid <- function(type) {
  
  if(missing(type)) { errormessage("Table name requiered") }
  
  if(type == "Taxa") { table <- "RDB_Taxa" }
  if(type == "Literature") { table <- "RDB_Ref_Literature" }
  if(type == "User") { table <- "RDB_Users" }
  if(type == "Abundance") { table <- "RDB_Abundance" }
  if(type == "Collection") { table <- "RDB_Sample_Collection" }
  if(type == "Repository") { table <- "RDB_Ref_Repository" }
  if(type == "Record") { table <- "MDB_GeoData" }
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  
  if( ! table %in% dbListTables(sql) ) { errormessage("table name not in Data Base") }
  
  last.fid <- dbGetQuery(sql,paste0("SELECT max(fid) FROM ", table)) 
  
  dbDisconnect(sql)
  
  last.fid <- as.numeric(unlist(last.fid))
  
  if( last.fid == -Inf | is.na(last.fid) ) { last.fid <- 0 }
  return(last.fid)
  
}

## ------------------------------------------------------------------------------------------------------------------

get.table.structure <- function(type) {
  
  if(missing(type)) { errormessage("Table name requiered") }
  
  if(type == "Taxa") { table <- "RDB_Taxa" }
  if(type == "Literature") { table <- "RDB_Ref_Literature" }
  if(type == "User") { table <- "RDB_Users" }
  if(type == "Repository") { table <- "RDB_Ref_Repository" }
  if(type == "Record") { table <- "MDB_GeoData" }
  if(type == "Abundance") { table <- "RDB_Abundance" }
  if(type == "Collection") { table <- "RDB_Sample_Collection" }
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  temp.data <- dbListFields(sql, table)
  dbDisconnect(sql)
  
  return(temp.data)
  
}

## ------------------------------------------------------------------------------------------------------------------

maintainence.sql <- function(pop.taxa.information,resolve.references,resolve.duplicated.references,resolve.coordinates,resolve.unused.references,resolve.species,resolve.duplicated.records,pop.missing.coordinates,backup.db) {
  
  if(missing(pop.taxa.information)) { pop.taxa.information <- FALSE }
  if(missing(resolve.references)) { resolve.references <- FALSE }
  if(missing(resolve.duplicated.references)) { resolve.duplicated.references <- FALSE }
  if(missing(resolve.unused.references)) { resolve.unused.references <- FALSE }
  if(missing(resolve.species)) { resolve.species <- FALSE }
  if(missing(resolve.coordinates)) { resolve.coordinates <- FALSE }
  if(missing(resolve.duplicated.records)) { resolve.duplicated.records <- FALSE }
  if(missing(pop.missing.coordinates)) { pop.missing.coordinates <- FALSE }
  
  species.id.all.records <- get.data.from.table("Record","SpeciesID")$SpeciesID
  species.id.all.taxa <- get.data.from.table("Taxa","SpeciesWormsID")$SpeciesWormsID
  
  repository.id.all.records <- get.data.from.table("Record","SourceRefID","SourceRefType","Repository")$Repository
  repository.id.all <- get.data.from.table("Repository","fid")$fid
  
  literature.id.all.records <- get.data.from.table("Record","SourceRefID","SourceRefType","Literature")$Literature
  literature.id.all <- get.data.from.table("Literature","fid")$fid
  
  if( FALSE %in% (unique(species.id.all.records) %in% unique(species.id.all.taxa)) ) { errormessage("In Records table there are species (ID) not included in Taxa Table") }
  if( FALSE %in% (unique(repository.id.all.records) %in% unique(repository.id.all)) ) { errormessage("In Records table there are sources (ID) not included in Repository Table") }
  if( FALSE %in% (unique(literature.id.all.records) %in% unique(literature.id.all)) ) { errormessage("In Records table there are sources (ID) not included in Literature Table") }
  
  if( resolve.duplicated.references ) { duplicated.references() }
  if( resolve.unused.references ) { unused.references() }
  if( resolve.references ) { references.crossref(force=TRUE) }
  
  if( resolve.coordinates ) { resolve.coordinates() }
  if( resolve.species ) { duplicated.taxa() ; species.by.worms() }
  if( pop.taxa.information ) { pop.taxa.information(taxa = as.vector(unlist(get.data.from.table("Taxa","SpeciesName")))) }
  
  if( pop.missing.coordinates ) { populate.missing.coordinates() }
  if( resolve.duplicated.records ) { duplicated.records() }
  
  if( backup.db ) {  backup.sql(100,export.raw=TRUE)  }
  
}

## ------------------------------------------------------------------------------------------------------------------

resolve.coordinates <- function() { 
  
  coordinates <- get.data.from.table("Record","fid, Lon, Lat")  
  fids <- coordinates[which( is.na(as.numeric(coordinates$Lon) ) & is.na(coordinates$Lon) | is.na(as.numeric(coordinates$Lat) ) & is.na(coordinates$Lat)),1]
  
  for(fid in fids) {
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbExecute(sql,paste0("UPDATE MDB_GeoData SET Lon = NULL WHERE fid = ", fid)) 
    dbExecute(sql,paste0("UPDATE MDB_GeoData SET Lat = NULL WHERE fid = ", fid)) 
    dbDisconnect(sql)
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

backup.sql <- function(max.versions,export.raw) {
  
  files.sql <- list.files(paste0(main.directory,backups.path),pattern = "sql",full.names = TRUE)
  version <- format(Sys.time(), "%d %b %Y %H.%M")
  
  if( length(files.sql) >= max.versions) {
    
    file.remove(files.sql[1])
    
  }
  
  file.copy( paste0(main.directory,sql.path) , paste0(main.directory,backups.path,"Main.biodiversity.database.v",version,".sql") )
  
  if( export.raw ) {  
    
    for( t in source.tables) {
      
      write.csv2( get.data.from.table(t) , paste0(main.directory,backups.path,"_ Source Tables/",t,".csv"), quote = FALSE)
      
    }
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

list.as.shiny <- function(table) {
  
  if(missing(table)) { errormessage("Table name requiered") }
  
  if( class(table) != "character" & class(table) != "data.frame"  ) { errormessage("Object should be character or data.frame") }
  
  if( class(table) == "data.frame"  ) { temp.data <- table }
  
  if( class(table) == "character" ) {
    
    if(table == "Taxa") { table <- "RDB_Taxa" }
    if(table == "Literature") { table <- "RDB_Ref_Literature" }
    if(table == "User") { table <- "RDB_Users" }
    if(table == "Repository") { table <- "RDB_Ref_Repository" }
    if(table == "Record") { table <- "MDB_GeoData" }
    if(table == "Abundance") { table <- "RDB_Abundance" }
    if(table == "Collection") { table <- "RDB_Sample_Collection" }
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    list.of.DB.tables <- dbListTables(sql)
    dbDisconnect(sql)
    
    if( ! table %in% list.of.DB.tables ) { errormessage("Table name not valid") }
    
    if( table %in% list.of.DB.tables ) {
      
      sql.parse.query <- paste0("SELECT * FROM ",table)
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      temp.data <- dbGetQuery(sql, sql.parse.query )
      dbDisconnect(sql)
    } 
    
    
  }
  
  # https://rstudio.github.io/DT/
  
  shinyApp(
    ui = fluidPage(
      title = "Examples of DataTables",
      DT::dataTableOutput('tbl')),
    server = function(input, output) {
      output$tbl = DT::renderDataTable(
        temp.data, options = list( pageLength = 10, lengthMenu = c(10, 50, 100, 1000)) ) # , filter = 'top'
    } )
}

## ------------------------------------------------------------------------------------------------------------------

add.manual.source <- function(type) {
  
  if( missing(type) | ! type %in% c(source.tables) ) { errormessage("Type of source is missing") }
  
  if(type == "Taxa") { table <- "RDB_Taxa" }
  if(type == "Literature") { table <- "RDB_Ref_Literature" }
  if(type == "User") { table <- "RDB_Users" }
  if(type == "Repository") { table <- "RDB_Ref_Repository" }
  if(type == "Record") { table <- "MDB_GeoData" }
  if(type == "Abundance") { table <- "RDB_Abundance" }
  if(type == "Collection") { table <- "RDB_Sample_Collection" }
  
  data.to.edit.fields <- get.table.structure(type)
  class.fields <- get.data.from.table(type,"*","fid",get.last.fid(type))
  class.fields <- sapply(class.fields,function(x) { class(x) })
  
  repeat {
    
    last.fid <- as.numeric(get.last.fid(type))
    
    cat('\014')
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat(paste0("Adding new entry to ",type))
    cat('\n')
    cat('\n')
    
    fields.to.data.frame <- data.frame(t(rep(NA,length(data.to.edit.fields))),stringsAsFactors=FALSE)
    colnames(fields.to.data.frame) <- data.to.edit.fields
    fields.to.data.frame[,colnames(fields.to.data.frame) == "fid"] <- last.fid + 1
    
    for( k in 2:length(data.to.edit.fields)) {
      
      field.to.edit <- data.to.edit.fields[k]
      new.value <- ""
      
      if( ! field.to.edit %in% non.edited.fields ) { 
        
        new.value <- readline(paste0("Value for ",data.to.edit.fields[k]," of type '", class.fields[k][[1]],"': "))
      }
      
      if( class.fields[k][[1]] == "numeric") { new.value <- as.numeric(new.value) }
      if( new.value == "" & class.fields[k][[1]] == "character") { new.value <- "" }
      
      fields.to.data.frame[,k] <- new.value
      
    }
    
    cat('\014')
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    print(fields.to.data.frame)
    cat('\n')
    
    comite.changes <- ""
    while( comite.changes !="Y" & comite.changes !="n" ) {   comite.changes <- readline(paste0("Add record? (Y/n) ")) }
    
    if(comite.changes == "Y") {
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbWriteTable(sql, table, fields.to.data.frame , append=TRUE, overwrite=FALSE )
      dbDisconnect(sql)
      updated.version.control()
    }
    
    continue <- ""
    while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Add extra record? (Y/n) ")) }
    if( continue == "n") { break }
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

edit.manual.source <- function(type,fid) {
  
  if( missing(type) | ! type %in% source.tables ) { errormessage("Type of source is missing") }
  if( missing(fid) ) { errormessage("fid is missing") }
  
  if( class(fid) != "numeric") { errormessage("fid must be numeric") }
  
  if(type == "Taxa") { table <- "RDB_Taxa" }
  if(type == "Literature") { table <- "RDB_Ref_Literature" }
  if(type == "User") { table <- "RDB_Users" }
  if(type == "Repository") { table <- "RDB_Ref_Repository" }
  if(type == "Record") { table <- "MDB_GeoData" }
  if(type == "Abundance") { table <- "RDB_Abundance" }
  if(type == "Collection") { table <- "RDB_Sample_Collection" }
  
  data.to.edit <- get.data.from.table(type,"*","fid",fid)
  data.to.edit.fields <- colnames(data.to.edit)
  class.fields <- sapply(data.to.edit,function(x) { class(x) })
  
  repeat {
    
    cat('\014')
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat(paste0("Editing fid #",fid," in ",type))
    cat('\n')
    cat('\n')
    print(data.to.edit)
    cat('\n')
    cat('\n')
    
    field.to.edit <- ""
    while( ! field.to.edit %in% data.to.edit.fields ) {  field.to.edit <- readline(paste0("Name the field to edit: ")) }
    
    field.to.edit.location <- which(data.to.edit.fields == field.to.edit)
    cat('\n')
    cat('Current value: ',unlist(data.to.edit[field.to.edit.location]))
    new.value <- ""
    while( new.value == "" ) {  new.value <- readline(paste0("New value of type '", class.fields[field.to.edit.location][[1]] ,"': ")) }
    
    cat('\n')
    
    if( new.value != "" & class.fields[field.to.edit.location][[1]] == "numeric") { new.value <- as.numeric(new.value) }
    
    comite.changes <- ""
    while( comite.changes !="Y" & comite.changes !="n" ) {   comite.changes <- readline(paste0("Change value to '", new.value , "' (Y/n) ")) }
    
    if(comite.changes == "Y") {
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbExecute(sql,paste0("UPDATE ",table," SET ",field.to.edit,"='",new.value,"' WHERE fid = ", fid)) 
      dbDisconnect(sql)   
      updated.version.control()
      
    }
    
    cat('\n')
    continue <- ""
    while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Edit extra field? (Y/n) ")) }
    if( continue == "n") { break }
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

edit.bulk.source <- function(type,fid,field,value,pass) {
  
  if( missing(pass) ) { pass <- FALSE }
  if( missing(type) | ! type %in% source.tables ) { errormessage("Type of source is missing") }
  if( missing(fid) ) { errormessage("Fid is missing") }
  if( missing(field) ) { errormessage("Field is missing") }
  if( missing(value) ) { errormessage("Value is missing") }
  
  if( class(fid) != "numeric") { errormessage("fid must be numeric") }
  
  if(type == "Taxa") { table <- "RDB_Taxa" }
  if(type == "Literature") { table <- "RDB_Ref_Literature" }
  if(type == "User") { table <- "RDB_Users" }
  if(type == "Repository") { table <- "RDB_Ref_Repository" }
  if(type == "Record") { table <- "MDB_GeoData" }
  if(type == "Abundance") { table <- "RDB_Abundance" }
  if(type == "Collection") { table <- "RDB_Sample_Collection" }
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  if( ! field %in% dbListFields(sql,table) ) { errormessage("Invalid field name") }
  dbDisconnect(sql)   
  
  data.to.edit <- get.data.from.table(type,"*","fid",fid)
  data.to.edit.fields <- colnames(data.to.edit)
  class.fields <- sapply(data.to.edit,function(x) { class(x) } )
  
  cat('\014')
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------')
  cat('\n')
  cat(paste0("Bulk editing ", length(fid), " records in table: ",type))
  cat('\n')
  cat('\n')
  
  if( ! pass ) { 
    
    ask <- ""
    while( ask != "Y" & ask != "n" ) {  ask <- readline(paste0("Are you sure? (Y/n) ")) }
    
  } else { ask <- "Y" }
  
  if( ask == "Y") {
    
    fid.to.change <- paste0("'", fid, "'", collapse=", ")
    fid.to.change <- paste0("(", fid.to.change, ")")
    
    saved.success <- FALSE
    
    while( ! saved.success ) {
      
      repeat {
        
        tryCatch( Unloked <- as.character(unlist(read.table(file = paste0(main.directory,"/Data/Main.biodiversity.database.Lk")))) == "Unloked" , finally = break )
        
      }
      
      if( Unloked ) {
        
        write("loked",file=paste0(main.directory,"/Data/Main.biodiversity.database.Lk"),append=FALSE)
        
        sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
        dbExecute(sql,paste0("UPDATE ",table," SET ",field,"='",value,"' WHERE fid IN ", fid.to.change)) 
        dbDisconnect(sql)   
        updated.version.control()
        
        write("Unloked",file=paste0(main.directory,"/Data/Main.biodiversity.database.Lk"),append=FALSE)
        saved.success <- TRUE
        
      }
    }
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

import.db.records <- function(taxa,rank,add.obis,add.gbif,ask) {
  
  if( rank != "species" ) {
    
    list.of.species <- wormsbynames(taxon_names=taxa)$AphiaID
    final.list <- data.frame()
    
    if( ! is.numeric(list.of.species) ) { stop("Error, Higher Taxa not found in Worms") }
    
    higher.taxa <- TRUE
    
    while(higher.taxa){
      
      list.of.childen <- data.frame()
      
      for(sp in 1:length(list.of.species)) {
        
        for( offset in seq(1,1000000,50) ) {
          
          tryCatch( list.of.childen.t <- wm_children(id=list.of.species[sp],offset=offset) , error=function(e) { error <- TRUE })
          
          if( ! exists("list.of.childen.t") ) {  break  }
          
          if( exists("list.of.childen.t") ) { 
            
            list.of.childen.t <- list.of.childen.t[list.of.childen.t$phylum == phylum,]
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
      
      list.of.species <- list.of.species[-to.include]
      
    }
    
  }
  
  # ------------------ 
  
  taxa <- final.list$scientificname
  
  taxa <- unique(taxa)
  taxa <- taxa[! grepl("sp.", taxa)]
  
  new.species <- data.frame()
  new.records.f <- data.frame()
  species.taxonomy.to.import <- data.frame()
  
  for( t in 1:length(taxa) ) {
    
    taxa.to.search <- taxa[t]
    
    queryied.data <- get.data.from.table("Taxa","*","SpeciesName",taxa.to.search)
    
    if( nrow(queryied.data) == 0 ) {
      
      new.species.t <- get.species.taxonomy(taxa.to.search,AllYes=TRUE)
      new.species <- rbind(new.species,new.species.t$worms.taxonomy)
      species.taxonomy.to.import <- rbind(species.taxonomy.to.import,new.species.t$worms.taxonomy)
      
    }
    
    if(nrow(queryied.data) != 0) {
      
      species.taxonomy.to.import <- rbind(species.taxonomy.to.import,queryied.data)
      
    }
    
  }
  
  # ------------------ 
  
  for( sp in 1:nrow(species.taxonomy.to.import)) {
    
    taxaname <- as.character(species.taxonomy.to.import[sp,]$SpeciesName)
    taxaid <- as.character(species.taxonomy.to.import[sp,]$SpeciesWormsID)
    
    cat("\n",taxaname)
    
    subset <- get.external.db(taxaname,taxaid,add.obis=add.obis,add.gbif=add.gbif)
    
    if( nrow(subset) > 0 ) {
      
      new.records.f <- rbind(new.records.f,
                             data.frame(  fid = rep(NA,nrow(subset)),
                                          SpeciesID = as.numeric(as.character(subset$SpeciesID)) ,
                                          Lon = as.numeric(as.character(subset$Lon)) ,
                                          Lat = as.numeric(as.character(subset$Lat)) ,
                                          CoordinateType = rep("",nrow(subset)) ,
                                          CoordinateUncertaintyInMeters = as.numeric(as.character(subset$CoordinateUncertaintyInMeters)) ,
                                          Country = as.character(subset$Country) ,
                                          Site = as.character(subset$SiteName) ,
                                          Depth = as.numeric(as.character(subset$Depth)) ,
                                          DepthRangeMin = as.numeric(as.character(subset$DepthRangeMin)) ,
                                          DepthRangeMax = as.numeric(as.character(subset$DepthRangeMax)) ,
                                          DateYear = as.numeric(as.character(subset$DateYear)) ,
                                          DateMonth = rep(NA,nrow(subset)) ,
                                          DateDay = rep(NA,nrow(subset)) ,
                                          AbundanceID = rep(NA,nrow(subset)) ,
                                          SampleCollectionID = rep(NA,nrow(subset)) ,
                                          SourceRefType = as.character(subset$SourceRefType) ,
                                          SourceRefID = as.character(subset$SourceRefID) ,
                                          RecordUserID = rep(NA,nrow(subset)) ,
                                          RecordDateYear = as.numeric(as.character(rep(format(Sys.time(), "%Y"),nrow(subset)))) ,
                                          RecordDateMonth = as.numeric(as.character(rep(format(Sys.time(), "%m"),nrow(subset)))) ,
                                          RecordDateDay = as.numeric(as.character(rep(format(Sys.time(), "%d"),nrow(subset)))) ,
                                          RecordPublic = rep(NA,nrow(subset)) ,
                                          RecordNotes = as.character(subset$Notes) , 
                                          OriginalSource = as.character(subset$OriginalSource) ,
                                          Flag = rep("",nrow(subset)),
                                          OriginalSourceID = as.character(subset$OriginalSourceID) ,
                                          FlagOnLand = rep("",nrow(subset)),
                                          FlagLightSuitable = rep("",nrow(subset)),
                                          FlagO2Suitable = rep("",nrow(subset)),
                                          FlagEcoRegion = rep("",nrow(subset)),
                                          
                                          stringsAsFactors=FALSE) )
    }
    
  }
  
  # ------------------
  
  new.records.f <- new.records.f[ new.records.f$Lon >= -180 & new.records.f$Lon <= 180 & new.records.f$Lat >= -90 & new.records.f$Lat <= 90 ,]
  
  new.species <- new.species[ new.species$SpeciesWormsID %in% new.records.f$SpeciesID, ]
  
  if( ! exists("new.species") | nrow(new.species) == 0) { new.species <- NULL }
  
  new.source.references <- NULL
  
  # ------------------ 
  
  cat('\014')
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------')
  cat('\n')
  cat('Number of records:', nrow(new.records.f))
  cat('\n')
  cat('Number of species:', length(unique(new.records.f$SpeciesID)))
  cat('\n')
  cat('\n')
  
  if(nrow(new.records.f) > 0) {
    
    continue <- ""
    
    if(!ask) { continue <- "Y"}
    
    while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Plot new records? (Y/n) ")) }
    if( continue == "Y") { 
      
      return.plot <- plot.coords.google(new.records.f[,3:4],new.records.f[,1],4,"black")
      print(return.plot) 
      
    }
    
    continue <- ""
    
    if(!ask) { continue <- "Y"}
    
    while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Inject new records to Database? (Y/n) ")) }
    
    if( continue == "Y") { 
      
      if( sum( ! is.na(new.records.f$CoordinateUncertaintyInMeters) ) > 0 ) {
        
        continue.2 <- ""
        while( continue.2 !="Y" & continue.2 !="n" ) { continue.2 <- readline(paste0("Some records have information about coordinate uncertainty. Do you wish to trim records? (Y/n) ")) }
        
        if( continue.2 == "Y") { 
          
          threhold.meters <- readline(paste0("Give a threshold in meters: "))
          new.records.f <- new.records.f[ which(is.na(new.records.f$CoordinateUncertaintyInMeters) | new.records.f$CoordinateUncertaintyInMeters <= as.numeric(threhold.meters)) , ]
          
        }
      }
      
      inject.to.sql(new.records.f,new.species,new.source.references) 
      
    }
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

get.species.taxonomy <- function(taxa,AllYes) {
  
  if( missing(taxa) ) { errormessage("Taxa is missing") }
  
  if( length(taxa) == 0 ) { errormessage("Taxa is missing") }
  
  if( length(taxa) > 0 ) {
    
    species.taxonomy <- data.frame()
    
    last.non.vallid.worms.id <- as.numeric(unlist(get.data.from.table("Taxa","SpeciesWormsID")))
    last.non.vallid.worms.id <- last.non.vallid.worms.id[last.non.vallid.worms.id < 1000]
    
    if( length(last.non.vallid.worms.id) == 0) { last.non.vallid.worms.id <- 0 }
    
    for( taxa.i in 1:length(taxa) ) { 
      
      taxa.to.search <- taxa[taxa.i]
      
      cat('\014')
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat("\n")
      cat("Getting taxonomy information for",taxa.to.search, " (",taxa.i,"out of",length(taxa),")")
      cat("\n")
      cat("\n")
      
      repeat {
        
        tryCatch( worms.query <- wormsbymatchnames( taxa.to.search , verbose = FALSE ), error=function(e) { } )
        
        if( ! exists("worms.query") ) { 
          
          use.worms <- FALSE
          
          cat('\n')
          cat('\n')
          cat(paste0("Species name not included in WORMS: ",taxa.to.search))
          
          if(AllYes) { x <- "Y" } else { x <- "" }
          
          while( x != "Y" & x != "n" ) { 
            x <- readline(paste0("Use current name? (Y/n): "))
          }
          
          if ( x == "n" ) { 
            
            taxa.to.search <- readline(paste0("Type new name: "))
            cat('\n')
            
          }
          
          if ( x == "Y" ) { 
            
            last.non.vallid.worms.id <- max(last.non.vallid.worms.id) + 1
            
            species.taxonomy.temp <-   data.frame(   fid = NA ,
                                                     SpeciesWormsID = last.non.vallid.worms.id ,
                                                     SpeciesName = taxa.to.search ,
                                                     Authority = "",
                                                     Status = "",
                                                     TaxKingdom = "",
                                                     TaxPhylum = "",
                                                     TaxClass = "",
                                                     TaxOrder = "",
                                                     TaxFamily = "",
                                                     TaxGenus = "",
                                                     ReviewedByWorms = 0,
                                                     RevisionByWormsDate = "",
                                                     AcceptedWormsID = NA,
                                                     AcceptedSpeciesName= "" ,
                                                     Description = "" , 
                                                     DescriptionSource = "" , 
                                                     DescriptionSourceLink = "" , 
                                                     Distribution = "" , 
                                                     DistributionSource = "" , 
                                                     CommonName = "" , stringsAsFactors = FALSE)
            break;
            
          }
        }
        
        if( exists("worms.query") ) {  
          
          use.worms <- FALSE
          
          if( taxa.to.search != worms.query$scientificname ) {  
            
            cat('\n')
            cat('\n')
            cat(paste0("Error in species name: ",taxa.to.search , " differs from, ",worms.query$scientificname, " in WORMS"))
            cat('\n')
            
            if(AllYes) { x <- "Y" } else { x <- "" }
            
            while( x != "Y" & x != "n" ) {  x <- readline(paste0("Use WORMS suggested name? (Y/n): ")) }
            
            if ( x == "Y" ) {  use.worms <- TRUE }
            
            if ( x == "n" ) { 
              
              use.worms <- FALSE
              
              last.non.vallid.worms.id <- last.non.vallid.worms.id + 1
              species.taxonomy.temp <-  data.frame(   fid = NA ,
                                                      SpeciesWormsID = last.non.vallid.worms.id ,
                                                      SpeciesName = taxa.to.search ,
                                                      Authority = "",
                                                      Status = "",
                                                      TaxKingdom = "",
                                                      TaxPhylum = "",
                                                      TaxClass = "",
                                                      TaxOrder = "",
                                                      TaxFamily = "",
                                                      TaxGenus = "",
                                                      ReviewedByWorms = 0,
                                                      RevisionByWormsDate = "",
                                                      AcceptedWormsID = NA,
                                                      AcceptedSpeciesName= "" ,
                                                      Description = "" , 
                                                      DescriptionSource = "" , 
                                                      DescriptionSourceLink = "" , 
                                                      Distribution = "" , 
                                                      DistributionSource = "" , 
                                                      CommonName = "" , stringsAsFactors = FALSE)
              break;
            }
          }
          
          if( taxa.to.search == worms.query$scientificname | use.worms ) {  
            
            species.taxonomy.temp <-  data.frame(   fid = NA ,
                                                    SpeciesWormsID = worms.query$AphiaID ,
                                                    SpeciesName = worms.query$scientificname ,
                                                    Authority = worms.query$authority,
                                                    Status = worms.query$status,
                                                    TaxKingdom = worms.query$kingdom,
                                                    TaxPhylum = worms.query$phylum,
                                                    TaxClass = worms.query$class,
                                                    TaxOrder = worms.query$order,
                                                    TaxFamily = worms.query$family,
                                                    TaxGenus = worms.query$genus,
                                                    ReviewedByWorms = 1,
                                                    RevisionByWormsDate = as.character(Sys.Date()),
                                                    AcceptedWormsID = worms.query$valid_AphiaID,
                                                    AcceptedSpeciesName= worms.query$valid_name,
                                                    Description = "" , 
                                                    DescriptionSource = "" , 
                                                    DescriptionSourceLink = "" , 
                                                    Distribution = "" , 
                                                    DistributionSource = "" , 
                                                    CommonName = "" , stringsAsFactors = FALSE)
            
            rm(worms.query)
            break;
          }
          
        }
        
        
      }
      
      species.taxonomy <- rbind(species.taxonomy , species.taxonomy.temp)
      
      if(exists("worms.query")) { rm(worms.query) }
      
    }
    
  }
  
  if( ! exists("species.taxonomy") ) { new.species <- NULL }
  
  return(list(original.taxa=taxa,worms.taxonomy=species.taxonomy))
  
}


## ------------------------------------------------------------------------------------------------------------------

remove.unused.species <- function() {
  
  internet.available()
  
  taxa.db <- get.data.from.table("Taxa","fid, SpeciesWormsID, SpeciesName")
  record.db <- get.data.from.table("Record","fid, SpeciesID")
  
  not.included <- taxa.db[ which( ! taxa.db$SpeciesWormsID %in% record.db$SpeciesID ) , ]
  
  delete.records("Taxa", not.included$fid ,FALSE)
  
}

# Check if Records have proper Taxa, if not add to Taxa table

add.missing.species.by.worms <- function() {
  
  internet.available()
  
  taxa.db <- get.data.from.table("Taxa","fid, SpeciesWormsID, SpeciesName")
  record.db <- get.data.from.table("Record","fid, SpeciesID")
  
  not.included <- unique(record.db$SpeciesID[which( ! record.db$SpeciesID %in% taxa.db$SpeciesWormsID )])
  not.included <- not.included[not.included > 100]
  not.included <- not.included[!is.na(not.included)]
  
  if( length(not.included) > 0 ) {
    
    for(g in 1:length(not.included) ) {
      
      tryCatch( worms.query <- wormsbyid( not.included[g] , verbose = FALSE )
                , error=function(e) { } )
      
      if( exists("worms.query") ) { 
        
        cat('\014')
        cat('\n')
        cat('\n')
        cat('----------------------------------------------------------------------')
        cat('\n')
        cat('Processing Species', g, "out of", length(not.included))
        cat('\n')
        cat('\n')
        cat('\n')
        
        append.result <- data.frame(
          
          fid = get.last.fid("Taxa") + 1,
          SpeciesWormsID = worms.query$AphiaID,
          SpeciesName = worms.query$scientificname,
          Authority = worms.query$authority,
          Status = worms.query$status,
          TaxKingdom = worms.query$kingdom,
          TaxPhylum = worms.query$phylum,
          TaxClass = worms.query$class,
          TaxOrder = worms.query$order,
          TaxFamily = worms.query$family,
          TaxGenus = worms.query$genus,
          ReviewedByWorms = 1,
          RevisionByWormsDate = Sys.Date(),
          AcceptedWormsID = worms.query$valid_AphiaID,
          AcceptedSpeciesName = worms.query$valid_name,
          Description = "",
          DescriptionSource = "",
          DescriptionSourceLink = "",
          Distribution = "",
          DistributionSource = "",
          CommonName = ""
        )
        
        
        sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
        dbWriteTable(sql, "RDB_Taxa", append.result , append=TRUE, overwrite=FALSE )
        dbDisconnect(sql)
        rm(worms.query)
        
      }
      
    }
  }  
  
}

## ------------------------------------------------------------------------------------------------------------------

# Populate Taxa Table with worms information

species.by.worms <- function(force.review) {
  
  if( missing(force.review) ) { force.review <- FALSE }
  
  internet.available()
  
  review.by.worms.db <- get.data.from.table("Taxa","fid, SpeciesWormsID, SpeciesName, ReviewedByWorms, AcceptedWormsID")
  review.by.worms.db <- review.by.worms.db[review.by.worms.db$ReviewedByWorms == 0 ,]
  
  if( force.review ) { review.by.worms.db <- get.data.from.table("Taxa","fid, SpeciesWormsID, SpeciesName, ReviewedByWorms, AcceptedWormsID") }
  
  if( nrow(review.by.worms.db) > 0 ) {
    
    for(g in 1:nrow(review.by.worms.db) ) {
      
      cat('\014')
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat('Processing', g, "out of", nrow(review.by.worms.db))
      cat('\n')
      cat('\n')
      cat('\n')
      
      species.name <- review.by.worms.db$SpeciesName[g]
      fid.to.change <- review.by.worms.db$fid[g]
      OriginalSpeciesWormsID <- review.by.worms.db$SpeciesWormsID[g]
      
      if( exists("worms.query") ) { rm(worms.query) }
      
      tryCatch( worms.query <- wormsbymatchnames( species.name , verbose = FALSE )
                , error=function(e) { } )
      
      if( exists("worms.query") ) { 
        
        SpeciesWormsID = worms.query$AphiaID
        SpeciesName = worms.query$scientificname
        Authority = worms.query$authority
        Status = worms.query$status
        TaxKingdom = worms.query$kingdom
        TaxPhylum = worms.query$phylum
        TaxClass = worms.query$class
        TaxOrder = worms.query$order
        TaxFamily = worms.query$family
        TaxGenus = worms.query$genus
        AcceptedWormsID = worms.query$valid_AphiaID
        AcceptedSpeciesName = worms.query$valid_name
        
        Authority <- gsub("'", "", Authority)
        
        sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET SpeciesWormsID='",SpeciesWormsID,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET SpeciesName='",SpeciesName,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET Authority='",Authority,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET Status='",Status,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET TaxKingdom='",TaxKingdom,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET TaxPhylum='",TaxPhylum,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET TaxClass='",TaxClass,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET TaxOrder = '",TaxOrder,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET TaxFamily='",TaxFamily,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET TaxGenus='",TaxGenus,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET AcceptedWormsID='",AcceptedWormsID,"' WHERE fid = ", fid.to.change)) 
        dbExecute(sql,paste0("UPDATE RDB_Taxa SET AcceptedSpeciesName='",AcceptedSpeciesName,"' WHERE fid = ", fid.to.change)) 
        dbDisconnect(sql)
        
        if( SpeciesWormsID != OriginalSpeciesWormsID ) {
          
          sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
          dbExecute(sql,paste0("UPDATE MDB_GeoData SET SpeciesID='",SpeciesWormsID,"' WHERE SpeciesID = ", OriginalSpeciesWormsID )) 
          dbDisconnect(sql)
          
        }
        
        rm(worms.query)
        
        
      }
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbExecute(sql,paste0("UPDATE RDB_Taxa SET ReviewedByWorms = 1 WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Taxa SET RevisionByWormsDate='",as.character(Sys.Date()),"' WHERE fid = ", fid.to.change)) 
      dbDisconnect(sql)
      updated.version.control()
      
    }
    
  }
  
  # -------------------------------
  
  review.by.worms.db <- get.data.from.table("Taxa","fid, SpeciesWormsID, SpeciesName, ReviewedByWorms, AcceptedWormsID, AcceptedSpeciesName")
  
  accepted.taxa <- review.by.worms.db$AcceptedWormsID
  accepted.taxa.to.include <- review.by.worms.db$AcceptedWormsID[which( ! accepted.taxa %in% review.by.worms.db$SpeciesWormsID )]
  accepted.taxa.to.include <- accepted.taxa.to.include[accepted.taxa.to.include != "" ]
  accepted.taxa.to.include <- accepted.taxa.to.include[ ! is.na(accepted.taxa.to.include) ]
  accepted.taxa.to.include <- accepted.taxa.to.include[accepted.taxa.to.include != "NA" ]
  accepted.taxa.to.include <- accepted.taxa.to.include[accepted.taxa.to.include != "0" ]
  accepted.taxa.to.include <- unique(accepted.taxa.to.include)
  
  if( length(accepted.taxa.to.include) > 0 ) {
    
    for(g in 1:length(accepted.taxa.to.include) ) {
      
      cat('\014')
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat('Processing', g, "out of", length(accepted.taxa.to.include))
      cat('\n')
      cat('\n')
      cat('\n')
      
      species.name <- accepted.taxa.to.include[g]
      
      if( exists("worms.query") ) { rm(worms.query) }
      
      tryCatch( worms.query <- wormsbyid( species.name , verbose = FALSE )
                , error=function(e) { } )
      
      if( exists("worms.query") ) { 
        
        fields.to.data.frame <- data.frame(   fid = get.last.fid("Taxa") + 1,
                                              SpeciesWormsID = worms.query$AphiaID,
                                              SpeciesName = worms.query$scientificname,
                                              Authority = worms.query$authority,
                                              Status = worms.query$status,
                                              TaxKingdom = worms.query$kingdom,
                                              TaxPhylum = worms.query$phylum,
                                              TaxClass = worms.query$class,
                                              TaxOrder = worms.query$order,
                                              TaxFamily = worms.query$family,
                                              TaxGenus = worms.query$genus,
                                              ReviewedByWorms = 1,
                                              RevisionByWormsDate = as.character(Sys.Date()),
                                              AcceptedWormsID = worms.query$valid_AphiaID,
                                              AcceptedSpeciesName = worms.query$valid_name,
                                              Description = "",
                                              DescriptionSource = "",
                                              DescriptionSourceLink = "",
                                              Distribution = "",
                                              DistributionSource = "",
                                              CommonName = "" , stringsAsFactors = FALSE)
        
        sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
        dbWriteTable(sql, "RDB_Taxa", fields.to.data.frame , append=TRUE, overwrite=FALSE )
        dbDisconnect(sql)
        updated.version.control()
        
        rm(worms.query)
        
      }
      
    }
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

duplicated.taxa <- function() {
  
  # Remove duplicates
  
  replicated.worms.id.db <- get.data.from.table("Taxa","fid, SpeciesWormsID")
  replicated.worms.id.db <- replicated.worms.id.db[replicated.worms.id.db$SpeciesWormsID != "",]
  
  replicated.worms.id <- which(duplicated(replicated.worms.id.db[,2]))
  unique.replicated.worms.id <- unique(replicated.worms.id.db[replicated.worms.id,2])
  
  if( length(unique.replicated.worms.id) > 0) {
    for(z in 1:length(unique.replicated.worms.id) ) {
      
      fid.to.delete <- replicated.worms.id.db[which( replicated.worms.id.db[,2] == unique.replicated.worms.id[z]),1]
      fid.to.keep <- fid.to.delete[1]
      fid.to.delete <- fid.to.delete[-1]
      
      fid.to.delete <- paste0("'", fid.to.delete, "'", collapse=", ")
      fid.to.delete <- paste0("(", fid.to.delete, ")")
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbExecute(sql,paste0("DELETE FROM RDB_Taxa WHERE fid IN ", fid.to.delete)) 
      dbDisconnect(sql)
      updated.version.control()
      
    }
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

references.crossref <- function(force) {
  
  if(missing(force)) { force <- FALSE }
  tolerance.consequtive.words <- 6
  internet.available()
  
  gaps.db <- get.data.from.table("Literature","*")
  
  if( ! force ) {
    
    gaps <- which(gaps.db$ReviewedByCrossref == 0 | is.na(gaps.db$ReviewedByCrossref) | gaps.db$ReviewedByCrossref == "" )
    
  }
  if( force ) {
    
    gaps <- which(gaps.db$Type == "" )
    
  }
  
  
  if( length(gaps) > 0) {
    
    gaps.id <- gaps.db$fid[gaps]
    main.source.list <- gaps.db$FullCitation[gaps]
    
    original.references <- get.data.from.table("Literature","*","fid","1")[-1,]
    crossref.references <- get.data.from.table("Literature","*","fid","1")[-1,]
    
    for( e in 1:length(main.source.list)) {
      
      cat('\014')
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat('Processing reference', e, "out of", length(main.source.list))
      cat('\n')
      cat('\n')
      cat('\n')
      
      query <- main.source.list[e]
      
      sink("/dev/null") 
      tryCatch( titles <- cr_types(works=TRUE, types = cross.ref.types, query = query) 
                , error=function(e) { } )
      sink()
      
      there.is.a.match <- FALSE
      
      if( exists("titles") ) {  
        
        list.up.to <- 100
        list.titles <- data.frame(Sources = paste0(0," : None"))
        
        for(j in 1:list.up.to ) {
          
          if( !is.null( titles$data[1:list.up.to,][j,]$title ) & !is.na(titles$data[1:list.up.to,][j,]$title) ) {
            list.titles <- rbind(list.titles,data.frame( Sources = paste0(j," : ", as.character(substr( titles$data[1:list.up.to,][j,]$title , 1, 100))) ))
          }
        }
        
        for( j in 1:list.up.to ) {
          
          title.to.test <- titles$data[1:list.up.to,][j,]$title
          
          clean.string.1 <- gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", tolower(title.to.test))
          clean.string.2 <- gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", tolower(query))
          
          clean.string.1 <- gsub("[[:punct:]]", "", clean.string.1)
          clean.string.2 <- gsub("[[:punct:]]", "", clean.string.2)
          
          clean.string.1 <- encoding.correction(clean.string.1)
          clean.string.2 <- encoding.correction(clean.string.2)
          
          clean.string.1 <- gsub("`|\\'", "", iconv(clean.string.1, to="ASCII//TRANSLIT"))
          clean.string.2 <- gsub("`|\\'", "", iconv(clean.string.2, to="ASCII//TRANSLIT"))
          
          consequtive.words <- consequtive.words.string.match(clean.string.1 , clean.string.2)
          
          if(consequtive.words >= tolerance.consequtive.words) {  
            
            there.is.a.match <- TRUE
            matcher <- j
            break
            
          } else {  there.is.a.match <- FALSE  }
          
          
        }
        
      }
      
      # --- 
      
      if( ! there.is.a.match ) {    
        
        DOI = "" 
        Type = "" 
        Authors = "" 
        Title = "" 
        Year = "" 
        Container = "" 
        Volume = "" 
        Issue = "" 
        Page = "" 
        FullCitation = simpleCap( tolower( encoding.correction(query) )) 
        
      }
      
      if ( there.is.a.match ) { 
        
        x <- as.numeric(matcher)
        x.names <- names(titles$data[1:list.up.to,][x,])
        
        if( "DOI" %in% x.names ) { DOI <- titles$data[1:list.up.to,][x,]$DOI } else { DOI <- "" }
        if( "type" %in% x.names ) { Type <- titles$data[1:list.up.to,][x,]$type } else { Type <- "" }
        if( "title" %in% x.names ) { Title <- titles$data[1:list.up.to,][x,]$title } else { Title <- "" }
        if( "created" %in% x.names ) { Year <- year(titles$data[1:list.up.to,][x,]$created) } else { Year <- "" }
        if( "container.title" %in% x.names ) { Container <- titles$data[1:list.up.to,][x,]$container.title } else { Container <- "" }
        if( "volume" %in% x.names ) { Volume <- titles$data[1:list.up.to,][x,]$volume } else { Volume <- "" }
        if( "issue" %in% x.names ) { Issue <- titles$data[1:list.up.to,][x,]$issue } else { Issue <- "" }
        if( "page" %in% x.names ) { Page <- titles$data[1:list.up.to,][x,]$page } else { Page <- "" }
        
        Authors <- cbind(c( titles$data[1:list.up.to,][x,]$author[[1]])$given , c( titles$data[1:list.up.to,][x,]$author[[1]])$family)
        if( is.null(Authors) ) { Authors <- "" } else { Authors <- paste(apply(Authors,1,function(x) {   paste0(  simpleCap(x[2])," ", simpleCap(x[1]) )   } ), collapse = "; ")  }
        
        Title <- simpleCap( tolower( Title ) )
        if( substr( Title , nchar(Title), nchar(Title) ) == " " ) { Title <- substr( Title , 1, nchar(Title)-1 ) }
        
        if( Authors == "" ) { 
          
          FullCitation <- paste0( "(" , Year , ") " , Title , ". " ,Container, ", " ,Volume, ": " ,Issue, ", " , Page )
          
        } else {
          
          FullCitation <- paste0( Authors , " (" , Year , ") " , Title , ". " ,Container, ", " ,Volume, ": " ,Issue, ", " , Page )
          
        }
        
        FullCitation <- gsub("NA", "", FullCitation)
        FullCitation <- gsub(": ,", ":", FullCitation)
        FullCitation <- gsub('\"', '', FullCitation)
        FullCitation <- gsub("'", "", FullCitation)
        rm(titles)
        
      }
      
      fid.to.change <- gaps.id[e]
      
      original.references <- rbind(original.references,data.frame(fid = fid.to.change, 
                                                                  Type = "", 
                                                                  Authors = "", 
                                                                  Title = "", 
                                                                  Year = "", 
                                                                  Container = "",
                                                                  Volume = "", 
                                                                  Issue = "", 
                                                                  Page = "", 
                                                                  FullCitation = simpleCap( tolower( query )) , 
                                                                  DOI = "", 
                                                                  ReviewedByCrossref = 1, 
                                                                  RevisionByCrossrefDate = as.character(Sys.time()) , stringsAsFactors = FALSE  ))
      
      crossref.references <- rbind(crossref.references,data.frame(fid = fid.to.change, 
                                                                  Type = Type, 
                                                                  Authors = Authors, 
                                                                  Title = Title, 
                                                                  Year = Year, 
                                                                  Container = Container,
                                                                  Volume = Volume, 
                                                                  Issue = Issue, 
                                                                  Page = Page, 
                                                                  FullCitation = FullCitation , 
                                                                  DOI = DOI, 
                                                                  ReviewedByCrossref = 1, 
                                                                  RevisionByCrossrefDate = as.character(Sys.time()) , stringsAsFactors = FALSE  ))
      
    }
    
    for( e in 1:length(main.source.list)) {
      
      if( crossref.references[e,]$Type == "" ) {
        
        data.to.use <- original.references[e,]
        
      }
      
      if( crossref.references[e,]$Type != "" ) {
        
        cat('\014')
        cat('\n')
        cat('\n')
        cat('----------------------------------------------------------------------')
        cat('\n')
        cat('Original reference:')
        cat('\n')
        cat(as.character(original.references[e,]$FullCitation))
        cat('\n')
        cat('\n')
        cat('Reference found in crossref:')
        cat('\n')
        cat(as.character(crossref.references[e,]$FullCitation))
        cat('\n')
        cat('\n')
        
        comite.changes <- "Key"
        while( comite.changes !="n" & comite.changes !="" ) {   comite.changes <- readline(paste0("Does the references match? (Return/n) ")) }
        
        if(comite.changes =="" ) {
          data.to.use <- crossref.references[e,]
        }
        if(comite.changes == "n" ) {
          data.to.use <- original.references[e,]
        }
        
        
      }
      
      fid.to.change <- data.to.use$fid
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET Type='",data.to.use$Type,"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET Authors='",gsub("'","",data.to.use$Authors),"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET Title='",gsub("'","",data.to.use$Title),"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET Year='",data.to.use$Year,"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET Container='",data.to.use$Container,"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET Volume='",data.to.use$Volume,"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET Issue='",data.to.use$Issue,"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET Page='",data.to.use$Page,"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET FullCitation='",gsub("'","",data.to.use$FullCitation), "' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET DOI='",data.to.use$DOI,"' WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET ReviewedByCrossref=1 WHERE fid = ", fid.to.change)) 
      dbExecute(sql,paste0("UPDATE RDB_Ref_Literature SET RevisionByCrossrefDate='",as.character(Sys.time()),"' WHERE fid = ", fid.to.change)) 
      dbDisconnect(sql)
      
    }
    
    updated.version.control()
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

duplicated.references <- function(force.review) {
  
  if( missing(force.review)) { force.review <- FALSE }
  
  # ----------------------------------------------------------
  
  # Absolute FullCitation duplicated
  
  replicated.FullCitation.db <- get.data.from.table("Literature","fid, FullCitation")
  replicated.FullCitation.db <- replicated.FullCitation.db[replicated.FullCitation.db$FullCitation != "",]
  
  replicated.FullCitation <- which(duplicated(replicated.FullCitation.db$FullCitation))
  unique.replicated.FullCitation <- unique(replicated.FullCitation.db[replicated.FullCitation,2])
  
  if( length(unique.replicated.FullCitation) > 0) {
    for(z in 1:length(unique.replicated.FullCitation) ) {
      
      fid.to.delete <- replicated.FullCitation.db[which( replicated.FullCitation.db[,2] == unique.replicated.FullCitation[z]),1]
      fid.to.keep <- fid.to.delete[1]
      fid.to.delete <- fid.to.delete[-1]
      
      fid.to.delete <- paste0("'", fid.to.delete, "'", collapse=", ")
      fid.to.delete <- paste0("(", fid.to.delete, ")")
      
      fid.to.keep <- paste0("'", fid.to.keep, "'", collapse=", ")
      fid.to.keep <- paste0("(", fid.to.keep, ")")
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbExecute(sql,paste0("DELETE FROM RDB_Ref_Literature WHERE fid IN ", fid.to.delete)) 
      dbExecute(sql,paste0("UPDATE MDB_GeoData SET SourceRefID=",fid.to.keep," WHERE SourceRefID IN ", fid.to.delete)) 
      dbDisconnect(sql)
      updated.version.control()
      
    }
    
  }
  
  # ----------------------------------------------------------
  
  # Relative FullCitation duplicated
  
  internet.available()
  tolerance.consequtive.words <- 10
  
  gaps.db <- get.data.from.table("Literature","fid, ReviewedByCrossref")
  gaps <- which(gaps.db$ReviewedByCrossref == 0 | is.na(gaps.db$ReviewedByCrossref) | gaps.db$ReviewedByCrossref == "" )
  gaps.id <- gaps.db[gaps,]$fid
  
  if( force.review ) {  gaps.id <- as.numeric(unlist( get.data.from.table("Literature","fid") ))    }
  
  if( length(gaps) > 0 ) {
    
    replicated.FullCitation.db <- get.data.from.table("Literature","fid, FullCitation","fid",gaps.id)
    replicated.FullCitation.db <- replicated.FullCitation.db[replicated.FullCitation.db$FullCitation != "",]
    
    replicated.FullCitation <- numeric(0)
    
    for( i in 1:(nrow(replicated.FullCitation.db)-1)) {
      
      if( i == 0) { next }
      
      for( j in (i+1):nrow(replicated.FullCitation.db) ) {
        
        if( i == j ) { next }
        
        consequtive.words <- consequtive.words.string.match(replicated.FullCitation.db$FullCitation[i] , replicated.FullCitation.db$FullCitation[j])
        
        if(consequtive.words >= tolerance.consequtive.words) {  replicated.FullCitation <- c(j,replicated.FullCitation)  }
        
      }
      
    }
    
    unique.replicated.FullCitation <- unique(replicated.FullCitation.db[replicated.FullCitation,2])
    
    if( length(unique.replicated.FullCitation) > 0) {
      for(z in 1:length(unique.replicated.FullCitation) ) {
        
        fid.to.delete <- replicated.FullCitation.db[which( replicated.FullCitation.db[,2] == unique.replicated.FullCitation[z]),1]
        fid.to.keep <- fid.to.delete[1]
        fid.to.delete <- fid.to.delete[-1]
        
        fid.to.delete <- paste0("'", fid.to.delete, "'", collapse=", ")
        fid.to.delete <- paste0("(", fid.to.delete, ")")
        
        fid.to.keep <- paste0("'", fid.to.keep, "'", collapse=", ")
        fid.to.keep <- paste0("(", fid.to.keep, ")")
        
        sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
        dbExecute(sql,paste0("DELETE FROM RDB_Ref_Literature WHERE fid IN ", fid.to.delete)) 
        dbExecute(sql,paste0("UPDATE MDB_GeoData SET SourceRefID=",fid.to.keep," WHERE SourceRefID IN ", fid.to.delete)) 
        dbDisconnect(sql)
        updated.version.control()
        
      }
      
    }
    
  }
}

## ------------------------------------------------------------------------------------------------------------------

unused.references <- function() {
  
  citations.fid <- as.numeric(unlist(get.data.from.table("Literature","fid")))
  records.citations.fid <- as.numeric(unlist(get.data.from.table("Record","SourceRefID","SourceRefType","Literature")))
  unused.citations <- citations.fid[ !  citations.fid %in% records.citations.fid]
  
  if( length(unused.citations) > 0 ) {
    
    cat('\014')
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------------------')
    cat('\n')
    cat('Found',length(unused.citations), 'unused citations')
    cat('\n')
    
    x <- ""
    while( x != "Y" & x != "n" ) { 
      x <- readline(paste0("Delete citations? (Y/n): "))
    }
    
    if ( x == "Y" ) { 
      
      unused.citations <- paste0("'", unused.citations, "'", collapse=", ")
      unused.citations <- paste0("(", unused.citations, ")")
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbExecute(sql,paste0("DELETE FROM RDB_Ref_Literature WHERE fid IN ", unused.citations)) 
      dbDisconnect(sql)
      updated.version.control()
      
    }
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

inject.species.from.file <- function(file,type,decimal,xls.sheet) {
  
  internet.available()
  
  if(type == "xlsx") {
    
    
    tryCatch( new.records <- gdata::read.xls(file,sheet=xls.sheet)
              , error=function(e) { } )
    
  }
  if(type == "csv") {
    tryCatch( new.records <- read.csv(file,dec=decimal, sep=";",header=TRUE)
              , error=function(e) { } )
    
  }
  
  if( ! exists("new.records") ) { errormessage("Error while reading file") }
  
  new.records <- new.records[!apply(is.na(new.records) | new.records == "", 1, all),]
  new.records <- as.character(new.records)
  
  if( length(new.records) == 0 ) { errormessage("Empty file") }
  
  # Correct species name
  
  species.in.new.records <- character(length(new.records))
  
  for(sp.n in 1:length(new.records)) {
    
    genus <- strsplit(new.records[sp.n], " ")[[1]]
    genus <- paste0( toupper(substring(genus, 1, 1 )) , tolower(substring(genus, 2, nchar(genus) )))[1]
    species <- tolower(strsplit(new.records[sp.n], " ")[[1]][2])
    
    if( !is.na(species) & species != "sp." ) {
      
      species.in.new.records[sp.n] <- as.character(paste0(genus," ",species))
      
    } else {
      
      species.in.new.records[sp.n] <- as.character(paste0(genus))
      
    }
    
  }
  
  # Which species not in DB
  
  species.in.new.records <- unique(species.in.new.records)
  
  species.in.db <- get.data.from.table("Taxa","SpeciesName")
  new.taxa <- which( ! species.in.new.records %in% as.character(species.in.db$SpeciesName))
  new.taxa <- species.in.new.records[new.taxa]
  
  if( length(new.taxa) > 0) {
    
    new.species <- get.species.taxonomy(new.taxa,AllYes=FALSE)
    new.species.worms <- new.species$worms.taxonomy
    new.species.worms <- new.species.worms[!duplicated(new.species.worms),]
    
    continue <- ""
    while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Inject ",nrow(new.species.worms)," species to Database? (Y/n) ")) }
    
    if( continue == "Y") { 
      
      max.fid <- get.last.fid("Taxa")
      new.fids <- (max.fid+1):(max.fid+nrow(new.species.worms))
      new.species.worms[,1] <- new.fids
      dbWriteTable(sql, "RDB_Taxa", new.species.worms , append=TRUE, overwrite=FALSE )
      
    }
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

inject.from.file <- function(file,type,decimal,xls.sheet,ask) {
  
  if(ask==FALSE) { AllYes = TRUE}
  if(ask==TRUE) { AllYes = FALSE}
  
  internet.available()
  
  options(warn=-1)
  
  if( grepl("xlsx", file) | grepl("xls", file) ) {
    
    tryCatch( new.records <- gdata::read.xls(file,sheet=1)
              , error=function(e) { } )
    
    
  }
  
  if( grepl("csv", file) )  {
    
    tryCatch( new.records <- read.csv(file,sep=";",dec=decimal,header=TRUE)
              , error=function(e) { } )
    
  }
  
  if( ! exists("new.records") ) { errormessage("Error while reading XLS file") }
  
  new.records <- as.data.frame(new.records,stringsAsFactors = FALSE)
  
  new.records <- new.records[ !apply(is.na(new.records) | new.records == "", 1, all) , ]
  
  is.letter <- function(x) grepl("[[:alpha:]]", x)
  
  to.remove.no.geo.data <- intersect(which(is.na(new.records$Lon)),which(!is.letter(new.records$Site)))
  
  if( length(to.remove.no.geo.data) > 0) {  
    
    new.records <- new.records[-to.remove.no.geo.data,]
    
  }
  
  new.records <- new.records[ !is.na(new.records$Lon) | !is.na(new.records$Site) , ]
  new.records <- new.records[ !is.na(new.records$Lon) | new.records$Site != "" , ]
  
  new.records[,which( colnames(new.records) == "Lon" )] <- sapply(new.records$Lon,function(x) { magic.decimals(x) })
  new.records[,which( colnames(new.records) == "Lat" )] <- sapply(new.records$Lat,function(x) { magic.decimals(x) })
  
  if( nrow(new.records) == 0 ) { errormessage("Empty file") }
  
  new.records <- new.records[new.records$Species != "",]
  
  new.records <- transform(new.records, Lon = as.character(Lon))
  new.records <- transform(new.records, Lon = as.numeric(Lon))
  
  new.records <- transform(new.records, Lat = as.character(Lat))
  new.records <- transform(new.records, Lat = as.numeric(Lat))
  
  sql.fields <- get.table.structure("Record")
  
  if( sum( ! required.data.fields %in% colnames(new.records) ) > 0 ) { 
    
    fields.missing <- required.data.fields[which(! required.data.fields %in% colnames(new.records))]
    errormessage(paste0("Missing required fields in XLS file: ",paste(fields.missing,collapse=", "))) 
    
  } 
  
  if( sum(is.na(new.records$Reference)) > 0 | sum(is.na(new.records$Species)) > 0 | sum( is.na(new.records$ReferenceType)) > 0 ) { errormessage("Missing important information") } 
  
  is.factor.columns <- which(sapply(new.records, is.factor))
  
  for(f in is.factor.columns) {
    
    new.records[f] <- unlist(lapply(new.records[f], as.character))
    
  }
  
  new.records[,which(colnames(new.records) ==  "Species" )] <- trimws( new.records[,which(colnames(new.records) ==  "Species" )] ,which="right")
  new.records[,which(colnames(new.records) ==  "Species" )] <- trimws( new.records[,which(colnames(new.records) ==  "Species" )] ,which="left")
  
  # Test coordinates
  test.coords <- data.frame(Lon=as.numeric(new.records$Lon),Lat=as.numeric(new.records$Lat))
  test.coords <- test.coords[ which( ! is.na(test.coords[,1]) & ! is.na(test.coords[,2]) ) , ]
  
  if( max(test.coords[,1], na.rm=T) > 180 | min(test.coords[,1], na.rm=T) < -180 | 
      max(test.coords[,2], na.rm=T) > 90 | min(test.coords[,2], na.rm=T) < -90 ) { errormessage("Coordinates out of WGS84 projection") }
  
  if( class(new.records$Lon) != "numeric" ) { new.records[,which(colnames(new.records) ==  "Lon" )] <- as.numeric(new.records[,which(colnames(new.records) ==  "Lon" )])  }
  if( class(new.records$Lat) != "numeric" ) { new.records[,which(colnames(new.records) ==  "Lat" )] <- as.numeric(new.records[,which(colnames(new.records) ==  "Lat" )]) }
  
  # Correct species name
  
  name.of.species.to.correct <- data.frame(oldname=unique(new.records$Species), stringsAsFactors = FALSE)
  
  for(sp.n in 1:nrow(name.of.species.to.correct)) {
    
    genus <- strsplit(name.of.species.to.correct$oldname[sp.n], " ")[[1]][1]
    genus <- paste0( toupper(substring(genus, 1, 1 )) , tolower(substring(genus, 2, nchar(genus) )))
    species <- tolower(strsplit(name.of.species.to.correct$oldname[sp.n], " ")[[1]][2])
    
    if( grepl("(", species, fixed=TRUE) ) { species <- tolower(strsplit(name.of.species.to.correct$oldname[sp.n], " ")[[1]][3]) }
    
    genus <- gsub(" ", "", genus, fixed = TRUE)
    species <- gsub(" ", "", species, fixed = TRUE)
    
    if( !is.na(species) & species != "sp" & species != "spp" & species != "indet" & species != "sp." & species != "spp." & species != "indet." ) {
      
      new.records[ which(new.records$Species == name.of.species.to.correct$oldname[sp.n]) , which(colnames(new.records) == "Species") ] <- as.character(paste0(genus," ",species))
      
    } else {
      
      new.records[ which(new.records$Species == name.of.species.to.correct$oldname[sp.n]) , which(colnames(new.records) == "Species") ] <- as.character(paste0(genus))
      
    }
    
  }
  
  # Which species not in DB
  
  species.in.new.records <- unique(new.records$Species)
  # species.in.db <- get.data.from.table("Taxa","SpeciesName,AcceptedSpeciesName")
  # new.taxa.1 <- which( ! species.in.new.records %in% as.character(species.in.db$SpeciesName))
  # new.taxa.2 <- which( ! species.in.new.records %in% as.character(species.in.db$AcceptedSpeciesName))
  # new.taxa <- intersect(species.in.new.records[new.taxa.1],species.in.new.records[new.taxa.2])
  
  species.in.db <- get.data.from.table("Taxa","SpeciesName")
  new.taxa <- which( ! species.in.new.records %in% as.character(species.in.db$SpeciesName))
  new.taxa <- species.in.new.records[new.taxa]
  
  if( length(new.taxa) > 0) {
    
    new.species <- get.species.taxonomy(new.taxa,AllYes=AllYes)
    new.species.worms <- new.species$worms.taxonomy
    new.species.worms <- new.species.worms[!duplicated(new.species.worms),]
    
    # Apply changes to new species vector
    
    # taxa.changed <- new.species$original.taxa[which( ! new.species$original.taxa %in% as.character(new.species.worms$SpeciesName) )]
    
    if( nrow(new.species.worms) > 0) {
      
      for( k in 1:nrow(new.species.worms) ) {
        corrected.name.worms <- as.character(new.species$worms.taxonomy[k,]$SpeciesName)
        new.records[ which( new.records$Species == new.species$original.taxa[k] ) , which(colnames(new.records) == "Species") ] <- corrected.name.worms
        
      } }
  }
  
  taxa <- new.records$Species
  
  # DF with all species and Worms ID (DB + New)
  # 
  
  species.data <- get.data.from.table("Taxa")
  
  if( length(new.taxa) > 0) { species.data.temp <- rbind(species.data,new.species.worms)
  
  } else { species.data.temp <- species.data }
  
  new.records.worms.id <- sapply(taxa,function(x) { species.data.temp$SpeciesWormsID[ which(species.data.temp$SpeciesName == x)][1] } )
  new.records.worms.id <- as.numeric( unlist( new.records.worms.id) )
  new.records[ , which(colnames(new.records) == "Species" )] <- new.records.worms.id
  
  # Resolve Source References
  
  if( sum (! new.records$ReferenceType %in% source.tables ) > 0 ) { errormessage("Corrupted reference type") }
  
  # Test if new.records.repository in db.repository
  
  new.records.repository <- new.records$Reference[which(new.records$ReferenceType == "Repository")]
  db.repository <- unlist(get.data.from.table("Repository","fid"))
  if( sum( ! new.records.repository %in% db.repository ) > 0 ) { errormessage("Unavailable repository id") }
  
  # Resolve literature
  
  new.records[,which(colnames(new.records) == "Reference")] <- encoding.correction( as.character(new.records$Reference ) )
  new.records.literature <- unique(as.character(new.records$Reference[which(new.records$ReferenceType == "Literature")]))
  
  if( length(new.records.literature) > 0) {
    
    max.fid <- get.last.fid("Literature")
    new.fids <- as.numeric(max.fid+1):as.numeric(max.fid+length(new.records.literature))
    
    new.source.references <- data.frame(    fid = new.fids ,
                                            Type = "" ,
                                            Authors = "" ,
                                            Title = "",
                                            Year = "" ,
                                            Container = "" ,
                                            Volume = "" ,
                                            Issue = "" ,
                                            Page = "" ,
                                            FullCitation=new.records.literature,
                                            DOI = "" ,
                                            ReviewedByCrossref = 0 ,
                                            RevisionByCrossrefDate="" , stringsAsFactors = FALSE)
    
    ## ---------------
    
    new.record.f.source.reference.id <- numeric(nrow(new.records))
    
    for(e in 1:length(new.records.literature) ) {
      
      index <- which( as.character(new.records$Reference) == new.records.literature[e])
      new.record.f.source.reference.id[index] <- new.fids[e]
      
    }
    
    index.of.non.ref.type <- which( ! as.character(new.records$Reference) %in% new.records.literature )
    non.ref.type <- as.numeric(as.character(new.records$Reference[index.of.non.ref.type]))
    new.record.f.source.reference.id[index.of.non.ref.type] <- non.ref.type
    new.records$Reference <- new.record.f.source.reference.id
    
  }
  
  # Resolve multiple years
  
  new.records[ new.records$DateYear == 0 & !is.na(new.records$DateYear) , which(colnames(new.records) == "DateYear")] <- NA
  new.records[ nchar(new.records$DateYear) < 4 & !is.na(new.records$DateYear) | nchar(new.records$DateYear) > 9 & !is.na(new.records$DateYear) , which(colnames(new.records) == "DateYear")] <- NA
  test.years <- nchar(new.records$DateYear)
  if( max(test.years,na.rm=T) > 9 | min(test.years,na.rm=T) < 4 & min(test.years,na.rm=T) > 0 ) { errormessage("Corrupted data for Year") }
  
  multiple.records <- which( nchar(new.records$DateYear) > 4)
  
  if(length(multiple.records) > 0) {
    
    for(rec in 1:length(multiple.records)) {
      
      years.to.partitioning <- new.records$DateYear[multiple.records[rec]]
      
      start.years.to.partitioning <- as.numeric(substr(years.to.partitioning, 1, 4))
      stop.years.to.partitioning <- as.numeric(substr(years.to.partitioning, 6, 9))
      
      length.years.to.partitioning <- start.years.to.partitioning:stop.years.to.partitioning
      
      new.records.rep <- do.call("rbind", replicate(length(length.years.to.partitioning),  new.records[multiple.records[rec],], simplify = FALSE))
      new.records.rep[, which(colnames(new.records.rep) == "DateYear") ] <- length.years.to.partitioning
      new.records <- rbind(new.records,new.records.rep)
      
    }
    
    new.records <- new.records[-multiple.records,]
    
  }
  
  if( "Species" %in% colnames(new.records) ) { Species <- new.records$Species } else { Species <- "" }
  if( "CoordinateType" %in% colnames(new.records) ) { CoordinateType <- encoding.correction(new.records$CoordinateType) } else { CoordinateType <- "" }
  if( "CoordinateUncertaintyInMeters" %in% colnames(new.records) ) { CoordinateUncertaintyInMeters <- new.records$CoordinateUncertaintyInMeters } else { CoordinateUncertaintyInMeters <- "" }
  if( "Country" %in% colnames(new.records) ) { Country <- encoding.correction(new.records$Country) } else { Country <- "" }
  if( "Site" %in% colnames(new.records) ) { Site <- encoding.correction(new.records$Site) } else { Site <- "" }
  
  if( "Depth" %in% colnames(new.records) ) { Depth <- new.records$Depth } else { Depth <- NA }
  if( "DepthRangeMin" %in% colnames(new.records) ) { DepthRangeMin <- new.records$DepthRangeMin } else { DepthRangeMin <- NA }
  if( "DepthRangeMax" %in% colnames(new.records) ) { DepthRangeMax <- new.records$DepthRangeMax } else { DepthRangeMax <- NA }
  
  depth.range <- do.call("rbind" , lapply(new.records$Depth,function(x) { magic.depths(x) } ))
  with.depth.range.data <- which(!is.na(depth.range[,1]))
  
  if( length(with.depth.range.data) > 0 ) {
    
    DepthRangeMin[ with.depth.range.data ] <- depth.range[,1]
    DepthRangeMax[ with.depth.range.data ] <- depth.range[,2]
    Depth[ with.depth.range.data ] <- NA
    
  }
  
  if( "DateYear" %in% colnames(new.records) ) { DateYear <- new.records$DateYear } else { DateYear <- NA }
  if( "DateMonth" %in% colnames(new.records) ) { DateMonth <- new.records$DateMonth } else { DateMonth <- NA }
  if( "DateDay" %in% colnames(new.records) ) { DateDay <- new.records$DateDay } else { DateDay <- NA }
  
  if( "AbundanceID" %in% colnames(new.records) ) { AbundanceID <- new.records$AbundanceID } else { AbundanceID <- "" }
  if( "SampleCollectionID" %in% colnames(new.records) ) { SampleCollectionID <- new.records$SampleCollectionID } else { SampleCollectionID <- "" }
  if( "RecordUserID" %in% colnames(new.records) ) { RecordUserID <- new.records$RecordUserID } else { RecordUserID <- UserIDSession }
  
  RecordDateYear <- as.numeric(format(Sys.time(), "%Y"))
  RecordDateMonth <- as.numeric(format(Sys.time(), "%m"))
  RecordDateDay <- as.numeric(format(Sys.time(), "%d"))
  
  if( "RecordPublic" %in% colnames(new.records) ) { RecordPublic <- new.records$RecordPublic } else { RecordPublic <- "" }
  if( "RecordNotes" %in% colnames(new.records) ) { RecordNotes <- encoding.correction(new.records$RecordNotes) } else { RecordNotes <- "" }
  if( "OriginalSource" %in% colnames(new.records) ) { OriginalSource <- encoding.correction(new.records$OriginalSource) } else { OriginalSource <- "" }
  if( "OriginalSourceID" %in% colnames(new.records) ) { OriginalSourceID <- encoding.correction(new.records$OriginalSourceID) } else { OriginalSourceID <- "" }
  
  new.records.f <- data.frame(  fid = NA,
                                SpeciesID = as.numeric(as.character(new.records$Species)) ,
                                Lon = as.numeric(as.character(new.records$Lon)) ,
                                Lat = as.numeric(as.character(new.records$Lat)) ,
                                CoordinateType = as.character(CoordinateType) ,
                                CoordinateUncertaintyInMeters = as.numeric(as.character(CoordinateUncertaintyInMeters)) ,
                                Country = as.character(Country) ,
                                Site = as.character(Site) ,
                                Depth = as.numeric(as.character(Depth)) ,
                                DepthRangeMin = as.numeric(as.character(DepthRangeMin)) ,
                                DepthRangeMax = as.numeric(as.character(DepthRangeMax)) ,
                                DateYear = as.numeric(as.character(DateYear)) ,
                                DateMonth = as.numeric(as.character(DateMonth)) ,
                                DateDay = as.numeric(as.character(DateDay)) ,
                                AbundanceID = as.numeric(as.character( rep(AbundanceID,nrow(new.records)) )) ,
                                SampleCollectionID = as.numeric(as.character( rep(SampleCollectionID,nrow(new.records)) )) ,
                                SourceRefType = as.character(new.records$ReferenceType) ,
                                SourceRefID = as.character(new.records$Reference) ,
                                RecordUserID = as.numeric(as.character(RecordUserID)) ,
                                RecordDateYear = as.numeric(as.character( rep(RecordDateYear,nrow(new.records)))) ,
                                RecordDateMonth = as.numeric(as.character(rep(RecordDateMonth,nrow(new.records)) )) ,
                                RecordDateDay = as.numeric(as.character( rep(RecordDateDay,nrow(new.records)))) ,
                                RecordPublic = as.numeric(as.character(RecordPublic)) ,
                                RecordNotes = as.character(RecordNotes) , 
                                OriginalSource = as.character(OriginalSource) ,
                                Flag = "0",
                                OriginalSourceID = as.character(OriginalSourceID) ,
                                FlagOnLand = "0",
                                FlagLightSuitable = "0",
                                FlagO2Suitable = "0",
                                FlagEcoRegion = "0",
                                stringsAsFactors=FALSE)
  
  if( ! exists("new.source.references") ) { new.source.references <- NULL }
  if( ! exists("new.species") ) { new.species <- NULL } else { new.species <- new.species.worms }
  
  if (sum(!is.na(new.records.f$AbundanceID)) > 0 ) { errormessage("Uncoded situation, AbundanceID") }
  if (sum(!is.na(new.records.f$SampleCollectionID)) > 0 ) { errormessage("Uncoded situation, SampleCollectionID") }
  
  cat('\014')
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------')
  cat('\n')
  cat('Number of records:', nrow(new.records.f))
  cat('\n')
  cat('Number of species:', length(unique(new.records.f$SpeciesID)))
  cat('\n')
  cat('\n')
  
  continue <- ""
  while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Show new records? (Y/n) ")) }
  if( continue == "Y") { 
    
    to.show <- new.records.f[1:50,]
    to.show <- to.show[complete.cases(to.show[,2]),]
    print(to.show,row.names=FALSE)
    
  }
  
  if( !is.null(new.species)) {
    
    continue <- ""
    while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Show new species? (Y/n) ")) }
    if( continue == "Y") { 
      
      to.show <- new.species[1:50,]
      to.show <- to.show[complete.cases(to.show[,2]),]
      print(to.show,row.names=FALSE)
      
    }
  }
  
  
  
  
  continue <- ""
  while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Plot new records? (Y/n) ")) }
  if( continue == "Y") { 
    
    return.plot <- plot.coords.google(new.records.f[,3:4],new.records.f[,1],4,"black")
    print(return.plot) 
    
  }
  
  continue <- ""
  while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Inject new records to Database? (Y/n) ")) }
  
  if( continue == "Y") { 
    
    if( sum( ! is.na(new.records.f$CoordinateUncertaintyInMeters) ) > 0 ) {
      
      continue.2 <- ""
      while( continue.2 !="Y" & continue.2 !="n" ) { continue.2 <- readline(paste0("Some records have information about coordinate uncertainty. Do you which to trim records? (Y/n) ")) }
      
      if( continue.2 == "Y") { 
        
        threhold.meters <- readline(paste0("Give a threshold in meters: "))
        new.records.f <- new.records.f[ which(is.na(new.records.f$CoordinateUncertaintyInMeters) | new.records.f$CoordinateUncertaintyInMeters <= as.numeric(threhold.meters)) , ]
        
      }
    }
    
    inject.to.sql(new.records.f,new.species,new.source.references) 
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

inject.to.sql <- function(new.records,new.taxa,new.reference) {
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  
  if( ! is.null(new.records)  ) {
    
    max.fid <- get.last.fid("Record")
    new.fids <- (max.fid+1):(max.fid+nrow(new.records))
    new.records[,1] <- new.fids
    dbWriteTable(sql, "MDB_GeoData", new.records , append=TRUE, overwrite=FALSE )
    
  }
  
  if( ! is.null(new.taxa) ) {
    
    max.fid <- get.last.fid("Taxa")
    new.fids <- (max.fid+1):(max.fid+nrow(new.taxa))
    new.taxa[,1] <- new.fids
    dbWriteTable(sql, "RDB_Taxa", new.taxa , append=TRUE, overwrite=FALSE )
    
  }
  
  if( ! is.null(new.reference) ) {
    
    dbWriteTable(sql, "RDB_Ref_Literature", new.reference , append=TRUE, overwrite=FALSE )
    
  }
  
  dbDisconnect(sql)
  updated.version.control()
  
}

## ------------------------------------------------------------------------------------------------------------------

duplicated.records <- function() {
  
  duplicated.db <- get.data.from.table("Record","fid, SpeciesID, Lon, Lat, DateYear, DateMonth, SourceRefID, SourceRefType, Site")
  duplicated.db <- duplicated.db[!is.na(duplicated.db$DateYear) & !is.na(duplicated.db$DateMonth),]
  
  # duplicated.db <- get.data.from.table("Record","fid, SpeciesID, Lon, Lat, DateYear")
  # duplicated.db <- duplicated.db[!is.na(duplicated.db$DateYear),]
  
  duplicated.db.fid <- duplicated.db$fid
  duplicated <- which(duplicated(duplicated.db[,-1]))
  
  if( length(duplicated) > 0 ) {
    
    cat('\014')
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat( length(duplicated), "duplicated records")
    cat('\n')
    cat('\n')
    
    fid.to.delete <- as.numeric(duplicated.db.fid[duplicated])
    delete.records("Record",fid.to.delete)
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

populate.missing.coordinates <- function() {
  
  missing.data.db <- get.data.from.table("Record","fid, Lon, Lat, Country, Site")
  missing.data <- which( missing.data.db$Lon == 0 | is.na(missing.data.db$Lat) | missing.data.db$Lat == 0 | is.na(missing.data.db$Lon) )
  
  if( length(missing.data) > 0) {
    
    internet.available()
    
    nas <- numeric(0)
    
    for( i in 1:length(missing.data)) {
      
      cat('\014')
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat('Geocoding',i,'out of',length(missing.data) , 'record(s) with no coordinate information')
      cat('\n')
      cat('\n')
      
      record <- missing.data.db[missing.data[i],]
      
      if( nchar(paste0(record$Site," ",record$Country)) > 1 ) {
        
        new.coords <- Geocode.address(paste0(record$Site," ",record$Country))
        
        if( is.na(new.coords[1]) ) { new.coords <- Geocode.address(record$Site) }
        
      }
      
      if( ! NA %in% new.coords ) {
        
        sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
        dbExecute(sql,paste0("UPDATE MDB_GeoData SET Lon = ",new.coords[2]," WHERE fid = ", record$fid)) 
        dbExecute(sql,paste0("UPDATE MDB_GeoData SET Lat = ",new.coords[1]," WHERE fid = ", record$fid)) 
        dbExecute(sql,paste0("UPDATE MDB_GeoData SET CoordinateType = 'Geocoded' WHERE fid = ", record$fid)) 
        dbDisconnect(sql)
        updated.version.control()
        
      }
      
      if( NA %in% new.coords ) { nas <- c(nas,1) }
      
    }
    
    if(length(nas) > 0 ) {  
      
      cat(length(nas) , 'record(s) for which no coordinate information were retrieved')
      
    }
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

summary.all <- function(type,exclude,level) {
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
  n.records <- as.numeric(dbGetQuery(sql,paste0("SELECT Count(fid) FROM MDB_GeoData")))
  n.taxa <- as.numeric(dbGetQuery(sql,paste0("SELECT Count(fid) FROM RDB_Taxa")))
  n.references <- as.numeric(dbGetQuery(sql,paste0("SELECT Count(fid) FROM RDB_Ref_Literature")))
  n.repository <- as.numeric(dbGetQuery(sql,paste0("SELECT Count(fid) FROM RDB_Ref_Repository")))
  n.users <- as.numeric( dbGetQuery(sql,paste0("SELECT Count(fid) FROM RDB_Users")))
  dbDisconnect(sql)
  
  summary.db <- get.data.from.table("Taxa","*")
  
  if(  type == "original" ) {
    
    SpeciesWormsID <- summary.db[ , ]$SpeciesWormsID
    SpeciesName <- summary.db[ , ]$SpeciesName
    
    results <- data.frame()
    
    for(sp in 1:length(SpeciesName)) {
      
      results.t <- summary.taxa(SpeciesName[sp],type,plot=FALSE)
      results.t.summary <- summary.db[which(summary.db$SpeciesName == results.t$Taxa),c("TaxKingdom","TaxPhylum","TaxClass","TaxOrder","TaxFamily","TaxGenus")]
      
      if( sum(sapply(exclude,function(x) { x %in% results.t.summary })) == 0 ) {
        
        results.t <- data.frame(results.t,results.t.summary, stringsAsFactors = FALSE)
        results <- rbind(results,results.t)
        
      }
      
    }
    
  }
  
  if(  type == "accepted" ) {
    
    SpeciesWormsID <- summary.db[summary.db$Status == "accepted", ]$SpeciesWormsID
    SpeciesName <- summary.db[summary.db$Status == "accepted", ]$AcceptedSpeciesName
    
    results <- data.frame()
    
    for(sp in 1:length(SpeciesName)) {
      
      results.t <- summary.taxa(SpeciesName[sp],type,plot=FALSE)
      results.t <- results.t[which(results.t$Status == "accepted" | results.t$Status == "unaccepted"),]
      results.t.sum <- sum(results.t$Records)
      results.t <- results.t[results.t$Status == "accepted",]
      results.t[,"Records"] <- results.t.sum
      
      results.t.summary <- summary.db[which(summary.db$SpeciesName == as.character(results.t$Taxa)),c("TaxKingdom","TaxPhylum","TaxClass","TaxOrder","TaxFamily","TaxGenus")]
      results.t <- data.frame(results.t,results.t.summary, stringsAsFactors = FALSE)
      
      if( level == "species" & sapply(strsplit(as.character(results.t$Taxa), " "), length) < 2 ) { 
        
        results.t <- NULL
        
      }
      
      if( sum(sapply(exclude,function(x) { x %in% results.t.summary })) == 0 ) {
        
        results <- rbind(results,results.t)
        
      }
    }
    
  }
  
  results <- results[ order( results[,"TaxKingdom"] , results[,"TaxPhylum"] , results[,"TaxClass"] , results[,"TaxOrder"] , results[,"TaxFamily"] , results[,"TaxGenus"] , results[,"Taxa"] ),]
  
  cat('\014')
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------')
  cat('\n')
  cat('\n')
  cat('Number of taxa: ', length(SpeciesName))
  cat('\n')
  cat('Number of records: ', sum(results$Records))
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------------------')
  cat('\n')
  cat('\n')
  cat('Sources (n)')
  cat('\n')
  cat('Literature: ', n.references)
  cat('\n')
  cat('Repository: ', n.repository)
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------------------')
  cat('\n')
  cat('\n')
  cat('Users: ',n.users)
  cat('\n')
  cat('\n')
  cat('DB version (Created): ',get.data.from.table("Version","Created")$Created)
  cat('\n')
  cat('DB version (Updated): ',get.data.from.table("Version","Updated")$Updated)
  cat('\n')
  cat('----------------------------------------------------------------------------------')
  cat('\n')
  cat('\n')
  
  return(results)
  
}

## ------------------------------------------------------------------------------------------------------------------

summary.taxa <- function(taxa,type,plot) {
  
  if( missing(plot)) { plot <- FALSE }
  if( missing(taxa)) { errormessage("missing taxa") }
  if( missing(type)) { type <- "original" }
  
  if(  type == "original" ) {
    
    summary.db <- get.data.from.table("Taxa","fid, SpeciesName, SpeciesWormsID, Status")
    
  }
  
  if(  type == "accepted" ) {
    
    summary.db <- get.data.from.table("Taxa","fid, AcceptedSpeciesName, SpeciesWormsID, Status")
    colnames(summary.db) <- c("fid","SpeciesName","SpeciesWormsID","Status")
    
  }
  
  taxa.in.db <- which(grepl(taxa,summary.db$SpeciesName, fixed = TRUE) & nchar(taxa) == nchar(summary.db$SpeciesName))
  
  summary.db <- summary.db[taxa.in.db,]
  
  if(nrow(summary.db) == 0) { errormessage("Taxa not present in the database.") }
  
  summary.t <- data.frame()
  
  for( t in 1:length(taxa.in.db) ) {
    
    name <- summary.db[t,]$SpeciesName
    status <- summary.db[t,]$Status
    wormsid <- summary.db[t,]$SpeciesWormsID
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    records <- as.numeric(dbGetQuery(sql,paste0("SELECT Count(fid) FROM MDB_GeoData WHERE SpeciesID= ",wormsid)))
    dbDisconnect(sql)
    summary.t <- rbind(summary.t,data.frame(Taxa=as.character(name),
                                            Status=as.character(status),
                                            WormsID=wormsid,
                                            Records=records))
    
  }
  
  summary.t[,2] <- as.character(summary.t[,2])
  summary.t[summary.t$Status == "",2] <- as.character("invalid")
  summary.t <- summary.t[sort( as.character(summary.t$Taxa), index.return = TRUE)$ix,]
  
  if( plot ) {
    
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    print(summary.t, row.names = FALSE)
    cat('\n')
    cat('\n')
    
    continue <- ""
    while( continue !="Y" & continue !="n" ) {  continue <- readline(paste0("Plot records? (Y/n) ")) }
    
    if( continue == "Y") { 
      
      if(nrow(summary.t) > 1) {
        
        continue.2 <- ""
        while( continue.2 !="Y" & continue.2 !="n" ) {  continue.2 <- readline(paste0("Plot species by species? (Y/n) ")) }
        
      } else { continue.2 <- "n" }
      
      if( continue.2 == "n") { 
        
        fids.to.plot <- data.frame()
        
        for( m in 1:nrow(summary.t)) {
          
          fids.to.plot <- rbind(fids.to.plot,
                                get.data.from.table("Record","fid","SpeciesID",summary.t$WormsID[m]) )
          
        }
        
        fids.to.plot <- as.numeric(unlist(fids.to.plot))
        print( plot.fid.google(fids.to.plot,4,"black") )
      }
      
      if( continue.2 == "Y") { 
        
        species.in.database <- unlist(summary.t$WormsID)
        
        for(p in 1:length(species.in.database)) {
          
          cat('\n')
          
          status <- get.data.from.table("Taxa","Status","SpeciesWormsID",as.numeric(species.in.database[p]))
          original.species.name <- get.data.from.table("Taxa","SpeciesName","SpeciesWormsID",as.numeric(species.in.database[p]))
          accepted.species.name <- get.data.from.table("Taxa","AcceptedSpeciesName","SpeciesWormsID",as.numeric(species.in.database[p]))
          
          if( status == "accepted" ) { 
            
            cat(paste0(original.species.name))
            
          }
          if( status != "accepted" ) { 
            
            cat(paste0(original.species.name , " with unaccepted status (Accepted: ",accepted.species.name,")" ))
            
          }
          
          print( plot.fid.google( get.data.from.table("Record","fid","SpeciesID",species.in.database[p]),4,"black") )
          continue.3 <- readline()
          
        }
        
      }
      
    }
    
  }
  
  return(summary.t)
}

## ------------------------------------------------------------------------------------------------------------------

delete.records <- function(type,fid,overpass) {
  
  if( missing(fid)) { errormessage("no fid number introduced.") }
  if( missing(overpass)) { overpass <- FALSE }
  
  if( missing(type) | ! type %in% source.tables ) { errormessage("Type of source is missing") }
  
  if(type == "Taxa") { table <- "RDB_Taxa" }
  if(type == "Literature") { table <- "RDB_Ref_Literature" }
  if(type == "User") { table <- "RDB_Users" }
  if(type == "Repository") { table <- "RDB_Ref_Repository" }
  if(type == "Record") { table <- "MDB_GeoData" }
  if(type == "Abundance") { table <- "RDB_Abundance" }
  if(type == "Collection") { table <- "RDB_Sample_Collection" }
  
  if( length(fid) > 0 & overpass) {  
    
    fid <- paste0("'", fid, "'", collapse=", ")
    fid <- paste0("(", fid, ")")
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbExecute(sql,paste0("DELETE FROM ",table," WHERE fid IN ", fid)) 
    dbDisconnect(sql)
    updated.version.control()
    
  }
  
  if( length(fid) > 0 & ! overpass ) {
    
    cat('\n')
    
    x <- ""
    
    while( x != "Y" & x != "n" ) { 
      x <- readline(paste0("Delete ",length(fid)," records from ",type," (Y/n): "))
    }
    
    cat('\n')
    cat('\n')
    
    if ( x == "Y" ) { 
      
      fid <- paste0("'", fid, "'", collapse=", ")
      fid <- paste0("(", fid, ")")
      
      sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
      dbExecute(sql,paste0("DELETE FROM ",table," WHERE fid IN ", fid)) 
      dbDisconnect(sql)
      updated.version.control()
      
    }
  } }

## ------------------------------------------------------------------------------------------------------------------

get.external.db <- function(taxaname,taxa.id,add.obis,add.gbif) {
  
  internet.available()
  
  if( missing(taxaname)) { errormessage("no taxa (Worms name) introduced.") }
  if( missing(taxa.id)) { errormessage("no taxa (Worms ID) introduced.") }
  if( missing(add.obis)) { add.obis <- FALSE }
  if( missing(add.gbif)) { add.gbif <- FALSE }
  
  temp.results <- data.frame()
  
  
  if( add.obis ) { 
    
    tryCatch( my_occs <- occurrence(scientificname = taxaname ) , error=function(e) { error <- TRUE })
    
    if( exists("my_occs") ) { if( nrow(my_occs) == 0 ) { my_occs <- data.frame() } }
    
    if( ! exists("my_occs") ) { my_occs <- data.frame() }
    
    if( nrow(my_occs) > 0) {
      
      my_occs <- subset(my_occs, my_occs$decimalLongitude !=0 & my_occs$decimalLatitude !=0)
      
      if( ! is.null(my_occs$rightsHolder) ) {
        
        my_occs <- my_occs[ is.na(my_occs$rightsHolder) ,]
        
      }
      
      if( nrow(my_occs) > 0) {
        
        OriginalSource <- ""
        if( is.null(my_occs$occurrenceID)) { OriginalSourceID <- ""} else { OriginalSourceID <- my_occs$occurrenceID }
        if( is.null(my_occs$depth)) { depths <- NA} else { depths <- my_occs$depth }
        if( is.null(my_occs$minimumDepthInMeters)) { minimumDepthInMeters <- ""} else { minimumDepthInMeters <- my_occs$minimumDepthInMeters }
        if( is.null(my_occs$maximumDepthInMeters)) { maximumDepthInMeters <- ""} else { maximumDepthInMeters <- my_occs$maximumDepthInMeters }
        if( is.null(my_occs$coordinateUncertaintyInMeters)) { coordinateUncertaintyInMeters <- ""} else { coordinateUncertaintyInMeters <- my_occs$coordinateUncertaintyInMeters }
        if( is.null(my_occs$occurrenceRemarks)) { Notes <- ""} else { Notes <- my_occs$occurrenceRemarks }
        if( is.null(my_occs$locality)) { SiteName <- ""} else { SiteName <- my_occs$locality }
        if( is.null(my_occs$yearcollected)) { years <- ""} else { years <- my_occs$yearcollected }
        
        temp.results <- rbind(temp.results,
                              data.frame(Species= taxaname,
                                         SpeciesID= taxa.id,
                                         Lon=my_occs$decimalLongitude,
                                         Lat=my_occs$decimalLatitude,
                                         CoordinateUncertaintyInMeters=coordinateUncertaintyInMeters,
                                         SiteName=SiteName,
                                         Country=rep("",nrow(my_occs)),
                                         Depth=depths,
                                         DepthRangeMin=minimumDepthInMeters,
                                         DepthRangeMax=maximumDepthInMeters,
                                         DateYear=years,
                                         OriginalSource=OriginalSource,
                                         Flag = rep("0",nrow(my_occs)),
                                         OriginalSourceID=OriginalSourceID,
                                         FlagOnLand = rep("0",nrow(my_occs)),
                                         FlagLightSuitable = rep("0",nrow(my_occs)),
                                         FlagO2Suitable = rep("0",nrow(my_occs)),
                                         FlagEcoRegion = rep("0",nrow(my_occs)),
                                         Notes=Notes,
                                         SourceRefType="Repository" ,
                                         SourceRefID=2 , stringsAsFactors=FALSE) )
        
      }
      
    }
    
    if( exists("my_occs") ) { rm(my_occs) }
    
  }
  
  if( add.gbif ) {  
    
    tryCatch( my_occs <- gbif(strsplit(as.character(taxaname), " ")[[1]][1], strsplit(as.character(taxaname), " ")[[1]][2], geo=T, removeZeros=T ) , error=function(e) { error <- TRUE })
    
    if( exists("my_occs") ) { if( is.null(my_occs) ) { my_occs <- data.frame() } }
    
    if( ! exists("my_occs") ) { my_occs <- data.frame() }
    
    if( ! is.null(my_occs$lat) & nrow(my_occs) > 0 ) {
      
      my_occs <- subset(my_occs, lat !=0 & lon !=0)
      
      if( !is.null(my_occs$accessRights) ) {
        
        my_occs <- my_occs[is.na(my_occs$accessRights) | my_occs$accessRights == "Free usage",]
        
      }
      
      if( nrow(my_occs) > 0 ) {
        
        if( is.null(my_occs$depth)) { depths <- NA} else { depths <- my_occs$depth }
        if( is.null(my_occs$year)) { years <- NA} else { years <- my_occs$year }
        if( is.null(my_occs$coordinateUncertaintyInMeters)) { coordinateUncertaintyInMeters <- ""} else { coordinateUncertaintyInMeters <- my_occs$coordinateUncertaintyInMeters }
        if( is.null(my_occs$gbifID)) { OriginalSource <- ""} else { OriginalSource <- paste0("https://www.gbif.org/occurrence/",my_occs$gbifID) }
        if( is.null(my_occs$gbifID)) { OriginalSourceID <- ""} else { OriginalSourceID <- my_occs$gbifID }
        if( is.null(my_occs$occurrenceRemarks)) { Notes <- ""} else { Notes <- my_occs$occurrenceRemarks }
        if( is.null(my_occs$locality)) { SiteName <- ""} else { SiteName <- my_occs$locality }
        if( is.null(my_occs$fullCountry)) { Country <- ""} else { Country <- my_occs$fullCountry }
        if( is.null(my_occs$year)) { years <- ""} else { years <- my_occs$year }
        
        temp.results <- rbind(temp.results,
                              data.frame(Species= taxaname,
                                         SpeciesID= taxa.id,
                                         Lon=my_occs$lon,
                                         Lat=my_occs$lat,
                                         CoordinateUncertaintyInMeters=coordinateUncertaintyInMeters,
                                         SiteName=SiteName,
                                         Country=Country,
                                         Depth= depths ,
                                         DepthRangeMin=rep(NA,nrow(my_occs)) ,
                                         DepthRangeMax=rep(NA,nrow(my_occs)) ,
                                         DateYear=years,
                                         OriginalSource=OriginalSource,
                                         Flag = rep("0",nrow(my_occs)),
                                         OriginalSourceID=OriginalSourceID,
                                         FlagOnLand = rep("0",nrow(my_occs)),
                                         FlagLightSuitable = rep("0",nrow(my_occs)),
                                         FlagO2Suitable = rep("0",nrow(my_occs)),
                                         FlagEcoRegion = rep("0",nrow(my_occs)),
                                         Notes=Notes,
                                         SourceRefType="Repository" ,
                                         SourceRefID=1 , stringsAsFactors=FALSE))
        
      }
      
      
    }
    
    if( exists("my_occs") ) { rm(my_occs) }
    
  }
  
  return(temp.results)
  
}

## ------------------------------------------------------------------------------------------------------------------

to.first.capital <- function(string) {
  
  string.1 <- toupper(substr(string, 1, 1))
  string.2 <- tolower(substr(string, 2, nchar(string)))
  
  return(paste0(string.1,string.2))  
  
}

## ------------------------------------------------------------------------------------------------------------------

export.records <- function(taxa.id,type,save.to.file,file.name,restrict.region) {
  
  if( missing(restrict.region) ) { restrict.region <- c(-180,180,-90,90) }
  if( is.null(restrict.region) ) { restrict.region <- c(-180,180,-90,90) }
  if( missing(type)) { type <- "original" }
  
  if(  type == "accepted" ) {
    
    taxa.id <- get.data.from.table("Taxa","SpeciesWormsID, AcceptedWormsID","AcceptedWormsID",taxa.id)$SpeciesWormsID   
    
  }
  
  if( length(taxa.id) > 0 ) { 
    
    taxa.id <- unique(taxa.id)
    temp.table <- get.data.from.table("Taxa","SpeciesWormsID, SpeciesName , Authority, Status , TaxKingdom , TaxPhylum , TaxClass , TaxOrder , TaxFamily , TaxGenus , AcceptedWormsID , AcceptedSpeciesName")    
    temp.table <- temp.table[temp.table$SpeciesWormsID %in% taxa.id ,]
    
    records.to.export <- data.frame()
    
    for( t in 1:length(taxa.id)) {   
      
      taxonomic.info <- temp.table[t,]
      taxa.to.import <- taxonomic.info$SpeciesWormsID
      
      contents <- get.data.from.table("Record","fid , Lon , Lat , CoordinateUncertaintyInMeters ,Country , Site , Depth , DepthRangeMin , DepthRangeMax , DateYear , DateMonth , DateDay , SourceRefType , SourceRefID , OriginalSource , OriginalSourceID ,  Flag , FlagOnLand , FlagLightSuitable , FlagO2Suitable , FlagEcoRegion , FlagEcoRegion2 , FlagEcoRegion3 , RecordNotes","SpeciesID",taxa.to.import)
      
      if(nrow(contents) > 0) {
        
        results.temp <-  data.frame(taxonomic.info, 
                                    contents)
        
        colnames(results.temp) <- c( "aphiaID" , "speciesName" , "authority" ,"status" , "taxKingdom", "taxPhylum", "taxClass", "taxOrder", "taxFamily", "taxGenus",  "acceptedAphiaID" , "acceptedSpeciesName" , "marineForestsID" ,   "decimalLongitude" , "decimalLatitude" , "coordinateUncertaintyInMeters" , "country" , "locality" , "depth" , "minimumDepthInMeters" , "maximumDepthInMeters" , "year" , "month" , "day" , "originalSourceType" , "SourceRefID" , "originalSource" , "originalSourceDOI" , "flagHumanCuratedDistribution" , "flagMachineOnLand" , "flagMachineSuitableLightBottom" , "flagMachineSuitableO2Bottom" , "flagMachineEcoRegionL1", "flagMachineEcoRegionL2", "flagMachineEcoRegionL3" , "RecordNotes")
        records.to.export <- rbind(records.to.export,results.temp)
        
      }
      
    }
    
    records.to.export <- records.to.export[  records.to.export$decimalLongitude >= restrict.region[1] & 
                                               records.to.export$decimalLongitude <= restrict.region[2] & 
                                               records.to.export$decimalLatitude >= restrict.region[3] &
                                               records.to.export$decimalLatitude <= restrict.region[4] , ]
    
    #  if(nrow(records.to.export) == 0) { errormessage("No records in restricted region")  }
    
    if( "Literature" %in% unique(records.to.export$originalSourceType) ) {
      
      for( l in which(records.to.export$originalSourceType == "Literature" ) ) {
        
        literature.included.id <- records.to.export[l,]$SourceRefID
        literature.included <- get.data.from.table("Literature","*","fid",literature.included.id)$FullCitation
        literature.included.doi <- get.data.from.table("Literature","*","fid",literature.included.id)$DOI
        literature.included.type <- get.data.from.table("Literature","*","fid",literature.included.id)$Type
        
        if(is.null(literature.included.type) | is.na(literature.included.type) | literature.included.type == "") { literature.included.type <- "journal-article" }
        
        if(length(literature.included.doi) == 0) { literature.included.doi <- "" }
        
        records.to.export[l, which(colnames(records.to.export) == "originalSource")] <- literature.included
        records.to.export[l, which(colnames(records.to.export) == "originalSourceDOI")] <- literature.included.doi
        records.to.export[l, which(colnames(records.to.export) == "originalSourceType")] <- literature.included.type
        
      }
      
    }
    
    if( "Repository" %in% unique(records.to.export$originalSourceType) ) {
      
      for( l in which(records.to.export$originalSourceType == "Repository" ) ) {
        
        repository.included <- records.to.export[l,]$SourceRefID
        repository.included <- get.data.from.table("Repository","*","fid",repository.included)
        repository.type <- repository.included$RepositoryType
        
        repository.included.name <- character()
        
        for( field in c("RepositoryShortName","RepositoryName","Institution","Responsible","ResponsibleContact","URL")) {
          
          if(repository.included[,field] != "") { repository.included.name <- paste0(repository.included.name,repository.included[,field],". ") }
          
        }
        
        repository.included.name <- substr(repository.included.name, 1, nchar(repository.included.name)-2)
        
        if( ! is.na(records.to.export[l,]$originalSourceDOI) ) { 
          
          if( records.to.export[l,]$originalSourceDOI != "" ) {
            
            repository.included.name <- paste0(repository.included.name," (Original record identifier: " , records.to.export[l,]$originalSourceDOI , ")" )
            
          }
        }
        
        if( ! is.na(records.to.export[l,]$originalSource) ) { 
          
          if( records.to.export[l,]$originalSource != "" ) {
            
            repository.included.doi <- records.to.export[l,]$originalSource 
            
          } 
          
        }
        
        if( ! is.na(records.to.export[l,]$originalSource) ) { 
          
          if( records.to.export[l,]$originalSource == "" ) {
            
            repository.included.doi <- repository.included$URL
            
          } 
          
        }
        
        if( is.na(records.to.export[l,]$originalSource) ) { 
          
          repository.included.doi <- repository.included$URL
          
        }
        
        records.to.export[l, which(colnames(records.to.export) == "originalSource")] <- repository.included.name
        records.to.export[l, which(colnames(records.to.export) == "originalSourceDOI")] <- repository.included.doi
        records.to.export[l, which(colnames(records.to.export) == "originalSourceType")] <- repository.type
        
      }
      
    }
    
    records.to.export <- records.to.export[,-which(colnames(records.to.export) %in% c("SourceRefID"))]
    
    records.to.export[,which(colnames(records.to.export) == "depth")] <- abs(records.to.export[,which(colnames(records.to.export) == "depth")])
    records.to.export[,which(colnames(records.to.export) == "minimumDepthInMeters")] <- abs(records.to.export[,which(colnames(records.to.export) == "minimumDepthInMeters")])
    records.to.export[,which(colnames(records.to.export) == "maximumDepthInMeters")] <- abs(records.to.export[,which(colnames(records.to.export) == "maximumDepthInMeters")])
    
  }
  
  records.to.export <- data.frame(lapply(records.to.export, function(x) { gsub(";", ",", x) }))
  
  if(save.to.file) {
    
    records.to.export <- data.frame(lapply(records.to.export, function(x) { gsub(";", ",", x) }))
    write.table(x=records.to.export,file=paste0(main.directory,"Data/Export/",file.name,".txt"),quote=FALSE,sep=";",dec=".",row.names=FALSE,na = "")
    
  }
  
  if(!save.to.file) {
    
    return(records.to.export)
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

plot.fid.google <- function(fid,radius,color,zoom) {
  
  zoom.define <- TRUE
  
  if( missing(radius)) { radius <- 2 }
  if( missing(color)) { color <- "black" }
  if( missing(zoom)) { zoom.define <- FALSE }
  
  fid <- as.numeric(unlist(fid))
  
  data.to.plot <- get.data.from.table("Record","*","fid", fid )
  
  # data.to.plot[,1] <- as.numeric(as.character(data.to.plot[,1]))
  # data.to.plot[,2] <- as.numeric(as.character(data.to.plot[,2]))
  # data.to.plot[,3] <- as.numeric(as.character(data.to.plot[,3]))
  
  if( max(data.to.plot[,"Lon"], na.rm=T) > 180 | min(data.to.plot[,"Lon"], na.rm=T) < -180 | 
      max(data.to.plot[,"Lat"], na.rm=T) > 90 | min(data.to.plot[,"Lat"], na.rm=T) < -90 ) { errormessage("Coordinates out of WGS84 projection") }
  
  na.records <- which(is.na(data.to.plot[,"Lon"]) )
  
  if( length(na.records) > 0 ) {
    
    cat('\n')
    cat('\n')
    cat(paste0(length(na.records)," records with no coordinate information"))
    cat('\n')
    cat('\n')
    
    data.to.plot <- data.to.plot[-na.records,]
    
  }
  
  temp.record.fid <- data.to.plot$fid
  temp.record.taxa <- get.data.from.table("Taxa","SpeciesWormsID, SpeciesName, Status","SpeciesWormsID",data.to.plot$SpeciesID)
  temp.record.taxa <- as.data.frame(t(sapply(data.to.plot$SpeciesID , function(x) {  temp.record.taxa[temp.record.taxa$SpeciesWormsID == x,c("SpeciesName","Status")]  } )))
  
  species.name <- unlist(temp.record.taxa$SpeciesName)
  species.status <- unlist(temp.record.taxa$Status)
  species.wormsid <- data.to.plot$SpeciesID
  
  temp.record.site <- data.to.plot$Site
  temp.record.country <- data.to.plot$Country
  temp.record.year <- data.to.plot$DateYear
  temp.record.depth <- data.to.plot$Depth
  temp.record.reference <- data.to.plot$SourceRefType
  temp.record.reference.id <- data.to.plot$SourceRefID
  
  # temp.record.site <- as.character(unlist(get.data.from.table("Record","Site","fid", fid  )))
  # temp.record.country <- as.character(unlist(get.data.from.table("Record","Country","fid", fid  )))
  # temp.record.year <- as.numeric(unlist(get.data.from.table("Record","DateYear","fid", fid  )))
  # temp.record.depth <- as.numeric(unlist(get.data.from.table("Record","Depth","fid", fid  )))
  # temp.record.reference <- as.character(unlist(get.data.from.table("Record","SourceRefType","fid", fid  )))
  # temp.record.reference.id <- as.numeric(unlist(get.data.from.table("Record","SourceRefID","fid", fid  )))
  # 
  
  popup = paste0(paste0("<b> Record (fid: ", fid ,")</b>" ),
                 paste0("<hr noshade size='1'/>"),
                 paste0("Species: <i>", species.name,"</i><br>"),
                 paste0("WORMS.ID: ", temp.record.fid,"<br>"),
                 paste0("Status: ", species.status,"<br>"),
                 ifelse( temp.record.site != "", paste0("Site: ", temp.record.site , " (",temp.record.country ,")","<br>"), "" ),
                 ifelse( !is.na(temp.record.year), paste0("Year: ", temp.record.year , "<br>"), "" ),
                 ifelse( !is.na(temp.record.depth), paste0("Depth: ", temp.record.depth , "<br>"), "" ),
                 paste0(temp.record.reference , " ID: ",temp.record.reference.id ,"","<br>")
                 
  )
  
  set.seed(42)
  
  epsg3006 <- leafletCRS(crsClass = "L.Proj.CRS", code = "EPSG:4326",
                         proj4def = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                         resolutions = 2^(13:-1), # 8192 down to 0.5
                         origin = c(0, 0)
  )
  
  if( zoom.define ) { m <- leaflet() %>% setView(lng=mean(data.to.plot[,"Lon"]),lat=mean(data.to.plot[,"Lat"]),zoom=zoom)  }
  if(! zoom.define ) { m <- leaflet()  }
  
  m <- addTiles(m)
  m <- addCircleMarkers(m, 
                        lng=data.to.plot[,"Lon"], 
                        lat=data.to.plot[,"Lat"], 
                        popup= popup , 
                        radius = radius, 
                        color = color , 
                        stroke = FALSE, 
                        fillOpacity = 0.5 )
  
  if(zoom.define) {   setView(m,lng=mean(data.to.plot[,"Lon"]),lat=mean(data.to.plot[,"Lat"]),zoom=zoom )  }
  
  return(m)
  
}

## ------------------------------------------------------------------------------------------------------------------

plot.coords.google <- function(coordinates,info,radius,color) {
  
  if( missing(info)) { info <- "" }
  if( missing(radius)) { radius <- 4 }
  if( missing(color)) { color <- "Black" }
  
  if( max(coordinates[,1], na.rm=T) > 180 | min(coordinates[,1], na.rm=T) < -180 | 
      max(coordinates[,2], na.rm=T) > 90 | min(coordinates[,2], na.rm=T) < -90 ) { errormessage("Coordinates out of WGS84 projection") }
  
  na.records <- which(is.na(coordinates[,1]) )
  
  if( length(na.records) > 0 ) {
    
    cat('\014')
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat(paste0(length(na.records)," records with no coordinate information"))
    cat('\n')
    cat('\n')
    cat('\n') 
    
    coordinates <- coordinates[-na.records,]
    info <- info[-na.records]
    
  }
  
  popup = paste0(paste0("<b> New Record </b>" ),
                 paste0("<hr noshade size='1'/>"),
                 paste0("Info: ",info)
                 
  )
  
  set.seed(42)
  m <- leaflet()
  m <- addTiles(m)
  m <- addCircleMarkers(m, 
                        lng=coordinates[,1], 
                        lat=coordinates[,2], 
                        popup= popup , 
                        radius = radius, 
                        color = color , 
                        stroke = FALSE, 
                        fillOpacity = 0.5 )
  
  return(m)
  
}

## ------------------------------------------------------------------------------------------------------------------

construct.geocode.url <- function(address, return.call = "json", sensor = "false") {
  
  root <- "http://www.mapquestapi.com/geocoding/v1/"
  u <- paste(root, "address?key=mpIx2AWq4Lj9R0mDbW1hNWrPe1Jju4X9&location=", address, sep = "")
  return(URLencode(u))
}

## ------------------------------------------------------------------------------------------------------------------

Geocode.address <- function(address,verbose=FALSE) {
  
  if(verbose) cat(address,"\n")
  
  u <- construct.geocode.url(address)
  doc <- getURL(u)
  x <- fromJSON(doc,simplify = FALSE)
  
  #  if(x$status=="OK") {
  lat <- x$results[[1]]$locations[[1]]$latLng$lat
  lng <- x$results[[1]]$locations[[1]]$latLng$lng
  return(c(lat, lng))
  # } else {
  #  return(c(NA,NA))
  #}
}

## ------------------------------------------------------------------------------------------------------------------

exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}

## ------------------------------------------------------------------------------------------------------------------

testInteger <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if(test == TRUE){ return(TRUE) }
  else { return(FALSE) }
}

## ------------------------------------------------------------------------------------------------------------------

to.valid.numeric <- function(vector.raw,type) { 
  
  t.vect <- numeric(length(vector.raw))
  
  for( i in 1:length(vector.raw)) {
    
    t.val <- NA
    tryCatch( t.val <- as.numeric(as.character(vector.raw[i])) , error=function(e) { error <- TRUE })
    
    if( ! is.na(t.val) & type == "absolute" ) { 
      
      if( ! testInteger(t.val) ) { t.val <- NA }
      
    }
    
    if( is.na(t.val) & type == "presence" ) { 
      
      t.val <- 1
      
    }
    
    if( ! is.na(t.val) & type == "presence" ) { 
      
      if( ! testInteger(t.val) ) { t.val <- NA }
      
    }
    
    t.vect[i] <- t.val
    
  }
  
  return(t.vect)
}

## ------------------------------------------------------------------------------------------------------------------

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

## ------------------------------------------------------------------------------------------------------------------

encoding.correction <- function(string) {
  
  results <- character(0)
  
  for(e in 1:length(string)) {
    
    string.test <- as.character(string[e])
    
    if( ! is.na(string.test)) {
      
      string.test <- enc2native(string.test)
      string.test <- sub( "\n" , " " , string.test)
      
      for( enc in 1:nrow(enconfing)) {
        
        special.character <- paste0("<",tolower(enconfing[enc,3]),">")
        encoded.special.character <- tolower(enconfing[enc,2])
        if ( gregexpr(special.character,string.test)[[1]][1] != -1 ) {
          string.test <- gsub( special.character , encoded.special.character , string.test)
        } }
      
    }
    
    results <- c(results,string.test)
    
  }
  
  results <- gsub("\n"," ",results)
  return(results)
  
}

## ------------------------------------------------------------------------------------------------------------------

consequtive.words.string.match <- function(x,y) { 
  
  x <- encoding.correction(x)
  y <- encoding.correction(y)
  
  x <- gsub('\\(', "", x)
  y <- gsub("\\(", "", y)
  
  x <- gsub('\\)', "", x, fixed=TRUE)
  y <- gsub("\\)", "", y, fixed=TRUE)
  
  x <- tolower(x)
  y <- tolower(y)
  
  char.x <- strsplit(x, " ")[[1]]
  char.y <- strsplit(y, " ")[[1]]
  
  r <- 0
  
  for(i in 1:length(char.x)) {
    length.vect <- i:(length(char.x))
    for(j in 1:length(length.vect)) {
      test.string <- paste( char.x[length.vect[1]:length.vect[j]] , collapse = ' ')
      test.string <- gsub( "\\", "", test.string, fixed=TRUE)
      test.string <- gsub("[","",test.string,fixed=TRUE)
      test.string <- gsub("]","",test.string,fixed=TRUE)
      if( length( grep(test.string,y)) > 0) { r <- max(r,length(strsplit(test.string, " ")[[1]])) }
    } }
  
  return(r)
  
}


## ------------------------------------------------------------------------------------------------------------------

records.ecoregions.climatic <- function(fid, level) {
  
  if( level == 1 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$RLM_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  if( level == 2 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$PROV_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  if( level == 3 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$ECO_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  
  data.points <- get.data.from.table("Record","fid, Lon, Lat","fid",fid)
  data.points <- data.points[!is.na(data.points$Lon),]
  data.points <- data.points[!is.na(data.points$Lat),]
  
  df <- rasters.mahalanobis
  df <- extract(df,data.points[,2:3])
  
  df[ df[,1] > mean(df[,1],na.rm=T),1] <- mean(df[,1],na.rm=T)
  df[ df[,2] < mean(df[,2],na.rm=T),2] <- mean(df[,2],na.rm=T)
  
  m_dist <- tryCatch( m_dist <- mahalanobis(df[, 1:ncol(df)], colMeans(df[, 1:ncol(df)],na.rm=T), cov(df[, 1:ncol(df)],use="pairwise.complete.obs")) , error=function(e) { return("e") })
  
  if( TRUE %in% ( m_dist != "e" ) ) { 
    
    df <- cbind(df,round(m_dist, 2))
    df <- cbind(df, df[,ncol(df)] >= as.vector(quantile(m_dist, probs = c(0.95),na.rm=T)))
    
    # plot(df[,1:2],col=c("Red","Grey")[df[,ncol(df)]+1])
    
    climatic.outliers <- data.points[which(df[,ncol(df)] == 1),1]
    climatic.insiders <- data.points[which(df[,ncol(df)] == 0),1]
    
    coordinates(data.points) <- c("Lon","Lat")
    crs(data.points) <- crs(ecoregions.shp)
    
    results <- cbind(extract(ecoregions.shp,data.points),data.frame(data.points))
    
    na <- results[is.na(results$poly.ID),]$fid
    
    polys.outliers <- unique(results[results$fid %in% climatic.outliers,]$poly.ID)
    polys.insiders <- unique(results[results$fid %in% climatic.insiders,]$poly.ID)
    
    polys.outliers <- polys.outliers[!is.na(polys.outliers)]
    polys.insiders <- polys.insiders[!is.na(polys.insiders)]
    
    library(rgeos)
    library(rgdal)
    
    adjacent.pairs <- gTouches(ecoregions.shp, byid=TRUE)
    
    outliers <- numeric()
    insiders <- numeric()
    
    insiders <- results[ results$poly.ID %in% polys.insiders,]$fid
    
    for( poly in polys.outliers) {
      
      if( (length(results[results$poly.ID %in% poly,]$fid) < 10 ) & ! TRUE %in% ( c(eco.names[poly],names(which(adjacent.pairs[poly,]))) %in% eco.names[polys.insiders]) ) {  
        
        outliers <- c(outliers,results[results$poly.ID %in% poly,]$fid)
        
      }
    }
    
    outliers <- unique(outliers)
    insiders <- unique(insiders)
    
    outliers <- outliers[!outliers %in% na]
    insiders <- insiders[!insiders %in% na]
  } else { outliers <- NULL; insiders <- NULL ;na <- NULL  }
  
  return(list(unsuitable=outliers,ids.suitable=insiders,undetermined=na))
  
}

## ------------------------------------------------------------------------------------------------------------------

records.ecoregions.adjancent <- function(fid, level,manual.threshold) {
  
  if(missing(manual.threshold)) { manual.threshold <- NULL }
  
  if( level == 1 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$RLM_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  if( level == 2 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$PROV_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  if( level == 3 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$ECO_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  
  data.points <- get.data.from.table("Record","fid, Lon, Lat","fid",fid)
  data.points <- data.points[!is.na(data.points$Lon),]
  data.points <- data.points[!is.na(data.points$Lat),]
  
  unique.data.points <- data.points[!duplicated(data.points[,2:3]),]
  coordinates(unique.data.points) <- c("Lon","Lat")
  crs(unique.data.points) <- crs(ecoregions.shp)
  
  results <- cbind(extract(ecoregions.shp,unique.data.points),data.frame(unique.data.points))
  eco.regions.count <- aggregate(rep(1,length(results$poly.ID)), by=list(results$poly.ID), FUN=sum)
  
  if( is.null(manual.threshold) ) {
    
    threhold.value <- 0
    
    repeat {
      
      threhold.value <- threhold.value + 1
      temp <- round(threhold.value / sum(eco.regions.count$x),digits=3) <= 0.01
      if( ! temp ) { break }
      
    }
    
    if( threhold.value == 0 ) { threhold.value <- 1 }
    
    
  }
  
  if( ! is.null(manual.threshold) ) { threhold.value <- manual.threshold }
  
  if( TRUE %in% ( eco.regions.count[,2] > threhold.value) & TRUE %in% (eco.regions.count[,2] <= threhold.value) ) {
    
    library(rgeos)
    library(rgdal)
    
    adjacent.pairs <- gTouches(ecoregions.shp, byid=TRUE)
    polys <- eco.regions.count[which(eco.regions.count[,2] <= threhold.value),1]
    ids.flaged <- numeric(0)
    
    for( poly in polys) {
      
      if( ! TRUE %in% ( eco.names[adjacent.pairs[ poly ,]] %in% eco.names[ eco.regions.count[which(eco.regions.count[,2] > threhold.value),1] ]  )) { 
        
        ids.flaged <- c(ids.flaged,
                        
                        data.points[data.points$Lon == results[which(results$poly.ID == poly),"Lon"] & data.points$Lat == results[which(results$poly.ID == poly),"Lat"]  , "fid"]
                        
        )
        
      }
      
    }
    
  } else { ids.flaged <- NULL }
  
  return(ids.flaged)
  
}

## ------------------------------------------------------------------------------------------------------------------

records.ecoregions <- function(fid, level,manual.threshold) {
  
  if(missing(manual.threshold)) { manual.threshold <- NULL }
  
  if( level == 1 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$RLM_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  if( level == 2 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$PROV_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  if( level == 3 ) {
    
    ecoregions.shp <- unionSpatialPolygons(ecoregions,ecoregions$ECO_CODE)
    eco.names <- names(ecoregions.shp)
    
  }
  
  data.points <- get.data.from.table("Record","fid, Lon, Lat","fid",fid)
  data.points <- data.points[!is.na(data.points$Lon),]
  data.points <- data.points[!is.na(data.points$Lat),]
  
  unique.data.points <- data.points[!duplicated(data.points[,2:3]),]
  coordinates(unique.data.points) <- c("Lon","Lat")
  crs(unique.data.points) <- crs(ecoregions.shp)
  
  results <- cbind(extract(ecoregions.shp,unique.data.points),data.frame(unique.data.points))
  eco.regions.count <- aggregate(rep(1,length(results$poly.ID)), by=list(results$poly.ID), FUN=sum)
  
  if(nrow(eco.regions.count) > 0 ) {
    
    
    
    if( is.null(manual.threshold) ) {
      
      threhold.value <- 0
      
      repeat {
        
        threhold.value <- threhold.value + 1
        temp <- round(threhold.value / sum(eco.regions.count$x),digits=3) <= 0.01
        if( ! temp ) { break }
        
      }
      
      if( threhold.value == 0 ) { threhold.value <- 1 }
      
      
    }
    
    if( ! is.null(manual.threshold) ) { threhold.value <- manual.threshold }
    
    if( TRUE %in% ( eco.regions.count[,2] > threhold.value) & TRUE %in% (eco.regions.count[,2] <= threhold.value) ) {
      
      ids.flaged <- results$fid[which(results$poly.ID %in% eco.regions.count[which(eco.regions.count$x <= threhold.value),1])]
      ids.flaged.suitable <- results$fid[which(results$poly.ID %in% eco.regions.count[which(eco.regions.count$x > threhold.value),1])]
      
    } else { ids.flaged <- NULL }
    
  } else { ids.flaged <- NULL }
  
  return(list(unsuitable=ids.flaged,ids.suitable=flaged.suitable))
  
}

## ------------------------------------------------------------------------------------------------------------------

magic.depths <- function(content) {
  
  content.i <- content
  
  if( grepl("-", content) ) {
    
    content.i <- gsub( " "  , "", content.i )
    content.i <- gsub( ","  , ".", content.i )
    content.i <- gsub( "m"  , "", content.i )
    
    cont.spacer <- unlist(gregexpr("-", content.i))
    cont.numeric <- unlist(gregexpr("[^-]", content.i))
    
    result.i.1 <- as.numeric(substr(content.i, 1, cont.spacer[1]-1)) 
    result.i.2 <- as.numeric(substr(content.i, cont.spacer[1] + 1, nchar(content.i)))
    
    return(c( min(result.i.1,result.i.2) , max(result.i.1,result.i.2)   ))
    
  }
  
  else {
    
    return(c(NA,NA))
    
  }
  
}

## ------------------------------------------------------------------------------------------------------------------

magic.decimals <- function(content) {
  
  if( grepl("[^0-9.-]", content) ) {
    
    is.negative <- grepl("S",content) | grepl("W",content) | grepl("-",content)
    
    content.i <- gsub( ","  , ".", content )
    content.i <- gsub( " "  , "", content.i )
    content.i <- gsub( "E"  , "", content.i )
    content.i <- gsub( "W"  , "", content.i )
    content.i <- gsub( "N"  , "", content.i )
    content.i <- gsub( "S"  , "", content.i )
    content.i <- gsub( "[^0-9.]"  , " ", content.i )
    
    cont.spacer <- unlist(gregexpr(" ", content.i))
    cont.numeric <- unlist(gregexpr("[^ ]", content.i))
    
    if( TRUE %in% (cont.spacer < cont.numeric[1]) ) { 
      
      substr(content.i, cont.numeric[1], nchar(content.i) )
      cont.spacer <- cont.spacer[! cont.spacer < cont.numeric[1]]
      
    }
    
    result.i.1 <- as.numeric(substr(content.i, 1, cont.spacer[1]-1)) 
    result.i.2 <- as.numeric(substr(content.i, cont.spacer[1] + 1, cont.spacer[2]-1)) / 60
    result.i.3 <- as.numeric(substr(content.i, cont.spacer[2] + 1, nchar(content.i))) / 3600
    result.i <- sum(c(result.i.1,result.i.2,result.i.3),na.rm=TRUE)
    
    if(floor(result.i.1) < floor(result.i)) { stop() }
    if( is.na(result.i)) { stop() }
    
    if(is.negative) { result.i <- result.i * (-1) }
    if(result.i == 0) { result.i <- NA }
    if(result.i == "") { result.i <- NA }
    
    return(result.i)
    
  } else {
    
    return(content)
    
  }
  
  
  
}

## ------------------------------------------------------------------------------------------------------------------

records.env.threshold <- function(fid) {
  
  data.points <- get.data.from.table("Record","fid, Lon, Lat","fid",fid)
  coordinates(data.points) <- c("Lon","Lat")
  results <- extract(rasters,data.points)
  
  return(results)
  
}

## ------------------------------------------------------------------------------------------------------------------

records.over.land <- function(fid,process.world,process.coastline) {
  
  options(warn=-1)
  
  data.points <- get.data.from.table("Record","fid, Lon, Lat","fid",fid)
  data.points <- data.points[!is.na(data.points$Lon) & !is.na(data.points$Lon) ,]
  
  data.points.geo <- data.points
  coordinates(data.points.geo) <- c("Lon","Lat")
  crs(data.points.geo) <- crs(world)
  
  if( process.world ) { 
    
    land.t <- crop(world,extent(data.points.geo) + c(-1,+1,-1,+1)) } else { land.t <- world }
  
  if( process.coastline ) { 
    
    tryCatch( coast.line.t <- crop(coast.line,extent(data.points.geo) + c(-1,+1,-1,+1)) , error=function(e) { } , finally = { coast.line.t <- coast.line } )
    coast.line.t <- as.data.frame(as(coast.line, "SpatialPointsDataFrame"))
    
  }
  
  if( ! is.null(land.t) ) { 
    
    over.land <- sp::over(data.points.geo, land.t)
    over.land <- as.vector(sapply(over.land[,1],function(x) { ifelse(is.na(x),1,-1)}))
    
    data.points.dist <- as.matrix(data.points,ncol=2)[,c("Lon","Lat")]
    coast.line.dist <- as.matrix(data.frame(Lon=coast.line.t[,5],Lat=coast.line.t[,6]))
    
    distances <- data.frame()
    
    for (d in 1:nrow(matrix(data.points.dist,ncol=2))) {
      
      if( over.land[d] != 1 ) {
        
        distance <- spDists(matrix(matrix(data.points.dist,ncol=2)[d,],ncol=2),coast.line.dist , longlat=TRUE)
        distances <- c(distances,apply(distance,1,min))
        
      } else {
        
        distances <- c(distances,0)
        
      }
      
      
    }
    
    results <- data.frame(fid = data.points$fid,
                          OnLand = over.land,
                          Distance = as.vector(unlist(distances)),
                          stringsAsFactors = FALSE)
    
    rm(coast.line.t)
    
  } else {
    
    results <- data.frame(fid = data.points$fid,
                          OnLand = rep(1,length(data.points$fid)),
                          Distance = rep(0,length(data.points$fid)),
                          stringsAsFactors = FALSE) 
    
  }
  
  return(results)
  
}

## ------------------------------------------------------------------------------------------------------------------

reset.database <- function() {
  
  cat('\014')
  cat('\n')
  cat('\n')
  cat('----------------------------------------------------------------------')
  cat('\n')
  cat("Overwriting current information with new clean DataBase.")
  cat('\n')
  cat('\n')
  
  x <- ""
  while( x != "Y" & x != "n" ) { 
    x <- readline("Are you Sure (Y/n): ")
  }
  
  cat('\n')
  
  if ( x == "Y" ) {  
    
    file.to.remove <- list.files(paste0(main.directory,"/Data"),pattern="sql",full.names = TRUE)
    if(length(file.to.remove) > 0) { file.remove(file.to.remove) }
    
    # First record
    
    new.record <- data.frame(     fid = 1 ,
                                  SpeciesID = as.numeric("") ,
                                  Lon = as.numeric("") ,
                                  Lat = as.numeric("") ,
                                  CoordinateType = "",
                                  CoordinateUncertaintyInMeters = as.numeric(""),
                                  Country = "" ,
                                  Site = "" ,
                                  Depth = as.numeric("") ,
                                  DepthRangeMin = as.numeric("") ,
                                  DepthRangeMax = as.numeric("") ,
                                  DateYear = as.numeric("") ,
                                  DateMonth = as.numeric("") ,
                                  DateDay = as.numeric("") ,
                                  AbundanceID = as.numeric("") ,            
                                  SampleCollectionID = as.numeric("") ,
                                  SourceRefType = "" ,       
                                  SourceRefID = as.numeric("") ,
                                  RecordUserID = as.numeric("") ,      
                                  RecordDateYear = as.numeric("") ,             
                                  RecordDateMonth = as.numeric("") ,            
                                  RecordDateDay = as.numeric("") ,
                                  RecordPublic = as.numeric("") ,
                                  RecordNotes="" , 
                                  OriginalSource = "" ,
                                  Flag = "0",
                                  OriginalSourceID = "" ,
                                  FlagOnLand = "0",
                                  FlagLightSuitable = "0",
                                  FlagO2Suitable = "0",
                                  FlagEcoRegion = "0",
                                  stringsAsFactors = FALSE )
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbWriteTable(sql, "MDB_GeoData", new.record , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
    ## ------------------------------
    # First Abundance
    
    new.abundance <- data.frame(  fid = 1 ,
                                  AbundanceVal = as.numeric("") ,
                                  AbundanceUnit = "" ,
                                  AbundanceStatistic = "" ,
                                  AbundanceMethod = "",
                                  AbundanceReplicates = as.numeric("") ,
                                  AbundanceNotes = "", stringsAsFactors = FALSE)
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbWriteTable(sql, "RDB_Abundance", new.abundance , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
    ## ------------------------------
    # First SampleCollection
    
    new.sample <- data.frame(     fid = 1 ,
                                  SampleType = "" , # HerbariumSpecimen Picture
                                  SampleDesign = "" ,
                                  SampleN = as.numeric("") ,
                                  SampleCollector = "" ,
                                  SampleStorageLabel = "" ,
                                  SampleStorageLocation = "" ,
                                  SampleConservationMethod = "" ,
                                  GeneticExtractionType = "" ,
                                  GeneticExtractionLocation = "" ,
                                  SampleNotes = "" , stringsAsFactors = FALSE)
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbWriteTable(sql, "RDB_Sample_Collection", new.sample , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
    ## ------------------------------
    # First Species
    
    new.species <- data.frame(    fid = 1 ,
                                  SpeciesWormsID = as.numeric("") ,
                                  SpeciesName = "" ,
                                  Authority = "",
                                  Status = "",
                                  TaxKingdom = "",
                                  TaxPhylum = "",
                                  TaxClass = "",
                                  TaxOrder = "",
                                  TaxFamily = "",
                                  TaxGenus = "",
                                  ReviewedByWorms = as.numeric("") ,
                                  RevisionByWormsDate = "",
                                  AcceptedWormsID = as.numeric("") ,
                                  AcceptedSpeciesName = "" , 
                                  Description = "" , 
                                  DescriptionSource = "" , 
                                  DescriptionSourceLink = "" , 
                                  Distribution = "" , 
                                  DistributionSource = "" , 
                                  CommonName = "" , stringsAsFactors = FALSE)
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbWriteTable(sql, "RDB_Taxa", new.species , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
    ## ------------
    # First reference
    
    new.reference <- data.frame(    fid = 1 ,
                                    Type = "" ,
                                    Authors = "" ,
                                    Title = "",
                                    Year = "" ,
                                    Container = "" ,
                                    Volume = "" ,
                                    Issue = "" ,
                                    Page = "" ,
                                    FullCitation="",
                                    DOI = "" ,
                                    ReviewedByCrossref = as.numeric("") ,
                                    RevisionByCrossrefDate="" , stringsAsFactors = FALSE)
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbWriteTable(sql, "RDB_Ref_Literature", new.reference , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
    ## ------------
    # First Repository
    
    new.database <- data.frame(     fid = 1 ,
                                    RepositoryType = "" ,
                                    RepositoryName = "" ,
                                    RepositoryShortName = "" ,
                                    Institution = "" ,
                                    Responsible = "",
                                    ResponsibleContact = "",
                                    DOI = "" ,
                                    URL = "" , stringsAsFactors = FALSE )
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbWriteTable(sql, "RDB_Ref_Repository", new.database , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
    ## --------------------------------
    # First User
    
    new.user <- data.frame( fid = 1 ,
                            Name = "Jorge Assis" ,
                            Login = "jorgeassis" ,
                            Password = "cruzeiro" ,
                            LastLogDate = "" ,
                            DBLevel = 1,
                            Institution = "Centre for Marine Sciences" ,
                            URL = "http://.ccmar.ualg.pt/" , stringsAsFactors = FALSE )
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbWriteTable(sql, "RDB_Users", new.user , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
    ## --------------------------------
    # First version
    
    version.control <- data.frame(  Created = as.character(Sys.time()) ,
                                    Created.by.user = UserIDSession , 
                                    Updated = as.character(Sys.time()) , 
                                    Last.update.by.user = UserIDSession , stringsAsFactors = FALSE )
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(main.directory,sql.path))
    dbWriteTable(sql, "Version", version.control , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
    ## --------------------------------
    
    delete.records("Abundance",1,TRUE)
    delete.records("Collection",1,TRUE)
    delete.records("Taxa",1,TRUE)
    delete.records("Literature",1,TRUE)
    delete.records("Repository",1,TRUE)
    delete.records("Record",1,TRUE)
    
  }
}

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------