

maxent_map <- function(sciName, minXlim, maxXlim, minYlim, maxYlim, modelStudyAreaOnly, mapStudyAreaOnly, futureModel) {
  

  ####### RUN FROM HERE #######
    
      # Load packages
     {
      Sys.setenv(NOAWT=TRUE)
      library(rgbif)
      suppressMessages(library(rgdal))
      suppressMessages(library(sf))
      suppressMessages(library(dplyr))
      library(dismo)
      library(rJava)
      library(rinat)
      library(ridigbio)
      library(CoordinateCleaner)
      suppressMessages(library(spThin))
      }
    # Download GBIF, iNaturalist, and iDigBio records for species
    
    
      print("Downloading GBIF records...")
      try(
      gbifRecords <- occ_search(scientificName = sciName, hasGeospatialIssue = FALSE, hasCoordinate = TRUE, limit = 5000)
      )
        if (typeof(gbifRecords[1]) == "character") {
          print("ERROR: No GBIF records found. Check correct spelling and taxonomic synonymy at GBIF.org")
        }
        else{
        gbifRecords <- dplyr::select(gbifRecords$data, decimalLongitude, decimalLatitude)
        
        print(nrow(gbifRecords))
        }
      
      
      print("Downloading iNaturalist records...")
      try(
        inatRecords <- get_inat_obs(query = sciName, geo = TRUE, maxresults = 5000, quality = "research")
      )
        if(exists("inatRecords")){
          inatRecords <- inatRecords %>% dplyr::select(longitude, latitude)
          names(inatRecords)[1] <- "decimalLongitude"
          names(inatRecords)[2] <- "decimalLatitude"
          print(nrow(inatRecords))
        }
      
      print("Downloading iDigBio records...")
      try(
      idigbioRecords <- idig_search_records(rq=list(scientificname = sciName, geopoint=list(type="exists")), fields=c("geopoint"), limit=5000)
      )
        if (nrow(idigbioRecords) == 0) {
          print("ERROR: No iDigBio records found. Check correct spelling and taxonomic synonymy at iDigBio.org")
        }
        else {
          names(idigbioRecords)[1] <- "decimalLongitude"
          names(idigbioRecords)[2] <- "decimalLatitude"
          print(nrow(idigbioRecords))
        }
      
      
      # Check for which coords exist and set to coords variable
      if (exists("gbifRecords") & exists("inatRecords") & exists("idigbioRecords")){
        coords <- rbind(gbifRecords, inatRecords, idigbioRecords)
      }
       else if (exists("gbifRecords") & exists("inatRecords")){
          coords <- rbind(gbifRecords, inatRecords)
       }
        else if (exists("inatRecords") & exists("idigbioRecords")){
            coords <- rbind(inatRecords, idigbioRecords)
        }
          else if (exists("gbifRecords") & exists("idigbioRecords")){
            coords <- rbind(gbifRecords, idigbioRecords)
          }
            else if (exists("gbifRecords")){
              coords <- gbifRecords
            }
              else if (exists("inatRecords")){
                coords <- inatRecords
              }
                else if (exists("idigbioRecords")){
                  coords <- idigbioRecords
                }
        
        
      
        
      # Check if any coords were found, if not end script
      if (exists("coords")){
        
      # Clean coordinates with 'Coordinate Cleaner' package
      print ("Cleaning coordinates...")
      coords <- cbind(sciName, coords)
      flags <- clean_coordinates(coords, species = "sciName", lon = "decimalLongitude", lat = "decimalLatitude", 
                    tests = c("capitals", "centroids", "duplicates", "equal", "gbif", "institutions", "seas", "zeros"))
      
      # Exclude problematic records
      coords <- coords[flags$.summary,]
      print("Number of coordinates passing tests:")
      print(nrow(coords))
      coords <- coords[,2:3]
      thin.coords <- thin.algorithm(coords, 8, 10)
      coords <- SpatialPoints(thin.coords[[10]], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      print("Number of coordinates after thinning:")
      print(nrow(coords@coords))
        
      # Get min/max values for latitude and longitude of sample points
      minLong <- min(coords@coords[,1])
      maxLong <- max(coords@coords[,1])
      minLat <- min(coords@coords[,2])
      maxLat <- max(coords@coords[,2])
      
      
      # Check if sample coordinates cover Carolina region, or if boundaries need to be manually enforced
      if (minXlim == 0 & maxXlim == 0 & minYlim == 0 & maxYlim == 0){
        print("Modeling entire area with coords...")
      
        # Adjust spatial extent to make square
        diffLong <- maxLong-minLong
        diffLat <- maxLat-minLat
        if (diffLong > diffLat){
          diffDim <- diffLong-diffLat
          spatialExtent <- c(minLong-2, minLat-(0.5*diffDim)-1, maxLong+2, maxLat+(0.5*diffDim)+1) 
        } else {
          diffDim <- diffLat-diffLong
          spatialExtent <- c(minLong-(0.5*diffDim)-2, minLat-1, maxLong+(0.5*diffDim)+2, maxLat+1) 
          if (spatialExtent[1] < -180){
            spatialExtent <- -180
            }
          if (spatialExtent[3] > 180){
            spatialExtent[3] <- 180
            }
          if (spatialExtent[2] < -90){
            spatialExtent[2] <- -90
            }
          if (spatialExtent[4] > 90){
            spatialExtent[4] <- 90
            }
          }
      }
      else if (modelStudyAreaOnly == "yes"){
        print("Modeling only study area...")
        spatialExtent <- c(minXlim-2,minYlim-1,maxXlim+2,maxYlim+1)
      } else {
        print("Modeling combination of area with coords and study area...")
        spatialExtent <- c(min(minXlim,minLong)-2,min(minYlim,minLat)-1,max(maxXlim,maxLong)+2,max(maxYlim,maxLat)+1)
      }
      spatialExtent <- matrix(spatialExtent, nrow = 2, byrow = TRUE)
      spatialExtent <- SpatialPoints(spatialExtent)
      
      print("Cropping predictor layers...")
      predictorsCropped <- crop(bioclimPredictors,spatialExtent)
      if (futureModel == "yes"){
      futurePredictorsCropped <- crop(futureBioclimPredictors, spatialExtent)
      }
    
    # Run MAXENT model
      print("Running MaxEnt model...this may take a long time!")
      me <- maxent(predictorsCropped, coords, removeDuplicates = TRUE)
      if (futureModel == "no"){
      mePredict <- predict(me, predictorsCropped, progress = 'text')
      } else if (futureModel == "yes"){
        futureMePredict <- predict(me, futurePredictorsCropped, progress = 'text')
      }
    
    # Plot MAXENT results
      print("Plotting results...")
      if (minXlim == 0 & maxXlim == 0 & minYlim == 0 & maxYlim == 0){
        print("Mapping entire area with coords...")
        mapXlim <- c(minLong-2, maxLong+2)
        mapYlim <- c(minLat-1, maxLat+1)
      }
      else if (mapStudyAreaOnly == "yes"){
        print("Mapping only study area...")
        mapXlim <- c(minXlim, maxXlim)
        mapYlim <- c(minYlim, maxYlim)
      } else {
        print("Mapping combination of area with coords and study area...")
        mapXlim <- c(min(minXlim,minLong-2), max(maxXlim,maxLong+2))
        mapYlim <- c(min(minYlim,minLat-1), max(maxYlim,maxLat+1))
      }
      par(mar = c(2,2,2,5))
      cexAdj <- (10/(mapXlim[2]-mapXlim[1]))
      if (futureModel == "yes"){
        plot(futureMePredict, axes = TRUE, main = paste(sciName, "Predicted", sep = " "), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1), col = rev(grDevices::terrain.colors(20)),xlim = mapXlim, ylim = mapYlim )# Uncomment to limit map extent to Carolinas and adjacent states
      } else {
        plot(mePredict, axes = TRUE, main = sciName, breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1), col = rev(grDevices::terrain.colors(20)),xlim = mapXlim, ylim = mapYlim )# Uncomment to limit map extent to Carolinas and adjacent states
        
      }
      plot(counties.simple, border = "black", lwd = 0.05, add = TRUE)
      plot(us.states.simple, add = TRUE)
      points(coords, pch = 21, col = "black", bg = "red", cex = cexAdj )

      
   }####### END RUN #######
      
  else {print("ERROR: No GBIF or iNaturalist records returned. Check spelling and synonymy of name")}

}
  