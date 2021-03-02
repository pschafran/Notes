setwd("/Volumes/Samsung_T5/GIS/Isoetes")
{
  Sys.setenv(NOAWT=TRUE)
  library(rgbif)
  library(maptools)
  library(rgeos)
  library(raster)
  library(rgdal)
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(dismo)
  library(rJava)
  library(rinat)
  library(maps)
}

# Prepare state boundaries
{
  states <- c("North Carolina", "South Carolina", "Georgia", "Virginia", "Tennessee", "Alabama", "West Virginia", "Kentucky")
  us.states <- usaHD[usaHD$NAME_1 %in% states,]
  us.states.simple <- gSimplify(us.states, tol = 0.01, topologyPreserve = TRUE)
}
# Prepare predictor rasters
{
  path <- file.path("/Volumes/Samsung_T5/GIS/CM10_1975H_Bio_ASCII_V1.2/CM10_1975H_Bio_V1.2")
  climondFiles <- list.files(path, pattern = 'txt$', full.names = TRUE)
  climondPredictors <- stack(climondFiles)
  spatialExtent <- c(-85,31,-75,38)
  spatialExtent <- matrix(spatialExtent, nrow = 2, byrow = TRUE)
  spatialExtent <- SpatialPoints(spatialExtent)
  predictorsCropped <- crop(climondPredictors,spatialExtent)
}

# Prepare Carolina and bordering cities
data("world.cities")
usa.cities <- world.cities[ which(world.cities$country.etc == "USA"),]
carolina.cities <- usa.cities[which(usa.cities$name == "Charlotte" | usa.cities$name == "Fayetteville" | usa.cities$name == "Asheville" |usa.cities$name == "Winston-Salem" | usa.cities$name == "Raleigh" | usa.cities$name == "Wilmington" | usa.cities$name == "Johnson City" | usa.cities$name == "Columbia" | usa.cities$name == "Charleston" | usa.cities$name == "Myrtle Beach" | usa.cities$name == "Greenville" | usa.cities$name == "Norfolk" | usa.cities$name == "Richmond" | usa.cities$name == "Danville" | usa.cities$name == "Savannah" | usa.cities$name == "Atlanta" | usa.cities$name == "Macon" | usa.cities$name == "Augusta-Richmond" | usa.cities$name == "Knoxville" | usa.cities$name == "Roanoke"),]
carolina.cities.labels <- carolina.cities$name
carolina.cities.coords <- SpatialPoints(carolina.cities[5:4], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

{####### RUN FROM HERE #######
  
# Input species name  
  sciName = "Asimina triloba"
  
# Download GBIF and iNaturalist records for species

{
  gbifRecords <- occ_search(scientificName = sciName, hasGeospatialIssue = FALSE, hasCoordinate = TRUE, return = "data", limit = 5000)
  gbifCoords <- gbifRecords %>% dplyr::select(decimalLongitude, decimalLatitude)  
  gbifCoords <- SpatialPoints(gbifCoords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  inatRecords <- get_inat_obs(query = sciName, geo = TRUE, maxresults = 5000, quality = "research")
  inatCoords <- inatRecords %>% dplyr::select(longitude, latitude)  
  inatCoords <- SpatialPoints(inatCoords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  coords <- rbind(gbifCoords, inatCoords)
}

# Run MAXENT model
{
  me <- maxent(predictorsCropped, coords, removeDuplicates = TRUE)
  mePredict <- predict(me, predictorsCropped)
}

# Plot MAXENT results
{
  par(mfrow = c(1,1))
  par(mar = c(2,2,2,5))
  plot(mePredict, axes = TRUE, main = sciName)
  plot(us.states.simple, add = TRUE)
  points(coords, pch = 21, col = "black", bg = "red")
  points(carolina.cities.coords, pch = 20, cex = 1)
  text(carolina.cities.coords, labels = carolina.cities.labels, halo = TRUE, hw = 0.2, hc = "white", pos = 1, cex = 0.75)
}

} ####### END RUN #######
