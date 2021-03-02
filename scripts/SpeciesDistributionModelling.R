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
}
usaHD <- getData(name="GADM", download = TRUE, level = 1, country = "USA")

#### Carolinas

states <- c("North Carolina", "South Carolina")
us.states <- usaHD[usaHD$NAME_1 %in% states,]
us.states.simple <- gSimplify(us.states, tol = 0.01, topologyPreserve = TRUE)

records <- occ_search(scientificName = "Vaccinium corymbosum", country = "US", hasGeospatialIssue = FALSE, hasCoordinate = TRUE, return = "data", limit = 1000)
coords <- records %>% dplyr::select(decimalLongitude, decimalLatitude)  
coords <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

path <- file.path(system.file(package="dismo"), 'ex')
files <- list.files(path, pattern='grd$', full.names=TRUE )
predictors <- stack(files)

plot(us.states.simple, add = TRUE)
points(coords[us.states.simple,], pch = 21, col = "black", bg = "red")
plot(predictors$bio1)
{
  getData(name = "worldclim", download = TRUE, var = "tmin", res = 0.5, lon = 35, lat = -76)
  getData(name = "worldclim", download = TRUE, var = "tmax", res = 0.5, lon = 35, lat = -76)
  getData(name = "worldclim", download = TRUE, var = "prec", res = 0.5, lon = 35, lat = -76)
  getData(name = "worldclim", download = TRUE, var = "bio", res = 0.5, lon = 35, lat = -76)
}
bc <- bioclim(predictors, coords)
bcPredict <- predict(bc, predictors)

dm <- domain(predictors, coords)
dmPredict <- predict(dm, predictors)

m <- mahal(predictors, coords)
mPredict <- predict(m, predictors)

me <- maxent(predictors, coords, factors='biome')
mePredict <- predict(me, predictors)

plot(us.states.simple)
plot(bcPredict, add = TRUE)
plot(us.states.simple, add = TRUE)
points(coords[us.states.simple,], pch = 21, col = "black", bg = "red")
  
plot(us.states.simple)
plot(dmPredict, add = TRUE)
plot(us.states.simple, add = TRUE)
points(coords[us.states.simple,], pch = 21, col = "black", bg = "red")

plot(us.states.simple)
plot(mPredict, add = TRUE)
plot(us.states.simple, add = TRUE)
points(coords[us.states.simple,], pch = 21, col = "black", bg = "red")

{
  sciName <- "Isoetes melanopoda"
  
  states <- c("North Carolina", "South Carolina")
  us.states <- usaHD[usaHD$NAME_1 %in% states,]
  us.states.simple <- gSimplify(us.states, tol = 0.01, topologyPreserve = TRUE)
  
  records <- occ_search(scientificName = sciName, country = "US", hasGeospatialIssue = FALSE, hasCoordinate = TRUE, return = "data", limit = 5000)
  coords <- records %>% dplyr::select(decimalLongitude, decimalLatitude)  
  coords <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  me <- maxent(climondPredictors, coords, removeDuplicates = TRUE)
  mePredict <- predict(me, climondPredictors)
  par(mfrow = c(1,1))
  par(mar = c(2,2,2,5))
  plot(us.states.simple, axes = TRUE, main = sciName)
  plot(mePredict, add = TRUE)
  plot(us.states.simple, add = TRUE)
  points(coords[us.states.simple,], pch = 21, col = "black", bg = "red")
}




#### Eastern US Isoetes collections


path <- file.path("/Volumes/Samsung_T5/GIS/CM10_1975H_Bio_ASCII_V1.2/CM10_1975H_Bio_V1.2")
climondFiles <- list.files(path, pattern = 'txt$', full.names = TRUE)
climondPredictors <- stack(climondFiles)
View(climondPredictors)

usaHD <- getData(name="GADM", download = TRUE, level = 1, country = "USA")
canadaHD <- getData(name="GADM", download = TRUE, level = 1, country = "CAN")
greatLakes <- st_read("./Great_Lakes/Great_Lakes.shp")

states <- c("Texas", "Oklahoma", "Nebraska", "Kansas", "North Dakota", "South Dakota", "Minnesota", "Iowa", "Missouri", "Arkansas", "Louisiana", "Wisconsin", "Illinois", "Kentucky", "Tennessee", "Mississippi", "Michigan", "Indiana", "Alabama", "Ohio", "Georgia", "Florida", "New York", "Pennsylvania", "West Virginia" ,"Maryland", "Virginia", "North Carolina", "South Carolina", "Delaware", "New Jersey", "Connecticut", "Rhode Island","Massachusetts", "New Hampshire", "Vermont","Maine")
provinces <- c("Ontario", "Quebec", "New Brunswick", "Nova Scotia")

us.states <- usaHD[usaHD$NAME_1 %in% states,]
us.states.simple <- gSimplify(us.states, tol = 0.01, topologyPreserve = TRUE)
can.provinces <- canadaHD[canadaHD$NAME_1 %in% provinces,]
can.provinces.simple <- gSimplify(can.provinces, tol = 0.01, topologyPreserve = TRUE)

allIsoetes <- read.csv("All_Isoetes.csv", header = TRUE, sep = "\t")
coords <- allIsoetes %>% dplyr::select(Longitude, Latitude)  
coords <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

piedNCpts <- read.csv("sch86_141-144_GPS.csv", header = FALSE, sep = ",")
coords <- piedNCpts %>% dplyr::select(V2, V1)  
coords <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
View(coords)


me <- maxent(climondPredictors, coords, removeDuplicates = TRUE)
mePredict <- predict(me, climondPredictors)
par(mfrow = c(1,1))
par(mar = c(2,2,2,5))
plot(us.states.simple)
plot(mePredict, add = TRUE)
plot(us.states.simple, add = TRUE)
points(coords[us.states.simple,], pch = 21, col = "black", bg = "red")
points(piedNCpts$V2, piedNCpts$V1)



par(mfrow = c(2,2))
### Run bioclim model on coords
{ #run all models
{
  bc <- bioclim(climondPredictors, coords)
  bcPredict <- predict(bc, climondPredictors)
  par(mar = c(5,5,2,2))
  plot(us.states.simple, axes = TRUE, xlim = c(-95,-65), ylim = c(27,47))
  plot(can.provinces.simple, add = TRUE)
  plot(bcPredict, add = TRUE)
  plot(us.states.simple, add = TRUE)
  plot(can.provinces.simple, add = TRUE)
  points(coords, pch = 21, col = "black", bg = "red")
}

### Run DOMAIN model on coords
{
  dm <- domain(climondPredictors, coords)
  dmPredict <- predict(dm, climondPredictors) 
  par(mar = c(5,5,2,2))
  plot(us.states.simple, axes = TRUE, xlim = c(-95,-65), ylim = c(27,47))
  plot(can.provinces.simple, add = TRUE)
  plot(dmPredict, add = TRUE)
  plot(us.states.simple, add = TRUE)
  plot(can.provinces.simple, add = TRUE)
  points(coords, pch = 21, col = "black", bg = "red")
}

### Run MAXENT model on coords
{
  me <- maxent(climondPredictors, coords, removeDuplicates = TRUE)
  mePredict <- predict(me, climondPredictors) 
  par(mar = c(5,5,2,2))
  plot(us.states.simple, axes = TRUE, xlim = c(-95,-65), ylim = c(27,47))
  plot(can.provinces.simple, add = TRUE)
  plot(mePredict, add = TRUE)
  plot(us.states.simple, add = TRUE)
  plot(can.provinces.simple, add = TRUE)
  points(coords, pch = 21, col = "black", bg = "red")
}
}


##### Isoetes GBIF

sciName <- "Isoetes riparia"

{ 
  states <- c("Texas", "Oklahoma", "Nebraska", "Kansas", "North Dakota", "South Dakota", "Minnesota", "Iowa", "Missouri", "Arkansas", "Louisiana", "Wisconsin", "Illinois", "Kentucky", "Tennessee", "Mississippi", "Michigan", "Indiana", "Alabama", "Ohio", "Georgia", "Florida", "New York", "Pennsylvania", "West Virginia" ,"Maryland", "Virginia", "North Carolina", "South Carolina", "Delaware", "New Jersey", "Connecticut", "Rhode Island","Massachusetts", "New Hampshire", "Vermont","Maine")
  provinces <- c("Ontario", "Quebec", "New Brunswick", "Nova Scotia")
  
  us.states <- usaHD[usaHD$NAME_1 %in% states,]
  us.states.simple <- gSimplify(us.states, tol = 0.01, topologyPreserve = TRUE)
  can.provinces <- canadaHD[canadaHD$NAME_1 %in% provinces,]
  can.provinces.simple <- gSimplify(can.provinces, tol = 0.01, topologyPreserve = TRUE)
  
  gbifRecords <- occ_search(scientificName = sciName, hasGeospatialIssue = FALSE, hasCoordinate = TRUE, return = "data", limit = 5000)
  coords <- records %>% dplyr::select(decimalLongitude, decimalLatitude)  
  coords <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  me <- maxent(climondPredictors, coords, removeDuplicates = TRUE)
  mePredict <- predict(me, climondPredictors)
  par(mfrow = c(1,1))
  par(mar = c(2,2,2,5))
  plot(us.states.simple, axes = TRUE, xlim = c(-100,-65), ylim = c(27,47), main = sciName)
  plot(mePredict, add = TRUE)
  plot(us.states.simple, add = TRUE)
  plot(can.provinces.simple, add = TRUE)
  plot(greatLakes, add = TRUE, col = "white")
  points(coords, pch = 21, col = "black", bg = "red")
}

### Isoetes GBIF and iNaturalist combined

sciName <- "Typha latifolia"

{ 
  states <- c("Texas", "Oklahoma", "Nebraska", "Kansas", "North Dakota", "South Dakota", "Minnesota", "Iowa", "Missouri", "Arkansas", "Louisiana", "Wisconsin", "Illinois", "Kentucky", "Tennessee", "Mississippi", "Michigan", "Indiana", "Alabama", "Ohio", "Georgia", "Florida", "New York", "Pennsylvania", "West Virginia" ,"Maryland", "Virginia", "North Carolina", "South Carolina", "Delaware", "New Jersey", "Connecticut", "Rhode Island","Massachusetts", "New Hampshire", "Vermont","Maine")
  provinces <- c("Ontario", "Québec", "New Brunswick", "Nova Scotia", "Newfoundland", "Prince Edward Island")
  
  us.states <- usaHD[usaHD$NAME_1 %in% states,]
  us.states.simple <- gSimplify(us.states, tol = 0.01, topologyPreserve = TRUE)
  can.provinces <- canadaHD[canadaHD$NAME_1 %in% provinces,]
  can.provinces.simple <- gSimplify(can.provinces, tol = 0.01, topologyPreserve = TRUE)
  
  gbifRecords <- occ_search(scientificName = sciName, hasGeospatialIssue = FALSE, hasCoordinate = TRUE, return = "data", limit = 5000)
  gbifCoords <- gbifRecords %>% dplyr::select(decimalLongitude, decimalLatitude)  
  gbifCoords <- SpatialPoints(gbifCoords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  inatRecords <- get_inat_obs(query = sciName, geo = TRUE, maxresults = 5000, quality = "research")
  inatCoords <- inatRecords %>% dplyr::select(longitude, latitude)  
  inatCoords <- SpatialPoints(inatCoords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  coords <- rbind(gbifCoords, inatCoords)
  
  predictorsCropped <- crop(climondPredictors, coords)
  
  me <- maxent(predictorsCropped, coords, removeDuplicates = TRUE)
  mePredict <- predict(me, predictorsCropped)
  par(mfrow = c(1,1))
  par(mar = c(2,2,2,5))
 # plot(us.states.simple, axes = TRUE, xlim = c(-100,-65), ylim = c(27,47), main = sciName)
  plot(mePredict, axes = TRUE, main = sciName)
  plot(us.states.simple, add = TRUE)
  plot(can.provinces.simple, add = TRUE)
  plot(greatLakes, add = TRUE, col = "white")
  points(coords, pch = 21, col = "black", bg = "red")
}
