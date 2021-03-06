#!/usr/bin/env Rscript

# Set study area limits and raster resolution
minXlim <- -85 # Longitude in decimal degrees
maxXlim <- -75 # Longitude in decimal degrees
minYlim <- 31 # Latitude in decimal degrees
maxYlim <- 38 # Latitude in decimal degrees
rasterResolution <- "2.5m" # Must be 10m, 5m, 2.5m, or 30s
path30s <- "GIS/WorldClim/wc2.0_30s_bio/"
path2.5m <- "GIS/WorldClim/wc2.0_2.5m_bio/"
path5m <- "GIS/WorldClim/wc2.0_5m_bio/"
path10m <- "GIS/WorldClim/wc2.0_10m_bio/"
setwd("~/Documents/Website/MaxEntMachine")
modelStudyAreaOnly <- "no" # Must be yes/no. Determines whether to use predictor layers and coordinates outside the study area to model the species distribution
                       # or only the area inside the study area
mapStudyAreaOnly <- "yes" # Must be yes/no. Determines whether to map only study area defined above or entire area containing sample coordinates
futureModel = "no" # Select whether to model future climate model or present day

#Shouldn't need to change anything below this line
################################################################################################################################
{
  Sys.setenv(NOAWT=TRUE)
  suppressMessages(library(rgeos))
  suppressMessages(library(raster))
  library(sp)
  library(maps)
  source("maxentHelper.R")
  library(mapdata)
}

{ ##### Run prep
  
  # Parse scientific names passed from terminal
  sciName = commandArgs(trailingOnly=TRUE)
  if (length(sciName) == 2) {
    sciName = paste(sciName[1], sciName[2], sep = " ")
  } else if (length(sciName) == 3){
    sciName = paste(sciName[1], sciName[2], sciName[3], sep = " ")
  } else if (length(sciName) == 4){
    sciName = paste(sciName[1], sciName[2], sciName[3], sciName[4], sep = " ")
  }
  print("**********************START***********************")
  print(sciName)
  
  # Prepare state & county boundaries
  usaHD <- getData(name="GADM", download = TRUE, level = 1, country = "USA")
  countiesHD <- getData(name="GADM", download = TRUE, country="USA", level=2)
  states <- c("North Carolina", "South Carolina", "Georgia", "Virginia", "Tennessee", "Alabama", "West Virginia", "Kentucky")
  us.states <- usaHD[usaHD$NAME_1 %in% states,]
  us.states.simple <- gSimplify(us.states, tol = 0.01, topologyPreserve = TRUE)
  counties <- countiesHD[countiesHD$NAME_1 %in% states,]
  counties.simple <- gSimplify(counties, tol = 0.01, topologyPreserve = TRUE)
  
  # Prepare predictor rasters
  if (rasterResolution == "30s"){
    path <- file.path(path30s)
  } 
  else if (rasterResolution == "2.5m") {
    path <- file.path(path2.5m)
  }
  else if (rasterResolution == "5m") {
    path <- file.path(path5m)
  }
  else {
    path <- file.path(path10m)
  }
  
  bioclimFiles <- list.files(path, pattern = 'tif$', full.names = TRUE)
  bioclimPredictors <- stack(bioclimFiles)
  
  
  # Prepare Carolina and bordering cities
  data("world.cities")
  usa.cities <- world.cities[ which(world.cities$country.etc == "USA"),]
  carolina.cities <- usa.cities[which(usa.cities$name == "Charlotte" | usa.cities$name == "Fayetteville" | usa.cities$name == "Asheville" |usa.cities$name == "Winston-Salem" | usa.cities$name == "Raleigh" | usa.cities$name == "Wilmington" | usa.cities$name == "Johnson City" | usa.cities$name == "Columbia" | usa.cities$name == "Charleston" | usa.cities$name == "Myrtle Beach" | usa.cities$name == "Greenville" | usa.cities$name == "Norfolk" | usa.cities$name == "Richmond" | usa.cities$name == "Danville" | usa.cities$name == "Savannah" | usa.cities$name == "Atlanta" | usa.cities$name == "Macon" | usa.cities$name == "Augusta-Richmond" | usa.cities$name == "Knoxville" | usa.cities$name == "Roanoke"),]
  carolina.cities.labels <- carolina.cities$name
  carolina.cities.coords <- SpatialPoints(carolina.cities[5:4], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
}

###### Run batch from file

#speciesList <- read.delim(speciesList, header = FALSE, sep = "\n")
quartzFonts(avenir = c("Avenir Book", "Avenir Black", "Avenir Book Oblique","Avenir Black Oblique"))
par(mfrow = c(1,1), family = 'avenir')
#for (sciName in speciesList$V1){
  tryCatch({
  pdf(file = paste(sciName, "_", rasterResolution, ".pdf", sep = ""), width = 8.97, height = 6.7)
  print(maxent_map(sciName, minXlim, maxXlim, minYlim, maxYlim, modelStudyAreaOnly, mapStudyAreaOnly, futureModel))
  points(carolina.cities.coords, pch = 20, cex = 1)
  text(carolina.cities.coords, font = 2, labels = carolina.cities.labels, halo = TRUE, hw = 0.2, hc = "white", pos = 1, cex = 0.75)
  }, error=function(e){cat("ERROR: ", conditionMessage(e), "\n")})
  dev.off()
  #}

  print(sciName)
  print("**********************END***********************")
  rm(list = ls())
  quit()
