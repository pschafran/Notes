setwd("/Volumes/Samsung_T5/GIS/Isoetes")
library(maptools)
library(rgeos)
library(raster)
library(rgdal)
library(sf)

data("wrld_simpl")
usaHD <- getData(name="GADM", download = TRUE, level = 1, country = "USA")
canadaHD <- getData(name="GADM", download = TRUE, level = 1, country = "CAN")
greatLakes <- st_read("./Great_Lakes/Great_Lakes.shp")

states <- c("Texas", "Oklahoma", "Nebraska", "Kansas", "South Dakota", "Minnesota", "Iowa", "Missouri", "Arkansas", "Louisiana", "Wisconsin", "Illinois", "Kentucky", "Tennessee", "Mississippi", "Michigan", "Indiana", "Alabama", "Ohio", "Georgia", "Florida", "New York", "Pennsylvania", "West Virginia" ,"Maryland", "Virginia", "North Carolina", "South Carolina", "Delaware", "New Jersey", "Connecticut", "Rhode Island","Massachusetts", "New Hampshire", "Vermont","Maine")
provinces <- c("Ontario", "Quebec", "New Brunswick", "Nova Scotia")

us.states <- usaHD[usaHD$NAME_1 %in% states,]
us.states.simple <- gSimplify(us.states, tol = 0.01, topologyPreserve = TRUE)
can.provinces <- canadaHD[canadaHD$NAME_1 %in% provinces,]
can.provinces.simple <- gSimplify(can.provinces, tol = 0.01, topologyPreserve = TRUE)

{
  wrld_alt <- getData(name="SRTM", download = TRUE, lon = -81, lat = 34)
  ga_alt <- getData(name="SRTM", download = TRUE, lon = -88, lat = 34)
  tn_alt <- getData(name="SRTM", download = TRUE, lon = -88, lat = 36)
  va_alt <- getData(name="SRTM", download = TRUE, lon = -81, lat = 36)
  nc_alt <- getData(name="SRTM", download = TRUE, lon = -81, lat = 34)
  vae_alt <- getData(name="SRTM", download = TRUE, lon = -76, lat = 36)
  nce_alt <- getData(name="SRTM", download = TRUE, lon = -76, lat = 34)
  ms_alt <- getData(name="SRTM", download = TRUE, lon = -92, lat = 34)
  tx_alt <- getData(name="SRTM", download = TRUE, lon = -92, lat = 26)
  tx2_alt <- getData(name="SRTM", download = TRUE, lon = -98, lat = 26)
  tx3_alt <- getData(name="SRTM", download = TRUE, lon = -98, lat = 34)
  ok_alt <- getData(name="SRTM", download = TRUE, lon = -92, lat = 38)
  ok2_alt <- getData(name="SRTM", download = TRUE, lon = -98, lat = 38)
  fl_alt <- getData(name="SRTM", download = TRUE, lon = -82, lat = 26)
  pa_alt <- getData(name="SRTM", download = TRUE, lon = -76, lat = 42)
  ny_alt <- getData(name="SRTM", download = TRUE, lon = -76, lat = 50)
  on_alt <- getData(name="SRTM", download = TRUE, lon = -84, lat = 50)
  qu_alt <- getData(name="SRTM", download = TRUE, lon = -72, lat = 50)
  vt_alt <- getData(name="SRTM", download = TRUE, lon = -72, lat = 42)
  me_alt <- getData(name="SRTM", download = TRUE, lon = -68, lat = 42)
  nb_alt <- getData(name="SRTM", download = TRUE, lon = -68, lat = 50)
  ns_alt <- getData(name="SRTM", download = TRUE, lon = -62, lat = 48)
  ns2_alt <- getData(name="SRTM", download = TRUE, lon = -62, lat = 42)
  oh_alt <- getData(name="SRTM", download = TRUE, lon = -82, lat = 43)
  il_alt <- getData(name="SRTM", download = TRUE, lon = -88, lat = 43)
  wi_alt <- getData(name="SRTM", download = TRUE, lon = -88, lat = 47)
  ia_alt <- getData(name="SRTM", download = TRUE, lon = -93, lat = 43)
  mn_alt <- getData(name="SRTM", download = TRUE, lon = -93, lat = 47)
  neb_alt <- getData(name="SRTM", download = TRUE, lon = -97, lat = 43)
  nd_alt <- getData(name="SRTM", download = TRUE, lon = -97, lat = 47)
  }

{
  piedNCpts <- read.csv("sch86_141-144_GPS.csv", header = FALSE, sep = ",")
  sandhills <- read.csv("sch_Sandhills.csv", header = FALSE, sep = ",")
  coast <- read.csv("sch_coast.csv", header = FALSE, sep = ",")
  georgianaSL <- read.csv("sch_georgianaSL.txt", header = FALSE, sep = ",")
  leary <- read.csv("sch_leary.csv", header = FALSE, sep = ",")
  allIsoetes <- read.csv("All_Isoetes.csv", header = TRUE, sep = "\t")
  appaN <- read.csv("appaN.csv", header = FALSE, sep = ",")
  appaS <- read.csv("appaS.csv", header = FALSE, sep = ",")
  sept <- read.csv("sept.csv", header = FALSE, sep = "\t")
}

grayscale <- gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2)
southeastAlt <- mosaic(wrld_alt, ga_alt, tn_alt, va_alt, vae_alt, nce_alt, ms_alt, tx_alt, fl_alt, tx2_alt, tx3_alt, ok_alt, ok2_alt, fun=mean)
plot(southeastAlt, col = grayscale)

{
eastCoastAlt <- mosaic(wrld_alt, ga_alt, tn_alt, va_alt, vae_alt, nce_alt, ms_alt, pa_alt, ny_alt, tx_alt, fl_alt, tx2_alt, tx3_alt, ok_alt, ok2_alt, qu_alt, on_alt, vt_alt, nb_alt, me_alt, ns_alt, ns2_alt,oh_alt, il_alt, wi_alt,ia_alt,neb_alt,nd_alt, mn_alt, fun=mean)
plot(eastCoastAlt, col = grayscale)
plot(wrld_simpl, add = TRUE)
points(allIsoetes$Longitude, allIsoetes$Latitude, col = "red", pch = 20, cex = 0.75)
points(piedNCpts$V2 , piedNCpts$V1, col = "blue", pch = 20 )
points(sandhills$V2 , sandhills$V1, col = "dark green", pch = 20 )
points(coast$V2, coast$V1, col = "purple" , pch = 20)
points(georgianaSL$V2, georgianaSL$V1, col = "light blue" , pch = 20)
points(leary$V2, leary$V1, col = "green" , pch = 20)
points(appaN$V2, appaS$V1, col = "gold" , pch = 20)
points(appaS$V2, appaN$V1, col = "gold" , pch = 20)
points(sept$V2, sept$V1, col = "purple" , pch = 20)
}

plot(us.states)
plot(can.provinces, add = TRUE)


## Mapping sample points with 'sf' and 'spData' packages
setwd("/Volumes/Samsung_T5/GIS/Isoetes")
library(sf)
library(spData)
library(raster)
library(spDataLarge)

data("us_states", package = "spData")
data("world", package = "spData")

# Load points from above

plot(st_geometry(south[south$NAME == "Alabama", ]), expandBB = c(0.5,1.25,2,1.5), axes = TRUE)
plot(st_geometry(us_states), add = TRUE)


{
  points(allIsoetes$Longitude, allIsoetes$Latitude, col = "black", cex = 0.75)
  points(piedNCpts$V2 , piedNCpts$V1, col = "black")
  points(sandhills$V2 , sandhills$V1, col = "black")
  points(coast$V2, coast$V1, col = "black")
  points(georgianaSL$V2, georgianaSL$V1, col = "black")
  points(leary$V2, leary$V1, col = "black")
  points(appaN$V2, appaN$V1, col = "black")
  points(appaS$V2, appaS$V1, col = "black")
  points(sept$V2, sept$V1, col = "black")
  
}
{
  points(allIsoetes$Longitude, allIsoetes$Latitude, col = "red", pch = 20, cex = 0.75)
  points(piedNCpts$V2 , piedNCpts$V1, col = "blue", pch = 20 )
  points(sandhills$V2 , sandhills$V1, col = "dark green", pch = 20 )
  points(coast$V2, coast$V1, col = "purple" , pch = 20)
  points(georgianaSL$V2, georgianaSL$V1, col = "light blue" , pch = 20)
  points(leary$V2, leary$V1, col = "green" , pch = 20)
  points(appaN$V2, appaN$V1, col = "gold" , pch = 20)
  points(appaS$V2, appaS$V1, col = "gold" , pch = 20)
  points(sept$V2, sept$V1, col = "pink" , pch = 20)

}
{
piedNcpts_hull <- chull(piedNCpts)
sandhills_hull <- chull(sandhills)
coast_hull <- chull(coast)
georgiana_hull <- chull(georgianaSL)
leary_hull <- chull(leary)
sept_hull <- chull(sept)
appaN_hull <- chull(appaN)
appaS_hull <- chull(appaS)
}


# Run whole plot at once
{
  par(mar = c(5,5,2,2))
  plot(us.states.simple, axes = TRUE, xlim = c(-95,-65), ylim = c(27,47))
  plot(can.provinces.simple, add = TRUE)
  
  points(allIsoetes$Longitude, allIsoetes$Latitude, col = "black", bg = "red", pch = 21, cex = 0.5)
  points(piedNCpts$V2 , piedNCpts$V1, col = "black", bg = "blue", pch = 21 )
  points(sandhills$V2 , sandhills$V1, col = "black", bg = "dark green", pch = 21 )
  points(coast$V2, coast$V1, col = "black", bg = "purple" , pch = 21)
  points(georgianaSL$V2, georgianaSL$V1, col = "black", bg = "light blue" , pch = 21)
  points(leary$V2, leary$V1, col = "black", bg = "green" , pch = 21)
  points(appaN$V2, appaN$V1, col = "black", bg = "gold" , pch = 21)
  points(appaS$V2, appaS$V1, col = "black", bg = "gold" , pch = 21)
  points(sept$V2, sept$V1, col = "black", bg = "brown1" , pch = 21)
  
  legend("bottomright",
         legend = c("I. appalachiana", "I. georgiana s.l.", "I. hyemalis 'NC Piedmont'", "I. hyemalis 'NC Sandhills'", "I. 'Leary'", "I. microvela + I. riparia s.l.", "I. septentrionalis", "Other taxa"),
         col = "black",
         pt.bg = c("gold", "light blue", "blue", "dark green", "green", "purple", "brown1", "red"),
         bty = "n",
         pch = c(21, 21, 21, 21, 21, 21, 21, 21),
         pt.cex = c(1,1,1,1,1,1,1,0.5)
         )
  plot(greatLakes, add = TRUE, col = "white")
}

# Plot sampling locations
{
  par(mar = c(5,5,2,2))
  plot(us.states.simple, axes = TRUE, xlim = c(-95,-65), ylim = c(27,47))
  plot(can.provinces.simple, add = TRUE)
  
  points(allIsoetes$Longitude, allIsoetes$Latitude, col = "black", bg = "red", pch = 21, cex = 1)
  plot(greatLakes, add = TRUE, col = "white")
  
}













###### Mantel Tests comparing genetic distance and geographic distance
{
   setwd("/Volumes/Samsung_T5/GIS/Isoetes")
  library(sp)
  library(vegan)
  library(ape)
}
polyploidIsoetes <- as.matrix(read.csv("Polyploid_Isoetes_GPS.csv", header = TRUE, sep = ",", row.names = 1))
polyploidIsoetes <- as.matrix(polyploidIsoetes)
polyploidIsoetesMatrix <- SpatialPoints(polyploidIsoetes)
isoetesDists <- spDists(polyploidIsoetes, longlat = TRUE)
View(isoetesDists)

geneticDist <- as.matrix(read.csv("Polyploid_reduced_for_mantel.csv", header = TRUE, sep = ",", row.names = 1))

mantel(isoetesDists, geneticDist, permutations = 1000)
mantel.test(isoetesDists, geneticDist, nperm = 1000, graph = TRUE)


polyploidIsoetesLatitude <- polyploidIsoetes[,1]
polyploidIsoetesLatitude <- cbind(polyploidIsoetesLatitude, 0)
polyploidIsoetesLatitudePoints <- SpatialPoints(polyploidIsoetesLatitude)
polyploidLatitudeDists <- spDists(polyploidIsoetesLatitude, longlat = TRUE)
 
mantel(polyploidLatitudeDists, geneticDist, permutations = 1000)
mantel.test(polyploidLatitudeDists, geneticDist, nperm = 1000, graph = TRUE)
polyploidLatitudes.correlog <- mantel.correlog(geneticDist, polyploidLatitudeDists, nperm = 10000, n.class = 30)
plot(polyploidLatitudes.correlog, alpha = 0.05, main = "Mantel test of polyploid latitude")

##### Southeast Tetraploids
seTetraploidGeographicDist <- as.matrix(read.csv("SE_tetraploids_Isoetes_GPS.csv",header = TRUE, sep = ",", row.names = 1 ))
seTetraploidGeographicDist <- spDists(seTetraploidGeographicDist, longlat = TRUE)

seTetraploidGeneticDist <- as.matrix(read.csv("SE_tetraploids_geneticDist.csv", header = TRUE, sep = ",", row.names = 1))

mantel.test(seTetraploidGeneticDist, seTetraploidGeographicDist, nperm = 10000, graph = TRUE)
mantel(seTetraploidGeographicDist, seTetraploidGeneticDist, permutations = 10000)
seTetraploids.correlog <- mantel.correlog(seTetraploidGeneticDist, seTetraploidGeographicDist, nperm = 10000)
plot(seTetraploids.correlog, alpha = 0.05, main = "Mantel test of tetraploids in SE US")


##### All Southeast Polyploids
sePolyploidGeographicDF <- as.data.frame(read.csv("SE_polyploids_Isoetes_GPS.csv",header = TRUE, sep = ",", row.names = 1 ))
sePolyploidGeographicDist <- spDists(as.matrix(sePolyploidGeographicDF), longlat = TRUE)

# To remove duplicate coordinates
sePolyploidGeographicDist <- cbind("Isoetes", sePolyploidGeographicDF)
flags <- clean_coordinates(sePolyploidGeographicDist, species = "Isoetes", lon = "Longitude", lat = "Latitude", 
                           tests = c("capitals", "centroids", "duplicates", "equal", "gbif", "institutions", "seas", "zeros"))
sePolyploidGeographicDistClean <- sePolyploidGeographicDist[flags$.dpl,]
sePolyploidGeographicDistClean <- sePolyploidGeographicDistClean[,2:3]
sePolyploidGeographicDistClean <- spDists(as.matrix(sePolyploidGeographicDistClean), longlat = TRUE)

sePolyploidsGeneticDist <- as.matrix(read.csv("SE_polyploids_geneticDist.csv", header = TRUE, sep = ",", row.names = 1))
sePolyploidsGeneticDist <- as.matrix(read.csv("SE_polyploids_geneticDist_noDups.csv", header = TRUE, sep = ",", row.names = 1))
dim(sePolyploidsGeneticDist)
dim(sePolyploidGeographicDistClean)

mantel.test(sePolyploidsGeneticDist, sePolyploidGeographicDistClean, nperm = 10000, graph = TRUE)
mantel(sePolyploidGeographicDistClean, sePolyploidsGeneticDist, permutations = 10000)
sePolyploids.correlog <- mantel.correlog(sePolyploidsGeneticDist, sePolyploidGeographicDistClean, nperm = 10000, n.class = 30)
plot(sePolyploids.correlog, alpha = 0.05, main = "Mantel test of polyploid in SE US")


View(sePolyploidsGeneticDist)

fit <- cmdscale(sePolyploidsGeneticDist,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", pch = 20)

fit2 <- cmdscale(sePolyploidGeographicDist, eig = TRUE, k = 2)
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
plot(x2, y2, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", pch = 20)
heatmap(sePolyploidsGeneticDist)
NMDS <- metaMDS(sePolyploidsGeneticDist, k = 3, try = 1000, trymax = 1000)
NMDS <- metaMDS(geneticDist, k = 2, try = 1000, trymax = 1000)
# NMDS with points colored by clades in NJ tree for southeastern dataset
plot(NMDS$points[,1], NMDS$points[,2])
points(NMDS$points[1:8,1],NMDS$points[1:8,2], pch = 21, bg = "gold", col = "black")
points(NMDS$points[73,1],NMDS$points[73,2], pch = 21, bg = "gold", col = "black")
points(NMDS$points[75,1],NMDS$points[75,2], pch = 21, bg = "gold", col = "black")
samples <- cbind(samples, "appalachiana")
points(NMDS$points[17:20,1],NMDS$points[17:20,2], pch = 21, bg = "orange", col = "black")
points(NMDS$points[37,1],NMDS$points[37,2], pch = 21, bg = "green", col = "black")
points(NMDS$points[39,1],NMDS$points[39,2], pch = 21, bg = "green", col = "black")
points(NMDS$points[63:67,1],NMDS$points[63:67,2], pch = 21, bg = "green", col = "black")
points(NMDS$points[55:58,1],NMDS$points[55:58,2], pch = 21, bg = "blue", col = "black")
points(NMDS$points[38,1],NMDS$points[38,2], pch = 21, bg = "dark green", col = "black")
points(NMDS$points[45:46,1],NMDS$points[45:46,2], pch = 21, bg = "dark green", col = "black")
points(NMDS$points[48,1],NMDS$points[48,2], pch = 21, bg = "dark green", col = "black")
points(NMDS$points[125:126,1],NMDS$points[125:126,2], pch = 21, bg = "dark green", col = "black")
points(NMDS$points[36,1],NMDS$points[36,2], pch = 21, bg = "brown", col = "black")
points(NMDS$points[23,1],NMDS$points[23,2], pch = 21, bg = "brown", col = "black")
points(NMDS$points[28,1],NMDS$points[28,2], pch = 21, bg = "brown", col = "black")
points(NMDS$points[21:22,1],NMDS$points[21:22,2], pch = 21, bg = "light blue", col = "black")
points(NMDS$points[23:27,1],NMDS$points[23:27,2], pch = 21, bg = "light blue", col = "black")
points(NMDS$points[29:31,1],NMDS$points[29:31,2], pch = 21, bg = "light blue", col = "black")
points(NMDS$points[10:11,1],NMDS$points[10:11,2], pch = 21, bg = "light blue", col = "black")
points(NMDS$points[92,1],NMDS$points[92,2], pch = 21, bg = "light blue", col = "black")
points(NMDS$points[49:52,1],NMDS$points[49:52,2], pch = 21, bg = "purple", col = "black")
points(NMDS$points[54,1],NMDS$points[54,2], pch = 21, bg = "purple", col = "black")
points(NMDS$points[86:91,1],NMDS$points[86:91,2], pch = 21, bg = "purple", col = "black")
points(NMDS$points[9,1],NMDS$points[9,2], pch = 21, bg = "purple", col = "black")
points(NMDS$points[60:62,1],NMDS$points[60:62,2], pch = 21, bg = "red", col = "black")
points(NMDS$points[97:102,1],NMDS$points[97:102,2], pch = 21, bg = "red", col = "black")
points(NMDS$points[67,1],NMDS$points[67,2], pch = 21, bg = "red", col = "black")
points(NMDS$points[110:118,1],NMDS$points[110:118,2], pch = 21, bg = "red", col = "black")
points(NMDS$points[77:81,1],NMDS$points[77:81,2], pch = 21, bg = "red", col = "black")
points(NMDS$points[85,1],NMDS$points[85,2], pch = 21, bg = "red", col = "black")
points(NMDS$points[121:123,1],NMDS$points[121:123,2], pch = 21, bg = "red", col = "black")

# NMDS with points colored by taxonomy for southeastern dataset
plot(NMDS$points[,1], NMDS$points[,2])
points(NMDS$points[1:8,1],NMDS$points[1:8,2], pch = 21, bg = "gold", col = "black")
points(NMDS$points[73,1],NMDS$points[73,2], pch = 21, bg = "gold", col = "black")
points(NMDS$points[75,1],NMDS$points[75,2], pch = 21, bg = "gold", col = "black")
points(NMDS$points[10:13,1],NMDS$points[10:13,2], pch = 21, bg = "orange", col = "black")
points(NMDS$points[21:31,1],NMDS$points[21:31,2], pch = 21, bg = "red", col = "black")
points(NMDS$points[33:59,1],NMDS$points[33:59,2], pch = 21, bg = "green", col = "black")
points(NMDS$points[60:62,1],NMDS$points[60:62,2], pch = 21, bg = "blue", col = "black")
points(NMDS$points[77:85,1],NMDS$points[77:85,2], pch = 21, bg = "purple", col = "black")
points(NMDS$points[86:91,1],NMDS$points[86:91,2], pch = 21, bg = "magenta", col = "black")
points(NMDS$points[95:107,1],NMDS$points[95:107,2], pch = 21, bg = "gray", col = "black")
points(NMDS$points[110:118,1],NMDS$points[110:118,2], pch = 21, bg = "black", col = "black")


# NMDS with points colored by NJ clades for whole dataset
plot(NMDS$points[,1], NMDS$points[,2])
points(NMDS$points[appalachiana,1],NMDS$points[appalachiana,2], pch = 21, bg = "gold", col = "black")
points(NMDS$points[septentrionalis,1],NMDS$points[septentrionalis,2], pch = 21, bg = "red", col = "black")
points(NMDS$points[tuckermanii,1],NMDS$points[tuckermanii,2], pch = 21, bg = "blue", col = "black")
points(NMDS$points["I_acadiensis_Schafran175-pop",1],NMDS$points["I_acadiensis_Schafran175-pop",2], pch = 21, bg = "white", col = "black")
points(NMDS$points[snowii,1],NMDS$points[snowii,2], pch = 21, bg = "magenta", col = "black")
points(NMDS$points[maritima,1],NMDS$points[maritima,2], pch = 21, bg = "purple", col = "black")
points(NMDS$points[microvela,1],NMDS$points[microvela,2], pch = 21, bg = "green", col = "black")
points(NMDS$points[georgiana2,1],NMDS$points[georgiana2,2], pch = 21, bg = "light blue", col = "black")
points(NMDS$points[georgiana,1],NMDS$points[georgiana,2], pch = 21, bg = "light green", col = "black")
points(NMDS$points[sandhills,1],NMDS$points[sandhills,2], pch = 21, bg = "tan", col = "black")
points(NMDS$points[piedmont,1],NMDS$points[piedmont,2], pch = 21, bg = "dark green", col = "black")
text(x = NMDS$points[,1], y = NMDS$points[,2], labels = rownames(NMDS$points), cex = 0.5)

appalachiana <- c("I_appalachiana_Cressler8-pop","I_appalachiana_Schafran105-pop", "I_appalachiana_Schafran148-pop" ,"I_appalachiana_Schafran150", "I_appalachiana_Schafran178" ,"I_appalachiana_Schafran199-pop", "I_appalachiana_Schafran200","I_appalachiana_Schafran201", "I_louisianensis_Brunton17581", "I_louisianensis_Schafran106", "I_septentrionalis_Schafran152-pop", "I_septentrionalis_Schafran153")
septentrionalis <- c("I_septentrionalis_Brunton15341","I_septentrionalis_Brunton19142","I_septentrionalis_Schafran151-pop","I_septentrionalis_Schafran170","I_septentrionalis_Schafran171-1","I_septentrionalis_Schafran172","I_septentrionalis_Schafran173-pop","I_tuckermanii_Schafran160-pop","I_riparia_Schafran159","I_riparia_Schafran161-pop","I_riparia_Taylor6706","I_X_dodgei_Brunton19143","I_X_eatonii_Taylor6750")
tuckermanii <- c("I_tuckermanii_Schafran166-pop","I_tuckermanii_Schafran168","I_tuckermanii_Schafran176-pop","I_tuckermanii_Taylor6707","I_riparia_Taylor6675","I_X_harveyi_Taylor6677","I_X_heterospora_Taylor6676","I_acadiensis_Schafran175-pop")
snowii <-c("I_snowii_Schafran78-2","I_snowii_Schafran80-10","I_snowii_Schafran80-11","I_snowii_Schafran80-7","I_snowii_Schafran80-9","I_snowii_SchafranGA02","I_snowii_SchafranGA05","I_snowii_SchafranGA06_rep3","I_snowii_SchafranGA12","I_melanopoda_Ciafre256-1","I_melanopoda_Ciafre256-2","I_melanopoda_Ciafre728-1","I_melanopoda_Ciafre728-2","I_melanopoda_Ciafre728-3","I_junciformis_Bolin","I_junciformis_Brunton17608","I_junciformis_Schafran104","I_Leary_Schafran114","I_sp_UnknownChickahominy2_rep1","I_sp_UnknownChickahominy2_rep2","I_sp_UnknownChickahominy5","I_piedmontana_Schafran101-1","I_piedmontana_Schafran101-2","I_piedmontana_Schafran101-3","I_piedmontana_Schafran102-3","I_piedmontana_Schafran103-1","I_piedmontana_Schafran103-2")   
maritima <- c("I_sp_Taylor-1","I_sp_Taylor-3","I_bolanderi_X_occidentalis_Taylor6756","I_X_herbwagneri_TaylorSN-pop","I_maritima_Taylor6983-pop","I_maritima_WoodbridgeSN2","I_maritimaXechinospora_Taylor6988-2_2","I_maritimaXechinospora_Taylor6988-3","I_splayedLeaves_Taylor6990-1","I_straightLeaves_Taylor6989-3","I_curledLeaves_Taylor6991-1_2","I_occidentalis_Taylor6755","I_occidentalis_WoodbridgeSN1-pop")
microvela <- c("I_microvela_BolinJB40NC","I_microvela_BolinJBNC199EO2","I_microvela_BolinJBNC200EO3-pop","I_microvela_BolinJBNC201EO4-all","I_microvela_BolinJBNC202EO5-all","I_microvela_MatthewsI09-35","I_hyemalis_Schafran134-pop","I_hyemalis_Schafran135","I_hyemalis_Schafran136","I_hyemalis_Schafran137-pop","I_hyemalis_Schafran138","I_hyemalis_Schafran139-pop","I_appalachiana_X_hyemalis_Brunton19011B","I_Laurentiana_Brunton20077","I_Laurentiana_Brunton20087","I_Laurentiana_Brunton20101-2","I_riparia_PotomacCreek")
georgiana2 <- c("I_georgiana_Cressler10","I_georgiana_Cressler11-1","I_georgiana_Schafran111-pop","I_georgiana_Schafran112","I_georgiana_SchafranGA17","I_georgiana_SchafranGA18","I_boomii_Leonard12408","I_georgiana_Taylor6769")
georgiana <- c("I_georgiana_Matthews3-pop","I_georgiana_SchafranGA16","I_hyemalis_Schafran107")
sandhills <- c("I_Uwharrie_Schafran76-2","I_Uwharrie_Schafran76-3","I_Uwharrie_SchafranSN","I_hyemalis_Schafran118-pop","I_hyemalis_Schafran133-pop","I_hyemalis_Schafran129-pop","I_hyemalis_Schafran130-1")
piedmont <- c("I_hyemalis_Schafran141-pop","I_hyemalis_Schafran142-pop","I_hyemalis_Schafran143","I_hyemalis_Schafran144","I_piedmontana_BolinJBNC17-3","I_piedmontana_Schafran86-2")




