#---
#title: "Tutorial Part I Basics - R goes spatial"
#author: Stephanie Kramer-Schadt
#date: Biomove GIS Course, 02/2016
#output: html_document
#
#---
#
### Thanks...
#...to Rafael Tietz and Jürgen Niedballa for providing data
#and support in compiling this manual.
#
#---
#
### Why R for Spatial Analysis?
#
#- free software
#- number of contributed packages
#  for spatial data handling and analysis
#- capacity to analyse and visualize large data
#- easily reproducable results (no endless clicking!)
#- can be used in conjunction with a GIS
#
#---
#
### **Packages**
#
#In **R** you can use some standard commands like  *apply*, *length* or *sqrt* included in the *base* package.
#
#But if you want to work with geospatial data, it is necessary to install and
#load **packages** with *classes* and *methods* for spatial data.
#
#---
#
### Some important packages I
#
#1. **sp** or  **sf** (spatial overlay) ## sf update 2019
#  * for plotting data as maps, spatial selection, retrieving coordinates
#    or subsetting (points, lines, polygons and rasters)
#2. **raster** (raster manipulation)
#  * reading, writing, manipulating, analyzing and moddeling
#    of gridded spatial data
#  * supports processing of very large files
#3. **rgdal** - *R Geospatial Data Abstraction Library*
#  * GDAL raster and OGR vector map data can be imported and exported
#  * variety of useful commandline utilities for data translation and processing
#  * useful for projections
#4. **tmap** ## sf update 2019
#  * To create static and interactive publication ready graphics  
#
#---
#
### Some important packages II
#
#4. **maptools**
#  * set of tools for manipulating and reading geographic data
#5. **rgeos**
#  * topology operations on geometries (e.g., Union, Intersection,...)
#  * contains some typical GIS operations
#6. **dismo** or **biomod2**
#  * functions for species distribution modeling
#7. **rgl**
#  * provides medium to high-level functions for 3D interactive graphics
#8. **geosphere**: spherical trigonometry
#
#---
#
#### How to install...
#  install.packages("sp")
#  install.packages("dismo")
#  install.packages("raster")
#  install.packages("GISTools")
#  install.packages("rgdal")
#  install.packages("maptools")
#  install.packages("rgeos")
#  install.packages("rgl")
#  install.packages("rasterVis")
#  install.packages("ellipse")
#  install.packages(pkgs=c("CircStats", "deSolve", "coda", "deldir",
#  "igraph", "RandomFields", "ks"))
#  install.packages("sf") ## sf update 2019
#  install.packages("tmap") ## sf update 2019

#### ...and load packages
   library(sp)
   library(dismo)
   library(raster)
   library(GISTools)
   library(rgdal)
   library(maptools)
   library(rgeos)
   library(rgl)
   library(rasterVis)
   library(sf) ## sf update 2019
   library(tmap) ## sf update 2019
   
#
#### How to call help
#- Information about a function or a package
   ?dismo
#- search for an item
   ??SpatialPolygons
#- show the instructions of a package
   vignette(package="dismo") # vignette(package="sp")
#- load the instruction of a package (embedded pdf)
 vignette("sdm") #vignette("over")
#
### Useful (web)sites
#
#1. R news and tutorials
#* http://www.r-bloggers.com/
#* http://www.inside-r.org/
#
#2. Quick introduction to spatial R
#* http://pakillo.github.io/R-GIS-tutorial/
#* http://rstudio-pubs-static.s3.amazonaws.com/7993_6b081819ba184047802a508a7f3187cb.html
#
#3. Spatial analysis, visualisation and resources
#* Introduction to Spatial Data and ggplot2
#* http://spatial.ly/
#
#4. overview over spatial packages
#* http://cran.r-project.org/web/views/Spatial.html
#* http://cran.r-project.org/view=Spatial
#* https://www.r-spatial.org/  ## sf update 2019
#
# 
#5. R-cheatsheets in your folder
#
#6. Infos about simple features in R (sf package) https://r-spatial.github.io/sf/index.html   ## sf update 2019
#7. Use EPSG codes as unique and specific identifies of your coordinate reference systsem instead of writing projection details. 
#   See http://spatialreference.org/ for EPSG codes   ## sf update 2019
#---

 #
### Further reading
#
#1. Check Roger Bivand's publications:
#* http://cran.r-project.org/web/views/Spatial.html
#* http://www.csiss.org/learning_resources/index.html
#* Roger Bivand et al. 2008, Appied Spatial Data Analysis in R, Springer
#
#2. Check Edzer Pebesma's course materials and publications:
#* http://ifgi.uni-muenster.de/~epebe_01/
#* https://edzer.github.io/UseR2017/  ## sf update 2019
#* https://www.rstudio.com/resources/videos/tidy-spatial-data-analysis/  ## sf update 2019
#
#3. Geocomputation with R a book of Robin Lovelace, Jakub Nowosad, Jannes Muenchow
#* https://geocompr.robinlovelace.net
 
# ##############################################################################
## Let us start
################################################################################

### Assign the workspace

# adjust the working directory [wd]
work_wd <- "D:/OneDrive/_BioMoveGISKurs/_BioMove_GISIntro_2019/Geodata"

# relative to wd
raster_wd <- paste(work_wd,"/","rasterfiles",sep='')
shapes_wd <- paste(work_wd,"/","shapefiles",sep='')
output_wd <- paste(work_wd,"/","output",sep='')

# assign the workspace and the file
setwd(raster_wd)
myfilenames <- list.files() #store filenames in object
head(myfilenames) #have a look at the first 6

### Some Vector data to work with (shapefiles folder)
#### Polygons - area/ perimeter + attribute table [a.t.]
#- borneo_admin.shp - country borders on Borneo
#- Bor_PA.shp - protected areas on Borneo
#
#### Spatial Line data - length + [a.t.]
#- sn_100000.shp       - main river net Borneo
#
#### Spatial Point data - georeferenced XY coordinate + [a.t.]
#- various, e.g. DHOsim.csv  -  species occurrence records
#- KettleHoles.txt - textfile with
#
#---
#
### Some Raster data (rasterfiles folder) to work with
#
#- Landsat_GLS2010_example - satellite image GeoTiff
#- Borneo_MAT  - mean annual temp
#- Borneo_DEM  - DGM
#- Borneo_LC  - land use categories (17)
#
###### Land use categories I (Borneo_LC)
#1. Lowland forest 0 - 500m
#2. Upland forest 501 - 1000m
#3. Lower montane forest 1001 - 1500m
#4. Upper montane forest > 1500m
#5. Lowland forest (fragmented or degraded)
#
#---
#
#### Raster data (Borneo_LC)
###### Land use categories II continued
#6. Upland forest (fragmented or degraded)
#7. Lower montane forest (fragmented or degraded)
#8. Upper montane forest (fragmented or degraded)
#9. Swamp forest
#10. Mangrove
#11. Old plantations
#12. Young plantation and agriculture
#13. Burned forest area
#14. Mixed crops
#15. Water and fishponds
#16. Water bodies
#17. No data
#

################################################################################
# Handling and visualising raster data
################################################################################
## Data import - Load the Raster from ascii

# transform the ascii data to a raster
Bor_mat <- raster(paste(raster_wd,'Borneo_MAT.asc',sep='/'))
#   Bor_mat   

# assign the projection (coordinate reference system)
Bor_mat@crs <- CRS("+proj=longlat +datum=WGS84")
Bor_mat@crs <- CRS("+init=epsg:4326") ## different way to specify. See http://spatialreference.org/ for EPSG codes   ## sf update 2019

Bor_dem <-raster(paste(raster_wd,'Borneo_DEM.asc',sep='/'))
Bor_dem@crs <- crs(Bor_mat)

### Before we start: clip a small area
extent(Bor_mat)
clip_extent <- extent(117.2,117.45,5.4,5.5)
Bor_mat_cr <- crop(Bor_mat, clip_extent)

### working with the 'raster' class I
#Have a look at the internal data structure
#of the raster object:
Bor_mat_cr #str(Bor_mat_cr) #attributes(Bor_mat_cr)
   Bor_mat_cr@data@values
   coordinates(Bor_mat_cr)

### working with the 'raster' class II
#Retrieve internal data and access single bits of information
bbox(Bor_mat_cr)
Bor_mat_cr@extent@xmin # = xmin(Bor_mat_cr) = xmin(extent(Bor_mat))
#Bor_mat_cr@extent = extent(Bor_mat_cr) = attr(Bor_mat_cr,'extent')
Bor_mat_cr@ncols        # = ncol(Bor_mat_cr) !n.b. ncols vs. ncol!
Bor_mat_cr@data@values  # = values(Bor_mat_cr)
Bor_mat_cr@crs #=crs(Bor_mat_cr) =projection(Bor_mat_cr) !=CRS(Bor_mat_cr)

# crs = coordinate reference system
# CRS creates projection and sets arguments for crs!
# e.g. CRS("+proj=longlat +datum=WGS84") OR CRS("+init=epsg:4326") ## sf update 2019


### working with the 'raster' class III
#Retrieve internal data cont.:
nrow(Bor_mat_cr) #ncol(Bor_mat_cr) # ncell(Bor_mat_cr)
dim(Bor_mat_cr) # 12 rows, 30 columns, 1 z-dimension
res(Bor_mat_cr) #resolution = cell size
coordinates(Bor_mat_cr) #command from package 'sp'

## Visualising (plotting) raster data
### Difference between *plot* and *image*
# if you *plot* a raster, the proportions will always be constant
plot(Bor_mat)
#   zoom(Bor_mat)
#plot(Bor_mat_cr)

## sf update 2019 - start

library(viridis) ## the colour map viridis provides clear patterns for all known perception problems 
#* https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

## plot map as interactive web map
tmap_mode("view")
(i.m <- tm_shape(Bor_mat) +
  tm_raster(palette = "viridis",  title = "Global Land Cover"))
save_tmap(i.m, paste0(output_wd,"/BorneoMap_4326.html"))

## plot map for print 
tmap_mode("plot")
(m <- tm_shape(Bor_mat) +
  tm_raster(palette = "viridis",  title = "Global Land Cover"))
save_tmap(m, paste0(output_wd,"/BorneoMap_4326.png"), units = "mm", width = 70, height = 70, dpi = 150 )

## sf update 2019 - end

#### Image (raster)
# if you use *image*, the proportions will be changed
image(Bor_mat)
# image(Bor_mat_cr)

### perspective: 3D plot
#- Cool 3D plots with rgl library, e.g. 'rgl.surface'
#- Why is it so spiny??? Check units!
persp(Bor_dem, xlab="Easting", ylab="Northing",
      zlab="elevation", main="Elevation model of Borneo",
      r=1, d=1.5, expand=0.7, ticktype="detailed")
# Solution: hillshade (see below)

### Project the Grid  (the small one)
proj_moll=CRS("+proj=moll +lat_0=65 +lon_0=10") # Mollweide
Bor_dem_moll <-projectRaster(Bor_mat_cr,crs=proj_moll)
persp(Bor_dem_moll, xlab="Easting", ylab="Northing",
      zlab="elevation", main="Elevation model of Borneo",
      r=1, d=5.5, expand=0.1, ticktype="detailed")


### Working with stacks I
#A raster stack is a collection of RasterLayer objects
#with the same spatial extent and resolution
myfiles <- list.files(path= raster_wd, pattern='.asc$',
                       full.names=TRUE )
myfiles
predictors <- stack(myfiles)

### Working with stacks II
# We can assign each single raster a projection with just one command!
# n.b.! - predictors@proj4string exchanged  by projection(predictors)
projection(predictors)
projection(predictors) <- CRS("+proj=longlat +datum=WGS84")
# n.b. example with capital CRS
cellStats(predictors, 'mean')
round(cellStats(predictors, 'mean'), digits=2)

### plot stack
plot(predictors)

### Bw-plot from raster stack
bwplot(predictors[[c(1,3)]]) #bwplot(predictors)
#bwplot(predictors[[c(3,3)]]) #bwplot(predictors)
   
### Create a hillshade (3D look)
#- simulates a 3D surface
#- computes shaded relief values for a raster surface
slope <- terrain(x = Bor_dem, opt = "slope",
                    unit = "radians", neighbors = 8)
aspect <- terrain(x = Bor_dem, opt = "aspect",
                     unit = "radians", neighbors = 8)
Bor_hs <- hillShade(slope, aspect, angle = 45,
                      direction = 270)

### Plot the hillshade
plot(Bor_hs, col=grey(0:100/100), legend=FALSE)
#   zoom(Bor_hs)
#and colors on top: alpha value gives semi-transparency
plot(Bor_dem, col=terrain.colors(25,alpha=0.3),add=T)
points(coordinates(Bor_mat_cr),cex=0.1,pch='+')
plot(extent(Bor_mat_cr),add=T, col='red')

### Solution to plotting problems
#- stretch the plot to see problem
#- problems with add = TRUE:  http://stackoverflow.com/questions/24213453/overlay-raster-plot-using-plot-add-t-leads-to-arbitrary-misalignment-of-fin
image(Bor_hs, col = gray(seq(0.1, 1, length.out = 100)), asp =  nrow(Bor_hs) / ncol(Bor_hs))
plot(Bor_dem, alpha = 0.6, col = terrain.colors(n = 100), add = TRUE)

## Data export
### Save the hillshade raster (not the plot...)
#- define the workspace where to save with setwd() - but remember that you changed the working directory with that!
#- Better: store path within the file name
#save as ascii-format
writeRaster(x=Bor_hs, filename=paste(output_wd, "my_Bor_hs", sep = "/"),
            format="ascii",overwrite=TRUE,NAflag = -9999)
#save as GeoTiff
writeRaster(x = Bor_hs,filename = paste(output_wd, "my_hillsh", sep = "/"),
            format = "GTiff")

## Manipulating rasters
vignette('Raster') #vignette(package='raster')
#Do any kind of raster algebra
new_ras <- Bor_mat + Bor_dem + 100
plot(new_ras,col=rev(rainbow(25)))

### Convert temperature to °C (divide by 10)
range(values(Bor_mat), na.rm = TRUE)
Bor_mat_cel <- Bor_mat / 10
range(values(Bor_mat_cel), na.rm = TRUE)
# find areas with mean annual temp >= 25°C
mean.t.c.25 <- Bor_mat_cel >= 25
plot(mean.t.c.25)

#### Tip of the day: Read the 'Raster' Manual
#Read chapters 5, 9 and 10 of vignette('Raster') carefully. You'll find a lot of useful commands, such as
#
#- aggregate, resample, crop, extend (!=extent), merge, trim, shift, flip, rotate
#- mask, cover
#- reclassify
#- focal
#- distance, gridDistance, pointDistance,distanceFromPoints
#- clump, edge, area
#- rasterize, rasterToPoints, rasterToPolygons
#- zonal, freq, crosstab
#

## Accessing cell values of a raster I
#Chapters 9 and 10 of vignette('Raster') example
#- getValues, extract, ...
plot(Bor_mat_cr)
cells <- cellFromRowCol(Bor_mat_cr, 5, 1:3)
cells #gives cell ID number!!!
extract(Bor_mat_cr, cells) ## gives cell values!
points(coordinates(Bor_mat_cr)[cells,], col='blue')

## Accessing cell values of a raster II
xy=xyFromCell(Bor_mat_cr, cells)#=coordinates(Bor_mat_cr)[cells,]
extract(Bor_mat_cr, xy) # = extract(Bor_mat_cr, cells)

# Change values e.g. for adding corridor
# take care! irrevesible! better work on a copy!
copy_ras <- Bor_mat_cr
copy_ras[cells] <- 250
plot(copy_ras)

## Compute distance to points
# calculate distance
my_dist <- distanceFromPoints(Bor_mat_cr,xy)
plot(my_dist) #units?

   
### many points distributed
cells1 <- c(cells, 250, 360)
xy=xyFromCell(Bor_mat_cr, cells1)#=coordinates(Bor_mat_cr)[cells,]
my_dist <- distanceFromPoints(Bor_mat_cr,xy)
plot(my_dist) #units?
   
   
##### extract rastervalues from all rasters at once
pred_dat <-   extract(predictors, xy) 
test <- data.frame(xy,pred_dat)
   
################################################################################
## Remote Sensing Data (example from JN)
################################################################################

#- load Landsat Stack (for NDVI calculation)
#-  example for multi-layer dataset
#
#scene:   LT51170562009223
#part of Global Land Survey 2010 dataset
#Landsat 5, Image taken 2009, day 223 (11 August)
#WRS path 117, row 56
#
ls <- raster(paste(raster_wd,"Landsat_GLS2010_example.tif",sep='/'))
ls
nlayers(ls)    # only 1 band out of 6
plot(ls, col= rainbow(200))

ls <- brick(paste(raster_wd,"Landsat_GLS2010_example.tif",sep='/'))
      # stack() also works
ls
nlayers(ls)    # 6 out of 6
plot(ls, col= rainbow(200))

### Plot Landsat image

#band combination 321 (Red, Green, Blue: true colour)
plotRGB(ls, r = 3, g = 2, b = 1,stretch = "lin",maxpixels = 1e6)

#band combination 352 (water)
plotRGB(ls, r = 3, g = 5, b = 2,stretch = "lin",maxpixels = 1e6)
   
#band combination 432 (NIR, Red, Green) - vegetation
plotRGB(ls, r = 4, g = 3, b = 2,stretch = "lin",maxpixels = 1e6)

#band combination 743 (IR, NIR, Red) to increase contrast 
plotRGB(ls, r = 6, g = 4, b = 3,stretch = "lin",maxpixels = 1e6)

   # Note: band 6 of original dataset was removed
#       (it is a Thermal IR band with 60m resolution)
# Band designations are sensor-specific (e.g. Landsat7 != Landsat 8)
# more information on band combinations
#      http://web.pdx.edu/~emch/ip1/bandcombinations.html

### Calculate NDVI
ndvi <- (ls[[4]] - ls[[3]]) / (ls[[4]] + ls[[3]])
range(values(ndvi), na.rm = T)

# set values < 0 to NA
ndvi2 <- ndvi
values(ndvi2)[which(values(ndvi2) <= 0)] <- NA
plot(ndvi2)

################################################################################
# Working with shapefiles
################################################################################
## Import shapefiles (polygons and lines)
#- dsn sets wd! no need to setwd().
Admin_shp <- readOGR(dsn=shapes_wd, layer="borneo_admin",
              stringsAsFactors=FALSE)[,c(1:3,5,7,17,18)]
PA_shp <-  readOGR(dsn=shapes_wd, layer="Bor_PA",
                   stringsAsFactors=FALSE)[,c(1:4)]
River_shp <- readOGR(dsn=shapes_wd, layer="sn_100000", 
                     stringsAsFactors=FALSE) ## sf update 2019 - hier ist kein subset 

test2 <- readShapePoly(fn=paste(shapes_wd,"borneo_admin.shp", sep='/'))[,c(1:3,5,7,17,18)]

## sf update 2019 - start

## Directly read as sf object - quicker for large files
Admin_sf <- read_sf(dsn = shapes_wd, 
                    layer = "borneo_admin",
                    stringsAsFactors = FALSE)[,c(1:3,5,7,17,18)]
## transformations
Admin_sf <- as(Admin_shp, "sf") ## from sp object ro sf object 
Admin_sp <- as(Admin_sf, "Spatial") ## from sf to sp object

PA_sf <- as(PA_shp, "sf") ## from sp object ro sf object 
River_sf <- as(River_shp, "sf") ## from sp object ro sf object 
test2_sf <- as(test2, "sf") ## from sp object ro sf object 

## sf update 2019 - end


### Retrieving general information of spatial layer
#Please note the similarity to accessing info from rasters.
PA_shp
str(extent(PA_shp))
extent(PA_shp)@xmin

## sf update 2019 - start
# for sf
Admin_sf
st_bbox(Admin_sf)
st_bbox(Admin_sf)$xmin
## sf update 2019 - end


# accessing d.f./ a.t.   
#PA_shp@data 
names(PA_shp) #returns column names of a.t.
names(PA_sf) ## sf update 2019 


### Retrieving summary for each column of attribute table [a.t.]
summary(PA_shp) # str(PA_shp) # PA_shp
print(PA_sf, n=1) ## sf update 2019 
summary(PA_sf) ## sf update 2019 
# not run: attributes(PA_shp)

### Retrieving information of a.t. I
head(PA_shp)   # tail(PA_shp)
head(PA_sf)   # tail(PA_sf) ## sf update 2019


### Retrieving information of a.t. II
PA_shp[1,] #returns first entry (row) of all 4 columns
PA_sf[1,] ## sf update 2019


### Retrieving information of a.t. III
PA_shp[,2] #returns summary of col only!
PA_sf[,2] # returns the second column of the data  ## sf update 2019 

PA_shp$NAME_ENG #returns a vector of comlumn elements # head(PA_shp$NAME_ENG)
# [1] "Kinabalu"
# [2] "Mulu"
# [3] "Niah"

PA_sf$NAME_ENG #returns a vector of comlumn elements # head(PA_shp$NAME_ENG) ## sf update 2019 

names(PA_shp)
PA_shp[1,1] # = PA_shp$SITE_ID[1]
PA_shp[2,3] # = PA_shp$COUNTRY[2]

names(PA_sf) ## sf update 2019
PA_sf[1,1] # = PA_shp$SITE_ID[1] ## sf update 2019
PA_sf[2,3] # = PA_shp$COUNTRY[2] ## sf update 2019


which(PA_shp$COUNTRY == 'Malaysia')
which(PA_sf$COUNTRY == 'Malaysia') ## sf update 2019

### Plotting the shapefile
spplot(PA_shp[,1]) #just plot first col with Site_ID
plot(PA_sf[,1]) #just plot first col with Site_ID ## sf update 2019

### Plot a detail on the map
#Another (the usual) plot command:
plot( Admin_shp, border = "deepskyblue4")
plot( PA_shp[,1], add = T, col="thistle2")
plot( PA_shp[1,], add = T, col="red")

## sf update 2019 - start

# base plot
plot( Admin_sf[,2], col = "transparent", border = "deepskyblue4")
plot( PA_sf[,1], add = T, col="thistle2", border = "black") ## be careful some few polygones are missing in sf object! 
plot( PA_sf[1,], add = T, col="red")

# tmap version of plot
tmap_mode("plot")
tm_shape(Admin_sf) + tm_borders( ) +
tm_shape(PA_sf[,1]) + tm_polygons( ) +
tm_shape(PA_sf[1,]) + tm_polygons( col = "red")

## sf update 2019 - end

## Simplify your life! Work with subsets
### Select Malaysia, plot and save as shapefile
Mal_PA_shp <- subset(PA_shp, PA_shp$COUNTRY == 'Malaysia')
plot( Admin_shp, border="deepskyblue4")
plot(Mal_PA_shp, add = T, col="green")
writeOGR(obj=Mal_PA_shp, dsn=output_wd, layer='test',
         driver ='ESRI Shapefile',overwrite=TRUE)

## sf update 2019 - start
# if not worling install dplyr package: install.packages(dplyr) ; library(dplyr)
Mal_PA_sf <- dplyr::filter(PA_sf, COUNTRY == 'Malaysia')
plot( Admin_sf[,1], col = NA, border = "deepskyblue4")
plot(Mal_PA_sf[,1], add = T, col="green")

write_sf(obj = Mal_PA_sf, 
         dsn = paste0(output_wd, "/Mal_PA_sf.shp"), 
         layer='test',
         driver ='ESRI Shapefile', 
         delete_layer = TRUE)
## sf update 2019 - end



### Save shapefiles in one go! [Write a loop]  
## no sf version. Loops work the same way in sf

ze_country <- unique(PA_shp$COUNTRY)
for (i in 1: length(ze_country))
{
  writeOGR(obj=PA_shp[PA_shp$COUNTRY == ze_country[i],],
  dsn=output_wd, layer=paste('PA_',toupper(substr(
  ze_country[i],1,3)),sep=''),driver ='ESRI Shapefile',overwrite=TRUE)
}

### Plot the four subsets in one go
image(Bor_mat)
for (i in 1: length(ze_country))
{ #plot.new()
  plot(PA_shp[PA_shp$COUNTRY==ze_country[i],],
       col= i+1,add=T)
#dev.off()
}

## Geospatial calculation and projection I
### command *gArea*
gArea(Admin_shp)# what does the warning message mean?
sum(st_area(Admin_sf))  ## in this case sf does not check projections but trasforms automatically. However, you should specify the projection. 
## sf update 2019

## Geospatial calculation and projection II
### Spatial transformation
Admin_shp_moll <-  spTransform(Admin_shp,
               CRS("+proj=moll +datum=WGS84"))

Admin_sf_moll <-  st_transform(Admin_sf, 4326) ## use EPSG codes to specify WGS84 lon/lat CRS ## sf update 2019
Admin_sf_moll <-  st_transform(Admin_sf,
                               "+proj=moll +datum=WGS84") ## sf update 2019

# plot(Admin_shp_moll) #do you see the difference?
gArea(Admin_shp_moll) # units?#gArea(Admin_shp_moll)/1e6
sum(st_area(Admin_sf_moll)) ## sf update 2019

### Example: Area of Malaysia in Borneo
Mal_Admin_shp <- subset(Admin_shp_moll,
                         Admin_shp_moll$NAME_0=='Malaysia')
gArea(Mal_Admin_shp)/1000000

Mal_Admin_sf <- subset(Admin_sf_moll,
                        Admin_shp_moll$NAME_0=='Malaysia') ## sf update 2019
sum(st_area(Mal_Admin_sf))/1000000 ## sf update 2019



## Area of a polygon - add new column to a.t.
gArea(Admin_shp_moll[3,])/1e6 # for a single polygon
# gArea(Admin_shp_moll, byid=TRUE) #for all polygons
plot(Admin_shp_moll)
plot(Admin_shp_moll[3,],col='red',add=T) 

st_area(Admin_sf_moll[3,])/1e6 # for a single polygon ## sf update 2019
# st_area(Admin_sf_moll)/1e6 #for all polygons ## sf update 2019



area_km2 <- gArea(Admin_shp_moll, byid=TRUE)/1e6
area_km2_sf <- st_area(Admin_sf_moll)/1e6 ## sf update 2019
head(area_km2)

Admin_shp_moll$my_area_km2 <- area_km2
Admin_sf_moll$my_area_km2 <- area_km2_sf ## sf update 2019
head(Admin_shp_moll)

## Extract elevation per Malaysian PA and save to a.t.
tmp <- extract(x = Bor_dem, y = Mal_PA_sf)   # gives a table ## same with sf
str(tmp)
lapply(X = tmp, FUN = summary)
lapply(X = tmp, FUN = mean)
   
Mal_PA_shp$mean_elevation <- round(unlist(lapply(tmp, FUN = mean)))
   plot(Admin_shp,col='transparent',borders='black',add=T)
View(Mal_PA_shp@data)

### Spatial queries
#### Vector data (Spatial Objects) (library 'rgeos')
#- gBuffer
#- gUnion
#- gIntersection
#- ...
#### library sf
#- st_buffer
#- st_union
#- st_intersection
#- ...
#
#### Raster data
#- buffer
#- ...and see before 'Raster' Manual
#


################################################################################
## Spatial Point Data
################################################################################
### Read spatial points - shp and csv/ txt
# setwd(recs_wd) #better include wd in the file name:
sp_shp <- readShapePoints(paste(shapes_wd,"FCsim.shp", sep='/'))
sp_shp <- rgdal::readOGR(paste(shapes_wd,"FCsim.shp", sep='/')) ## sf update 2019
sp_sf <- sf::st_read(paste(shapes_wd,"FCsim.shp", sep='/')) ## sf update 2019


pt_file <- paste(shapes_wd,'MyNewSpecies.csv',sep='/')
sp_csv <- read.table(pt_file, header=TRUE, sep=',')
class(sp_csv)
head(sp_csv)
class(sp_shp)
head(sp_shp)
class(sp_sf) ## - advantage is a data frame!
head(sp_sf)

### Plot the points (data frame!)
plot(Admin_shp,col='grey',border='white')
points(sp_csv$long,sp_csv$lat,cex=0.5,pch=15)
plot(sp_shp,col='blue',add=T) # easy handling of shapefiles

### Create a spatial point object
#Check package 'sp', e.g. vignette('over')
coordinates(sp_csv) <- ~long + lat
crs(sp_csv)  <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
class(sp_csv) # plot(sp_csv) #plotting command is simple
head(sp_csv) # this is @data slot; str(sp_csv) for more!

## create spatial object using sf
sp_sf <- st_as_sf(sp_csv,  
                  coords = c('long', 'lat'), 
                  crs = 4326) ## sf update 2019

## Spatial overlay
# retrieve the geometry (location) indices of PA_shp at
# the locations of sp_csv: which points are in PA_shp
insidePA <- sp_csv[PA_shp,] #returns SpatialPointsDataFrame
plot(insidePA, add=T, col='red',pch=17)
head(insidePA@data)

insidePA.sf <- sp_sf[PA_sf,] #returns sf object and warning that CRS 4326  may be dificult for intersections ## sf update 2019
plot(insidePA.sf, add=T, col='red',pch=17) ## sf update 2019
head(insidePA.sf)## sf update 2019


# with 'over' the information of the PA_shp can be added
insidePAdf <- over(insidePA, PA_shp)# returns data.frame
insidePAdf   
insidePA@data <- data.frame(insidePA@data,insidePAdf) #to add > 1 column


insidePAsf <- st_intersection(insidePA.sf, PA_sf)# returns sf and data.frame object  ## sf update 2019
insidePAsf    ## sf update 2019

# Repetition: for a raster: extract mean ann. temp. from Bor_mat
# and add it to a.t.
sp_csv$mean_t <- extract(Bor_mat,sp_csv ) ## same for sf  ## sf update 2019
mean(sp_csv$mean_t) #hist(sp_csv$mean_t,breaks=4)

################################################################################
# Least Cost Path LCP for connectivity 
################################################################################
## - atm no sf version here added  ## sf update 2019

library(gdistance)
start <- gCentroid(PA_shp[1,1], byid=FALSE, id = NULL)
end <- gCentroid(PA_shp[2,1], byid=FALSE, id = NULL)
# Several commands in one line when separated by ';':
plot(Bor_mat); points(start); points(end)
extent(end)

### Clip raster to save computation time
cr_extent <- c(113,117,3.5,6.5)
me_cr <- crop(Bor_mat,cr_extent)
plot(me_cr); points(start,pch=15); 
points(end,pch=15)

## Necessary calculations
### Calculate the transition layer
trans <- transition((1/me_cr), transitionFunction=mean,
                    directions=4, symm=F)
trans1 <- geoCorrection(trans)

### Calculate the shortest weighted connection
sPath <- shortestPath(trans1, start, end,output="SpatialLines")
### Calculate the length of the path
costDistance(trans, start, end) #units?

## Make a plot
plot(1/me_cr); points(start,pch=15); points(end,pch=15)
lines(sPath,lwd=3); plot(PA_shp, border='blue', lwd=2,add=T)

### A last cool interactive command:
zoom(me_cr) # click twice on the plot to define zoom extent
lines(sPath,lwd=3); plot(PA_shp, border='blue', lwd=2,add=T)

### Which connections intersect the PAs?
PA_inter1 <- gIntersection(spgeom1=PA_shp,spgeom2=sPath,
                          byid=TRUE)
plot(1/me_cr)
plot(sPath,lwd=3,add=T); plot(PA_shp,border='blue',lwd=2,add=T)
lines(PA_inter1,col='red',lwd=4); box()

identicalCRS(PA_shp,sPath)

### Save new path as shapefile
crs(sPath) <- CRS('+proj=longlat +datum=WGS84')
class(sPath)
cp_line <- SpatialLinesDataFrame(sl=sPath,
          data = data.frame(name = c(1:length(sPath@lines))))

writeOGR(obj=cp_line,dsn=output_wd,layer='costpath_line',
         driver = 'ESRI Shapefile', overwrite = TRUE)

### Final end of session plot
image(Bor_hs,col=grey(0:100/100))
plot(sp_csv,cex=0.5,pch=15, add=T)
plot(insidePA, col= 'red',add=T)
plot(PA_shp, border='white',add=T)
plot(River_shp, col="dodgerblue3", add=T)
text(114, 0, 'END OF SESSION', cex=3, col= 'grey')

writeOGR(obj=insidePA,dsn=output_wd, layer='inPA',
         driver ='ESRI Shapefile',overwrite=TRUE)

################################################################################
# Additional thingies
################################################################################

## Other plot commands
### Sp-plot
spplot(Bor_mat)

### Levelplot (package 'rasterVis')
levelplot(Bor_mat, FUN.margin=max)

### Densityplot
densityplot(Bor_dem)

### Use colors
#1. general information about colors in R
#   *http://www.stat.tamu.edu/~jkim/Rcolorstyle.pdf*
#
#2. names of defined colors in R
#   *http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf*
#
#3. RColorBrewer- package provides palettes for drawing nice maps
#   *http://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf*
#
#4. If possible use the viridis colour scheme (R library Viridis) ## sf update 2019



### Plot raster and vector data with colors
image(Bor_hs,col=grey.colors(10))
plot(Admin_shp,border='white',
     col=rgb(red=Admin_shp$ID_0/max(Admin_shp$ID_0),
             green=0,blue=1,alpha=0.25),add=T)

### Hexadecimal, RGB, color names
image(Bor_hs,col=grey.colors(50))
plot(Admin_shp[Admin_shp$NAME_0 == "Malaysia",],
        border='white',col="#6699CC",add=T)
plot(Admin_shp[Admin_shp$NAME_0 == "Indonesia",],
        border='white',col=rgb(0,1,0,0.5),add=T)
plot(Admin_shp[Admin_shp$NAME_0 == "Brunei",],
        border='white',col="red",add=T)

### Predefine colors in a vector...
ze_country <- unique(Admin_shp$NAME_0)
ze_country
mycol <- c('darkseagreen3','lightsteelblue1','red')
### ...and assign them [Write a loop]
image(Bor_hs,col=grey.colors(25))
for(i in 1: (length(ze_country))){
  plot(Admin_shp[Admin_shp$NAME_0 == ze_country[i],],
       col= mycol[i],add=T)
}

## Pimp up your map

## sf update 2019 - start 
## tmap examples: 
##  - https://geocompr.robinlovelace.net/adv-map.html
##  - https://cran.r-project.org/web/packages/tmap/vignettes/tmap-getstarted.html
##  - https://rpubs.com/chrisbrunsdon/UQ_tmap
## sf update 2019 - end

### Gridlines
image(Bor_mat, xlab="longitude", ylab="latitude")
plot(Admin_shp,border='gray',add=T)
grd <- gridlines(Admin_shp,
                    easts = pretty(bbox(Admin_shp)[1,]),
                    norths = pretty(bbox(Admin_shp)[2,]),
                    ndiscr = 20)
plot(grd, add=T, col='pink', lty=2)

###North Arrow and Scale Bar
image(Bor_mat, xlab="longitude", ylab="latitude")
plot(Admin_shp,border='gray',add=T)
north.arrow(110, 4.5, 0.25, cex.lab=0.7, col="grey80")
map.scale(118.4, -2, 1.8,"Km",2, 100, sfcol='red')

#########################   THE END ############################################
