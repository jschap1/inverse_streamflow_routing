# SWOT virtual gages
#
# 6/6/2023 JRS
# Creating synthetic SWOT observations for ISR
#
# INPUTS
# swot
# grwl
# bbox
# rivpoints
#
# OUTPUTS
# overpass_table
# crds_overpasses

library(rgdal)
library(sf)

setwd("/home/jschap/Documents/ISR/Ohio/")

# Assign river widths to each point on the ISR model river channel ------------

crds.sf <- st_read("/home/jschap/Documents/ISR/inverse_streamflow_routing/ohio_data/river_pixels.shp")
grwl <- st_read("/home/jschap/Documents/ISR/inverse_streamflow_routing/ohio_data/GRWL_summaryStats_V01.01/GRWL_summaryStats.shp")
bbox <- st_read("/home/jschap/Documents/ISR/inverse_streamflow_routing/ohio_data/boundingbox_for_cropping.shp")
grwl_crop <- st_crop(grwl, bbox)

# Find the nearest feature in riv for each point in crds
nearest_feature <- st_nearest_feature(crds.sf, grwl_crop)

# Join the width values from riv to crds based on the nearest feature
crds.sf <- st_join(crds.sf, grwl_crop, join = st_nearest_feature)

# Extract the width values to the "width" column in crds
crds.sf$width <- crds.sf$width_mean

# Extract longitude and latitude from the geometry column
lon <- st_coordinates(crds.sf)[,1]
lat <- st_coordinates(crds.sf)[,2]

plot(crds.sf["width"], col = crds.sf$width, 
     main = "Map of Widths", 
     pch = 19,
)

# Find SWOT overpasses ---------------------------------------------------

swot <- st_read("/home/jschap/Documents/ISR/inverse_streamflow_routing/ohio_data/swot_science_hr_Aug2021-v05_shapefile_swath/swot_science_hr_2.0s_4.0s_Aug2021-v5_swath.shp")

plot(crds.sf$geometry, pch = 19, cex = 0.5, main = "river channel and swot")
plot(swot$geometry, add  = TRUE)

crds.overpasses <- st_intersection(crds.sf, swot)

plot(crds.overpasses["START_TIME"], pch=19)

crds.overpasses$START_TIME[crds.overpasses$width_mean<100] <- NA

plot(crds.overpasses["START_TIME"], pch=19)

st_write(crds.overpasses, dsn = "./swot_overpasses_ohio.shp",
         driver = "ESRI Shapefile")

# There are more than n_riv_pixels features in crds.overpasses bc 
# SWOT observes some points more than once in one orbit cycle

# Create table of overpasses for use in MATLAB -------------------------

lon <- st_coordinates(crds.overpasses)[,1]
lat <- st_coordinates(crds.overpasses)[,2]
start_time <- crds.overpasses$START_TIME
orbit_day <- as.integer(substr(start_time, 5, 6))

T <- data.frame(lon=lon, lat=lat, orbit_day = orbit_day)
T1 <- na.omit(T) # remove "observations" where river is too narrow

write.table(T1, "./overpass_table.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE
            )


