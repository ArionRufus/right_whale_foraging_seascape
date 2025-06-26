#########################################################################
## Creates density maps from outputs of eularian_grids.R               ##
##                                                                     ##
## Andeol Bourgouin                                                    ##
## andeol.bourgouin.1@ulaval.ca                                        ##
#########################################################################


#### Initialisation ####

library(dplyr)       # for data manipulation
library(terra)       # to deal with rasters, for the bathymetry map
library(ggnewscale)  # to have two scales on the same plot
library(ggplot2)
library(tidyterra)   # to plot the bathymetry map
library(sf)          # to read whale sightings data
setwd("~/Desktop/travail/these/analyses/git/MINDZ/output_analyses")





#### Prepare data ####

#-- Get all file names :

folder_path <- "./results/2018/to_compute" # location of the input files
all_files <- list.files(folder_path)
#all_files = all_files[-c(1:5)] # if we want to select just specific files



#-- Get data for the map background : 

world <- map_data("world") #the borders
basin = rast("gebco_grand.tif") #bathymetry, downloaded before from gebco
# Just keep thethe ocean basin as a qualitative layer :
basin[basin >= -600] <- NA
basin[basin <  -600] <- 1
basin <- as.factor(basin)



#-- Get whale sightings :

# get opportunistic data : 
opport = readRDS("Sightings_8days_Opportunistic_Global_2018.RDS")
# translate coordinates from sf object to numbers :
ouca <- st_coordinates(opport$geometry) 
opport$longitude <- ouca[, 1]
opport$latitude  <- ouca[, 2]
# just keep important variables :
opport <- opport |>
  select(c(Date8days, longitude, latitude))

# Get survey data :
bibopeloula = readRDS("Sightings_8days_Surveys_Global_2018.RDS")
bibopeloula = bibopeloula |>
  filter(SPECCODE == 1) #just keep presence data

# Merge the 2 datasets :
all_whales <- data.frame(
  date = c(bibopeloula$Date8days, opport$Date8days), 
  lon  = c(bibopeloula$Long,      opport$longitude), 
  lat  = c(bibopeloula$Lat,       opport$latitude), 
  cat  = c(rep("sighting", nrow(bibopeloula)) , rep("opportunistic", nrow(opport)) )
)
all_whales$date = as.character(all_whales$date)



#-- There is a miss match of 1 day in the dates, so equi_dates is a variable
#that make the equivalence : 

equi_dates = c("2018-05-08")

# For LH :
# equi_dates = c("2018-05-09", "2018-05-17", "2018-05-25", "2018-06-02",
#                "2018-06-10", "2018-06-18", "2018-06-26", "2018-07-04", 
#                "2018-07-12", "2018-07-20", "2018-07-28", "2018-08-05", 
#                "2018-08-13", "2018-08-21", "2018-08-29")


# For YF :
# equi_dates = c("2018-04-07", "2018-04-15", "2018-04-23", "2018-05-01", 
#                "2018-05-09", "2018-05-17", "2018-05-25", "2018-06-02", 
#                "2018-06-10", "2018-06-18", "2018-06-26", "2018-07-04", 
#                "2018-07-12", "2018-07-20", "2018-07-28", "2018-08-05", 
#                "2018-08-13", "2018-08-21", "2018-08-29", "2018-09-06")






#### Create the maps ####


#bwa = 1 #to debug
for (bwa in 1 : length(all_files)) { #for each file, 
  
  # Load the file  (grid_test):
  current_file = all_files[bwa] 
  load( paste0(folder_path, '/', current_file) )
  
  
  # Sum the densities : 
  grid_test$tot <- NA
  for (i in 1:nrow(grid_test)) {
    grid_test$tot[i] <- sum(grid_test[i, -c(1:2)], na.rm = T)
  }
  
  #show histogram of densities :
  # bop = ggplot() +
  #   geom_histogram(data = grid_test, aes(x = tot)) +
  #   scale_y_continuous(trans='sqrt')
  # 
  # print(bop)
  
  # Apply a treshold to the densities :
  the.max   <- 1500
  grid_test$proj <- grid_test$tot
  grid_test$proj[which(grid_test$proj > the.max)] <- the.max
  #grid_test$proj <- grid_test$proj / max(grid_test$proj) #if we want realtive density
  
  # Find the sightings of the right date :
  whales <- all_whales |>
    filter(date == equi_dates[bwa])
  
  
  # aThe plot :

  the_plot <- ggplot( data = grid_test, aes( x = longitude, y = latitude ) ) +
    coord_quickmap() +

    # density of particles :
    geom_tile(data = grid_test, aes(x = longitude,   y = latitude,
                               fill = proj, color = proj), color = NA) +
    scale_fill_viridis_c(option = 'G', begin = 0, limits = c(0,the.max)) +

    # coast limit and countries borders :
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "black", fill = "#B38A3E", size = 0.1) +

    # the basin :
    new_scale_fill() +
    geom_spatraster(data = basin, maxcell = 5e+06) +
    scale_fill_manual(values = "white", na.value = NA) +
    guides(fill = "none") +
    
    # the whales :
    new_scale_color() +
    geom_point(data = whales, aes(x = lon, y = lat, color = cat),
               alpha = 0.3, size = 1.5) +
    scale_color_manual(values = c('green', 'red')) +

    # size :
    coord_sf(xlim=c(-70.1, -50), ylim=c(40.5, 50.9)) +

    # look of the plot :
    theme_classic() +
    theme(
      axis.line   = element_line(size = 0.3, colour = "black"),
      axis.text.y = element_text(size=6),
      axis.text.x = element_text(size=6),
      plot.title  = element_text(size=9)
    )
  
  # Save the plot : 
  ggsave(paste0("map_", current_file, ".pdf"),
         the_plot, width = 11, height = 11)
}





