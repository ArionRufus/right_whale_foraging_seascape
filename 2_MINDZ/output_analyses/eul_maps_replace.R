#########################################################################
## Creates density maps from outputs of eularian_grids.R               ##
##                                                                     ##
## Andeol Bourgouin                                                    ##
## andeol.bourgouin.1@ulaval.ca                                        ##
#########################################################################


#### Initialisation ####
setwd("~/projects/def-frmap5/MINDZ_andeol")
library(renv)
renv::load()
renv::activate()

library(dplyr)       # for data manipulation
library(terra)       # to deal with rasters, for the bathymetry map
#library(ggnewscale)  # to have two scales on the same plot
library(ggplot2)
library(tidyterra)   # to plot the bathymetry map
library(sf)          # to read whale sightings data


#setwd("~/Desktop/travail/these/analyses/git/MINDZ/output_analyses")
setwd("~/projects/def-frmap5/MINDZ_andeol/MINDZ/output_analyses")



#### Prepare data ####

#-- Get all file names :

folder_path <- "./results/to_compute" # location of the input files
all_files <- list.files(folder_path)
#all_files = all_files[-c(1:5)] # if we want to select just specific files



#-- Get data for the map background : 

world <- map_data("world") #the borders




#### Create the maps ####


bwa = 3 #to debug
for (bwa in 1 : length(all_files)) { #for each file, 
  
  # Load the file  (grid_test):
  current_file = all_files[bwa] 
  agreg_map = rast(  paste0(folder_path, '/', current_file) )
  themax <- as.integer(global(agreg_map, quantile, probs=0.9999, na.rm=TRUE))
  agreg_map[agreg_map > themax] <- themax
  # aThe plot :
  
  the_plot <- ggplot( ) +
    geom_spatraster(data = agreg_map) +
    
    # density of particles :
    # geom_tile(data = grid_test, aes(x = longitude,   y = latitude,
    #                                 fill = proj, color = proj), color = NA) +
    scale_fill_viridis_c(option = 'G', begin = 0, na.value = 'gray20') +
    
    # coast limit and countries borders :
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "black", fill = "grey70", size = 0.1) +

    
    # size :
    coord_sf(xlim=c(-70.25, -53), ylim=c(41.5, 51.4)) +
    
    # look of the plot :
    theme_classic() +
    theme(
      axis.line   = element_line(size = 0.3, colour = "black"),
      axis.text.y = element_text(size=6),
      axis.text.x = element_text(size=6),
      plot.title  = element_text(size=9)
    )
  
  the_plot
  
  # Save the plot : 
   ggsave(
    paste0("map_", current_file, ".png"),
    dpi    = 100,     # the resolution
    width  = 10,
    height = 10  )
   
  }





