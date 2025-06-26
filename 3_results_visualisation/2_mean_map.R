#########################################################################
## Creates density maps from outputs of the MINDZ model                ##
##                                                                     ##
#########################################################################

##
#-- initialisation : 

library(stringr  )   # manipulate characters to extract date from file name
library(lubridate)   # to deal with date formats
library(ggplot2  )
library(dplyr    )
library(terra    )   # to deal with rasters, for the bathymetry map
library(tidyterra)   # to plot the bathymetry map
library(sf       )   # to read whale sightings data

setwd(".../codes_article/1_depth_distribution")

load("./other_data/all_sightings.Rdata") # load sightings
world <- map_data("world")               # the map borders



#### plot maps ####

# specify copepod class (LH0, LF0, LF1, YF1) :
type <- "YF1"                       

##
#-- find all map files : 
all_files <- list.files(path = paste0("./other_data/fingrid_", type)) # find all maps of the year and copepod type

#j = "fingrid_LH0_2017_07_19_traj.tif" #if we want a specific map
for (j in all_files) { # for each file, 
  
  ##
  #-- extract the sightings of the corresponding dates : 
  date_temp <- str_extract(j, "\\d{4}_\\d{2}_\\d{2}")      # extract date
  date_temp <- ymd(str_replace_all(date_temp, "_", "-"))   # convert into date format
  
  sightings_temp <- all_sightings |>
    filter(Date == date_temp)           # select sightings inside these dates
  
  if (nrow(sightings_temp) == 0){       # if there is no sighting, look at a next map
    print(paste0(" !! no sightings found at at date : ", date_temp, " !!"))
    next
  }                     
  
  ##
  #-- Read map file : 
  agreg_map = rast(  paste0("./other_data/fingrid_", type, "/", j) )  # open the file
  
  # set the max value of the map depending on copepod class :
  if(type == "LH0"){themax <- 50} else{themax <- 35}      
  agreg_map[agreg_map > themax] <- themax
  
  #agreg_map <- aggregate(agreg_map, fact=2, fun=mean) # if we want to decrease resolution
  
  # set the color palette :
  cols3 <- c("#263B50", "#577590", "#4d908e", "#43aa8b", "#90be6d", "#f9c74f", "#f8961e", "#f3722c", "#f94144")
  
  ##
  #-- plot the map : 
  the_plot <- ggplot( ) +
    # aggregation values :
    geom_spatraster(data = agreg_map) +
    scale_fill_gradientn(colors = cols3) +
    labs(fill = "aggreg. value") +
    
    # coast limit and countries borders :
    geom_map(
      data = world, map = world,
      aes(map_id = region),
      color = "black", fill = "grey60", linewidth = 0.1) +
    
    # presence and absence data :
    geom_point(
      data = sightings_temp |> filter(obs == "surv", Sp_code == 0),
      aes(Long, Lat),
      size = 0.3, stroke = 0.1, alpha = 0.6, color = "grey90" #grey30
    ) +
    geom_point(
      data = sightings_temp |> filter(obs == "surv", Sp_code == 1),
      aes(Long, Lat),
      size = 1.4, alpha = 1, color = "black"
    ) +
    geom_point(
      data = sightings_temp |> filter(obs == "oport"),
      aes(Long, Lat),
      size = 1.4, alpha = 1, color = "black" #"ED7676"
    ) +
    
    # limits :
    coord_sf(xlim=c(-70.25, -56), ylim=c(41.9, 51.4)) +
    theme_void()
  
  #the_plot
  
  # Save the plot : 
  fact <- 1.1
  ggsave(
    paste0("./results/maps/", type, "/", type, "_", date_temp, ".png"),
    dpi    = 300,     # the resolution
    width  = 7 *fact,
    height = 5 *fact,
    bg     = "white"  # the background color
  )
}








#### -- variance/mean -- ####
# Mean and variance of all aggregation maps of each copepod class

# specify copepod class (LH0, LF0, LF1, YF1) :
type        <- "YF1"        
all_files   <- list.files(path = paste0("./other_data/fingrid_", type), full.names = TRUE) # find all maps of the year and copepod type
all_rasters <- rast(all_files)

# set the max value of the map depending on copepod class :
if(type == "LH0"){themax <- 50} else{themax <- 35}
all_rasters[all_rasters > themax] <- themax

# determine standard deviation and mean :
aggreg_sd    <- app(all_rasters, sd,   na.rm = TRUE)
aggreg_mean  <- app(all_rasters, mean, na.rm = TRUE)


# Export in tiff for correlations : 
# Used in script 1_correlations.R
# !! do not run the themax line before !!
writeRaster(aggreg_mean, paste0("./other_data/mean_aggreg_", type, ".tiff"))


##
#-- plot maps :

# color palette :
cols3 <- c("#263B50", "#577590", "#4d908e", "#43aa8b", "#90be6d", "#f9c74f", "#f8961e", "#f3722c", "#f94144")

the_plot <- ggplot( ) +
  geom_spatraster(data = aggreg_sd) +
  
  # density of particles :
  scale_fill_gradientn(colors = cols3) +
  labs(fill = "aggreg. value") +
  
  # coast limit and countries borders :
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "black", fill = "grey60", linewidth = 0.1) +
  
  # size :
  coord_sf(xlim=c(-70.25, -56), ylim=c(41.9, 51.4)) +
  
  theme_void()


#the_plot

# Save the plot : 
fact <- 1.5
ggsave(
  paste0("./results/maps/", type, "_SD_map.png"),
  dpi    = 300,     # the resolution
  width  = 7 *fact,
  height = 5 *fact, 
  bg     = "white"  # the background color
)



the_plot <- ggplot( ) +
  geom_spatraster(data = aggreg_mean) +
  
  # density of particles :
  scale_fill_gradientn(colors = cols3) +
  labs(fill = "aggreg. value") +
  
  # coast limit and countries borders :
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "black", fill = "grey60", linewidth = 0.1) +
  
  # size :
  coord_sf(xlim=c(-70.25, -56), ylim=c(41.9, 51.4)) +
  
  theme_void()


#the_plot

# Save the plot : 
fact <- 1.5
ggsave(
  paste0("./results/maps/", type, "_mean_map.png"),
  dpi    = 300,     # the resolution
  width  = 7 *fact,
  height = 5 *fact, 
  bg     = "white"  # the background color
)

