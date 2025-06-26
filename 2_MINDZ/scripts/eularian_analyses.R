#########################################################################
## Creates eularian analyses from the MINDZ outputs                    ##
## (works only with CIOPS inputs for now)                              ##
##                                                                     ##
## Andeol Bourgouin                                                    ##
## andeol.bourgouin.1@ulaval.ca                                        ##
#########################################################################


#### initialisation ####

library(dplyr)
library(ggplot2)
library(terra)        #to deal with rasters, for the bathymetry map
library(tidyterra)    #to plot the bathymetry map
library(ggnewscale)   #to put 2 scales in a ggplot
library(gganimate)    #animated plots


setwd("/Users/andy/Desktop/travail/these/analyses/git/MINDZ/scripts")

# load MINDZ output file modified in another R code :
load("traj_part.Rdata") #file with the position of each particle at each time




#### spatial extent of the analyses ####

# choose resolution of the eularian analyses : 
resolution <- 0.1

# set the boudaries of the analyses :
# from the extreme positions of the particles, 
# (a buffer can be set to be sure that the particles in the border
# are taken into account)
min_lon <- min(traj.part$lon) #-resolution
max_lon <- max(traj.part$lon) #+resolution
min_lat <- min(traj.part$lat) #-resolution
max_lat <- max(traj.part$lat) #+resolution





#### Select particles ####

#-- Select only certain times
unique(traj.part$time) # (the time here is the time step)
particles <- traj.part[which(traj.part$time >= 24),]


#-- Remove particles at bottom depth deeper than 600m :
# Because we don't want to see the agregations over the ocean basin.
# to do so, 
# extract bottom depth from the bathymetry file
kartras <- rast("gebco_grand.tif") 
coords  <- cbind(particles$lon, particles$lat)
values  <- extract(kartras, coords)

# add it to the 'particles' dataset
particles$bdep <- values

#and keep the particles above the choosen depth
particles <- particles[which(particles$bdep > -600),]


#-- Remove particles too close to the border :
# those particles interact with the border and therefore do wrong patterns
#we don't do it for min latitude, as the ocean basin already removed the border
particles <- particles[-which(particles$lon < min(traj.part$lon)+0.5),]
particles <- particles[-which(particles$lon > max(traj.part$lon)-1.5),]
particles <- particles[-which(particles$lat > max(traj.part$lat)-1.5),]




#### set the grid ####

# Create a sequence of latitudes and longitudes to set the grid :
latitudes  <- seq(min_lat, max_lat, by = resolution)
longitudes <- seq(min_lon, max_lon, by = resolution)

# Create the grid :
grid <- expand.grid(latitude = latitudes, longitude = longitudes)


#-- Count the nb of parts in each grid cell, for each selected time step :
for (i in min(particles$time): max(particles$time)) {
  print(i) # To see where we are
  
  # Select the particles of the time step :
  parts <- particles[which(particles$time == i),]
  
  # Create a new column and count the nb of particles for each cell :
  grid <- grid %>%
    mutate(tme = sapply(1:nrow(grid), function(i) { #for each cell,
      sum(parts$lon >= grid[i, "longitude"] & parts$lon < grid[i, "longitude"] + resolution &
          parts$lat >= grid[i, "latitude"]  & parts$lat < grid[i, "latitude"]  + resolution)
    }))
  
  # Set the name of the new column :
  colnames(grid)[length(colnames(grid))] = paste0('time_', i)
}


#max density : 
max(grid$time_24)


#-- Mean night - day :
# de 10h a 14h
# de 20h a 24h


# Save the gridded data : 
save(grid, file = paste0("grid_3days.Rdata") )


#### plots ####

#-- Just a simple one at a time step : 
ggplot(data = grid, aes(x = longitude, y = latitude, fill = time_24)) +
  geom_tile(aes(alpha = time_24)) +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Nombre de particles par cellule de grid") +
  theme_minimal()


#-- Reshape bathymetry data :

# Remove positive altitude data
kartras[kartras > 0.99] <- NA

# Fix bathym of the ocean basin (to have more contrast in the shelf)
kartras[kartras <= -600] <- -650


# There is a mismatch with the projetction that needs to be fixed...
# kartras2 <- kartras
# kartras <- kartras2
# crsuggest::suggest_crs(kartras)
# dest_proj <- "EPSG:32183"
# kartras <- project(kartras2, dest_proj)

# Download cartography for the ground :
world <- map_data("world")



#-- The big plot at a time step :

ggplot()           +
  # bathymetry :
  geom_spatraster(data = kartras)+
  scale_fill_gradientn(colors = gray.colors(200, #nb de couleurs
                                            start = 0.03, #valeur la + sombre
                                            end = .92)) + #valeur la + claire
  
  # coast limit and countries borders :
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "#a18f7c", size = 0.1) +
  
  # density of particles : 
  new_scale_fill() + #because we use a second scale after bathymetry
  geom_tile(data = grid, aes(x = longitude,   y = latitude, 
                             fill = time_24, alpha = time_24)) +
  scale_fill_gradient(low = "white", high = "red") +
  
  # size :
  coord_sf(xlim=c(min_lon, max_lon), ylim=c(min_lat, max_lat)) +
  
  # look of the plot :
  theme_classic() +
  theme(
    axis.line   = element_line(size = 0.3, colour = "black"),
    axis.text.y = element_text(size=5),
    axis.text.x = element_text(size=5), 
    plot.title  = element_text(size=9) 
    #legend.position = "none" 
  )                 



#-- The same big plot animated for all time steps :

#d'abord il faut changer les données pour avoir un dataframe avec dune colonne 
#density, et une colonne time step 
library(tidyr)

# Gather columns into key-value pairs
grid_animate <- grid |>
  gather(key = "time", value = "density", -c(latitude, longitude)) |>
  filter(density > 0)


p <- ggplot(data = grid_animate)           +
  # bathymetry :
  geom_spatraster(data = kartras)+
  scale_fill_gradientn(colors = gray.colors(200, #nb de couleurs
                                            start = 0.03, #valeur la + sombre
                                            end = .92),  #valeur la + claire
                       guide = guide_none()) +  #no legend
  
  # coast limit and countries borders :
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "#a18f7c", size = 0.1) +
  
  # density of particles : 
  new_scale_fill() + #because we use a second scale after bathymetry
  geom_tile(aes(x = longitude,   y = latitude, alpha = density,
                             fill = density)) +
  guides(alpha = "none") + #remove alpha legend
  scale_fill_gradient('density', low = "white", high = "red") +
  
  # size :
  coord_sf(xlim=c(min_lon, max_lon), ylim=c(min_lat, max_lat)) +
  
  # animate :
  transition_states(time) +
  
  # look of the plot :
  labs( title = "Time : {closest_state}",
        x     = "",
        y     = "" ) +
  #guides(fill = guide_legend(title = "density")) +
  theme_classic() +
  theme(
    axis.line   = element_line(size = 0.3, colour = "black"),
    axis.text.y = element_text(size=5),
    axis.text.x = element_text(size=5), 
    plot.title  = element_text(size=9),
    legend.key.size = unit(0.1, 'cm'), 
    legend.text = element_text(size=3), 
    legend.title = element_text(size=5)
  )                 

gganimate :: animate( p,
                      height = 8,
                      width  = 10,
                      units  = "cm",
                      res    = 600,
                      renderer = gifski_renderer("density_evolution.gif") )


 
 #What are the next moves ?
 #- Save the new created data
 #- To mean the data of half a day
 #- To plot that mean
 #- To parallelise the count of particles