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
library(parallel)   #to make parallel computations
library(doParallel) #to make parallel computations
library(tictoc)

#setwd("/Users/andy/Desktop/travail/these/analyses/git/MINDZ/output_analyses")
setwd("~/Disque_20TO_DATA/MINDZ/output_analyses")

# Model characteristics :
# used for the name of inputs and outputs
part         = 'LH'  # particle class : LH=late_hyperboreus
init_distrib = 'all' # initial distribution


# load MINDZ output file modified in the "main_figures" R code :
# file with the position of each particle at each time
load(paste0("traj_part", part, init_distrib, ".Rdata")) 




#### spatial extent of the analyses ####

# choose resolution of the eularian analyses : 
resolution <- 0.05

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
particles <- traj.part[which(traj.part$time > 24 &
                               traj.part$time < 121),] #121
rm(traj.part)

#-- Remove particles at bottom depth deeper than 600m :
# Because we don't want to see the agregations over the ocean basin.
# to do so, 
# extract bottom depth from the bathymetry file
kartras <- rast("gebco_grand.tif") 
coords  <- cbind(particles$lon, particles$lat)
values  <- terra::extract(kartras, coords)
rm(coords)

# add it to the 'particles' dataset
particles$bdep <- unlist(as.vector(values))
rm(values)

#and keep the particles above the choosen depth
particles <- particles[which(particles$bdep > -600),]



#-- Remove particles too close to the border :
# those particles interact with the border and therefore do wrong patterns
#we don't do it for min latitude, as the ocean basin already removed the border
particles <- particles[-which(particles$lon < min_lon+0.5),]
particles <- particles[-which(particles$lon > max_lon-1.5),]
particles <- particles[-which(particles$lat > max_lat-1.5),]

save(particles, file = paste0("particles", part, init_distrib, "_res0.5.Rdata") )
load(paste0("particles", part, init_distrib, "_res0.5.Rdata") )



#### set the grid ####

# Create a sequence of latitudes and longitudes to set the grid :
min_lat_grid <- round(min_lat*2, 1) / 2
max_lat_grid <- round(max_lat*2, 1) / 2 - 1.5
min_lon_grid <- round(min_lon*2, 1) / 2 + 0.5
max_lon_grid <- round(max_lon*2, 1) / 2 - 1.5


latitudes  <- seq(min_lat_grid, max_lat_grid, by = resolution)
longitudes <- seq(min_lon_grid, max_lon_grid, by = resolution)

# Create the grid :
grid <- expand.grid(latitude = latitudes, longitude = longitudes)
grid_dpt <- grid #this grid will stock the bottom depth values

#-- Count the nb of parts in each grid cell, for each selected time step :
tmax <- max(particles$time)



for (i in min(particles$time): tmax) {
  print(paste0(i, " / ", tmax)) # To see where we are
  
  # Select the particles of the time step :
  parts <- particles[which(particles$time == i),]
  
  # Create a new column and count the nb of particles for each cell :
  grid <- grid %>%
    mutate(tme = sapply(1:nrow(grid), function(i) { #for each cell,
      sum(parts$lon >= grid[i, "longitude"] & parts$lon < grid[i, "longitude"] + resolution &
          parts$lat >= grid[i, "latitude"]  & parts$lat < grid[i, "latitude"]  + resolution)
    }))
  
  grid$tme      <- 0
  #grid_dpt$bdep <- 0
  
  for (j in 1:nrow(grid)) {
    bop <- parts[ which(
      parts$lon >= grid[i, "longitude"] & parts$lon < grid[i, "longitude"] + resolution &
        parts$lat >= grid[i, "latitude"]  & parts$lat < grid[i, "latitude"]  + resolution), ]
    
    grid$tme[j]      <- nrow(bop)
    
    #grid_dpt$bdep[j] <- mean(bop$bdep)
  }
  
  # Set the name of the new column :
  colnames(grid)    [length(colnames(grid))]     = paste0('time_', i, '_res0.5')
  #colnames(grid_dpt)[length(colnames(grid_dpt))] = paste0('time_', i, '_res0.5')
}





#__________________________________________________________________________




#mode parallélisation : 
detectCores()
parCluster <- makeCluster( 10, methods = FALSE)
registerDoParallel( parCluster )

# Parallelized call to pricefit()
eularian.parallel <- function( i, particles, grid, resolution ){
  
  
  parts <- particles[which(particles$time == i),]
  parts$lon <- round(parts$lon*2, 1) / 2
  parts$lat <- round(parts$lat*2, 1) / 2
  
  # Create a new column and count the nb of particles for each cell :
  
    # tme = sapply(1:nrow(grid), function(k) { #for each cell,
    #   sum(parts$lon >= grid[k, "longitude"] & parts$lon < grid[k, "longitude"] + resolution &
    #         parts$lat >= grid[k, "latitude"]  & parts$lat < grid[k, "latitude"]  + resolution)
    # })
  
  
  #tme  <- rep(0, nrow(grid)) 
  dpt  <- rep(0, nrow(grid))
                                  
  for (j in 1:nrow(grid)) {
   bop <- which(
     near(grid[j, "longitude"], parts$lon, tol = 0.04) &
     near(grid[j, "latitude"] , parts$lat, tol = 0.04) )

   if (length(bop)>0){
     parts  <- parts[-bop,]
     bop    <- parts[bop,]
     #tme[j] <- nrow(bop)
     dpt[j] <- mean(bop$bdep, na.rm = T)
   }
  }
  
  out <- data.frame(
    #a = tme , 
    b = dpt
  )
  #!!!colnames(out) [1] = paste0('dens_time_', i, '_res0.5')
  colnames(out) [1] = paste0('dpt_time_', i, '_res0.5')
  return(out)
  
}

#ptime <- system.time({ # this is to check the time of execution
tic()
z.out <- foreach( i=min(particles$time): tmax , .packages = 'dplyr') %dopar% { 
  eularian.parallel( i, particles, grid, resolution )
}

stopCluster( parCluster )
toc()


grid_prof <- bind_cols(grid, z.out)
#grid_test <- bind_cols(grid, z.out)


head(z.out[[3]])


#__________________________________________________________________________







#max density : 
max(grid_test$time_25_res0.5, na.rm = T)
save(grid_test, file = paste0("grid_60days", part, init_distrib, "_res0.05.Rdata") )
load(paste0("grid_60days", part, init_distrib, "_res0.05.Rdata") )
grid <- grid_test
rm(grid_test)
#-- Mean night - day :
# de 10h a 14h
# de 20h a 24h


# Save the grided data : 
#save(grid, file = paste0("grid_10days", part, init_distrib, "_res0.5.Rdata") )
load( paste0("grid_10days", part, init_distrib, "_res0.5.Rdata") )


#### plots ####

#-- Just a simple one at a time step : 
ggplot(data = grid, aes(x = longitude, y = latitude, fill = time_120)) +
  geom_tile(aes(alpha = time_84)) +
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



#-- The big plot at a time step or total:
#mettre le temps en h, +2
grid$td30 <- NA
for (i in 1:nrow(grid)) {
  grid$td30[i] <- sum(grid[i, c(3:362)], na.rm = T)
}

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
                             fill = time_25_res0.5, alpha = time_25_res0.5)) +
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



#-- For the poster : 
# Déjà il faut que je gère le décalage entre les couches d'informations. 
# Ensuite il faut que je m'inspire de la carte de main figures
grid$proj = grid$td30 / max(grid$td30, na.rm = T)

# grid <- grid  %>%
#   filter( time_120 >= 1)

library(ggnewscale)

#add a layer to put over the basin : 
kart = rast("gebco_grand.tif") #downloaded before from gebco
kart[kart > 0.99] <- NA     #remove terrestrial data
kart[kart < -450 & kart > -550] <- -450
basin <- kart
kart[kart <= -550] <- NA
basin[basin >=  -450] <- NA
basin[basin <  -450] <- 1
basin <- as.factor(basin)

ggplot( data = grid, aes( x = longitude, y = latitude ) ) +
  coord_quickmap() +   
  
  # bathymetry :
  geom_spatraster(data = kart, maxcell = 5e+06)+
  scale_fill_gradientn(na.value = "#155176", 
                       colors = gray.colors(200, #nb de couleurs
                                            start = 0.03, #valeur la + sombre
                                            end = .92)) + #valeur la + claire

    # density of particles : 
  new_scale_fill() + #because we use a second scale after bathymetry
  # geom_tile(data = grid, aes(x = longitude,   y = latitude, 
  #                            fill = proj, alpha = proj)) +
  # scale_fill_gradient(low = "white", high = "red") +
  geom_tile(data = grid, aes(x = longitude,   y = latitude, 
                             fill = proj), color = NA) +
  scale_fill_viridis_c(option = 'G', begin = 0.1) +
  
  # coast limit and countries borders :
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "#B38A3E", size = 0.1) +
  
  
  # the basin
  new_scale_fill() +
  geom_spatraster(data = basin, maxcell = 5e+06) +
  scale_fill_manual(values = "#155176", na.value = NA) +
  
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
  #transition_states (time) +
  transition_manual(time) +
  # look of the plot :
  # labs( title = "Time : {frame_along}",
  #       x     = "",
  #       y     = "" ) +
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
                      renderer = gifski_renderer(
                        paste0("./results/density_evolution", part, init_distrib, ".gif")) )







#### look at depth distribution ####
particles <- traj.part[which(traj.part$time >= 50),]
coords  <- cbind(particles$lon, particles$lat)
values  <- terra::extract(kartras, coords)
particles$bdep <- values

traj_100 <- particles[which( round(particles$bdep) == 100 ),]

ggplot()+
  geom_histogram(data = traj_100, aes(x = dep), alpha = 1, 
                 bins = 50, fill = '#0A7999') +
  coord_flip() +
  scale_x_reverse() + 
  theme_minimal() +
  xlab("depth") + ylab("")
