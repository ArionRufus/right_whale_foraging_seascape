#########################################################################
## Creates eularian analyses from the MINDZ outputs                    ##
## Determine for each time step in a grid of a specific resolution the ##
## number of particles for each pixel of the grid.                     ##
## Inputs : netcdf files, each presenting position of particles        ##
## for several time steps.                                             ##
## Outputs : R files, each presenting densities for several time steps.##
## (works only with CIOPS inputs for now)                              ##
##                                                                     ##
## Andeol Bourgouin                                                    ##
## andeol.bourgouin.1@ulaval.ca                                        ##
#########################################################################



#### initialisation ####

library(ncdf4)       # to read netcdf files
library(dplyr)       # for data manipulation
library(terra)       # to deal with rasters, for the bathymetry map
library(tictoc)      # to know the computation time
library(ggnewscale)  # to have two scales on the same plot
library(parallel)    # to make parallel computations
library(doParallel)  # to make parallel computations
library(ggplot2)
library(tidyterra)   # to plot the bathymetry map


#setwd("~/Desktop/travail/these/analyses/git/MINDZ/output_analyses")
setwd("/Ext_20T_andromede/0_ARCTUS_Projects/15_SmartWhales/DATA/MINDZ/output_analyses")


# Model characteristics :
part         = 'LH'  # particle class : LH=late_hyperboreus
init_distrib = 'all' # initial distribution


# Set the path where the ncdf MINDZ output data are located
folder_path <- "../outputs/to_compute"

# Get a list of files in the folder
all_files <- list.files(folder_path)


# Loop for all files : 

tic() # Optional, to have an idea of the time taken
for (current_file in all_files) {
  
  # Print the actual file analysed :
  print( paste('TIME : ', which(all_files == current_file), ' / ', length(all_files), ', : ', current_file))
  
  
  
  #### extract data ####
  
  # open ncdf output file
  ncid.out <-  nc_open( paste0(folder_path, '/', current_file) )
  
  
  # Define the first and last time steps of the simulation 
  # to do the extraction
  time.1 <- 1 
  time.n <- ncid.out$dim$time$len
  
  # Define the first and last particle of the simulation 
  # to do the extraction
  part.1 <- 1
  part.n <- ncid.out$dim$part$len
  
  # Counter for data extraction
  start  <- c( part.1, 1, time.1 )
  count  <- c( part.n, 1, time.n )
  
  
  
  #--- Get 3D position variables
  # Each row is a particle, each column is a time step.
  
  # Particles' longitude
  p.lon  <- ncvar_get( ncid.out, 
                       varid = 'lon',
                       start,
                       count )
  
  # Particles' latitude
  p.lat  <- ncvar_get( ncid.out,
                       varid = 'lat',
                       start,
                       count )
  # Particles' depth
  p.dep  <- ncvar_get( ncid.out,
                       varid = 'zpo',
                       start,
                       count )
  
  
  
  
  #### reshape data ####
  
  # this reshaping version gets rid of the first time step, 
  #that has weird long lat data :
  traj.part <- data.frame( time  = rep( (time.1+1):time.n,     
                                        times = part.n-part.1+1 ),
                           part  = rep( part.1:part.n, 
                                        each  = time.n-time.1 )  ,
                           lon   = as.vector( t(p.lon[,-1])  )   ,
                           lat   = as.vector( t(p.lat[,-1])  )   ,
                           dep   = as.vector( t(p.dep[,-1]) )
  )
  
  # Remove variables used to build the data.frame :
  rm(ncid.out)
  rm(p.dep)
  rm(p.lon)
  rm(p.lat)
  

  
  
  #### spatial extent of the analyses ####
  
  # Choose resolution of the grid : 
  resolution <- 0.05
  
  # set the boudaries of the analyses :
  min_lon <- min(traj.part$lon) 
  max_lon <- max(traj.part$lon) 
  min_lat <- min(traj.part$lat) 
  max_lat <- max(traj.part$lat) 
  
  
  
  
  
  #### Select particles ####
  
  #-- Select only certain times / longitude
  # we don't want the first 12 hours -> the first 6 time steps
  # we don't want the particles close to the border
  # those particles interact with the border and therefore do wrong patterns
  # we don't do it for min latitude, as the ocean basin already removed the border
  #(time.1+1):time.n  # (the time here is the time step)
  
  particles <- traj.part |>
    filter(time > 12           & 
             lon  > min_lon+0.5  & 
             lon  < max_lon-1.5  &
             lat  < max_lat-1.5)
  
  # particles <- traj.part |>
  #   filter(time > 24           & 
  #            lon  > min_lon+0.5  & 
  #            lon  < max_lon-1.5  &
  #            lat  < max_lat-1.5)
  
  rm(traj.part)
  
  
  #-- Remove particles at bottom depth deeper than 600m :
  # Because we don't want to see the aggregations over the ocean basin.
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
  
  
  
  
  #### set the grid ####
  
  # Create a sequence of latitudes and longitudes to set the grid :
  # we approximate the values at 0.05
  min_lat_grid <- round(min_lat*2, 1) / 2
  max_lat_grid <- round(max_lat*2, 1) / 2 - 1.5
  min_lon_grid <- round(min_lon*2, 1) / 2 + 0.5
  max_lon_grid <- round(max_lon*2, 1) / 2 - 1.5
  
  latitudes  <- seq(min_lat_grid, max_lat_grid, by = resolution)
  longitudes <- seq(min_lon_grid, max_lon_grid, by = resolution)
  
  # Create the grid :
  grid <- expand.grid(latitude = latitudes, longitude = longitudes)
  
  
  
  
  #### count the density ####
  
  # parallélisation : 
  detectCores()
  parCluster <- makeCluster( 10, methods = FALSE)
  registerDoParallel( parCluster )
  
  # function that will count the density at each time step : 
  eularian.parallel <- function( i, particles, grid, resolution ){
    
    parts     <- particles[which(particles$time == i),]
    # we approximate the values at 0.05 :
    parts$lon <- round(parts$lon*2, 1) / 2
    parts$lat <- round(parts$lat*2, 1) / 2
    
    
    
    # the variable that will host the values of density
    #parts$pos <- 1 : nrow(parts)
    tme  <- rep(0, nrow(grid)) 
    
    for (j in 1:nrow(grid)) { #for each position, 
      # find the particles at that position 
      # (we use the "near" function bcause the "which" function struggles 
      # with small values)
      bop <- which(
        near(grid[j, "longitude"], parts$lon, tol = 0.04) &
          near(grid[j, "latitude"] , parts$lat, tol = 0.04))
      
      # Another way of doing the computation, 
      # that take longer for now but that may be usefull later 
      # if we take depth into account :
      # bop <- parts |>
      #   filter(lon == grid[j, "longitude"] &
      #            lat == grid[j, "latitude"])
      
      
      if (length(bop)>0){ #if there are particles in this part of the grid, 
        parts  <- parts[-bop,] # remove these lines, so that the computation 
        # time is shorter
        tme[j] <- length(bop)  # keep the nb of particles
      }
    }
    
    out <- data.frame(
      a = tme#, 
      #b = dpt
    )
    colnames(out) [1] = paste0('time_', i, '_res0.5')
    #colnames(out) [2] = paste0('time_', i, '_res0.5')
    
    return(out)
    
  }
  
  # --COMPUTATION-- #
  #tic()
  bep = unique(particles$time)
  z.out <- foreach( i= bep, .packages = 'dplyr') %dopar% { 
    eularian.parallel( i, particles, grid, resolution )
  }
  
  stopCluster( parCluster )
  #toc()
  grid_test <- bind_cols(grid, z.out)
  
  # Save the grided data : 
  save(grid_test, file = paste0("grid_", tools::file_path_sans_ext(current_file)) )
  #load( paste0("grid_", tools::file_path_sans_ext(current_file)) )
  
  grid <- grid_test
  rm(grid_test)
  rm(z.out)

  }
toc()



