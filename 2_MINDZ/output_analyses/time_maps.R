
#setwd("/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol")
#setwd("~/projects/def-frmap5/MINDZ_andeol")
#renv::load()
#renv::activate()

#### initialisation ####
library(dplyr)
library(terra)
library(tidync)
library(tictoc)
library(stringr)
library(ggplot2)
library(tidyterra)   # to plot the bathymetry map
library(sf)          # to read whale sightings data


setwd("/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/MINDZ/output_analyses")

world <- map_data("world") #the borders


raster_template <- rast(
  xmin       = -71.1,
  xmax       = -50.7,
  ymin       = 41,
  ymax       = 52.5,
  resolution = 0.05,
  crs        = "EPSG:4326"
)




taxon        <- 'YF'
current_file <- 'YF1_20j_traj.nc'
folder_path  <- '../outputs'
nc_file      <- paste0(folder_path, '/', current_file)




# get the time steps :
time_step <- tidync(nc_file) |> # open the file
  tidync::activate("D2")     |> # activate the "time" variable
  hyper_tibble()             |> # change into a tibble with 1 column
  pull()                        # select the first column -> change into a vector

days = c(1, 2, 5, 10, 20)
ts   = days * 24

for (tlaps in 1:length(ts)) {
  print( paste('  inside time : ',tlaps, ' / ', length(ts)))
  
  #### extract data ####
  loop_time <- time_step[2 : ts[tlaps]] # get rid of first time that has wrong long lat
  
  for (t in loop_time) {

    # open nc file and transform it into tibble :
    pts <- tidync(nc_file)                              |> # open the file
      tidync::activate("D2")                            |> # activate the "time" variable
      hyper_filter(time = time == t)                    |> # just select data of the desired time
      tidync::activate("D0,D1,D2")                      |> # activate the 3 dims (D0 is the particle, D1 is useless, D2 time)
      hyper_tibble(select_var = c("lon", "lat", "zpo")) |> # select the desired variables
      filter(zpo <= 400)                                   # remove depth data under 400m
    

    
    # Assign to each particle a depth bin :
    breaks <- seq(0, 400, by = 50) # the desired depth bin (50m wide)
    
    pts <- pts |>
      mutate(depth_bin = cut(zpo, breaks = breaks))        |> # assign the depth bin
      mutate(depth_bin = str_remove(depth_bin, "\\("))     |> # remove undesired strings
      mutate(depth_bin = str_remove(depth_bin, "\\]"))     |> # remove undesired strings
      mutate(depth_bin = str_replace(depth_bin, ",", "_")) |> # remove undesired strings
      vect(geom = c("lon", "lat"), crs = "EPSG:4326")         # spatialise the tibble
    
    # Transform the tibble into a raster, counting nb a particles for each lon lat dpth_bin :
    if (t == loop_time[1]) {
      the_r <- rasterize(pts, raster_template, fun = "count", by = "depth_bin")
      names(the_r) <- str_extract(names(the_r), "\\d+_\\d+")
    } else{
      r <- rasterize(pts, raster_template, fun = "count", by = "depth_bin")
      # I am trying to fix what seems to be a bug when using the "by" argument.
      # This hack seems to work for now.
      # Issue opened here: https://github.com/rspatial/terra/issues/1435
      names(r) <- str_extract(names(r), "\\d+_\\d+")
      
      # Add the values of the new raster to the former one :
      the_r <- the_r + r
      rm(r)
    }
    rm(pts)
  }
  
  
  # function that takes the biggest density of all dpth bins :
  bloup <- function(x){ 
    if (length(x) == length(which(is.na(x)))){
      NA
    } else {
      max(x, na.rm = TRUE)
    }
  }
  
  #if we want to have the depth istead of the density : 
  # bloup <- function(x){
  #   if (length(x) == length(which(is.na(x)))){
  #     NA
  #   } else {
  #     which.max(x)
  #   }
  # }
  
  # Aggregate the raster to find the maximum value across the 'nlayer' dimension :
  max_raster <- aggregate(the_r, fact=c(1, 1, 8), fun=bloup)
  rm(the_r)
  
  
  
  
  
  
  
  #____________________________________________________________________________________
  
  
  
  # on a notre raster, maintenant on veut en faire une carte qu'on va exporter. 
  
  
  ggplot( ) +
    geom_spatraster(data = max_raster) +
    
    # density of particles :
    # geom_tile(data = grid_test, aes(x = longitude,   y = latitude,
    #                                 fill = proj, color = proj), color = NA) +
    scale_fill_viridis_c(option = 'G', begin = 0, na.value = 'gray20') +
    
    # coast limit and countries borders :
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "black", fill = "#B38A3E", size = 0.1) +
    
    
    # size :
    coord_sf(xlim=c(-70.1, -50), ylim=c(40.5, 50.9)) +
    
    # look of the plot :
    theme_classic() +
    theme(
      axis.line   = element_line(linewidth = 0.3, colour = "black"),
      axis.text.y = element_text(size=6),
      axis.text.x = element_text(size=6),
      plot.title  = element_text(size=9)
    )
  
  # Save the plot : 
  ggsave(
    paste0("mapdenstraj_", taxon, '_', days[tlaps], "j.png"),
    dpi    = 100,     # the resolution
    width  = 10,
    height = 10  )
  
}



#____________________________________________________________________________________



