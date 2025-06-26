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

# setwd("/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol")
# setwd("~/projects/def-frmap5/MINDZ_andeol")
# renv::load()
# renv::activate()

#### initialisation ####
library(dplyr)
library(terra)
library(tidync)
library(tictoc)
library(stringr)
#renv::snapshot()

#setwd("~/Desktop/travail/these/analyses/git/MINDZ/output_analyses")
#setwd("/Ext_20T_andromede/0_ARCTUS_Projects/15_SmartWhales/DATA/MINDZ/output_analyses")
setwd("/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/MINDZ/output_analyses")


# Set the path where the ncdf MINDZ output data are located
folder_path <- "../outputs/to_compute"
#folder_path <- "../outputs/valer_transfer"

# Get a list of files in the folder
all_files <- list.files(folder_path)
YF = F
if (YF == T){max_dep = 150}else{max_dep = 500}

# Loop for all files : 
current_file = all_files[1]
tic() # Optional, to have an idea of the time taken

raster_template <- rast(
  xmin       = -71.1,
  xmax       = -50.7,
  ymin       = 41,
  ymax       = 52.5,
  resolution = 0.05,
  crs        = "EPSG:4326"
)

for (current_file in all_files) {
  
  # Print the actual file analysed :
  print( paste('||TIME : ', which(all_files == current_file), ' / ', length(all_files), ', : ', current_file))
  
  
  
  #### extract data ####
  
  # file path :
  nc_file <-  paste0(folder_path, '/', current_file)
  
  # get the time steps :
  time_step <- tidync(nc_file) |> # open the file
    tidync::activate("D2")     |> # activate the "time" variable
    hyper_tibble()             |> # change into a tibble with 1 column
    pull()                        # select the first column -> change into a vector
  
  ts <- time_step[time_step > time_step[12]]  # enlever les données des 12 premiers TS
  rm(time_step)
  t = ts[7]
  for (t in ts) {
    print( paste('  inside time : ', which(ts == t), ' / ', length(ts)))
    
    # open nc file and transform it into tibble :
    pts <- tidync(nc_file)                              |> # open the file
      tidync::activate("D2")                            |> # activate the "time" variable
      hyper_filter(time = time == t)                    |> # just select data of the desired time
      tidync::activate("D0,D1,D2")                      |> # activate the 3 dims (D0 is the particle, D1 is useless, D2 time)
      hyper_tibble(select_var = c("lon", "lat", "zpo", "bdep")) |> # select the desired variables
      filter(zpo <= max_dep)                                   # remove depth data under 400m
    
    
    # Assign to each particle a depth bin :
    breaks <- seq(0, max_dep, by = 50) # the desired depth bin (50m wide)
    
    pts <- pts |>
      filter(bdep < 600)                                   |>                                 
      mutate(depth_bin = cut(zpo, breaks = breaks))        |> # assign the depth bin
      mutate(depth_bin = str_remove(depth_bin, "\\("))     |> # remove undesired strings
      mutate(depth_bin = str_remove(depth_bin, "\\]"))     |> # remove undesired strings
      mutate(depth_bin = str_replace(depth_bin, ",", "_")) |> # remove undesired strings
      vect(geom = c("lon", "lat"), crs = "EPSG:4326")         # spatialise the tibble
    
    # Transform the tibble into a raster, counting nb a particles for each lon lat dpth_bin :
    if (t == ts[1]) {
      the_r <- rasterize(pts, raster_template, fun = "count", by = "depth_bin")
      names(the_r) <- str_extract(names(the_r), "\\d+_\\d+")
      the_r <- subst(the_r, NA, 0) # replae NAs by 0 (otherwise they will cumulate)
    } else{
      r <- rasterize(pts, raster_template, fun = "count", by = "depth_bin")
      # I am trying to fix what seems to be a bug when using the "by" argument.
      # This hack seems to work for now.
      # Issue opened here: https://github.com/rspatial/terra/issues/1435
      names(r) <- str_extract(names(r), "\\d+_\\d+")
      r <- subst(r, NA, 0)  # replae NAs by 0 (otherwise they will cumulate)
      # Add the values of the new raster to the former one :
      the_r <- the_r + r
      #table(is.na(r[]))
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
  max_raster <- aggregate(the_r, fact=c(1, 1, max_dep/50), fun=bloup)
  
  rm(the_r)
  
  #-- Export the raster
  
  filename <- paste0("grid_3d_", tools::file_path_sans_ext(current_file), ".tif")
  writeRaster(max_raster, filename, overwrite = TRUE)
  
  toc()
}
