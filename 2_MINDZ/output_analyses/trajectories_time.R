#### Initialisation ####
# setwd("~/projects/def-frmap5/MINDZ_andeol")
# library(renv)
# renv::load()
# renv::activate()
# 

# mode plus simple avec tidync : 
library(tidync )
library(dplyr  )       #for data manipulation
library(ggplot2)       #for the plots
library(sf)
setwd("/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/MINDZ/output_analyses")
#setwd("~/projects/def-frmap5/MINDZ_andeol/MINDZ/output_analyses")


world        <- map_data("world") #the borders

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
#days = c(1, 2, 5)
#days = c(2)
ts   = days * 24

for (t in 1:length(ts)) {
   print( paste('  inside time : ',t, ' / ', length(ts)))
  
  # open nc file and transform it into tibble :
  pts <- tidync(nc_file)                              |> # open the file
    tidync::activate("D2")                            |> # activate the "time" variable
    hyper_filter(time = time <= time_step[ts[t]])     |> # just select data of the desired time
    tidync::activate("D0,D1,D2")                      |> # activate the 3 dims (D0 is the particle, D1 is useless, D2 time)
    hyper_tibble(select_var = c("lon", "lat"))           # select the desired variables
  
  #remove the first ts because it has a bug with it longlat
  pts <- pts |>
    filter(time != time_step[1])

  # only plot 0.2% of the particles, select them randomly
  traj <- pts  |>
    filter( part %in% sample( max(pts$part), round(max(pts$part) * 0.002) ) ) |>
    arrange(part, time)
  
  
  
  ggplot( data = traj, aes( x = lon, y = lat ) ) +
    geom_path(    aes( group = part ),
                  linewidth = 0.2) +
                  
       # coast limit and countries borders :
      geom_map(
        data = world, map = world,
        aes(long, lat, map_id = region),
        color = "black", fill = "#B38A3E", linewidth = 0.1) +
      
      # size :
      coord_sf(xlim=c(-70.1, -55), ylim=c(42, 50.9)) +
      
    theme_classic()                                   +
    theme(
      axis.line       = element_line(linewidth = 0.3, colour = "black"),
      axis.text.y     = element_text(size=5),
      axis.text.x     = element_text(size=5), 
      plot.title      = element_text(size=9),
      legend.position = "none" )
  
    ggsave(
      paste0("maptraj_", taxon, '_', days[t], "j.png"),
      dpi    = 200,     # the resolution
      width  = 10,
      height = 7  )
    
}





