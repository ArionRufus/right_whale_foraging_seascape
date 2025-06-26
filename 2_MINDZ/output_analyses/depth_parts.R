
# setwd("/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol")
# setwd("~/projects/def-frmap5/MINDZ_andeol")
# renv::load()
# renv::activate()

#### initialisation ####
library(dplyr)
library(terra)
library(tidync)
#library(stringr)
library(ggplot2)
library(tidyterra)   # to plot the bathymetry map
library(sf)          # to read whale sightings data


setwd("/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/MINDZ/output_analyses")



#-- Select the depth parameters
bdepmax      = 600                    # max bottom depth of the graph
dpt_interval = 10                     # depth resolution of the graph
nb_dpths     = bdepmax / dpt_interval # nb of rows and columns of the graph

# Create the data_frame of the graph :
dpt_profile = matrix(0, nrow = nb_dpths, ncol = nb_dpths)


taxon        <- 'LF'
diapause     <- "0"
current_file <- paste0('test_', taxon, diapause, '_traj.nc')
current_file <- paste0(taxon, diapause, '_4j_traj.nc')
folder_path  <- '../outputs'
nc_file      <- paste0(folder_path, '/', current_file)


# get the time steps of the file :
time_step <- tidync(nc_file) |> # open the file
  tidync::activate("D2")     |> # activate the "time" variable
  hyper_tibble()             |> # change into a tibble with 1 column
  pull()                        # select the first column -> change into a vector

ts <- time_step[time_step > time_step[12]]  # remove the 12 first time steps
rm(time_step)

t = ts[30]
for (t in ts) { # for each time step : 
  print( paste('  inside time : ', which(ts == t), ' / ', length(ts)))
  
  # open nc file and transform it into tibble :
  pts <- tidync(nc_file)                              |> # open the file
    tidync::activate("D2")                            |> # activate the "time" variable
    hyper_filter(time = time == t)                    |> # just select data of the desired time
    tidync::activate("D0,D1,D2")                      |> # activate the 3 dims (D0 is the particle, D1 is useless, D2 time)
    hyper_tibble(select_var = c("zpo", 'bdep'))       |> # select the depth and bottom depth variables
    filter(zpo <= bdepmax)                                   # remove depth data under 400m
  
  # hist(pts$zpo)
  # hist(pts$bdep)
  # mean(pts$zpo/pts$bdep)
  # probs <- pts   |> filter(zpo > bdep)
  # hist(probs$zpo / probs$bdep)
  
  # transform  and bdepth into data frame ids :
  # ex, if a bdepth = 25 and dpt_interval = 20, then this bdepth should be 
  # of the second column. 25/20 = 1.25, round(1.25)+1 = 2
  pts$zpo  = round(pts$zpo  / dpt_interval) + 1
  pts$bdep = round(pts$bdep / dpt_interval) + 1
  
  
  # for each particle, add the count in the data frame : 
  for (i in 1:nrow(pts)) {
    dep_part  = pts$zpo[i] 
    bdep_part = pts$bdep[i] 
    if (dep_part <= nb_dpths & bdep_part <= nb_dpths){ #if the depth and bdepth are in the data frame boundaries
      dpt_profile[dep_part, bdep_part] = dpt_profile[dep_part, bdep_part] + 1
    }
  }    
}


# for each bottom depth, compute a relactive count 
# pour i de 1 a ncol, ,les donnees = les donnes sur la donnee max

for (i in 1:nb_dpths) {
  if ( sum(dpt_profile[, i]) > 0 )
    dpt_profile[, i] = dpt_profile[, i] / max(dpt_profile[, i])
}



tab_dpt_profile = tibble(
  bdep_min = rep( seq(0, bdepmax-dpt_interval, by = dpt_interval), each =  nb_dpths), 
  bdep_max = rep( seq(dpt_interval, bdepmax  , by = dpt_interval), each =  nb_dpths), 
  dep_min  = rep( seq(0, bdepmax-dpt_interval, by = dpt_interval),         nb_dpths), 
  dep_max  = rep( seq(dpt_interval, bdepmax  , by = dpt_interval),         nb_dpths), 
  dens = as.vector(dpt_profile) )



#-- The final plot : 

bop <- ggplot() +
  geom_abline(slope = -1, intercept = 0, linewidth = 1, col = 'grey40') +
  geom_rect(data = tab_dpt_profile, aes(xmin  = bdep_min, xmax = bdep_max, 
                                        ymin  = dep_min,  ymax = dep_max , 
                                        alpha = dens ), fill = '#ee502a') + 
  scale_alpha_continuous(range = c(0, 1)) +
  scale_y_reverse() + 
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        # axis.text.y = element_text (size=14),
        # axis.text.x = element_text (size=14)
        aspect.ratio=1/4) +
  xlab('bottom depth') +
  ylab('depth')

bop
ggsave(paste0("pred_dens_bdep_", taxon, diapause, ".png"),
       bop, width = 6, height = 2)
