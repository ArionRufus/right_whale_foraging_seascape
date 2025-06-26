#########################################################################
## Initialise the position of particles                                ##
## To be used as imput in the MINDZ model.                             ##
##                                                                     ##
## Andeol Bourgouin                                                    ##
## andeol.bourgouin.1@ulaval.ca                                        ##
#########################################################################


#### Initialisation ####

setwd("~/Desktop/travail/these/analyses/git/MINDZ")
library("ncdf4")
library("raster")
library('dplyr')



#### Set variables ####

# Take a physical imput file to get the boudaries : 
nc_data <- nc_open( "./Data/2022040306_003.nc")

# Extract the boundaries and resolution : 
lon <- ncvar_get(nc_data, "nav_lon")
lat <- ncvar_get(nc_data, "nav_lat")

# Extract a variable with data in 3 dimensions, 
# to have access to the depth :
sal <- ncvar_get(nc_data, "so")

#dpth <- ncvar_get(nc_data, "depth")
#dpth[bop]

# Create the matrix that will stock the particles position :
particles <- lon
particles[] <- 0


# Select borders for the particles release : 
min(lon)
max(lon)
minlon <- 360 - 61.0259 #from degree east to degree west
maxlon <- 360 - 59.8656
minlat <- 46.6424
maxlat <- 47.1928


for (i in 1:nrow(lon)) {
  
  longs <- unique(lon[i,])
  lats  <- unique(lat[i,])
  
  saldep <- sal[i,,]
  result <- saldep |>
    t() |>
    as_tibble() |>
    summarise_all( ~ which( is.na(.))[1] )
  #summarise_all(~ length( na.omit(.) ))
  
  
  butwhere <- which(longs >= minlon & longs <= maxlon & 
                     lats >= minlat &  lats <= maxlat &
                      result > 0)
  particles[i, butwhere] <- 1
}
sum(particles)


#### Special to have a file covering the entire area ####
particles[] <- 2

for (i in 1:nrow(lon)) {

  saldep <- sal[i,,]
  result <- saldep |>
    t() |>
    as_tibble() |>
    summarise_all(~ length( na.omit(.) ))


  butwhere <- which(result == 0)
  particles[i, butwhere] <- 0
}

sum(particles)

output_file <- "initial_positions_allarea_2.csv"



#### Export file ####

output_file <- "initial_positions.csv"
write.table(particles, file = paste0('./run/data/', output_file), 
            sep = ',', row.names = FALSE, col.names = FALSE)

