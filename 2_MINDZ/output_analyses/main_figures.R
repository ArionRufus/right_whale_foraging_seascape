# A ameliorer : 
#- Inscrire les nom de ce qu'on fait tourner dès le debut. 
# a partir de ce nom, les outpus seront écrits tout seuls.
#- La première partie, lecture des données, est un peu compliquée à comprendre. La ranger.
#- Pour l'instant il y a 2 cartes de bathym, il faut que j'en garde qu'une. 
#- En parlant de ces cartes de bathym, il va falloir que je crée un système de projection, 
# pour pouvoir faire correspondre au mieux les données CIOPS. 




#########################################################################
## Reads the output file of the MINDZ model, and do fancy plots.       ##
## (works with CIOPS or NEMO inputs)                                   ##
##                                                                     ##
## Andeol Bourgouin                                                    ##
## andeol.bourgouin.1@ulaval.ca                                        ##
## From Frederic Maps' code                                            ##
#########################################################################


#### initialisation ####
require(gifski)
require(ncdf4)       #to read netcdf files
require(terra)       #for the bathymetry map
library(tidyterra)   #for the bathymetry map
require(ggplot2)     #for the plots
require(gganimate)   #to animate the plots
require(dplyr)       #for data manipulation

setwd("~/Desktop/travail/these/analyses/git/MINDZ/output_analyses")

# Model characteristics :
# used in the code to choose the right code, and to set the names of the outputs
ciops        = TRUE  # type of imput file used by the MINDZ model. if false, imput = NEMO
part         = 'YF'  # particle class : LH=late_hyperboreus
init_distrib = 'all' # initial distribution



#### Load input/output MINDZ files ####


# netCDF input file
if( ciops ){
  name.in <- '2022040100_001.nc'
} else {
  name.in <- 'L_2011_M2_10d.nc'
}

path.in <- '../Data/'
ncid.in <-  nc_open( paste0( path.in, name.in ) )

# netCDF output file
if( ciops ){
  name.out <- 'CIOPSE_LH_traj.nc'
} else {
  name.out <- 'CIOPSE_LH_traj.nc'
}
path.out <- '../outputs/'
ncid.out <-  nc_open( paste0( path.out, name.out ) )





#### Extract bathymetry ####

# Get some information from input file to build the model's bathymetry 

# Get vertical layers upper depth
# (MINDZ imputs work with vertical layers, and each layer is  associated to a depth)
if( ciops ){
  v.depth <- ncvar_get( ncid.in, varid = "depth" )
} else {
  v.depth <- ncid.in$dim$depthw$vals
}

# Check whether the "vomecrty" meridional currents variable is present
# and if so, get it. 
# These data will be used to determine the bathymetry : 
if( ciops ){
  v.check <- as.logical( sum( agrepl( "vo",       names(ncid.in$var) ) ) )
} else {
  v.check <- as.logical( sum( agrepl( "vomecrty", names(ncid.in$var) ) ) )
}

if( v.check ) {
  
  if( ciops ){
    v.ndims <- ncid.in$var$vo$ndims #number of dimensions associated to vo
  } else {
    v.ndims <- ncid.in$var$vomecrty$ndims
  }
  
  v.start <- rep( 1, v.ndims )
  
  if( ciops ){
    v.count <- ncid.in$var$vo$size #size of each dimension
    
    v <- ncvar_get( ncid.in,
                    varid = "vo",
                    start = v.start,
                    count = v.count )
  } else {
    v.count <- ncid.in$var$vomecrty$size
    
    v <- ncvar_get( ncid.in,
                    varid = "vomecrty",
                    start = v.start,
                    count = v.count )
  }
  v[v==0] <- NA

} else {
  stop( paste( "Problem with the velocity 4D variable read from the netCDF input file", name.in ) )
}


# Function to get the bottom depth of last wet layer : 
# the first NA value is the layer of bottom depth. Then we look within the 
# 'depth' data which is the depth associated to that layer. 
bathy.fun <- function( val, depth ) {
  bath.i <- which( is.na(val) )[1]
  bathy  <- NA
  if( is.finite(bath.i) && bath.i>1 ) bathy <- depth[bath.i]
  return(bathy)
}

# Get model's bathymetry :
bathy <- apply( X = v, MARGIN = c(1,2), bathy.fun, depth = v.depth )
rm(v) #remove the data on currents, it's large and useless





#### extract data from output file ####

# => Need to know the structure of the file; check with "ncdump -h outfile.nc"

#--- Information for extracting variables

# Number of particles to extract data for
part.1 <- 1
part.n <- ncid.out$dim$part$len

# Sequence dimension to use 
seq    <- 1

# Define the first time step of the simulation and 
# the numbers of time steps we want for data extraction
time.1 <- 1 
time.n <- ncid.out$dim$time$len

# Counter for data extraction
start  <- c( part.1, seq, time.1 )
count  <- c( part.n, seq, time.n )

# Variables names
v.names <- names(ncid.out$var)


#--- Get time

# Timestep for figures
tstep <- ncvar_get( ncid.out,
                    varid = 'time' )

# Time origin from file attributes
if( ciops ){
  time.origin <- ncatt_get( ncid.out, 
                          varid   = "time", 
                          attname = "time_origin" )$value %>%
               as.POSIXct(.)
} else {
  time.origin <- ncatt_get( ncid.out, 
                            varid   = "time", 
                            attname = "time_origin" )$value %>%
    as.POSIXct(.,
               format  = "%Y-%b-%d %H:%M:%S" )
}

# Time format conversion
tstep <- as.POSIXct( tstep,
                     origin = time.origin )

# Change into integer because of gganimate issue with transition_reveal
tstep.plot <- as.character(tstep) %>% 
               sub(":00:00","",.) %>%
              gsub("-",     "",.) %>%
              gsub(" ",     "",.) %>%
              as.integer(.)



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

# Some other variables, for CIOPS analyses we don't need them : 

# Particles' env. temperature
#p.tp   <- ncvar_get( ncid.out,
#                     varid = 'tp',
#                     start,
#                     count )

# Particles' mass
# if( "mass" %in% v.names) {
#   p.mass <- ncvar_get( ncid.out,
#                        varid = 'mass',
#                        start,
#                        count )
# }

# Particles' development
# if( "stage" %in% v.names) {
#   p.stage <- ncvar_get( ncid.out,
#                         varid = 'stage',
#                         start,
#                         count )
#}

# Particles' nighttime depth
# if( "znight" %in% v.names) {
#   p.znight <- ncvar_get( ncid.out,
#                          varid = 'znight',
#                          start,
#                          count )
# }




#### Histograms ####

#--- Histogram of bathymetry & particles' depth

# As the histogramm has the bathymetry AND the particles, the 2 ones 
#won't have the same nb of observations. So we create a histogramm data_frame 
#for both and then we transform in porportion. 
#bathymetry is in proportion 1.5 to be a bit larger.
depth.plot = c(v.depth, max(p.dep, na.rm=T))

#histogram data-frames : 
bathy.plot <- hist( bathy,
                    breaks = c(0, depth.plot),
                    plot   = FALSE )$density
bathy.plot <- bathy.plot * 1.5 / sum(bathy.plot)

pdep.plot <- hist( p.dep,
                    breaks = c(0, depth.plot),
                    plot   = FALSE )$density
pdep.plot <- pdep.plot * 1 / sum(pdep.plot)

#put everything in a data-frame ready to plot :
hist_data <- data.frame(
  dpth_min = c(0, depth.plot[-length(depth.plot)]), 
  dpth_max = depth.plot,
  bathy = bathy.plot, 
  parts = pdep.plot
)

#plot it :
if( ciops ){
  png(paste0("./results/CIOPS-E_depth_frequencies", part, init_distrib, ".png"))
} else {
  png("./results/ZOO_CANOPA_depth_frequencies.png")
}
ggplot(data = hist_data) + 
  # bathy :
  geom_rect(aes(xmin = dpth_max, xmax = dpth_min, ymin = 0, ymax = bathy, 
                fill = "Model's bathymetry"), alpha = 0.4) +
  # particles :
  geom_rect(aes(xmin = dpth_max, xmax = dpth_min, ymin = 0, ymax = parts, 
                fill = 'particles'), alpha = 0.4) +
  #turn the plot :
  coord_flip() +
  scale_x_reverse() +
  #Theme :
  labs(x ="Depth", y = "Density") +
  scale_fill_manual(name = '', breaks=c("Model's bathymetry", 'particles'),
                    values=c("Model's bathymetry"='#0000CC', 'particles'='#CC0000')) +
  theme(aspect.ratio = 1.3, 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line   = element_line(size = 0.3, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10), 
        legend.position = c(0.8, 0.15) )

dev.off()


#--- Histograms of particles' mass

#plot.new()
#png("CALANUS_mass_frequencies.png")

#hist( p.mass,
#      border = NA,
#      col    = rgb(0,0.5,0,0.7),
#      xlab   = "Mass (g)",
#      ylab   = "Prob. density",
#      main   = "" )

#dev.off()


#--- Histograms of particles' development

#plot.new()
#png("CALANUS_stage_frequencies.png")

#hist( p.stage,
#      border = NA,
#      col    = rgb(0,0.5,0,0.7),
#      xlab   = "Development stage",
#      ylab   = "Prob. density",
#      main   = "" )

#dev.off()





#### horizontal particles path ####

#-- Reshape data to use gganimate

#this reshaping version gets rid of the first time step, 
#that has weird long lat data :
traj.part <- data.frame( time  = rep( (time.1+1):time.n,     
                                      times = part.n-part.1+1 ),
                         date  = rep( tstep[-1],
                                       times = part.n-part.1+1 ),
                         part  = rep( part.1:part.n, 
                                       each  = time.n-time.1 ),
                         lon   = as.vector( t(p.lon[,-1])  ),
                         lat   = as.vector( t(p.lat[,-1])  ),
                         dep   = as.vector( t(p.dep[,-1]) )
                         ) %>%
  filter( lat>0 )

#this one is when there is no problem with the first time step :
# traj.part <- data.frame( time  = rep( time.1:time.n,     
#                                       times = part.n-part.1+1 ),
#                          date  = rep( tstep,
#                                       times = part.n-part.1+1 ),
#                          part  = rep( part.1:part.n, 
#                                       each  = time.n-time.1+1 ),
#                          lon   = as.vector( t(p.lon)  ),
#                          lat   = as.vector( t(p.lat)  ),
#                          dep   = as.vector( t(p.dep)  ) ) %>%
#              filter( lat>0 )




#-- Prepare bathymetry data for the bacground map :
kartras = rast("gebco_grand.tif") #downloaded before from gebco
kartras[kartras > 0.99] <- NA     #remove terrestrial data

#to have more contrast on the map, put a limit to bathymetry that will 
#make the ocean basins all black :
kartras[kartras <= -600] <- -650

#world data, to plot borders :
world <- map_data("world")


#-- Subsample data.frame for faster ploting :

if ( ciops ) {
  traj.plot <- traj.part  %>%
    filter( part %in% sample( part.n, round(part.n*0.03) ) ) #0..03
} else {
  traj.plot <- traj.part  %>%
    filter( part %in% sample( part.n, round(part.n*0.05) ) )
}


#-- Plot with gganimate :

gg <- ggplot( data = traj.plot, aes( x = lon, y = lat ) ) +
  coord_quickmap()                                  +
  # bathymetry map :
  geom_spatraster(data = kartras)+
  scale_fill_gradientn(colors = gray.colors(200,            #nb of colours
                                            start = 0.03,   #darkest value
                                            end   = .92)) + #whitest value
  # coast limit and countries borders :
  geom_map( data = world, map = world, aes(long, lat, map_id = region),
    color = "black", fill = "#a18f7c", size = 0.1) +
  # limits :
  coord_sf(xlim=c(-72, -49),
           ylim=c(38, 54)) +
  # particles positions :
  geom_point( aes( group = part ),
                stroke     = 0,       #trick to have points shorter than 0.1
                size       = 0.15,    #point size
                colour     = '#CC0000')             +
  theme_classic()                                   +
  theme(
    axis.line       = element_line(size = 0.3, colour = "black"),
    axis.text.y     = element_text(size=5),
    axis.text.x     = element_text(size=5), 
    plot.title      = element_text(size=9),
    legend.position = "none" )                      +
  # gganimate specific bits
  transition_reveal(date)                           + 
  shadow_wake(wake_length = 0.3, 
              alpha       = TRUE,
              size        = FALSE, 
              wrap        = FALSE) +
  labs( title = "Time : {frame_along}",
        x     = "",
        y     = "" )
  

if ( ciops ) {
  gganimate :: animate( gg,
           height = 8,
           width  = 10,
           units  = "cm",
           res    = 600,
           renderer = gifski_renderer(
             paste0("./results/CIOPS-E_circulation_", part, init_distrib, ".gif")) )
  
} else {
  gganimate :: animate( gg,
           height = 8,
           width  = 10,
           units  = "cm",
           res    = 300,
           renderer = gifski_renderer("ZOO_CANOPA_trajectories.gif") )
}





#-- Special plot for the poster :
library(ggnewscale)

#less particles to plot : 
traj.poster <- traj.part  %>%
  filter( time >= 121,
          time <= 132,
    part %in% sample( part.n, round(part.n*0.2) ) ) #0..03

#add a layer to put over the basin : 
kart = rast("gebco_grand.tif") #downloaded before from gebco
kart[kart > 0.99] <- NA     #remove terrestrial data
kart[kart < -450 & kart > -550] <- -450
basin <- kart
kart[kart <= -550] <- NA
basin[basin >=  -450] <- NA
basin[basin <  -450] <- 1
basin <- as.factor(basin)

# ggplot()+
#   geom_spatraster(data = basin)

ggplot( data = traj.poster, aes( x = lon, y = lat ) ) +
  coord_quickmap()                                  +
  # bathymetry map :
  geom_spatraster(data = kart, maxcell = 9e+06)+
  scale_fill_gradientn(colors = gray.colors(200,            #nb of colours
                                            start = 0.03,   #darkest value
                                            end   = .92)) + #whitest value
  
  # coast limit and countries borders :
  geom_map( data = world, map = world, aes(long, lat, map_id = region),
            color = "black", fill = "#B38A3E", size = 0.1) +
  # particles positions :
  # geom_point( aes( group = part ),
  #             stroke     = 0,       #trick to have points shorter than 0.1
  #             size       = 0.15,    #point size
  #             colour     = '#CC0000')             +
  geom_path(    aes( group = part ),
                linewidth = 0.1,
                colour = "#78325D")          +
  
  # the basin
  new_scale_fill() +
  geom_spatraster(data = basin) +
  scale_fill_manual(values = "#155176", na.value = NA) +
  
  # limits :
  coord_sf(xlim=c(-72, -49),
           ylim=c(38, 54)) +

  theme_classic()                                   +
  theme(
    axis.line       = element_line(size = 0.3, colour = "black"),
    axis.text.y     = element_text(size=5),
    axis.text.x     = element_text(size=5), 
    plot.title      = element_text(size=9),
    legend.position = "none" )




#### vertical particles path ####

# Plot with gganimate
gg <- ggplot(     data = traj.plot, 
                  aes( x = date, y = -dep ) )       +
      geom_line(  aes( group = part ),
                  size       = 0.1,    #point size
                  colour = "darkblue" )            +
      geom_point( aes( group = part ),
                  size       = 0.1,    #point size
                  pch = "." )                      +
      theme_classic()                              +
      theme( legend.position = "none" )            +
      # gganimate specific bits
      transition_reveal(time)                      +
      labs( title = "time step: {frame_along}",
            x     = "Time step",
            y     = "Individual depth (m)" )

if ( ciops ) {
  gganimate :: animate( gg,
           height = 10,
           width  = 10,
           units  = "cm",
           res    = 600,
           renderer = gifski_renderer(
             paste0("./results/CIOPS-E_vertical_distrib_", part, init_distrib, ".gif")) )
} else {
  gganimate :: animate( gg,
           height = 10,
           width  = 10,
           units  = "cm",
           res    = 300,
           renderer = gifski_renderer("ZOO_CANOPA_DVM.gif") )
}






#### Other plots ####

# We are not interested in those plots with CIOPS-E outputs

# Plot individual mass with gganimate
# gg <- ggplot(     data = traj.plot, 
#                   aes( x = date, y = mass ) )      +
#       geom_line(  aes( group = part ),
#                   colour = "black" )               +
#       geom_point( aes( group = part ),
#                   pch = "." )                      +
#       theme_classic()                              +
#       theme( legend.position = "none" )            +
#       # gganimate specific bits
#       transition_reveal(date)                      +
#       labs( title = "Part. growth: {frame_along}",
#             x     = "Date",
#             y     = "Individual mass (g)" )

#animate( gg,
#         height = 10,
#         width  = 10,
#         units  = "cm",
#         res    = 300,
#         renderer = gifski_renderer("CALANUS_growth.gif") )



# Plot individual stage with gganimate
# gg <- ggplot( data = traj.plot, 
#               aes( x = date, y = stage ) )     +
#   geom_line(  aes( group = part ),
#               colour = "black" )               +
#   geom_point( aes( group = part ),
#               pch = "." )                      +
#   theme_classic()                              +
#   theme( legend.position = "none" )            +
#   # gganimate specific bits
#   transition_reveal(date)                      +
#   labs( title = "Part. growth: {frame_along}",
#         x     = "Date",
#         y     = "Individual stage" )

#animate( gg,
#         height = 10,
#         width  = 10,
#         units  = "cm",
#         res    = 300,
#         renderer = gifski_renderer("CALANUS_development.gif") )




#### export data ####
save(traj.part, file = paste0("traj_part", part, init_distrib, ".Rdata") )
