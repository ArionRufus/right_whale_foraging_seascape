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
require(marmap)      #to find bathym data
require(ggplot2)     #for the plots
require(gganimate)   #to animate the plots
require(dplyr)       #for data manipulation

setwd("~/Desktop/travail/these/analyses/git/MINDZ/scripts")



#### Load input/output files ####
#chose type of imput file used by the MINDZ model. if false, imput = NEMO
ciops = TRUE

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
  name.out <- 'ZOO_CIOPSE_001_traj.nc'
} else {
  name.out <- 'ZOO_particles_CANOPA_traj.nc'
}
path.out <- '../NETCDF/'
ncid.out <-  nc_open( paste0( path.out, name.out ) )



#--- Extract variables

#### model's bathymetry ####
# Get some information from input file
# to build the model's bathymetry 

# Get vertical layers upper depth
# (MINDZ imputs work with vertical layers, and each layer is  associated to a depth)
if( ciops ){
  v.depth <- ncvar_get( ncid.in, varid = "depth" )
} else {
  v.depth <- ncid.in$dim$depthw$vals
}

# Check whether the "vomecrty" meridional currents variable is present
# and if so, get it
if( ciops ){
  v.check <- as.logical( sum( agrepl( "vo", names(ncid.in$var) ) ) )
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


# Function to get the bottom depth of last wet layer
bathy.fun <- function( val, depth ) {
  bath.i <- which( is.na(val) )[1]
  bathy  <- NA
  if( is.finite(bath.i) && bath.i>1 ) bathy <- depth[bath.i]
  return(bathy)
}

# Get model's bathymetry
bathy <- apply( X = v, MARGIN = c(1,2), bathy.fun, depth = v.depth )
rm(v)



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

#--- Particles variables

# Variables names
v.names <- names(ncid.out$var)

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


# 3D position variables

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

# Particles' env. temperature
#p.tp   <- ncvar_get( ncid.out,
#                     varid = 'tp',
#                     start,
#                     count )

# Particles' mass
if( "mass" %in% v.names) {
  p.mass <- ncvar_get( ncid.out,
                       varid = 'mass',
                       start,
                       count )
}

# Particles' development
if( "stage" %in% v.names) {
  p.stage <- ncvar_get( ncid.out,
                        varid = 'stage',
                        start,
                        count )
}

# Particles' nighttime depth
if( "znight" %in% v.names) {
  p.znight <- ncvar_get( ncid.out,
                         varid = 'znight',
                         start,
                         count )
}



#-----------------------------------------------------------------------------------------------------------------
# Make rapid diagnostic figures !

#--- Histograms of bathymetry & particles' depth

# Model's bathymetry in the GSL area
bathy.max  <- 600  #max( bathy, na.rm=TRUE )

if( ciops ){
  bathy.plot <- bathy
  bathy.plot[bathy.plot>bathy.max] = NA
} else {
  bathy.plot <- bathy[30:170,115:220]
}

depth.plot <-  v.depth[ v.depth <= bathy.max ]

bathy.plot <- hist( bathy.plot,
                    breaks = c(0, depth.plot),
                    plot   = FALSE )$density %>% rev()

z.width <- diff(depth.plot) %>% rev()

# Histograms in probability densities for intercomparison
if( ciops ){
  png("ZOO_CIOPS-E_depth_frequencies.png")
} else {
  png("ZOO_CANOPA_depth_frequencies.png")
}
barplot( bathy.plot,
         width  = z.width,
         xlim   = c(0,0.1),
         ylim   = c(0,bathy.max),
         space  = 0,
         border = NA,
         horiz  = TRUE,
         col    = rgb(1,0,0,0.7),
         axes   = FALSE,
         xlab   = "Prob. density",
         ylab   = "Depth (m)",
         main   = "" )
         
# Trick to get the ticks from 0 to ... by 50 m increments
# needed because of the "reverse" vertical barplot
axis( side   = 2,
      at     = seq( bathy.max, 0, by = -50 ),
      labels = seq( 0, bathy.max, by =  50 ) )

# Particles depth distribution
p.dep.plot <- hist( p.dep,
                    breaks = c(0, depth.plot),
                    plot   = FALSE )$density %>% rev()

barplot( p.dep.plot,
         width  = z.width,
         space  = 0,
         border = NA,
         horiz  = TRUE,
         col    = rgb(0,0,1,0.5),
         axes   = FALSE,
         add    = TRUE )

legend( "bottomright",
        legend = c("Particles","Model's bathymetry"),
        pch    = c(15,15),
        col    = c( rgb(0,0,1,0.5),
                    rgb(1,0,0,0.7) ),
        cex    = 1.2,
        bty    = "n" )

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


#--- Generate animated horizontal path of sub-sampled particles
# Reshape data to use gganimate

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

#why not to use "tstep.plot" ?

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

# Get bathymetric data for ggplot2
if( !exists("gsl.map") ) {
  if( !file.exists("CIOPS-E.bathy.Rdata") ) {
    gsl.map <- getNOAA.bathy( lon1 = -73,
                              lon2 = -49,
                              lat1 =  39,
                              lat2 =  53,
                              resolution = 4 ) %>%
               fortify.bathy()
  
    gsl.bathy <- gsl.map[gsl.map$z<=0,]
    
    save( list = c("gsl.map",
                   "gsl.bathy"),
          file = "CIOPS-E.bathy.Rdata" )
  } else {
    load( file = "CIOPS-E.bathy.Rdata" )
  }
}





#___________________________________________________________________________________


#recuperation des donnees : 
library(terra)
#on veut une carte plus large que d'habitude, pour pouvoir loger des graphiques par dessus le background
kartras = rast("gebco_grand.tif")
kartras[kartras > 0.99] <- NA
#kartras = classify(kartras, cbind(0.99, Inf, NA)) #toutes les valeurs >0.99 deviennent NA

#pour gagner en contraste dans la carte, on limite la profondeur à 4000
#(qui est plus profond qui la profondeur max de notre zone d'étude):')
kartras[kartras <= -600] <- -650

#données monde, pour avoir les frontières :
world <- map_data("world")

#limites de la carte:
lonmin = -73
lonmax = -49
latmin = 39
latmax = 53


#___________________________________________________________________________________




# Subsample data.frame for faster ploting
if ( ciops ) {
  traj.plot <- traj.part  %>%
    filter( part %in% sample( part.n, round(part.n*0.3) ) ) #0..03
} else {
  traj.plot <- traj.part  %>%
    filter( part %in% sample( part.n, round(part.n*0.05) ) )
}

# Plot with gganimate
gg <- ggplot(         data = traj.plot, 
                aes( x = lon, y = lat ) )           +
  coord_quickmap()                                  +
  # geom_raster(  data   = gsl.bathy,
  #               aes( x = x, y = y, fill = z ),
  #               alpha  = 0.5 )                      +
  # scale_fill_gradientn( colours = etopo.colors(4) ) +
  # geom_contour( data   = gsl.map,
  #               aes( x = x, y = y, z = z ),
  #               breaks     = c(0, -100, -500),
  #               colour     = "black",
  #               linewidth  = 0.1 )                  +
  #depth map :
  geom_spatraster(data = kartras)+
  scale_fill_gradientn(colors = gray.colors(200, #nb de couleurs
                                            start = 0.03, #valeur la + sombre
                                            end = .92)) + #valeur la + claire
  #coast limit and countries borders :
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "#a18f7c", size = 0.1)+
  #size :
  coord_sf(xlim=c(-72, -49),
           ylim=c(38, 54)) +
  
  # geom_spatraster(data = kartras)+
  # geom_map(
  #   data = world, map = world,
  #   aes(long, lat, map_id = region),
  #   color = "black", fill = "#a18f7c", size = 0.1)+
  # coord_sf(xlim=c(lonmin, lonmax),ylim=c(latmin, latmax)) +
  
  # geom_path(    aes( group = part ),
  #               colour     = "red", 
  #               alpha      = 0.5)           +
  geom_point(   aes( group = part ),
                #pch        = ".",
                stroke     = 0, 
                size       = 0.15, 
                colour     = '#CC0000')                   +
  theme_classic()                                   +
  theme(
    axis.line   = element_line(size = 0.3, colour = "black"),
    axis.text.y = element_text(size=5),
    axis.text.x = element_text(size=5), 
    plot.title  = element_text(size=9))            +
  theme( legend.position = "none" )                 +
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
           renderer = gifski_renderer("CIOPS-E_trajectories_initial.gif") )
  
} else {
  gganimate :: animate( gg,
           height = 8,
           width  = 10,
           units  = "cm",
           res    = 300,
           renderer = gifski_renderer("ZOO_CANOPA_trajectories.gif") )
}


#--- Generate animated vertical trajectories of sub-sampled particles

# Plot with gganimate
gg <- ggplot(     data = traj.plot, 
                  aes( x = date, y = -dep ) )       +
      geom_line(  aes( group = part ),
                  colour = "darkblue" )            +
      geom_point( aes( group = part ),
                  pch = "." )                      +
      theme_classic()                              +
      theme( legend.position = "none" )            +
      # gganimate specific bits
      transition_reveal(time)                      +
      labs( title = "Output time step: {frame_along}",
            x     = "Time step",
            y     = "Individual depth (m)" )

if ( ciops ) {
  gganimate :: animate( gg,
           height = 10,
           width  = 10,
           units  = "cm",
           res    = 600,
           renderer = gifski_renderer("ZOO_CIOPS-E_DVM.gif") )
} else {
  gganimate :: animate( gg,
           height = 10,
           width  = 10,
           units  = "cm",
           res    = 300,
           renderer = gifski_renderer("ZOO_CANOPA_DVM.gif") )
}


#--- Generate animated growth trajectories of sub-sampled particles

# Plot individual mass with gganimate
gg <- ggplot(     data = traj.plot, 
                  aes( x = date, y = mass ) )      +
      geom_line(  aes( group = part ),
                  colour = "black" )               +
      geom_point( aes( group = part ),
                  pch = "." )                      +
      theme_classic()                              +
      theme( legend.position = "none" )            +
      # gganimate specific bits
      transition_reveal(date)                      +
      labs( title = "Part. growth: {frame_along}",
            x     = "Date",
            y     = "Individual mass (g)" )

#animate( gg,
#         height = 10,
#         width  = 10,
#         units  = "cm",
#         res    = 300,
#         renderer = gifski_renderer("CALANUS_growth.gif") )

# Plot individual stage with gganimate
gg <- ggplot( data = traj.plot, 
              aes( x = date, y = stage ) )     +
  geom_line(  aes( group = part ),
              colour = "black" )               +
  geom_point( aes( group = part ),
              pch = "." )                      +
  theme_classic()                              +
  theme( legend.position = "none" )            +
  # gganimate specific bits
  transition_reveal(date)                      +
  labs( title = "Part. growth: {frame_along}",
        x     = "Date",
        y     = "Individual stage" )

#animate( gg,
#         height = 10,
#         width  = 10,
#         units  = "cm",
#         res    = 300,
#         renderer = gifski_renderer("CALANUS_development.gif") )

