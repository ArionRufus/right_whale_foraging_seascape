#########################################################################
## The Orsnstein - Ulhenbek depth distribution model,                  ##
## applied to data_krumh depth profiles                                ##
##                                                                     ##
## Andeol Bourgouin                                                    ##
## andeol.bourgouin.1@ulaval.ca                                        ##
## From Frederic Maps' code                                            ##
#########################################################################

#How it works :
#- import the data
#- transform it into a list of profiles
#- transform them into an array of profiles too
#- create the functions
#- test them
#- run pricefit to find theoretical profiles that fit the best, and store them in a list
#- transform the results into an array
#- merge obs and sim arrays
#- create a table of environmental and simulation parameters
#- export z.opt (for security), and the last 2 tables.



#### initialisation ####

library(tidyverse)  # for data manipulation
library(ecolMod)    # for the pricefit (= gradient descent) function
library(parallel)   # to make parallel computations
library(doParallel) # to make parallel computations

setwd("/Users/andeolbourgouin/Library/CloudStorage/OneDrive-UniversitéLaval/travail/0_these/1_foraging_areas/codes_article/1_depth_distribution")



#### get and select data ####
#--- Open the large vertical profiles spreadsheet
# (choose the copepod group)
taxon <- "YF"

if (taxon == "YF"){load("./data/data_Cfin_young.Rdata")}
if (taxon == "LF"){load("./data/data_Cfin_late.Rdata" )}
if (taxon == "LH"){load("./data/data_Chyp_late.Rdata" )}



#--- Split this data.frame into a list of distinct profiles
#    ID is a unique identifier of profiles
data_vert <- hop |> arrange(ID)
profiles  <- split( data_vert, as.factor(data_vert$ID) )
n_prof    <- length(profiles)

print( paste0( "We have detected ", n_prof, " profiles, taxon : ", taxon) )




#### functions ####

#--- Ornstein - Uhlenbeck

# set a depth position to each particle with the Ornstein - Uhlenbeck equation
OU_module <- function( Z, sigma, theta, mu, dt ){
  with( as.list(Z, sigma, theta, mu, dt), {
    
    # Z     : the variable, in our case depth                                  [m]
    # sigma : > 0; volatility => average magnitude of the random fluctuations  [time^-0.5]
    # theta : > 0; rate of mean reversion => parameterize the "swimming speed" [time^-1]
    # mu    : long-term mean of the process, in our case the target depth      [m]
    # dt    : time-step                                                        [time]
    
    M  <- mu + ( Z - mu ) * exp(-theta*dt)
    
    SD <- sigma * sqrt( ( 1 - exp(-2*theta*dt) ) * 0.5 / theta )
    
    Z  <- rnorm( length(Z), mean=M, sd=SD )
    
    return( Z )
  })
}



#--- Vertical binning and counts

# We use the `hist()` function to compute the frequencies (counts)
# in depth bin identical to those sampled
bin.fun <-  function( new_interval, sim ) {
  with( as.list( c( new_interval, sim ) ), {
    
    # build the breaks vector for hist()
    ret <- hist( sim,
                 breaks = new_interval,
                 plot   = FALSE )$counts
    
    # Output: counts for the same stratas as sampled
    return(ret)
  })
}



#--- New boundaries

# modify the boundaries of the simulations and of the observations,
# so that they can be compared.
# only used outside the optimisation.
bound.fun <-  function( z.min, z.max, p.obs, z.out.1, z.out.2 ) {
  with( as.list( c( z.min, z.max, p.obs, z.out.1, z.out.2  ) ), {
    
    new_interval <- c(z.min, max(z.max))
    new_p.obs    <- p.obs
    
    # if the min depth of the data is deeper than min predicted depth,
    # add the 0 - min pred depth bin to the data :
    if( min(new_interval) > min(c(z.out.1, z.out.2)) ) {
      new_interval <- c( min(c(z.out.1, z.out.2)) , new_interval )
    }
    
    # if the max depth of the data is shallower than max predicted depth,
    # add the max depth data - max pred depth bin to the data :
    if( max(new_interval) < max(c(z.out.1, z.out.2)) ) {
      new_interval <- c( new_interval , max(c(z.out.1, z.out.2)) )
    }
    
    # Output
    return(list(new_interval = new_interval,
                new_p.obs    = new_p.obs))
  })
}



#--- Mean absolute error

# alternative to Root mean square deviation.
# adaptated with inegal depth bins inside the dataset.
# deals only with two datasets with the same depth bins.
mae.fun <-  function( obs, sim, interval ) {
  with( as.list( c( obs, sim, interval ) ), {
    
    lgth <- interval[-1] - interval[-length(interval)]
    ret  <- sum( abs( sim - obs ) * lgth ) / length(obs)
    
    # Output
    return(ret)
  })
}



#--- Cost function for later use with pricefit()

#- take the values of OU parameters (for 2 OU distributions, to consider
#  the bimodal distribution) ,
#- compute the two distribution with the OU function ,
#- put it in the right depth bin ,
#- apply the cost function.

cost.fun <- function( pars, p.obs, z.min, z.max, bot.depth ) {
  
  with( as.list( c( pars, p.obs, z.min, z.max, bot.depth ) ), {
    
    #in case we want to debug :
    # print(mu1)
    # print(mu2)
    # print(sigma1)
    # print(sigma2)
    # print(dens1)
    # print(dens2)
    # print("________________________")
    
    
    #--Simulate the profiles :
    
    dens <- dens1 + dens2        # the cumulated density of the simulated profiles
    p    <- dens1 / dens         # the proportion of density
    p.simul <- 10000             # the nb of particles to simulate
    p.1     <-    p  * p.simul   # the nb for the first mode
    p.2     <- (1-p) * p.simul   # and for the second mode
    
    # first mode simulation :
    z.out.1 <- OU_module( Z     = rep(mu1,ceiling(p.1)),
                          mu    = mu1,
                          sigma = sigma1,
                          theta = 3e-4,
                          dt    = 600 )
    # we use ceiling in case p.1 is less than 1
    
    # second mode simulation :
    z.out.2 <- OU_module( Z     = rep(mu2,ceiling(p.2)),
                          mu    = mu2,
                          sigma = sigma2,
                          theta = 3e-4,
                          dt    = 600 )
    
    
    
    #-- Remove all particles outside the surface - bottom area :
    
    # change depth of particles over the surface :
    if (min(z.out.1) < 0){
      patate      <- length(which(z.out.1 < 0))
      prof_buffer <- 1 : 20
      youpi       <- sample(prof_buffer, size=patate, replace = T)
      z.out.1[which(z.out.1 < 0)] <- youpi }
    
    if (min(z.out.2) < 0){
      patate      <- length(which(z.out.2 < 0))
      prof_buffer <- 1 : 20
      youpi       <- sample(prof_buffer, size=patate, replace = T)
      z.out.2[which(z.out.2 < 0)] <- youpi }
    
    
    # change depth of particles under bottom depth :
    if (max(z.out.1) > bot.depth){
      patate      <- length(which(z.out.1 > bot.depth))
      prof_buffer <- (0.9 * bot.depth) : bot.depth
      youpi       <- sample(prof_buffer, size=patate, replace = T)
      z.out.1[which(z.out.1 > bot.depth)] <- youpi }
    
    if (max(z.out.2) > bot.depth){
      patate      <- length(which(z.out.2 > bot.depth))
      prof_buffer <- (0.9 * bot.depth) : bot.depth
      youpi       <- sample(prof_buffer, size=patate, replace = T)
      z.out.2[which(z.out.2 > bot.depth)] <- youpi }
    
    
    #-- If the max observed depth is over the bottom depth,
    # remove the predicted depth over this depth,
    # and update the density (because it is not used to determine the cost,
    # but it still exists) :
    nb_parts_before <- length(c(z.out.1, z.out.2))
    if (max(z.out.1) > max(z.max)){ z.out.1 <- z.out.1[-which(z.out.1 > max(z.max))] }
    if (max(z.out.2) > max(z.max)){ z.out.2 <- z.out.2[-which(z.out.2 > max(z.max))] }
    
    # same for surface :
    if (min(z.out.1) < min(z.min)){ z.out.1 <- z.out.1[-which(z.out.1 < min(z.min))] }
    if (min(z.out.2) < min(z.min)){ z.out.2 <- z.out.2[-which(z.out.2 < min(z.min))] }
    
    nb_parts_after <- length(c(z.out.1, z.out.2))
    prop_ext <- (nb_parts_before - nb_parts_after) / nb_parts_before
    
    # if there are still some particles to determine the cost,
    # then determine the proportion of particles to remove
    if (prop_ext != 1) {dens <- dens - (prop_ext * dens)}
    
    #-- Bin the profiles, so they have the same depth bins than the observed ones :
    
    sim      <- c(z.out.1, z.out.2)
    interval <- c(z.min, max(z.max))
    p.sim    <- bin.fun(interval, sim)
    
    # put the number of simulated partiles equal to the total nb of wanted simulated particles :
    if (prop_ext != 1) {
      nb_parts <- length(sim)
      p.sim    <- p.sim * dens / nb_parts }
    
    
    #-- Determine the cost :
    
    # here we use the nrmsd, but it doesn't work for inegal rectangles
    cost <- mae.fun( p.obs, p.sim, interval )
    
    
    return(cost)
  })
}




#### test pricefit function for one profile ####
#select a profile
i <- 37

p.obs     <-     profiles[[i]]$dens_zoo
z.min     <-     profiles[[i]]$Z_MIN
z.max     <-     profiles[[i]]$Z_MAX
bot.depth <- max(profiles[[i]]$Z_STATION)


# Initial parameters vector with limits
OU.pars <- c( mu1     = 20,
              mu2     = max(z.max)-10,
              sigma1  = 0.15,
              sigma2  = 0.15,
              dens1   = sum(p.obs) / 2,
              dens2   = sum(p.obs) / 2,
              theta   = 3e-4,
              dt      = 600)

OU.min  <- c( min(z.min),
              100,
              0.12,
              0.12,
              1,
              0,
              OU.pars["theta"],
              OU.pars["dt"] )

OU.max  <- c( 80,
              min(c( bot.depth/2 , max(z.max) )),
              1.5,
              1.5,
              1.5 * sum(p.obs),
              1.5 * sum(p.obs),
              OU.pars["theta"],
              OU.pars["dt"] )


# run pricefit
# pricefit will try to find the best parameters by using a lot of random values
# of these variables inbetween the bin minpar - maxpar.
# we run the pricefit function several times, stock the cost and just keep 
# the run with the lowest cost

temp_opts   <- rep( list(NULL), 4 )
costs = NULL

for (b in 1:4) {
  #On fait le pricefit
  #on loge les resultats dans une liste
  #et on garde le cost
  temp_opts[[b]] <- pricefit( par       = OU.pars,
                              minpar    = OU.min,
                              maxpar    = OU.max,
                              func      = cost.fun,
                              p.obs     = p.obs,
                              z.min     = z.min,
                              z.max     = z.max,
                              bot.depth = bot.depth,
                              numiter   = 10000)

  costs = c(costs, temp_opts[[b]]$cost)

}

# just keep the run with the lowest cost
z.opt <- temp_opts[[ which.min(costs) ]]

# Show "best" set of parameters
optpar <- round( z.opt$par, 4 )
optpar

mu1    <- optpar["mu1"]
mu2    <- optpar["mu2"]
sigma1 <- optpar["sigma1"]
sigma2 <- optpar["sigma2"]
dens1  <- optpar["dens1"]
dens2  <- optpar["dens2"]
theta  <- optpar["theta"]
dt     <- optpar["dt"]


#--Simule the profiles :

dens    <- dens1 + dens2     #the cumulated density of the simulated profiles
p       <- dens1 / dens      #the proportion of density
p.simul <- 10000             #the nb of particles to simulate
p.1     <-    p  * p.simul   #the nb for the first mode
p.2     <- (1-p) * p.simul   #and for the second mode

#first mode simulation :
z.out.1 <- OU_module( Z     = rep(mu1,ceiling(p.1)),
                      mu    = mu1,
                      sigma = sigma1,
                      theta = 3e-4,
                      dt    = 600 )
#we use ceiling in case p.1 is less than 1

#second mode simulation :
z.out.2 <- OU_module( Z     = rep(mu2,ceiling(p.2)),
                      mu    = mu2,
                      sigma = sigma2,
                      theta = 3e-4,
                      dt    = 600 )



#-- change all particles outside the surface - bottom area :


#change depth of particles over the surface :
if (min(z.out.1) < 0){
  patate      <- length(which(z.out.1 < 0))
  prof_buffer <- 1 : 20
  youpi       <- sample(prof_buffer, size=patate, replace = T)
  z.out.1[which(z.out.1 < 0)] <- youpi }

if (min(z.out.2) < 0){
  patate      <- length(which(z.out.2 < 0))
  prof_buffer <- 1 : 20
  youpi       <- sample(prof_buffer, size=patate, replace = T)
  z.out.2[which(z.out.2 < 0)] <- youpi }


#change depth of particles under bottom depth :
if (max(z.out.1) > bot.depth){
  patate      <- length(which(z.out.1 > bot.depth))
  prof_buffer <- (0.9 * bot.depth) : bot.depth
  youpi       <- sample(prof_buffer, size=patate, replace = T)
  z.out.1[which(z.out.1 > bot.depth)] <- youpi }

if (max(z.out.2) > bot.depth){
  patate      <- length(which(z.out.2 > bot.depth))
  prof_buffer <- (0.9 * bot.depth) : bot.depth
  youpi       <- sample(prof_buffer, size=patate, replace = T)
  z.out.2[which(z.out.2 > bot.depth)] <- youpi }




#-- Bin the profiles, so they have the same depth bins than the observed ones :

sim      <- c(z.out.1, z.out.2)
interval <- c(z.min, max(z.max))

if (  max(z.max) < bot.depth ){
  interval <- c(interval, bot.depth)
}

if (  min(z.min) > 0 ){
  interval <- c(0, interval)
}

p.sim    <- bin.fun(interval, sim)

#put the number of simulated partiles equal to the total nb of wanted simulated particles :
nb_parts <- length(sim)
p.sim    <- p.sim * dens / nb_parts



#plot the sim profile
opbinatgg <- data.frame(
  StartX <- interval[-length(interval)],
  EndX   <- interval[-1],
  StartY <- c(0),
  EndY   <- p.sim)

colnames(opbinatgg) = c('StartX', 'EndX', 'StartY', 'EndY')

#the observed profile
attest <- data.frame(
  StartX <- z.min,
  EndX   <- z.max,
  StartY <- c(0),
  EndY   <- p.obs)

colnames(attest) = c('StartX', 'EndX', 'StartY', 'EndY')


ggplot() +
  #observed data :
  geom_rect(data = attest, aes(xmin = StartX, xmax = EndX, ymin = StartY, ymax = EndY),
            fill ="#f79378") +

  #pedicted data :
  geom_rect(data = opbinatgg, aes(xmin = StartX, xmax = EndX, ymin = StartY, ymax = EndY),
            fill ="blue", color = 'blue',linewidth = 0.3,  alpha = 0.3) +
  coord_flip() +
  scale_x_reverse() +
  ggtitle(paste0("Bined simulated and observed profiles : ", profiles[[i]]$profile[1])) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        plot.title       = element_text(size=13),
        axis.text.y      = element_text(size=15),
        axis.text.x      = element_text(size=15))










#### Run the Price algorithm calibration on all valid profiles ####

# Initial parameters vector with limits
OU.pars <- c( mu1     = 20,   # m
              mu2     = NA,
              sigma1  = 0.15,
              sigma2  = 0.15,
              dens1   = NA,
              dens2   = NA,
              theta   = 3e-4,
              dt      = 600)  # time

OU.min  <- c( mu1    = NA,
              mu2    = NA,
              sigma1 = 0.01,
              sigma2 = 0.01,
              dens1  = 0.01, #to get rid of a double 0
              dens2  = 0.01,
              theta  = OU.pars["theta"],
              dt     = OU.pars["dt"] )

OU.max  <- c( mu1     = NA,
              mu2     = NA,
              sigma1  = 0.8,
              sigma2  = NA,
              dens1   = NA,
              dens2   = NA,
              theta   = OU.pars["theta"],
              dt      = OU.pars["dt"] )



z.opt   <- rep( list(NULL), n_prof )
prof.id <- 1:length(profiles) #for the loop


# !!! the following loop will run for a long time and necessitate memory !!!
# !!! make sure to know what it does before to run it                    !!!

if( !file.exists( paste0("./data/calibration_krumh_dia_", taxon, ".Rdata") ) ) {
  
  #--- Parallelization
  #    Create a cluster (a group of R processes):
  #    * make ncores-1 of them
  #    * their communication type is "SOCK" (compatible with Windows, Linux, MacOS)
  #    * Bonus: skip loading the `methods` package for some extra speed :)
  
  ncores     <- detectCores() # find the nb of cores of the computer
    
  parCluster <- makeCluster( ncores-1, methods = FALSE)
  registerDoParallel( parCluster )
  
  # Parallelized call to pricefit()
  price.parallel <- function( i, prof.id, OU.pars, OU.min, OU.max, profiles ){
    with( as.list( c(i, prof.id, OU.pars, OU.min, OU.max, profiles) ), {
      
      set.seed(42)
      
      # set the boundaries of mu - mu2:
      bot.depth <- max(profiles[[i]]$Z_STATION)
      z.min <- profiles[[i]]$Z_MIN
      z.max <- profiles[[i]]$Z_MAX
      p.obs <- profiles[[i]]$dens_zoo
      
      
      
      OU.min["mu1"]    <- min(z.min)
      OU.max["mu1"]    <- min(c(80, bot.depth-50))
      OU.min["mu2"]    <- min(c(100, bot.depth-50))
      OU.max["mu2"]    <- max(c( bot.depth/2 , max(z.max) ))
      OU.max["dens1"]  <- 1.5 * sum(p.obs)
      OU.max["dens2"]  <- 1.5 * sum(p.obs)
      OU.max["sigma2"] <- 0.01 * bot.depth
      OU.pars["mu2"]   <- max(z.max) - 5
      OU.pars["dens1"] <- sum(p.obs) / 2
      OU.pars["dens2"] <- sum(p.obs) / 2
    
      
      # if we want to do several departures of the optim function,
      # to prevent fits in local minimums :
      temp_opts   <- rep( list(NULL), 30 )
      costs = NULL
      
      for (b in 1:30) {
        temp_opts[[b]] <- pricefit( par       = OU.pars,
                                    minpar    = OU.min,
                                    maxpar    = OU.max,
                                    func      = cost.fun,
                                    p.obs     = p.obs,
                                    z.min     = z.min,
                                    z.max     = z.max,
                                    bot.depth = bot.depth,
                                    numiter   = 10000)
        
        costs = c(costs, temp_opts[[b]]$cost)
        
      }
      
      z.opt <- temp_opts[[ which.min(costs) ]]
      return(z.opt)
    } )
  }
  
  z.out <- foreach( i=prof.id, .packages=c('dplyr','ecolMod')) %dopar% {
    print(i)
    price.parallel( i, prof.id, OU.pars, OU.min, OU.max, profiles )
  }
  
  stopCluster( parCluster )
  
  for( i in 1:length(prof.id) ) {
    j          <- prof.id[i]
    z.opt[[j]] <- z.out[[i]]
  }
  rm(z.out)
    
  save(z.opt, file = paste0("./data/calibration_krumh_dia_", taxon, ".Rdata") )
} else {
  load(       file = paste0("./data/calibration_krumh_dia_", taxon, ".Rdata") )
}









#### Compute the profiles ####


# dataframe stocking the predicted profiles :
pred_profiles <- data.frame(
  profile    = NULL,
  order      = NULL,
  sim_startX = NULL,
  sim_endX   = NULL,
  sim_startY = NULL,
  sim_endY   = NULL
)

# and the predicted parameters :
pred_params <- data.frame(
  profile =  NULL,
  order   =  NULL,
  p       =  NULL,
  mu1     =  NULL,
  mu2     =  NULL,
  sigma1  =  NULL,
  sigma2  =  NULL,
  nrmsd   =  NULL
)



#--- Get observations
for( i in 1:length(z.opt)) {
  
  # If the model didn't converge :
  if( is.null( z.opt[[i]] ) ) {
    print(paste0(i, ", model didn't converge"))
    
    lgt <- nrow(profiles[[i]])
    
    pred_profiles <- rbind(pred_profiles,
                           data.frame(
                             profile    = profiles[[i]]$ID[1],
                             order      = i,
                             sim_startX = NA,
                             sim_endX   = NA,
                             sim_startY = NA,
                             sim_endY   = NA
                           ))
    
    pred_params <- rbind(pred_params,
                         data.frame(
                           profile =  profiles[[i]]$ID[1],
                           order   =  i,
                           p       =  NA,
                           mu1     =  NA,
                           mu2     =  NA,
                           sigma1  =  NA,
                           sigma2  =  NA,
                           nrmsd   =  NA
                         ))
    next()
    
  } else { # if the model converged :
    bot.depth <- max(profiles[[i]]$Z_STATION)
    z.min     <- profiles[[i]]$Z_MIN
    z.max     <- profiles[[i]]$Z_MAX
    p.obs     <- profiles[[i]]$dens_zoo
    
    
    # extract the optimised parameters
    mu1    <- z.opt[[i]]$par["mu1"]
    mu2    <- z.opt[[i]]$par["mu2"]
    dens1  <- z.opt[[i]]$par["dens1"]
    dens2  <- z.opt[[i]]$par["dens2"]
    sigma1 <- z.opt[[i]]$par["sigma1"]
    sigma2 <- z.opt[[i]]$par["sigma2"]
    theta  <- z.opt[[i]]$par["theta"]
    dt     <- z.opt[[i]]$par["dt"]
    
    
    # if max depth <= half of bottom depth (or inverse), dens2 = 0
    if(max(z.max) <= (bot.depth/2)){
      dens2 = 0}
    if(min(z.min) >= (bot.depth/2)){
      dens1 = 0}
    
    #--Simulate the profiles :
    
    dens    <- dens1 + dens2     # the cumulated density of the simulated profiles
    p       <- dens1 / dens      # the proportion of density
    
    p.simul <- 10000             # the nb of particles to simulate
    p.1     <-    p  * p.simul   # the nb for the first mode
    p.2     <- (1-p) * p.simul   # and for the second mode
    
    # first mode simulation :
    z.out.1 <- OU_module( Z     = rep(mu1,ceiling(p.1)),
                          mu    = mu1,
                          sigma = sigma1,
                          theta = 3e-4,
                          dt    = 600 )
    # we use ceiling in case p.1 is less than 1
    
    # second mode simulation :
    z.out.2 <- OU_module( Z     = rep(mu2,ceiling(p.2)),
                          mu    = mu2,
                          sigma = sigma2,
                          theta = 3e-4,
                          dt    = 600 )
    
    
    
    #-- Remove all particles outside the surface - bottom area :
    
    # remove particles uover the surface :
    if (min(z.out.1) < 0){ z.out.1 <- z.out.1[-which(z.out.1 < 0)] }
    if (min(z.out.2) < 0){ z.out.2 <- z.out.2[-which(z.out.2 < 0)] }
    
    # remove particles under bottom depth :
    if (max(z.out.1) > bot.depth){ z.out.1 <- z.out.1[-which(z.out.1 > bot.depth)] }
    if (max(z.out.2) > bot.depth){ z.out.2 <- z.out.2[-which(z.out.2 > bot.depth)] }
    
    
    
    #-- Bin the profiles, so they have the same depth bins than the observed ones :
    
    sim      <- c(z.out.1, z.out.2)
    interval <- c(z.min, max(z.max))
    
    if (  max(z.max) < bot.depth ){
      interval <- c(interval, bot.depth)
    }
    
    if (  min(z.min) > 0 ){
      interval <- c(0, interval)
    }
    
    p.sim    <- bin.fun(interval, sim)
    
    # put the number of simulated particles equal to the total nb of
    # wanted simulated particles :
    nb_parts <- length(sim)
    p.sim    <- p.sim * dens / nb_parts
    
    
    
    # plot the sim profile
    pred_profiles <- rbind(pred_profiles,
                           data.frame(
                             profile = profiles[[i]]$ID[1],
                             order   = i,
                             StartX  = interval[-length(interval)],
                             EndX    = interval[-1],
                             StartY  = c(0),
                             EndY    = p.sim) )
    
    colnames(pred_profiles) = c('profile', 'order', 'StartX', 'EndX', 'StartY', 'EndY')
    
    
    # add parameters to the dataset
    pred_params <- rbind(pred_params,
                         data.frame(
                           profile =  profiles[[i]]$ID[1],
                           order   =  i,
                           p       =  p,
                           mu1     =  mu1,
                           mu2     =  mu2,
                           sigma1  =  sigma1,
                           sigma2  =  sigma2,
                           nrmsd   =  z.opt[[i]]$cost
                         ))
    
  }
}






#### add env variables pred_params ####

#--- select the explanatory variables to extract
pred_params$LONGITUDE            <- NA
pred_params$LATITUDE             <- NA
pred_params$Z_STATION            <- NA
pred_params$T_MEAN               <- NA
pred_params$T_MIN                <- NA
pred_params$T_MAX                <- NA
pred_params$T_0_50               <- NA
pred_params$T_NEAR_BOTTOM        <- NA
pred_params$PSAL_MEAN            <- NA
pred_params$PSAL_MIN             <- NA
pred_params$PSAL_MAX             <- NA
pred_params$PSAL_0_50            <- NA
pred_params$PSAL_NEAR_BOTTOM     <- NA
pred_params$DENSITY_MEAN         <- NA
pred_params$MIN_SUB_T            <- NA
pred_params$DENSITY_NEAR_BOTTOM  <- NA
pred_params$DEPTH_MIN_SUB_T      <- NA
pred_params$CIL_UP               <- NA
pred_params$CIL_LOW              <- NA
pred_params$CIL_THICK            <- NA
pred_params$T_B_CIL              <- NA
pred_params$MEAN_KD490           <- NA
pred_params$SPAN                 <- NA
pred_params$BLOOM_SEASON         <- NA
pred_params$date                 <- NA
pred_params$time                 <- NA
pred_params$hour                 <- NA
pred_params$REGION               <- NA
pred_params$diapause             <- NA




#--- extract the explanatory variables

for( i in prof.id) {
  pred_params$LONGITUDE [i]           <- profiles[[i]]$LONGITUDE[1]
  pred_params$LATITUDE [i]            <- profiles[[i]]$LATITUDE[1]
  pred_params$Z_STATION [i]           <- profiles[[i]]$Z_STATION[1]
  pred_params$T_MEAN [i]              <- profiles[[i]]$T_MEAN[1]
  pred_params$T_MIN [i]               <- profiles[[i]]$T_MIN[1]
  pred_params$T_MAX [i]               <- profiles[[i]]$T_MAX[1]
  pred_params$T_0_50 [i]              <- profiles[[i]]$T_0_50[1]
  pred_params$T_NEAR_BOTTOM [i]       <- profiles[[i]]$T_NEAR_BOTTOM[1]
  pred_params$PSAL_MEAN [i]           <- profiles[[i]]$PSAL_MEAN[1]
  pred_params$PSAL_MIN [i]            <- profiles[[i]]$PSAL_MIN[1]
  pred_params$PSAL_MAX  [i]           <- profiles[[i]]$PSAL_MAX[1]
  pred_params$PSAL_0_50 [i]           <- profiles[[i]]$PSAL_0_50[1]
  pred_params$PSAL_NEAR_BOTTOM  [i]   <- profiles[[i]]$PSAL_NEAR_BOTTOM[1]
  pred_params$DENSITY_MEAN  [i]       <- profiles[[i]]$DENSITY_MEAN[1]
  pred_params$MIN_SUB_T [i]           <- profiles[[i]]$MIN_SUB_T[1]
  pred_params$DENSITY_NEAR_BOTTOM [i] <- profiles[[i]]$DENSITY_NEAR_BOTTOM[1]
  pred_params$DEPTH_MIN_SUB_T [i]     <- profiles[[i]]$DEPTH_MIN_SUB_T[1]
  pred_params$CIL_UP [i]              <- profiles[[i]]$CIL_UP[1]
  pred_params$CIL_LOW   [i]           <- profiles[[i]]$CIL_LOW[1]
  pred_params$CIL_THICK  [i]          <- profiles[[i]]$CIL_THICK[1]
  pred_params$T_B_CIL  [i]            <- profiles[[i]]$T_B_CIL[1]
  pred_params$MEAN_KD490  [i]         <- profiles[[i]]$MEAN_KD490[1]
  pred_params$SPAN  [i]               <- profiles[[i]]$SPAN[1]
  pred_params$BLOOM_SEASON  [i]       <- profiles[[i]]$BLOOM_SEASON[1]
  pred_params$date  [i]               <- profiles[[i]]$date[1]
  pred_params$time [i]                <- profiles[[i]]$time[1]
  pred_params$hour  [i]               <- profiles[[i]]$hour[1]
  pred_params$REGION [i]              <- profiles[[i]]$REGION[1]
  pred_params$diapause [i]            <- profiles[[i]]$diapause[1]
}





#### export data ####
save(pred_profiles, file = paste0("./data/pred_profiles_krumh_", taxon, ".Rdata") )
save(pred_params,   file = paste0("./data/pred_params_krumh_", taxon, ".Rdata"  ) )
load(paste0("./data/pred_profiles_krumh_", taxon, ".Rdata") )
load(paste0("./data/pred_params_krumh_", taxon, ".Rdata"  ) )





#### Plot the profiles ####
# Plot all the profiles in a long pdf

# Select if we want to plot active stages, diapause or all

d_state <- "D" # one of D for diapause, A for active, or all

if (taxon   == 'YF'){d_state        = "all"}
if (d_state == 'D' ){diapause_state = "YES"}
if (d_state == 'A' ){diapause_state = "NO" }

if (d_state == 'all'){
  # params :
  params_diapause <- pred_params          |>
    select( c(mu1, mu2, p, profile) )
  
  # obs profiles :
  obs_diapause    <- data_vert            |>
    mutate( profile    = ID)              |>
    select( c(Z_MIN, Z_MAX, dens_zoo, Z_STATION, profile) )
  
  # sim profiles :
  pred_diapause  <- pred_profiles         |>
    select( c(StartX, EndX, StartY, EndY, profile) )
} else {
  
  # params :
  params_diapause <- pred_params          |>
    filter( diapause  == diapause_state)        |>
    select( c(mu1, mu2, p, profile) )
  
  # obs profiles :
  obs_diapause    <- data_vert            |>
    filter( diapause  == diapause_state)        |>
    mutate( profile    = ID)              |>
    select( c(Z_MIN, Z_MAX, dens_zoo, Z_STATION, profile) )
  
  # sim profiles :
  pred_diapause  <- pred_profiles         |>
    filter( profile %in% obs_diapause$profile) |>
    select( c(StartX, EndX, StartY, EndY, profile) )
}



# title variation :
plot_size <- n_distinct(params_diapause$profile) / 1.72
title     <- paste0("./results_krumh_optim/simprofs_krumh_dia_",
                    taxon, d_state, ".pdf")


# the plot :

simprofiles_plot = ggplot() +
  #observed profiles :
  geom_rect(data = obs_diapause, aes(xmin = Z_MIN, xmax = Z_MAX, ymin = 0, ymax = dens_zoo),
            fill ="#f79378") +
  #sim profiles :
  geom_rect(data = pred_diapause, aes(xmin = StartX, xmax = EndX,
                                      ymin = StartY, ymax = EndY),
            fill ="#78dcf7", col = "#78dcf7", linewidth = 0.1, alpha = 0.4) +
  
  #bottom depth line :
  geom_vline(data = obs_diapause, aes(xintercept=Z_STATION), linetype="dashed",
             color = "#386087", linewidth = 0.2) +
  
  #mu depth lines :
  geom_vline(data = params_diapause, aes(xintercept=mu1), linetype="dashed",
             color = "#78dcf7", linewidth = 0.3, alpha = 0.8) +
  geom_vline(data = params_diapause, aes(xintercept=mu2), linetype="dashed",
             color = "#78dcf7", linewidth = 0.3, alpha = 0.8) +
  
  #p value :
  geom_text(data = params_diapause, x=0,  y=Inf,
            aes(label = paste0("p = ", round(p, 2))),
            size = 2, hjust = "inward") +
  
  #turn the plot :
  coord_flip() +
  scale_x_reverse(limits = c(NA, 0)) +
  
  theme(aspect.ratio = 1.5,
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line   = element_line(size = 0.3, colour = "black"),
        axis.text.y = element_text(size=5),
        axis.text.x = element_text(size=5),
        strip.text  = element_text(size=3)
  ) +
  facet_wrap(~profile, ncol = 4, scales = "free")


# save the big plot in a pdf :
ggsave(title, simprofiles_plot, device = "pdf",
       height = plot_size, width = 6, scale = 1, limitsize = F)

