################################################################
## Shape the dataset from Krumhansl et al.                    ##
##                                                            ##
## Andeol Bourgouin                                           ##
## andeol.bourgouin.1@ulaval.ca                               ##
################################################################




#### initialisation ####
setwd("/Users/andeolbourgouin/Library/CloudStorage/OneDrive-UniversitéLaval/travail/0_these/1_foraging_areas/codes_article/1_depth_distribution")
library(tidyverse)    # for data manipulation

# Read dataset
data_vert <- read.csv("./data/Calanus_Vertical_Distribution_Dataset_m3.csv")
  


#### Manage time ####
# Extract months : 
split_dates     <- strsplit(data_vert$DATE_TIME, "-")
data_vert$month <- sapply(split_dates, function(x) x[[2]]) |>
  as.numeric()

# Change the way of representing hour
split_dates    <- strsplit(data_vert$DATE_TIME, " ")
data_vert$date <- sapply(split_dates, function(x) ifelse(length(x) >= 1, x[[1]], NA))
data_vert$time <- sapply(split_dates, function(x) ifelse(length(x) == 2, x[[2]], NA))
head(data_vert$time)

# Extract hours
split_hours    <- strsplit(data_vert$time, ":")
data_vert$hour <- sapply(split_hours, function(x) ifelse(length(x) >= 1, x[[1]], NA))

rm(split_hours, split_dates)


table(data_vert$time)
table(data_vert$hour)



#### create global calanus profiles ####
#Cfin mean Calanus finmarchicus, Chyp means Calanus hyperboreus

data_vert$Cfin_tot   <- rowSums(data_vert[, c("Cfin_CI", "Cfin_CII", "Cfin_CIII", 
                                            "Cfin_CIV", "Cfin_CV", "Cfin_CVIf", 
                                            "Cfin_CVIm", "Cfin_CVIu")], 
                                na.rm = TRUE)

data_vert$Cfin_young <- rowSums(data_vert[, c("Cfin_CI", "Cfin_CII", "Cfin_CIII")], 
                                na.rm = TRUE)

data_vert$Cfin_late  <- rowSums(data_vert[, c("Cfin_CIV", "Cfin_CV", "Cfin_CVIf", 
                                             "Cfin_CVIm", "Cfin_CVIu")], 
                                na.rm = TRUE)

data_vert$Chyp_tot   <- rowSums(data_vert[, c("Chyp_CI", "Chyp_CII", "Chyp_CIII", 
                                            "Chyp_CIV", "Chyp_CV", "Chyp_CVIf", 
                                            "Chyp_CVIm", "Chyp_CVIu")], 
                                na.rm = TRUE)

data_vert$Chyp_young <- rowSums(data_vert[, c("Chyp_CI", "Chyp_CII", "Chyp_CIII")], 
                                na.rm = TRUE)

data_vert$Chyp_late <- rowSums(data_vert[, c("Chyp_CIV", "Chyp_CV", "Chyp_CVIf", 
                                             "Chyp_CVIm", "Chyp_CVIu")], 
                               na.rm = TRUE)




#### Split dat. frame into a list of distinct profiles ####
#    TID is a unique identifier of profiles

data_vert$ID <- paste(data_vert$DATE_TIME, data_vert$LONGITUDE, data_vert$LATITUDE)
profiles     <- split( data_vert, as.character(data_vert$ID) )
n_prof       <- length(profiles)
n_prof
# 
# 
# profiles_hyp <- profiles
# profiles_YH  <- profiles
# profiles_LH  <- profiles
# profiles_fin <- profiles
# profiles_YF  <- profiles
# profiles_LF  <- profiles





#### Diapause ####

# Create a diapause variable : 
#data_vert$diapause <- "YES"







#### exclude bad profiles ####

prof.id <- 1:n_prof #for the loop



# for each copepod group 
cop.groups <- c("Cfin_young", "Cfin_late", "Chyp_young", "Chyp_late")

for (actual.group in cop.groups) {
  hop <- data_vert
  prof.id.check <- prof.id # will note as NA the excluded profiles
  
  for( i in prof.id) { # for each profile
    # copepod densities :
    p.obs <- eval(parse(text = paste0("profiles[[i]]$" , actual.group))) 
    

    # Check for profiles with NO copepods | copepods in a SINGLE layer
    if( sum(p.obs,na.rm=TRUE)==0 | sum(p.obs,na.rm=TRUE)==max(p.obs,na.rm=TRUE) ){
      prof.id.check[i] <- NA
      hop <- hop[-which(hop$ID == profiles[[i]]$ID[1]),]
      next()
    }
    
    # Remove profiles with bottom depth > 600m
    profiles[[i]]$Z_STATION[1]
    if( profiles[[i]]$Z_STATION[1] > 600 ){
      prof.id.check[i] <- NA
      hop <- hop[-which(hop$ID == profiles[[i]]$ID[1]),]
      next()
    }
    
    # Remove profiles with a bin of more than 200m : 
    bins <- profiles[[i]]$Z_MAX - profiles[[i]]$Z_MIN
    if( max(bins) > 200 ){
      prof.id.check[i] <- NA
      hop <- hop[-which(hop$ID == profiles[[i]]$ID[1]),]
      next()
    }
    
    # Remove profiles less than 4 bins : 
    nb.bins <- nrow( profiles[[i]] )
    if( nb.bins < 4 ){
      prof.id.check[i] <- NA
      hop <- hop[-which(hop$ID == profiles[[i]]$ID[1]),]
      next()
    }
    
    # Remove profiles with max depth more than 50m shallower than bottom depth : 
    bottom.distance <- profiles[[i]]$Z_STATION[1] - max(profiles[[i]]$Z_MAX)
    if( bottom.distance > 50 ){
      prof.id.check[i] <- NA
      hop <- hop[-which(hop$ID == profiles[[i]]$ID[1]),]
      next()
    }
    
    
    #Check for OVERLAPPING depth strata
    #ZMin : upper depth of each bin of the profile
    #diff : the difference between the values of depth between 
    #each bin and it deeper neighboor.
    zmin.check <- any( diff( profiles[[i]]$Z_MIN ) <= 0 )
    zmax.check <- any( diff( profiles[[i]]$Z_MAX ) <= 0 )
    
    if( zmin.check | zmax.check ){
      prof.id.check[i] <- NA
      hop <- hop[-which(hop$ID == profiles[[i]]$ID[1]),]
      next()
    }
  }
  
  print( paste0(actual.group, ", on the ",n_prof,
                 " profiles, ",
                 length(prof.id.check[-which(is.na(prof.id.check))]), 
                 " profiles are selected") )
  
  # Delete all the unused density data :
  density <- eval(parse( text = paste0("hop$", actual.group) ))
  hop <- hop |>
    select(- c(Cfin_CI, Cfin_CII, Cfin_CIII, Cfin_CIV, Cfin_CV, Cfin_CVIf, 
               Cfin_CVIm, Cfin_CVIu, Cfin_tot, Cfin_late, Cfin_young, 
               Chyp_CI, Chyp_CII, Chyp_CIII, Chyp_CIV, Chyp_CV, Chyp_CVIf, 
               Chyp_CVIm, Chyp_CVIu, Chyp_tot, Chyp_late, Chyp_young, 
               Cglac_CI, Cglac_CII, Cglac_CIII, Cglac_CIV, Cglac_CV, Cglac_CVIf, 
               Cglac_CVIm, Cglac_CVIu) )
  hop$dens_zoo <- density
  
  cop.groups <- c("Cfin_young", "Cfin_late", "Chyp_young", "Chyp_late")
  
  #-- Add diapause variable : 
  hop$diapause = "YES"
  if (actual.group %in% c("Cfin_young", "Cfin_late")) {
    hop$diapause[which(
      hop$REGION == "GSL" &
        hop$month > 3 &
        hop$month < 8
    )] = "NO" 
    
    hop$diapause[which(
      hop$REGION == "SS" &
        hop$month > 2 &
        hop$month < 7
    )] = "NO" 
    
    hop$diapause[which(
      hop$REGION == "GOM" &
        hop$month > 4 &
        hop$month < 10
    )] = "NO" 
    
    hop$diapause[which(
      hop$REGION == "NL" &
        hop$month > 3 &
        hop$month < 10  )] = NA
  }
  
  if (actual.group %in% c("Chyp_young", "Chyp_late")) {
    hop$diapause[which(
      hop$REGION == "GSL" &
        hop$month > 3 &
        hop$month < 7
    )] = "NO" 
    
    hop$diapause[which(
      hop$REGION == "SS" &
        hop$month > 2 &
        hop$month < 7
    )] = "NO" 
    
    hop$diapause[which(
      hop$REGION == "GOM" &
        hop$month > 1 &
        hop$month < 6
    )] = "NO" 
    
    hop$diapause[which(
      hop$REGION == "NL" &
        hop$month > 1 &
        hop$month < 6
    )] = NA
  }
  
  
  
  # Save :
  eval(parse(text = paste0("data_" , actual.group, " <- hop"))) 
  save(hop, file = paste0("./data/data_", actual.group, ".Rdata") )
  }






