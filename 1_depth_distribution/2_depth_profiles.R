##########################################################################
## Creates mean depth distribution profiles of Calanus abundances       ##
## for finmarchicus ans hyperboreus at diff stages active / in diapause ##   
##                                                                      ##
##########################################################################

# Creates a mean depth distribution profile

#-- Initialisation

setwd(".../codes_article/1_depth_distribution")
library(ggplot2)
library(tidyverse)

# Choose the copepod group : 
# Chyp_late, 
# Cfin_late, 
# Chyp_young, 
# Cfin_young 
load("./data/data_Chyp_late.Rdata") 

# In case we just want diapausing / no diapausing data only :
hop <- hop |>
  filter(diapause == "YES")

# if we want just the Late calanus actives :
# ( = Chyp_late active + Cfin_late active)
# hop1 <- hop
# load("./data/data_Cfin_late.Rdata")
# hop <- hop |>
#   filter(diapause == "NO")
# hop = rbind(hop, hop1)

bdmax   <- max(hop$Z_STATION)        # max depth of the graph = maximum. bottom depth  

vertres <- 10                        # vertical resolution
bin     <- 10                        # bottom depth resolution
bds     <- seq(0, bdmax+bin-1, bin)  # bottom depth bins


#-- Create the initial dataframe, 1 line = 1 pixel : 

posits <- NULL
for (i in 1:(length(bds)-1)) {
  vert      <- seq(0, bds[i+1]+vertres-1, vertres)
  temposits <- data.frame(
    StartX = bds [i],                # x = horizontal
    EndX   = bds [i+1],
    StartY = vert[-length(vert)],    # y = vertical
    EndY   = vert[-1], 
    transp = NA,                     # transparency will be related to the density
    ncase  = NA                      # ncase = the number of profiles that we have
  )
  posits <- rbind(posits, temposits)
}
head(posits)


#-- Determine a relative density for each profile
hop$rel_dens <- NA
for (i in unique(hop$ID)) {                                            # For each depth profile,
  temp_profile <- hop$dens_zoo[which(hop$ID == i)]                     # stock the profile
  hop$rel_dens[which(hop$ID == i)] = temp_profile / max(temp_profile)  # rel dens = values of the profile / maximum value of the profile
}



#-- Determine the tranparancy of each pixel from density :

for (i in 1:nrow(posits)) { # For each depth,
  
  # take the upper and lower depth boundary
  do = posits[i,]$StartY
  up = posits[i,]$EndY
  
  # we select relative density values of hop that are included in bottom depth bin, 
  # and if their depth range respect 1 of the 3 conditions :
  # -it is included in the targeted depth range
  # -it is inferior but overlaps more that half of the targeted depth range 
  # -same with superior
  good_pix <- hop |>
    filter(Z_STATION <  posits[i,]$EndX & Z_STATION >= posits[i,]$StartX, # inside bottom depth bin
           Z_MIN     >= do  & Z_MAX <= up |
             Z_MIN   <  do  & Z_MAX >  do + 0.5*(up - do) |
             Z_MAX   >  up  & Z_MIN <  do + 0.5*(up - do) ) |>
    dplyr::select(rel_dens)
  
  #save the mean value and the nb of values : 
  posits[i,]$transp <-  mean(good_pix$rel_dens)
  posits[i,]$ncase  <-  nrow(good_pix)
}

# Set at 0 transparency NA values :
posits$transp[which(is.na(posits$transp))] <- 0



#-- Determine the nb of profiles of each bottom depth :
profiles_bdep <- data.frame(
  bdep = bds[-length(bds)], 
  nb_prof = 0 
)

for (i in 1:nrow(profiles_bdep)) {
  profiles_bdep$nb_prof[i] = n_distinct(
    hop |>
      filter(Z_STATION  <  bds[i+1], 
             Z_STATION  >= bds[i]) |>
      dplyr::select(ID)
  )
}

profiles_bdep <- profiles_bdep |>
  filter(nb_prof != 0)



#-- plot :

# proposed colors for each copepod class, in hex :
# colors : Chyp_young, EDBC4E
#          Cfin_young, 6BAE68
#          Cfin_late,  ED7676
#          Chyp_late,  3293AE

bibop <- ggplot() +
  geom_line(aes(x = 0:600, y = 0:600)) +                                 # the line delimitates bottom depth
  geom_rect(data = posits,                                               # the dataset
            aes(xmin = StartX, xmax = EndX, ymin = StartY, ymax = EndY),
            alpha = posits$transp,                                       # the transparency depending on the density
            fill = '#3293AE') +
  # If desired to write the nb of profiles for each bottom depth :
  # geom_text(data = profiles_bdep,
  #           aes(bdep+2.5, bdep+30, label = nb_prof),
  #           size = 1.5) +
  scale_y_reverse() +                                                    # reverse scale because it is depth and not altitude
  xlim(0, 600) +
  ylim(600, 0) +
  theme_void() +                                                         # theme without anything
  theme(aspect.ratio = 0.25,                                             # specify the ratio to have the same one between all profiles
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.text.y = element_text (size=14),                            # text size of the depth and bottom depth
        axis.text.x = element_text (size=14)
  ) +
  xlab('') +                                                             # no lab text
  ylab('')

bibop

# Save the graph
ggsave("./results_krumh_visu/dens_bdep_color_lhd.png",
       bibop, width = 6, height = 2)

