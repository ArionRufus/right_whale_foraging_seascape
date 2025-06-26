#########################################################################
## Analyses and Regressions from the outputs of the O-U model,         ##
## Applied to data_krumh depth profiles                                ##
##                                                                     ##
## Andeol Bourgouin                                                    ##
## andeol.bourgouin.1@ulaval.ca                                        ##
#########################################################################




#### initialisation ####
setwd("/Users/andeolbourgouin/Library/CloudStorage/OneDrive-UniversitéLaval/travail/0_these/1_foraging_areas/analyses/profondeur/dfo")
library(ggplot2)
library(scales)       # to put scale with commas in plots
library(RColorBrewer) # for color palettes
library(tidyverse)



#-- Get and select data 
# load the observed and simulated profiles (atgg / atpred) :
# load( file = paste0("./results_krumh_optim/atgg.Rdata") )
# load( file = paste0("./results_krumh_optim/atpred.Rdata") )

# (CAN CHANGE TAXA & STAGE HERE)
tax   <- "Chyp"
sta   <- "CIV_VI"
grp   <- "LH"

#load the parameters data (pred_params) :
load( file = paste0("./data/pred_params_krumh_", grp, ".Rdata") )


# If we want to have active stages of finmarchicus + hyperboreus :
# tax   <- "Chyp"
# load( file = paste0("./data/pred_params_krumh_", "LH", ".Rdata") )
# # faire tourner le morceau 'diapause' du script
# LHdat <- pred_params
# tax   <- "Cfin"
# load( file = paste0("./data/pred_params_krumh_", "LF", ".Rdata") )
# # faire tourner le morceau 'diapause' du script
# pred_params <- rbind(LHdat, pred_params)
# tax   <- "Calanus"
# sta   <- "CIV_VI"
# rm(LHdat)
# grp   <- "Calanus"

#change values of some mu ton NA when p si close to 1 or 0
#remove mus when p close to 1 or 0 :
pred_params$mu2[ which(pred_params$p >= 0.95) ] = NA
pred_params$mu1[ which(pred_params$p <= 0.05) ] = NA

#change names of the variables mu and sigma :
names(pred_params)[names(pred_params) == 'mu1']    <- 'mu_upper'
names(pred_params)[names(pred_params) == 'mu2']    <- 'mu_lower'
names(pred_params)[names(pred_params) == 'sigma1'] <- 'sigma_upper'
names(pred_params)[names(pred_params) == 'sigma2'] <- 'sigma_lower'

# Extract months : 
split_dates       <- strsplit(pred_params$date, "-")
pred_params$month <- sapply(split_dates, function(x) x[[2]]) |>
  as.numeric()

# Create a diapause variable : 
pred_params$diapause <- "YES"

if (tax == "Cfin") {
  pred_params$diapause[which(
    pred_params$REGION == "GSL" &
      pred_params$month > 3 &
      pred_params$month < 8
  )] = "NO" 
  
  pred_params$diapause[which(
    pred_params$REGION == "SS" &
      pred_params$month > 0 &
      pred_params$month < 7
  )] = "NO" 
  
  pred_params$diapause[which(
    pred_params$REGION == "GOM" &
      pred_params$month > 4 &
      pred_params$month < 10
  )] = "NO" 
  
  pred_params$diapause[which(
    pred_params$REGION == "NL" &
      pred_params$month > 2 &
      pred_params$month < 10  )] = NA
}

if (tax == "Chyp") {
  pred_params$diapause[which(
    pred_params$REGION == "GSL" &
      pred_params$month > 3 &
      pred_params$month < 7
  )] = "NO" 
  
  pred_params$diapause[which(
    pred_params$REGION == "SS" &
      pred_params$month > 0 &
      pred_params$month < 6
  )] = "NO" 
  
  pred_params$diapause[which(
    pred_params$REGION == "GOM"
  )] = NA
  
  pred_params$diapause[which(
    pred_params$REGION == "NL"
  )] = NA
  
  pred_params$diapause[which(
    pred_params$REGION == "NL" &
      pred_params$month > 4 &
      pred_params$month < 7
  )] = "NO"
  
  pred_params$diapause[which(
    pred_params$REGION == "NL"&
      pred_params$month > 9 &
      pred_params$month < 12
  )] = "YES"
}



##





#### Diapause ####


#mu upper :
ggplot(data = pred_params, aes(x=diapause, y = mu_upper)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color='#54AEFF', size=0.7, alpha=0.8) +
  scale_y_reverse() +
  ggtitle("mu_upper depth and Diapause state", 
          subtitle = paste0(tax, ", stage ", tax)) +
  xlab("Diapause") + ylab("depth") +
  theme_minimal()


#mu lower :
ggplot(data = pred_params, aes(x=diapause, y = mu_lower)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color='#FF5459', size=0.7, alpha=0.8) +
  scale_y_reverse() +
  ggtitle("mu_lower depth and Diapause state", 
          subtitle = paste0(tax, ", stage ", tax)) +
  xlab("Diapause") + ylab("depth") +
  theme_minimal()


#p :
# for P, I remove the NA of the plot
ggplot(data = pred_params |> filter(!is.na(diapause)), aes(x=diapause, y = p)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color='black', size=0.7, alpha=0.8) +
  ggtitle("p and Diapause state", 
          subtitle = paste0(tax, ", stage ", sta)) +
  xlab("Diapause") + ylab("p") +
  theme_minimal()

ggsave(paste0("./results_krumh_optim/p_diapause", grp, ".jpeg"), 
       height = 5, width = 5)


#bottom depth :
ggplot(data = pred_params, aes(x=diapause, y = Z_STATION)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color='black', size=0.7, alpha=0.8) +
  ggtitle("bottom depth and Diapause state", 
          subtitle = paste0(tax, ", stage ", tax)) +
  xlab("Diapause") + ylab("bottom depth") +
  scale_y_reverse() +
  theme_minimal()





#### model mu ####

# I use weighted least squares regression, to integrate the increase of 
# variance with increasing bottom depth


#-- mu upper

# Find data with NA :
pos_mu_up <- which(is.na(pred_params$mu_upper) | is.na(pred_params$diapause))

# Define weights to use :
reg_mu_up <-  lm(mu_upper[-pos_mu_up] ~ 0 + Z_STATION[-pos_mu_up], data = pred_params)
wt_mu_up <- 1 / lm(abs(reg_mu_up$residuals) ~ reg_mu_up$fitted.values)$fitted.values^2

# Regression with diapause effect
WLM_mu_up_D <-  lm(mu_upper[-pos_mu_up] ~ Z_STATION[-pos_mu_up] + 
                     diapause[-pos_mu_up]:Z_STATION[-pos_mu_up] , 
                   data = pred_params, weights=wt_mu_up)
AIC(WLM_mu_up_D)

# Regression without diapause effect
WLM_mu_up <-  lm(mu_upper[-pos_mu_up] ~ Z_STATION[-pos_mu_up], 
                 data = pred_params, weights=wt_mu_up)
AIC(WLM_mu_up)

# keep the coefficients 
intercept_up <- WLM_mu_up_D$coefficients[1]

if (AIC(WLM_mu_up_D) < AIC(WLM_mu_up)) {  #if the diapause model fit better
  bdep_up     <- WLM_mu_up_D$coefficients[3]
  YES_up       <- WLM_mu_up_D$coefficients[2]  #coef rectification if not in diapause
} else {
  bdep_up     <- WLM_mu_up$coefficients[2]
  YES_up       <- 0
}


#-- mu lower

# Find data with NA :
pos_mu_lo <- which(is.na(pred_params$mu_lower) | is.na(pred_params$diapause))

# Define weights to use :
reg_mu_lo <-  lm(mu_lower[-pos_mu_lo] ~ 0 + Z_STATION[-pos_mu_lo], data = pred_params)
wt_mu_lo <- 1 / lm(abs(reg_mu_lo$residuals) ~ reg_mu_lo$fitted.values)$fitted.values^2

# Regression with diapause effect
WLM_mu_lo_D <-  lm(mu_lower[-pos_mu_lo] ~ diapause[-pos_mu_lo]:Z_STATION[-pos_mu_lo] +
                     Z_STATION[-pos_mu_lo], data = pred_params, weights=wt_mu_lo)
AIC(WLM_mu_lo_D)

# Regression without diapause effect
WLM_mu_lo <-  lm(mu_lower[-pos_mu_lo] ~ Z_STATION[-pos_mu_lo], 
                 data = pred_params, weights=wt_mu_lo)
AIC(WLM_mu_lo)

# keep the coefficients 
intercept_lo <- WLM_mu_lo_D$coefficients[1]

if (AIC(WLM_mu_lo_D) < AIC(WLM_mu_lo)) {  #if the diapause model fit better
  bdep_lo     <- WLM_mu_lo_D$coefficients[2]
  YES_lo       <- WLM_mu_lo_D$coefficients[3] #coef rectification if not in diapause
} else {
  bdep_lo     <- WLM_mu_lo$coefficients[2]
  YES_lo       <- 0
}




#-- graphs
#the dashed lines of the middle boundaries : 
mid_sep_up <- data.frame(
  bdp = c(0, 130, 600), 
  dpt = c(0, 80, 80)
)
mid_sep_lo <- data.frame(
  bdp = c(0, 150, 600), 
  dpt = c(0, 100, 100)
)


# Regression in diapause :
# get the 3 points x y of the corresation path :
changing_up <- (intercept_up + 10*(bdep_up+YES_up)) / (1-(bdep_up+YES_up))
upper_line <- data.frame(
  bdp = c(0, changing_up+10, 600), 
  dpt = c(0, changing_up, intercept_up + 600*(bdep_up+YES_up))
)

changing_lo <- (intercept_lo + 10*(bdep_lo+YES_lo)) / (1-(bdep_lo+YES_lo))
lower_line <- data.frame(
  bdp = c(0, changing_lo+10, 600), 
  dpt = c(0, changing_lo, intercept_lo + 600*(bdep_lo+YES_lo))
)

ggplot(data = pred_params[which(pred_params$diapause == "YES"),]) +
  # mu upper
  geom_point(aes(x = Z_STATION, y = mu_upper), 
             size = 1, 
             color = '#54AEFF') +
  # regression up
  geom_line(data = upper_line, aes(x = bdp, y =dpt),
            color='#54AEFF',
            size=0.5,
            alpha = 0.4) +
  #bottom depth line
  geom_abline(intercept = 0, slope = -1, 
              color="#FF5459", 
              linetype="dashed", 
              size=0.2) +
  # mu lower
  geom_point(aes(x = Z_STATION, y = mu_lower), 
             size = 1, 
             color = '#FF5459') +
  # regression low
  geom_line(data = lower_line, aes(x = bdp, y =dpt),
            color='#FF5459',
            size=0.5,
            alpha = 0.4) +
  # split lower line of mu upper
  geom_line(data = mid_sep_up, aes(x = bdp, y =dpt),
            color="#54AEFF", 
            linetype="dashed", 
            size=0.2) +
  # split upper line of mu lower
  geom_line(data = mid_sep_lo, aes(x = bdp, y =dpt),
            color="#FF5459", 
            linetype="dashed", 
            size=0.2) +
  # surface line
  geom_abline(intercept = 0, slope = 0, 
              color="#54AEFF", 
              linetype="dashed", 
              size=0.2) +
  scale_y_reverse() +
  xlim(0, 600) +
  ylim(600, 0) +
  ggtitle("position of mu upper and lower in the water column", 
          subtitle = paste0(tax, " ", sta, " diapause")) +
  xlab("bottom depth") + ylab("mu") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/mu_bdep_diapause", grp, ".jpeg"), 
       height = 5, width = 7)



# Regression out of diapause :
# get the 3 points x y of the corresation path :
changing_up <- (intercept_up + 10*(bdep_up)) / (1-(bdep_up))
upper_line <- data.frame(
  bdp = c(0, changing_up+10, 600), 
  dpt = c(0, changing_up, intercept_up + 600*(bdep_up))
)

changing_lo <- (intercept_lo + 10*(bdep_lo)) / (1-(bdep_lo))
lower_line <- data.frame(
  bdp = c(0, changing_lo+10, 600), 
  dpt = c(0, changing_lo, intercept_lo + 600*(bdep_lo))
)

ggplot(data = pred_params[which(pred_params$diapause == "NO"),]) +
  # mu upper
  geom_point(aes(x = Z_STATION, y = mu_upper), 
             size = 1, 
             color = '#54AEFF') +
  # regression up
  geom_line(data = upper_line, aes(x = bdp, y =dpt),
            color='#54AEFF',
            size=0.5,
            alpha = 0.4) +
  #bottom depth line
  geom_abline(intercept = 0, slope = -1, 
              color="#FF5459", 
              linetype="dashed", 
              size=0.2) +
  # mu lower
  geom_point(aes(x = Z_STATION, y = mu_lower), 
             size = 1, 
             color = '#FF5459') +
  # regression low
  geom_line(data = lower_line, aes(x = bdp, y =dpt),
            color='#FF5459',
            size=0.5,
            alpha = 0.4) +
  # split lower line of mu upper
  geom_line(data = mid_sep_up, aes(x = bdp, y =dpt),
            color="#54AEFF", 
            linetype="dashed", 
            size=0.2) +
  # split upper line of mu lower
  geom_line(data = mid_sep_lo, aes(x = bdp, y =dpt),
            color="#FF5459", 
            linetype="dashed", 
            size=0.2) +
  # surface line
  geom_abline(intercept = 0, slope = 0, 
              color="#54AEFF", 
              linetype="dashed", 
              size=0.2) +
  scale_y_reverse() +
  xlim(0, 600) +
  ylim(600, 0) +
  ggtitle("position of mu upper and lower in the water column", 
          subtitle = paste0(tax, " ", sta, " active")) +
  xlab("bottom depth") + ylab("mu") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/mu_bdep_active", grp, ".jpeg"), 
       height = 5, width = 7)



# Regression without diapause effect (in case no effect was the best model):
# get the 3 points x y of the corresation path :
changing_up <- (intercept_up + 10*bdep_up) / (1-bdep_up)
upper_line <- data.frame(
  bdp = c(0, changing_up+10, 600), 
  dpt = c(0, changing_up, intercept_up + 600*bdep_up)
)

changing_lo <- (intercept_lo + 10*bdep_lo) / (1-bdep_lo)
lower_line <- data.frame(
  bdp = c(0, changing_lo+10, 600), 
  dpt = c(0, changing_lo, intercept_lo + 600*bdep_lo)
)

ggplot(data = pred_params) +
  # mu upper
  geom_point(aes(x = Z_STATION, y = mu_upper), 
             size = 1, 
             color = '#54AEFF') +
  # regression up
  geom_line(data = upper_line, aes(x = bdp, y =dpt),
            color='#54AEFF',
            size=0.5,
            alpha = 0.4) +
  #bottom depth line
  geom_abline(intercept = 0, slope = -1, 
              color="#FF5459", 
              linetype="dashed", 
              size=0.2) +
  # mu lower
  geom_point(aes(x = Z_STATION, y = mu_lower), 
             size = 1, 
             color = '#FF5459') +
  # regression low
  geom_line(data = lower_line, aes(x = bdp, y =dpt),
            color='#FF5459',
            size=0.5,
            alpha = 0.4) +
  # split lower line of mu upper
  geom_line(data = mid_sep_up, aes(x = bdp, y =dpt),
            color="#54AEFF", 
            linetype="dashed", 
            size=0.2) +
  # split upper line of mu lower
  geom_line(data = mid_sep_lo, aes(x = bdp, y =dpt),
            color="#FF5459", 
            linetype="dashed", 
            size=0.2) +
  # surface line
  geom_abline(intercept = 0, slope = 0, 
              color="#54AEFF", 
              linetype="dashed", 
              size=0.2) +
  scale_y_reverse() +
  xlim(0, 600) +
  ylim(600, 0) +
  ggtitle("position of mu upper and lower in the water column", 
          subtitle = paste0(tax, " ", sta)) +
  xlab("bottom depth") + ylab("mu") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/mu_bdep", grp, ".jpeg"), 
       height = 5, width = 7)



#-- Params
mean(pred_params$mu_upper, na.rm = T)
bdep_up
bdep_up + YES_up
bdep_lo
bdep_lo + YES_lo
intercept_up
intercept_lo




#### model P ####

#-- Graphs 

# P in relation to bottom depth :
ggplot(data = pred_params) +
  geom_point(aes(x = Z_STATION, y = p), 
             size = 1, 
             color = 'black') +
  ggtitle("p linked to bottom depth", 
          subtitle = paste0(tax, " ", sta)) +
  xlab("bottom depth") + ylab("p") +
  theme_minimal() +
  theme(plot.title=element_text(size=10))

# P linked to diapause (without NA) :
ggplot(data = pred_params |> filter(!is.na(diapause)), aes(x=diapause, y = p)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color='black', size=0.7, alpha=0.8) +
  ggtitle("p and Diapause state", 
          subtitle = paste0(tax, ", stage ", tax)) +
  xlab("Diapause") + ylab("p") +
  theme_minimal()

ggsave(paste0("./results_krumh_optim/p_diapause", grp, ".jpeg"), 
       height = 7, width = 5)


# Value of P, and relation to diapause : 

# Value of P
mean(pred_params$p)                                        # without diapause effect
mean(pred_params$p[which(pred_params$diapause == "YES")])  # in diapause
mean(pred_params$p[which(pred_params$diapause == "NO" )])  # active

# statistical difference 
# (wilcoxon test because the variance differ)
wilcox.test(pred_params$p[which(pred_params$diapause == "YES")], 
            pred_params$p[which(pred_params$diapause == "NO")] )





#### model sigma ####

# Just to know the mean :
mean(pred_params$sigma_upper)
mean(pred_params$sigma_lower)


# I use weighted least squares regression, to integrate the increase of 
# variance with increasing bottom depth


#-- sigma upper

# Find data with NA :
# pos_sigma_up <- which(is.na(pred_params$sigma_upper) | is.na(pred_params$diapause))
# 
# # Define weights to use :
# reg_sigma_up <-  lm(sigma_upper[-pos_sigma_up] ~ 0 + Z_STATION[-pos_sigma_up], data = pred_params)
# wt_sigma_up <- 1 / lm(abs(reg_sigma_up$residuals) ~ reg_sigma_up$fitted.values)$fitted.values^2
# 
# # Regression with diapause effect
# WLM_sigma_up_D <-  lm(sigma_upper[-pos_sigma_up] ~ Z_STATION[-pos_sigma_up] + 
#                         diapause[-pos_sigma_up]:Z_STATION[-pos_sigma_up] , 
#                       data = pred_params, weights=wt_sigma_up)
# AIC(WLM_sigma_up_D)
# 
# # Regression without diapause effect
# WLM_sigma_up <-  lm(sigma_upper[-pos_sigma_up] ~ Z_STATION[-pos_sigma_up], 
#                     data = pred_params, weights=wt_sigma_up)
# AIC(WLM_sigma_up)
# 
# # keep the coefficients 
# if (round(AIC(WLM_sigma_up_D)) < round(AIC(WLM_sigma_up))) {  #if the diapause model fit better
#   coef_sig_up     <- WLM_sigma_up_D$coefficients[2]
#   YES_sig_up      <- WLM_sigma_up_D$coefficients[3]  #coef rectification if not in diapause
# } else {
#   coef_sig_up     <- WLM_sigma_up$coefficients[1]
#   YES_sig_up      <- 0
# }
# 


#-- sigma lower

# Find data with NA :
pos_sigma_lo <- which(is.na(pred_params$sigma_lower) | is.na(pred_params$diapause))

# Define weights to use :
reg_sigma_lo <-  lm(sigma_lower[-pos_sigma_lo] ~ 0 + Z_STATION[-pos_sigma_lo], data = pred_params)
wt_sigma_lo <- 1 / lm(abs(reg_sigma_lo$residuals) ~ reg_sigma_lo$fitted.values)$fitted.values^2

# Regression with diapause effect
WLM_sigma_lo_D <-  lm(sigma_lower[-pos_sigma_lo] ~ 0 + Z_STATION[-pos_sigma_lo] + 
                        diapause[-pos_sigma_lo]:Z_STATION[-pos_sigma_lo] , 
                      data = pred_params, weights=wt_sigma_lo)
AIC(WLM_sigma_lo_D)

# Regression without diapause effect
WLM_sigma_lo <-  lm(sigma_lower[-pos_sigma_lo] ~ 0 + Z_STATION[-pos_sigma_lo], 
                    data = pred_params, weights=wt_sigma_lo)
AIC(WLM_sigma_lo)

# keep the coefficients 
if (round(AIC(WLM_sigma_lo_D)) < round(AIC(WLM_sigma_lo))) {  #if the diapause model fit better
  coef_sig_lo     <- WLM_sigma_lo_D$coefficients[1]
  NO_sig_lo       <- WLM_sigma_lo_D$coefficients[2]  #coef rectification if not in diapause
} else {
  coef_sig_lo     <- WLM_sigma_lo$coefficients[1]
  NO_sig_lo       <- 0
}




#-- Graphs

# Sigma upper in diapause : 
ggplot(data = pred_params[which(pred_params$diapause == "YES"),]) +
  # sigma upper
  geom_point(aes(x = Z_STATION, y = sigma_upper), 
             size = 1, 
             color = '#54AEFF') +
  # regression
  # geom_abline(intercept = 0, slope = coef_sig_up,
  #             color='#54AEFF',
  #             size=0.5,
  #             alpha = 0.4) +
  # upper limit
  geom_abline(intercept = 1.5, slope = 0,
              color='black',
              size=0.3,
              alpha = 0.4) +
  xlim(0, 600) +
  ylim(0, 1.6) +
  ggtitle("position of sigma upper in the water column", 
          subtitle = paste0(tax, " ", sta, " diapause")) +
  xlab("bottom depth") + ylab("sigma") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/sigma_bdep_upper", grp, "_diapause", ".jpeg"), 
       height = 5, width = 7)


# Sigma upper active : 
ggplot(data = pred_params[which(pred_params$diapause == "NO"),]) +
  # sigma upper
  geom_point(aes(x = Z_STATION, y = sigma_upper), 
             size = 1, 
             color = '#54AEFF') +
  # regression
  # geom_abline(intercept = 0, slope = coef_sig_up + NO_sig_up,
  #             color='#54AEFF',
  #             size=0.5,
  #             alpha = 0.4) +
  # upper limit
  geom_abline(intercept = 1.5, slope = 0,
              color='black',
              size=0.3,
              alpha = 0.4) +
  xlim(0, 600) +
  ylim(0, 1.6) +
  ggtitle("position of sigma lower in the water column", 
          subtitle = paste0(tax, " ", sta, "_active")) +
  xlab("bottom depth") + ylab("sigma") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/sigma_bdep_upper", grp, "_active", ".jpeg"), 
       height = 5, width = 7)


# Sigma lower in diapause : 
ggplot(data = pred_params[which(pred_params$diapause == "YES"),]) +
  # sigma lower
  geom_point(aes(x = Z_STATION, y = sigma_lower), 
             size = 1, 
             color = '#FF5459') +
  # regression
  geom_abline(intercept = 0, slope = coef_sig_lo,
              color='#FF5459',
              size=0.5,
              alpha = 0.4) +
  # upper limit
  geom_abline(intercept = 0, slope = 0.01,
              color='black',
              size=0.3,
              alpha = 0.4) +
  xlim(0, 600) +
  ylim(0, (600*0.01)) +
  ggtitle("position of sigma lower in the water column", 
          subtitle = paste0(tax, " ", sta, " diapause")) +
  xlab("bottom depth") + ylab("sigma") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/sigma_bdep_lower", grp, "_diapause", ".jpeg"), 
       height = 5, width = 7)


# Sigma lower active : 
ggplot(data = pred_params[which(pred_params$diapause == "NO"),]) +
  # sigma lower
  geom_point(aes(x = Z_STATION, y = sigma_lower), 
             size = 1, 
             color = '#FF5459') +
  # regression
  geom_abline(intercept = 0, slope = coef_sig_lo + NO_sig_lo,
              color='#FF5459',
              size=0.5,
              alpha = 0.4) +
  # upper limit
  geom_abline(intercept = 0, slope = 0.01,
              color='black',
              size=0.3,
              alpha = 0.4) +
  xlim(0, 600) +
  ylim(0, (600*0.01)) +
  ggtitle("position of sigma lower in the water column", 
          subtitle = paste0(tax, " ", sta, " active")) +
  xlab("bottom depth") + ylab("sigma") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/sigma_bdep_lower", grp, " active", ".jpeg"), 
       height = 5, width = 7)


# Sigma upper without diapause effect (if no effect is the best model): 
ggplot(data = pred_params) +
  # sigma upper
  geom_point(aes(x = Z_STATION, y = sigma_upper), 
             size = 1, 
             color = '#54AEFF') +
  # regression
  # geom_abline(intercept = 0, slope = coef_sig_up,
  #             color='#54AEFF',
  #             size=0.5,
  #             alpha = 0.4) +
  # upper limit
  geom_abline(intercept = 1.5, slope = 0,
              color='black',
              size=0.3,
              alpha = 0.4) +
  xlim(0, 600) +
  ylim(0, 1.6) +
  ggtitle("position of sigma upper in the water column", 
          subtitle = paste0(tax, " ", sta)) +
  xlab("bottom depth") + ylab("sigma") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/sigma_bdep_upper", grp, ".jpeg"), 
       height = 5, width = 7)


# Sigma lower without diapause effect (if no effect is the best model): 
ggplot(data = pred_params) +
  # sigma lower
  geom_point(aes(x = Z_STATION, y = sigma_lower), 
             size = 1, 
             color = '#FF5459') +
  # regression
  geom_abline(intercept = 0, slope = coef_sig_lo,
              color='#FF5459',
              size=0.5,
              alpha = 0.4) +
  # upper limit
  geom_abline(intercept = 0, slope = 0.01,
              color='black',
              size=0.3,
              alpha = 0.4) +
  xlim(0, 600) +
  ylim(0, (600*0.01)) +
  ggtitle("position of sigma lower in the water column", 
          subtitle = paste0(tax, " ", sta)) +
  xlab("bottom depth") + ylab("sigma") +
  theme(aspect.ratio = 0.5,
        plot.title=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white")
  )

ggsave(paste0("./results_krumh_optim/sigma_bdep_lower", grp, ".jpeg"), 
       height = 5, width = 7)



#-- Params

mean(pred_params$sigma_upper)
mean(pred_params$sigma_upper[which(pred_params$diapause == 'YES')])
coef_sig_up
coef_sig_up + NO_sig_up       
coef_sig_lo
coef_sig_lo + NO_sig_lo




