# Nouveau code pour faire le graphe de variance 

library("dplyr")
library("ggplot2")
library(terra)

setwd("~/projects/def-frmap5/MINDZ_andeol/MINDZ/output_analyses")

folder_path <- "./results/to_compute" # location of the input files

agreg_values = c(1, 2, 3, 4, 5, 10, 15)

output_variances = NULL
loop = 0
for (i in agreg_values) {
  loop = loop + 1
  print(paste0('loop ', loop, ' of ', length(agreg_values)))
  
  # Lister les fichiers commençant par "grid_3d_LF_3"
  fichiers <- list.files(folder_path, 
                         pattern = paste0("^grid_3d_LF_", i, "_"), 
                         full.names = TRUE)
  
  
  # Lire les rasters et les stocker dans une liste
  rasters <- lapply(fichiers, rast)
  
  # Combiner les rasters en un seul objet SpatRaster avec plusieurs couches
  combined_rasters <- do.call(c, rasters)
  
  # Normaliser les valeurs de chaque couche individuellement
  for (j in 1:nlyr(combined_rasters)) {
    layer <- combined_rasters[[j]]
    max_val <- global(layer, fun=max, na.rm=TRUE)[[1]]
    combined_rasters[[j]] <- layer / max_val
  }
  
  # Calculer la variance pixel par pixel
  variance_raster <- app(combined_rasters, var)
  
  # Calculer la variance totale sur tous les pixels
  total_variance <- global(variance_raster, var, na.rm = T)
  
  # add variance value to the data frame :
  output_variances = rbind(output_variances, 
                           tibble(density  = i, 
                                  variance =as.numeric(total_variance)))
}



ggplot() +
  geom_point(data = output_variances, aes(x = density, y = variance))





