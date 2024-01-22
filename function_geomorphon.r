###############################################################################
############ Traitement géomorphone sur un / des rasters ######################
###############################################################################
############################### Wojciech MEDYK ################################
################################ LP Pro SIGDAT ################################
###############################################################################
### Script basé sur : https://rzine.fr/docs/20230425_geomorphon/index.html#####
###############################################################################

################################## Packages ###################################
# Noms de packages nécessaires
my_packages <- c("sf", "terra", "vegan", "tidyverse", "qgisprocess", "ggplot2", "gtools")

# Vérifier si ces packages sont installés
missing_packages <- my_packages[!(my_packages %in% installed.packages()[,"Package"])]

# Installation des packages manquants depuis le CRAN
if(length(missing_packages)) install.packages(missing_packages, 
                                              repos = "http://cran.us.r-project.org")

# Chargement des packages nécessaires
lapply(my_packages, library, character.only = TRUE)

####### Configuration QGIS (nécessaire pour fonctionnement du script) ##########

qgis_configure()
qgis_enable_plugins()

################## Function création raster pour traitement ####################
## Exemple path_rasters = "/data/rast/"
## Exemple nom_raster = "nomachoisir"
## Exemple limites (!limites en format vecteur!) = "data/ariege.gpkg"

creation_raster <- function(path_rasters, nom_raster, limites){

  ## Charger les rasters 
    path_rasters <- paste(getwd(),toString(path_rasters), sep="")
    path_final <- c("data/", toString(nom_raster),".vrt")
    path_final <- paste(path_final, collapse = "") 
    ras_lst <- list.files(path_rasters, full.names =TRUE, pattern = ".asc$", recursive=TRUE)
    ras_lst <- qgis_list_input(!!! ras_lst)
  ## Assembler les rasters
    qgis_run_algorithm("gdal:buildvirtualraster",
                   INPUT=ras_lst,
                   SEPARATE=FALSE,
                   ASSIGN_CRS='EPSG:2154',
                   OUTPUT=path_final)
    path_final <- rast(path_final)
  
  ## Limites zone études 
    limites <- vect(toString(limites))
  
  ## Découpage de raster 
    zone_rast <- crop(path_final, limites)
    zone_rast <- mask(path_final, limites)

  ## Nom de raster pour l'export
    nom_raster <- c("data/", toString(nom_raster), ".tif")
    nom_raster <- paste(nom_raster, collapse = "")

  ## Export raster
    terra::writeRaster(zone_rast,
                   nom_raster,
                   filetype="GTiff", 
                   overwrite=TRUE)
}
######################### Traitement Geomorphon ################################
## Exemple zone_rast (résultat de création de raster) = "/data/nom_raster.tif
## Exemple limites (!limites en format vecteur!) = "data/ariege.gpkg"
## Exemple nom_ = "nomachoisir"
## Exemple rayon_recherche = 100 (taille de la focale en pixels)
## Exemple res = 5000 (taille de la maille en m)

traitement_geomorphon <- function(zone_rast, limites, nom_, rayon_recherche, res){

  ## Création de géomorphon par rayon de recherche
    for (rayon in rayon_recherche) {
      geom_ <- c("data/geomorphons_",toString(nom_),"_rayon_",toString(rayon*100/4),"m.tif")
      geom_ <- paste(geom_, collapse = "")
      qgis_run_algorithm("grass7:r.geomorphon",
                     elevation = zone_rast,
                     search = rayon, #rayon de recherche
                     flat = 1,
                     forms = geom_)
  ## Création de grille sur plusieurs resolutions 
    for (maille in res) {
      limites <- vect(toString(limites))
      grid <- c("data/",toString(nom_),"_maille_", toString(maille), ".gpkg")
      grid <- paste(grid, collapse = "")
      creation_maille <- st_make_grid(limites, maille)
      creation_maille <- vect(creation_maille)
      grid_final <- crop(creation_maille, limites)
      terra::writeVector(grid_final, 
                       grid,
                       filetype = "GPKG",
                       overwrite = T)
    ## Calculation de nombre de pixel
      path_nb_pixels <- c("data/",toString(nom_),"_nb_pixels.gpkg")
      path_nb_pixels <- paste(path_nb_pixels, collapse = "")
    
      qgis_run_algorithm("native:zonalhistogram",
                       INPUT_RASTER = geom_,
                       RASTER_BAND = qgisprocess:::qgis_default_value(),
                       INPUT_VECTOR = grid_final,
                       COLUMN_PREFIX = "nb_pixels_",
                       OUTPUT = path_nb_pixels,
    )
    ## Calculation de surface de la maille
      rast_final <- vect(path_nb_pixels)
      df <- data.frame(rast_final)
      df <-df %>% add_column(surface_maille = terra::expanse(rast_final))
    ## Calculation de surface de pixel 
      df2 <- df
      for(i in 1:ncol(df)) {   
        df2[ , i] <- df[ , i] * 25 / df$surface_maille
        names(df2)[i] <- paste("surface_pixels_", toString(i), sep ="")
      }
      df2 <- df2[,-ncol(df2)]
      df3 <- data.frame(c(df, df2))
    ## Indice de Shannon par maille
      for (j in 1:nrow(df3)){
        df3$shannon[j]<- diversity(df3[j,1:10])
      }
    ## Class indice Shannon
      df3$shannon_class <- quantcut(df3$shannon, 5) 
      levels(df3$shannon_class) <- c("1","2","3", "4", "5")
    ## Indice de Piélou
      S <- rowSums(df3[,1:10])
      S <- data.frame(S)
      df3$pielou <- (df3$shannon/base::log(S))
      df3$pielou <- as.numeric(unlist(df3$pielou))
    ## Class indice de Piélou
      df3$pielou_class <- quantcut(df3$pielou, 5)
      levels(df3$pielou_class) <- c("1","2","3", "4","5") 
    ## Sauvegarde du résultat final
      resultat_final <- c("data/geomorphon_", toString(nom_), "_rayon_", toString(rayon), "_zone_", toString(maille),"m.csv")
      final <- paste(resultat_final, collapse = "")
      write_csv(df3, final)
    }
  }
}
