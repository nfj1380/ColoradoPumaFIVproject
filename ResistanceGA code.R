###################################################
#                                                 #
#   Resistance GA code    #
# Adaptation by Nick Fountain-Jones           #
#                                                 #
###################################################

tinytex::install_tinytex()

devtools::install_github("wpeterman/ResistanceGA", build_vignettes = F) # Download package. Can be a slow process
options(tinytex.verbose = TRUE)

library(ResistanceGA) # Installs package and the other required packages needed

remotes::install_github("wpeterman/ResistanceGA") # when installing need old versions of of sp(1.2-5)/raster(2.8.19)

 
library(ResistanceGA)