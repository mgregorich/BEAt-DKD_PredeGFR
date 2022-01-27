#################################
# Author:MG
# Date: 07/10/2021
# Info: app file - shiny
################################

# Create Shiny app ----
pacman::p_load(shiny, shinyjs, shinythemes, nlme, ggplot2, reshape2, dplyr, tidyverse, png)
source("functions_shiny.R")

shinyApp(ui = ui, server = server)