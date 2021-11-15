###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Linear mixed effects model for the prediction of eGFR (BEAt-DKD)
###################################################

rm(list=ls())

# Initialise global values

source("scr/setup.R")


# ------- Start the fun
# Data preparation
source("scr/01_dataprep.R", print.eval=F)

# Initial data analysis
source("scr/02_IDA.R", print.eval=F)

# Model building and validation
source("scr/03_model.R", print.eval=F)

# Model evaluation (figures, tables)
source("scr/04_eval.R", print.eval=F)

# External validation
source("scr/05_extval.R", print.eval=F)

# Run shiny app
source("scr/06_shinyapp.R")


