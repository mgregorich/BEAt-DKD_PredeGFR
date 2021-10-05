###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Linear mixed effects model for the prediction of eGFR (BEAt-DKD)
###################################################

rm(list=ls())

# Initialise global values

source("setup.R")


# ------- Start the fun
# Data preparation
source("01_dataprep.R", print.eval=F)

# Initial data analysis
source("02_IDA.R", print.eval=F)

# Model building and validation
source("03_model.R", print.eval=F)

# Model evaluation (figures, tables)
source("04_eval.R", print.eval=F)

# Shiny web implementation
# source("05_extval.R", print.eval=F)

# Shiny web implementation
source("06_shiny.R", print.eval=F)

