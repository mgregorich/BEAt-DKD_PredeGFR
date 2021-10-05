###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Linear mixed effects model for the prediction of eGFR (BEAt-DKD)
###################################################


rm(list=ls())

pacman::p_load(tidyr, plyr, reshape2, ggplot2, openxlsx, stringr, transplantr, skimr,
               lme4, readxl, purrr, janitor, tableone, dplyr, splitstackshape,
               nlme, lmerTest, MuMIn, JMbayes, splines, rms, Hmisc, concreg, caret, MASS, performance,
               pmsampsize)
 

# ------ Initialization 
set.seed(12345)

# Data and Paths
data.path = "../Data/"
GCKD.path = paste0(data.path, "GCKD/")
PROVALID.path = paste0(data.path, "PROVALID/")
DIACORE.path = paste0(data.path, "DIACORE/")

out.path = "../Output/"

slope_cutpoint=-3

# Load auxiliary functions 
source("functions_aux.R")


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

