###################################################
# Date: 30/06/2021
# Info: Initialisation and global values
###################################################


pacman::p_load(tidyr, plyr, reshape2, ggplot2, openxlsx, stringr, transplantr, skimr,shiny,
               lme4, readxl, purrr, janitor, tableone, dplyr, splitstackshape,r2glmm,
               nlme, lmerTest, MuMIn, splines, rms, Hmisc, concreg, caret, MASS, performance,
               pmsampsize, broom, broom.mixed, gridExtra, future.apply, parallel, future,
               sjmisc, sjPlot, forcats, pmsampsize, jsonlite, 
               shinyvalidate, shinyalert)


# ------ Initialization 
set.seed(12345)

# Data and Paths
data.path = "../Data/"
GCKD.path = paste0(data.path, "GCKD/")
PROVALID.path = paste0(data.path, "PROVALID/")
DIACORE.path = paste0(data.path, "DIACORE/")

out.path = "../Output/"
shiny.path = "./shiny"
slope_cutpoint=-3
nboot = 1000

# Load auxiliary functions 
source("scr/functions_aux.R")

# Set of predictors
pred.vars <- c("BL_age","BL_sex","BL_bmi", "BL_smoking", "BL_map",
               "BL_hba1c", "BL_serumchol", "BL_hemo", "BL_uacr_log2",
               "BL_med_dm", "BL_med_bp", "BL_med_lipid")
