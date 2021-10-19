###################################################
# Date: 30/06/2021
# Info: Initialisation and global values
###################################################


pacman::p_load(tidyr, plyr, reshape2, ggplot2, openxlsx, stringr, transplantr, skimr,
               lme4, readxl, purrr, janitor, tableone, dplyr, splitstackshape,
               nlme, lmerTest, MuMIn, JMbayes, splines, rms, Hmisc, concreg, caret, MASS, performance,
               pmsampsize, broom)


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
source("scr/functions_aux.R")