###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Initial data analysis: data screening & cleaning
###################################################

source("01_dataprep.R")

# --------- Missingness
data.full <- data.all[complete.cases(data.all),]
data.full$Country <- factor(data.full$Country, levels=c("AT", "GE", "HU", "PL", "UK", "NL"))
data.rem <- data.all[!complete.cases(data.all),]

length(unique(data.all$PatID))
length(unique(data.full$PatID))

table(data.all[data.all$Time==0,]$Country)
table(data.full[data.full$Time==0,]$Country)
table(data.rem[data.rem$Time==0,]$Country)

pred.vars <- c("Cohort","Time","BL_age", "BL_gender", "BL_smoking", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c", 
               "BL_serumchol", "BL_hemo", "BL_uacr_log2", "BL_med_dm", "BL_med_bp", "BL_med_lipid")

# ----- Descriptive statistics

# Number of subjects per country and cohort
table(data.full[!duplicated(data.full$PatID),]$Cohort)
table(data.full[!duplicated(data.full$PatID),]$Country)


# Distribution of baseline variables between PROVALID and GCKD
data.all %>% 
  filter(Time==0 & Cohort %in% "PROVALID") %>%
  skim()

d <- Hmisc::describe(data.all)

data.all %>% 
  filter(Time==0 & Cohort %in% "GCKD") %>%
  skim()

CreateTableOne(vars=c(colnames(data.full)[str_detect(colnames(data.full), "BL_")], "FU_eGFR_epi"), strata="Cohort", data=data.full[data.full$Time==0,], test=F)

# Histogram of all baseline variables
par(mfrow = c(3,5))
hist(data.full[data.full$Time==0,]$BL_age, col="grey", xlab=NULL, main = "Age [years]")
barplot(table(data.full[data.full$Time==0,]$BL_gender), ylab = "Frequency", main = "Gender [1=female]")
hist(data.full[data.full$Time==0,]$BL_bmi, col="grey", xlab=NULL, main = "BMI [kg/m]")
barplot(table(data.full[data.full$Time==0,]$BL_smoking), ylab = "Frequency", main = "Smoking [1=ever]")
hist(data.full[data.full$Time==0,]$BL_bpsys, col="grey", xlab=NULL, main = "Systolic blood pressure [mmHg]")
hist(data.full[data.full$Time==0,]$BL_bpdia, col="grey", xlab=NULL, main = "Diastolic blood pressure [mmHg]")
hist(data.full[data.full$Time==0,]$BL_serumchol, col="grey", xlab=NULL, main = "Serum Cholesterin")
hist(data.full[data.full$Time==0,]$BL_hemo, col="grey", xlab=NULL, main = "Hemoglobin [g/dL]")
hist(data.full[data.full$Time==0,]$BL_hba1c, col="grey", xlab=NULL, main = "HbA1c [mmol/mol]")
hist(data.full[data.full$Time==0,]$BL_hba1c_perc, col="grey", xlab=NULL, main = "HbA1c [%]")
hist(data.full[data.full$Time==0,]$BL_uacr, col="grey", xlab=NULL, main = "UACR")
hist(data.full[data.full$Time==0,]$BL_uacr_log2, col="grey", xlab=NULL, main = "Log2 UACR")
barplot(table(data.full[data.full$Time==0,]$BL_med_bp), ylab = "Frequency", main = "Blood pressure-lowering meds [1=yes]")
barplot(table(data.full[data.full$Time==0,]$BL_med_dm), ylab = "Frequency", main = "Glucose-lowering meds [1=yes]")
barplot(table(data.full[data.full$Time==0,]$BL_med_lipid), ylab = "Frequency", main = "Lipid-lowering meds [1=yes]")


# Histograms per time point for eGFR
par(mfrow = c(2,4))
data.full %>% filter(Time==0) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=0", xlab=NULL) 
data.full %>% filter(Time==1) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=1", xlab=NULL) 
data.full %>% filter(Time==2) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=2", xlab=NULL) 
data.full %>% filter(Time==3) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=3", xlab=NULL) 
data.full %>% filter(Time==4) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=4", xlab=NULL) 
data.full %>% filter(Time==5) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=5", xlab=NULL) 
data.full %>% filter(Time==6) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=6", xlab=NULL) 
par(mfrow = c(1,1))


# ---------- Correlation matrix
x_tmp <- model.matrix(FU_eGFR_epi~BL_age + BL_gender + BL_smoking + BL_bmi + BL_bpsys + BL_bpdia + BL_hba1c 
                      + BL_serumchol + BL_hemo + BL_uacr + BL_med_dm + BL_med_bp + BL_med_lipid, 
                      data=data.full[data.full$Time==0,])[,-1]
corr <- abs(cor(x_tmp))
(corr>0.5)


# ------ Bivariate distribution
plot_fun <- function( x, y, time=0 ){
  df <- data.frame(dep=as.numeric(unlist(data.full[data.full$Time==time, x])),indep=as.numeric(unlist(data.full[data.full$Time==time, y])))
  ggplot(df, aes(x=dep, y=indep)) + 
    geom_point() +
    scale_x_continuous(name = x)+
    scale_y_continuous(name = y) +
    theme_bw()
}

cont.vars <- c("BL_age", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c", 
               "BL_serumchol", "BL_hemo", "BL_uacr_log2")
pairw.int <- expand.grid(cont.vars, cont.vars)
pairw.int <- pairw.int[!pairw.int[,1]==pairw.int[,2],]
plots <- lapply(1:nrow(pairw.int), function(x) plot_fun(x=pairw.int[x,1], y=pairw.int[x,2], time=0))
plots


# ---- Check non-linearities
plots <- lapply(cont.vars, function(x) plot_fun(x=x, y="FU_eGFR_epi", time=1))
plots

# Non-linearity for BL_age
plots[[3]]

