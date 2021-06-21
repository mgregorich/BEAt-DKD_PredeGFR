rm(list=ls())

pacman::p_load(tidyr, dplyr, reshape2, ggplot2, openxlsx, stringr, transplantr, skimr,
               lme4)


# ------ Initialization 
# Data and Paths
data.path = "../Data/"
GCKD.path = paste0(data.path, "GCKD/")
PROVALID.path = paste0(data.path, "PROVALID/")

out.path = "../Output/"


# ------ Get data 
data.gckd <- read.xlsx(paste0(GCKD.path, "20200203_GCKD_BEAtDKD_data_EOS.xlsx"), sheet=2)
data.provalid <- as.data.frame(readRDS(paste0(PROVALID.path, "PROVALID_20170921.rds")))


# ----- Data cleaning 
data.gckd$BL_creavalue[which(data.gckd$BL_creavalue<0)] <- NA
data.gckd$BL_ein_visdat <- convertToDate(data.gckd$BL_ein_visdat)
data.gckd$FU6_fu_visdat <- convertToDate(data.gckd$FU6_fu_visdat)
years <- (data.gckd$FU6_fu_visdat-data.gckd$BL_ein_visdat)/365.25
data.gckd$FU6_age <- as.numeric(data.gckd$BL_age + round(years,1))

# ----- Data Pre-Processing
data.gckd$eth <- "non-black"
data.gckd <- data.gckd %>%
  mutate(BL_eGFR_epi = ckd_epi(creat = BL_creavalue, age = BL_age, sex = BL_gender, eth = eth, units = "US"),
         FU2_eGFR_epi = ckd_epi(creat = FU2_creavalue, age = FU2_fu2_age, sex = BL_gender, eth = eth, units = "US"),
         FU4_eGFR_epi = ckd_epi(creat = FU4_creavalue, age = FU4_fu4_age, sex = BL_gender, eth = eth, units = "US"),
         FU6_eGFR_epi = ckd_epi(creat = FU6_creavalue, age = FU6_age, sex = BL_gender, eth = eth, units = "US"),
         BL_eGFR_mdrd = mdrd(creat = BL_creavalue, age = BL_age, sex = BL_gender, eth = eth, units = "US"),
         FU2_eGFR_mdrd = mdrd(creat = FU2_creavalue, age = FU2_fu2_age, sex = BL_gender, eth = eth, units = "US"),
         FU4_eGFR_mdrd = mdrd(creat = FU4_creavalue, age = FU4_fu4_age, sex = BL_gender, eth = eth, units = "US"),
         FU6_eGFR_mdrd = mdrd(creat = FU6_creavalue, age = FU6_age, sex = BL_gender, eth = eth, units = "US"))
  

data.gckd %>%
  select(subjid, BL_eGFR_epi, BL_eGFR_mdrd, FU2_eGFR_epi, FU2_eGFR_mdrd, FU4_eGFR_epi, FU4_eGFR_mdrd, FU6_eGFR_epi, FU6_eGFR_mdrd) %>%
  skim()




# ------ Merge datasets
# DO WHEN BL DATA OF GCKD IS AVAILABLE


# ------ Mixed effects model for slope estimation (adjusted by study cohort)

# GCKD + PROVAILD
data.gckd.long <- data.gckd %>%
  select(subjid, BL_eGFR_epi, FU2_eGFR_epi, FU4_eGFR_epi,  FU6_eGFR_epi) %>%
  `colnames<-`(c("ID","0","1","2","3")) %>%
  pivot_longer(cols=2:5, names_to = "Time", values_to = "eGFR_epi") %>%
  mutate(Cohort=1)

data.provalid <- data.provalid %>%
  mutate(Time=round(as.numeric(data.provalid$Time_FU)/365.25,0),
         Cohort=2) %>%
  rename(ID=PID,
         eGFR_epi = eGFR_EPI,
         eGFR_mdrd = eGFR_MDRD)

data.eGFR <- rbind(data.gckd.long[,c("ID","Cohort","Time","eGFR_epi")],data.provalid[,c("ID","Cohort","Time","eGFR_epi")] )
data.eGFR$ID <- as.factor(data.eGFR$ID)
data.eGFR$Time <- as.numeric(data.eGFR$Time)
data.eGFR$eGFR_epi <- as.numeric(data.eGFR$eGFR_epi)

lme.eGFR <- lmer(eGFR_epi ~ Time + Cohort + (Time|ID), data=data.eGFR)
eGFR.slopes <- data.frame("ID"=cbind(rownames(coef(lme.eGFR)$ID), "Intercept"=coef(lme.eGFR)$ID[,2]), "Slope"=coef(lme.eGFR)$ID[,2])

table(cut(eGFR.slopes$Slope, quantile(eGFR.slopes$Slope, seq(0,1,0.33))))


