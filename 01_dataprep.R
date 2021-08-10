###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Data preparation of PROVALID & GCKD
###################################################


rm(list=ls())

pacman::p_load(tidyr, plyr, reshape2, ggplot2, openxlsx, stringr, transplantr, skimr,
               lme4, readxl, purrr, janitor, tableone, Hmisc, dplyr)


# ------ Initialization 
set.seed(12345)

# Data and Paths
data.path = "../Data/"
GCKD.path = paste0(data.path, "GCKD/")
PROVALID.path = paste0(data.path, "PROVALID/")
out.path = "../Output/"



## ------------------- FUNCTIONS ----------

getProvalidFU <- function(sheet, path){
  df <- read_excel(path, sheet)
  df <- clean_names(df)
  rel.cols <- c("patient_reference_id", "date_yyyy_mm_dd_10","serum_creatinine_mg_dl", "serum_creatinine_mmol_l")
  df <- df[,rel.cols]
  df$source <- sheet
  df$Time <- ifelse(str_detect(sheet, "FOLLOW-UP"), str_match(sheet, "FOLLOW-UP\\s*(.*?)\\s*LOCAL LABORATORY")[,2], "0")
  df$FU_serumcrea <- ifelse(is.na(df$serum_creatinine_mg_dl), df$serum_creatinine_mmol_l/100, df$serum_creatinine_mg_dl)
  df <- df[!is.na(df$FU_serumcrea),]
  colnames(df)[1] <- "Patient reference [ID]"

  return(df)
}

calcBMI <- function(height_m, weight_kg){
  BMI <-(weight_kg/height_m^2)
  return(BMI)
}

calcUACR <- function(df){
  # df= data.provalid
  
  df <- df %>%
    mutate(urin_alb1st_mg_dl = ifelse(is.na(x1st_urinary_albumin_mg_l), x1st_urinary_albumin_g_l*100, x1st_urinary_albumin_mg_l/10),
           urin_alb2nd_mg_dl = ifelse(is.na(x2nd_urinary_albumin_mg_l), x2nd_urinary_albumin_g_l*100, x2nd_urinary_albumin_mg_l/10),
           urin_alb3rd_mg_dl = ifelse(is.na(x3rd_urinary_albumin_mg_l), x3rd_urinary_albumin_g_l*100, x3rd_urinary_albumin_mg_l/10),
           urin_crea1st_g_dl = ifelse(is.na(x1st_urinary_creatinine_mg_dl), x1st_urinary_creatinine_mmol_l*0.01132, x1st_urinary_creatinine_mg_dl/1000),
           urin_crea2nd_g_dl = ifelse(is.na(x2nd_urinary_creatinine_mg_dl), x2nd_urinary_creatinine_mmol_l*0.01132, x2nd_urinary_creatinine_mg_dl/1000),
           urin_crea3rd_g_dl = ifelse(is.na(x3rd_urinary_creatinine_mg_dl), x3rd_urinary_creatinine_mmol_l*0.01132, x3rd_urinary_creatinine_mg_dl/1000))
  
  urin_alb_mg_dl <- df %>% 
    dplyr::select(urin_alb1st_mg_dl, urin_alb2nd_mg_dl, urin_alb3rd_mg_dl) %>%
    mutate(urin_alb_mg_dl = apply(., 1, function(x) first(na.omit(x)))) %>%
    dplyr::select(urin_alb_mg_dl)
  
  
  urin_crea_g_dl <- df %>% 
    dplyr::select(urin_crea1st_g_dl, urin_crea2nd_g_dl, urin_crea3rd_g_dl) %>%
    mutate(urin_crea_g_dl = apply(., 1, function(x) first(na.omit(x)))) %>%
    dplyr::select(urin_crea_g_dl)
  
  UACR_mgg <- (urin_alb_mg_dl/urin_crea_g_dl)[[1]]
  
  # Set bl uacr to 0.01 to allow log2 transformation due to skewness
  UACR_mgg_adj <- ifelse(as.numeric(UACR_mgg) == 0, 0.05, as.numeric(UACR_mgg))

  return(UACR_mgg_adj)
}


# ------------- Get data 
## GCKD
data.gckd <- read_excel(paste0(GCKD.path, "Datenanfrage2_export_20210625.xlsx"))

## PROVALID
file = paste0(PROVALID.path, "PROVALID BASE Final Data Export_Clinical_LocalLab_Endpoints_20200119.xls")
sheets <- excel_sheets(file)
fu.sheets <- sheets[str_detect(sheets,"LABORATORY")]

# Select relevant columns
clin.cols <- c("ID", "Patient reference", "Patient reference [ID]", "Date of visit [yyyy/mm/dd]", "Year of birth", "Country","Gender","Smoking","Body weight [kg]","Height [cm]",
               "Diastolic [mm Hg]", "Systolic [mm Hg]")
lab.cols <- c("Patient reference [ID]","Serum cholesterol (total) [mg/dl]", "Serum cholesterol (total) [mmol/l]", 
              "Serum creatinine [mg/dl]", "Hemoglobin [g/l]", "Hemoglobin [g/dl]", "HbA1C [%]", "HbA1C [mmol/l]",
              "1st Urinary creatinine [mg/dl]", "1st Urinary creatinine [mmol/l]", "1st Urinary albumin [mg/l]","1st Urinary albumin [g/l]", 
              "2nd Urinary creatinine [mg/dl]", "2nd Urinary creatinine [mmol/l]", "2nd Urinary albumin [mg/l]","2nd Urinary albumin [g/l]",
              "3rd Urinary creatinine [mg/dl]", "3rd Urinary creatinine [mmol/l]", "3rd Urinary albumin [mg/l]","3rd Urinary albumin [g/l]")
gluc.meds <- c("Biguanides (metformin)", "Insulines","Sulfonylureas", "DPPIV inhibitors or GLP1 analogs", "Meglitinides (glinides)", "Thiazolinediones (glitazones)",
               "Alpha-Glucosidase-inhibitors","SGLT2 inhibitors (gliflozins)")
bp.meds <- c("ACE inhibitors", "Renin inhibitors","Angiotensin-II-receptor blockers","Beta-receptor blockers", "Calcium-antagonists", 
             "Centrally acting antihypertensives", "Alpha-receptor blockers","Direct vasodilators","Loop diuretics", "Thiazides",
             "Potassium saving diuretics", "Aldosterone antagonists")
lip.meds <- c("Clofibric acid derivative", "Statins", "Others (ezetimibe, omega 3 acid)")

data.tmp1 <- read_excel(file, sheet=1)[,c(clin.cols, gluc.meds, bp.meds, lip.meds)]
data.tmp1$BL_med_dm <- apply(data.tmp1[,gluc.meds],1, function(x) any(!is.na(x)))*1
data.tmp1$BL_med_bp <- apply(data.tmp1[,bp.meds],1, function(x) any(!is.na(x)))*1
data.tmp1$BL_med_lipid <- apply(data.tmp1[,lip.meds],1, function(x) any(!is.na(x)))*1

data.tmp2 <- read_excel(file, sheet=2)[,lab.cols]
data.tmp3 <- lapply(fu.sheets, function(x) getProvalidFU(x, path=file))
data.tmp3 <- do.call(rbind, data.tmp3)

data.provalid <- left_join(left_join(data.tmp1, data.tmp2, by="Patient reference [ID]"), 
                           data.tmp3[,c("Patient reference [ID]", "Time", "FU_serumcrea")], by="Patient reference [ID]")
rm(data.tmp1,data.tmp2,data.tmp3)



# ----- Data Pre-Processing
data.gckd$eth <- "non-black"
data.gckd <- data.gckd %>%
  mutate(FU2_age = as.numeric(BL_age + round((FU2_fu_visdat-BL_ein_visdat)/365.25,1)),
         FU3_age = as.numeric(BL_age + round((FU3_fu_visdat-BL_ein_visdat)/365.25,1)),
         FU4_age = as.numeric(BL_age + round((FU4_fu_visdat-BL_ein_visdat)/365.25,1)),
         FU5_age = as.numeric(BL_age + round((FU5_fu_visdat-BL_ein_visdat)/365.25,1)),
         FU6_age = as.numeric(BL_age + round((FU6_fu_visdat-BL_ein_visdat)/365.25,1))) %>%
  mutate(BL_eGFR_epi = ckd_epi(creat = BL_creavalue, age = BL_age, sex = BL_gender, eth = eth, units = "US"),
         FU2_eGFR_epi = ckd_epi(creat = FU2_creavalue, age = FU2_age, sex = BL_gender, eth = eth, units = "US"),
         FU3_eGFR_epi = ckd_epi(creat = FU3_creavalue, age = FU3_age, sex = BL_gender, eth = eth, units = "US"),
         FU4_eGFR_epi = ckd_epi(creat = FU4_creavalue, age = FU4_age, sex = BL_gender, eth = eth, units = "US"),
         FU5_eGFR_epi = ckd_epi(creat = FU5_creavalue, age = FU5_age, sex = BL_gender, eth = eth, units = "US"),
         FU6_eGFR_epi = ckd_epi(creat = FU6_creavalue, age = FU6_age, sex = BL_gender, eth = eth, units = "US"),
         BL_eGFR_mdrd = mdrd(creat = BL_creavalue, age = BL_age, sex = BL_gender, eth = eth, units = "US"),
         FU2_eGFR_mdrd = mdrd(creat = FU2_creavalue, age = FU2_age, sex = BL_gender, eth = eth, units = "US"),
         FU4_eGFR_mdrd = mdrd(creat = FU4_creavalue, age = FU4_age, sex = BL_gender, eth = eth, units = "US"),
         FU6_eGFR_mdrd = mdrd(creat = FU6_creavalue, age = FU6_age, sex = BL_gender, eth = eth, units = "US"),
         PatID=as.character(subjid),
         Country="GE",
         BL_smoking = ifelse(BL_smoking == 0, 0, 1)) %>%
  dplyr::rename(BL_bpsys = BL_bloodpr_sys,
         BL_bpdia = BL_bloodpr_dias,
         BL_hemo = BL_hemovalue,
         BL_hba1c = BL_phbavalue,
         BL_hba1c_perc = BL_hbavalue,
         BL_serumchol = BL_cholvalue1,
         BL_med_lipid = BL_med_lipidsenker) 

data.gckd.long <- data.gckd %>%
  dplyr::select(PatID, BL_eGFR_epi, FU2_eGFR_epi, FU3_eGFR_epi, FU4_eGFR_epi, FU5_eGFR_epi, FU6_eGFR_epi) %>%
  `colnames<-`(c("PatID","0","2","3","4","5","6")) %>%
  pivot_longer(cols=2:5, names_to = "Time", values_to = "FU_eGFR_epi") %>%
  mutate(Cohort="GCKD") %>%
  filter(!is.na(FU_eGFR_epi)) %>%
  mutate(Time=as.numeric(Time))

data.gckd <- right_join(data.gckd, data.gckd.long, by="PatID")


data.provalid <- data.provalid  %>%
  clean_names()
data.provalid <- data.provalid  %>%
  mutate(BL_age = round(as.numeric(as.Date(date_of_visit_yyyy_mm_dd)-as.Date(ISOdate(year_of_birth, 1, 1)))/365.25,1),
         BL_gender = ifelse(gender=="Female",1,0),
         BL_smoking = ifelse(smoking=="never",0,1),
         BL_bmi = calcBMI(height_m=height_cm/100, body_weight_kg),
         BL_diabetic=1,
         BL_serumchol = ifelse(is.na(serum_cholesterol_total_mg_dl), serum_cholesterol_total_mmol_l*38.67 , serum_cholesterol_total_mg_dl),
         BL_hemo = ifelse(is.na(hemoglobin_g_dl), hemoglobin_g_l/10, hemoglobin_g_dl),
         BL_uacr = calcUACR(data.provalid),
         BL_hba1c_perc = ifelse(is.na(hb_a1c_mmol_l), ((hb_a1c_percent-2.15) * 10.929), hb_a1c_mmol_l),
         BL_hba1c = ifelse(is.na(hb_a1c_percent), (hb_a1c_mmol_l/10.929)+2.15, hb_a1c_percent),
         eth = "non-black",
         FU_age = BL_age + as.numeric(data.provalid$time),
         FU_eGFR_epi = ckd_epi(creat = data.provalid$fu_serumcrea, age = FU_age, sex = BL_gender, eth = eth, units = "US"),
         BL_eGFR_epi = ckd_epi(creat = serum_creatinine_mg_dl, age = BL_age, sex=BL_gender, ethnicity = eth, units="US"),
         Cohort="PROVALID",
         PatID=as.character(patient_reference_id),
         Country= str_sub(data.provalid$country, -2,-1),
         Time=as.numeric(time)) %>%
  dplyr::rename(BL_bpsys = systolic_mm_hg,
                BL_bpdia = diastolic_mm_hg,
                BL_med_lipid = bl_med_lipid,
                BL_med_bp = bl_med_bp,
                BL_med_dm = bl_med_dm)
      

# ------------ Merge datasets
cols <- c("PatID", "Cohort","Country","Time","FU_eGFR_epi","BL_age", "BL_gender", "BL_smoking",
          "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc", "BL_hba1c", "BL_serumchol", "BL_hemo", "BL_uacr", "BL_med_dm", "BL_med_bp", "BL_med_lipid")

data.all <- rbind(data.provalid[,cols], data.gckd[,cols])

data.all <- data.all %>%
  mutate_at(c("Time","FU_eGFR_epi", "BL_age", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc", "BL_hba1c", "BL_serumchol", "BL_hemo", "BL_uacr"),as.numeric) %>%
  mutate_at(c("Cohort", "Country", "BL_gender", "BL_smoking", "BL_med_dm", "BL_med_bp", "BL_med_lipid"), as.factor)



# ------------ Inclusion & Exclusion criteria
data.all <- data.all %>%
  filter(!(Time==0 & FU_eGFR_epi<=30)) %>%
  mutate(BL_uacr_log2=log(BL_uacr,2),
         Cohort=ifelse(Cohort %in% "PROVALID", 1, 0)) %>%  # Code PROVALID = 1, GCKD = 0
  filter(BL_hemo > 2)                                     # Sorts out values in wrong unit column (g/l instead of g/dL)

rec.visits <- apply(as.data.frame.matrix(table(data.all$PatID, data.all$Time)),1, sum)
data.all <- data.all[data.all$PatID %in% names(which(rec.visits >=3)),]


  
