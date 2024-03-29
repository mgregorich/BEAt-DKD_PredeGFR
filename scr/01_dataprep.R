#===============================================================================#
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Data preparation of PROVALID & GCKD
#===============================================================================#



# ============= TRIPOD flowchart setup=============================
tripod_flowchart <- data.frame("Step"=c("all", 
                                        "lack_var", 
                                        "excl2_baseunder30",  
                                        "excl3_nobase", 
                                        "excl1_not3eGFR",
                                        "post_excl",
                                        "completecases",
                                        "DevCohort"), 
                        "GCKD"=rep(NA,8), 
                        "PROVALID"=rep(NA,8),
                        "DIACORE"=rep(NA,8))



# ================== Data Pre-processing =======================================

## -----  DIACORE -----
recodeDates <- function(x){
  x=as.character(x)
  x=str_replace(x, "Mrz", "März")
  x.new=if_else(str_detect(x, "^[:digit:]+$"), as.character(as.Date(as.numeric(x), origin = "1899-12-30")), as.character(as.Date(x, "%d%b%Y")))
  x.new = as.Date(x.new)
  return(x.new)
}

data.diacore <-  read_excel(paste0(DIACORE.path, "BeatDKD_DIACORE_V1_V2R_V3R_13Mar2019_korr_13Aug2019_all_variables_nopw.xlsx"))
data.diacore$visitDate_char_V1=do.call("c",lapply(data.diacore$visitDate_char_V1, recodeDates))
data.diacore$visitDate_char_V2=do.call("c",lapply(data.diacore$visitDate_char_V2, recodeDates))
data.diacore$visitDate_char_V3=do.call("c",lapply(data.diacore$visitDate_char_V3, recodeDates))

data.diacore <- data.diacore %>%
  dplyr::select(patNr, AGE_V1, genderCode, smoke_ever_V1, BMI_V1, RRsys_mean_V1, RRdia_mean_V1,UACR_V1, HbA1c_percent_V1,  CHOL_V1, Hb_V1,
                Med_Lipid_V1, Med_BPlowering_V1, Med_Antidiabetika_V1,
                visitDate_char_V1, eGFR_CKDEPI_V1, visitDate_char_V2, eGFR_CKDEPI_V2, visitDate_char_V3, eGFR_CKDEPI_V3) %>%
  mutate(BL_map = (RRsys_mean_V1+ 2*RRdia_mean_V1)/3,
         Country="Unknown",
         Cohort=3,
         genderCode = as.factor(genderCode-1),
         smoke_ever_V1 = as.factor(ifelse(smoke_ever_V1==1, 0, 1)),
         AGE_V1 = round(AGE_V1, 0),
         patNr = as.character(patNr),
         Hb_V1 = as.numeric(sub(",", ".", Hb_V1, fixed = TRUE))) %>%
  relocate(c(Cohort, Country), .before = AGE_V1,) %>%
  relocate(BL_map, .before=RRsys_mean_V1) 

data.diacore <- merged.stack(data.diacore, id.vars="patNr", var.stubs=c("eGFR_CKDEPI"), sep = "_V") 
data.diacore$.time_1 <- as.numeric(data.diacore$.time_1)-1
colnames(data.diacore) <- c("PatID", "Time_cat","FU_eGFR_epi_2021", "Cohort", "Country","BL_age", "BL_sex", "BL_smoking","BL_bmi",
                            "BL_map", "BL_bpsys", "BL_bpdia","BL_uacr", "BL_hba1c_perc", "BL_serumchol","BL_hemo", 
                            "BL_med_lipid", "BL_med_bp", "BL_med_dm", "FU0_date", "FU1_date", "FU2_date")
data.diacore <- data.diacore %>%
  mutate(BL_uacr_log2 =log(data.diacore$BL_uacr, 2),
         BL_hba1c = (BL_hba1c_perc-2.15) * 10.929) %>%
  mutate_at(vars(BL_med_bp, BL_med_dm, BL_med_lipid), ~as.factor(.)) %>%
  filter(!(is.na(FU_eGFR_epi_2021) & (Time_cat %in% c(1,2)))) %>%
  mutate(Time_cont=0)

data.diacore$Time_cont=ifelse((data.diacore$Time_cat==1), round(difftime(data.diacore$FU1_date,data.diacore$FU0_date,units="days")/365.25,3), data.diacore$Time_cont)
data.diacore$Time_cont=ifelse((data.diacore$Time_cat==2), round(difftime(data.diacore$FU2_date,data.diacore$FU0_date,units="days")/365.25,3), data.diacore$Time_cont)
data.diacore$Time_cat= round(data.diacore$Time_cont)

data.diacore <- data.diacore %>%
  dplyr::select(PatID, Time_cat, Time_cont, FU_eGFR_epi_2021, Cohort, Country, BL_age, BL_sex, BL_smoking, BL_bmi, BL_map, BL_bpsys,     
               BL_bpdia, BL_uacr, BL_hba1c_perc, BL_serumchol, BL_hemo, BL_med_lipid, BL_med_bp, BL_med_dm, BL_uacr_log2, BL_hba1c)


## ----- GCKD ------
data.gckd <- read_excel(paste0(GCKD.path, "Datenanfrage2_export_20210625.xlsx"))

data.gckd$eth <- "non-black"
data.gckd <- data.gckd %>%
  mutate(FU2_age = as.numeric(BL_age + round((FU2_fu_visdat-BL_ein_visdat)/365.25,0)),
         FU3_age = as.numeric(BL_age + round((FU3_fu_visdat-BL_ein_visdat)/365.25,0)),
         FU4_age = as.numeric(BL_age + round((FU4_fu_visdat-BL_ein_visdat)/365.25,0)),
         FU5_age = as.numeric(BL_age + round((FU5_fu_visdat-BL_ein_visdat)/365.25,0)),
         FU6_age = as.numeric(BL_age + round((FU6_fu_visdat-BL_ein_visdat)/365.25,0)),
         BL_sex = BL_gender,
         BL_gender = ifelse(BL_gender==1, "F", "M")) %>%
  mutate(BL_eGFR_epi_2009 = ckd_epi(creat = BL_creavalue, age = BL_age, sex = BL_gender, eth = eth, units = "US"),
         FU2_eGFR_epi_2009 = ckd_epi(creat = FU2_creavalue, age = FU2_age, sex = BL_gender, eth = eth, units = "US"),
         FU3_eGFR_epi_2009 = ckd_epi(creat = FU3_creavalue, age = FU3_age, sex = BL_gender, eth = eth, units = "US"),
         FU4_eGFR_epi_2009 = ckd_epi(creat = FU4_creavalue, age = FU4_age, sex = BL_gender, eth = eth, units = "US"),
         FU5_eGFR_epi_2009 = ckd_epi(creat = FU5_creavalue, age = FU5_age, sex = BL_gender, eth = eth, units = "US"),
         FU6_eGFR_epi_2009 = ckd_epi(creat = FU6_creavalue, age = FU6_age, sex = BL_gender, eth = eth, units = "US"),
         BL_eGFR_epi_2021 = ckd_epi_2021(creat = BL_creavalue, age = BL_age, sex = BL_gender, units = "US"),
         FU2_eGFR_epi_2021 = ckd_epi_2021(creat = FU2_creavalue, age = FU2_age, sex = BL_gender, units = "US"),
         FU3_eGFR_epi_2021 = ckd_epi_2021(creat = FU3_creavalue, age = FU3_age, sex = BL_gender, units = "US"),
         FU4_eGFR_epi_2021 = ckd_epi_2021(creat = FU4_creavalue, age = FU4_age, sex = BL_gender, units = "US"),
         FU5_eGFR_epi_2021 = ckd_epi_2021(creat = FU5_creavalue, age = FU5_age, sex = BL_gender, units = "US"),
         FU6_eGFR_epi_2021 = ckd_epi_2021(creat = FU6_creavalue, age = FU6_age, sex = BL_gender, units = "US"),
         BL_eGFR_mdrd = mdrd(creat = BL_creavalue, age = BL_age, sex = BL_gender, eth = eth, units = "US"),
         FU2_eGFR_mdrd = mdrd(creat = FU2_creavalue, age = FU2_age, sex = BL_gender, eth = eth, units = "US"),
         FU4_eGFR_mdrd = mdrd(creat = FU4_creavalue, age = FU4_age, sex = BL_gender, eth = eth, units = "US"),
         FU6_eGFR_mdrd = mdrd(creat = FU6_creavalue, age = FU6_age, sex = BL_gender, eth = eth, units = "US"),
         PatID=as.character(subjid),
         Country="GE",
         BL_smoking = ifelse(BL_smoking == 0, 0, 1),
         BL_crea = BL_creavalue) %>%
  dplyr::rename(BL_bpsys = BL_bloodpr_sys,
                BL_bpdia = BL_bloodpr_dias,
                BL_hemo = BL_hemovalue,
                BL_hba1c = BL_hbavalue,
                BL_hba1c_perc = BL_phbavalue,
                BL_serumchol = BL_cholvalue1,
                BL_med_lipid = BL_med_lipidsenker) 

data.gckd.long.2021 <- data.gckd %>%
  dplyr::select(PatID, BL_eGFR_epi_2021, FU2_eGFR_epi_2021, FU3_eGFR_epi_2021, FU4_eGFR_epi_2021, FU5_eGFR_epi_2021, FU6_eGFR_epi_2021) %>%
  `colnames<-`(c("PatID","0","2","3","4","5","6")) %>%
  pivot_longer(cols=2:7, names_to = "Time_cat", values_to = "FU_eGFR_epi_2021") %>%
  mutate(Cohort="GCKD") %>%
  filter(!is.na(FU_eGFR_epi_2021)) %>%
  mutate(Time_cat=as.numeric(Time_cat))

data.gckd.long.2009 <- data.gckd %>%
  dplyr::select(PatID, BL_eGFR_epi_2009, FU2_eGFR_epi_2009, FU3_eGFR_epi_2009, FU4_eGFR_epi_2009, FU5_eGFR_epi_2009, FU6_eGFR_epi_2009) %>%
  `colnames<-`(c("PatID","0","2","3","4","5","6")) %>%
  pivot_longer(cols=2:7, names_to = "Time_cat", values_to = "FU_eGFR_epi_2009") %>%
  mutate(Cohort="GCKD") %>%
  filter(!is.na(FU_eGFR_epi_2009)) %>%
  mutate(Time_cat=as.numeric(Time_cat))

data.gckd.date <- data.gckd %>%
  dplyr::select(PatID, BL_ein_visdat, FU2_fu_visdat, FU3_fu_visdat, FU4_fu_visdat, FU5_fu_visdat, FU6_fu_visdat) %>%
  `colnames<-`(c("PatID","0","2","3","4","5","6")) %>%
  pivot_longer(cols=2:7, names_to = "Time_cat", values_to = "FU_date") %>%
  filter(!is.na(FU_date)) %>%
  mutate(Time_cat=as.numeric(Time_cat))


data.gckd.long <- full_join(data.gckd.long.2009, data.gckd.long.2021, by=c("PatID", "Time_cat", "Cohort"))
# ggplot(data = data.gckd.long[data.gckd.long$Time_cat==0,], aes(x=FU_eGFR_epi_2021, y=FU_eGFR_epi_2009)) +
#   geom_point() +
#   geom_abline(xintercept = 0, col="red")+
#   theme_bw()



data.gckd <- right_join(data.gckd, left_join(data.gckd.long.2021, data.gckd.date, by=c("PatID","Time_cat")), by="PatID")
data.gckd <- data.gckd %>% 
  mutate(Time_cont = round(difftime(FU_date,BL_ein_visdat,units="days")/365.25,3))
data.gckd$Time_cat <- round(data.gckd$Time_cont)



## -----  PROVALID -----
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
raas <- c("ACE inhibitors", "Renin inhibitors","Angiotensin-II-receptor blockers")

data.tmp1 <- read_excel(file, sheet=1)[,c(clin.cols, gluc.meds, bp.meds, lip.meds)]
data.tmp1$BL_med_dm <- apply(data.tmp1[,gluc.meds],1, function(x) any(!is.na(x)))*1
data.tmp1$BL_med_bp <- apply(data.tmp1[,bp.meds],1, function(x) any(!is.na(x)))*1
data.tmp1$BL_med_lipid <- apply(data.tmp1[,lip.meds],1, function(x) any(!is.na(x)))*1
data.tmp1$BL_raas <- apply(data.tmp1[,raas],1, function(x) any(!is.na(x)))*1

data.tmp1$baseline_date <- data.tmp1$`Date of visit [yyyy/mm/dd]`

data.tmp2 <- read_excel(file, sheet=2)[,lab.cols]
data.tmp3 <- lapply(fu.sheets, function(x) getProvalidFU(x, path=file))
data.tmp3 <- do.call(rbind, data.tmp3)

# Change dates with wrongly have a 3 at start e.g. 3012
data.tmp3$date_yyyy_mm_dd_10[str_detect(data.tmp3$date_yyyy_mm_dd_10, "^3") & !is.na(data.tmp3$date_yyyy_mm_dd_10)] <- na.omit(sub("^.", "2",data.tmp3$date_yyyy_mm_dd_10[str_detect(data.tmp3$date_yyyy_mm_dd_10, "^3")]))
# Change dates with wrongly 9 at second postion e.g. 2912
data.tmp3$date_yyyy_mm_dd_10[substr(data.tmp3$date_yyyy_mm_dd_10,2,2) %in% "9"]   <- "2014-10-29 UTC"
data.tmp3$date_yyyy_mm_dd_10[substr(data.tmp3$date_yyyy_mm_dd_10,3,3) %in% "4"]   <- "2016-04-11 UTC"

data.provalid <- left_join(left_join(data.tmp1, data.tmp2, by="Patient reference [ID]"), 
                           data.tmp3[,c("Patient reference [ID]", "time", "date_yyyy_mm_dd_10","FU_serumcrea")], by="Patient reference [ID]")
rm(data.tmp1,data.tmp2,data.tmp3)

data.provalid <- data.provalid  %>%
  clean_names()
data.provalid <- data.provalid  %>%
  mutate(BL_age = round(as.numeric(as.Date(baseline_date)-as.Date(ISOdate(year_of_birth, 1, 1)))/365.25,0),
         BL_gender = ifelse(gender=="Female","F","M"),
         BL_sex = ifelse(gender=="Female", 1, 0),
         BL_smoking = ifelse(smoking=="never",0,1),
         BL_bmi = calcBMI(height_m=height_cm/100, body_weight_kg),
         BL_diabetic=1,
         BL_serumchol = ifelse(is.na(serum_cholesterol_total_mg_dl), serum_cholesterol_total_mmol_l*38.67 , serum_cholesterol_total_mg_dl),
         BL_hemo = ifelse(is.na(hemoglobin_g_dl), hemoglobin_g_l/10, hemoglobin_g_dl),
         BL_uacr = calcUACR(data.provalid),
         BL_hba1c_perc = ifelse(is.na(hb_a1c_percent), (hb_a1c_mmol_l/10.929)+2.15, hb_a1c_percent ),
         BL_hba1c = ifelse(is.na(hb_a1c_mmol_l), ((hb_a1c_percent-2.15) * 10.929), hb_a1c_mmol_l),
         eth = "non-black",
         FU_age = round(BL_age + as.numeric(data.provalid$time),0),
         FU_eGFR_epi_2009 = ckd_epi(creat = data.provalid$fu_serumcrea, age = FU_age, sex = BL_gender, eth = eth, units = "US"),
         BL_eGFR_epi_2009 = ckd_epi(creat = serum_creatinine_mg_dl, age = BL_age, sex=BL_gender, ethnicity = eth, units="US"),
         FU_eGFR_epi_2021 = ckd_epi_2021(creat = data.provalid$fu_serumcrea, age = FU_age, sex = BL_gender, units = "US"),
         BL_eGFR_epi_2021 = ckd_epi_2021(creat = serum_creatinine_mg_dl, age = BL_age, sex=BL_gender, units="US"),
         BL_crea =serum_creatinine_mg_dl,
         Cohort="PROVALID",
         PatID=as.character(patient_reference_id),
         Country= str_sub(data.provalid$country, -2,-1),
         Time_cont =as.numeric(round(difftime(date_yyyy_mm_dd_10,baseline_date,units="days")/365.25,3)),
         Time_cat=as.numeric(round(Time_cont))) %>%
  dplyr::rename(BL_bpsys = systolic_mm_hg,
                BL_bpdia = diastolic_mm_hg,
                BL_med_lipid = bl_med_lipid,
                BL_med_bp = bl_med_bp,
                BL_med_dm = bl_med_dm,
                BL_raas = bl_raas,
                FU_date=date_yyyy_mm_dd_10) %>%
  dplyr::select(PatID, Time_cat, Time_cont, Country, Cohort, FU_date, baseline_date, FU_eGFR_epi_2009,FU_eGFR_epi_2021, BL_crea,
                BL_age, BL_sex, BL_smoking, BL_bmi, 
         BL_uacr, BL_hemo, BL_serumchol, BL_hba1c_perc, BL_hba1c, BL_bpsys, BL_bpdia,
         BL_med_bp, BL_med_lipid, BL_med_dm, BL_raas) %>%
  filter(!Time_cat < 0)

# Delete duplicated rows with Time_cat == 0 (only one baseline!) - take closest to 0
tmp <- data.provalid[data.provalid$Time_cat==0,]
tmp <- tmp[tmp$PatID %in% tmp[duplicated(tmp[,c("PatID","Time_cat")]),]$PatID,]
tmp1 <- tmp %>%
  group_by(PatID) %>%
  filter(abs(Time_cont) == min(abs(Time_cont)))
tmp1$Time_cont <- 0

data.provalid <- data.frame(rbind(tmp1, anti_join(data.provalid, tmp)))
tripod_flowchart$PROVALID[1] <- length(unique(data.provalid$PatID))
data.provalid <- data.provalid[!data.provalid$Country %in% "NL",]
data.provalid.raas <- data.provalid[,c("PatID", "Cohort","Country","Time_cat","Time_cont","FU_date","BL_raas")]


# =========================== Merge datasets  ==================================
cols <- c("PatID", "Cohort","Country","Time_cat","Time_cont","FU_date", "FU_eGFR_epi_2021", "BL_age", "BL_sex", "BL_smoking",
          "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c", "BL_hba1c_perc", "BL_serumchol", "BL_hemo", "BL_uacr", "BL_med_dm", "BL_med_bp", "BL_med_lipid")
  
data.all <- rbind(data.provalid[,cols], data.gckd[,cols])
data.all$Time_cont[data.all$Time_cont <0] <-0
data.all$Time_cat[data.all$Time_cat <=0] <-0
data.all$Time_cont[data.all$Time_cat ==0] <-0

tripod_flowchart$GCKD[1] <- length(unique(data.all[data.all$Cohort %in% "GCKD",]$PatID))
tripod_flowchart$DIACORE[1] <- length(unique(data.diacore$PatID))

# Lack of data in the NL cohort
tripod_flowchart$PROVALID[2] <- tripod_flowchart$PROVALID[1]-length(unique(data.all[data.all$Cohort %in% "PROVALID",]$PatID))
tripod_flowchart$GCKD[2] <- tripod_flowchart$DIACORE[2] <- 0

data.all <- data.all %>%
  mutate(BL_uacr_log2=log(BL_uacr,2)) %>%
  mutate_at(c("Time_cat","Time_cont","FU_eGFR_epi_2021", "BL_age", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc", "BL_hba1c", "BL_serumchol", "BL_hemo", "BL_uacr"),as.numeric) %>%
  mutate_at(c("Cohort", "Country", "BL_sex", "BL_smoking", "BL_med_dm", "BL_med_bp", "BL_med_lipid"), as.factor) %>%
  mutate_at(c("BL_age", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c", "BL_serumchol", "BL_hemo", "BL_uacr",  "BL_uacr_log2"),
            ~pmin(pmax(.x, quantile(.x, .01, na.rm=T)), quantile(.x, .99, na.rm=T)))




# ======================== Inclusion & Exclusion criteria ======================
data.all <- data.all %>%
  mutate(BL_uacr_log2=log(BL_uacr,2),
         Cohort=ifelse(Cohort %in% "PROVALID", 1, 0),     # Code PROVALID = 1, GCKD = 0
         BL_map = (BL_bpsys+ 2*BL_bpdia)/3,
         BL_hemo=to_numeric(BL_hemo))              # Sorts out values in wrong unit column (g/l instead of g/dL)

# Exclude Patients that fall below 30 at baseline
excl_patid_1 <- data.all[data.all$Time_cat==0 & data.all$FU_eGFR_epi_2021 <31,]$PatID
tripod_flowchart$PROVALID[3] <- length(unique(data.provalid[data.provalid$PatID %in% excl_patid_1,]$PatID))
tripod_flowchart$GCKD[3] <- length(unique(data.gckd[data.gckd$PatID %in% excl_patid_1,]$PatID))

drop.egfr30 <- unique(data.diacore[(data.diacore$Time_cat==0 & data.diacore$FU_eGFR_epi_2021 < 31),]$PatID)
tripod_flowchart$DIACORE[3] <- length(drop.egfr30)

# Exclude patients with no baseline value
df.time <- as.data.frame.matrix(table(data.all$PatID, data.all$Time_cat))
excl_patid_2 <- rownames(df.time[df.time$`0` == 0,])
tripod_flowchart$PROVALID[4] <- length(unique(data.provalid[data.provalid$PatID %in% excl_patid_2,]$PatID))
tripod_flowchart$GCKD[4] <- length(unique(data.gckd[data.gckd$PatID %in% excl_patid_2,]$PatID))

drop.nobaseline <- data.diacore[data.diacore$Time_cat == 0 & is.na(data.diacore$FU_eGFR_epi_2021),]$PatID
tripod_flowchart$DIACORE[4] <- length(drop.nobaseline)

# Only consider patients with more than 2 visits
rec.visits <- apply(as.data.frame.matrix(table(data.all$PatID, data.all$Time_cat)),1, sum)
excl_patid_3 <- unique(names(which(rec.visits <3)))
tripod_flowchart$PROVALID[5] <- length(unique(data.provalid[data.provalid$PatID %in% excl_patid_3,]$PatID))
tripod_flowchart$GCKD[5] <- length(unique(data.gckd[data.gckd$PatID %in% excl_patid_3,]$PatID))

drop.fu3 <- names(which(table(data.diacore$PatID) < 3))
tripod_flowchart$DIACORE[5] <- length(drop.fu3)


# Remove
data.all <- data.all[!data.all$PatID %in% c(excl_patid_1, excl_patid_2,excl_patid_3),]
data.diacore <- data.diacore %>%
  filter(!(PatID %in% drop.egfr30 | PatID %in% drop.fu3))


tripod_flowchart$PROVALID[6] <- length(unique(data.all[data.all$Cohort==1,]$PatID))
tripod_flowchart$GCKD[6] <- length(unique(data.all[data.all$Cohort==0,]$PatID))
tripod_flowchart$DIACORE[6] <- length(unique(data.diacore$PatID))



# =========================== Complete-cases ===================================



# Patient characteristics prior to complete case analysis

df.tmp <- data.frame(data.diacore[data.diacore$Time_cat == 0, ])[c(pred.vars, "BL_uacr", "FU_eGFR_epi_2021", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc")]
df.tmp$Cohort <- "Diacore"
df.tmp$Cohort <- ifelse(df.tmp$Cohort=="0", "GCKD", ifelse(df.tmp$Cohort=="1", "PROVALID", "DIACORE"))

data.tmp <- rbind(data.all[data.all$Time_cat==0,c(pred.vars, "Cohort", "BL_uacr", "FU_eGFR_epi_2021", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc")], 
                  df.tmp)

table1 <- as.data.frame.matrix(print(CreateTableOne(data=data.tmp, vars= c(pred.vars, "BL_uacr", "FU_eGFR_epi_2021", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc"), strata="Cohort", test=F), nonnormal="BL_uacr"))
tbl_table1 <- data.frame("Variable"=rownames(table1), table1) %>%
  `colnames<-`(c("Variable", "GCKD", "PROVALID", "DIACORE"))
rownames(tbl_table1) <- NULL

list.str <- strsplit(gsub("[()]", "", str_trim(tbl_table1[which(str_detect(tbl_table1$Variable, "mean")),]$GCKD, "left")), "\\s+")
tbl_table1[which(str_detect(tbl_table1$Variable, "mean")),]$GCKD <- unlist(lapply(list.str, function(x) paste0(x[1], " \U00B1 ", x[2])))
list.str <- strsplit(gsub("[()]", "", str_trim(tbl_table1[which(str_detect(tbl_table1$Variable, "mean")),]$PROVALID, "left")), "\\s+")
tbl_table1[which(str_detect(tbl_table1$Variable, "mean")),]$PROVALID <- unlist(lapply(list.str, function(x) paste0(x[1], " \U00B1 ", x[2])))
list.str <- strsplit(gsub("[()]", "", str_trim(tbl_table1[which(str_detect(tbl_table1$Variable, "mean")),]$DIACORE, "left")), "\\s+")
tbl_table1[which(str_detect(tbl_table1$Variable, "mean")),]$DIACORE <- unlist(lapply(list.str, function(x) paste0(x[1], " \U00B1 ", x[2])))

write.xlsx(tbl_table1, here::here(out.path, "model_main","tbl_tableone_before_completecaseanalysis.xlsx"), overwrite=T)

tmp.diacore <- data.diacore %>%
  filter(Time_cat==0) %>%
  skim()

tmp.gckd <- data.all %>%
  filter(Time_cat==0 & Cohort == 0) %>%
  skim()

tmp.provalid <- data.all %>%
  filter(Time_cat==0 & Cohort == 1) %>%
  skim()

data.full <- data.all[complete.cases(data.all),]
data.full$Country <- factor(data.full$Country, levels=c("AT", "GE", "HU", "PL", "UK", "NL"))
data.rem <- data.all[!complete.cases(data.all),]
data.diacore <- data.diacore[complete.cases(data.diacore),]

length(unique(data.all$PatID))
length(unique(data.full$PatID))

table(data.all[data.all$Time_cat==0,]$Country)
table(data.full[data.full$Time_cat==0,]$Country)
table(data.rem[data.rem$Time_cat==0,]$Country)

tripod_flowchart$PROVALID[7] <- length(unique(data.full[data.full$Cohort==1,]$PatID))
tripod_flowchart$GCKD[7] <- length(unique(data.full[data.full$Cohort==0,]$PatID))
tripod_flowchart$DIACORE[7] <- length(unique(data.diacore$PatID))

tripod_flowchart$PROVALID[8] <- length(unique(data.full$PatID))
tripod_flowchart$GCKD[8] <- length(unique(data.full$PatID))
tripod_flowchart$DIACORE[8] <- length(unique(data.diacore$PatID))

tripod_flowchart


# ==================== True progression slope estimation ============================

# Development cohort
df.tmp <- data.full %>%
  group_by(PatID) %>%
  summarise(eGFRslope=summary(lm(FU_eGFR_epi_2021 ~ Time_cont))$coefficients[2]) %>%
  data.frame() %>%
  `colnames<-`(c("PatID", "true.slope"))
data.full <- left_join(data.full, df.tmp, by="PatID")
data.full$true.prob <- (data.full$true.slope <= slope_cutpoint)*1

# Validation cohort
df.tmp <- data.diacore %>%
  group_by(PatID) %>%
  summarise(eGFRslope=summary(lm(FU_eGFR_epi_2021 ~ Time_cont))$coefficients[2]) %>%
  data.frame() %>%
  `colnames<-`(c("PatID", "true.slope"))
data.diacore <- left_join(data.diacore, df.tmp, by="PatID")
data.diacore$true.prob <- (data.diacore$true.slope <= slope_cutpoint)*1

data.model <- list("development"=data.full, "validation"=data.diacore)
saveRDS(data.model, here::here(out.path, "data", "cohort_data_model.rds"))

# ==== Save PROVALID with RAAS ====

data.provalid.withRAAS <- left_join(data.full[data.full$Cohort==1,], data.provalid.raas, by=c("PatID", "Country","Time_cat","Time_cont")) 
data.provalid.withRAAS <- data.provalid.withRAAS[complete.cases(data.provalid.withRAAS$BL_raas),]

saveRDS(data.provalid.withRAAS, here::here(out.path, "data", "cohort_data_provalid_raas.rds"))
