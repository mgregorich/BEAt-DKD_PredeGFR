#===============================================================================#
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Initial data analysis: data screening & cleaning
#===============================================================================#



# =================== Descriptive statistics & data inspection =================

# Number of subjects per country and cohort
table(data.full[!duplicated(data.full$PatID),]$Cohort)
table(data.full[!duplicated(data.full$PatID),]$Country)

# Number of subjects stratified by country and follow-up
table(data.full$Country,data.full$Time_cat)

# Number of time points per cohort
table(data.full[data.full$Cohort==0,]$Time_cat)
table(data.full[data.full$Cohort==1,]$Time_cat)
table(data.diacore$Time_cat)


# Distribution of baseline variables between PROVALID and GCKD
data.all %>% 
  filter(Time_cont==0 & Cohort == 1) %>%
  skim()

data.all %>% 
  filter(Time_cont==0 & Cohort == 0) %>%
  skim()

data.all %>% 
  filter(Time_cont==0) %>%
  skim()


# Distribution of follow-up time
data.full %>% group_by(PatID) %>% summarise(maxT=max(Time_cont, na.rm=T)) %>% pull(maxT) %>% summary()

# Histogram of all baseline variables
par(mfrow=c(1,1))
hist(data.full[data.full$Time_cont==0,]$BL_age, col="grey", xlab=NULL, main = "Age [years]")
barplot(table(data.full[data.full$Time_cont==0,]$BL_sex), ylab = "Frequency", main = "Gender [1=female]")
hist(data.full[data.full$Time_cont==0,]$BL_bmi, col="grey", xlab=NULL, main = "BMI [kg/m]")
barplot(table(data.full[data.full$Time_cont==0,]$BL_smoking), ylab = "Frequency", main = "Smoking [1=ever]")
hist(data.full[data.full$Time_cont==0,]$BL_bpsys, col="grey", xlab=NULL, main = "Systolic blood pressure [mmHg]")
hist(data.full[data.full$Time_cont==0,]$BL_bpdia, col="grey", xlab=NULL, main = "Diastolic blood pressure [mmHg]")
hist(data.full[data.full$Time_cont==0,]$BL_serumchol, col="grey", xlab=NULL, main = "Serum Cholesterin")
hist(data.full[data.full$Time_cont==0,]$BL_hemo, col="grey", xlab=NULL, main = "Hemoglobin [g/dL]")
hist(data.full[data.full$Time_cont==0,]$BL_hba1c, col="grey", xlab=NULL, main = "HbA1c [mmol/mol]")
hist(data.full[data.full$Time_cont==0,]$BL_hba1c_perc, col="grey", xlab=NULL, main = "HbA1c [%]")
hist(data.full[data.full$Time_cont==0,]$BL_uacr, col="grey", xlab=NULL, main = "UACR")
hist(data.full[data.full$Time_cont==0,]$BL_uacr_log2, col="grey", xlab=NULL, main = "Log2 UACR")
barplot(table(data.full[data.full$Time_cont==0,]$BL_med_bp), ylab = "Frequency", main = "Blood pressure-lowering meds [1=yes]")
barplot(table(data.full[data.full$Time_cont==0,]$BL_med_dm), ylab = "Frequency", main = "Glucose-lowering meds [1=yes]")
barplot(table(data.full[data.full$Time_cont==0,]$BL_med_lipid), ylab = "Frequency", main = "Lipid-lowering meds [1=yes]")


# Histograms per time point for eGFR
par(mfrow=c(1,1))
data.full %>% filter(Time_cat==0) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=0", xlab=NULL) 
data.full %>% filter(Time_cat==1) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=1", xlab=NULL) 
data.full %>% filter(Time_cat==2) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=2", xlab=NULL) 
data.full %>% filter(Time_cat==3) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=3", xlab=NULL) 
data.full %>% filter(Time_cat==4) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=4", xlab=NULL) 
data.full %>% filter(Time_cat==5) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=5", xlab=NULL) 
data.full %>% filter(Time_cat==6) %>%  pull(FU_eGFR_epi) %>% hist(main = "eGFR - time=6", xlab=NULL) 


# ======================= Correlation heatmap ===================================
x_tmp <- data.frame(model.matrix(FU_eGFR_epi~BL_age + BL_sex + BL_smoking + BL_bmi + BL_map + BL_hba1c 
                      + BL_serumchol + BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid, 
                      data=data.full[data.full$Time_cat==0,])[,-1])
colnames(x_tmp) <- c("Age", "Sex", "Smoking", "BMI", "MAP", "Hba1C", "Serum chol.", "Hemoglobin", "log2UACR", "Glucose-low. Med.", "BP-low. Med.", "Lipid-low. Med.")
corr <- cor(x_tmp, method = "spearman")
(abs(corr)>0.5)

melted_cormat <- melt(round(corr,4))
corr[lower.tri(corr)]<-NA
diag(corr) <- NA
melted_cormat$value_label <- round_0(c(corr),2)
melted_cormat$value_label <- ifelse(melted_cormat$value_label %in% "NA", NA, melted_cormat$value_label)

ggplot(data = melted_cormat, aes(x=Var1, y=reorder(Var2, desc(Var2)), fill=value)) + 
  geom_tile() +
  geom_text(aes(label = value_label), size = 4) +
  scale_x_discrete("") +
  scale_y_discrete("") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_bw() +
  theme(legend.position = "right", text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust=1), 
        legend.text = element_text(size=12), legend.title = element_text(size=14))
ggsave(paste0(out.path, "fig_correlation.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')



# ================================ Bivariate distribution ======================

cont.vars <- c("BL_age", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc", 
               "BL_serumchol", "BL_hemo", "BL_uacr_log2")
pairw.int <- expand.grid(cont.vars, cont.vars)
pairw.int <- pairw.int[!pairw.int[,1]==pairw.int[,2],]
plots <- lapply(1:nrow(pairw.int), function(x) plot_fun(x=pairw.int[x,1], y=pairw.int[x,2], time=0))



# ========================== Check non-linearities =============================
plots <- lapply(cont.vars, function(x) plot_fun(x=x, y="FU_eGFR_epi", time=1))
#plots

# Non-linearity for BL_age
plots[[3]]


# ======================== Longitudinal measurement analysis =========
summary(data.full[data.full$Time_cat==0,]$Time_cont)
summary(data.full[data.full$Time_cat==1,]$Time_cont)
summary(data.full[data.full$Time_cat==2,]$Time_cont)
summary(data.full[data.full$Time_cat==3,]$Time_cont)
summary(data.full[data.full$Time_cat==4,]$Time_cont)
summary(data.full[data.full$Time_cat==5,]$Time_cont)
summary(data.full[data.full$Time_cat==6,]$Time_cont)
summary(data.full[data.full$Time_cat==7,]$Time_cont)


# --- Median and IQR of eGFR decline per year & true progression estimation

summary(data.full[data.full$Cohort==0,]$true.slope)
summary(data.full[data.full$Cohort==1,]$true.slope)
summary(data.diacore$true.slope)


# --- Average years of follow-up ----

tmp <- data.full %>%
  group_by(PatID) %>%
  summarise(max_time=max(Time_cont))
mean(tmp$max_time)
hist(tmp$max_time)

tmp <- data.diacore %>%
  group_by(PatID) %>%
  summarise(max_time=max(Time_cont))
mean(tmp$max_time)
hist(tmp$max_time)

# ------------- Table: longitudinal eGFR measurements stratified by cohort
table(data.full$Time_cat, data.full$Cohort)
table(data.diacore$Time_cat)


# =========================== CKD stage at baseline ============================
data.tmp <- data.full %>%
  filter(Time_cat==0) %>%
  mutate(CKDstage = as.factor(ifelse(FU_eGFR_epi >=90, 1, ifelse(FU_eGFR_epi >=60, 2, 3))))

ggplot(data.tmp, aes(x=CKDstage, fill=Cohort)) +
  geom_bar(stat="count") +
  theme_bw() +
  facet_wrap(~Cohort)

Nobs.CKDstages <- list("CKD1"=(sum(data.tmp$CKDstage==1)/nrow(data.tmp))*100,
     "CKD2"=(sum(data.tmp$CKDstage==2)/nrow(data.tmp))*100,
     "CKD3"=(sum(data.tmp$CKDstage==3)/nrow(data.tmp))*100)

data.tmp <- data.diacore %>%
  filter(Time_cat==0) %>%
  mutate(CKDstage = as.factor(ifelse(FU_eGFR_epi >=90, 1, ifelse(FU_eGFR_epi >=60, 2, 3))))
Nobs.CKDstages <- list("CKD1"=(sum(data.tmp$CKDstage==1)/nrow(data.tmp))*100,
                       "CKD2"=(sum(data.tmp$CKDstage==2)/nrow(data.tmp))*100,
                       "CKD3"=(sum(data.tmp$CKDstage==3)/nrow(data.tmp))*100)




# ================================ TABLE 1: Patient characteristics at baseline ====================================

# Strafied by cohort
df.tmp <- data.frame(data.diacore[data.diacore$Time_cat == 0, ])[c(pred.vars, "FU_eGFR_epi", "BL_bpsys", "BL_bpdia")]
df.tmp$Cohort <- "Diacore"
df.tmp$Cohort <- ifelse(df.tmp$Cohort=="0", "GCKD", ifelse(df.tmp$Cohort=="1", "PROVALID", "DIACORE"))

data.tmp <- rbind(data.full[data.full$Time_cat==0,c(pred.vars, "Cohort","FU_eGFR_epi", "BL_bpsys", "BL_bpdia")], 
                  df.tmp)

table1 <- as.data.frame.matrix(print(CreateTableOne(data=data.tmp, vars= c(pred.vars, "FU_eGFR_epi", "BL_bpsys", "BL_bpdia"), strata="Cohort", test=F)))

tbl_table1 <- data.frame("Variable"=rownames(table1), table1) %>%
  `colnames<-`(c("Variable", "GCKD", "PROVALID", "DIACORE"))
rownames(tbl_table1) <- NULL

na_count <-sapply(data.tmp[,c(pred.vars, "FU_eGFR_epi", "BL_bpsys", "BL_bpdia")], 
                  function(y) sum(length(which(is.na(y)))))
tbl_table1$NA.count <- c("",na_count)

write.xlsx(tbl_table1, paste0(out.path, "tbl_tableone_all.xlsx"), overwrite=T)

# For the development cohort in total
data.tmp <- data.full[data.full$Time_cat==0,]
CreateTableOne(data=data.tmp, vars= c(pred.vars, "FU_eGFR_epi"), test=F)


# ========================= Median rate of eGFR decline ===========================

summary(data.full[data.full$Time_cat==0 & data.full$Cohort==1,]$true.slope)  #PROVALID
summary(data.full[data.full$Time_cat==0 & data.full$Cohort==0,]$true.slope) #GCKD
summary(data.diacore[data.diacore$Time_cat==0,]$true.slope)


# =============================== Sample size calculation ======================
pmsampsize(type = "c", rsquared = 0.7, parameters = 13, intercept = 78.4, sd = 21.4, shrinkage = 0.99)
