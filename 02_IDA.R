###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Initial data analysis: data screening & cleaning
###################################################



# ----- Descriptive statistics

# Number of subjects per country and cohort
table(data.full[!duplicated(data.full$PatID),]$Cohort)
table(data.full[!duplicated(data.full$PatID),]$Country)


# Distribution of baseline variables between PROVALID and GCKD
data.all %>% 
  filter(Time==0 & Cohort == 1) %>%
  skim()

data.all %>% 
  filter(Time==0 & Cohort == 0) %>%
  skim()

data.all %>% 
  filter(Time==0) %>%
  skim()


# Distribution of follow-up time
data.full %>% group_by(PatID) %>% summarise(maxT=max(Time, na.rm=T)) %>% pull(maxT) %>% summary()

# Histogram of all baseline variables
par(mfrow = c(3,5))
hist(data.full[data.full$Time==0,]$BL_age, col="grey", xlab=NULL, main = "Age [years]")
barplot(table(data.full[data.full$Time==0,]$BL_sex), ylab = "Frequency", main = "Gender [1=female]")
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
x_tmp <- model.matrix(FU_eGFR_epi~BL_age + BL_sex + BL_smoking + BL_bmi + BL_bpsys + BL_bpdia + BL_hba1c_perc 
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

cont.vars <- c("BL_age", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc", 
               "BL_serumchol", "BL_hemo", "BL_uacr_log2")
pairw.int <- expand.grid(cont.vars, cont.vars)
pairw.int <- pairw.int[!pairw.int[,1]==pairw.int[,2],]
plots <- lapply(1:nrow(pairw.int), function(x) plot_fun(x=pairw.int[x,1], y=pairw.int[x,2], time=0))
#plots


# ---- Check non-linearities
plots <- lapply(cont.vars, function(x) plot_fun(x=x, y="FU_eGFR_epi", time=1))
#plots

# Non-linearity for BL_age
plots[[3]]



# -------------- TABLE 1
data.tmp <- data.full[data.full$Time==0,]
table1 <- CreateTableOne(data=data.tmp, vars= c("BL_age", "BL_sex", "BL_smoking", "BL_bmi", "BL_map","BL_bpsys", "BL_bpdia", "BL_hba1c_perc", 
                                          "BL_serumchol", "BL_hemo", "BL_uacr", "FU_eGFR_epi","BL_med_dm", "BL_med_bp", "BL_med_lipid"), strata="Cohort", test=F)
write.xlsx(as.data.frame.matrix(print(table1)), paste0(out.path, "tbl_tableone_dev.xlsx"), overwrite=T)
CreateTableOne(data=data.tmp, vars= c("BL_age", "BL_sex", "BL_smoking", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c_perc", 
                                      "BL_serumchol", "BL_hemo", "BL_uacr", "BL_med_dm", "BL_med_bp", "BL_med_lipid"),  test=F)

data.tmp <- data.diacore[data.diacore$Time == 0,]
table1 <- CreateTableOne(data=data.tmp, vars= c("BL_age", "BL_sex", "BL_smoking", "BL_bmi", "BL_map","BL_bpsys", "BL_bpdia", "BL_hba1c_perc", 
                                                "BL_serumchol", "BL_hemo", "BL_uacr", "FU_eGFR_epi","BL_med_dm", "BL_med_bp", "BL_med_lipid"), test=F)
write.xlsx(as.data.frame.matrix(print(table1)), paste0(out.path, "tbl_tableone_val.xlsx"), overwrite=T)


# ------------- Table: longitudinal eGFR measurements stratified by cohort
table(data.full$Time, data.full$Cohort)
table(data.diacore$Time)


########################################################################################################
# --------------- Sample size calculation -------------------------------------------------------
########################################################################################################
library(pmsampsize)

pmsampsize(type = "c", rsquared = 0.7, parameters = 13, intercept = 78.4, sd = 21.4, shrinkage = 0.99)

#########################################################################################################
# ----------- INITIAL MODEL BUILDING (Apparent)- Non-linearities? Interactions? -------------------------
##########################################################################################################

# # ---------- Non-linear effects
# # Fit model and check for non-linear behavior in residuals vs independent variables
# fit.init <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_sex + BL_smoking + BL_bpsys + BL_bpdia + 
#                                                       BL_hba1c + BL_serumchol + BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid)), 
#                 random=list(~1|Country, ~1+Time|PatID), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
# ggplot(data.frame(x1=data.full$BL_uacr_log2,residual=resid(fit.init, scaled=TRUE), time=data.full$Time),
#        aes(x=x1,y=residual, color=time)) +
#   geom_smooth(se=T) +
#   theme_bw() +
#   facet_wrap(~time)
# ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.init, scaled=TRUE), time=data.full$Time),
#        aes(x=x1,y=residual, color=time)) +
#   geom_smooth(se=T) +
#   theme_bw() +
#   facet_wrap(~time)
# 
# # Non-linear term for BL_uacr_log2 according to plots above
# fit.tmp1 <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_sex + BL_smoking + BL_bpsys + BL_bpdia + 
#                                                       BL_hba1c + BL_serumchol + BL_hemo + rcs(BL_uacr_log2,df=2) + BL_med_dm + BL_med_bp + BL_med_lipid)), 
#                 random=list(~1+Time|PatID, ~1|Country), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
# anova(fit.init, fit.tmp1)
# ggplot(data.frame(x1=data.full$BL_uacr_log2,residual=resid(fit.tmp1, scaled=TRUE), time=data.full$Time),
#        aes(x=x1,y=residual, color=time)) +
#   geom_smooth(se=T) +
#   theme_bw() +
#   facet_wrap(~time)
# ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.tmp1, scaled=TRUE), time=data.full$Time),
#        aes(x=x1,y=residual, color=time)) +
#   geom_smooth(se=T) +
#   theme_bw() +
#   facet_wrap(~time)
# 
# # Non-linear term for BL_serumchol according to plots above
# fit.tmp2 <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_sex + BL_smoking + BL_bpsys + BL_bpdia + 
#                                                       BL_hba1c + rcs(BL_serumchol,df=2) + BL_hemo + rcs(BL_uacr_log2,df=2) + BL_med_dm + BL_med_bp + BL_med_lipid)), 
#                 random=list(~1+Time|PatID, ~1|Country), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
# anova(fit.tmp1, fit.tmp2)
# ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.tmp2, scaled=TRUE), time=data.full$Time),
#        aes(x=x1,y=residual, color=time)) +
#   geom_smooth(se=T) +
#   theme_bw() +
#   facet_wrap(~time)
# ggplot(data.frame(x1=data.full$BL_BL_uacr_log2,residual=resid(fit.tmp2, scaled=TRUE), time=data.full$Time),
#        aes(x=x1,y=residual, color=time)) +
#   geom_smooth(se=T) +
#   theme_bw() +
#   facet_wrap(~time)
# 
# 
# 
# # --------------------------- Pairwise interactions -------------------------------------------
# out.step <- stepAIC(fit.tmp2, 
#                     direction = 'forward', 
#                     scope=list(upper = ~ (Time + Time * (BL_bmi + BL_age + BL_sex + BL_smoking + 
#                                                            (BL_bpsys + BL_bpdia)^2  + rcs(BL_serumchol, df = 2) + 
#                                                            rcs(BL_uacr_log2, df = 2) + (BL_hba1c + BL_med_dm)^2 + (BL_hemo +BL_med_bp)^2 + 
#                                                            BL_med_lipid))))
# anova(out.step)
# 
# 
# 
# # -------------------  Mixed model - Substantial increase in R2 of the more complex model ? ---------
# # Linear model
# fit.final <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_sex + BL_smoking + BL_bpsys + BL_bpdia + BL_hba1c + BL_serumchol +
#                                                        BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid)),
#                  random=list(~1|Country,~1+Time|PatID), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
# 
# data.full.t0 <- data.full[data.full$Time==0,]
# res <- LongPred_ByBase(lmeObject=fit.final, newdata = data.full.t0, timeVar = "Time", idVar="PatID", idVar2="Country",times = unique(data.full$Time), all_times=F)
# data.preds <- full_join(data.full, res$Pred[,c("PatID", "Time","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")], by=c("PatID", "Time"))
# 
# plot_calibration(yobs=data.preds[data.preds$Time==1,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==1,]$pred, time=1)
# plot_calibration(yobs=data.preds[data.preds$Time==2,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==2,]$pred, time=2)
# plot_calibration(yobs=data.preds[data.preds$Time==3,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==3,]$pred, time=3)
# plot_calibration(yobs=data.preds[data.preds$Time==4,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==4,]$pred, time=4)
# plot_calibration(yobs=data.preds[data.preds$Time==5,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==5,]$pred, time=5)
# 
# 
# # Complex model
# fit.final <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_sex + BL_smoking + 
#                                                        (BL_bpsys + BL_bpdia)^2  + rcs(BL_serumchol, df = 2) + 
#                                                        rcs(BL_uacr_log2, df = 2) + (BL_hba1c + BL_med_dm)^2 + (BL_hemo +BL_med_bp)^2 + 
#                                                        BL_med_lipid)),
#                  random=list(~1|Country,~1+Time|PatID), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
# 
# data.full.t0 <- data.full[data.full$Time==0,]
# res <- LongPred_ByBase(lmeObject=fit.final, newdata = data.full.t0, timeVar = "Time", idVar="PatID", idVar2="Country",times = unique(data.full$Time), all_times=F)
# data.preds <- full_join(data.full, res$Pred[,c("PatID", "Time","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")], by=c("PatID", "Time"))
# 
# plot_calibration(yobs=data.preds[data.preds$Time==1,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==1,]$pred, time=1)
# plot_calibration(yobs=data.preds[data.preds$Time==2,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==2,]$pred, time=2)
# plot_calibration(yobs=data.preds[data.preds$Time==3,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==3,]$pred, time=3)
# plot_calibration(yobs=data.preds[data.preds$Time==4,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==4,]$pred, time=4)
# plot_calibration(yobs=data.preds[data.preds$Time==5,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==5,]$pred, time=5)


# 16.08: the more complex model shows no significant increase (~ 10%) in R2 compared to the linear model, 
#        which is why the linear model is continued.