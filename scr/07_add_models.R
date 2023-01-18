#===============================================================================#
# Author: Mariella Gregorich
# Date: 15.12.2022
# Info: Comparison of models as suggested by reviewers' comments
#===============================================================================#



# ============ Data preparation, IDA and general setup ==========================

rm(list = ls())
source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)


# ============== (1) Model 1: Clinical covariates only ==============

data.full$Time <- data.full$Time_cont
lmer_formula_model_clin <- as.formula("FU_eGFR_epi_2021 ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID))")
fit.clin <- lmer(lmer_formula_model_clin,
                 data=data.full, REML=F, control=lmerControl(optimizer="bobyqa"))


# UPDATE prediction with BL value
data.full.t0 <- data.full[data.full$Time==0,]
data.full.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=fit.clin, newdata = data.full.t0, timeVar = "Time", idVar="PatID", idVar2="Country",
                         times = unique(data.full$Time_cat)[-1], all_times=F)
data.full$Time <- round(data.full$Time,0)
data.sens <- full_join(data.full, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))

# ================== Internal validation  =====================================

# ----  Save final model and extract components for shiny ----
saveRDS(fit.clin, here::here(out.path, "model_clin","predmodel_lmerObject_model_clin.rds"))

# ----  Residual check -----
check_residuals(fit.clin, path = here::here(out.path, "model_clin","fig_residual_analysis_model_clin.tiff"))


# ====================== Internal-external validation =========================

set.seed(123)
j=1; b=nboot

# Cross-validated predictions
data.full$fold <- as.numeric(data.full$Country)
res <- intext_crossvalidate(data.full, NA, mod_formula=lmer_formula_model_clin, return_preds = T, sa=T)
df.preds <- res$pred

# Cross-validated performance measures
plan(multisession, gc=T, workers=detectCores()*.6)
res_cv_boot <- future_lapply(1:b, function(x)  intext_crossvalidate(draw_bootstrap_sample(data.full),x, mod_formula=lmer_formula_model_clin, return_preds = F, sa=T), future.seed = T)
plan(sequential)
df.stats <- data.frame(do.call(rbind, lapply(res_cv_boot, `[[`, 1)))


# =============== Table: Model discrimination, precision and fit =====================
tbl_performance <- list()
tbl_performance[[1]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], mean, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[2]]<- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], quantile, 0.05, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[3]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], quantile, 0.95, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))
df.tmp <- data.frame(do.call(cbind, tbl_performance))
tbl_performance[[4]] <-df.tmp %>%
  data.frame() %>%
  mutate(Time=Time,
         N=round_0(Nobs,3),
         N.ci=paste0("(",round_0(Nobs.1,3),", ",round_0(Nobs.2,3),")"),
         R2=round_0(R2,3),
         R2.ci = paste0("(",round_0(R2.1,3),", ",round_0(R2.2,3),")"),
         C=round_0(C,3),
         C.ci = paste0("(",round_0(C.1,3),", ",round_0(C.2,3),")"),
         CS=round_0(CalbSlope,3),
         CS.ci = paste0("(",round_0(CalbSlope.1,3),", ",round_0(CalbSlope.2,3),")")) %>%
  dplyr::select(Time, N, N.ci, R2, R2.ci, C, C.ci, CS, CS.ci) 

names(tbl_performance) <- c("avg", "ci.lo", "ci.up", "reported")
write.xlsx(tbl_performance, here::here(out.path, "model_clin", "tbl_perform_val_full_sa.xlsx"), overwrite = T)

# =================== External validation ======================================

rm(list=ls())
source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)
risk_model <- readRDS(here::here(out.path, "model_clin","predmodel_lmerObject_model_clin.rds"))

data.diacore$Time <- data.diacore$Time_cont
data.diacore.t0 <- data.frame(data.diacore[data.diacore$Time_cat==0,])
data.diacore.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=risk_model, 
                         newdata = data.diacore.t0, 
                         cutpoint = slope_cutpoint,
                         timeVar = "Time", idVar="PatID", idVar2="Country",
                         times =sort(unique(data.diacore$Time_cat)), 
                         all_times=F)

# Summarize and prepare output
data.diacore$Time <- round(data.diacore$Time,0)
data.diacore.new <- full_join(data.diacore, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))


# -------------------- Calibration plot ---------------------------------------

for(t in 2:7){
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi_2021, yhat=data.diacore.new[data.diacore.new$Time==t,]$prior.pred,
                        cohort="extval", time=t, save=T, folder="model_clin", out.path = out.path, type="preUp_model_clin")
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi_2021, yhat=data.diacore.new[data.diacore.new$Time==t,]$pred,
                        cohort="extval", time=t, save=T, folder="model_clin", out.path = out.path, type="postUp_model_clin")
}

ggplot(data.diacore.new[data.diacore.new$Time_cat <=5,], aes(x=pred, y=FU_eGFR_epi_2021, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Predicted eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = expression(paste("Observed eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=22)) 
ggsave(here::here(out.path, "model_clin", "fig_calibration_all_diacore_post_model_clin.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')


df <- melt(data.diacore.new[data.diacore.new$Time_cat <=5 ,c("Time_cont","FU_eGFR_epi_2021", "pred","prior.pred")], id.vars = c("Time_cont", "FU_eGFR_epi_2021"))
levels(df$variable) <- c("Post-update", "Pre-update")
df$variable <- relevel(df$variable, "Pre-update")
df <- df[df$Time_cont>0,]

ggplot(df, aes(x=value, y=FU_eGFR_epi_2021, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Predicted eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = expression(paste("Observed eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=16), strip.text = element_text(size = 18)) +
  facet_wrap(~variable)
ggsave(here::here(out.path, "model_clin", "fig_calibration_all_diacore_prepost_model_clin.tiff"),  width=12, height=6, device='tiff', dpi=350, compression = 'lzw')



# ------------------ Table: Performance measures -------------------------------

b=nboot
plan(multisession, gc=T, workers=detectCores()*.6)
res_ext_boot <- future_lapply(1:b, function(x){
  data.boot <- data.frame(draw_bootstrap_sample(data.diacore, cv=F))
  data.boot.t0 <- data.boot[data.boot$Time_cat==0,]
  res <- update_PredByBase(lmerObject=risk_model,
                           newdata=data.boot.t0,
                           cutpoint = slope_cutpoint,
                           timeVar = "Time", idVar="PatID", idVar2="Country",
                           times =seq(0,8,1), 
                           all_times=F)
  data.boot$Time <- data.boot$Time_cat
  data.boot.new <- full_join(data.boot, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], 
                             by=c("PatID", "Time"))
  
  data.boot.list <- split(data.boot.new, as.factor(data.boot.new$Time)) 
  res <- lapply(data.boot.list, function(x) eval_preds(pred=x$pred, obs=x$FU_eGFR_epi_2021, lmerObject=risk_model))
  df.res <- data.frame("Boot"=x,"Time"=as.numeric(names(res)), do.call(rbind, res))
  
  return(df.res)
}
, future.seed = T)
plan(sequential)

df.stats <- data.frame(do.call(rbind, res_ext_boot))



# Summarize and prepare output
tbl_performance <- list()
tbl_performance[[1]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], mean, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[2]]<- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], quantile, 0.05, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[3]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], quantile, 0.95, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))
df.tmp <- data.frame(do.call(cbind, tbl_performance))
tbl_performance[[4]] <-df.tmp %>%
  data.frame() %>%
  mutate(Time=Time,
         N=round_0(Nobs,3),
         N.ci=paste0("(",round_0(Nobs.1,3),", ",round_0(Nobs.2,3),")"),
         R2=round_0(R2,3),
         R2.ci = paste0("(",round_0(R2.1,3),", ",round_0(R2.2,3),")"),
         C=round_0(C,3),
         C.ci = paste0("(",round_0(C.1,3),", ",round_0(C.2,3),")"),
         CS=round_0(CalbSlope,3),
         CS.ci = paste0("(",round_0(CalbSlope.1,3),", ",round_0(CalbSlope.2,3),")")) %>%
  dplyr::select(Time, N, N.ci, R2, R2.ci, C, C.ci, CS, CS.ci) 

names(tbl_performance) <- c("avg", "ci.lo", "ci.up", "reported")
write.xlsx(tbl_performance, here::here(out.path, "model_clin", "tbl_perform_extval_full_model_clin.xlsx"), overwrite = T)

# ================= Table 2: Model coefficients ================================
df <- data.frame(variable=rownames(summary(risk_model)$coefficients),
                 effect=summary(risk_model)$coefficients[,1],
                 lower=summary(risk_model)$coefficients[,1] - 1.96*summary(risk_model)$coefficients[,2],
                 upper=summary(risk_model)$coefficients[,1] + 1.96* summary(risk_model)$coefficients[,2]) 
df$variable <- factor(df$variable, levels = df$variable[length(df$variable):1])
df$Group <- ifelse(str_detect(df$variable, "Time:"), "Main effect", "Interaction effect")
df <- mutate(df, Group = fct_rev(Group)) 
df$variable <- as.factor(str_replace(df$variable, "Time:", ""))
df$variable <- factor(df$variable, levels=c("(Intercept)","Time","BL_bmi","BL_smoking1","BL_map","BL_hba1c","BL_serumchol", "BL_hemo",
                                            "BL_uacr_log2", "BL_med_dm1","BL_med_bp1","BL_med_lipid1"))
levels(df$variable) <- list("Intercept"="(Intercept)","Time"="Time", BMI="BL_bmi", "Smoking (ever)"="BL_smoking1", MAP="BL_map", HbA1c="BL_hba1c",
                            "Serum Chol."="BL_serumchol","Hemoglobin"="BL_hemo", "log2 UACR"="BL_uacr_log2", "Glucose-low. Med."="BL_med_dm1", "Blood pressure-low. Med."="BL_med_bp1", "Lipid-low. Med."="BL_med_lipid1")


df[,2:4] <- apply(df[,2:4],2, round, digits=3)
df$name <- rownames(df)
df$Variable <- c(rep("Constant",2), rep(c("Age","Sex","BMI", "Smoking", "MAP", "GL Med.", "BPL Med.","LL Med."),2))
df$CI <- paste0("(",df$lower,", " ,df$upper,")")
tbl_fixeff <- cbind(df[!str_detect(df$name, "Time"),c(7,2,8)], df[str_detect(df$name, "Time"),c(2,8)])
colnames(tbl_fixeff) <- c("Variable", "Effect.Baseline", "95% CI", "Effect.Slope", "95% CI")
write.xlsx(tbl_fixeff, here::here(out.path, "model_clin", "tbl_coeff_unstd_model_clin.xlsx"), overwrite = TRUE)


# ===================== Forest plot of standardized coefficients ===============
data.full$Time <- data.full$Time_cont
data.scaled <- data.full %>%mutate_at(vars(Time, BL_bmi, BL_map, BL_hba1c, 
                                           BL_serumchol, BL_hemo, BL_uacr_log2, FU_eGFR_epi_2021), scale)
fit.scaled <- lmer(FU_eGFR_epi_2021 ~ (Time + Time  *(BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +                                                  
                                                        BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                   data=data.scaled, REML=F, control=lmerControl(optimizer="bobyqa"))

df <- data.frame(variable=rownames(summary(fit.scaled)$coefficients),
                 effect=summary(fit.scaled)$coefficients[,1],
                 lower=summary(fit.scaled)$coefficients[,1] - 1.96*summary(fit.scaled)$coefficients[,2],
                 upper=summary(fit.scaled)$coefficients[,1] + 1.96* summary(fit.scaled)$coefficients[,2]) 
df$variable <- factor(df$variable, levels = df$variable[length(df$variable):1])
df$Group <- ifelse(str_detect(df$variable, "Time:"), "Interaction effect", "Main effect")
df <- mutate(df, Group = fct_rev(Group)) 
df$variable <- as.factor(str_replace(df$variable, "Time:", ""))
df$variable <- factor(df$variable, levels=c("(Intercept)","Time","BL_bmi","BL_smoking1","BL_map","BL_hba1c","BL_serumchol", "BL_hemo",
                                            "BL_uacr_log2", "BL_med_dm1","BL_med_bp1","BL_med_lipid1"))
levels(df$variable) <- list("Intercept"="(Intercept)","Time"="Time", BMI="BL_bmi", "Smoking (ever)"="BL_smoking1", MAP="BL_map", HbA1c="BL_hba1c",
                            "Serum Chol."="BL_serumchol","Hemoglobin"="BL_hemo", "log2 UACR"="BL_uacr_log2", "Glucose-low. Med."="BL_med_dm1", "Blood pressure-low. Med."="BL_med_bp1", "Lipid-low. Med."="BL_med_lipid1")

df[,2:4] <- apply(df[,2:4],2, round, digits=3)
df$name <- rownames(df)
df$Variable <- c(rep("Constant",2), rep(c("BMI", "Smoking", "MAP", "HbA1c","Serum chol.",
                                          "Hemoglobin", "log2UACR", "GL Med.", "BPL Med.","LL Med."),2))
df$CI <- paste0("(",df$lower,", " ,df$upper,")")
tbl_fixeff <- cbind(df[!str_detect(df$name, "Time"),c(7,2,8)], df[str_detect(df$name, "Time"),c(2,8)])
colnames(tbl_fixeff) <- c("Variable", "Effect.Baseline", "95% CI", "Effect.Slope", "95% CI")
write.xlsx(tbl_fixeff, here::here(out.path, "model_clin", "tbl_coeff_std_model_clin.xlsx"), overwrite = TRUE)

ggplot(data=df[-1,], aes(y=reorder(variable,desc(variable)), x=effect, xmin=lower, xmax=upper)) +
  geom_point(size=2, shape=1) + 
  geom_errorbarh(height=.25) +
  scale_y_discrete("", labels=c("HbA1c" = expression(HbA[1][c]), "log2 UACR"=expression(paste(log[2], ' UACR')))) +
  scale_x_continuous("Standardized Effect", limits = c(-0.4,0.4), breaks=seq(-0.4,0.4,0.1)) +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  facet_wrap(~Group) +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave(here::here(out.path, "model_clin", "plot_forest_standardized_model_clin.tiff"),  width=8, height=4, device='tiff', dpi=350, compression = 'lzw')







# ============== (2) Model 2: UACR only ==============
rm(list = ls())
source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)

data.full$Time <- data.full$Time_cont
lmer_formula_model_uacr <- as.formula("FU_eGFR_epi_2021 ~ (Time + Time  *(BL_uacr_log2) + (1|Country) + (1+Time|PatID))")
fit.uacr <- lmer(lmer_formula_model_uacr,
                 data=data.full, REML=F, control=lmerControl(optimizer="bobyqa"))


# UPDATE prediction with BL value
data.full.t0 <- data.full[data.full$Time==0,]
data.full.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=fit.uacr, newdata = data.full.t0, timeVar = "Time", idVar="PatID", idVar2="Country",
                         times = unique(data.full$Time_cat)[-1], all_times=F)
data.full$Time <- round(data.full$Time,0)
data.sens <- full_join(data.full, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))

# ================== Internal validation  =====================================

# ----  Save final model and extract components for shiny ----
saveRDS(fit.uacr, here::here(out.path, "model_uacr","predmodel_lmerObject_model_uacr.rds"))

# ----  Residual check -----
check_residuals(fit.uacr, path = here::here(out.path, "model_uacr","fig_residual_analysis_model_uacr.tiff"))


# ====================== Internal-external validation =========================

set.seed(123)
j=1; b=nboot

# Cross-validated predictions
data.full$fold <- as.numeric(data.full$Country)
res <- intext_crossvalidate(data.full, NA, mod_formula=lmer_formula_model_uacr, return_preds = T, sa=T)
df.preds <- res$pred

# Cross-validated performance measures
plan(multisession, gc=T, workers=detectCores()*.6)
res_cv_boot <- future_lapply(1:b, function(x)  intext_crossvalidate(draw_bootstrap_sample(data.full),x, mod_formula=lmer_formula_model_uacr, return_preds = F, sa=T), future.seed = T)
plan(sequential)
df.stats <- data.frame(do.call(rbind, lapply(res_cv_boot, `[[`, 1)))


# =============== Table: Model discrimination, precision and fit =====================
tbl_performance <- list()
tbl_performance[[1]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], mean, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[2]]<- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], quantile, 0.05, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[3]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], quantile, 0.95, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))
df.tmp <- data.frame(do.call(cbind, tbl_performance))
tbl_performance[[4]] <-df.tmp %>%
  data.frame() %>%
  mutate(Time=Time,
         N=round_0(Nobs,3),
         N.ci=paste0("(",round_0(Nobs.1,3),", ",round_0(Nobs.2,3),")"),
         R2=round_0(R2,3),
         R2.ci = paste0("(",round_0(R2.1,3),", ",round_0(R2.2,3),")"),
         C=round_0(C,3),
         C.ci = paste0("(",round_0(C.1,3),", ",round_0(C.2,3),")"),
         CS=round_0(CalbSlope,3),
         CS.ci = paste0("(",round_0(CalbSlope.1,3),", ",round_0(CalbSlope.2,3),")")) %>%
  dplyr::select(Time, N, N.ci, R2, R2.ci, C, C.ci, CS, CS.ci) 

names(tbl_performance) <- c("avg", "ci.lo", "ci.up", "reported")
write.xlsx(tbl_performance, here::here(out.path, "model_uacr", "tbl_perform_val_full_uacr.xlsx"), overwrite = T)

# =================== External validation ======================================

rm(list=ls())
source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)
risk_model <- readRDS(here::here(out.path, "model_uacr","predmodel_lmerObject_model_uacr.rds"))

data.diacore$Time <- data.diacore$Time_cont
data.diacore.t0 <- data.frame(data.diacore[data.diacore$Time_cat==0,])
data.diacore.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=risk_model, 
                         newdata = data.diacore.t0, 
                         cutpoint = slope_cutpoint,
                         timeVar = "Time", idVar="PatID", idVar2="Country",
                         times =sort(unique(data.diacore$Time_cat)), 
                         all_times=F)

# Summarize and prepare output
data.diacore$Time <- round(data.diacore$Time,0)
data.diacore.new <- full_join(data.diacore, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))


# -------------------- Calibration plot ---------------------------------------

for(t in 2:7){
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi_2021, yhat=data.diacore.new[data.diacore.new$Time==t,]$prior.pred,
                        cohort="extval", time=t, save=T, folder="model_uacr", out.path = out.path, type="preUp_model_uacr")
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi_2021, yhat=data.diacore.new[data.diacore.new$Time==t,]$pred,
                        cohort="extval", time=t, save=T, folder="model_uacr", out.path = out.path, type="postUp_model_uacr")
}

ggplot(data.diacore.new[data.diacore.new$Time_cat <=5,], aes(x=pred, y=FU_eGFR_epi_2021, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Predicted eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = expression(paste("Observed eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=22)) 
ggsave(here::here(out.path, "model_uacr", "fig_calibration_all_diacore_post_model_uacr.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')


df <- melt(data.diacore.new[data.diacore.new$Time_cat <=5 ,c("Time_cont","FU_eGFR_epi_2021", "pred","prior.pred")], id.vars = c("Time_cont", "FU_eGFR_epi_2021"))
levels(df$variable) <- c("Post-update", "Pre-update")
df$variable <- relevel(df$variable, "Pre-update")
df <- df[df$Time_cont>0,]

ggplot(df, aes(x=value, y=FU_eGFR_epi_2021, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Predicted eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = expression(paste("Observed eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=16), strip.text = element_text(size = 18)) +
  facet_wrap(~variable)
ggsave(here::here(out.path, "model_uacr", "fig_calibration_all_diacore_prepost_model_uacr.tiff"),  width=12, height=6, device='tiff', dpi=350, compression = 'lzw')



# ------------------ Table: Performance measures -------------------------------

b=nboot
plan(multisession, gc=T, workers=detectCores()*.6)
res_ext_boot <- future_lapply(1:b, function(x){
  data.boot <- data.frame(draw_bootstrap_sample(data.diacore, cv=F))
  data.boot.t0 <- data.boot[data.boot$Time_cat==0,]
  res <- update_PredByBase(lmerObject=risk_model,
                           newdata=data.boot.t0,
                           cutpoint = slope_cutpoint,
                           timeVar = "Time", idVar="PatID", idVar2="Country",
                           times =seq(0,8,1), 
                           all_times=F)
  data.boot$Time <- data.boot$Time_cat
  data.boot.new <- full_join(data.boot, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], 
                             by=c("PatID", "Time"))
  
  data.boot.list <- split(data.boot.new, as.factor(data.boot.new$Time)) 
  res <- lapply(data.boot.list, function(x) eval_preds(pred=x$pred, obs=x$FU_eGFR_epi_2021, lmerObject=risk_model))
  df.res <- data.frame("Boot"=x,"Time"=as.numeric(names(res)), do.call(rbind, res))
  
  return(df.res)
}
, future.seed = T)
plan(sequential)

df.stats <- data.frame(do.call(rbind, res_ext_boot))



# Summarize and prepare output
tbl_performance <- list()
tbl_performance[[1]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], mean, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[2]]<- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], quantile, 0.05, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[3]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], quantile, 0.95, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))
df.tmp <- data.frame(do.call(cbind, tbl_performance))
tbl_performance[[4]] <-df.tmp %>%
  data.frame() %>%
  mutate(Time=Time,
         N=round_0(Nobs,3),
         N.ci=paste0("(",round_0(Nobs.1,3),", ",round_0(Nobs.2,3),")"),
         R2=round_0(R2,3),
         R2.ci = paste0("(",round_0(R2.1,3),", ",round_0(R2.2,3),")"),
         C=round_0(C,3),
         C.ci = paste0("(",round_0(C.1,3),", ",round_0(C.2,3),")"),
         CS=round_0(CalbSlope,3),
         CS.ci = paste0("(",round_0(CalbSlope.1,3),", ",round_0(CalbSlope.2,3),")")) %>%
  dplyr::select(Time, N, N.ci, R2, R2.ci, C, C.ci, CS, CS.ci) 

names(tbl_performance) <- c("avg", "ci.lo", "ci.up", "reported")
write.xlsx(tbl_performance, here::here(out.path, "model_uacr", "tbl_perform_extval_full_model_uacr.xlsx"), overwrite = T)

# ================= Table 2: Model coefficients ================================
df <- data.frame(variable=rownames(summary(risk_model)$coefficients),
                 effect=summary(risk_model)$coefficients[,1],
                 lower=summary(risk_model)$coefficients[,1] - 1.96*summary(risk_model)$coefficients[,2],
                 upper=summary(risk_model)$coefficients[,1] + 1.96* summary(risk_model)$coefficients[,2]) 
df$variable <- factor(df$variable, levels = df$variable[length(df$variable):1])
df$Group <- ifelse(str_detect(df$variable, "Time:"), "Main effect", "Interaction effect")
df <- mutate(df, Group = fct_rev(Group)) 
df$variable <- as.factor(str_replace(df$variable, "Time:", ""))
df$variable <- factor(df$variable, levels=c("(Intercept)","Time","BL_uacr_log2"))
levels(df$variable) <- list("Intercept"="(Intercept)","Time"="Time", BMI="BL_bmi", "Smoking (ever)"="BL_smoking1", MAP="BL_map", HbA1c="BL_hba1c",
                            "Serum Chol."="BL_serumchol","Hemoglobin"="BL_hemo", "log2 UACR"="BL_uacr_log2", "Glucose-low. Med."="BL_med_dm1", "Blood pressure-low. Med."="BL_med_bp1", "Lipid-low. Med."="BL_med_lipid1")


df[,2:4] <- apply(df[,2:4],2, round, digits=3)
df$name <- rownames(df)
df$Variable <- c(rep("Constant",2), rep(c("log2UACR"),2))
df$CI <- paste0("(",df$lower,", " ,df$upper,")")
tbl_fixeff <- cbind(df[!str_detect(df$name, "Time"),c(7,2,8)], df[str_detect(df$name, "Time"),c(2,8)])
colnames(tbl_fixeff) <- c("Variable", "Effect.Baseline", "95% CI", "Effect.Slope", "95% CI")
write.xlsx(tbl_fixeff, here::here(out.path, "model_uacr", "tbl_coeff_unstd_model_uacr.xlsx"), overwrite = TRUE)


# ===================== Forest plot of standardized coefficients ===============
data.full$Time <- data.full$Time_cont
data.scaled <- data.full %>%mutate_at(vars(Time, BL_bmi, BL_map, BL_hba1c, 
                                           BL_serumchol, BL_hemo, BL_uacr_log2, FU_eGFR_epi_2021), scale)
fit.scaled <- lmer(FU_eGFR_epi_2021 ~ (Time + Time  *(BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +                                                  
                                                        BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                   data=data.scaled, REML=F, control=lmerControl(optimizer="bobyqa"))

df <- data.frame(variable=rownames(summary(fit.scaled)$coefficients),
                 effect=summary(fit.scaled)$coefficients[,1],
                 lower=summary(fit.scaled)$coefficients[,1] - 1.96*summary(fit.scaled)$coefficients[,2],
                 upper=summary(fit.scaled)$coefficients[,1] + 1.96* summary(fit.scaled)$coefficients[,2]) 
df$variable <- factor(df$variable, levels = df$variable[length(df$variable):1])
df$Group <- ifelse(str_detect(df$variable, "Time:"), "Interaction effect", "Main effect")
df <- mutate(df, Group = fct_rev(Group)) 
df$variable <- as.factor(str_replace(df$variable, "Time:", ""))
df$variable <- factor(df$variable, levels=c("(Intercept)","Time","BL_bmi","BL_smoking1","BL_map","BL_hba1c","BL_serumchol", "BL_hemo",
                                            "BL_uacr_log2", "BL_med_dm1","BL_med_bp1","BL_med_lipid1"))
levels(df$variable) <- list("Intercept"="(Intercept)","Time"="Time", BMI="BL_bmi", "Smoking (ever)"="BL_smoking1", MAP="BL_map", HbA1c="BL_hba1c",
                            "Serum Chol."="BL_serumchol","Hemoglobin"="BL_hemo", "log2 UACR"="BL_uacr_log2", "Glucose-low. Med."="BL_med_dm1", "Blood pressure-low. Med."="BL_med_bp1", "Lipid-low. Med."="BL_med_lipid1")

df[,2:4] <- apply(df[,2:4],2, round, digits=3)
df$name <- rownames(df)
df$Variable <- c(rep("Constant",2), rep(c("BMI", "Smoking", "MAP", "HbA1c","Serum chol.",
                                          "Hemoglobin", "log2UACR", "GL Med.", "BPL Med.","LL Med."),2))
df$CI <- paste0("(",df$lower,", " ,df$upper,")")
tbl_fixeff <- cbind(df[!str_detect(df$name, "Time"),c(7,2,8)], df[str_detect(df$name, "Time"),c(2,8)])
colnames(tbl_fixeff) <- c("Variable", "Effect.Baseline", "95% CI", "Effect.Slope", "95% CI")
write.xlsx(tbl_fixeff, here::here(out.path, "model_uacr", "tbl_coeff_std_model_uacr.xlsx"), overwrite = TRUE)

ggplot(data=df[-1,], aes(y=reorder(variable,desc(variable)), x=effect, xmin=lower, xmax=upper)) +
  geom_point(size=2, shape=1) + 
  geom_errorbarh(height=.25) +
  scale_y_discrete("", labels=c("HbA1c" = expression(HbA[1][c]), "log2 UACR"=expression(paste(log[2], ' UACR')))) +
  scale_x_continuous("Standardized Effect", limits = c(-0.4,0.4), breaks=seq(-0.4,0.4,0.1)) +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  facet_wrap(~Group) +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave(here::here(out.path, "model_uacr", "plot_forest_standardized_model_uacr.tiff"),  width=8, height=4, device='tiff', dpi=350, compression = 'lzw')

## ----------- (3) model comparison ----------------


risk_model_main <- readRDS(here::here(out.path, "model_main","predmodel_lmerObject.rds"))
risk_model_clin <- readRDS(here::here(out.path, "model_clin","predmodel_lmerObject_model_clin.rds"))
risk_model_uacr <- readRDS(here::here(out.path, "model_uacr","predmodel_lmerObject_model_uacr.rds"))


anova(risk_model_clin,risk_model_uacr,risk_model_main)


## ------------ (4) RAAS therapy with PROVALID ------------

rm(list = ls())
source("scr/setup.R")
data.prov <- readRDS(here::here(out.path, "data", "cohort_data_provalid_raas.rds"))

# Model development and validation
data.prov$Time <- data.prov$Time_cont
lmer_formula_model_raas <- as.formula("FU_eGFR_epi_2021 ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_med_dm + BL_med_bp + BL_med_lipid + BL_raas) + (1|Country) + (1+Time|PatID))")
fit.raas <- lmer(lmer_formula_model_raas,
                 data=data.prov, REML=F, control=lmerControl(optimizer="bobyqa"))


# UPDATE prediction with BL value
data.prov.t0 <- data.prov[data.prov$Time==0,]
data.prov.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=fit.raas, newdata = data.prov.t0, timeVar = "Time", idVar="PatID", idVar2="Country",
                         times = unique(data.prov$Time_cat)[-1], all_times=F)
data.prov$Time <- round(data.prov$Time,0)
data.sens <- full_join(data.prov, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))

# ================== Internal validation  =====================================

# ----  Save final model and extract components for shiny ----
saveRDS(fit.raas, here::here(out.path, "model_raas","predmodel_lmerObject_model_raas.rds"))

# ----  Residual check -----
check_residuals(fit.raas, path = here::here(out.path, "model_raas","fig_residual_analysis_model_raas.tiff"))


# ====================== Internal-external validation =========================

set.seed(123)
j=1; b=nboot

# Cross-validated predictions
data.prov$Country <- as.factor(data.prov$Country)
data.prov$fold <- as.numeric(as.factor(data.prov$Country))
res <- intext_crossvalidate(data.prov, NA, mod_formula=lmer_formula_model_raas, return_preds = T, sa=T)
df.preds <- res$pred

# Cross-validated performance measures
plan(multisession, gc=T, workers=detectCores()*.6)
res_cv_boot <- future_lapply(1:b, function(x)  intext_crossvalidate(draw_bootstrap_sample(data.prov),x, mod_formula=lmer_formula_model_raas, return_preds = F, sa=T), future.seed = T)
plan(sequential)
df.stats <- data.frame(do.call(rbind, lapply(res_cv_boot, `[[`, 1)))


# =============== Table: Model discrimination, precision and fit =====================
tbl_performance <- list()
tbl_performance[[1]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], mean, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[2]]<- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], quantile, 0.05, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[3]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[5:13], quantile, 0.95, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))
df.tmp <- data.frame(do.call(cbind, tbl_performance))
tbl_performance[[4]] <-df.tmp %>%
  data.frame() %>%
  mutate(Time=Time,
         N=round_0(Nobs,3),
         N.ci=paste0("(",round_0(Nobs.1,3),", ",round_0(Nobs.2,3),")"),
         R2=round_0(R2,3),
         R2.ci = paste0("(",round_0(R2.1,3),", ",round_0(R2.2,3),")"),
         C=round_0(C,3),
         C.ci = paste0("(",round_0(C.1,3),", ",round_0(C.2,3),")"),
         CS=round_0(CalbSlope,3),
         CS.ci = paste0("(",round_0(CalbSlope.1,3),", ",round_0(CalbSlope.2,3),")")) %>%
  dplyr::select(Time, N, N.ci, R2, R2.ci, C, C.ci, CS, CS.ci) 

names(tbl_performance) <- c("avg", "ci.lo", "ci.up", "reported")
write.xlsx(tbl_performance, here::here(out.path, "model_raas", "tbl_perform_val_full_sa.xlsx"), overwrite = T)


# ================= Table 2: Model coefficients ================================
risk_model <- readRDS(here::here(out.path, "model_raas","predmodel_lmerObject_model_raas.rds"))
df <- data.frame(variable=rownames(summary(risk_model)$coefficients),
                 effect=summary(risk_model)$coefficients[,1],
                 lower=summary(risk_model)$coefficients[,1] - 1.96*summary(risk_model)$coefficients[,2],
                 upper=summary(risk_model)$coefficients[,1] + 1.96* summary(risk_model)$coefficients[,2]) 
df$variable <- factor(df$variable, levels = df$variable[length(df$variable):1])
df$Group <- ifelse(str_detect(df$variable, "Time:"), "Main effect", "Interaction effect")
df <- mutate(df, Group = fct_rev(Group)) 
df$variable <- as.factor(str_replace(df$variable, "Time:", ""))
df$variable <- factor(df$variable, levels=c("(Intercept)","Time","BL_bmi","BL_smoking1","BL_map","BL_hba1c","BL_serumchol", "BL_hemo",
                                            "BL_uacr_log2", "BL_med_dm1","BL_med_bp1","BL_med_lipid1", "BL_raas"))
levels(df$variable) <- list("Intercept"="(Intercept)","Time"="Time", BMI="BL_bmi", "Smoking (ever)"="BL_smoking1", MAP="BL_map", HbA1c="BL_hba1c",
                            "Serum Chol."="BL_serumchol","Hemoglobin"="BL_hemo", "log2 UACR"="BL_uacr_log2", "Glucose-low. Med."="BL_med_dm1", "Blood pressure-low. Med."="BL_med_bp1", "Lipid-low. Med."="BL_med_lipid1")


df[,2:4] <- apply(df[,2:4],2, round, digits=3)
df$name <- rownames(df)
df$Variable <- c(rep("Constant",2), rep(c("Age","Sex","BMI", "Smoking", "MAP", "GL Med.", "BPL Med.","LL Med.", "RAAS"),2))
df$CI <- paste0("(",df$lower,", " ,df$upper,")")
tbl_fixeff <- cbind(df[!str_detect(df$name, "Time"),c(7,2,8)], df[str_detect(df$name, "Time"),c(2,8)])
colnames(tbl_fixeff) <- c("Variable", "Effect.Baseline", "95% CI", "Effect.Slope", "95% CI")
write.xlsx(tbl_fixeff, here::here(out.path, "model_raas", "tbl_coeff_unstd_model_raas.xlsx"), overwrite = TRUE)

