#===============================================================================#
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Model building and validation
#===============================================================================#





# ============ Data preparation, IDA and general setup ==========================

# rm(list = ls())
# source("scr/setup.R")
# source("scr/01_dataprep.R", print.eval=F)


# ====================== Model development =====================================

data.full$Time <- data.full$Time_cont
fit.final <- lmer(FU_eGFR_epi ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +                                                  
                BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                  data=data.full, REML=F, control=lmerControl(optimizer="bobyqa"))


# UPDATE prediction with BL value
data.full.t0 <- data.full[data.full$Time==0,]
data.full.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=fit.final, newdata = data.full.t0, timeVar = "Time", idVar="PatID", idVar2="Country",
                       times = unique(data.full$Time_cat)[-1], all_times=F)
data.full$Time <- round(data.full$Time,0)
data.preds <- full_join(data.full, res$Pred[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))


# ================== Internal validation  ======================

# ---- (1) Save final model and extract components for shiny ----
saveRDS(fit.final, paste0(out.path,"predmodel_lmerObject.rds"))

pred_components <- extractLMER(fit.final, idVar = "PatID")
saveRDS(pred_components, paste0(shiny.path,"/predmodel_shinyObject.rds"))


# ---- (2) Model parameter ----
# summary(fit.final)
# coef(fit.final)$Country[,1]
# print(coef(fit.final)$Country[,1])
# hist(coef(fit.final)$PatID[,"(Intercept)"])
# summary(coef(fit.final)$PatID[,"(Intercept)"])
# hist(coef(fit.final)$PatID[,"Time"])
# summary(coef(fit.final)$PatID[,"Time"])
# summary(fit.final)$coefficients
# data.full$resid <- residuals(fit.final)


# ---- (3) Calibration -----
for(t in 1:7){
  #Before update
  plot_calibration_cont(yobs=data.preds[data.preds$Time==t,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==t,]$prior.pred,
                        cohort="dev", time=t, save=F, out.path = out.path, type="preUp")
  #After update
  plot_calibration_cont(yobs=data.preds[data.preds$Time==t,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==t,]$pred, 
                        cohort="dev", time=t, save=F, out.path = out.path, type="postUp")}


# ---- (4) Residual check -----
df_model <- augment(fit.final)
df_model[".stdresid"] <- resid(fit.final, type = "pearson")

p1 <- ggplot(df_model, aes(.fitted, .resid)) + 
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_smooth(method = "gam",se=T, formula = y ~ x) +
  scale_x_continuous("Fitted Values") +
  scale_y_continuous("Standardized Residuals") +
  theme_bw() +
  theme(text=element_text(size=16))

p2 <- ggplot(df_model, aes(sample = .stdresid)) +
  geom_qq() +
  geom_qq_line() +
  scale_x_continuous("Sample Quantiles") +
  scale_y_continuous("Theoretical Quantiles") +
  theme_bw() +
  theme(text=element_text(size=16))

p3<-grid.arrange(p1,p2)
ggsave(paste0(out.path, "fig_residual_analysis.tiff"), plot=p3  ,width=10, height=6, device='tiff', dpi=350, compression = 'lzw')


# ---- (5) Linearity in each variables -----
ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.final, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(method = "gam",se=T, formula = y ~ x) +
  theme_bw() +
  facet_wrap(~time)

ggplot(data.frame(x1=data.full$BL_uacr_log2,residual=resid(fit.final, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(method = "gam",se=T, formula = y ~ x) +
  theme_bw() +
  facet_wrap(~time)

ggplot(data.frame(x1=as.factor(data.full$Time),residual=resid(fit.final, scaled=TRUE)),
       aes(x=x1,y=residual)) +
  geom_point() +
  geom_boxplot() +
  theme_bw() 



# ====================== Internal-external validation =========================

set.seed(123)
j=1; b=nboot

# Cross-validated predictions
data.full$fold <- as.numeric(data.full$Country)
res <- intext_crossvalidate(data.full,NA, return_preds = T)
df.preds <- res$pred

# Cross-validated performance measures
plan(multisession, gc=T, workers=detectCores()*.6)
res_cv_boot <- future_lapply(1:b, function(x)  intext_crossvalidate(draw_bootstrap_sample(data.full),x, return_preds = F), future.seed = T)
plan(sequential)
df.stats <- data.frame(do.call(rbind, lapply(res_cv_boot, `[[`, 1)))




