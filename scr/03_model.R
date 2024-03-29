#===============================================================================#
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Model building and validation
#===============================================================================#





# ============ Data preparation, IDA and general setup ==========================

# rm(list = ls())
# source("scr/setup.R")
# source("scr/01_dataprep.R", print.eval=F)
# 

# ====================== Model development =====================================

data.full$Time <- data.full$Time_cont
fit.final <- lmer(FU_eGFR_epi_2021 ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +                                                  
                BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                  data=data.full, REML=F, control=lmerControl(optimizer="bobyqa"))


# UPDATE prediction with BL value
data.full.t0 <- data.full[data.full$Time==0,]
data.full.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=fit.final, newdata = data.full.t0, timeVar = "Time", idVar="PatID", idVar2="Country",
                       times = unique(data.full$Time_cat)[-1], all_times=F)
data.full$Time <- round(data.full$Time,0)
data.preds <- full_join(data.full, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))


# ================ Predicted eGFR slope distribution ==========================
pred.slopes <- data.preds[duplicated(data.preds[,c("PatID", "pred.slope")]),]$pred.slope
distr.slopes <- quantile(pred.slopes, seq(0,1,0.0025))
dens.slopes <- data.frame("x"=density(pred.slopes)$x, "y"=density(pred.slopes)$y)
dens.slopes <- dens.slopes[seq(1,nrow(dens.slopes),10),]
distr.slopes <- data.frame("quantil"=names(distr.slopes), "Value"=distr.slopes, row.names=NULL) %>%
  mutate(quantil=as.numeric(str_remove(quantil, "%")))


# ================== Internal validation  =====================================

# ---- (1) Save final model and extract components for shiny ----
saveRDS(fit.final, here::here(out.path,"model_main","predmodel_lmerObject.rds"))

pred_components <- extractLMER(fit.final, idVar = "PatID")
pred_components$distr.slopes <- distr.slopes
pred_components$dens.slopes <- dens.slopes

saveRDS(pred_components, 
        file.path(shiny.path, "www/predmodel_shinyObject.rds"))
saveJSON(pred_components, 
         file.path(shiny.path,"www/predmodel_shinyObject.json"))

# Model equation 
# equatiomatic::extract_eq(fit.final)

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
  plot_calibration_cont(yobs=data.preds[data.preds$Time==t,]$FU_eGFR_epi_2021, yhat=data.preds[data.preds$Time==t,]$prior.pred,
                        cohort="dev", time=t, save=F, out.path = out.path, type="preUp")
  #After update
  plot_calibration_cont(yobs=data.preds[data.preds$Time==t,]$FU_eGFR_epi_2021, yhat=data.preds[data.preds$Time==t,]$pred, 
                        cohort="dev", time=t, save=F, out.path = out.path, type="postUp")
}


# ---- (4) Residual check -----

check_residuals(fit.final, path = here::here(out.path, "model_main","fig_residual_analysis.tiff"))


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
res <- intext_crossvalidate(data.full, NA, return_preds = T)
df.preds <- res$pred

# Cross-validated performance measures
plan(multisession, gc=T, workers=detectCores()*.6)
res_cv_boot <- future_lapply(1:b, function(x)  intext_crossvalidate(draw_bootstrap_sample(data.full),x, return_preds = F), future.seed = T)
plan(sequential)
df.stats <- data.frame(do.call(rbind, lapply(res_cv_boot, `[[`, 1)))




