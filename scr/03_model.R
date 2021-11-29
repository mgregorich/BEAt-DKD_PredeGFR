###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Model building and validation
###################################################



# -------- Data preparation, initial data analysis and general functions
# 
# rm(list = ls())
# source("scr/setup.R")
# source("scr/01_dataprep.R", print.eval=F)
# 

################################################################################
# ----------------------------- Model development ------------------------------
################################################################################

data.full$Time <- data.full$Time_cont
fit.final <- lmer(FU_eGFR_epi ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +                                                  
                BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                  data=data.full, REML=F, control=lmerControl(optimizer="bobyqa"))

# fit.final <- lme(fixed=FU_eGFR_epi ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +
#                                                        BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid)),
#                random=list(~1|Country,~1+Time|PatID), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
# summary(fit.final)$coefficients$fixed

# UPDATE prediction with BL value
data.full.t0 <- data.full[data.full$Time==0,]
data.full.t0$Country <- "Unknown"
res <- LongPred_ByBase_lmer(lmerObject=fit.final, newdata = data.full.t0, timeVar = "Time", idVar="PatID", idVar2="Country",
                       times = unique(data.full$Time_cat)[-1], all_times=F)
data.full$Time <- round(data.full$Time,0)
data.preds <- full_join(data.full, res$Pred[,c("PatID", "Time","prior.pred","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")], by=c("PatID", "Time"))

for(t in 1:7){
  #Before update
  plot_calibration_cont(yobs=data.preds[data.preds$Time==t,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==t,]$prior.pred, fit = fit.final,
                        cohort="dev", time=t, save=F, out.path = out.path, type="preUp")
  #After update
  plot_calibration_cont(yobs=data.preds[data.preds$Time==t,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==t,]$pred, fit = fit.final,
                        cohort="dev", time=t, save=F, out.path = out.path, type="postUp")}


# ----- Final model
# summary(fit.final)
# coef(fit.final)$Country[,1]
# print(coef(fit.final)$Country[,1])
# hist(coef(fit.final)$PatID[,"(Intercept)"])
# summary(coef(fit.final)$PatID[,"(Intercept)"])
# hist(coef(fit.final)$PatID[,"Time"])
# summary(coef(fit.final)$PatID[,"Time"])
# summary(fit.final)$coefficients
# data.full$resid <- residuals(fit.final)
#saveRDS(fit.final, paste0(out.path,"riskpred_model.rds"))
#saveRDS(fit.final, paste0(shiny.path,"/riskpred_model.rds"))



################################################################################
# ------------------ Internal validation (Model accuracy) ----------------------
################################################################################

# ---- Residual check
df_model[".stdresid"] <- resid(fit.final, type = "pearson")

# Fitted values vs standardized residuals
p1 <- ggplot(df_model, aes(.fitted, .resid)) + 
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_smooth(se=FALSE) +
  scale_x_continuous("Fitted Values") +
  scale_y_continuous("Standardized Residuals") +
  theme_bw() +
  theme(text=element_text(size=16))

# QQ Plot - Normality of residuals: qqnorm(residuals(fit.final))
p2 <- ggplot(df_model, aes(sample = .stdresid)) +
  geom_qq() +
  geom_qq_line() +
  scale_x_continuous("Sample Quantiles") +
  scale_y_continuous("Theoretical Quantiles") +
  theme_bw() +
  theme(text=element_text(size=16))

p3<-grid.arrange(p1,p2)
ggsave(paste0(out.path, "fig_residual_analysis.tiff"), plot=p3  ,width=10, height=6, device='tiff', dpi=350, compression = 'lzw')


# ---- Linearity in each variables
ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.final, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(se=T) +
  theme_bw() +
  facet_wrap(~time)

ggplot(data.frame(x1=data.full$BL_uacr_log2,residual=resid(fit.final, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(se=T) +
  theme_bw() +
  facet_wrap(~time)

ggplot(data.frame(x1=as.factor(data.full$Time),residual=resid(fit.final, scaled=TRUE)),
       aes(x=x1,y=residual)) +
  geom_point() +
  geom_boxplot() +
  # geom_smooth(se=T) +
  theme_bw() 



################################################################################
# -------------------- Internal-external validation ----------------------------
################################################################################

set.seed(123)
data.full$fold <- as.numeric(data.full$Country)
df.pred <- df.res <- df.re <- list()
f = length(unique(data.full$Country))
i=1

for(i in 1:f){
  print(i)
  data.train <- data.full[data.full$fold != i, ] 
  data.test <- data.full[data.full$fold == i, ] 
 
  fit.lmer <- lmer(FU_eGFR_epi ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +
                                                  BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                  data=data.train, REML=F, control=lmerControl(optimizer="bobyqa"))
  
  # Update prediction with BL value
  data.test.t0 <- data.test[data.test$Time==0,]
  data.test.t0$Country <- "Unknown"
  res <- LongPred_ByBase_lmer(lmerObject=fit.lmer, 
                         newdata = data.test.t0, 
                         cutpoint = slope_cutpoint,
                         timeVar = "Time", idVar="PatID", idVar2="Country",
                         times =unique(data.full$Time_cat)[-1], 
                         all_times=F)
  
  # Summarize and prepare output
  data.test$Time <- round(data.test$Time,0)
  data.test.new <- full_join(data.test, res$Pred[,c("PatID", "Time","prior.pred","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")], by=c("PatID", "Time"))
  data.test.new$fold <- i
  data.test.new$Country <- data.test$Country[1]
  data.test.new$Cohort <- data.test$Cohort[1]
  df.pred[[i]] <- data.test.new
  df.re[[i]] <- data.frame("PatID"= data.test.t0$PatID, "Fold"= i, res$RE.est)
  
  data.test.list <- split(data.test.new, as.factor(data.test.new$Time)) 
  res <- lapply(data.test.list, function(x) eval_preds(pred=x$pred, obs=x$FU_eGFR_epi, lmerObject=fit.lmer))
  df.res[[i]] <- data.frame("Fold"=i,"Time"=unique(data.full$Time),do.call(rbind, res))
}


# --- Concatenate
df.preds <- do.call(rbind, df.pred)
df.stats <- do.call(rbind, df.res)
df.reeff <- do.call(rbind, df.re)

