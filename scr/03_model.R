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

# Before update
plot_calibration_cont(yobs=data.preds[data.preds$Time==1,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==1,]$prior.pred, fit=fit.final, time=1)
plot_calibration_cont(yobs=data.preds[data.preds$Time==2,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==2,]$prior.pred, fit=fit.final, time=2)
plot_calibration_cont(yobs=data.preds[data.preds$Time==3,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==3,]$prior.pred, fit=fit.final, time=3)
plot_calibration_cont(yobs=data.preds[data.preds$Time==4,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==4,]$prior.pred, fit=fit.final, time=4)
plot_calibration_cont(yobs=data.preds[data.preds$Time==5,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==5,]$prior.pred, fit=fit.final, time=5)

# After update
plot_calibration_cont(yobs=data.preds[data.preds$Time==1,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==1,]$pred, fit=fit.final, time=1)
plot_calibration_cont(yobs=data.preds[data.preds$Time==2,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==2,]$pred, fit=fit.final, time=2)
plot_calibration_cont(yobs=data.preds[data.preds$Time==3,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==3,]$pred, fit=fit.final, time=3)
plot_calibration_cont(yobs=data.preds[data.preds$Time==4,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==4,]$pred, fit=fit.final, time=4)
plot_calibration_cont(yobs=data.preds[data.preds$Time==5,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==5,]$pred, fit=fit.final, time=5)
plot_calibration_cont(yobs=data.preds[data.preds$Time==6,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==6,]$pred, fit=fit.final, time=6)
plot_calibration_cont(yobs=data.preds[data.preds$Time==7,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==7,]$pred, fit=fit.final, time=7)
plot_calibration_cont(yobs=data.preds[data.preds$Time==8,]$FU_eGFR_epi, yhat=data.preds[data.preds$Time==8,]$pred, fit=fit.final, time=8)


# ----- Final model
summary(fit.final)
coef(fit.final)$Country[,1]
print(coef(fit.final)$Country[,1])
hist(coef(fit.final)$PatID[,"(Intercept)"])
summary(coef(fit.final)$PatID[,"(Intercept)"])
hist(coef(fit.final)$PatID[,"Time"])
summary(coef(fit.final)$PatID[,"Time"])
summary(fit.final)$coefficients
data.full$resid <- residuals(fit.final)
#saveRDS(fit.final, paste0(out.path,"riskpred_model.rds"))
#saveRDS(fit.final, paste0(shiny.path,"/riskpred_model.rds"))



################################################################################
# ------------------ Internal validation (Model accuracy) ----------------------
################################################################################
# Fitted values vs standardized residuals
plot(fit.final, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

# Linearity in each variables
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

# Normality of residuals
plot(qqnorm(residuals(fit.final)))
points(qqline(residuals(fit.final)))


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
  data.test.new <- full_join(data.test, res$Pred[,c("PatID", "Time","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")], by=c("PatID", "Time"))
  data.test.new$fold <- i
  data.test.new$Country <- data.test$Country[1]
  data.test.new$Cohort <- data.test$Cohort[1]
  df.pred[[i]] <- data.test.new
  df.re[[i]] <- data.frame("PatID"= data.test.t0$PatID, "Fold"= i, res$RE.est)
  
  data.test.list <- split(data.test.new, as.factor(data.test.new$Time)) 
  res <- lapply(data.test.list, function(x) eval_preds(pred=x$pred, obs=x$FU_eGFR_epi, lmerObject=fit.lmer))
  df.res[[i]] <- data.frame("Fold"=i,"Time"=unique(data.full$Time),do.call(rbind, res))
}

# Concatenate
df.preds <- do.call(rbind, df.pred)
df.stats <- do.call(rbind, df.res)
df.reeff <- do.call(rbind, df.re)



################################################################################
# ---------------------------- "True" probability of progression ---------------
################################################################################

true.prog <- df.preds %>%
  mutate(Time=as.numeric(Time)) %>%
  group_by(PatID) %>%
  do(broom::tidy(lm(FU_eGFR_epi ~ Time,  data = .))) %>%
  filter(term %in% "Time") %>%
  dplyr::select(PatID, estimate) %>%
  `colnames<-`(c("PatID", "true.slope")) %>%
  mutate(true.prob = (true.slope <= slope_cutpoint)*1) %>%
  data.frame()

df.preds <- left_join(df.preds, true.prog, by="PatID")

