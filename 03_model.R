##################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Initial data analysis: data screening & cleaning
###################################################

rm(list=ls())
source("01_dataprep.R", print.eval=F)
source("02_IDA.R", print.eval=F)


# ----------- FUNCTIONS ------------
source("functions.R")


pacman::p_load(nlme, lmerTest, MuMIn, JMbayes, splines,rms, Hmisc, concreg, caret, MASS, performance)


# ----- Model building
pred.vars <- c("Time","BL_age", "BL_gender", "BL_smoking", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c", 
               "BL_serumchol", "BL_hemo", "BL_uacr_log2", "BL_med_dm", "BL_med_bp", "BL_med_lipid", "BL_med_lipid")


# ----- Non-linear interactions
# Fit model and check for non-linear behaviour in residuals vs independent variables
fit.init <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
                                                      BL_hba1c + BL_serumchol + BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid)), 
                random=list(~1+Time|PatID, ~1|Country), data=data.full, control=lmeControl(opt = "optim"), method = "ML")


ggplot(data.frame(x1=data.full$BL_uacr_log2,residual=resid(fit.init, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(se=T) +
  theme_bw() +
  facet_wrap(~time)

ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.init, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(se=T) +
  theme_bw() +
  facet_wrap(~time)

# Non-linear term for BL_uacr_log2
fit.tmp1 <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
                                                      BL_hba1c + BL_serumchol + BL_hemo + rcs(BL_uacr_log2,df=2) + BL_med_dm + BL_med_bp + BL_med_lipid)), 
                random=list(~1+Time|PatID, ~1|Country), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
anova(fit.init, fit.tmp1)


ggplot(data.frame(x1=data.full$BL_uacr_log2,residual=resid(fit.tmp1, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(se=T) +
  theme_bw() +
  facet_wrap(~time)

ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.tmp1, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(se=T) +
  theme_bw() +
  facet_wrap(~time)

# Non-linear term for BL_serumchol
fit.tmp2 <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
                                                     BL_hba1c + rcs(BL_serumchol,df=2) + BL_hemo + rcs(BL_uacr_log2,df=2) + BL_med_dm + BL_med_bp + BL_med_lipid)), 
               random=list(~1+Time|PatID, ~1|Country), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
anova(fit.tmp1, fit.tmp2)

ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.tmp2, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(se=T) +
  theme_bw() +
  facet_wrap(~time)

ggplot(data.frame(x1=data.full$BL_BL_uacr_log2,residual=resid(fit.tmp2, scaled=TRUE), time=data.full$Time),
       aes(x=x1,y=residual, color=time)) +
  geom_smooth(se=T) +
  theme_bw() +
  facet_wrap(~time)



# ----- Pairwise interactions
out.step <- stepAIC(fit.tmp2, 
                    direction = 'forward', 
                    scope=list(upper = ~ (Time + Time * (BL_bmi + BL_age + BL_gender + BL_smoking + 
                                                           (BL_bpsys + BL_bpdia)^2 + BL_hba1c + rcs(BL_serumchol, df = 2) + 
                                                           BL_hemo + rcs(BL_uacr_log2, df = 2) + BL_med_dm + BL_med_bp + 
                                                           BL_med_lipid))))

# saveRDS(out.step, file = "out.step.rds")


# ----- Final model
fit.final <- out.step
summary(fit.final)
hist(fit.final$coefficients$random$PatID[,2])
summary(fit.final$coefficients$random$PatID[,2])
fit.final$coefficients$fixed



# ------- Internal validation (Model accuracy)

# Fitted values vs standardized residuals
plot(fit.final, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

# Linearity in each variables

ggplot(data.frame(x1=data.full$BL_serumchol,residual=resid(fit.final, scaled=TRUE), time=data.full$Time),
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
qqnorm(residuals(fit.final))




# -------------------- Internal-external validation --------------------------------
set.seed(123)
data.full$fold <- as.numeric(data.full$Country)
df.pred <- df.res <- df.re <- list()
k = length(unique(data.full$Country))
i=1

for(i in 1:k){
  print(i)
  data.train <- data.full[data.full$fold != i, ] 
  data.test <- data.full[data.full$fold == i, ] 
  
  # fit.lme <- lme(fixed=formula(fit.final), 
  #                random=list(~1+Time|PatID, ~1|Country), data=data.train, control=lmeControl(opt = "optim"), method = "ML")
  # fit.lme <- lme(fixed=formula(fit.final), 
  #                random=~1+Time|PatID, data=data.train, control=lmeControl(opt = "optim"), method = "ML")
  
  fit.lme <- lme(fixed=FU_eGFR_epi ~ Time + Time * BL_bmi, 
                 random=list(~1+Time|PatID, ~1|Country), data=data.train, control=lmeControl(opt = "optim"), method = "ML")
  
  # Update prediction with BL value
  data.test.t0 <- data.test[data.test$Time==0,]
  res <- LongPred_ByBase(fit.lme, newdata = data.test.t0, timeVar = "Time", idVar="PatID", times = unique(data.test$Time))
  
  #pred <- IndvPred_lme(fit.lme, newdata = data.test.t0, timeVar = "Time", times = unique(data.test$Time), M = 500, return_data = T, interval = "prediction")
  
  pred <- res$Pred
  data.test.new <- full_join(data.test, pred[,c("PatID", "Time","prob.prog" ,"pred", "low", "upp")], by=c("PatID", "Time"))
  df.pred[[i]] <- data.test.new
  df.re[[i]] <- data.frame("PatID"= data.test.t0$PatID, "Fold"= i, res$RE.est)
  
  data.test.list <- split(data.test.new, as.factor(data.test.new$Time)) 
  res <- lapply(data.test.list, function(x) eval_preds(pred=x$pred, obs=x$FU_eGFR_epi, fit = fit.lme))
  df.res[[i]] <- data.frame("Fold"=i,"Time"=unique(data.test$Time),do.call(rbind, res))
}

df.preds <- do.call(rbind, df.pred)
df.stats <- do.call(rbind, df.res)
df.reeff <- do.call(rbind, df.re)


# Calibration plot per time t
plot_calibration(yobs=df.preds[df.preds$Time==1,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==1,]$pred, time=1)
plot_calibration(yobs=df.preds[df.preds$Time==2,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==2,]$pred, time=2)
plot_calibration(yobs=df.preds[df.preds$Time==3,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==3,]$pred, time=3)
plot_calibration(yobs=df.preds[df.preds$Time==4,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==4,]$pred, time=4)
plot_calibration(yobs=df.preds[df.preds$Time==5,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==5,]$pred, time=5)


# Plot subject-specific trajectories
which.max(abs(fit.final$coefficients$random$PatID[,2]))

set.seed(666)
df.melt <- melt(df.preds[,c("PatID","Time", "FU_eGFR_epi", "pred", "low", "upp")], id.vars = c("PatID","Time"))
df.melt.small <- df.melt[df.melt$PatID %in% c(sample(unique(df.melt$PatID),8),"PV18562"),]

ggplot(data =df.melt.small, aes(x = Time, y = value, group = variable, col=variable)) +
  geom_point() +
  geom_line(data=df.melt.small[!is.na(df.melt.small$value),], aes(linetype=variable, color=variable)) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed")) +
  scale_color_manual(values=c("black", "blue", "red", "red")) +
  scale_x_continuous(breaks = seq(1,5,1), limits = c(0,5)) +
  ggtitle(df.melt$PatID[1]) +
  theme_bw() +
  facet_wrap(~PatID)



