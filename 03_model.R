##################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Initial data analysis: data screening & cleaning
###################################################

rm(list=ls())
source("01_dataprep.R", print.eval=F)
source("02_IDA.R", print.eval=F)


pacman::p_load(nlme, lmerTest, MuMIn, JMbayes, splines,rms, Hmisc, concreg, caret, MASS, performance)


# ---------------- FUNCTIONS --------------------------

c_index <- function(pred, obs){
  c.model <- concreg(data=data.frame(predicted=pred, observed=obs), observed~predicted, npar=TRUE)
  return(1-cindex(c.model))
}

eval_preds <- function(pred, obs, fit){
  df <- data.frame(pred=pred, obs=obs)
  df <- df[complete.cases(df),]
  r2 <- r2_nakagawa(fit.final)
  
  res <- data.frame(R2cond = r2[1],
                    R2marg = r2[2],
              RMSE = RMSE(df$pred, df$obs),
              MAE = MAE(df$pred, df$obs),
              CalbinLarge = mean(df$obs)- mean(df$pred),
              CalbSlope=lm(df$pred ~ df$obs)$coefficients[2],
              C = c_index(pred=df$pred, obs=df$obs))
  return(res)
}

plot_calibration <- function(yobs, yhat, time="Not specified!"){
  df <- data.frame(yobs=yobs, yhat=yhat)
  res <- round(eval_preds(pred=yhat, obs = yobs),2)
  
  # Plot
  plot_calibration_curve_ps <- ggplot(df, aes(x=yhat, y=yobs)) +
    geom_point() +
    scale_x_continuous(expand = c(0, 0), name = "Predicted", limits = c(0,150)) + 
    scale_y_continuous(expand = c(0, 0), name = "Observed", limits = c(0,150)) +
    geom_abline(intercept = 0) + 
    ggtitle(paste0("FU time = ", time)) +
    theme_bw() +
    annotate("text", x = 25, y = 125, label = paste0("R2 cond. = ",res$R2_conditional ,
                                                     "\nR2 marg. = ", res$R2_marginal, 
                                                     "\nC-index = ", res$C, 
                                                     "\nCalLarge = ", res$CalbinLarge, 
                                                     "\nCalSlope = ", res$CalbSlope))
  plot_calibration_curve_ps
  
  # filename <- paste0("plot_calcurve_", mod.method)
  # png(file = paste(out.path, filename,".png", sep = ""), width = 5*600, height = 5*600, units = "px", res = 600)
  # plot_calibration_curve_ps
  # dev.off()
}


# ------------------------ END ------------------------------



# ----------------  Model building -------------


pred.vars <- c("Time","BL_age", "BL_gender", "BL_smoking", "BL_bmi", "BL_bpsys", "BL_bpdia", "BL_hba1c", 
               "BL_serumchol", "BL_hemo", "BL_uacr_log2", "BL_med_dm", "BL_med_bp", "BL_med_lipid", "BL_med_lipid")

# --------- Non-linear interactions: BL_age (detected in IDA)
init.fit1 <- lme(fixed=FU_eGFR_epi ~ Time + BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
                   BL_hba1c + BL_serumchol + BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid, 
                 random=~Time+Time|PatID, data=data.full, control=lmeControl(opt = "optim"), method = "ML")


init.fit2 <- lme(fixed=FU_eGFR_epi ~ Time + BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
                  BL_hba1c + BL_serumchol + BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid, 
                random=~Time+Time|PatID, data=data.full, control=lmeControl(opt = "optim"), method = "ML")
anova(init.fit1, init.fit2)


# ----------- Pairwise interactions
init.fit1 <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
                BL_hba1c + BL_serumchol + BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid)), 
                 random=list(~1+Time|PatID, ~1|Country), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
init.fit2 <- lme(fixed=FU_eGFR_epi ~ (Time + Time * (BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
                   BL_hba1c + rcs(BL_serumchol,df=3) + BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid)), 
                 random=list(~1+Time|PatID, ~1|Country), data=data.full, control=lmeControl(opt = "optim"), method = "ML")
anova(init.fit1, init.fit2)


# ------ All pairwise interactions ??
# init.fit <- lme(fixed=FU_eGFR_epi ~ (Time + BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
#                                        BL_hba1c + BL_serumchol + BL_hemo + BL_uacr + BL_med_dm + BL_med_bp + BL_med_lipid)^2, 
#                 random=~Time+Time|PatID, data=data.full, control=lmeControl(opt = "optim"), method = "ML")
# out.step <- stepAIC(init.fit, direction = 'backward')
# saveRDS(out.step, file = "out.step.rds")


# ----- Final model
fit.final <- init.fit2
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


ggplot(data.frame(lev=hatvalues(fit.final),pearson=residuals(fit.final,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()




# -------------------- Internal-external validation --------------------------------
set.seed(123)
data.full$fold <- as.numeric(data.full$Country)
df.pred <- list()
df.res <- list()
k = length(unique(data.full$Country))
i=1

for(i in 1:k){
  print(i)
  data.train <- data.full[data.full$fold != i, ] 
  data.test <- data.full[data.full$fold == i, ] 
  
  fit.lme <- lme(fixed=FU_eGFR_epi ~ Time + BL_bmi + BL_age + BL_gender + BL_smoking + BL_bpsys + BL_bpdia + 
                   BL_hba1c + BL_serumchol + BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid, 
                 random=~Time+Time|PatID, data=data.train, control=lmeControl(opt = "optim"), method = "ML")
  summary(fit.lme)
  
  # Update prediction with BL value
  data.test.t0 <- data.test[data.test$Time==0,]
  pred <- IndvPred_lme(fit.lme, newdata = data.test.t0, timeVar = "Time", times = unique(data.test$Time), M = 500, return_data = T, interval = "prediction")
  data.test.new <- full_join(data.test, pred[,c("PatID", "Time", "pred", "low", "upp")], by=c("PatID", "Time"))
  df.pred[[i]] <- data.test.new
  
  data.test.list <- split(data.test.new, as.factor(data.test.new$Time)) 
  res <- lapply(data.test.list, function(x) eval_preds(pred=x$pred, obs=x$FU_eGFR_epi))
  df.res[[i]] <- data.frame("Fold"=i,"Time"=unique(data.test$Time),do.call(rbind, res))
}

df.preds <- do.call(rbind, df.pred)
df.stats <- do.call(rbind, df.res)


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



