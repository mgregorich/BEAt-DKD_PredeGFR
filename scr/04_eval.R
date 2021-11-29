##################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Model evaluation
###################################################


################################################################################
# --------------- Model discrimination, precision and fit ----------------------
################################################################################

tbl_tmp <- lapply(unique(data.full$Time), function(x) eval_preds(pred=df.preds[df.preds$Time==x,]$pred, 
                                                                 obs=df.preds[df.preds$Time==x,]$FU_eGFR_epi, 
                                                                 lmerObject = fit.final))
time.perf <- data.frame("Time"=seq(0,8,1),round(do.call(rbind, tbl_tmp),3))
# Overall performance
ov.perf <- eval_preds(df.preds[!df.preds$Time==0,]$pred, df.preds[!df.preds$Time==0,]$FU_eGFR_epi, lmerObject = fit.final)
tbl_performance <- list("overall"=ov.perf, "time-specific"=time.perf)

write.xlsx(tbl_performance, paste0(out.path, "tbl_perform_val.xlsx"), 
           overwrite = TRUE)



################################################################################
# ----------------------------- Calibration ------------------------------------
################################################################################

for(t in 2:7){
  # Before update
  plot_calibration_cont(yobs=df.preds[df.preds$Time==t,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==t,]$prior.pred, fit = fit.final, time=t, save=T, out.path = out.path, cohort="dev", type="preUp")
  # After update
  plot_calibration_cont(yobs=df.preds[df.preds$Time==t,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==t,]$pred, fit = fit.final, time=t, save=T, out.path = out.path, cohort="dev", type="postUp")
}

ggplot(df.preds, aes(x=pred, y=FU_eGFR_epi, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = "Predicted eGFR", limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = "Observed EGFR", limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=22)) 
ggsave(paste0(out.path, "fig_calibration_all_dev_post.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')


# ------ Probability of progression
pred_prob = df.preds[df.preds$Time==0,]$prob.prog
true_prob = df.preds[df.preds$Time==0,]$true.prob
plot_calibration_bin(pred=pred_prob, true=as.factor(true_prob), out.path = out.path, save=F)
Brier <- mean((pred_prob-true_prob)^2)
C.index <- somers2(pred_prob, true_prob)["C"]
Deviance<- -2*sum(true_prob*log(pred_prob) + (1-true_prob)*log(1-pred_prob) )
CS <- lm(true_prob ~ pred_prob)$coefficients[2]
print(paste0("Brier = ", round(Brier,2), ", C statistic = ", round(C.index,2), ", calibration slope = ",round(CS,2)," and deviance = ", round(Deviance,2)))



################################################################################
# ---------------------------- Model fit ---------------------------------------
################################################################################

# ----- Plot subject-specific trajectories
which.max(abs(coef(fit.final)$PatID[,"Time"]))

set.seed(666)
df.melt <- melt(df.preds[,c("PatID","Time", "FU_eGFR_epi","pred", "pred.low", "pred.upp")], id.vars = c("PatID","Time", "pred.low", "pred.upp"))
df.melt.small <- df.melt[df.melt$PatID %in% c(sample(unique(df.melt$PatID),16)),]
df.melt.small$Time <- as.numeric(df.melt.small$Time)

ggplot(data =df.melt.small, aes(x = Time, y = value,  col=variable)) +
  geom_ribbon(aes(x=Time, ymax=pred.upp, ymin=pred.low), fill="blue", alpha=.05) +
  geom_point(size=1) +
  geom_line(data=df.melt.small[!is.na(df.melt.small$value),], aes(color=variable)) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed")) +
  scale_color_manual(values=c("black", "blue", "red", "red"), labels=c("Observed", "Predicted", "95% CI", "95% CI")) +
  scale_x_continuous(expand=c(0,0),breaks = seq(0,8,1), limits = c(0,8)) +
  scale_y_continuous(expand=c(0,0), "eGFR") +
  ggtitle("") +
  theme_bw() +
  facet_wrap(~PatID) +
  theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=16))
ggsave(paste0(out.path, "fig_indvPred_eGFR_dev.png"),width=10, height=6)

df.melt.small <- df.melt[df.melt$PatID %in% "PV14750",]
ggplot(data =df.melt.small, aes(x = Time, y = value,  col=variable)) +
  geom_ribbon(aes(x=Time, ymax=pred.upp, ymin=pred.low), fill="blue", alpha=.05) +
  geom_point(size=1) +
  geom_line(data=df.melt.small[!is.na(df.melt.small$value),], aes(color=variable)) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed")) +
  scale_color_manual(values=c("black", "blue", "red", "red"), labels=c("Observed", "Predicted", "95% CI", "95% CI")) +
  scale_x_continuous(expand=c(0,0),breaks = seq(0,8,1), limits = c(0,8)) +
  scale_y_continuous(expand=c(0,0), "eGFR") +
  ggtitle("") +
  theme_bw() +
  facet_wrap(~PatID) +
  theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=16))



# ---- subject specific trajectories, comparing updated and no-update predictions
df.melt_0 <- melt(df.preds[df.preds$Time==0,c("PatID", "Time", "pred")], id.vars = c("PatID","Time"))
df.melt.small_0 <- df.melt_0[df.melt_0$PatID %in% unique(df.melt.small$PatID), ]
df.melt.small_0$variable = "pred_0"
df.melt.small_0$Time <- as.numeric(df.melt.small_0$Time)
df.melt.small_0 <- rbind(df.melt.small[, c("PatID", "Time", "variable", "value")], 
                         df.melt.small_0)

ggplot(data =df.melt.small_0, aes(x = Time, y = value,  col=variable)) +
  geom_point(size=1) +
  geom_line(data=df.melt.small_0[!is.na(df.melt.small_0$value),], aes(color=variable)) +
  scale_linetype_manual(values=c("solid", "solid", "dashed")) +
  scale_color_manual(values=c("black", "blue", "red"), labels=c("Observed", "Predicted", "No Update")) +
  scale_x_continuous(expand=c(0,0),breaks = seq(0,7,1), limits = c(0,7)) +
  scale_y_continuous(expand=c(0,0), "eGFR") +
  ggtitle("") +
  theme_bw() +
  facet_wrap(~PatID) +
  theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=16))
ggsave(paste0(out.path, "fig_indvPred_eGFR_dev_noUpdate.png"),width=10, height=6)

# ---- Mean country trajectories
df.preds %>%
  group_by(Country, Time) %>%
  summarise(pm = median(FU_eGFR_epi, na.rm=T)) %>%
  ggplot(aes(x = Time, y = pm, group = Country, col=Country)) +
  geom_point() +
  geom_line(aes(color=Country)) +
  theme_bw()


# ---- Longitudinal trajectory of observed eGFR per country
df.preds$Time <- as.factor(df.preds$Time)
ggplot(df.preds,aes(x = Time, y = FU_eGFR_epi, fill=Country))  +
  geom_boxplot(aes(x=Time, fill=Country), outlier.shape = NA) +
  stat_boxplot(aes(fill=Country),geom ='errorbar') +
  scale_y_continuous("Observed eGFR") +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(paste0(out.path, "fig_longiObs_eGFR.png"),width=12, height=6)


# Longitudinal trajectory of predicted eGFR per cohort
df.tmp <- data.frame(rbind(df.preds[,c("Time_cat", "FU_eGFR_epi", "Cohort")], data.diacore[,c("Time_cat", "FU_eGFR_epi", "Cohort")]))
df.tmp$Cohort <- as.factor(df.tmp$Cohort)
df.tmp$Time_cat <- as.factor(df.tmp$Time_cat)
df.tmp <- df.tmp[!is.na(df.tmp$Time_cat),]
ggplot(df.tmp,aes(x = Time_cat, y = FU_eGFR_epi, fill=Cohort))  +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom ='errorbar') +
  scale_y_continuous(expand=c(0,0),"Observed eGFR", lim=c(0,150), breaks=seq(0,150,25)) +
  scale_x_discrete("Time") +
  scale_fill_manual(values=c("skyblue1", "royalblue1", "lightgreen"), labels=c("GCKD", "PROVALID", "DIACORE")) +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(paste0(out.path, "fig_longiobs_eGFR_cohort.png"),width=10, height=6)



# ---- Longitudinal trajectory of predicted eGFR per country
ggplot(df.preds,aes(x = Time, y = pred, fill=Country))  +
  geom_boxplot(aes(x=Time, fill=Country), outlier.shape = NA) +
  stat_boxplot(aes(fill=Country),geom ='errorbar') +
  scale_y_continuous("Predicted eGFR") +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(paste0(out.path, "fig_longipred_eGFR.png"),width=12, height=6)

# ---- Longitudinal trajectory of predicted eGFR per cohort
df.preds$Cohort <- as.factor(df.preds$Cohort)
ggplot(df.preds,aes(x = Time, y = pred, fill=Cohort))  +
  geom_boxplot(aes(x=Time, fill=Cohort), outlier.shape = NA) +
  stat_boxplot(aes(fill=Cohort),geom ='errorbar') +
  scale_y_continuous("Predicted eGFR") +
  scale_fill_manual(values=c("orange1", "royalblue1"), labels=c("GCKD", "PROVALID")) +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(paste0(out.path, "fig_longipred_eGFR_cohort_dev.png"),width=10, height=6)


# ---- Forest plot of model coefficients
df <- data.frame(variable=rownames(summary(fit.final)$coefficients),
                 effect=summary(fit.final)$coefficients[,1],
                 lower=summary(fit.final)$coefficients[,1] - 1.96*summary(fit.final)$coefficients[,2],
                 upper=summary(fit.final)$coefficients[,1] + 1.96* summary(fit.final)$coefficients[,2])
df$variable <- factor(df$variable, levels = df$variable[length(df$variable):1])

head(df)
ggplot(data=df[-1,], aes(y=variable, x=effect, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_discrete("", labels=c("(Intercept)"="Intercept","BL_age"="Age", "BL_sex1"="Sex", "BL_bmi"="BMI", "BL_smoking1"="Smoking", "BL_map"="MAP", "BL_hba1c_perc"="Hba1c","BL_serumchol"= "Serum chol.",
                                "BL_hemo"="Hemoglobin", "BL_uacr_log2"="log2UACR", "BL_med_dm1"="GL Med.", "BL_med_bp1"="BPL Med.", "BL_med_lipid1"="LL Med.",
                                "Time:BL_age"="Time:Age", "Time:BL_sex1"="Time:Sex", "Time:BL_bmi"="Time:BMI", "Time:BL_smoking1"="Time:Smoking", "Time:BL_map"="Time:MAP",
                                "Time:BL_hba1c_perc"="Time:Hba1c","Time:BL_serumchol"= "Time:Serum chol.", "Time:BL_hemo"="Time:Hemoglobin", "Time:BL_uacr_log2"="Time:log2UACR", 
                                "Time:BL_med_dm1"="Time:GL Med.", "Time:BL_med_bp1"="Time:BPL Med.","Time:BL_med_lipid1"="Time:LL Med." )) +
  scale_x_continuous("Effect") +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave(paste0(out.path, "forest_plot.png"),width=6, height=6)

# Table of model coefficients
df[,2:4] <- apply(df[,2:4],2, round, digits=3)
df$name <- c(rep("Constant",2), rep(c("Age", "Sex", "BMI", "Smoking", "MAP", "Hba1c","Serum chol.",
                                      "Hemoglobin", "log2UACR", "GL Med.", "BPL Med.","LL Med."),2))
df$CI <- paste0("(",df$lower,", " ,df$upper,")")
tbl_fixeff <- cbind(df[!str_detect(df$variable, "Time"),c(5,2,6)], df[str_detect(df$variable, "Time"),c(2,6)])
write.xlsx(tbl_fixeff, paste0(out.path, "tbl_fixeff.xlsx"), overwrite = TRUE)


# ---- Individual risk predictions
df.preds.t0 <- df.preds[df.preds$Time==0,]

ggplot(df.preds.t0, aes(x=prob.prog, fill=Cohort)) +
  geom_density(alpha=0.4) +
  scale_fill_manual(values=c("orange1", "royalblue1"), labels=c("GCKD", "PROVALID")) +
  scale_x_continuous("Probability of progression") +
  scale_y_continuous("Density") +
  theme_bw() +
  theme(text = element_text(size=18), legend.position = "bottom", legend.title = element_blank()) 
ggsave(paste0(out.path, "fig_prob_progression_",abs(slope_cutpoint),"_dev.png"),width=8, height=6)

summary(df.preds.t0$prob.prog)
summary(df.preds.t0[df.preds.t0$Cohort==0,]$prob.prog)
summary(df.preds.t0[df.preds.t0$Cohort==1,]$prob.prog)



##########################################################################
# ---------- Partial explained variance of individual predictors ---------
##########################################################################

# (1) With partR2 function
#remotes::install_github("mastoffel/partR2") 
library(partR2)
library(future)
library(furrr)

pred.vars <- c("BL_age","BL_sex","BL_bmi", "BL_smoking", "BL_map",
               "BL_hba1c", "BL_serumchol", "BL_hemo", "BL_uacr_log2",
               "BL_med_dm", "BL_med_bp", "BL_med_lipid")

# plan(multisession, workers=detectCores()*.60)
# partial.R2.cond <- partR2(fit.final, partvars = pred.vars,
#        R2_type = "conditional", nboot=10, parallel=T, 
#        max_level = 1,allow_neg_r2 = T)
# partial.R2.marg <- partR2(fit.final, partvars = pred.vars,
#                           R2_type = "marginal", nboot=10, parallel=T, 
#                           max_level = 1,allow_neg_r2 = T)
# plan(sequential)
# 
# partial.R2.cond <- partial.R2.cond$R2 %>% 
#   data.frame() %>% 
#   mutate(drop.R2=first(estimate)-estimate) %>%
#   mutate(scaled.drop = (drop.R2/sum(drop.R2))*first(estimate))
# 
# partial.R2.marg <- partial.R2.marg$R2 %>% 
#   data.frame() %>% 
#   mutate(drop.R2=first(estimate)-estimate) %>%
#   mutate(scaled.drop = (drop.R2/sum(drop.R2))*first(estimate))
# 
# partial.R2 <- list("Marginal"=partial.R2.marg, "Conditional"=partial.R2.cond)
# 
# write.xlsx(partial.R2, paste0(out.path, "tbl_partialR2.xlsx"), overwrite = TRUE)


# (1) Manual
r2.final <- r.squaredGLMM(fit.final)
out <- list()
for(i in 1:length(pred.vars)){
  model.vars <- pred.vars[-i]
  model.formula <- as.formula(paste0("FU_eGFR_epi ~ (Time + Time  *(",paste(model.vars,collapse="+") ,") + (1|Country) + (1+Time|PatID))"))
  fit.tmp <- lmer(model.formula,
                  data=data.full, REML=F, control=lmerControl(optimizer="bobyqa"))
  r2 <- r.squaredGLMM(fit.tmp)
  drop.r2 <- r2.final - r2
  
  out[[i]] <- c("pred.vars"=pred.vars[i], "R2m"=r2[1], "R2c"=r2[2], "drop.R2m"=drop.r2[1], "drop.R2c"=drop.r2[2])
}
res <- do.call(rbind, out)
tmp <- res %>%
  data.frame() %>%
  mutate_at(2:5, as.numeric) %>%
  mutate(scaled.dropR2m = (drop.R2m/sum(drop.R2m))*first(R2m),
         scaled.dropR2c = (drop.R2c/sum(drop.R2c))*first(R2c)) %>%
  mutate_at(2:7, round, 4)
tmp  

