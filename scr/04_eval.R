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
                                                                 N=length(unique(fit.final$data$PatID)), 
                                                                 k=sum(anova(fit.final)$numDF)))
tbl_performance <- data.frame("Time"=seq(0,8,1),round(do.call(rbind, tbl_tmp),3))
# tbl_performance$C <- round(sapply(0:7, function(x) mean(df.stats[df.stats$Time==x,]$C, na.rm=T)),3)
write.xlsx(tbl_performance, paste0(out.path, "tbl_perform_val.xlsx"), 
           overwrite = TRUE)


################################################################################
# ----------------------------- Calibration ------------------------------------
################################################################################
# Before update
# plot_calibration_cont(yobs=df.preds[df.preds$Time==1,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==1,]$prior.pred, time=1, save=T)
# plot_calibration_cont(yobs=df.preds[df.preds$Time==2,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==2,]$prior.pred, time=2, save=T)
# plot_calibration_cont(yobs=df.preds[df.preds$Time==3,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==3,]$prior.pred, time=3, save=T)
# plot_calibration_cont(yobs=df.preds[df.preds$Time==4,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==4,]$prior.pred, time=4, save=T)
# plot_calibration_cont(yobs=df.preds[df.preds$Time==5,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==5,]$prior.pred, time=5, save=T)
# plot_calibration_cont(yobs=df.preds[df.preds$Time==6,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==6,]$prior.pred, time=6, save=T)
# plot_calibration_cont(yobs=df.preds[df.preds$Time==7,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==7,]$prior.pred, time=7, save=T)

# After update
plot_calibration_cont(yobs=df.preds[df.preds$Time==1,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==1,]$pred, time=1, save=T, out.path = out.path)
plot_calibration_cont(yobs=df.preds[df.preds$Time==2,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==2,]$pred, time=2, save=T, out.path = out.path)
plot_calibration_cont(yobs=df.preds[df.preds$Time==3,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==3,]$pred, time=3, save=T, out.path = out.path)
plot_calibration_cont(yobs=df.preds[df.preds$Time==4,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==4,]$pred, time=4, save=T, out.path = out.path)
plot_calibration_cont(yobs=df.preds[df.preds$Time==5,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==5,]$pred, time=5, save=T, out.path = out.path)
plot_calibration_cont(yobs=df.preds[df.preds$Time==6,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==6,]$pred, time=6, save=T, out.path = out.path)
plot_calibration_cont(yobs=df.preds[df.preds$Time==7,]$FU_eGFR_epi, yhat=df.preds[df.preds$Time==7,]$pred, time=7, save=T, out.path = out.path)

# Probability of progression
pred_prob = df.preds[df.preds$Time==0,]$prob.prog
true_prob = df.preds[df.preds$Time==0,]$true.prob
plot_calibration_bin(pred=pred_prob, true=true_prob, out.path = out.path, save=T)
Brier <- mean((pred_prob-true_prob)^2)
C.index <- somers2(pred_prob, true_prob)["C"]
Deviance<- -2*sum(true_prob*log(pred_prob) + (1-true_prob)*log(1-pred_prob) )
print(paste0("Brier = ", round(Brier,2), ", C statistic = ", round(C.index,2), " and deviance = ", round(Deviance,2)))



################################################################################
# ---------------------------- Model fit ---------------------------------------
################################################################################

# ----- Plot subject-specific trajectories
which.max(abs(fit.final$coefficients$random$PatID[,2]))

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
df.preds$Cohort <- as.factor(df.preds$Cohort)
ggplot(df.preds,aes(x = Time, y = FU_eGFR_epi, fill=Cohort))  +
  geom_boxplot(aes(x=Time, fill=Cohort), outlier.shape = NA) +
  stat_boxplot(aes(fill=Cohort),geom ='errorbar') +
  scale_y_continuous("Observed eGFR") +
  scale_fill_manual(values=c("orange1", "royalblue1"), labels=c("GCKD", "PROVALID")) +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(paste0(out.path, "fig_longiobs_eGFR_cohort_dev.png"),width=10, height=6)



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
df <- data.frame(variable=rownames(intervals(fit.final)[[1]]),
                 effect=as.numeric(intervals(fit.final)[[1]][,2]),
                 lower=as.numeric(intervals(fit.final)[[1]][,1]),
                 upper=as.numeric(intervals(fit.final)[[1]][,3]))
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

