#===============================================================================#
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: Model evaluation of the cross-validated and bias-corrected model results
#===============================================================================#


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
write.xlsx(tbl_performance, here::here(out.path, "model_main","tbl_perform_val_full.xlsx"), overwrite = T)



# ===================== Calibration (continuous) ==============================

for(t in 2:7){
  # Before update
  plot_calibration_cont(yobs=df.preds[df.preds$Time==t,]$FU_eGFR_epi_2021, yhat=df.preds[df.preds$Time==t,]$prior.pred, 
                        time=t, save=T, folder="model_main", out.path = out.path, cohort="dev", type="preUp")
  # After update
  plot_calibration_cont(yobs=df.preds[df.preds$Time==t,]$FU_eGFR_epi_2021, yhat=df.preds[df.preds$Time==t,]$pred, 
                        time=t, save=T, folder="model_main", out.path = out.path, cohort="dev", type="postUp")
}

ggplot(df.preds, aes(x=pred, y=FU_eGFR_epi_2021, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = "Predicted eGFR", limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = "Observed EGFR", limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=22)) 
ggsave(here::here(out.path, "model_main", "fig_calibration_all_dev_post.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')


# ==================== Calibration: Probability of progression =================
pred_prob = df.preds[df.preds$Time==0,]$pred.prob
true_prob = df.preds[df.preds$Time==0,]$true.prob
plot_calibration_bin(pred=pred_prob, true=as.factor(true_prob), out.path = out.path, save=F)
Brier <- mean((pred_prob-true_prob)^2)
C.index <- somers2(pred_prob, true_prob)["C"]
Deviance<- -2*sum(true_prob*log(pred_prob) + (1-true_prob)*log(1-pred_prob) )
CS <- lm(true_prob ~ pred_prob)$coefficients[2]
cat(paste0("RISK ASSESSMENT: \nBrier = ", round(Brier,2), ", C statistic = ", round(C.index,2), ", calibration slope = ",round(CS,2)," and deviance = ", round(Deviance,2)),"\n")


par(mfrow=c(2,2))
tiff(here::here(out.path, "model_main", "hist_slopes.png"),
     res = 300,                                              
     width = 5, height = 4, units = 'in',
     compression = c( "lzw"))
hist(df.preds$pred.slope, main = "Predicted eGFR slope")
hist(df.preds$true.slope, main = "Observed eGFR slope (LM)")
dev.off()

ggplot(df.preds, aes(x=pred.slope, y=true.slope)) +
  geom_point() +
  geom_abline(intercept = 0, col="red") + 
  scale_x_continuous("Predicted eGFR slope ") +
  scale_y_continuous("Observed eGFR slope (LM)") +
  theme_bw() +
  theme(text=element_text(size=16))
ggsave(here::here(out.path, "model_main", "fig_prob_progression_full.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')


# Pick individual with largest deviation of treu and predicted
df.preds[which.max(df.preds$true.slope - df.preds$pred.slope),]$PatID


# ==========================  Estimated eGFR slopes  ===============

# ----- Plot individual-specific trajectories
which.max(abs(coef(fit.final)$PatID[,"Time"]))

set.seed(666)
df.melt <- melt(df.preds[,c("PatID","Time", "FU_eGFR_epi_2021","pred", "pred.lo", "pred.up", "true.slope", "pred.slope")], id.vars = c("PatID","Time", "pred.lo", "pred.up", "true.slope", "pred.slope"))
df.melt.small <- df.melt[df.melt$PatID %in% c(sample(unique(df.melt$PatID),16)),]
df.melt.small$Time <- as.numeric(df.melt.small$Time)
graphLabels <- df.melt.small[!duplicated(df.melt.small$PatID),]
graphLabels$label1 = paste0("Obs.(LM): ", round(graphLabels$true.slope,2), "\nPred: ", round(graphLabels$pred.slope,2))

ggplot(data =df.melt.small, aes(x = Time, y = value,  col=variable)) +
  geom_ribbon(aes(x=Time, ymax=pred.up, ymin=pred.lo), fill="blue", alpha=.05) +
  geom_point(size=1) +
  geom_line(data=df.melt.small[!is.na(df.melt.small$value),], aes(color=variable)) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed")) +
  scale_color_manual(values=c("black", "blue", "red", "red"), labels=c("Observed", "Predicted", "95% CI", "95% CI")) +
  scale_x_continuous(expand=c(0,0),breaks = seq(0,8,1), limits = c(0,8)) +
  scale_y_continuous(expand=c(0,0), "eGFR", limits=c(0,150)) +
  ggtitle("") +
  theme_bw() +  
  facet_wrap(~PatID) +
  theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=16)) +
  geom_text(data=graphLabels, aes(x = 4, y = 125, label =label1),  size=3) 
ggsave(here::here(out.path, "model_main", "fig_indvPred_eGFR_dev.png"),width=10, height=15)

df.melt.small <- df.melt[df.melt$PatID %in% "PV14750",]
df.melt.small$Time <- as.numeric(df.melt.small$Time)
ggplot(data =df.melt.small, aes(x = Time, y = value,  col=variable)) +
  geom_ribbon(aes(x=Time, ymax=pred.up, ymin=pred.lo), fill="blue", alpha=.05) +
  geom_point(size=1) +
  geom_line(data=df.melt.small[!is.na(df.melt.small$value),]) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed")) +
  scale_color_manual(values=c("black", "blue", "red", "red"), labels=c("Observed", "Predicted", "95% CI", "95% CI")) +
  scale_x_continuous(expand=c(0,0),breaks = seq(0,8,1), limits = c(0,8)) +
  scale_y_continuous(expand=c(0,0), "eGFR") +
  ggtitle("") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=16))



# ---- Individual specific trajectories, comparing updated and no-update predictions
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
ggsave(here::here(out.path, "model_main", "fig_indvPred_eGFR_dev_noUpdate.png"),width=10, height=6)


# ---- Mean country trajectories
df.preds %>%
  group_by(Country, Time) %>%
  summarise(pm = median(FU_eGFR_epi_2021, na.rm=T)) %>%
  ggplot(aes(x = Time, y = pm, group = Country, col=Country)) +
  geom_point() +
  geom_line(aes(color=Country)) +
  theme_bw()


# ---- Longitudinal trajectory of observed eGFR per country
df.preds$Time <- as.factor(df.preds$Time)
ggplot(df.preds,aes(x = Time, y = FU_eGFR_epi_2021, fill=Country))  +
  geom_boxplot(aes(x=Time, fill=Country), outlier.shape = NA) +
  stat_boxplot(aes(fill=Country),geom ='errorbar') +
  scale_y_continuous("Observed eGFR") +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(here::here(out.path, "model_main", "fig_longiObs_eGFR.png"),width=12, height=6)


# Longitudinal trajectory of predicted eGFR per cohort
data.diacore$Cohort <- as.factor(data.diacore$Cohort)
df.tmp <- data.frame(rbind(df.preds[,c("Time_cat", "FU_eGFR_epi_2021", "Cohort")], data.diacore[,c("Time_cat", "FU_eGFR_epi_2021", "Cohort")]))
df.tmp$Time_cat <- as.factor(df.tmp$Time_cat)
df.tmp <- df.tmp[!is.na(df.tmp$Time_cat),]
ggplot(df.tmp,aes(x = Time_cat, y = FU_eGFR_epi_2021, fill=Cohort))  +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom ='errorbar') +
  scale_y_continuous(expand=c(0,0),"Observed eGFR", lim=c(0,150), breaks=seq(0,150,25)) +
  scale_x_discrete("Time") +
  scale_fill_brewer(palette = "Set2", labels=c("GCKD","PROVALID","DIACORE")) +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(here::here(out.path, "model_main", "fig_longiobs_eGFR_cohort.png"),width=10, height=6)



# ---- Longitudinal trajectory of predicted eGFR per country
ggplot(df.preds,aes(x = Time, y = pred, fill=Country))  +
  geom_boxplot(aes(x=Time), outlier.shape = NA) +
  stat_boxplot(geom ='errorbar') +
  scale_y_continuous("Predicted eGFR") +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(here::here(out.path, "model_main", "fig_longipred_eGFR.png"),width=12, height=6)

# ---- Longitudinal trajectory of predicted eGFR per cohort
df.preds$Cohort <- as.factor(df.preds$Cohort)
ggplot(df.preds,aes(x = Time, y = pred, fill=Cohort))  +
  geom_boxplot(aes(x=Time), outlier.shape = NA) +
  stat_boxplot(geom ='errorbar') +
  scale_y_continuous("Predicted eGFR") +
  scale_fill_manual(values=c("orange1", "royalblue1"), labels=c("GCKD", "PROVALID")) +
  theme_bw() +
  theme(text = element_text(size=16))
ggsave(here::here(out.path, "model_main", "fig_longipred_eGFR_cohort_dev.png"),width=10, height=6)


# ================= Table: Model coefficients ================================
df <- data.frame(variable=rownames(summary(fit.final)$coefficients),
                 effect=summary(fit.final)$coefficients[,1],
                 lower=summary(fit.final)$coefficients[,1] - 1.96*summary(fit.final)$coefficients[,2],
                 upper=summary(fit.final)$coefficients[,1] + 1.96* summary(fit.final)$coefficients[,2]) 
df$variable <- factor(df$variable, levels = df$variable[length(df$variable):1])
df$Group <- ifelse(str_detect(df$variable, "Time:"), "Main effect", "Interaction effect")
df <- mutate(df, Group = fct_rev(Group)) 
df$variable <- as.factor(str_replace(df$variable, "Time:", ""))
df$variable <- factor(df$variable, levels=c("(Intercept)","Time","BL_age","BL_sex1","BL_bmi","BL_smoking1","BL_map","BL_hba1c","BL_serumchol", "BL_hemo",
                                            "BL_uacr_log2", "BL_med_dm1","BL_med_bp1","BL_med_lipid1"))
levels(df$variable) <- list("Intercept"="(Intercept)","Time"="Time", Age  = "BL_age", "Sex (female)" = "BL_sex1", BMI="BL_bmi", "Smoking (ever)"="BL_smoking1", MAP="BL_map", HbA1c="BL_hba1c",
                            "Serum Chol."="BL_serumchol","Hemoglobin"="BL_hemo", "log2 UACR"="BL_uacr_log2", "Glucose-low. Med."="BL_med_dm1", "Blood pressure-low. Med."="BL_med_bp1", "Lipid-low. Med."="BL_med_lipid1")


df[,2:4] <- apply(df[,2:4],2, round, digits=3)
df$name <- rownames(df)
df$Variable <- c(rep("Constant",2), rep(c("Age", "Sex", "BMI", "Smoking", "MAP", "HbA1c","Serum chol.",
                                          "Hemoglobin", "log2UACR", "GL Med.", "BPL Med.","LL Med."),2))
df$CI <- paste0("(",df$lower,", " ,df$upper,")")
tbl_fixeff <- cbind(df[!str_detect(df$name, "Time"),c(7,2,8)], df[str_detect(df$name, "Time"),c(2,8)])
colnames(tbl_fixeff) <- c("Variable", "Effect.Baseline", "95% CI", "Effect.Slope", "95% CI")
write.xlsx(tbl_fixeff, here::here(out.path, "model_main", "tbl_coeff_unstd.xlsx"), overwrite = TRUE)




# ===================== Forest plot of standardized coefficients ===============
data.scaled <- data.full %>%mutate_at(vars(Time, BL_age, BL_bmi, BL_map, BL_hba1c, 
                                           BL_serumchol, BL_hemo, BL_uacr_log2, FU_eGFR_epi_2021), scale)
fit.scaled <- lmer(FU_eGFR_epi_2021 ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +                                                  
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
df$variable <- factor(df$variable, levels=c("(Intercept)","Time","BL_age","BL_sex1","BL_bmi","BL_smoking1","BL_map","BL_hba1c","BL_serumchol", "BL_hemo",
                                            "BL_uacr_log2", "BL_med_dm1","BL_med_bp1","BL_med_lipid1"))
levels(df$variable) <- list("Intercept"="(Intercept)","Time"="Time", Age  = "BL_age", "Sex (female)" = "BL_sex1", BMI="BL_bmi", "Smoking (ever)"="BL_smoking1", MAP="BL_map", HbA1c="BL_hba1c",
                            "Serum Chol."="BL_serumchol","Hemoglobin"="BL_hemo", "log2 UACR"="BL_uacr_log2", "Glucose-low. Med."="BL_med_dm1", "Blood pressure-low. Med."="BL_med_bp1", "Lipid-low. Med."="BL_med_lipid1")

df[,2:4] <- apply(df[,2:4],2, round, digits=3)
df$name <- rownames(df)
df$Variable <- c(rep("Constant",2), rep(c("Age", "Sex", "BMI", "Smoking", "MAP", "HbA1c","Serum chol.",
                                          "Hemoglobin", "log2UACR", "GL Med.", "BPL Med.","LL Med."),2))
df$CI <- paste0("(",df$lower,", " ,df$upper,")")
tbl_fixeff <- cbind(df[!str_detect(df$name, "Time"),c(7,2,8)], df[str_detect(df$name, "Time"),c(2,8)])
colnames(tbl_fixeff) <- c("Variable", "Effect.Baseline", "95% CI", "Effect.Slope", "95% CI")
write.xlsx(tbl_fixeff, here::here(out.path, "model_main","tbl_coeff_std.xlsx"), overwrite = TRUE)


ggplot(data=df[-1,], aes(y=reorder(variable,desc(variable)), x=effect, xmin=lower, xmax=upper)) +
  geom_point(size=2, shape=1) + 
  geom_errorbarh(height=.25) +
  scale_y_discrete("", labels=c("HbA1c" = expression(HbA[1][c]), "log2 UACR"=expression(paste(log[2], ' UACR')))) +
  scale_x_continuous("Standardized Effect", limits = c(-0.4,0.4), breaks=seq(-0.4,0.4,0.1)) +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  facet_wrap(~Group) +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave(here::here(out.path, "model_main", "plot_forest_standardized.tiff"),  width=8, height=4, device='tiff', dpi=350, compression = 'lzw')
ggsave(here::here(out.path, "model_main", "plot_forest_standardized.wmf"),  width=8, height=4, dpi=350)


# ====================== Individual risk predictions ==========================
df.preds.t0 <- df.preds[df.preds$Time==0,]

ggplot(df.preds.t0, aes(x=pred.prob, fill=Cohort)) +
  geom_density(alpha=0.4) +
  scale_fill_manual(values=c("orange1", "royalblue1"), labels=c("GCKD", "PROVALID")) +
  scale_x_continuous("Probability of progression") +
  scale_y_continuous("Density") +
  theme_bw() +
  theme(text = element_text(size=18), legend.position = "bottom", legend.title = element_blank()) 
filename <- paste0("fig_prob_progression_",abs(slope_cutpoint),"_dev.png")
ggsave(here::here(out.path, "model_main", filename),width=8, height=6)

summary(df.preds.t0$pred.prob)
summary(df.preds.t0[df.preds.t0$Cohort==0,]$pred.prob)
summary(df.preds.t0[df.preds.t0$Cohort==1,]$pred.prob)



# ================= Relative partial importance of individual predictors ========

tbl_partialR2 <- data.frame(matrix(NA, nrow=13, ncol = 7))
colnames(tbl_partialR2) <- c("term", "R2m", "R2m.lo",  "R2m.up", "R2c",  "R2c.lo", "R2c.up")

# --- Full
i=1
plan(multisession, gc=T, workers=detectCores()*.6)
out <- t(future_sapply(1:nboot, function(x) {
  data.boot <- draw_bootstrap_sample(data.full)
  fit.full <- lmer(FU_eGFR_epi_2021 ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +
                                                   BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                   data=data.boot, REML=F, control=lmerControl(optimizer="bobyqa"))
  r2.full <- r.squaredGLMM(fit.full)
  return(r2.full)
  }, future.seed = T))
plan(sequential)
tbl_partialR2[i, ] <- c("Full", apply(out, 2, function(x) c(mean(x), quantile(x, 0.05), quantile(x, 0.95))))


# --- Partial
for(i in 1:length(pred.vars)){
  # Reduced
  model.vars <- pred.vars[-i]
  model.formula <- as.formula(paste0("FU_eGFR_epi_2021 ~ (Time + Time  *(",paste(model.vars,collapse="+") ,") + (1|Country) + (1+Time|PatID))"))
  plan(multisession, gc=T, workers=detectCores()*.6)
  out <- t(future_sapply(1:nboot, function(x) {
    data.boot <- draw_bootstrap_sample(data.full)
    fit.tmp <- lmer(model.formula,
                    data=data.boot, REML=F, control=lmerControl(optimizer="bobyqa"))
    r2.tmp <- r.squaredGLMM(fit.tmp)
    return(r2.tmp)
  }, future.seed = T))
  plan(sequential)

  tbl_partialR2[i+1, ] <- c(pred.vars[i], apply(out, 2, function(x) c(mean(x), quantile(x, 0.05), quantile(x, 0.95))))
}
write.xlsx(tbl_partialR2, here::here(out.path, "model_main", "tbl_partialR2.xlsx"), overwrite = T)


# --- Visualization
tbl_partR2 <- read.xlsx(here::here(out.path, "model_main", "tbl_partialR2.xlsx"))
tbl_partR2$term <- as.factor(tbl_partR2$term)
levels(tbl_partR2$term) <- list("Full Model"="Full",Age  = "BL_age", "Sex" = "BL_sex", BMI="BL_bmi", "Smoking"="BL_smoking", MAP="BL_map", HbA1c="BL_hba1c",
                                "Serum Chol."="BL_serumchol","Hemoglobin"="BL_hemo", "log2 UACR"="BL_uacr_log2", 
                                "Glucose-low. Med."="BL_med_dm", "Blood pressure-low. Med."="BL_med_bp", "Lipid-low. Med."="BL_med_lipid")
tmp1 <- tbl_partR2[,c("term","R2m", "R2m.lo", "R2m.up")]
tmp1$type <- "marginal"
tmp2 <- tbl_partR2[,c("term","R2c", "R2c.lo", "R2c.up")]
tmp2$type <- "conditional"
names(tmp1) <- names(tmp2) <- c("term", "est", "ci.lo", "ci.up", "type")

tbl_partR2 <- data.frame(rbind(tmp1, tmp2))
tbl_partR2[,2:4] <- apply(tbl_partR2[,2:4], 2, as.numeric)
dummy <- data.frame(type = c("marginal", "conditional"), 
                    val = c(tbl_partR2[tbl_partR2$term %in% "Full Model",]$est[1], 
                            tbl_partR2[tbl_partR2$term %in% "Full Model",]$est[2]),
                    cols=c("skyblue4", "skyblue2"))

ggplot(tbl_partR2, aes(x=est, y=term, col=type)) +
  geom_point() +
  geom_line() +
  geom_errorbarh(aes(xmin = ci.lo, xmax=  ci.up), height=0.25) +
  geom_vline(data=dummy, aes(xintercept=val, colour=cols), size=0.2, linetype="dashed") +
  scale_y_discrete("",limits=rev) +
  scale_x_continuous("Partial R2 and 95% CI") +
  scale_color_manual("", values=c("skyblue2", "skyblue4","skyblue2", "skyblue4")) +
  theme_bw() +
  theme(text=element_text(size=16), legend.position = "None",
        axis.text.y = element_text(face = c('plain', 'plain', 'plain', 'plain', 'plain', 'plain', 
                                            'plain', 'plain', 'plain', 'plain', 'plain', 'plain','bold'))) +
  facet_wrap(~type, scales="free_x")
ggsave(here::here(out.path, "model_main", "fig_partialR2_bold.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')

