###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: External validation with DIACORE
###################################################


# rm(list=ls())
# 
# source("scr/setup.R")
# source("scr/01_dataprep.R", print.eval=F)

# ------------ Model validation ---------------------------
risk_model <- readRDS(file.path(out.path, "riskpred_model.rds"))

data.diacore$Time <- data.diacore$Time_cont
data.diacore.t0 <- data.frame(data.diacore[data.diacore$Time_cat==0,])
data.diacore.t0$Country <- "Unknown"
res <- LongPred_ByBase_lmer(lmerObject=risk_model, 
                       newdata = data.diacore.t0, 
                       cutpoint = slope_cutpoint,
                       timeVar = "Time", idVar="PatID", idVar2="Country",
                       times =sort(unique(data.diacore$Time_cat)), 
                       all_times=F)

# Summarize and prepare output
data.diacore$Time <- round(data.diacore$Time,0)
data.diacore.new <- full_join(data.diacore, res$Pred[,c("PatID", "Time","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")], by=c("PatID", "Time"))


# Calibration plot

plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==2,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==2,]$pred, fit = risk_model,
                 cohort="extval", time=2, save=T, out.path = out.path)
plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==3,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==3,]$pred, fit = risk_model,
                      cohort="extval", time=3, save=T, out.path = out.path)
plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==4,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==4,]$pred, fit = risk_model, 
                      cohort="extval", time=4, save=T, out.path = out.path)
plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==5,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==5,]$pred, fit = risk_model, 
                      cohort="extval", time=5, save=T, out.path = out.path)
plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==6,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==6,]$pred, fit = risk_model, 
                      cohort="extval", time=6, save=T, out.path = out.path)
plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==7,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==7,]$pred, fit = risk_model, 
                      cohort="extval", time=7, save=T, out.path = out.path)


# Model performance and validation
tbl_tmp <- lapply(sort(unique(data.diacore.new$Time)), function(x) eval_preds(pred=data.diacore.new[data.diacore.new$Time==x,]$pred, 
                                                                              obs=data.diacore.new[data.diacore.new$Time==x,]$FU_eGFR_epi, 
                                                                              lmerObject = risk_model))
tbl_performance <- data.frame("Time"=sort(unique(data.diacore.new$Time)),round(do.call(rbind, tbl_tmp),3))
write.xlsx(tbl_performance, paste0(out.path, "tbl_perform_extval.xlsx"), 
           overwrite = TRUE)


# Probability of progression
data.diacore.new.t0 <- data.diacore.new[data.diacore.new$Time==0,]

ggplot(data.diacore.new.t0, aes(x=prob.prog)) +
  geom_density(alpha=0.4, fill="grey") +
  scale_x_continuous(expand=c(0,0),"Probability of progression", limits = c(0,1)) +
  scale_y_continuous(expand=c(0,0),"Density") +
  theme_bw() +
  theme(text = element_text(size=18), legend.position = "bottom", legend.title = element_blank()) 
ggsave(paste0(out.path, "fig_prob_progression_SC",abs(slope_cutpoint),"_extval.png"),
       width=8, height=6)


# Subject-specific trajectories
set.seed(666)
df.melt <- melt(data.diacore.new[,c("PatID","Time", "FU_eGFR_epi","pred", "pred.low", "pred.upp")], id.vars = c("PatID","Time", "pred.low", "pred.upp"))
df.melt.small <- df.melt[df.melt$PatID %in% c(sample(unique(df.melt$PatID),16)),]

ggplot(data =df.melt.small, aes(x = Time, y = value,  col=variable)) +
  geom_ribbon(aes(x=Time, ymax=pred.upp, ymin=pred.low), fill="blue", alpha=.05) +
  geom_point(size=1) +
  geom_line(data=df.melt.small[!is.na(df.melt.small$value),], aes(color=variable)) +
  scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed")) +
  scale_color_manual(values=c("black", "blue", "red", "red"), labels=c("Observed", "Predicted", "95% CI", "95% CI")) +
  scale_x_continuous(expand=c(0,0),breaks = seq(0,7,1), limits = c(0,7)) +
  scale_y_continuous(expand=c(0,0), "eGFR") +
  ggtitle("") +
  theme_bw() +
  facet_wrap(~PatID) +
  theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=16))
ggsave(paste0(out.path, "fig_indvPred_eGFR_diacore_extval.png"),width=10, height=6)

