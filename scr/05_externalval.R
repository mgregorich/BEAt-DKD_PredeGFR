###################################################
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: External validation with DIACORE
###################################################


rm(list=ls())

source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)

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
data.diacore.new <- full_join(data.diacore, res$Pred[,c("PatID", "Time","prior.pred","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")], by=c("PatID", "Time"))



# --- Calibration plot
for(t in 2:7){
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==t,]$prior.pred, fit = risk_model,
                        cohort="extval", time=t, save=T, out.path = out.path, type="preUp")
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==t,]$pred, fit = risk_model,
                        cohort="extval", time=t, save=T, out.path = out.path, type="postUp")}

ggplot(data.diacore.new, aes(x=pred, y=FU_eGFR_epi, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = "Predicted eGFR", limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = "Observed EGFR", limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=22)) 
ggsave(paste0(out.path, "fig_calibration_all_diacore_post.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')

df <- melt(data.diacore.new[,c("Time_cont","FU_eGFR_epi", "pred","prior.pred")], id.vars = c("Time_cont", "FU_eGFR_epi"))
levels(df$variable) <- c("Post-update", "Pre-update")
df$variable <- relevel(df$variable, "Pre-update")

ggplot(df, aes(x=value, y=FU_eGFR_epi, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = "Predicted eGFR", limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = "Observed EGFR", limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=16), strip.text = element_text(size = 18)) +
  facet_wrap(~variable)
ggsave(paste0(out.path, "fig_calibration_all_diacore_prepost.tiff"),  width=12, height=6, device='tiff', dpi=350, compression = 'lzw')



# --- Model performance and validation
tbl_tmp <- lapply(sort(unique(data.diacore.new$Time)), function(x) eval_preds(pred=data.diacore.new[data.diacore.new$Time==x,]$pred, 
                                                                              obs=data.diacore.new[data.diacore.new$Time==x,]$FU_eGFR_epi, 
                                                                              lmerObject = risk_model))
time.perf <- data.frame("Time"=sort(unique(data.diacore.new$Time)),round(do.call(rbind, tbl_tmp),3))
ov.perf <- eval_preds(data.diacore.new[!data.diacore.new$Time==0,]$pred, data.diacore.new[!data.diacore.new$Time==0,]$FU_eGFR_epi, lmerObject = risk_model)
tbl_performance <- list("overall"=ov.perf, "time-specific"=time.perf)
write.xlsx(tbl_performance, paste0(out.path, "tbl_perform_extval.xlsx"), overwrite = TRUE)



# --- Probability of progression
data.diacore.new.t0 <- data.diacore.new[data.diacore.new$Time==0,]

ggplot(data.diacore.new.t0, aes(x=prob.prog)) +
  geom_density(alpha=0.4, fill="grey") +
  scale_x_continuous(expand=c(0,0),"Probability of progression", limits = c(0,1)) +
  scale_y_continuous(expand=c(0,0),"Density") +
  theme_bw() +
  theme(text = element_text(size=18), legend.position = "bottom", legend.title = element_blank()) 
ggsave(paste0(out.path, "fig_prob_progression_SC",abs(slope_cutpoint),"_extval.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')



# --- Subject-specific trajectories
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
ggsave(paste0(out.path, "fig_indvPred_eGFR_diacore_extval.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')

