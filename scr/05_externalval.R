#===============================================================================#
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: External validation with DIACORE
#===============================================================================#



# ======================== Data preparation ====================================

rm(list=ls())
source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)


# =================== External model validation ======================================
risk_model <- readRDS(file.path(out.path, "riskpred_model.rds"))

data.diacore$Time <- data.diacore$Time_cont
data.diacore.t0 <- data.frame(data.diacore[data.diacore$Time_cat==0,])
data.diacore.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=risk_model, 
                            newdata = data.diacore.t0, 
                            cutpoint = slope_cutpoint,
                            timeVar = "Time", idVar="PatID", idVar2="Country",
                            times =sort(unique(data.diacore$Time_cat)), 
                            all_times=F)

# Summarize and prepare output
data.diacore$Time <- round(data.diacore$Time,0)
data.diacore.new <- full_join(data.diacore, res$Pred[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))



# ========================= Calibration plot ===================================
for(t in 2:7){
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==t,]$prior.pred,
                        cohort="extval", time=t, save=T, out.path = out.path, type="preUp")
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==t,]$pred,
                        cohort="extval", time=t, save=T, out.path = out.path, type="postUp")}

ggplot(data.diacore.new[data.diacore.new$Time_cat <=5,], aes(x=pred, y=FU_eGFR_epi, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = "Predicted eGFR", limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = "Observed EGFR", limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=22)) 
ggsave(paste0(out.path, "fig_calibration_all_diacore_post.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')

df <- melt(data.diacore.new[data.diacore.new$Time_cat <=5,c("Time_cont","FU_eGFR_epi", "pred","prior.pred")], id.vars = c("Time_cont", "FU_eGFR_epi"))
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



# ====================== Table: Performance measures ===========================
tbl_tmp <- lapply(sort(unique(data.diacore.new$Time)), function(x) eval_preds(pred=data.diacore.new[data.diacore.new$Time==x,]$pred, 
                                                                              obs=data.diacore.new[data.diacore.new$Time==x,]$FU_eGFR_epi, 
                                                                              lmerObject = risk_model))
time.perf <- data.frame("Time"=sort(unique(data.diacore.new$Time)),round(do.call(rbind, tbl_tmp),3))
ov.perf <- eval_preds(data.diacore.new[!data.diacore.new$Time==0,]$pred, data.diacore.new[!data.diacore.new$Time==0,]$FU_eGFR_epi, lmerObject = risk_model)
tbl_performance <- list("overall"=ov.perf, "time-specific"=time.perf)
write.xlsx(tbl_performance, paste0(out.path, "tbl_perform_extval.xlsx"), overwrite = TRUE)



# ========================== Probability of progression ========================
data.diacore.new.t0 <- data.diacore.new[data.diacore.new$Time==0,]

ggplot(data.diacore.new.t0, aes(x=pred.prob)) +
  geom_density(alpha=0.4, fill="grey") +
  scale_x_continuous(expand=c(0,0),"Probability of progression", limits = c(0,1)) +
  scale_y_continuous(expand=c(0,0),"Density") +
  theme_bw() +
  theme(text = element_text(size=18), legend.position = "bottom", legend.title = element_blank()) 
ggsave(paste0(out.path, "fig_prob_progression_SC",abs(slope_cutpoint),"_extval.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')



# ========================== Subject-specific trajectories =====================
set.seed(666)
df.melt <- melt(data.diacore.new[,c("PatID","Time", "FU_eGFR_epi","pred", "pred.lo", "pred.up")], id.vars = c("PatID","Time", "pred.lo", "pred.up"))
df.melt.small <- df.melt[df.melt$PatID %in% c(sample(unique(df.melt$PatID),16)),]

ggplot(data =df.melt.small, aes(x = Time, y = value,  col=variable)) +
  geom_ribbon(aes(x=Time, ymax=pred.up, ymin=pred.lo), fill="blue", alpha=.05) +
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



# ====================  Visualization of time-specific performance measures  ====
tbl_tmp.est <- read.xlsx(paste0(out.path, "tbl_perform_val_full.xlsx"), sheet=1) %>%
  dplyr::select(Time, R2, C, CalbSlope) %>%
  filter(!Time %in% "Overall") %>%
  mutate(type="est") %>%
  rename(CS=CalbSlope)

tbl_tmp.lo <- read.xlsx(paste0(out.path, "tbl_perform_val_full.xlsx"), sheet=2) %>%
  dplyr::select(Time, R2, C, CalbSlope) %>%
  filter(!Time %in% "Overall") %>%
  mutate(type="lo") %>%
  rename(CS=CalbSlope)

tbl_tmp.up <- read.xlsx(paste0(out.path, "tbl_perform_val_full.xlsx"), sheet=3)%>%
  dplyr::select(Time, R2, C, CalbSlope) %>%
  filter(!Time %in% "Overall") %>%
  mutate(type="up") %>%
  rename(CS=CalbSlope)


tbl_int <- rbind(tbl_tmp.est, tbl_tmp.lo, tbl_tmp.up) %>%
  mutate(Cohort="dev")

tbl_ext <- tbl_performance$`time-specific` %>%
  dplyr::select(Time, R2, C, CalbSlope) %>%
  rename(CS=CalbSlope) %>%
  mutate(type="est",
         Cohort="ext")

tbl_perf <- rbind(tbl_int, tbl_ext)
tbl_perf_melt <- melt(tbl_perf, id.vars = c("Time", "type", "Cohort"))
tbl_perf_melt$tmp <- paste0(tbl_perf_melt$variable,"+",tbl_perf_melt$Cohort)
tbl_new <- pivot_wider(tbl_perf_melt, names_from = type, values_from = value)
tbl_new <- tbl_new[tbl_new$Time %in% 0:5,]
tbl_new$variable <- fct_recode(tbl_new$variable, "Calibration slope"="CS", "C statistic"="C", "R2"="R2")
tbl_new$Cohort <- fct_recode(tbl_new$Cohort, "Development cohort"="dev", "Validation cohort"="ext")

ggplot(tbl_new, aes(x=Time, y=est, group=tmp, col=tmp, fill=variable)) +
  geom_point(size=2) +
  geom_line(aes(linetype = Cohort), size=0.5) +
  scale_y_continuous("Estimate") +
  geom_ribbon(aes(ymin = lo,
                  ymax = up, fill=variable), alpha=0.1) +
  scale_linetype_manual("",values=c("solid", "longdash"))+
  scale_fill_manual(values=c("red2","dodgerblue2", "mediumseagreen"), guide="none") +
  scale_color_manual(values=c("dodgerblue2","dodgerblue4", "mediumseagreen","seagreen", "red2","red4"),guide="none") +
  theme_bw() +
  theme(text = element_text(size=16), legend.position = "top")
ggsave(paste0(out.path, "fig_performance_dev_ext.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')

