#===============================================================================#
# Author: Mariella Gregorich
# Date: 30/06/2021
# Info: External validation with DIACORE
#===============================================================================#



# ======================== Data preparation ====================================

rm(list=ls())
source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)


# =================== External predictions ======================================
risk_model <- readRDS(paste0(out.path,"predmodel_lmerObject.rds"))

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
data.diacore.new <- full_join(data.diacore, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))




# ========================= Calibration plot ===================================
for(t in 2:7){
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==t,]$prior.pred,
                        cohort="extval", time=t, save=T, out.path = out.path, type="preUp")
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==t,]$pred,
                        cohort="extval", time=t, save=T, out.path = out.path, type="postUp")
  }

ggplot(data.diacore.new[data.diacore.new$Time_cat <=5,], aes(x=pred, y=FU_eGFR_epi, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Predicted eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = expression(paste("Observed eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=22)) 
ggsave(paste0(out.path, "fig_calibration_all_diacore_post.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')

df <- melt(data.diacore.new[data.diacore.new$Time_cat <=5 ,c("Time_cont","FU_eGFR_epi", "pred","prior.pred")], id.vars = c("Time_cont", "FU_eGFR_epi"))
levels(df$variable) <- c("Post-update", "Pre-update")
df$variable <- relevel(df$variable, "Pre-update")
df <- df[df$Time_cont>0,]

ggplot(df, aes(x=value, y=FU_eGFR_epi, fill=Time_cont, col=Time_cont)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Predicted eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) + 
  scale_y_continuous(expand = c(0, 0), name = expression(paste("Observed eGFR  [ml/min/1.73",m^2,"]")), limits = c(0,150)) +
  scale_color_continuous("Time") +
  scale_fill_continuous("Time") +
  geom_abline(intercept = 0) + 
  
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=16), strip.text = element_text(size = 18)) +
  facet_wrap(~variable)
ggsave(paste0(out.path, "fig_calibration_all_diacore_prepost.tiff"),  width=12, height=6, device='tiff', dpi=350, compression = 'lzw')



# ====================== Table: Performance measures ===========================
b=nboot
plan(multisession, gc=T, workers=detectCores()*.6)
res_ext_boot <- future_lapply(1:b, function(x){
  data.boot <- data.frame(draw_bootstrap_sample(data.diacore, cv=F))
  data.boot.t0 <- data.boot[data.boot$Time_cat==0,]
  res <- update_PredByBase(lmerObject=risk_model,
                           newdata=data.boot.t0,
                           cutpoint = slope_cutpoint,
                           timeVar = "Time", idVar="PatID", idVar2="Country",
                           times =seq(0,8,1), 
                           all_times=F)
  data.boot$Time <- data.boot$Time_cat
  data.boot.new <- full_join(data.boot, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))
  
  data.boot.list <- split(data.boot.new, as.factor(data.boot.new$Time)) 
  res <- lapply(data.boot.list, function(x) eval_preds(pred=x$pred, obs=x$FU_eGFR_epi, lmerObject=risk_model))
  df.res <- data.frame("Boot"=x,"Time"=as.numeric(names(res)), do.call(rbind, res))
  
  return(df.res)
}
, future.seed = T)
plan(sequential)

df.stats <- data.frame(do.call(rbind, res_ext_boot))



# Summarize and prepare output
tbl_performance <- list()
tbl_performance[[1]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], mean, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[2]]<- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], quantile, 0.05, na.rm=T) %>%
  mutate(Time=as.character(Time),
         Time=replace(Time, Time == "100", "Overall"))

tbl_performance[[3]] <- df.stats %>%
  group_by(Time) %>%
  summarize_at(.vars=names(.)[3:11], quantile, 0.95, na.rm=T) %>%
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
write.xlsx(tbl_performance, paste0(out.path, "tbl_perform_extval_full.xlsx"), overwrite = T)


# tbl_tmp <- lapply(sort(unique(data.diacore.new$Time)), function(x) eval_preds(pred=data.diacore.new[data.diacore.new$Time==x,]$pred, 
#                                                                               obs=data.diacore.new[data.diacore.new$Time==x,]$FU_eGFR_epi, 
#                                                                               lmerObject = risk_model))
# time.perf <- data.frame("Time"=sort(unique(data.diacore.new$Time)),round(do.call(rbind, tbl_tmp),3))
# ov.perf <- eval_preds(data.diacore.new[!data.diacore.new$Time==0,]$pred, data.diacore.new[!data.diacore.new$Time==0,]$FU_eGFR_epi, lmerObject = risk_model)
# tbl_performance <- list("overall"=ov.perf, "time-specific"=time.perf)
# write.xlsx(tbl_performance, paste0(out.path, "tbl_perform_extval.xlsx"), overwrite = TRUE)



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

tbl_ext <- tbl_performance$avg %>%
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
tbl_new[is.nan(tbl_new$est),] <-NA
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

