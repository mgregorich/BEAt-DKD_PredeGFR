#===============================================================================#
# Author: Mariella Gregorich
# Date: 24.05.2022
# Info: Sensitivity analysis
#===============================================================================#



# ============ Data preparation, IDA and general setup ==========================

rm(list = ls())
source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)

# ============== Sensitivity analysis - Model excluding sex & age ==============

data.full$Time <- data.full$Time_cont
fit.sens <- lmer(FU_eGFR_epi ~ (Time + Time  *(BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +                                                  
                                                 BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                 data=data.full, REML=F, control=lmerControl(optimizer="bobyqa"))


# UPDATE prediction with BL value
data.full.t0 <- data.full[data.full$Time==0,]
data.full.t0$Country <- "Unknown"
res <- update_PredByBase(lmerObject=fit.sens, newdata = data.full.t0, timeVar = "Time", idVar="PatID", idVar2="Country",
                         times = unique(data.full$Time_cat)[-1], all_times=F)
data.full$Time <- round(data.full$Time,0)
data.sens <- full_join(data.full, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))

# ================== Internal validation  =====================================

# ----  Save final model and extract components for shiny ----
saveRDS(fit.sens, paste0(out.path,"predmodel_lmerObject_sa.rds"))

# ----  Residual check -----
check_residuals(fit.sens, filename = "fig_residual_analysis_sa")


# ====================== Internal-external validation =========================

set.seed(123)
j=1; b=nboot

# Cross-validated predictions
data.full$fold <- as.numeric(data.full$Country)
res <- intext_crossvalidate(data.full, NA, return_preds = T, sa=T)
df.preds <- res$pred

# Cross-validated performance measures
plan(multisession, gc=T, workers=detectCores()*.6)
res_cv_boot <- future_lapply(1:b, function(x)  intext_crossvalidate(draw_bootstrap_sample(data.full),x, return_preds = F, sa=T), future.seed = T)
plan(sequential)
df.stats <- data.frame(do.call(rbind, lapply(res_cv_boot, `[[`, 1)))



# =================== External validation ======================================

rm(list=ls())
source("scr/setup.R")
source("scr/01_dataprep.R", print.eval=F)
risk_model <- readRDS(paste0(out.path,"predmodel_lmerObject_sa.rds"))

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


# -------------------- Calibration plot ---------------------------------------

for(t in 2:7){
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==t,]$prior.pred,
                        cohort="extval", time=t, save=T, out.path = out.path, type="preUp_sa")
  plot_calibration_cont(yobs=data.diacore.new[data.diacore.new$Time==t,]$FU_eGFR_epi, yhat=data.diacore.new[data.diacore.new$Time==t,]$pred,
                        cohort="extval", time=t, save=T, out.path = out.path, type="postUp_sa")
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
ggsave(paste0(out.path, "fig_calibration_all_diacore_post_sa.tiff"),  width=8, height=6, device='tiff', dpi=350, compression = 'lzw')


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
ggsave(paste0(out.path, "fig_calibration_all_diacore_prepost_sa.tiff"),  width=12, height=6, device='tiff', dpi=350, compression = 'lzw')



# ------------------ Table: Performance measures -------------------------------

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
  data.boot.new <- full_join(data.boot, res[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], 
                             by=c("PatID", "Time"))
  
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
write.xlsx(tbl_performance, paste0(out.path, "tbl_perform_extval_full_sa.xlsx"), overwrite = T)
