#===============================================================================#
# Author: MG
# Date: 03.08.2021
# Info: Update longitudinal predictions by baseline measurements
#===============================================================================#



# --------------------------- GENERAL ---------------------

recodeDates <- function(x){
  x=str_replace(x, "Mrz", "März")
  x.new=if_else(str_detect(x, "^[:digit:]+$"), as.character(as.Date(as.numeric(x), origin = "1899-12-30")), as.character(as.Date(x, "%d%b%Y")))
  x.new = as.Date(x.new)
  return(x.new)
}

right_rows <- function (data, times, ids, Q_points) {
  fids <- factor(ids, levels = unique(ids))
  if (!is.list(Q_points)){Q_points <- split(Q_points, row(Q_points))}
  ind <- mapply(findInterval, Q_points, split(times, fids))
  ind[ind < 1] <- 1
  rownams_id <- split(row.names(data), fids)
  ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
  data[c(ind), ]
}


round_0 <- function(x, digits){
  return(sprintf(paste0('%.',digits,'f'),x))
}
to_numeric <- function(x){
  return(as.numeric(as.character(x)))
}

# ======================== DATA CLEANING ======================================

getProvalidFU <- function(sheet, path){
  df <- read_excel(path, sheet)
  df <- clean_names(df)
  rel.cols <- c("patient_reference_id", "date_yyyy_mm_dd_10","serum_creatinine_mg_dl", "serum_creatinine_mmol_l")
  df <- df[,rel.cols]
  df$source <- sheet
  df$time <- ifelse(str_detect(sheet, "FOLLOW-UP"), str_match(sheet, "FOLLOW-UP\\s*(.*?)\\s*LOCAL LABORATORY")[,2], "0")
  df$FU_serumcrea <- ifelse(is.na(df$serum_creatinine_mg_dl), df$serum_creatinine_mmol_l/100, df$serum_creatinine_mg_dl)
  df <- df[!is.na(df$FU_serumcrea),]
  colnames(df)[1] <- "Patient reference [ID]"
  
  return(df)
}


calcBMI <- function(height_m, weight_kg){
  BMI <-(weight_kg/height_m^2)
  return(BMI)
}

calcUACR <- function(df){
  # df= data.provalid
  
  df <- df %>%
    mutate(urin_alb1st_mg_dl = ifelse(is.na(x1st_urinary_albumin_mg_l), x1st_urinary_albumin_g_l*100, x1st_urinary_albumin_mg_l/10),
           urin_alb2nd_mg_dl = ifelse(is.na(x2nd_urinary_albumin_mg_l), x2nd_urinary_albumin_g_l*100, x2nd_urinary_albumin_mg_l/10),
           urin_alb3rd_mg_dl = ifelse(is.na(x3rd_urinary_albumin_mg_l), x3rd_urinary_albumin_g_l*100, x3rd_urinary_albumin_mg_l/10),
           urin_crea1st_g_dl = ifelse(is.na(x1st_urinary_creatinine_mg_dl), x1st_urinary_creatinine_mmol_l*0.01132, x1st_urinary_creatinine_mg_dl/1000),
           urin_crea2nd_g_dl = ifelse(is.na(x2nd_urinary_creatinine_mg_dl), x2nd_urinary_creatinine_mmol_l*0.01132, x2nd_urinary_creatinine_mg_dl/1000),
           urin_crea3rd_g_dl = ifelse(is.na(x3rd_urinary_creatinine_mg_dl), x3rd_urinary_creatinine_mmol_l*0.01132, x3rd_urinary_creatinine_mg_dl/1000))
  
  urin_alb_mg_dl <- df %>% 
    dplyr::select(urin_alb1st_mg_dl, urin_alb2nd_mg_dl, urin_alb3rd_mg_dl) %>%
    mutate(urin_alb_mg_dl = apply(., 1, function(x) first(na.omit(x)))) %>%
    dplyr::select(urin_alb_mg_dl)
  
  
  urin_crea_g_dl <- df %>% 
    dplyr::select(urin_crea1st_g_dl, urin_crea2nd_g_dl, urin_crea3rd_g_dl) %>%
    mutate(urin_crea_g_dl = apply(., 1, function(x) first(na.omit(x)))) %>%
    dplyr::select(urin_crea_g_dl)
  
  UACR_mgg <- (urin_alb_mg_dl/urin_crea_g_dl)[[1]]
  
  # Set bl uacr to 0.01 to allow log2 transformation due to skewness
  UACR_mgg_adj <- ifelse(as.numeric(UACR_mgg) == 0, 0.05, as.numeric(UACR_mgg))
  
  return(UACR_mgg_adj)
}

# -------------------------- DATA STORAGE -----------------

#' @title Wrapper to store list object as .json file
#' 
#' @details
#' Uses `jsonlite` package for conversion to minified JSON format. 
#' 
#' Converts formulas to strings. 
#' 
#' Cannot deal with all R types, only those used in this project.
saveJSON <- function(object, file) {
  # convert unmapped R types
  for (i in 1:length(object)) {
    if (is.formula(object[[i]]))
      object[[i]] <- deparse1(object[[i]])
  }
  
  write(jsonlite::toJSON(object, digits = NA), file)
}

#' @title Wrapper to read stored .json files
readJSON <- function(file) {
  jsonlite::read_json(file, simplifyVector = TRUE)
}

# -------------------------- MODELLING --------------------

adjusted_R2 <- function(pred, obs, N, k){
  r2 <- cor(pred,obs, use="pairwise.complete.obs")^2
  r2 <- 1-(((1-r2)*(N-1))/(N-k-1))
  return(r2)
}

c_index <- function(pred, obs){
  c.model <- concreg(data=data.frame(predicted=pred, observed=obs), observed~predicted, npar=TRUE)
  return(1-cindex(c.model))
}

eval_preds <- function(pred, obs, lmerObject){
  # x=data.test.list[[8]]; pred=x$pred; obs=x$FU_eGFR_epi; N=length(unique(fit.lme$data$PatID)); k=sum(anova(fit.lme)$numDF)
  
  df <- data.frame(pred=pred, obs=obs)
  df <- df[complete.cases(df),]
  
  if(sum(!is.na(obs))>100){
    r2mm <- r.squaredGLMM(lmerObject)
    r2 <- cor(pred,obs, use="pairwise.complete.obs")^2
    CS <- lm(df$obs ~ df$pred)$coefficients[2]
    Cind <- c_index(pred=df$pred, obs=df$obs)
    rmse.val <- RMSE(df$pred, df$obs)
    mae.val <- MAE(df$pred, df$obs)
    cil.val <- mean(df$obs)- mean(df$pred)
  }else{r2 <- r2mm <- CS <- rmse.val <- mae.val <- cil.val <- Cind<-NA}

  res <- data.frame(Nobs = sum(!is.na(df$obs)),
                    R2marg = r2mm[1],
                    R2cond = r2mm[2],
                    R2 = r2,
                    RMSE = rmse.val,
                    MAE = mae.val,
                    CalbinLarge = cil.val,
                    CalbSlope=CS,
                    C = Cind)
  return(res)
}

extractLMER <- function(lmerObject, idVar){
  
  #Formula
  formYx <- nobars(formula(lmerObject))
  
  # covariance matrix of random effect
  V.re <- VarCorr(lmerObject, )[[idVar]]
  
  # random effect correlation matrix
  C.re <- attr(V.re, "correlation")
  
  # estimated standard deviation of the errors
  sigma <- lmerObject@sigma   
  
  # fixed effects covariance matrix
  V.fe <- as.data.frame(as.matrix(vcov(lmerObject)))
  
  # fixed effects coefficients
  betas <- fixef(lmerObject)
  
  out <- list("form"=formYx,
              "VarRE"=V.re, "CorrRE"=C.re, "VarFE"=V.fe, 
              "betas"=betas, "betas_names"=names(betas),
              "sigma"=sigma)
  return(out)
}

draw_bootstrap_sample <- function(data.tmp){
  pats.boot <- data.tmp %>%
    filter(Time==0) %>%
    group_by(Country) %>%
    sample_frac(2/3, replace=T) %>%
    dplyr::select(PatID, Country) %>%
    data.frame()
  data.boot <- data.tmp[data.tmp$PatID %in% pats.boot$PatID,]
  data.boot$fold <- as.numeric(data.boot$Country)
  
  return(data.boot)
}


intext_crossvalidate <- function(data.boot, b.nr, return_preds=F){
  # data.boot=data.full; b.nr=NA; return_preds=T
  
  df.pred <- df.res <- df.re <- list()
  c = length(unique(data.full$Country))
  i=1
  for(i in 1:c){
    data.train <- data.boot[data.boot$fold != i, ] 
    data.test <- data.boot[data.boot$fold == i, ] 
    
    fit.lmer <- lmer(FU_eGFR_epi ~ (Time + Time  *(BL_age + BL_sex + BL_bmi + BL_smoking + BL_map + BL_hba1c + BL_serumchol +
                                                     BL_hemo + BL_uacr_log2 + BL_med_dm + BL_med_bp + BL_med_lipid) + (1|Country) + (1+Time|PatID)),
                     data=data.train, REML=F, control=lmerControl(optimizer="bobyqa"))
    
    # Update prediction with BL value
    data.test.t0 <- data.test[data.test$Time==0,]
    data.test.t0$Country <- "Unknown"
    res <- update_PredByBase(lmerObject=fit.lmer, 
                                newdata = data.test.t0, 
                                cutpoint = slope_cutpoint,
                                timeVar = "Time", idVar="PatID", idVar2="Country",
                                times =unique(data.boot$Time_cat)[-1], 
                                all_times=F)
    
    # Summarize and prepare output
    data.test$Time <- round(data.test$Time,0)
    data.test.new <- full_join(data.test, res$Pred[,c("PatID", "Time","prior.pred","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")], by=c("PatID", "Time"))
    data.test.new$fold <- i
    data.test.new$Country <- data.test$Country[1]
    data.test.new$Cohort <- data.test$Cohort[1]
    df.pred[[i]] <- data.test.new
    df.re[[i]] <- data.frame("PatID"= data.test.t0$PatID, "Fold"= i, res$RE.est)
    
    data.test.list <- split(data.test.new, as.factor(data.test.new$Time)) 
    res <- lapply(data.test.list, function(x) eval_preds(pred=x$pred, obs=x$FU_eGFR_epi, lmerObject=fit.lmer))
    ov.perf <- eval_preds(data.test.new[!data.test.new$Time==0,]$pred, data.test.new[!data.test.new$Time==0,]$FU_eGFR_epi, lmerObject = fit.lmer)
    df.res[[i]] <- data.frame("Boot"=b.nr,"Fold"=i,"Time"=c(100,unique(data.boot$Time)),rbind( ov.perf, do.call(rbind, res)))
  }
  
  if(return_preds){
    out <- list("pred"=do.call(rbind, df.pred), "res"=do.call(rbind, df.res), "re"=do.call(rbind, df.re))
  }else{
    out <- list("res"=do.call(rbind, df.res))
  }
  
  return(out)
}



# ------------------------------ PLOTS ------------------------

plot_calibration_cont <- function(yobs, yhat, time="Not specified!", cohort="dev",save=F, type="postUp",
                             out.path = "."){
  df <- data.frame(yobs=yobs, yhat=yhat)

  # Plot
  p <- ggplot(df, aes(x=yhat, y=yobs)) +
    geom_point() +
    scale_x_continuous(expand = c(0, 0), name = "Predicted", limits = c(0,150)) + 
    scale_y_continuous(expand = c(0, 0), name = "Observed", limits = c(0,150)) +
    geom_abline(intercept = 0) + 
    ggtitle(paste0("Time of follow-up (years) = ", time)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=18)) 

  if(save){ggsave(paste0(out.path, "fig_cal_",cohort,"_t",time,"_",type,".png"),width=8, height=8)}
  print(p)  
}

plot_calibration_bin <- function(pred, true, out.path, save=F){
  # Generate calibration curve plot for binary outcome (0,1)
  data <- data.frame("pred"=pred,"true"=true)
  
  df <- data.frame(obs_prob=rep(0,10), std_obs=rep(0,10),pred_prob=rep(0,10)) 
  qPS <- quantile(data$pred, seq(0,0.9,0.1))
  deciles <- sapply(1:nrow(data), function(x) sum(data$pred[x] >=qPS))
  
  # Mean Predicted Probability in each decile
  df$pred_prob <- sapply(1:10, function(x) mean(data[deciles == x,]$pred))
  
  # Relative Frequency of receiving CT in each decile + St.D
  df$obs_prob <- sapply(1:10, function(x) nrow(data[deciles == x & data$true == 1,])/nrow(data[deciles == x,]))
  df$std_obs <- sapply(1:10,function(x) sqrt((df$obs_prob[x]*(1-df$obs_prob[x]))/nrow(data[deciles == x & data$true == 1,])))
  
  p <- ggplot(df, aes(x=pred_prob, y=obs_prob,group=pred_prob)) +
    geom_point() +
    stat_boxplot(geom ='errorbar') +
    scale_x_continuous(expand = c(0, 0), limits = c(0,1.05), breaks=seq(0,1,0.2), name = "Predicted Probability") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.05), breaks=seq(0,1,0.2), name = "Observed Probability") +
    geom_abline(intercept = 0) + 
    geom_errorbar(aes(x=pred_prob, ymin=obs_prob-std_obs, ymax=obs_prob+std_obs), width=0.025)+
    theme_bw()+ 
    theme(axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  print(p)
  if(save){ggsave(paste0(out.path, "fig_cal_probability.png"), width = 20, height = 20, dpi = 640, units="cm", limitsize = F)}
  return(p)
}

plot_fun <- function( x, y, time=0 ){
  df <- data.frame(dep=as.numeric(unlist(data.full[data.full$Time_cat==time, x])),indep=as.numeric(unlist(data.full[data.full$Time_cat==time, y])))
  ggplot(df, aes(x=dep, y=indep)) + 
    geom_point() +
    scale_x_continuous(name = x)+
    scale_y_continuous(name = y) +
    theme_bw()
}


# --------------------- PREDICTION UPDATE WITH BASELINE VALUES  ---------------------



#' @title update_PredByBase for lme4 model objects
#' 
#' @details 
#' Ported from implementation in `JMBayes::indvPred_lme`. 

update_PredByBase <- function (lmerObject, newdata, timeVar, idVar, idVar2=NULL, times,
                                 level = 0.95, cutpoint=-3, all_times) 
{
  # Specify to try the function
  # lmerObject = fit.final
  # newdata = data.full.t0
  # timeVar = "Time"
  # idVar <- "PatID"
  # idVar2="Country"
  # times = seq(1,8,1)
  # all_times = F
  # level = 0.95
  # cutpoint=-3

  
  # ---- Assign elements of lmer to objects
  data <- lmerObject@frame
  # get fixed effects part of formula
  formYx <- nobars(formula(lmerObject))
  mfX <- model.frame(terms(formYx), data = data)
  TermsX <- attr(mfX, "terms")
  
  # get ids for random effects
  rf_spec = findbars(formula(lmerObject)) %>%
    sapply(deparse1) %>% 
    sapply(str_split, pattern = "\\|")
  rf = map_chr(rf_spec, ~ str_trim(.x[1])) %>% 
    set_names(map_chr(rf_spec, ~ str_trim(.x[2])))
  
  # get random effects part of formula for idVar
  formYz <- as.formula(sprintf("~ %s", rf[[idVar]]))
  mfZ <- model.frame(terms(formYz), data = data)
  TermsZ <- attr(mfZ, "terms")
  
  # random intercepts per idVar2
  formYz_2 <- as.formula(sprintf("~ %s", rf[[idVar2]]))
  mfZ_2 <- model.frame(terms(formYz_2), data = data)
  TermsZ_2 <- attr(mfZ_2, "terms")
  Z2.icpt <- ranef(lmerObject)[[idVar2]] %>% 
    tibble::rownames_to_column("country") %>% 
    add_row(country="Unknown", "(Intercept)" = 0)
  
  betas <- fixef(lmerObject)
  # estimated standard deviation of the errors
  sigma <- lmerObject@sigma   
  # random effects covariance matrixfor idVar
  G <- VarCorr(lmerObject, )[[idVar]]
  # random effect correlation matrix
  C <- attr(G, "correlation")
  # fixed effects covariance matrix
  V.fe <- as.matrix(vcov(lmerObject))
  
  all_vars <- unique(c(all.vars(TermsX), all.vars(TermsZ)))
  ind_cc = complete.cases(newdata[all_vars])
  if (sum(ind_cc) < nrow(newdata)) {
    print(sprintf("Using %d out of %d observations from newdata with complete data.", 
                  sum(ind_cc), nrow(newdata)))  
  }
  newdata_nomiss <- newdata[ind_cc, ]
  mfX_new <- model.frame(TermsX, data = newdata_nomiss)
  X_new <- model.matrix(formYx, mfX_new)
  mfZ_new <- model.frame(TermsZ, data = newdata_nomiss)
  Z_new <- model.matrix(formYz, mfZ_new)
  
  mfZ_2_new <- model.frame(TermsZ_2, data = newdata_nomiss)
  Z_2_new <- sapply(model.frame(formula("~ Country"), data=newdata_nomiss)[,1], 
                    function(x) Z2.icpt[Z2.icpt$country %in% x,2])
  
  na_ind <- attr(mfX_new, "na.action")
  y_new <- model.response(mfX_new, "numeric")
  
  id_nomiss <- match(newdata_nomiss[[idVar]], unique(newdata_nomiss[[idVar]]))
  n <- length(unique(id_nomiss))
  
  # --- Determine times for which to predict for
  times_orig <- data[[timeVar]]
  times_orig <- times_orig[!is.na(times_orig)]
  if (is.null(times) || !is.numeric(times)) {
    times <- seq(min(times_orig), max(times_orig), length.out = 100)
  }
  id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
  last_time <- tapply(newdata[[timeVar]], id, max)
  times_to_pred <- lapply(last_time, function(t) if (all_times) 
    times
    else times[times > t])
  id_pred <- rep(seq_len(n), sapply(times_to_pred, length))
  
  # ------ Create set with new rows for preds
  newdata_pred <- right_rows(newdata, newdata[[timeVar]], id, 
                             times_to_pred)
  newdata_pred[[timeVar]] <- unlist(times_to_pred)
  mfX_new_pred <- model.frame(TermsX, data = newdata_pred, 
                              na.action = NULL)
  X_new_pred <- model.matrix(formYx, mfX_new_pred)
  mfZ_new_pred <- model.frame(TermsZ, data = newdata_pred, na.action = NULL)
  Z_new_pred <- model.matrix(formYz, mfZ_new_pred)
  
  mfZ_2_new_pred <- model.frame(TermsZ_2, data = newdata_pred, na.action = NULL)
  Z_2_new_pred <- sapply(model.frame(formula("~ Country"), data=newdata_pred)[,1], function(x) Z2.icpt[Z2.icpt$country %in% x,2])
  
  # ---- Initial predictions
  pred_y_i0 <- c(X_new %*% betas) + Z_2_new
  pred_y_it <- c(X_new_pred %*% betas) + Z_2_new_pred

  # --------  Update predictions
  outcomeVar <- formula(lmerObject)[[2]]
  obs_y <- newdata[[outcomeVar]]
  b0.post <- (G[1,1]/(G[1,1] + sigma^2))*(obs_y-pred_y_i0)
  b1.post <- (b0.post* G[1,2])/(G[1,1])
  b.new <- data.frame("Intercept"=b0.post, "Time"=b1.post)
  
  std.error <- function(x) sd(x)/sqrt(length(x))
  se.b0 <- std.error(b0.post)
  se.b1 <- std.error(b1.post)
  
  RE_estimates <- data.frame("b0.est"=b0.post, "b0.se"=se.b0, 
                             "b0.lo"=b0.post-qnorm((1-level)/2,lower.tail=FALSE)*se.b0, 
                             "b0.up"=b0.post+qnorm((1-level)/2,lower.tail=FALSE)*se.b0, 
                             "b1.est"=b1.post, "b1.se"=se.b1, 
                             "b1.lo"=b1.post-qnorm((1-level)/2,lower.tail=FALSE)*se.b1, 
                             "b0.up"=b1.post+qnorm((1-level)/2,lower.tail=FALSE)*se.b1)
  
  # -------- Compute E(Y(t=time)|y_i0) + 95% - CI
  y_hat_0 <- c(X_new %*% betas) + rowSums(Z_new * b.new) + Z_2_new
  y_hat_time <- c(X_new_pred %*% betas) + rowSums(Z_new_pred * b.new[id_pred, , drop = FALSE]) + Z_2_new_pred
  
  # ------ eGFR slope per individual
  dyit_hat <- betas[str_detect(names(betas), paste0(timeVar,"$"))] + 
    c(X_new[,!str_detect(colnames(X_new_pred), paste0(timeVar,"|Intercept"))] %*% 
        betas[str_detect(names(betas), paste0(timeVar,":"))]) +
    c(b.new[,str_detect(colnames(b.new), "Time")])
  
  # ------ Prediction interval
  V.fe.dev <- V.fe[str_detect(rownames(V.fe), "Time"), str_detect(colnames(V.fe), "Time")]
  SE.yit <- SE.dyit <- checkV0i <- checkGstar <- list()
  for (i in seq_len(n)) {
    # Variance of yit
    se <- list()
    id_i <- id_nomiss == i
    X_new_id <- X_new[id_i, , drop = FALSE]
    
    V0i= X_new_id %*% tcrossprod(V.fe,X_new_id)
    diag1<-(G[1,1]/(G[1,1] + sigma^2))^2*V0i
    diag2<-(1-(C[1,2])^2)*G[2,2]
    off12<-(G[1,2]^2/(G[1,1]))^2*diag1
    G_star<-cbind(c(diag1,off12),c(off12,diag2))
    
    ### compute Var(Y(t=time))
    for (time in times_to_pred[[i]]) {
      z <- c(1, time)
      Vit_star <- V0i + (z %*% G_star%*% z) + sigma^2
      se[[time]] <- sqrt(Vit_star)
    }
    SE.yit[[i]] <- unlist(se)
    
    # Variance of dyit
    X_new_id_dev <- X_new_id[,!str_detect(colnames(X_new_pred), paste0(timeVar))] 
    V0i= X_new_id_dev %*% V.fe.dev %*% X_new_id_dev
    diag2<-(1-(C[1,2])^2)*G[2,2]

    se <- sqrt(V0i + diag2)
    SE.dyit[[i]] <- se
  }
  
  se.pred <- unlist(SE.yit)
  low <-  y_hat_time - qnorm((1-level)/2,lower.tail=FALSE) * se.pred
  upp <- y_hat_time + qnorm((1-level)/2,lower.tail=FALSE) * se.pred
  
  se.pred.dev <- unlist(SE.dyit)
  low.dev <-   dyit_hat - qnorm((1-level)/2,lower.tail=FALSE) * se.pred.dev
  upp.dev <-  dyit_hat + qnorm((1-level)/2,lower.tail=FALSE) * se.pred.dev
  
  # ------ Probability of Progression:
  prob.prog <- round(pnorm((cutpoint-dyit_hat)/se.pred.dev, lower.tail=T),4)
  
  # ------Output
  # Predictions + CI
  out_data <- rbind(newdata, newdata_pred)
  out_data$pred <- c(y_hat_0, y_hat_time)
  out_data$Time_cat <- out_data$Time
  out_data$prior.pred <- c(pred_y_i0, pred_y_it)
  out_data$pred.lo <- c(rep(NA, length(pred_y_i0)), low)
  out_data$pred.up <- c(rep(NA, length(pred_y_i0)), upp)
  
  times_rep <- c(sapply(times_to_pred, length))
  out_data$pred.slope <- c(dyit_hat, rep(dyit_hat, times_rep))
  out_data$pred.slope.lo <- c(low.dev, rep(low.dev, times_rep))
  out_data$pred.slope.up <- c(upp.dev, rep(upp.dev, times_rep))
  out_data$pred.prob <- c(prob.prog, rep(prob.prog, times_rep))
  
  out_data <- out_data[order(out_data[[idVar]], out_data[[timeVar]]),]
  
  res <- list("Pred" = out_data, "RE.est" = RE_estimates)
  return(res)
}

