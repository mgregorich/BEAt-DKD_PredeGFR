# ----------------------------------------------------------------------
# Author: MG
# Date: 03.08.2021
# Info: Update longitudinal predictions by baseline measurements
# ---------------------------------------------------------------------



# --------------------------- GENERAL ---------------------


recodeDates <- function(x){
  x=str_replace(x, "Mrz", "März")
  x.new=if_else(str_detect(x, "^[:digit:]+$"), as.character(as.Date(as.numeric(x), origin = "1899-12-30")), as.character(as.Date(x, "%d%b%Y")))
  x.new = as.Date(x.new)
  return(x.new)
}

right_rows <- function (data, times, ids, Q_points) {
  fids <- factor(ids, levels = unique(ids))
  if (!is.list(Q_points))
    Q_points <- split(Q_points, row(Q_points))
  ind <- mapply(findInterval, Q_points, split(times, fids))
  ind[ind < 1] <- 1
  rownams_id <- split(row.names(data), fids)
  ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
  data[c(ind), ]
}


# -------------------------- MODELLING --------------------

c_index <- function(pred, obs){
  c.model <- concreg(data=data.frame(predicted=pred, observed=obs), observed~predicted, npar=TRUE)
  return(1-cindex(c.model))
}

eval_preds <- function(pred, obs, lmerObject){
  # x=data.test.list[[8]]; pred=x$pred; obs=x$FU_eGFR_epi; N=length(unique(fit.lme$data$PatID)); k=sum(anova(fit.lme)$numDF)
  
  df <- data.frame(pred=pred, obs=obs)
  df <- df[complete.cases(df),]
  
  r2 <- r.squaredGLMM(lmerObject)
  CS <- ifelse(all(is.na(df$obs)),NA, lm(df$obs ~ df$pred)$coefficients[2])
  
  if(sum(!is.na(obs))>25){
    Cind <- ifelse(all(is.na(df$obs)),NA,c_index(pred=df$pred, obs=df$obs))
  }else{Cind<-NA}

  res <- data.frame(Nobs = sum(!is.na(df$obs)),
                    R2marg = r2[1],
                    R2cond = r2[2],
                    RMSE = RMSE(df$pred, df$obs),
                    MAE = MAE(df$pred, df$obs),
                    CalbinLarge = mean(df$obs)- mean(df$pred),
                    CalbSlope=CS,
                    C = Cind)
  return(res)
}



# ------------------------------ PLOTS ------------------------

plot_calibration_cont <- function(yobs, yhat, fit=NULL, time="Not specified!", cohort="dev",save=F, 
                             out.path = "."){
  df <- data.frame(yobs=yobs, yhat=yhat)
  if(!is.null(fit)){res <- round(eval_preds(pred=yhat, obs = yobs, lmerObject=fit),3)}
  
  # Plot
  p <- ggplot(df, aes(x=yhat, y=yobs)) +
    geom_point() +
    scale_x_continuous(expand = c(0, 0), name = "Predicted", limits = c(0,150)) + 
    scale_y_continuous(expand = c(0, 0), name = "Observed", limits = c(0,150)) +
    geom_abline(intercept = 0) + 
    ggtitle(paste0("FU time = ", time)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=16)) 
  
  if(!save){
    p<- p +
    annotate("text", x = 25, y = 125, label = paste0("R2 = ",res$R2,
                                                     "\nC-index = ", res$C, 
                                                     "\nCalLarge = ", res$CalbinLarge, 
                                                     "\nCalSlope = ", res$CalbSlope))}
  print(p)
  if(save){ggsave(paste0(out.path, "fig_cal_",cohort,"_t",time,".png"),width=8, height=8)}
  # filename <- paste0("plot_calcurve_", mod.method)
  # png(file = paste(out.path, filename,".png", sep = ""), width = 5*600, height = 5*600, units = "px", res = 600)
  # plot_calibration_curve_ps
  # dev.off()
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



# --------------------- PREDICTION UPDATE WITH BASELINE VALUES  ---------------------



#' @title LongPred_ByBase for lme4 model objects
#' 
#' @details 
#' Ported from implementation in `JMBayes::indvPred_lme`. 
LongPred_ByBase_lmer <- function (lmerObject, newdata, timeVar, idVar, idVar2=NULL,  times,
                                 level = 0.95, cutpoint=-3, all_times, interval="prediction", M=500, seed=123) 
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
  # interval="prediction"
  # M=100
  # seed=123
  # 
  # lmeObject=risk_model; 
  # newdata = data.diacore.t0; 
  # cutpoint = slope_cutpoint;
  # timeVar = "Time"; idVar="PatID"; idVar2="Country";
  # times =sort(unique(data.diacore$Time_cat)); 
  # all_times=F
  
  
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
  D <- VarCorr(lmerObject)[[idVar]]
  # fixed effects covariance matrix
  V.fe <- as.matrix(vcov(lmerObject))
  
  all_vars <- unique(c(all.vars(TermsX), all.vars(TermsZ)))
  newdata_nomiss <- newdata[complete.cases(newdata[all_vars]), ]
  mfX_new <- model.frame(TermsX, data = newdata_nomiss)
  X_new <- model.matrix(formYx, mfX_new)
  mfZ_new <- model.frame(TermsZ, data = newdata_nomiss)
  Z_new <- model.matrix(formYz, mfZ_new)
  
  mfZ_2_new <- model.frame(TermsZ_2, data = newdata_nomiss)
  Z_2_new <- sapply(model.frame(formula("~ Country"), data=newdata_nomiss)[,1], function(x) Z2.icpt[Z2.icpt$country %in% x,2])
  
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
  
  # ------ Compute random coeffs with baseline value; assumes country is unknown
  b <- matrix(0.0, n, ncol(Z_new))
  post_vars <- DZtVinv <- vector("list", n)
  for (i in seq_len(n)) {
    id_i <- id_nomiss == i
    X_new_id <- X_new[id_i, , drop = FALSE]
    Z_new_id <- Z_new[id_i, , drop = FALSE]
    Z_2_new_id <- Z_2_new[id_i, drop = FALSE]
    
    Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) + sigma^2 * diag(sum(id_i)))
    DZtVinv[[i]] <- tcrossprod(D, Z_new_id) %*% Vi_inv
    b[i, ] <- c(DZtVinv[[i]] %*% (y_new[id_i] - (X_new_id %*% betas + Z_2_new_id)))
    t1 <- DZtVinv[[i]] %*% Z_new_id %*% D
    t2 <- DZtVinv[[i]] %*% X_new_id %*% V.fe %*% 
      crossprod(X_new_id, Vi_inv) %*% Z_new_id %*% D
    post_vars[[i]] <- D - t1 + t2
  }
  fitted_y <- c(X_new %*% betas) + rowSums(Z_new * b[id_nomiss, , drop = FALSE] + Z_2_new)
  
  
  # ---- Initial predictions
  pred_y_i0 <- c(X_new %*% betas) + Z_2_new
  pred_y_it <- c(X_new_pred %*% betas) + Z_2_new_pred
  
  # --------  Update predictions
  outcomeVar <- formula(lmerObject)[[2]]
  obs_y <- newdata[[outcomeVar]]
  b0.post <- (D[1,1]/(D[1,1] + sigma^2))*(obs_y-pred_y_i0)
  b1.post <- (b0.post* D[1,2])/(D[1,1])
  b.new <- data.frame("Intercept"=b0.post, "Time"=b1.post)
  
  # -------- Compute E(Y(t=time)|y_i0) + 95% - CI
  y_hat_time <- c(X_new_pred %*% betas) + rowSums(Z_new_pred * b.new[id_pred, , drop = FALSE]) + Z_2_new_pred
  
  std.error <- function(x) sd(x)/sqrt(length(x))
  se.b0 <- std.error(b0.post)
  se.b1 <- std.error(b1.post)
  
  RE_estimates <- data.frame("b0.est"=b0.post, "b0.se"=se.b0, 
                             "b0.lo"=b0.post-qnorm((1-level)/2,lower.tail=FALSE)*se.b0, 
                             "b0.up"=b0.post+qnorm((1-level)/2,lower.tail=FALSE)*se.b0, 
                             "b1.est"=b1.post, "b1.se"=se.b1, 
                             "b1.lo"=b1.post-qnorm((1-level)/2,lower.tail=FALSE)*se.b1, 
                             "b0.up"=b1.post+qnorm((1-level)/2,lower.tail=FALSE)*se.b1)
  
  
  # ------ eGFR slope per individual
  dyit_hat <- betas[str_detect(names(betas), paste0(timeVar,"$"))] + 
    c(X_new[,!str_detect(colnames(X_new_pred), paste0(timeVar,"|Intercept"))] %*% 
        betas[str_detect(names(betas), paste0(timeVar,":"))]) +
    c(b.new[,str_detect(colnames(b.new), "Time")])
  
  # ------- Confidence/Prediction interval for Y 
  if(interval=="confidence"){
    SE.yit <- list()
    for (i in seq_len(n)) {
      se <- list()
      id_i <- id_nomiss == i
      X_new_id <- X_new[id_i, , drop = FALSE]
      
      Vy_0i_hat = X_new_id %*% tcrossprod(V.fe, X_new_id)
      diag0 <- (D[1,1]/(D[1,1] + sigma^2)) * Vy_0i_hat
      diag1 <- (1-(D[1,2]/(D[1,1]*D[2,2]))) * D[2,2]
      off01 <- off10 <- (D[1,2]^2)/D[2,2]
      D_star <- matrix(c(diag0, off01, off10, diag1), ncol=2, byrow=T) # adjusted covariance matrix of random effects given known baseline value
      
      ### compute Var(Y(t=time))
      for (time in times_to_pred[[i]]) {
        z <- c(1, time)
        Vit_star <- ((X_new_id) %*% tcrossprod(V.fe, X_new_id)) + (z %*% (D_star %*% z)) + sigma^2
        se[[time]] <- sqrt(Vit_star)
      }
      SE.yit[[i]] <- unlist(se)
    }
    
    se.pred <- unlist(SE.yit)
    low <-  y_hat_time - qnorm((1-level)/2,lower.tail=FALSE) * se.pred
    upp <- y_hat_time + qnorm((1-level)/2,lower.tail=FALSE) * se.pred
  }else if(interval=="prediction"){
    set.seed(seed)
    betas_M <- MASS::mvrnorm(M, betas, V.fe)
    modes_fun <- function (betas) {
      t(mapply("%*%", DZtVinv, split(y_new - ((X_new %*% betas) + Z_2_new), id_nomiss)))
    }
    modes_M <- lapply(split(betas_M, row(betas_M)), modes_fun) # MxN
    matrix_row <- function (m, i) m[i, , drop = FALSE]
    modes_M <- lapply(seq_len(n), function (i) t(sapply(modes_M, matrix_row, i = i))) # NxM
    b_M <- modes_M
    for (i in seq_len(n)) {
      b_M[[i]] <- t(apply(modes_M[[i]], 1, MASS::mvrnorm, n = 1, Sigma = post_vars[[i]]))
    }
    n_pred <- length(y_hat_time)
    sampled_y <- sampled_dy <-matrix(0.0, n_pred, M)
    
    for (m in seq_len(M)) {
      betas_m <- betas_M[m, ]
      b_m <- t(sapply(b_M, function (x) x[m, ]))
      mean_m <- c(X_new_pred %*% betas_m) + 
        rowSums(Z_new_pred * b_m[id_pred, , drop = FALSE]) + Z_2_new_pred
      dmean_m <- betas_m[str_detect(names(betas_m), paste0(timeVar,"$"))] + 
        c(X_new_pred[,!str_detect(colnames(X_new_pred), paste0(timeVar,"|Intercept"))] %*% betas_m[str_detect(names(betas_m), paste0(timeVar,":"))]) + 
        c(b_m[id_pred,str_detect(colnames(b_m), timeVar)])
      
      sampled_y[, m] <- rnorm(n_pred, mean_m, sigma)
      sampled_dy[,m] <- dmean_m
    }
    low <- apply(sampled_y, 1, quantile, probs = (1 - level) / 2)
    upp <- apply(sampled_y, 1, quantile, probs = 1 - (1 - level) / 2)
    
    if(nrow(newdata)>1){
      low.dyi <- apply(sampled_dy[!duplicated(sampled_dy),],1, quantile, probs = (1 - level) / 2)
      upp.dyi <- apply(sampled_dy[!duplicated(sampled_dy),],1, quantile, probs = 1 - (1 - level) / 2)
      sd.dyi <- apply(sampled_dy[!duplicated(sampled_dy),],1, sd)
    }else{
      low.dyi <- apply(sampled_dy,1, quantile, probs = (1 - level) / 2)[1]
      upp.dyi <- apply(sampled_dy,1, quantile, probs = 1 - (1 - level) / 2)[1]
      sd.dyi <- apply(sampled_dy,1, sd)[1]
    }
  }
  
  # ------ Probability of Progression:
  
  prob.prog <- round(pnorm(cutpoint,dyit_hat,sd.dyi, lower.tail=T),4)
  
  
  # ------Output
  # Predictions + CI
  out_data <- rbind(newdata, newdata_pred)
  out_data$pred <- c(fitted_y, y_hat_time)
  out_data$Time_cat <- out_data$Time
  out_data$prior.pred <- c(pred_y_i0, pred_y_it)
  out_data$pred.low <- c(rep(NA, length(pred_y_i0)), low)
  out_data$pred.upp <- c(rep(NA, length(pred_y_i0)), upp)
  
  times_rep <- c(sapply(times_to_pred, length))
  out_data$slope <- c(dyit_hat, rep(dyit_hat, times_rep))
  out_data$slope.low <- c(low.dyi, rep(low.dyi, times_rep))
  out_data$slope.upp <- c(upp.dyi, rep(upp.dyi, times_rep))
  out_data$prob.prog <- c(prob.prog, rep(prob.prog, times_rep))
  
  out_data[order(out_data[[idVar]], out_data[[timeVar]]),]
  
  res <- list("Pred" = out_data, "RE.est" = RE_estimates)
  return(res)
}

LongPred_ByBase_lme <- function (lmeObject, newdata, timeVar, idVar, idVar2=NULL,  times,
                              level = 0.95, cutpoint=-3, all_times, interval="prediction", M=500, seed=123) 
{
  # Specify to try the function
  # lmeObject = fit.final
  # newdata = data.full.t0
  # timeVar = "Time"
  # idVar <- "PatID"
  # idVar2="Country"
  # times = seq(1,8,1)
  # all_times = F
  # level = 0.95
  # cutpoint=-3
  # interval="prediction"
  # M=100
  # seed=123
  # 
  # lmeObject=risk_model; 
  # newdata = data.diacore.t0; 
  # cutpoint = slope_cutpoint;
  # timeVar = "Time"; idVar="PatID"; idVar2="Country";
  # times =sort(unique(data.diacore$Time_cat)); 
  # all_times=F

  
  # ---- Assign elements of lme to objects
  data <- lmeObject$data
  formYx <- formula(lmeObject)
  mfX <- model.frame(terms(formYx), data = data)
  TermsX <- attr(mfX, "terms")
  
  formYz <- formula(lmeObject$modelStruct$reStruct[[idVar]])
  mfZ <- model.frame(terms(formYz), data = data)
  TermsZ <- attr(mfZ, "terms")
  
  formYz_2 <- formula(lmeObject$modelStruct$reStruct[[idVar2]])
  mfZ_2 <- model.frame(terms(formYz_2), data = data)
  TermsZ_2 <- attr(mfZ_2, "terms")
  Z2.icpt <- ranef(lmeObject)[[1]] %>% tibble::rownames_to_column("country") %>% add_row(country="Unknown", "(Intercept)"=0)

  betas <- fixef(lmeObject)
  sigma <- lmeObject$sigma   # estimated standard deviation of the errors
  D <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", sigma^2)[[idVar]]  # random effects matrix lmeObject$modelStruct$reStruct
  V.fe <- vcov(lmeObject) # fixed effects covariance matrix
  
  all_vars <- unique(c(all.vars(TermsX), all.vars(TermsZ)))
  newdata_nomiss <- newdata[complete.cases(newdata[all_vars]), ]
  mfX_new <- model.frame(TermsX, data = newdata_nomiss)
  X_new <- model.matrix(formYx, mfX_new)
  mfZ_new <- model.frame(TermsZ, data = newdata_nomiss)
  Z_new <- model.matrix(formYz, mfZ_new)
  
  mfZ_2_new <- model.frame(TermsZ_2, data = newdata_nomiss)
  Z_2_new <- sapply(model.frame(formula("~ Country"), data=newdata_nomiss)[,1], function(x) Z2.icpt[Z2.icpt$country %in% x,2])
  
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
  
  
  # ------ Compute random coeffs with baseline value; assumes country is unknown
  b <- matrix(0.0, n, ncol(Z_new))
  post_vars <- DZtVinv <- vector("list", n)
  for (i in seq_len(n)) {
    id_i <- id_nomiss == i
    X_new_id <- X_new[id_i, , drop = FALSE]
    Z_new_id <- Z_new[id_i, , drop = FALSE]
    Z_2_new_id <- Z_2_new[id_i, drop = FALSE]
    
    Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) + sigma^2 * diag(sum(id_i)))
    DZtVinv[[i]] <- tcrossprod(D, Z_new_id) %*% Vi_inv
    b[i, ] <- c(DZtVinv[[i]] %*% (y_new[id_i] - (X_new_id %*% betas + Z_2_new_id)))
    t1 <- DZtVinv[[i]] %*% Z_new_id %*% D
    t2 <- DZtVinv[[i]] %*% X_new_id %*% V.fe %*% 
      crossprod(X_new_id, Vi_inv) %*% Z_new_id %*% D
    post_vars[[i]] <- D - t1 + t2
  }
  fitted_y <- c(X_new %*% betas) + rowSums(Z_new * b[id_nomiss, , drop = FALSE] + Z_2_new)
  
  
  # ---- Initial predictions
  pred_y_i0 <- c(X_new %*% betas) + Z_2_new
  pred_y_it <- c(X_new_pred %*% betas) + Z_2_new_pred
  
  # ------  Update predictions 
  outcomeVar <- all.vars(lmeObject$terms)[1]
  obs_y <- newdata[[outcomeVar]]
  b0.post <- (D[1,1]/(D[1,1] + sigma^2))*(obs_y-pred_y_i0)
  b1.post <- (b0.post* D[1,2])/(D[1,1])
  b.new <- data.frame("Intercept"=b0.post, "Time"=b1.post)
  
  # -------- Compute E(Y(t=time)|y_i0) + 95% - CI
  y_hat_time <- c(X_new_pred %*% betas) + rowSums(Z_new_pred * b.new[id_pred, , drop = FALSE]) + Z_2_new_pred
  
  std.error <- function(x) sd(x)/sqrt(length(x))
  se.b0 <- std.error(b0.post)
  se.b1 <- std.error(b1.post)
  
  RE_estimates <- data.frame("b0.est"=b0.post, "b0.se"=se.b0, 
                             "b0.lo"=b0.post-qnorm((1-level)/2,lower.tail=FALSE)*se.b0, 
                             "b0.up"=b0.post+qnorm((1-level)/2,lower.tail=FALSE)*se.b0, 
                             "b1.est"=b1.post, "b1.se"=se.b1, 
                             "b1.lo"=b1.post-qnorm((1-level)/2,lower.tail=FALSE)*se.b1, 
                             "b0.up"=b1.post+qnorm((1-level)/2,lower.tail=FALSE)*se.b1)
  
  
  
  # ------ eGFR slope per individual
  dyit_hat <- betas[str_detect(names(betas), paste0(timeVar,"$"))] + c(X_new[,!str_detect(colnames(X_new_pred), paste0(timeVar,"|Intercept"))] %*% 
                                                              betas[str_detect(names(betas), paste0(timeVar,":"))]) + c(b.new[,str_detect(colnames(b.new), "Time")])
  
  
  # ------- Confidence/Prediction interval for Y 
  if(interval=="confidence"){
    SE.yit <- list()
    for (i in seq_len(n)) {
      se <- list()
      id_i <- id_nomiss == i
      X_new_id <- X_new[id_i, , drop = FALSE]
  
      Vy_0i_hat = X_new_id %*% tcrossprod(V.fe, X_new_id)
      diag0 <- (D[1,1]/(D[1,1] + sigma^2)) * Vy_0i_hat
      diag1 <- (1-(D[1,2]/(D[1,1]*D[2,2]))) * D[2,2]
      off01 <- off10 <- (D[1,2]^2)/D[2,2]
      D_star <- matrix(c(diag0, off01, off10, diag1), ncol=2, byrow=T) # adjusted covariance matrix of random effects given known baseline value
      
      ### compute Var(Y(t=time))
      for (time in times_to_pred[[i]]) {
        z <- c(1, time)
        Vit_star <- ((X_new_id) %*% tcrossprod(V.fe, X_new_id)) + (z %*% (D_star %*% z)) + sigma^2
        se[[time]] <- sqrt(Vit_star)
      }
      SE.yit[[i]] <- unlist(se)
    }
    
    se.pred <- unlist(SE.yit)
    low <-  y_hat_time - qnorm((1-level)/2,lower.tail=FALSE) * se.pred
    upp <- y_hat_time + qnorm((1-level)/2,lower.tail=FALSE) * se.pred
  }else if(interval=="prediction"){
    set.seed(seed)
    betas_M <- MASS::mvrnorm(M, betas, V.fe)
    modes_fun <- function (betas) {
      t(mapply("%*%", DZtVinv, split(y_new - ((X_new %*% betas) + Z_2_new), id_nomiss)))
    }
    modes_M <- lapply(split(betas_M, row(betas_M)), modes_fun) # MxN
    matrix_row <- function (m, i) m[i, , drop = FALSE]
    modes_M <- lapply(seq_len(n), function (i) t(sapply(modes_M, matrix_row, i = i))) # NxM
    b_M <- modes_M
    for (i in seq_len(n)) {
      b_M[[i]] <- t(apply(modes_M[[i]], 1, MASS::mvrnorm, n = 1, Sigma = post_vars[[i]]))
    }
    n_pred <- length(y_hat_time)
    sampled_y <- sampled_dy <-matrix(0.0, n_pred, M)

    for (m in seq_len(M)) {
      betas_m <- betas_M[m, ]
      b_m <- t(sapply(b_M, function (x) x[m, ]))
      mean_m <- c(X_new_pred %*% betas_m) + 
        rowSums(Z_new_pred * b_m[id_pred, , drop = FALSE]) + Z_2_new_pred
      dmean_m <- betas_m[str_detect(names(betas_m), paste0(timeVar,"$"))] + 
        c(X_new_pred[,!str_detect(colnames(X_new_pred), paste0(timeVar,"|Intercept"))] %*% betas_m[str_detect(names(betas_m), paste0(timeVar,":"))]) + 
        c(b_m[id_pred,str_detect(colnames(b_m), timeVar)])
      
      sampled_y[, m] <- rnorm(n_pred, mean_m, lmeObject$sigma)
      sampled_dy[,m] <- dmean_m
    }
    low <- apply(sampled_y, 1, quantile, probs = (1 - level) / 2)
    upp <- apply(sampled_y, 1, quantile, probs = 1 - (1 - level) / 2)
    
    if(nrow(newdata)>1){
      low.dyi <- apply(sampled_dy[!duplicated(sampled_dy),],1, quantile, probs = (1 - level) / 2)
      upp.dyi <- apply(sampled_dy[!duplicated(sampled_dy),],1, quantile, probs = 1 - (1 - level) / 2)
      sd.dyi <- apply(sampled_dy[!duplicated(sampled_dy),],1, sd)
    }else{
      low.dyi <- apply(sampled_dy,1, quantile, probs = (1 - level) / 2)[1]
      upp.dyi <- apply(sampled_dy,1, quantile, probs = 1 - (1 - level) / 2)[1]
      sd.dyi <- apply(sampled_dy,1, sd)[1]
    }

    
  }
  
  # ------ Probability of Progression:

  prob.prog <- round(pnorm(cutpoint,dyit_hat,sd.dyi, lower.tail=T),4)

  
  # ------Output
  # Predictions + CI
  out_data <- rbind(newdata, newdata_pred)
  out_data$pred <- c(fitted_y, y_hat_time)
  out_data$Time_cat <- out_data$Time
  out_data$prior.pred <- c(pred_y_i0, pred_y_it)
  out_data$pred.low <- c(rep(NA, length(pred_y_i0)), low)
  out_data$pred.upp <- c(rep(NA, length(pred_y_i0)), upp)
  
  times_rep <- c(sapply(times_to_pred, length))
  out_data$slope <- c(dyit_hat, rep(dyit_hat, times_rep))
  out_data$slope.low <- c(low.dyi, rep(low.dyi, times_rep))
  out_data$slope.upp <- c(upp.dyi, rep(upp.dyi, times_rep))
  out_data$prob.prog <- c(prob.prog, rep(prob.prog, times_rep))
  
  out_data[order(out_data[[idVar]], out_data[[timeVar]]),]
  
  res <- list("Pred" = out_data, "RE.est" = RE_estimates)
  return(res)
}
