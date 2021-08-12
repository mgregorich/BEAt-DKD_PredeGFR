# ----------------------------------------------------------------------
# Author: MG
# Date: 03.08.2021
# Info: Update longitudinal predictions by baseline measurements
# ---------------------------------------------------------------------


c_index <- function(pred, obs){
  c.model <- concreg(data=data.frame(predicted=pred, observed=obs), observed~predicted, npar=TRUE)
  return(1-cindex(c.model))
}

eval_preds <- function(pred, obs){
  df <- data.frame(pred=pred, obs=obs)
  df <- df[complete.cases(df),]
  r2 <- cor(pred,obs)^2
  res <- data.frame(R2 = r2,
                    RMSE = RMSE(df$pred, df$obs),
                    MAE = MAE(df$pred, df$obs),
                    CalbinLarge = mean(df$obs)- mean(df$pred),
                    CalbSlope=lm(df$pred ~ df$obs)$coefficients[2],
                    C = c_index(pred=df$pred, obs=df$obs))
  return(res)
}

plot_calibration <- function(yobs, yhat, fit=NULL, time="Not specified!"){
  df <- data.frame(yobs=yobs, yhat=yhat)
  res <- round(eval_preds(pred=yhat, obs = yobs),3)
  
  # Plot
  plot_calibration_curve_ps <- ggplot(df, aes(x=yhat, y=yobs)) +
    geom_point() +
    scale_x_continuous(expand = c(0, 0), name = "Predicted", limits = c(0,150)) + 
    scale_y_continuous(expand = c(0, 0), name = "Observed", limits = c(0,150)) +
    geom_abline(intercept = 0) + 
    ggtitle(paste0("FU time = ", time)) +
    theme_bw() +
    annotate("text", x = 25, y = 125, label = paste0("R2 = ",res$R2,
                                                     "\nC-index = ", res$C, 
                                                     "\nCalLarge = ", res$CalbinLarge, 
                                                     "\nCalSlope = ", res$CalbSlope))
  plot_calibration_curve_ps
  
  # filename <- paste0("plot_calcurve_", mod.method)
  # png(file = paste(out.path, filename,".png", sep = ""), width = 5*600, height = 5*600, units = "px", res = 600)
  # plot_calibration_curve_ps
  # dev.off()
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



########################################################################################################################

LongPred_ByBase <- function (lmeObject, newdata, timeVar, idVar, idVar2=NULL,  times,
                              level = 0.95, cutpoint=-3, all_times, interval="prediction", M=100, seed=123) 
{
  # Specify to try the function
  lmeObject = fit.final
  newdata = data.full.t0
  timeVar = "Time"
  idVar <- "PatID"
  idVar2="Country"
  times = unique(data.full$Time)
  all_times = FALSE
  level = 0.95
  cutpoint=-3
  interval="prediction"
  M=100
  seed=123

  # ---- Assign elements of lme to objects
  data <- lmeObject$data
  formYx <- formula(lmeObject)
  mfX <- model.frame(terms(formYx), data = data)
  TermsX <- attr(mfX, "terms")
  
  formYz <- formula(lmeObject$modelStruct$reStruct[[idVar]])
  mfZ <- model.frame(terms(formYz), data = data)
  TermsZ <- attr(mfZ, "terms")
  add.icpt <- ranef(lmeObject)[[1]] %>% tibble::rownames_to_column("country")

  betas <- fixef(lmeObject)
  sigma <- lmeObject$sigma   # estimated standard deviation of the errors
  D <- pdMatrix(lmeObject$modelStruct$reStruct[[idVar]]) * sigma^2  # random effects matrix lmeObject$modelStruct$reStruct
  V.fe <- vcov(lmeObject) # fixed effects covariance matrix
  
  all_vars <- unique(c(all.vars(TermsX), all.vars(TermsZ)))
  newdata_nomiss <- newdata[complete.cases(newdata[all_vars]), ]
  mfX_new <- model.frame(TermsX, data = newdata_nomiss)
  X_new <- model.matrix(formYx, mfX_new)
  mfZ_new <- model.frame(TermsZ, data = newdata_nomiss)
  Z_new <- model.matrix(formYz, mfZ_new)
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
  
  # ------ Compute random coeffs with baseline value
  b <- matrix(0.0, n, ncol(Z_new))
  post_vars <- DZtVinv <- vector("list", n)
  for (i in seq_len(n)) {
    id_i <- id_nomiss == i
    X_new_id <- X_new[id_i, , drop = FALSE]
    Z_new_id <- Z_new[id_i, , drop = FALSE]
    Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) + sigma^2 * diag(sum(id_i)))
    DZtVinv[[i]] <- tcrossprod(D, Z_new_id) %*% Vi_inv
    b[i, ] <- c(DZtVinv[[i]] %*% (y_new[id_i] - X_new_id %*% betas))
    t1 <- DZtVinv[[i]] %*% Z_new_id %*% D
    t2 <- DZtVinv[[i]] %*% X_new_id %*% V.fe %*% 
      crossprod(X_new_id, Vi_inv) %*% Z_new_id %*% D
    post_vars[[i]] <- D - t1 + t2
  }
  fitted_y <- c(X_new %*% betas) + rowSums(Z_new * b[id_nomiss, , drop = FALSE])
  
  # ---- Initial predictions
  pred_y_i0 <- c(X_new %*% betas)
  pred_y_it <- c(X_new_pred %*% betas)
  
  # ------  Update predictions ------
  outcomeVar <- all.vars(lmeObject$terms)[1]
  obs_y <- newdata[[outcomeVar]]
  b0.post <- (D[1,1]/(D[1,1] + sigma^2))*(obs_y-pred_y_i0)
  b1.post <- (b0.post* D[1,2])/(D[1,1])
  b.new <- data.frame("Intercept"=b0.post, "Time"=b1.post)
  
  # -------- Compute E(Y(t=time)|y_i0) + 95% - CI
  y_hat_time <- c(X_new_pred %*% betas) + rowSums(Z_new_pred * b.new[id_pred, , drop = FALSE])
  
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
  dyit_hat <- betas[str_detect(names(betas), paste0(timeVar,"$"))] + c(X_new[,!str_detect(colnames(X_new_pred), "Time|Intercept")] * 
                                                              betas[str_detect(names(betas), paste0(timeVar,":"))]) + c(b.new[,str_detect(colnames(b.new), timeVar)])
  
  
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
      t(mapply("%*%", DZtVinv, split(y_new - X_new %*% betas, id_nomiss)))
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
        rowSums(Z_new_pred * b_m[id_pred, , drop = FALSE])
      dmean_m <- betas_m[str_detect(names(betas_m), paste0(timeVar,"$"))] + 
        c(X_new_pred[,!str_detect(colnames(X_new_pred), "Time|Intercept")] * betas_m[str_detect(names(betas_m), paste0(timeVar,":"))]) + 
        c(b_m[id_pred,str_detect(colnames(b_m), timeVar)])
      
      sampled_y[, m] <- rnorm(n_pred, mean_m, lmeObject$sigma)
      sampled_dy[,m] <- dmean_m
    }
    low <- apply(sampled_y, 1, quantile, probs = (1 - level) / 2)
    upp <- apply(sampled_y, 1, quantile, probs = 1 - (1 - level) / 2)
    
    low.dy <- apply(sampled_dy[!duplicated(sampled_dy),],1, quantile, probs = (1 - level) / 2)
    upp.dy <- apply(sampled_dy[!duplicated(sampled_dy),],1, quantile, probs = 1 - (1 - level) / 2)
    se.dy <- apply(sampled_dy[!duplicated(sampled_dy),],1, std.error)
    
  }
  
  # ------ Probability of Progression:
  prob.prog <- round(pnorm((cutpoint-dyit_hat)/se.dy, lower.tail=FALSE),4)
  
  
  # ------Output
  # Predictions + CI
  out_data <- rbind(newdata, newdata_pred)
  out_data$pred <- c(fitted_y, y_hat_time)
  out_data$pred.low <- c(rep(NA, length(pred_y_i0)), low)
  out_data$pred.upp <- c(rep(NA, length(pred_y_i0)), upp)
  out_data$slope <- c(dyit_hat, rep(dyit_hat, sapply(times_to_pred, length)))
  out_data$slope.low <- c(low.dy, rep(low.dy, sapply(times_to_pred, length)))
  out_data$slope.upp <- c(upp.dy, rep(upp.dy, sapply(times_to_pred, length)))
  out_data$prob.prog <- c(prob.prog, rep(prob.prog, sapply(times_to_pred, length)))
  
  out_data[order(out_data[[idVar]], out_data[[timeVar]]),]
  
  res <- list("Pred" = out_data, "RE.est" = RE_estimates)
  return(res)
}
