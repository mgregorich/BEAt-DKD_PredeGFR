# ============================================================================ #
# Author: MG
# Date: 03.08.2021
# Info: Update longitudinal predictions by baseline measurements
# ============================================================================ #



# =============================== GENERAL ======================================

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

calc_map <- function(sys, dia){
  map <- (sys+ 2*dia)/3
  return(map)
}

# ============================== DATA STORAGE ==============================

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

# =============================== PLOTS ========================================

plot_eGFRslopeDistribution <- function(distr.slopes, dens.slopes, pred.slope){
  
  p <- ggplot(dens.slopes, aes(x=x, y=y)) +
    geom_line() +
    scale_y_continuous("Density", expand=c(0,0), lim=c(0,0.7), breaks = seq(0,0.7, 0.1)) +
    scale_x_continuous("eGFR slope",lim=c(-5,2.5), breaks = seq(-5,2.5,1)) +
    geom_area(alpha=0.1) +
    geom_vline(xintercept = pred.slope, col="red3", linetype="dashed") +
    geom_text(aes(x=pred.slope, label="Predicted slope", y=0.65), colour="red3", size=5) +
    geom_segment(x = 0, y = 0, xend = 0, yend = 0.7, color = 1,
                 arrow = arrow(length=unit(0.3,"cm"),)) +
    theme_minimal() +
    theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.2, "cm"),ends = "both")),
          text=element_text(size=18),
          panel.grid.minor=element_blank())
  p
  prob <- distr.slopes[which.min(abs(distr.slopes$Value-pred.slope)),]$quantil
  
  return(list("plot"=p, "prob"=prob))
}

plot_trajectory <- function(x){
  
  p1 <- ggplot2::ggplot(x, ggplot2::aes(x=Time, y=pred)) +
    ggplot2::geom_point(ggplot2::aes(colour="Predicted eGFR"),size=5, shape=8)  +
    ggplot2::geom_point(ggplot2::aes(x=0, y=BaseGFR[1], colour="Observed eGFR"), shape=8, size=5) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual("",values = c("Predicted eGFR" = "black", "Observed eGFR" = "red3")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pred.lo, ymax = pred.up), 
                         alpha = 0.1) +
    ggplot2::scale_x_continuous("Follow-up time (in years)", 
                                limits = c(0, 5), breaks = seq(0, 10, 1)) +
    ggplot2::scale_y_continuous(expression(paste("eGFR (mL/min/1.73 ", m^2, ")"))) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size=20), legend.position = "top")
  p1
  return(p1)
}


# ==================== PREDICTION UPDATE WITH BASELINE VALUES  =================

update_PredByBase <- function (lmerList, newdata, cutpoint=-3, level=0.95, baseline_eGFR=T) 
{
  # ---- Specify relevant components
  idVar = "PatID"
  idVar2 = "Country"
  timeVar = "Time"
  times_to_pred=list("1"=0:5)
  formYx = as.formula(lmerList$form)
  outcomeVar = stringr::str_extract(formYx, "[^~]+")[2]
  betas=lmerList$betas
  G=lmerList$VarRE
  C=lmerList$CorrRE
  sigma=lmerList$sigma
  V.fe=as.matrix(lmerList$VarFE)
  
  # ---- Assign elements of lmer to objects

  # get fixed effects part of formula
  mfX <- model.frame(terms(formYx), data = newdata)
  TermsX <- attr(mfX, "terms")
  X_new <- model.matrix(formYx, mfX)
  
  # ------ Create set with new rows for preds
  newdata_pred <- right_rows(newdata, newdata[[timeVar]], 1, times_to_pred)
  newdata_pred[[timeVar]] <- unlist(times_to_pred)
  mfX_new_pred <- model.frame(TermsX, data = newdata_pred, na.action = NULL)
  X_new_pred <- model.matrix(formYx, mfX_new_pred)
  
  
  formYz=as.formula(paste0("~1 + ", timeVar))
  mfZ <- model.frame(terms(formYz), data = newdata)
  TermsZ <- attr(mfZ, "terms")
  Z_new <- model.matrix(formYz, mfZ)
  
  mfZ_new_pred <- model.frame(TermsZ, data = newdata_pred, na.action = NULL)
  Z_new_pred <- model.matrix(formYz, mfZ_new_pred)
  
  
  # ---- Initial predictions
  pred_y_i0 <- c(X_new %*% betas) 
  pred_y_it <- c(X_new_pred %*% betas)
  
  # --------  Update predictions
  obs_y <- newdata[[outcomeVar]]
  b0.post <- (G[1,1]/(G[1,1] + sigma^2))*(obs_y-pred_y_i0)
  b1.post <- (b0.post* G[1,2])/(G[1,1])
  b.new <- data.frame("Intercept"=b0.post, "Time"=b1.post)
  
  # -------- Compute E(Y(t=time)|y_i0) + 95% - CI
  y_hat_0 <- c(X_new %*% betas) + rowSums(Z_new * b.new) 
  y_hat_time <- c(X_new_pred %*% betas) + rowSums(Z_new_pred * b.new) 
  
  # ------ eGFR slope per individual
  dyit_hat <- betas[stringr::str_detect(names(betas), paste0(timeVar,"$"))] + 
    c(X_new[,!stringr::str_detect(colnames(X_new_pred), paste0(timeVar,"|Intercept"))] %*% 
        betas[stringr::str_detect(names(betas), paste0(timeVar,":"))]) +
    c(b.new[,stringr::str_detect(colnames(b.new), "Time")])
  
  # ------ Prediction interval
  V.fe.dev <- V.fe[stringr::str_detect(rownames(V.fe), "Time"), stringr::str_detect(colnames(V.fe), "Time")]
  
  V0i= X_new %*% tcrossprod(V.fe,X_new)
  diag1<-(G[1,1]/(G[1,1] + sigma^2))^2*V0i
  diag2<-(1-(C[1,2])^2)*G[2,2]
  off12<-(G[1,2]^2/(G[1,1]))^2*diag1
  G_star<-cbind(c(diag1,off12),c(off12,diag2))
  
  ### compute Var(Y(t=time))
  se <- list()
  for (time in unlist(times_to_pred)) {
    z <- c(1, time)
    Vit_star <- V0i + (z %*% G_star%*% z) + sigma^2
    se[[time+1]] <- sqrt(Vit_star)
  }
  SE.yit <- unlist(se)
  
  # Variance of the deviation (dev) dyit
  X_new_dev <- X_new[,!stringr::str_detect(colnames(X_new_pred), paste0(timeVar))] 
  V0i= X_new_dev %*% V.fe.dev %*% X_new_dev
  diag2<-(1-(C[1,2])^2)*G[2,2]
  
  se <- sqrt(V0i + diag2)
  SE.dyit <- se
  
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
  out_data <- data.frame(newdata_pred, "prior.pred"=pred_y_it, "pred"=y_hat_time, "pred.lo"=low, "pred.up"=upp,
                         "pred.slope"=dyit_hat, "pred.slope.lo"=low.dev, "pred.slope.up"=upp.dev,
                         "pred.prob"=prob.prog)
  out_data$Time_cat <- out_data$Time
  out_data <- out_data[order(out_data[[idVar]], out_data[[timeVar]]),]
  return(out_data)
}

