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


# =============================== PLOTS ========================================

plot_trajectory <- function(x){

  p1 <- ggplot2::ggplot(x, ggplot2::aes(x=Time, y=pred)) +
    ggplot2::geom_point(size=3, shape=8) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous("Time", limits = c(0,7), breaks = seq(0,7,1)) +
    ggplot2::scale_y_continuous("Predicted eGFR") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size=16))
  
  return(p1)
}

smilegraph <- function(risk){
  happy<-png::readPNG("www/happy.png")
  sad<-png::readPNG("www/sad.png")
  plot(0:10,0:10,ty="n", xlab="", ylab="")
  
  num<-round(risk)
  
  if(num==0){
    for(i in 1:10){
      for(j in 1:10) rasterImage(happy,j-1,i-1, j,i)  
    }
    
  }else{
    numy<-ceiling(num/10)
    lastnum <- (num/10 - floor(num/10))*10
    lastnum[lastnum==0]<-10
    if(numy>1){ numx <- c(rep(10, numy-1),lastnum)}else{
      numx <- lastnum
    }
    apply(as.matrix(1:numy),1, function(i) {
      apply(as.matrix(1:numx[i]),1, function(j) rasterImage(sad,j-1,i-1,j,i))  
    })
    
    numsad<-round(100-risk)
    numsady<-ceiling(numsad/10)
    lastnumsad <- (numsad/10 - floor(numsad/10))*10
    lastnumsad[lastnumsad==0]<-10
    if(numsady>1){ numsadx <- c(rep(10, numsady-1),lastnumsad)}else{
      numsadx <- lastnumsad
    }
    apply(as.matrix(10:(11-numsady)),1, function(i){
      index<- 11-i
      apply(as.matrix(1:numsadx[index]),1,function(j) rasterImage(happy,10-j+1,i-1, 10-j,i) ) 
    })
  }
}

# ==================== PREDICTION UPDATE WITH BASELINE VALUES  =================

update_PredByBase <- function (lmerList, newdata, cutpoint=-3, times, level=0.95) 
{
  # ---- Specify relevant components
  idVar = "PatID"
  idVar2 = "Country"
  timeVar = "Time"
  times_to_pred=list("1"=times)
  formYx =as.formula(lmerList$form)
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
    se[[time]] <- sqrt(Vit_star)
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
  
  res <- out_data
  return(res)
}

