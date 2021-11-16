#################################
# Author:MG
# Date: 07/10/2021
# Info: server file - shiny
################################



# data.new <- data.frame(PatID="Pnew", Time=0, Country="Unknown", FU_eGFR_epi=55,
#                        BL_age=69, BL_sex=1, BL_bmi=24, BL_smoking=1,BL_hemo=14, BL_hba1c=62,BL_serumchol=130, BL_map=calc_map(sys=125,dia=85),
#                        BL_uacr_log2=log(10,2),BL_med_dm=1, BL_med_bp=1, BL_med_lipid=1)


pacman::p_load(shiny, shinyjs, shinythemes, nlme, ggplot2, reshape2, dplyr, tidyverse, png)
source("functions_aux.R")

# ---------------------------- SHINY Server --------------------------
risk_model <- readRDS("riskpred_model.rds")


plot_trajectory <- function(x){
  #df.melt <- melt(x[,c("PatID","Time", "FU_eGFR_epi","pred", "pred.low", "pred.upp")], id.vars = c("PatID","Time", "pred.low", "pred.upp"))

  p1 <- ggplot(x, aes(x=Time, y=pred)) +
    geom_point(size=3, shape=8) +
    geom_line() +
    scale_x_continuous("Time", limits = c(0,7), breaks = seq(0,7,1)) +
    scale_y_continuous("Predicted eGFR") +
    theme_bw() +
    theme(text = element_text(size=16))
  
  return(p1)
}

calc_map <- function(sys, dia){
  map <- (sys+ 2*dia)/3
  return(map)
}




smilegraph <- function(risk){
  happy<-readPNG("happy.png")
  sad<-readPNG("sad.png")
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

shinyServer(function(input, output, session) {
  

  observeEvent(input$add_pred, {
    if(input$add_pred == 2){
      shinyjs::show(id = "lab1")
    }else{
      shinyjs::hide(id = "lab1")}})
  
  
  inputdata <- eventReactive(input$goButton ,{
    if(input$add_pred==2){
      data.new <- data.frame(PatID="Pnew", Time=0, Country="Unknown", FU_eGFR_epi=input$BL_eGFR,
                             BL_age=input$BL_age, 
                             BL_sex=as.numeric(input$BL_sex), 
                             BL_bmi=input$BL_bmi, 
                             BL_smoking=as.numeric(input$BL_smoking),
                             BL_hemo=input$BL_hemo, BL_hba1c=input$BL_hba1c, BL_serumchol=input$BL_serumchol, BL_map=calc_map(sys=input$BL_sysbp,dia=input$BL_diabp),
                             BL_uacr_log2=log(input$BL_uacr,2),
                             BL_med_dm=input$BL_med_dm*1, BL_med_bp=input$BL_med_bp*1, BL_med_lipid=input$BL_med_lipid*1)
    }else{
      data.new <- data.frame(PatID="Pnew", Time=0, Country="Unknown", FU_eGFR_epi=input$BL_eGFR,
                             BL_age=input$BL_age, 
                             BL_sex=as.numeric(input$BL_sex), 
                             BL_bmi=input$BL_bmi, 
                             BL_smoking=as.numeric(input$BL_smoking),
                             BL_hemo=14.5, BL_hba1c=input$BL_hba1c, BL_serumchol=202, BL_map=calc_map(sys=139,dia=76),
                             BL_uacr_log2=log(41,2),
                             BL_med_dm=input$BL_med_dm*1, BL_med_bp=input$BL_med_bp*1, BL_med_lipid=input$BL_med_lipid*1)
    }
    
    
    data.new
  })
  
  
  calc_out <- eventReactive(input$goButton,{
    
    data.new=inputdata()
    slope.cutpoint <- input$cutpoint
    
    res <- LongPred_ByBase(risk_model, newdata = data.new,  
                           cutpoint = slope.cutpoint,
                           timeVar = "Time", idVar="PatID", idVar2="Country",
                           times = seq(1,7,1), 
                           all_times=T)
    data.long <- data.new[rep(1,8),]
    data.long$Time <- seq(0,7,1)
    data.pred <- full_join(data.long,
                           res$Pred[,c("PatID", "Time","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")],
                           by=c("PatID", "Time"))
    data.pred[data.pred$Time==0,]$pred <- input$BL_eGFR
    data.pred <- data.frame(data.pred[,c("PatID", "Time","pred", "pred.low", "pred.upp", "slope", "slope.low", "slope.upp","prob.prog")])
    data.pred
    
      })
  
  makeOutputTable <- eventReactive(input$goButton, {
    odata <- calc_out()
    odata <- data.frame(odata[,-1])
    
    odata$pred.CI <- paste0("(",round(odata$pred.low,2),", ",round(odata$pred.upp,2), ")")
    odata$slope.CI <- paste0("(",round(odata$slope.low,2),", ",round(odata$slope.upp,2), ")")
    
    pred <- data.frame("Time"=odata$Time, "Prediction"=odata$pred, "CI"=odata$pred.CI)
    pred$CI[1] <- ""
    colnames(pred) <- c("FU Time", "Prediction","95% Confidence interval")
   
    slope <- data.frame("Slope"=odata$slope[1], "CI"=odata$slope.CI[1])
    colnames(slope) <- c("eGFR slope", "95% Confidence interval")

     return(list("pred"=pred, "slope"=slope, "prob"=odata$prob.prog[1]))
  })
  
  makeInputTable <- eventReactive(input$goButton, {
    idata <- inputdata()
    idata <- data.frame(idata)[,-c(1,2,3)]
    idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")] <- apply(idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")],2, function(x) ifelse(x==1, "yes", "no"))
    idata$BL_smoking <- ifelse(idata$BL_smoking==1, "ever", "never")
    idata$BL_sex <- ifelse(idata$BL_sex==1, "female", "male")
    
    if(input$add_pred==2){
      idemo <- data.frame("General"=c("Baseline eGFR", "Age","Sex", "BMI", "Smoking"), "Values"=unlist(idata[,c("FU_eGFR_epi","BL_age", "BL_sex", "BL_bmi", "BL_smoking")]))
      imed <- data.frame("Medication"=c("Glucose-lowering med.", "Blood pressure-lowering med.", "Lipid-lowering med."), "Values"=unlist(idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")]))
      ilab <- data.frame("Laboratory"=c("Hemoglobin", "Hba1C", "Serum Cholesterol", "MAP", "log2UACR"), "Values"=unlist(idata[,c("BL_hemo", "BL_hba1c", "BL_serumchol", "BL_map", "BL_uacr_log2")]))

      return(list("idemo"=idemo,"ilab"=ilab,"imed"=imed))
    }
    else{
      idemo <- data.frame("General"=c("Baseline eGFR", "Age","Sex", "BMI", "Smoking"), "Values"=unlist(idata[,c("FU_eGFR_epi","BL_age", "BL_sex", "BL_bmi", "BL_smoking")]))
      imed <- data.frame("Medication"=c("Glucose-lowering med.", "Blood pressure-lowering med.", "Lipid-lowering med."), "Values"=unlist(idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")]))
      
      return(list("idemo"=idemo,"imed"=imed))
    }
    
  })
  
  # ---------- Execute --------------------------------------------------------
  text_risk1 <- eventReactive(input$goButton, {  
    paste("Based on the provided information and a slope cutpoint of ", input$cutpoint,"mL/min/1.73m2, 
          the probability of the subject to progress into the group of rapid decline in renal function is ",  
          round(calc_out()[1,"prob.prog"]*100,2), "%.", sep="" )
    })
  
  text_risk2 <- eventReactive(input$goButton, {  
    paste("This means that in 100 people with these characteristics ",round(calc_out()[1,"prob.prog"]*100,0)," are expected to exhibit fast progression of renal decline in the following years.", sep="" )
  })
  
  text_longitudinal <- eventReactive(input$goButton, {  
    paste("The figure below illustrates the estimated longitudinal trajectory of the patient's eGFR measurements given the their information.", sep="" )
  })
  
  plotsmile <- eventReactive(input$goButton, {smilegraph(round(calc_out()[1,"prob.prog"]*100,2))})
  
  text_data <- eventReactive(input$goButton,{
    if(input$add_pred==2){
      HTML(paste0("The tables below contain both the patient information given and the predictions resulting from the model, including 95% prediction intervals for the individual time points. The patient's probability of progression to fast renal decline is also given."))
      
    }else{
      HTML(paste("The tables below contain both the patient information given and the predictions resulting from the model, including 95% prediction intervals for the individual time points. The patient's probability of progression to fast renal decline is also given.", 
            "By selecting the simple model, the average laboratory measurements of the development cohort were used to calculate the prediction values and the risk estimate. This can skew the accuracy and precision of the estimates.", 
            sep="<br/>"))
    }
  })
  
  # ------- Output ------------------------------
  output$text_risk1 <- renderText({
    text_risk1()
  })
  
  output$text_risk2 <- renderText({
    text_risk2()
  })
  
  output$text_data <- renderText({
    text_data()
  })
  
  output$text_longitudinal <- renderText({
    text_longitudinal()
  })
  
  output$plot_risk <- renderPlot({ 
    plotsmile()})
  
  
  output$plot_trajectory <- renderPlot({
    data.pred <- calc_out()
    plot_trajectory(data.pred)
  })
  
  output$table_new1 <- renderTable({
    makeInputTable()$idemo
  })
  
  output$table_new2 <- renderTable({
    if(input$add_pred==2){
      makeInputTable()$ilab}
    else{return()}
  })
  
  output$table_new3 <- renderTable({
    makeInputTable()$imed
  })
  
  output$table_pred1 <- renderTable({
    makeOutputTable()$pred
  })
  
  output$table_pred2 <- renderTable({
    makeOutputTable()$slope
  })



  
  # Reset
  observeEvent(input$reset_input, {
    updateNumericInput(session, "BL_age", value = 67)
    updateRadioButtons(session, "BL_sex", selected = "female")
    updateNumericInput(session, "BL_bmi", value=25)
    updateRadioButtons(session, "BL_smoking", selected = "never")
    
    updateNumericInput(session, "BL_hba1c", value = 62)
    updateNumericInput(session, "BL_serumchol", value = 130)
    updateNumericInput(session, "BL_sysbp", value = 125)
    updateNumericInput(session, "BL_diabp", value = 85)
    updateNumericInput(session, "BL_hemo", value = 15)
    updateNumericInput(session, "BL_uacr", value = 10)
    updateNumericInput(session, "BL_eGFR", value = 70)

    updateCheckboxInput(session, "BL_med_bp", value=F)
    updateCheckboxInput(session, "BL_med_dm", value=T)
    updateCheckboxInput(session, "BL_med_lipid", value=F)
    
    
  })
})
  

