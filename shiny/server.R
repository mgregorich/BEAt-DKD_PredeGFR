#################################
# Author:MG
# Date: 07/10/2021
# Info: server file - shiny
################################


# ---------------------------- SHINY Server --------------------------
source("www/functions_shiny.R",local=TRUE)

# load model
pred_model <- readJSON("www/predmodel_shinyObject.json")
# input conversion
pred_model$form = as.formula(pred_model$form)
names(pred_model$betas) = pred_model$betas_names

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
    slope_cutpoint <- input$cutpoint
    
    res <- update_PredByBase(pred_model, 
                             newdata = data.new,  
                             cutpoint = slope_cutpoint,
                             times = seq(1,5,1))
    data.long <- data.new[rep(1,6),]
    data.long$Time <- seq(0,5,1)
    data.pred <- dplyr::full_join(data.long,
                           res[,c("PatID", "Time","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")],
                           by=c("PatID", "Time"))
    #data.pred[data.pred$Time==0,]$pred <- input$BL_eGFR
    data.pred <- data.frame(data.pred[,c("PatID", "Time","pred", "pred.lo", "pred.up", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")])
    data.pred
      })
  
  makeOutputTable <- eventReactive(input$goButton, {
    odata <- calc_out()
    odata <- data.frame(odata[,-1])
    
    odata$pred.CI <- paste0("(",round(odata$pred.lo,2),", ",round(odata$pred.up,2), ")")
    odata$slope.CI <- paste0("(",round(odata$pred.slope.lo,2),", ",round(odata$pred.slope.up,2), ")")
    
    pred <- data.frame("Time"=odata$Time, "Prediction"=odata$pred, "CI"=odata$pred.CI)
    pred$CI[1] <- ""
    colnames(pred) <- c("Year", "Prediction","95% Confidence interval")
   
    slope <- data.frame("Slope"=odata$pred.slope[1], "CI"=odata$slope.CI[1])
    colnames(slope) <- c("eGFR slope", "95% Confidence interval")

     return(list("pred"=pred, "slope"=slope, "prob"=odata$pred.prob[1]))
  })
  
  makeInputTable <- eventReactive(input$goButton, {
    idata <- inputdata()
    idata <- data.frame(idata)[,-c(1,2,3)]
    idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")] <- apply(idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")],2, function(x) ifelse(x==1, "yes", "no"))
    idata$BL_smoking <- ifelse(idata$BL_smoking==1, "ever", "never")
    idata$BL_sex <- ifelse(idata$BL_sex==1, "female", "male")
    
    if(input$add_pred==2){
      idemo <- data.frame("General"=c("Baseline eGFR", "Sex", "Age", "BMI", "Smoking"), 
                          "Values"=unlist(idata[,c("FU_eGFR_epi", "BL_sex", "BL_age", "BL_bmi", "BL_smoking")]))
      imed <- data.frame("Medication"=c("Glucose-lowering", "Blood pressure-lowering", "Lipid-lowering"), 
                         "Values"=unlist(idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")]))
      ilab <- data.frame("Laboratory"=c("Hemoglobin", "Mean arterial pressure",  "Hba1C", "Serum cholesterol", "log-2 UACR"), 
                         "Values"=unlist(idata[,c("BL_hemo", "BL_map", "BL_hba1c", "BL_serumchol", "BL_uacr_log2")]))

      return(list("idemo"=idemo,"ilab"=ilab,"imed"=imed))
    }
    else{
      idemo <- data.frame("General"=c("Baseline eGFR", "Sex", "Age", "BMI", "Smoking"), 
                          "Values"=unlist(idata[,c("FU_eGFR_epi","BL_sex", "BL_age", "BL_bmi", "BL_smoking")]))
      imed <- data.frame("Medication"=c("Glucose-lowering", "Blood pressure-lowering", "Lipid-lowering"), 
                         "Values"=unlist(idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")]))
      
      return(list("idemo"=idemo,"imed"=imed))
    }
    
  })
  
  # ---------- Execute --------------------------------------------------------
  text_risk1 <- eventReactive(input$goButton, {  
    paste("Based on the provided information and a slope cutpoint of ", input$cutpoint," mL/min/1.73m², 
          the patient's probability of rapid decline in renal function is ",  
          round(calc_out()[1,"pred.prob"]*100,2), "%.", sep="" )
    })
  
  text_risk2 <- eventReactive(input$goButton, {  
    paste("This means that in 100 people with these characteristics ",round(calc_out()[1,"pred.prob"]*100,0)," are expected to exhibit fast progression of renal decline in the following years.", sep="" )
  })
  
  text_longitudinal <- eventReactive(input$goButton, {  
    paste("The figure below illustrates the expected longitudinal trajectory of the patient's future eGFR measurements given the provided information.", sep="" )
  })
  
  plotsmile <- eventReactive(input$goButton, {smilegraph(round(calc_out()[1,"pred.prob"]*100,2))})
  
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
    updateNumericInput(session, "BL_age", value = 65)
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
  

