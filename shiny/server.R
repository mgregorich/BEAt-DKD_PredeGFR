#################################
# Author:MG
# Date: 07/10/2021
# Info: server file - shiny
################################



# ---------------------------- SHINY Server --------------------------
source("functions_shiny.R",local=TRUE)

# load model
pred_model <- readJSON("www/predmodel_shinyObject.json")
# input conversion
pred_model$form = as.formula(pred_model$form)
names(pred_model$betas) = pred_model$betas_names



shinyServer(function(input, output, session) {
  iv <- shinyvalidate::InputValidator$new()
  iv$add_rule("BL_age", shinyvalidate::sv_required())
  iv$add_rule("BL_age", shinyvalidate::sv_between(left=18, right=75))
  iv$add_rule("BL_bmi", shinyvalidate::sv_required())
  iv$add_rule("BL_bmi", shinyvalidate::sv_between(left=10, right=40))
  iv$add_rule("BL_hemo", shinyvalidate::sv_required())
  iv$add_rule("BL_hemo", shinyvalidate::sv_between(left=10, right=18))
  iv$add_rule("BL_diabp", shinyvalidate::sv_required())
  iv$add_rule("BL_diabp", shinyvalidate::sv_between(left=50, right=105))
  iv$add_rule("BL_sysbp", shinyvalidate::sv_required())
  iv$add_rule("BL_sysbp", shinyvalidate::sv_between(left=100, right=190))
  iv$add_rule("BL_hba1c", shinyvalidate::sv_required())
  iv$add_rule("BL_hba1c", shinyvalidate::sv_between(left=33, right=93))
  iv$add_rule("BL_serumchol", shinyvalidate::sv_required())
  iv$add_rule("BL_serumchol", shinyvalidate::sv_between(left=99, right=328))
  iv$add_rule("BL_uacr", shinyvalidate::sv_required())
  iv$add_rule("BL_uacr", shinyvalidate::sv_between(left=0.05, right=2549))
  iv$add_rule("BL_eGFR", shinyvalidate::sv_required())
  iv$add_rule("BL_eGFR", shinyvalidate::sv_between(left=30, right=145))
  iv$add_rule("cutpoint", shinyvalidate::sv_required())
  iv$add_rule("cutpoint", shinyvalidate::sv_between(left=-5, right=0))
  
  ## ==== Display message at start =====
  query_modal <- modalDialog(
    title = 'Note',
    HTML(
      '
        <p>
        This website is still in the development phase. In case of crashes or unexpected behavior please report to the main developer
        <a href="mailto:mariella.gregorich@meduniwien.ac.at">mariella.gregorich@meduniwien.ac.at</a>. 
        <br>
        Thank you for your understanding!
        </p>
        '
    ),
    easyClose = F, 
    footer = modalButton("Cancel")
  )
  showModal(query_modal)
  
  # ===== Display message when switching to model tab =====
  observeEvent(input$tabselected, {
    if (input$tabselected == "Model") 
      shinyalert::shinyalert("Info", "The presented results are only for informative purposes and do not replace medical consultation!", type = "warning")  
  })
  
  # ===== Check input values ======
  observeEvent(input$goButton, {
    if (!iv$is_valid()) {
      iv$enable() # Start showing validation feedback
      shinyalert::shinyalert("Oops!", "Please correct the errors indicated in red in the form and try again.", type = "error")
    }
  })
  
  ## ====== INPUT ================
  inputdata <- eventReactive(input$goButton ,{
    
  data.new <- data.frame(PatID="Pnew", Time=0, Country="Unknown", FU_eGFR_epi_2021=input$BL_eGFR,
                         BL_age=input$BL_age, 
                         BL_sex=as.numeric(input$BL_sex), 
                         BL_bmi=input$BL_bmi, 
                         BL_smoking=as.numeric(input$BL_smoking),
                         BL_hemo=input$BL_hemo, BL_hba1c=input$BL_hba1c, BL_serumchol=input$BL_serumchol, 
                         BL_map=calc_map(sys=input$BL_sysbp,dia=input$BL_diabp),
                         BL_uacr_log2=log(input$BL_uacr,2),
                         BL_med_dm=input$BL_med_dm*1, BL_med_bp=input$BL_med_bp*1, BL_med_lipid=input$BL_med_lipid*1)

    
    data.new
  })
  
  # ====== Compute predictions ================
  calc_out <- eventReactive(input$goButton,{
    req(iv$is_valid())
    
    data.new=inputdata()
    slope_cutpoint <- input$cutpoint
    
    res <- update_PredByBase(pred_model, 
                             newdata = data.new,  
                             cutpoint = slope_cutpoint)
    data.long <- data.new[rep(1,6),]
    data.long$Time <- seq(0,5,1)
    data.pred <- dplyr::full_join(data.long,
                                  res[,c("PatID", "Time","pred", "pred.lo.95", "pred.up.95", 
                                         "pred.lo.50", "pred.up.50","pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")],
                                  by=c("PatID", "Time"))
    #data.pred[data.pred$Time==0,]$pred <- input$BL_eGFR
    data.pred <- data.frame(data.pred[,c("PatID", "Time","pred", "pred.lo.95", "pred.up.95", 
                                         "pred.lo.50", "pred.up.50", "pred.slope", "pred.slope.lo", "pred.slope.up","pred.prob")])
    data.pred$BaseGFR <- c(data.new$FU_eGFR_epi_2021, rep(NA, nrow(data.pred)-1))
    data.pred[data.pred$Time==0,c("pred", "pred.lo.95", "pred.up.95", "pred.lo.50", "pred.up.50")] <-NA
    data.pred

      })
  
  # ==== Output table =====
  makeOutputTable <- eventReactive(input$goButton, {
    odata <- calc_out()
    odata <- data.frame(odata[,-1])
    
    odata$pred.CI <- paste0("(",round(odata$pred.lo.95,2),", ",round(odata$pred.up.95,2), ")")
    odata$slope.CI <- paste0("(",round(odata$pred.slope.lo,2),", ",round(odata$pred.slope.up,2), ")")
    
    pred <- data.frame("Time"=odata$Time, "Prediction"=odata$pred, "CI"=odata$pred.CI)[-1,]
    colnames(pred) <- c("Year", "Prediction","95% Prediction interval")
   
    slope <- data.frame("Slope"=odata$pred.slope[1], "CI"=odata$slope.CI[1])
    colnames(slope) <- c("eGFR slope", "95% Prediction interval")

     return(list("pred"=pred, "slope"=slope, "prob"=odata$pred.prob[1]))
  })
  
  # ===== Input table =======
  makeInputTable <- eventReactive(input$goButton, {
    idata <- inputdata()
    idata <- data.frame(idata)[,-c(1,2,3)]
    idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")] <- apply(idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")],2, function(x) ifelse(x==1, "yes", "no"))
    idata$BL_smoking <- ifelse(idata$BL_smoking==1, "ever", "never")
    idata$BL_sex <- ifelse(idata$BL_sex==1, "female", "male")
    
    idemo <- data.frame("General"=c("Baseline eGFR", "Sex", "Age", "BMI", "Smoking"), 
                        "Values"=unlist(idata[,c("FU_eGFR_epi_2021", "BL_sex", "BL_age", "BL_bmi", "BL_smoking")]))
    imed <- data.frame("Medication"=c("Glucose-lowering", "Blood pressure-lowering", "Lipid-lowering"), 
                       "Values"=unlist(idata[,c("BL_med_dm", "BL_med_bp", "BL_med_lipid")]))
    ilab <- data.frame("Laboratory"=c("Hemoglobin", "Mean arterial pressure",  "Hba1C", "Serum cholesterol", "log-2 UACR"), 
                       "Values"=unlist(idata[,c("BL_hemo", "BL_map", "BL_hba1c", "BL_serumchol", "BL_uacr_log2")]))

    return(list("idemo"=idemo,"ilab"=ilab,"imed"=imed))

    
  })
  
  plot_slopedist <- eventReactive(input$goButton, 
                                  {plot_eGFRslopeDistribution(distr.slopes=pred_model$distr.slopes, 
                                                              dens.slopes=pred_model$dens.slopes, 
                                                              pred.slope=calc_out()[1,"pred.slope"])})
  
  
  # ================= Text chunks ===================================
  text_risk1 <- eventReactive(input$goButton, {  
    paste("Based on the provided information and a slope cutpoint of ", input$cutpoint," mL/min/1.73m², 
          the patient's probability of rapid decline in kidney function is ",  
          round(calc_out()[1,"pred.prob"]*100,2), "%.", sep="" )
    })
  
  text_risk2 <- eventReactive(input$goButton, {  
    paste("This means that in 100 people with these characteristics ",round(calc_out()[1,"pred.prob"]*100,0)," are expected to exhibit fast progression of kidney function decline in the following years.", sep="" )
  })
  
  text_risk3 <- eventReactive(input$goButton, {
    paste("The graph below shows the distribution of the predicted eGFR curves in the development cohort with which the prediction model was developed. 
          The red line indicates the predicted eGFR slope associated with the input data and illustrates how the predicted eGFR slope of the new individual compares to the distribution seen in
          the development cohort.", sep="" )
  })
  
  text_risk4 <- eventReactive(input$goButton, {
    eGFR_prob = plot_slopedist()[["prob"]]
    if(eGFR_prob < 40){
      eGFR_pos <- "below average"
    }else if(eGFR_prob > 60){
      eGFR_pos <- "above average"
    }else{
      eGFR_pos <-"mid-range"
    }
    paste("The predicted eGFR curve has a slope of ",round(calc_out()[1,"pred.slope"],2)," (",round(calc_out()[1,"pred.slope.lo"],2),", ",round(calc_out()[1,"pred.slope.up"],2),"). This GFR slope is ",eGFR_pos," compared to the predicted eGFR slopes of the development cohort.
          ",round(eGFR_prob,1),"% of the included individuals in the development cohort were below the predicted eGFR curve of ",
          round(calc_out()[1,"pred.slope"],2),".", sep="" )
  })
  
  text_longitudinal <- eventReactive(input$goButton, {  
    paste("The figure below illustrates the expected longitudinal trajectory of the patient's future eGFR measurements and the corresponding 95% (light grey) and 50% prediction interval (dark grey) given the provided information. These intervals comprise both sampling variability as well as the variability of individual data points. The observed baseline eGFR of the patient is indicated in red.", sep="" )
  })
  
  
  text_data <- eventReactive(input$goButton,{
      HTML(paste0("The tables below contain both the patient information given and the predictions resulting from the model, including 95% prediction intervals for the individual time points. The patient's probability of progression to fast kidney function decline is also given."))
  })
  
  # ============= Output chunks ====================================
  output$text_risk1 <- renderText({
    text_risk1()
  })
  
  output$text_risk2 <- renderText({
    text_risk2()
  })
  
  output$text_risk3 <- renderText({
    text_risk3()
  })
  
  output$text_risk4 <- renderText({
    text_risk4()
  })
  
  output$text_data <- renderText({
    text_data()
  })
  
  output$text_longitudinal <- renderText({
    text_longitudinal()
  })
  
  output$plot_trajectory <- renderPlot({
    data.pred <- calc_out()
    plot_trajectory(data.pred)
  })
  
  output$plot_slopedist <- renderPlot({
    plot_slopedist()
  })
  
  output$table_new1 <- renderTable({
    makeInputTable()$idemo
  })
  
  output$table_new2 <- renderTable({
      makeInputTable()$ilab
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

  
  # =======  Reset =============
  observeEvent(input$reset_input, {
    updateNumericInput(session, "BL_age", value = 65)
    updateRadioButtons(session, "BL_sex", selected = 1)
    updateNumericInput(session, "BL_bmi", value = 25)
    updateRadioButtons(session, "BL_smoking", selected = 0)
    
    updateNumericInput(session, "BL_hba1c", value = 62)
    updateNumericInput(session, "BL_serumchol", value = 130)
    updateNumericInput(session, "BL_sysbp", value = 125)
    updateNumericInput(session, "BL_diabp", value = 85)
    updateNumericInput(session, "BL_hemo", value = 15)
    updateNumericInput(session, "BL_uacr", value = 10)
    updateNumericInput(session, "BL_eGFR", value = 70)
    updateNumericInput(session, "cutpoint", value = -3)

    updateCheckboxInput(session, "BL_med_bp", value = FALSE)
    updateCheckboxInput(session, "BL_med_dm", value = FALSE)
    updateCheckboxInput(session, "BL_med_lipid", value = FALSE)
  })
  session$onSessionEnded(stopApp) #automatically stop when closing browser
  
})
