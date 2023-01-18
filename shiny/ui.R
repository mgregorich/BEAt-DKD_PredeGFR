#################################
# Author:MG
# Date: 07/10/2021
# Info: user interface file - shiny
################################



# -------------------- SHINY User Interface  ---------------------------
corner_element = HTML(paste0('<a href=\"https://www.beat-dkd.eu/\">BEAt-DKD</a>'))

navbarPage(corner_element, 
           id = "tabselected",
           windowTitle = "BEAt-DKD | Prediction model",
           theme = shinythemes::shinytheme("spacelab"),
           tabPanel(title="Home",
                    fluidPage(
                      h1(id="big-heading" ,"Web implementation of the prediction model for kidney function in individuals with type 2 diabetes mellitus"),
                      
                      h4(id="text1", "Chronic kidney disease (CKD) is a prevalent issue in people with type 2 diabetes mellitus. In practice, the kidney function is essentially monitored by 
                         measuring the patient's estimated globular filtration rate (eGFR) over time. Diagnostic tools may facilitate continuous monitoring of disease progression and allow to identify people with 
                         an increased risk for kidney function decline over the next few years."),
                      h4(id="text1", "The prediction model for the future eGFR trajectory has been developed for people:"),
                      h4(id="text1", HTML("<ul><li>with diabetes mellitus type 2,</li><li>aged between 18 to 75 years, and</li><li>with a current eGFR above 30 mL/min/1.73m². </li></ul>")),
                      h4(id="text1", "    "),
                      h4(id="text1", "Users of this tool are urged to read the disclaimer carefully."),
                      h4(id="text1", "    "),
                      
                      
                      h4(id="text4","Disclaimer:"),
                      h4(id="text3", "Users of the prediction tool should not rely on information provided by the prediction tool for their own health problems. Questions should be addressed to the treating physician or other healthcare providers."),
                      h4(id="text3", "The authors make no warranties, nor expressed or implied representations whatsoever, regarding the accuracy, completeness, timeliness, comparative or controversial nature, or usefulness of any information contained or referenced in the prediction tool. The authors do not assume any risk whatsoever for the use of the prediction tool or the information contained herein. Health related information changes frequently and therefore information contained in the prediction tool may be outdated, incomplete or incorrect."),
                      h4(id="text3", "Use of the prediction tool does not create an expressed or implied physician-patient relationship. The authors do not endorse or claim validity for the prediction tool found on this webpage. This web tool does not record specific prediction tool user information and does not contact users of the prediction tools."),
                      h4(id="text3", "You are hereby advised to consult with a physician or other professional healthcare provider prior to making any decisions, or undertaking any actions or not undertaking any actions related to any healthcare problem or issue you might have at any time, now or in the future. In using the prediction tool, you agree that neither the authors nor any other party is or will be liable or otherwise responsible for any decision made or any action taken or any action not taken due to your use of any information presented in the prediction tool.")
                    ), 
           ),
           tabPanel(title="Model",
                    
                    fluidPage( 
                      tags$head(
                        tags$style(
                          HTML(".checkbox {margin: 0}
                                .checkbox p {margin: 0;}
                                .shiny-input-container {margin-bottom: 0;}" ),
                          HTML("hr {border-top: 1px solid #000000;}"),
                          HTML("#big-heading{color: #4181BB;
                               padding-bottom: 20px;}"),
                          HTML("#text1{color: black;
                               font-size: 16px;
                               line-height: 1.5em;
                               border: 2px double white;
                               background-color: white;
                               padding-right: 2.5%;
                               padding-left:2.5%
                               }"),
                          HTML("#text10{color: #1476CB;
                               font-size: 16px;
                               line-height: 1.5em;
                               font-weight: bold;
                               border: 2px double white;
                               background-color: white;
                               padding-right: 2.5%;
                               padding-left:2.5%
                               }"),
                          HTML("#text11{color: black;
                               font-size: 14px;
                               line-height: 1.5em;
                               border: 2px double white;
                               background-color: white;
                               padding-right: 2.5%;
                               padding-left:2.5%
                               }"),
                          HTML("#text2{color: #1476CB;
                               font-weight: bold;
                               font-size: 18px;
                               line-height: 1.5em;
                               margin-top: 5px;
                               white-space: pre-wrap;
                               }"),
                          HTML("#text5{color: #1476CB;
                               font-weight: bold;
                               font-size: 18px;
                               line-height: 1.5em;
                               margin-top: 5px;
                               padding-top: 7.5%;
                               white-space: pre-wrap;
                               }"),
                          HTML("#text0{color: #1476CB;
                               font-weight: bold;
                               font-size: 18px;
                               line-height: 1.5em;
                               margin-top: 5px;
                               white-space: pre-wrap;
                               }"),
                          HTML("#text3{color: black;
                               font-size: 14px;
                               line-height: 1.5em;
                               border: 2px double white;
                               background-color: white;
                                margin-top: 5px;
                               padding-right: 2.5%;
                               padding-left:2.5%
                               }"),
                          HTML("#text4{color: black;
                               font-weight: bold;
                               font-size: 14px;
                               line-height: 1.5em;
                               margin-top: 5px;
                               }"),
                          HTML("#basic1 {
                            border: 2px double #E5F1FC;
                            background-color: #E5F1FC;}"),
                          HTML("#lab1 {
                            border: 2px double #E5F1FC;
                            background-color: #E5F1FC;}"),
                          HTML("#border1 {
                          border: 2px solid #1165AE;
                          padding-top:2.55%;
                          padding-bottom:2.5%},
                          #border2 {
                          padding-top:15%;}")
                        )
                      ), 
                      
                      # App title ----
                      h1(id="big-heading" ,"Prediction model for kidney function in individuals with type 2 diabetes mellitus"),
                      # Sidebar layout with input and output definitions ----
                      
                      fluidRow(shinyjs::useShinyjs(),
                               tags$head(
                                 tags$style(type="text/css", 
                                 ".form-control.shiny-bound-input, 
                                 .selectize-input {height: 25px; margin-bottom:2 px; margin-top: 2px;}
                                 .checkbox-inline, .radio-inline {margin-bottom: 2px; margin-top: 2px;}")
                               ),
                               column(4,id="basic1",
                                      fluidRow(id="border1",align="center",
                                               h4(id="text0","Insert your information below.")),
                 
                                      h4(""),
                                      fluidRow(id="border1",
                                        column(6,id="basic1", 
                                               h4(id="text2","Demographics:"),
                                               # Input: Slider for the number of bins ----
                                               h5(shinyWidgets::prettyRadioButtons("BL_sex", "Sex:",
                                                                     c("female"=1, "male"=0), inline=T, shape = "square")),
                                               h5(shinyWidgets::numericInputIcon("BL_age", "Age (years):", value=65, min=17, max=76)), # min & max are a bit larger e.g. [18,75] but need to be wider here for shinyvalidate
                                               h5(shinyWidgets::numericInputIcon("BL_bmi", "BMI (kg/m²):", value=25, min=9, max=41)),
                                               h5(shinyWidgets::prettyRadioButtons("BL_smoking", "Smoking status:",
                                                                     choices = c("never"=0, "ever"=1),inline=T,shape = "square")),
                                               h4(""),
                                               h4(id="text5","Medication intake:"),
                                               shinyWidgets::prettyCheckbox("BL_med_bp", "Blood pressure-lowering", status = "primary"),
                                               shinyWidgets::prettyCheckbox("BL_med_lipid", "Lipid-lowering", status = "primary"),
                                               shinyWidgets::prettyCheckbox("BL_med_dm", "Glucose-lowering", status = "primary")),
                                        column(6,id="basic1",
                                                 fluidRow(id="lab1",
                                                          h4(id="text2", "Laboratory:"),
                                                          h5(div(style="margin:5px;",shinyWidgets::numericInputIcon("BL_hemo", "Hemoglobin (g/dL):",value=15, min=9, max=19))), 
                                                          h5(div(style="margin:5px;",shinyWidgets::numericInputIcon("BL_diabp", "Diastolic blood pressure (mmHg):",value=85, min=49, max=106))),
                                                          h5(div(style="margin:5px;",shinyWidgets::numericInputIcon("BL_sysbp", "Systolic blood pressure (mmHg):", value=125, min=99, max=191))),
                                                          h5(div(style="margin:5px;",shinyWidgets::numericInputIcon("BL_hba1c", "Hba1C (mmol/mol):", value=62, min=32, max=94))),
                                                          h5(div(style="margin:5px;",shinyWidgets::numericInputIcon("BL_serumchol", "Serum cholesterol (mg/dL):", value=130, min=98, max=329))),
                                                          h5(div(style="margin:5px;",shinyWidgets::numericInputIcon("BL_uacr", "Urinary albumin-creatinine ratio (mg/g):", value=10, min=0.04, max=2550)))
                                                 ))),
                                      h4(""),
                                      fluidRow(align="center", id="border1",
                                        column(5,h5(div(style="color: #1165AE;margin-top:-15px;width=400px", shinyWidgets::numericInputIcon("BL_eGFR", "Baseline eGFR (mL/min/1.73m²):",value=70, min=29, max=146, width="50%")))),
                                        column(6,h5(div(style="color: #1165AE;margin-top:-15px;width=400px", shinyWidgets::numericInputIcon("cutpoint", "eGFR cutpoint to classify stable and fast progression:", value=-3, min=-6, max=1)))),
                                        
                                      ), 
                                      fluidRow(
                                        h2(""),
                                        column(12,)
                                      ),
                                      fluidRow( 
                                        column(6, align="right", actionButton("goButton", "Compute")),
                                        column(6, align="left", actionButton("reset_input", "Reset"))
                                      ),
                                      fluidRow(
                                        h2(""),
                                        column(12,)
                                      ),
                                      ),
                               # Main panel for displaying outputs ----
                               column(8,
                                      h4(id="text0", "Results of the prediction model"),
                                      
                                      # Output: Histogram ----
                                      tabsetPanel(
                                        tabPanel("Longitudinal",
                                                 h4(id="text1", textOutput("text_longitudinal1")),
                                                 h4(id="text1", textOutput("text_longitudinal2")),
                                                 fluidRow(align="center",
                                                          plotOutput("plot_trajectory", height="500px", width="700px"))),
                                        tabPanel("Risk", 
                                                 h4(id="text1", textOutput("text_risk1")),
                                                 h4(id="text1", textOutput("text_risk2")),
                                                 h4(id="text1", textOutput("text_risk3")),
                                                 fluidRow(align="center",
                                                          plotOutput("plot_slopedist", height="400px", width="700px")
                                                          ),
                                                 h4(id="text1", textOutput("text_risk4")),
                                        ),
                                        tabPanel("Data", 
                                                 fluidRow(
                                          h4(id="text1", htmlOutput("text_data")),
                                          column(6, 
                                                 h4(id="text2", "Provided patient information"),
                                                 fluidRow(
                                                   column(6, h4(id="text3",tableOutput("table_new1"))),
                                                   column(6, h4(id="text3",tableOutput("table_new2")))),
                                                 h4(id="text3",tableOutput("table_new3"))),
                                          column(6,
                                                 h4(id="text2", "Predicted eGFR per time point"),
                                                 h4(id="text3",tableOutput("table_pred1")),
                                                 h4(id="text3",tableOutput("table_pred2")) )
                                          ))
                                        
                                      )
                               ),
                      )
                    )
           ),
           tabPanel(title="About",
                    fluidPage(
                      h3(id="text2", "Details"),
                      h4(id="text1", "The web tool is the online implementation of the prediction model for kidney function in individuals with type 2 diabetes mellitus presented in (1) and (2)."),
                      h3(id="text1", "   "),
                      h3(id="text2", "Contact"),
                      h4(id="text10", "For technical inquiries:"),
                      
                      h4(id="text1", "Mariella Gregorich (Section for Clinical Biometrics, Center for Medical Statistics, Informatics and Intelligent Systems, Medical University of Vienna)"),
                      h4(id="text1", "Email:" ,a("mariella.gregorich@meduniwien.ac.at", href="mailto:mariella.gregorich@meduniwien.ac.at")),
                      h3(id="text1", "   "),
                      h4(id="text10", "For content-related inquiries:"),
                      h4(id="text1", "Rainer Oberbauer (Division of Nephrology and Dialysis, Department of Medicine III, Medical University of Vienna, Vienna, Austria)"),
                      h4(id="text1", "Email:" ,a("rainer.oberbauer@meduniwien.ac.at", href="mailto:rainer.oberbauer@meduniwien.ac.at")),
                      h3(id="text1", "   "),
                  
                      h3(id="text2", "Funding"),
                      h4(id="text1", "The research leading to these results has received support from the Innovative Medicines Initiative Undertaking under grant agreement no. 115974 BEAt-DKD. "),
                      h3(id="text1", "   "),
                      h3(id="text2", "References"),
                      
                      h4(id="text1", "(1) Gregorich, M., Heinzel, A., Kammer, M., Meiselbach, H., Böger, C., Eckardt, K. U., Mayer, G., Heinze, G. & Oberbauer, R. (2022). Individual-specific prediction of future eGFR in people with type 2 diabetes mellitus: development and external validation (in submission)"),
                      h4(id="text1", "(2) Gregorich, M., Heinzel, A., Kammer, M., Meiselbach, H., Böger, C., Eckardt, K. U., Mayer, G., Heinze, G. & Oberbauer, R. (2021). A prediction model for the decline in renal function in people with type 2 diabetes mellitus: study protocol. Diagnostic and Prognostic Research, 5(1), 1-9.")
                      


 )
           )
)