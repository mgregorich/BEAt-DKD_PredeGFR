#################################
# Author:MG
# Date: 07/10/2021
# Info: user interface file - shiny
################################



# -------------------- SHINY User Interface  ---------------------------

navbarPage("BEAt-DKD", theme = shinythemes::shinytheme("spacelab"),
           tabPanel(title="Home",
                    
                    fluidPage(
                      h1(id="big-heading" ,"Web implementation of the prediction model for renal decline in individuals with type 2 diabetes mellitus"),
                      
                      h4(id="text1", "Renal decline is a prevalent issue in patients with type 2 diabetes mellitus. Early on detection through the use of diagnostic 
                         tools facilitating continuous monitoring which allow to identify subjects with an increased risk for kidney decline over the course 
                         of the next years."),
                      h4(id="text1", "The prediction model has been developed for people with diabetes mellitus type 2, in between 18 to 75 years and an eGFR above 30mL/min/1.73m2."),
                      h4(id="text1", "Users are urged to read the disclaimer carefully."),
                      
                      h4(id="text4","Disclaimer:"),
                      h4(id="text3", "Users of the prediction tool should not rely on information provided by the prediction tool for their own health problems. Questions should be addressed to your own physician or other healthcare provider. Moreover, the prediction tool does not predict, whether a single patient will suffer from hypothermia during surgery or not, because this is influenced by a large variety of genetic and acquired factors, some of which are still unknown."),
                      h4(id="text3", "CeMSIIS makes no warranties, nor express or implied representations whatsoever, regarding the accuracy, completeness, timeliness, comparative or controversial nature, or usefulness of any information contained or referenced in the prediction tool. CeMSIIS does not assume any risk whatsoever for your use of the prediction tool or the information contained herein. Health related information changes frequently and therefore information contained in the prediction tool may be outdated, incomplete or incorrect."),
                      h4(id="text3", "Use of the prediction tool does not create an express or implied physician-patient relationship. CeMSIIS does not endorse or claim validity for the prediction tools found on the meduniwien.ac.at website. The activities and products of CeMSIIS and its developers and agents are not endorsed by our past, present, or future employers. CeMSIIS does not record specific prediction tool user information and does not contact users of the prediction tools."),
                      h4(id="text3", "You are hereby advised to consult with a physician or other professional healthcare provider prior to making any decisions, or undertaking any actions or not undertaking any actions related to any healthcare problem or issue you might have at any time, now or in the future. In using the prediction tools, you agree that neither CeMSIIS nor any other party is or will be liable or otherwise responsible for any decision made or any action taken or any action not taken due to your use of any information presented in the prediction tool.")
                    )
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
                          HTML("#text2{color: #1476CB;
                               font-weight: bold;
                               font-size: 18px;
                               line-height: 1.5em;
                               margin-top: 5px;
                               white-space: pre-wrap;
                               }"),
                          HTML("#text0{color: #1165AE;
                               font-weight: bold;
                               font-size: 18px;
                               line-height: 1.5em;
                               margin-top: 5px;
                               white-space: pre-wrap;
                               }"),
                          HTML("#text3{color: black;
                               font-size: 12px;
                               line-height: 1.5em;
                               border: 2px double white;
                               background-color: white;
                                margin-top: 5px;
                               padding-right: 2.5%;
                               padding-left:2.5%
                               }"),
                          HTML("#text4{color: black;
                               font-weight: bold;
                               font-size: 12px;
                               line-height: 1.5em;
                               margin-top: 5px;
                               }"),
                          HTML("#basic1 {
                            border: 2px double #E5F1FC;
                            background-color: #E5F1FC;}"),
                          HTML("#lab1 {
                            border: 2px double #E5F1FC;
                            background-color: #E5F1FC;}")
                        )
                      ), 
                      
                      # App title ----
                      h1(id="big-heading" ,"Web implementation of the prediction model for renal decline in individuals with type 2 diabetes mellitus"),
                      # Sidebar layout with input and output definitions ----
                      
                      fluidRow(shinyjs::useShinyjs(),
                               tags$head(
                                 tags$style(type="text/css", 
                                 ".form-control.shiny-bound-input, 
                                 .selectize-input {height: 25px; margin-bottom:2 px; margin-top: 2px;}
                                 .checkbox-inline, .radio-inline {margin-bottom: 2px; margin-top: 2px;}")
                               ),
                               column(4,id="basic1",
                                      fluidRow(
                                               column(4,h4(id="text0","Choose Model: ")),
                                               column(5,h4(radioButtons("add_pred", "", 
                                                                     choices = c("Simple" = 1, 
                                                                                 "Extended" = 2),
                                                                     selected = 1, inline = T)))),
                                      h4(" ------------------------------------------------------------------------------------------------------ "),
                                      fluidRow(
                                        column(6,id="basic1", 
                                               h4(id="text2","Demographics:"),
                                               # Input: Slider for the number of bins ----
                                               h5(shinyWidgets::prettyRadioButtons("BL_sex", "Sex:",
                                                                     c("female"=1, "male"=0), inline=T, shape = "square")),
                                               h5(shinyWidgets::numericInputIcon("BL_age", "Age, years:", value=65, min=18, max=75)),
                                               h5(shinyWidgets::numericInputIcon("BL_bmi", "BMI (kg/m2):", value=25, min=11, max=40)),
                                               h5(shinyWidgets::prettyRadioButtons("BL_smoking", "Smoking status:",
                                                                     choices = c("never"=0, "ever"=1),inline=T,shape = "square")),
                                               
                                               h4(" ------------------------------------ "),
                                               h4(id="text2","Medication intake:"),
                                               shinyWidgets::prettyCheckbox("BL_med_bp", "Blood pressure-lowering"),
                                               shinyWidgets::prettyCheckbox("BL_med_lipid", "Lipid-lowering"),
                                               shinyWidgets::prettyCheckbox("BL_med_dm", "Glucose-lowering")),
                                        column(6,id="basic1", shinyjs::hidden(
                                                 fluidRow(id="lab1",
                                                          h4(id="text2", "Laboratory:"),
                                                          h5(div(style="margin:5px;",numericInput("BL_hemo", "Hemoglobin (g/dL):",value=15, min=10, max=18))), 
                                                          h5(div(style="margin:5px;",numericInput("BL_diabp", "Diastolic blood pressure:",value=85, min=50, max=105))),
                                                          h5(div(style="margin:5px;",numericInput("BL_sysbp", "Systolic blood pressure:", value=125, min=100, max=190))),
                                                          h5(div(style="margin:5px;",numericInput("BL_hba1c", "Hba1C (mmol/mol):", value=62, min=33.3, max=93))),
                                                          h5(div(style="margin:5px;",numericInput("BL_serumchol", "Serum Cholesterol (mg/dL)", value=130, min=99, max=328))),
                                                          h5(div(style="margin:5px;",numericInput("BL_uacr", "Urinary albumin-creatinine ratio (mg/g):", value=10, min=0.05, max=2549)))
                                                 )))),
                                      h4(" ----------------------------------------------------------------------------------------------------- "),
                                      fluidRow(align="center",
                                        column(5,h5(div(style="color: #1165AE;margin-top:-15px;width=400px", numericInput("BL_eGFR", "Baseline eGFR",value=70, min=30, max=145, width="50%")))),
                                        column(6,h5(div(style="color: #1165AE;margin-top:-15px;width=400px", numericInput("cutpoint", "Cutpoint for stable and fast progression", value=-3, min=-5, max=0)))),
                                        
                                      )
                                      ),
                               # Main panel for displaying outputs ----
                               column(8,
                                      h4(id="text0", "Results of the prediction model"),
                                      
                                      # Output: Histogram ----
                                      tabsetPanel(
                                        tabPanel("Risk", 
                                                 h4(id="text1", textOutput("text_risk1")),
                                                 h4(id="text1", textOutput("text_risk2")),
                                                 plotOutput("plot_risk", height="400px", width="400px")),
                                        tabPanel("Longitudinal",
                                                 h4(id="text1", textOutput("text_longitudinal")),
                                                 plotOutput("plot_trajectory", height="500px", width="700px")),
                                        tabPanel("Data", 
                                                 fluidRow(
                                          h4(id="text1", htmlOutput("text_data")),
                                          column(6, 
                                                 h4(id="text2", "Provided patient information:"),
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
                      ),
                      fluidRow(
                        column(4,align="center",
                               h5(div(style="margin-bottom:-5px;width=400px", numericInput("cutpoint", "Cutpoint for stable and fast progression", value=-3, min=-5, max=0)))),
                        column(8,)),
                      fluidRow(
                        column(2, align="right", actionButton("goButton", "Compute")),
                        column(2, align="left", actionButton("reset_input", "Reset")),
                        column(8,)),
                      fluidRow(
                        hr(),
                        column(12,)
                      )
                    )
           )
)