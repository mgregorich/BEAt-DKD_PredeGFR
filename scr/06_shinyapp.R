# Copy necessary data
file.copy(file.path(out.path, "riskpred_model.rds"), 
          file.path("./shiny", "riskpred_model.rds"), 
          overwrite = TRUE)

# Run app in Rstudio from correct directory
rm(list=ls())
runApp("./shiny")
