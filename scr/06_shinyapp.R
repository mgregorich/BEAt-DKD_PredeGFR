# Copy necessary data
file.copy(file.path(out.path, "riskpred_model.rds"), 
          file.path(shiny.path, "riskpred_model.rds"))

# Run app in Rstudio from correct directory
runApp(shiny.path)
