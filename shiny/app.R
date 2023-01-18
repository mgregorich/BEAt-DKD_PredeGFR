#################################
# Author:MG
# Date: 07/10/2021
# Info: app file - shiny
################################


# Create Shiny app ----
shinyApp(ui = ui, server = server)

# Upload to shinyapps.io
library(rsconnect)
# rsconnect::setAccountInfo(name='beatdkd',
#                           token='<TOKEN>',  # token and secret must be inserted
#                           secret='<SECRET>')
deployApp(appDir = here::here("shiny/"))
