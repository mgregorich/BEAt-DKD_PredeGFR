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
rsconnect::setAccountInfo(name='beatdkd',
                          token='5250922C5FF6DBA42E91AC641FE261BA',
                          secret='aTNlStV0f/F1UmmS1LuZ3N6nWVeWg8n9uTxhn0Ar')
deployApp(appDir = here::here("shiny/"))
