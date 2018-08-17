##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Main application launch
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

source('server.R') 
source('ui.R',local = TRUE)

shinyApp(
  ui = ui,
  server = server
)