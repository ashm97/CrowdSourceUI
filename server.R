##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Main Server Function
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

# --- Main Server Function
server <- function(input, output, session) {
  options(shiny.maxRequestSize=200*1024^2) # Increase max upload six to 200mb
  
  ## Calling modules
  callModule(homePage,"mainPage")
  current_dataSet_server_side <- callModule(dataPage,"datP",current_dataSet_server_side)#data input page containing 
  callModule(singleScatPage,"scat",current_dataSet_server_side)# single scatter display page
  callModule(hist, "hist", current_dataSet_server_side)#hist page
  callModule(cleavPage, "cleav", current_dataSet_server_side)#cleav page
  callModule(ptmPage,"ptmPage",current_dataSet_server_side)#PTM page
  callModule(decoyDis, "dec", current_dataSet_server_side)#decoy display page
  callModule(scorePageDisplay,"scoreDis",current_dataSet_server_side)# score display page
  
}