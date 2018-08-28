##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Contains the body for the UI function
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

source('support/tabs.R')

body <- dashboardBody(

  tags$style(HTML("

                  .box.box-solid.box-primary>.box-header {
                  
                  background:#D0D0D0
                  }
                  
                  .box.box-solid.box-primary{
                  border-bottom-color:#D0D0D0;
                  border-left-color:#D0D0D0;
                  border-right-color:#D0D0D0;
                  border-top-color:#D0D0D0;
                  background:navy
                  }
                  
                  ")),
  
  tags$style(HTML("

                  .box.box-solid.box-info>.box-header {
                  
                  background:#D0D0D0
                  }
                  
                  .box.box-solid.box-info{
                  border-bottom-color:#D0D0D0;
                  border-left-color:#D0D0D0;
                  border-right-color:#D0D0D0;
                  border-top-color:#D0D0D0
                  
                  }
                  
                  ")),
  tabItems(
    # Tab list
    main_tab,
    data_tab,
    scatter_tab,
    histograms_tab,
    cleavages_tab,
    ptm_tab,
    decoy_tab,
    score_tab,
    #protein_tab
    info_tab
    
  )
)
