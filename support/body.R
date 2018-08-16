##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Contains the body for the UI function
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

source('support/tabs.R')

body <- dashboardBody(
  tabItems(
    # Tab list
    data_tab,
    scatters_tab,
    histograms_tab,
    cleavages_tab,
    ptm_tab,
    decoy_tab,
    score_tab,
    #protein_tab
    info_tab
    
  )
)
