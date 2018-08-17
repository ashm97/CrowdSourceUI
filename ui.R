##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Main UI Function
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

library(ggplot2)
library(shinydashboard)
library(shiny)
library(plotly)
library(DT)
library(stringr)
library(mzR)
library(dplyr)
library(shinyalert)
library(leaflet)

source('support/support_functions.R')
source('support/plots.R')
source('support/modules.R')
source('support/body.R')



ui <- dashboardPage(skin = "black",
  dashboardHeader(title = "Results Analyser",disable = F),
  dashboardSidebar(
    useShinyalert(), 
    sidebarMenu(
      menuItem("Home", tabName = "main", icon = icon("home")),
      menuItem("Data", tabName = "data", icon = icon("table")),
      
      menuItem("Analysis", tabName = "analysis", icon = icon("bolt"),
               #menuItem("Scatters", tabName = "scatters"),
               menuItem("Scatter", tabName = "scatter"),
               menuItem("Histograms", tabName = "histograms"),
               menuItem("Cleavages", tabName = "cleavages"),
               menuItem("PTM", tabName = "ptm"),
               menuItem("Statistical", tabName = "statistical",icon=icon("line-chart"),
                        menuItem("Decoy", tabName = "decoy"),
                        menuItem("Score", tabName = "score"))
      ),
      #menuItem("Protein Level", tabName = "protein"),
      hr(),
      menuItem("Information", tabName = "info",icon=icon("info")),
      hr()#,
      #sidebarUserPanel(name = a("Prof Andy Jones", target = "_blank_",
      #                          href = "https://www.liverpool.ac.uk/integrative-biology/staff/andrew-jones/"), 
      #                 subtitle = "Professor of Bioinformatics",
      #                 image = "data:image/png;base64,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"),
      #sidebarUserPanel(name = a("Dr Andrew Collins", target = "_blank_",
      #                          href = "https://www.liverpool.ac.uk/integrative-biology/staff/a-collins/"), 
      #                 subtitle = "Research Associate",
      #                 image = "data:image/png;base64,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")
    )
  ),
  body
  
)

