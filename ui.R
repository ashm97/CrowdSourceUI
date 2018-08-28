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
library(splitstackshape)
library(readr)

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
      hr()
    )
  ),
  body
  
)

