##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Main UI Function
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

# --- Main User Interface Function
ui <- dashboardPage(skin = "black",
  dashboardHeader(title = "Results Analyser",disable = F), 
  dashboardSidebar(
    useShinyalert(), # For generating custom alert pop ups 
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
                        menuItem("Score", tabName = "score"))),
      hr(), # Line Break
      menuItem("Information", tabName = "info",icon=icon("info")),
      hr()
    )
  ),
  body # Callig the body for tabs 
)

