##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Page containing code for Shiny Modules
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

# -------------------------------------------------------------------

### Module for the Data input page

## UI Module

dataPageInput <- function(id) {
  ns <- NS(id)
  
  tagList(       #taglist for the UI output of this page
    fluidRow(
      column(width = 4,
             box(width = 12,
                 title = "File Upload",
                 # Input: Select a file ----
                 fileInput(ns("file1"), "Choose File",
                           multiple = FALSE,
                           accept = c(".csv",
                                      ".mzid",
                                      ".gz")),
                 
                 helpText("Upload a file. Accepted formats: .csv, .mzid, .gz")
             ),
             
             box(width = 12, title = "Scoring Column Selector",
                 
                 # Input: Select a column ----
                 uiOutput(ns("choose_columns")),
                 helpText("Select a column to be used as the Scoring variable")
             ),
             
             box(title = "Decoy Term",width = 12,
                 
                 # Input: Set Decoy Term ----
                 textInput(ns("decoy"), "Identifier", value = "DECOY"),
                 helpText("Enter the String term to distinguish
                                        a Decoy Entry from the Accession column")
             )
      ),
      
      column(width = 8,
             box(title= "Summary of Data Set",width = 12,
                 tabPanel('No filtering',style = 'overflow-x: scroll',       DT::dataTableOutput(ns('ex')))),
             box(h3("Data URL"),
                 verbatimTextOutput(ns("queryText")))
      )
    )
  )
}

## Server Module
dataPage <- function(input, output, session,current_dataSet_server_side) {
  ns <- session$ns
  
  # -------------------------------------------------------------------
  
  # Create the conductor for uploading a dataset
  current_dataSet <- reactive({
    get_current_dataSet(input$file1,passedUrlData(),queryLen())#Check the input type is valid - check if all checks are false
  })
  
  
  # Create a reactive conductor to take input from the score col selector and update the current data frame with it
  current_dataSet_server_side <- reactive({
    #return the server side version of the df with the stats calculated
    if(!is.null(current_dataSet()$pep)){
      returnCurrentServerDF(current_dataSet,input$column,input$decoy)
    }else{
      return(NULL)
    }
  })
  
  
  # -------------------------------------------------------------------
  
  ### Server Output section
  
  ## Drop down menu for scoring column choice
  
  output$choose_columns <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(current_dataSet()$pep)){
      validate(
        need(FALSE, "No data set uploaded ")
      )
      return()
    }
      
    
    colnames <- names(current_dataSet()$pep) # Get the data set with the appropriate name
    
    #get the selected option
    if("mzid.Scoring" %in% colnames){
      col_mat <- grepl("mzid.Scoring",colnames)
      selected_default <- colnames[which(col_mat)]
    }else if("X.10lgP" %in% colnames){
      col_mat <- grepl("X.10lgP",colnames)
      selected_default <- colnames[which(col_mat)]
    }else if("scr.PEAKS.peptideScore" %in% colnames){
      col_mat <- grepl("scr.PEAKS.peptideScore",colnames)
      selected_default <- colnames[which(col_mat)]
    }else{
      selected_default <- NULL
    }
    
    selectInput(ns("column"), "Choose column", as.list(colnames),selected = selected_default) #Create the drop down list
  })
  
  
  ## Table with turn off filtering (no searching boxes)
  
  output$ex <- DT::renderDataTable(
    if(is.null(current_dataSet()$pep)){
      validate(
        need(FALSE, "No data set uploaded ")
      )
    }else{
      DT::datatable(current_dataSet()$pep, options = list(searching = FALSE))
    }

  )
  
  # -------------------------------------------------------------------

  ### Validation Section for User selection
  
  ## Give user validation of score column selection
  observeEvent(input$column, {
    #warn user of column non numeric
    if(checkScoreColNum(current_dataSet()$pep,input$column)){
      shinyalert(title = "Warning!",text = "Score must be a numeric value", type = "warning")
    }
    
    showNotification("Score Column Set.",type = "message")
  })
  
  
  ## Give user validation of Decoy Term selection
  observeEvent(input$decoy, {
    showNotification("Decoy Term Set.",type = "message")
  })
  
  # -------------------------------------------------------------------
  
  ### Query String Parsing Section
  
  # Parse the GET query string
  output$queryText <- renderText({
    passedUrlData()
  })
  
  
  #reactive conductor to set the file pathway for 
  passedUrlData <- reactive({
    query <- session$clientData$url_search
    query <- returnDataUrl(query)
    # Return a string of data
    paste(query)
  })
  
  queryLen <- reactive(({
    query <- parseQueryString(session$clientData$url_search)
    paste(query)
  }))
  
  
  return(current_dataSet_server_side)
}



# -------------------------------------------------------------------

### Scat Page Display Modular Page


## Inner module for Scatter Count

innerGGPlotInput <- function(id){
  ns <- NS(id)
  tagList(
    plotOutput(ns("plot"))
  )
}


## Inner module for Scatter Count of Decoy
innerGGPlot <- function(input, output, session, plotting_func,current_dataSet_server_side,slider_range){
  output$plot <- renderPlot({
    
    validate(
      need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded ")
    )
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    plotting_func(current_dataSet_server_side(),slider_range())

  })
}



## Outer Module

scatDisInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  tagList(
    tabBox(
      title = "Retention Time vs Mass",
      tabPanel("Scatter Count Plot of Retention Time vs Mass Coloured by Score", title = "Score",
               innerGGPlotInput(ns("massRtbyP"))),
      tabPanel("Scatter Plot of Retention Time vs Mass, Coloured by Score & Decoy", title = "Score & Decoy",
               innerGGPlotInput(ns("massRtbyP&D")))
    ),
    tabBox(
      title = "Mass Error vs Mass/Charge",
      tabPanel("Scatter Count Plot of Mass Error vs Mass over Charge", title = "Score",
               innerGGPlotInput(ns("ppmMzbyP"))),
      tabPanel("Scatter Plot of Mass Error vs Mass over Charge, Coloured by Score & Decoy", title = "Score & Decoy",
               innerGGPlotInput(ns("ppmMzbyP&D")))
    ),
    box(width = 6,
        uiOutput(ns("choose_range"))
    ),
    box(width = 6, title = "Information",
        status = "info", solidHeader = TRUE,
        collapsible = TRUE,
        h5("Information Text"))
  )
}



## Module Server Function

scatDis <- function(input, output, session, current_dataSet_server_side) {
  ns <- session$ns

  
  
  # Drop down menu for scoring column choice
  output$choose_range <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(current_dataSet_server_side()))
      return()
    
    #Get Range
    range_to_use <- get_score_range(current_dataSet_server_side())
    #Create the Score slider
    sliderInput(ns("score_slider"), "Score Range:", range_to_use[1], range_to_use[2], range_to_use)
    
  })
  
  #Reactive slider input
  slider_range <- reactive({
    validate(
      need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded ")
    )
    input$score_slider
  })
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Making plot", value = 0.99)
  
  #Call the inner modules for ggploting
  callModule(innerGGPlot,"massRtbyP",plot_mass_rt_by_p,current_dataSet_server_side,slider_range)#Scatter plot of Mass vs RT coloured by score
  callModule(innerGGPlot,"massRtbyP&D",plot_scat_rt_mass_deocy,current_dataSet_server_side,slider_range)#Scatter plot of ppm vs m/z coloured by score differing with decoy
  callModule(innerGGPlot,"ppmMzbyP",plot_scatter_ppm_mz,current_dataSet_server_side,slider_range)#Scatter plot of ppm vs m/z coloured by score
  callModule(innerGGPlot,"ppmMzbyP&D",plot_scat_ppm_mz_deocy,current_dataSet_server_side,slider_range)

}



# -------------------------------------------------------------------

## Inner Module (per hist comb with slider box)

innerBarInput <- function(id){
  ns <- NS(id)
  tagList(
    box(width=12,
        plotlyOutput(ns("plot"))
    )
  )
}


innerBar <- function(input, output, session, current_dataSet_server_side, title, yAxisLab,col_id){
  output$plot <- renderPlotly({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    plot_bar(current_dataSet_server_side(),title,yAxisLab,col_id)
  })
}

innerHistInput <- function(id){
  ns <- NS(id)
  tagList(
    box(width=12,
        plotlyOutput(ns("plot")),
        sliderInput(ns("hist_slider"), "Number of Bins:", 20, 40, 2)
    )
  )
}


innerHist <- function(input, output, session, current_dataSet_server_side, title, yAxisLab,col_id){
  output$plot <- renderPlotly({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    plot_hist(current_dataSet_server_side(),input$hist_slider,title,yAxisLab,col_id)
  })
}

### Hist Outer Page Module

histInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  tagList(
    column(width = 6,
           innerHistInput(ns("ppm")),
           box(width = 12, title = "Information",
               status = "info", solidHeader = TRUE,
               collapsible = TRUE,
               h5("Information Text")),
           innerBarInput(ns("charge")),
           box(width = 12, title = "Information",
               status = "info", solidHeader = TRUE,
               collapsible = TRUE,
               h5("Information Text"))
    ),
    column(width = 6,innerHistInput(ns("mz")),
           box(width = 12, title = "Information",
               status = "info", solidHeader = TRUE,
               collapsible = TRUE,
               h5("Information Text")),
           innerHistInput(ns("RetT")),
           box(width = 12, title = "Information",
               status = "info", solidHeader = TRUE,
               collapsible = TRUE,
               h5("Information Text"))
    )
  )
}


hist <- function(input, output, session, current_dataSet_server_side) {
  callModule(innerHist,"ppm",current_dataSet_server_side,"Histogram of Mass Error","Mass Error","ppm")
  callModule(innerHist,"mz",current_dataSet_server_side,"Histogram of Mass over Charge","M/z","m.z")
  callModule(innerBar,"charge",current_dataSet_server_side,"Bar graph of Charge","Charge","z")
  callModule(innerHist,"RetT",current_dataSet_server_side,"Histogram of Retention Time","RT","RT")
}



# -------------------------------------------------------------------

### Cleavages Module Page

cleavPageInput <- function(id){
  ns <- NS(id)
  tagList(
    box(title = "Bar Graph of Cleavages",
        plotlyPlotInput(ns("cleavBar"),600)
    ),
    box(title = "Box Plot for Score Distrubution amongst cleavages",
        plotlyPlotInput(ns("cleavBox"),600)
    ),
    box(width = 9, title = "Information",
        status = "info", solidHeader = TRUE,
        collapsible = TRUE,
        h5("Information Text")),
    box(width = 3, title = "Toggle",status = "primary",
        checkboxInput(ns("checkbox"), label = "Show Decoys", value = TRUE))
  )
}

cleavPage <- function(input,output,session,current_dataSet_server_side){
  toggleDecoy <- reactive({input$checkbox})
  callModule(plotlyPlot2,"cleavBar",plot_bar_cleav,current_dataSet_server_side,toggleDecoy)# Bar of Cleavage distrubution
  callModule(plotlyPlot2,"cleavBox",plot_box_clev_by_score,current_dataSet_server_side,toggleDecoy)# Box plot of Cleavage distrubution amongst Score
  
}

# -------------------------------------------------------------------

### PTM page module: with a bar graph and a pie chart in plotly

ptmPageInput <- function(id){
  ns <- NS(id)
  tagList(
    
    column(width = 12,
           box(title = "Donut Chart of peptides with PTMs", width = 12,
               plotlyPlotInput(ns("ptmDon"),800)
           ),
           box(width = 12, title = "Information",
               status = "info", solidHeader = TRUE,
               collapsible = TRUE,
               h5("Information Text")))
  )
}


ptmPage <- function(input,output,session,current_dataSet_server_side){
  #callModule(plotlyPlot,"ptmBar",plot_bar_ptm,current_dataSet_server_side)# Plot: Bar of PTM
  callModule(plotlyPlot,"ptmDon",plot_donut_ptm,current_dataSet_server_side)# Plot: Donut of PTM
}

# -------------------------------------------------------------------

### Decoy Display

decoyDisInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  tagList(
    column(width = 8,
           box(title = "Scatter plot of Mass error (ppm) and score (-10lgP)",width = 12,
               plotlyPlotInput(ns("scorePpmScat"),400)),
           box(title = "Box Plot Marginal Denstity of Mass Error",width = 12,
               plotlyPlotHeightInput(ns("scoreBox"),200))
    ),
    column(width = 4,
           box(title = "Box Plot Marginal Denstity of Score", width = 12,
               plotlyPlotInput(ns("ppmBox"),200)),
           box(width = 12, title = "Information",
               status = "info", solidHeader = TRUE,
               collapsible = TRUE,
               h5("Information Text"))
    )
  )
}


decoyDis <- function(input, output, session, current_dataSet_server_side) {
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Making plot", value = 0.99)
  callModule(plotlyPlot,"scorePpmScat",plot_scat_score_ppm_by_decoy,current_dataSet_server_side)# Plot: Scatter Plot of score and ppm
  callModule(plotlyPlot,"scoreBox",plot_box_marg_score,current_dataSet_server_side)# Plot: Box plot Score
  callModule(plotlyPlot,"ppmBox",plot_box_marg_ppm,current_dataSet_server_side)# Plot: Box plot Score
}


# -------------------------------------------------------------------

### Module for the score page

scorePageDisplayInput <- function(id){
  ns <- NS(id)
  tagList(box(width = 12,plotlyPlotInput(ns("scorePep"))),
          box(width = 12, title = "Information",
              status = "info", solidHeader = TRUE,
              collapsible = TRUE,
              h5("The plot shows distubution of scores, by count, split by whether it is Target or Decoy peptide.")),
          column(width = 6,box(title = "FDR Curve",width = 12,plotlyPlotInput(ns("FDRcurvePep"))),
                 box(width = 12, title = "Information",
                     status = "info", solidHeader = TRUE,
                     collapsible = TRUE,
                     h5("Information Text"))),
          column(width = 6,box(title = "q Value Curve",width = 12,plotlyPlotInput(ns("QcurvePep"))),
                 box(width = 12, title = "Information",
                     status = "info", solidHeader = TRUE,
                     collapsible = TRUE,
                     h5("Information Text")))
  )
}

scorePageDisplay <- function(input,output,session,current_dataSet_server_side){
  
  callModule(plotlyPlot,"FDRcurvePep",plot_FDR_curve,current_dataSet_server_side)#false discovery rate of peptides
  callModule(plotlyPlot,"QcurvePep",plot_Q_curve,current_dataSet_server_side)#q val curve
  callModule(plotlyPlot,"scorePep",plot_hist_score,current_dataSet_server_side)#Score distrubution page
}

# -------------------------------------------------------------------

### Module to produce  plotly plot

plotlyPlotInput <- function(id,plotHeight){
  ns <- NS(id)
  tagList(
    plotlyOutput(ns("plot"))
  )
}

#Input UI module with a height arguement
plotlyPlotHeightInput<- function(id,plotHeight){
  ns <- NS(id)
  tagList(
    plotlyOutput(ns("plot"),height = plotHeight)
  )
}

#Standard plotly server function
plotlyPlot <- function(input, output, session, ploting_function,current_dataSet_server_side)  {
  

  output$plot <- renderPlotly({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    ploting_function(current_dataSet_server_side())
  })
  
}

#Server with an extra arguement
plotlyPlot2 <- function(input, output, session, ploting_function,current_dataSet_server_side,decoyToggle)  {
  

  
  output$plot <- renderPlotly({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    ploting_function(current_dataSet_server_side(),decoyToggle())
  })
  
}


# -------------------------------------------------------------------