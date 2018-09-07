##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Page containing code for Shiny Modules
## Date: 24.08.2018
## Author: Ashleigh Myall
##################################################



# -------------------------------------------------------------------

### Modules for the crowd source running statistics page

# First page with a background of the map, with points of user nodes. Also contains a side bar with server stats

## UI
homePageInput <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(column(width = 12,box(width = 12, leafletOutput(ns("mymap"), height = 500))),
             column(width = 8,box(title= "Top Users by location",width = 12,
                                  tabPanel('No filtering',style = 'overflow-x: scroll',
                                           DT::dataTableOutput(ns('ex'))))),
             column(width = 4,
                    box(width = 12,
                        plotOutput(ns("plot")),
                        verbatimTextOutput(ns("txtout"))
                        )
                    )
    )
  )
}

## Server
homePage <- function(input, output, session){
  
  # -------------------------------------------------------------------
  
  ### Query String Parsing Section
  
  ## reactive conductor to set the file pathway for 
  passedUrlData <- reactive({paste(returnServerDataCsvUrl(session$clientData$url_search))})
  
  queryLen <- reactive(({parseQueryString(session$clientData$url_search)}))
  
  # -------------------------------------------------------------------
  
  data <- reactive({  #Get the dataset from the URL
    if(!length(queryLen())<1){
      tryCatch({
        read.csv(url(passedUrlData()))
      },error = function(e) {return(NULL)})# return a safeError if a parsing error occurs
    }else{
      return(NULL)
    }
  })
  
  # -------------------------------------------------------------------
  
  output$mymap <- renderLeaflet({ #shiny output for the map with points from dataset
    if(!is.null(data())){
      df <- data()
      m <- leaflet(data = df) %>% addTiles() %>%
        addMarkers(lng = ~Longitude,
                   lat = ~Latitude,
                   popup = paste("City: ", df$City,"\n","Results: ", df$Results)) %>%
        addProviderTiles(providers$CartoDB.Positron)
      
    }else{
      m <- leaflet() %>%addTiles() %>% addProviderTiles(providers$CartoDB.Positron)
    }
    m
  })
  
  # -------------------------------------------------------------------
  
  output$plot <- renderPlot({ # Placeholder bar plot for server stats
    df <- data()
    plotDf <- data.frame("Type" = c("Requests","Results"),"Num" = c(sum(df$Requests),sum(df$Results)))
    ggplot(plotDf , aes(Type, Num)) +
      geom_bar(stat="identity")+ 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  })
  
  # -------------------------------------------------------------------
  
  output$txtout <- renderText({ #Output text stats of server info
    df <- data()
    paste("Total Processes per day:",sum(df$Results),"\n","Total Users per day: ",nrow(df),
          "\n","Top user of the day: ",df[which.max(df$Results),]$Users,"\n","Requests to results: : ",round(100*sum(df$Results)/sum(df$Requests),digits = 2),"%")
  })
  
  # -------------------------------------------------------------------
  
  ## Table with turn off filtering (no searching boxes)
  output$ex <- DT::renderDataTable({
    if(is.null(data())){validate(need(FALSE, "No data set uploaded "))}
    df <- data()
    #create a new dataframe for display
    dfDisplay <- data.frame("City" = df$City,"Country" = df$Country,"Requests" = df$Requests,"Results" = df$Results)
    #sort and head the top 10
    dfDisplay <- head(dfDisplay[order(dfDisplay$Results)])
    dfDisplay$Country <- NULL
    DT::datatable(dfDisplay, options = list(searching = FALSE))
  })
}

# -------------------------------------------------------------------

### Module for the Data input page

#Contains inputs for users selection. By default no data is uploaded so no table of data will be displayed

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
                 helpText("Accepted formats are : csv, mzid and gz. 
                          To be processed properly the file should contain a column labeled Peptide or Sequence."),
                 downloadButton(ns('downloadData'), 'Download')
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
             ),
             box(width = 12,
                 solidHeader = TRUE,
                 collapsible = TRUE,
                 h5("Users can upload a data set (of file types csv,mzid and gz) and view them in the table. Options are available to select the column for scoting and set a decoy marker term.
                    Note that the decoy marker term will be taken from the Accession column if no isDecoy column is present."),
                 h5("This application also filters input by a max rank of 1 and selects unique peptide sequences(selecting the peptide with the highest score. This column is the 2nd column present in the scoring file for mzIdentML 
                    and set as X.10lgP or the first column containing the string Score for a given csv. Please ensure a column contains this string otherwise filtering for unique peptides will not be done.)"))
      ),
      
      column(width = 8,
             box(title= "Summary of Data Set",width = 12,
                 tabPanel('No filtering',style = 'overflow-x: scroll',       DT::dataTableOutput(ns('ex'))))
      )
    )
  )
}

## Server Module
dataPage <- function(input, output, session,current_dataSet_server_side) {
  ns <- session$ns
  
  # -------------------------------------------------------------------
  
  ## Create the conductor for uploading a dataset
  current_dataSet <- reactive({
    get_current_dataSet(input$file1,passedUrlData(),queryLen())#Check the input type is valid - check if all checks are false
  })
  
  ## Create a reactive conductor to take input from the score col selector and update the current data frame with it
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
      validate(need(FALSE, "No data set uploaded "))
      return()
    }
    selectInput(ns("column"), "Choose column", as.list(colnames(current_dataSet()$pep)),selected = getInitScoreCol(current_dataSet()$pep)) #Create the drop down list
  })
  
  
  
  ## Table with turn off filtering (no searching boxes)
  output$ex <- DT::renderDataTable(
    if(is.null(current_dataSet()$pep)){validate(need(FALSE, "No data set uploaded "))}
    else{DT::datatable(current_dataSet()$pep, options = list(searching = FALSE))}
  )
  
  # -------------------------------------------------------------------

  ### Validation Section for User selections
  
  ## Give user validation of score column selection 
  observeEvent(input$column, {  #warn user of column non numeric
    if(checkScoreColNum(current_dataSet()$pep,input$column)){
      shinyalert(title = "Warning!",text = "Score must be a numeric value",type = "warning")
    }else{
      showNotification("Score Column Set.",type = "message")
    }
  })
  
  ## Give user validation of Decoy Term selection
  observeEvent(input$decoy, {showNotification("Decoy Term Set.",type = "message")})
  
  # -------------------------------------------------------------------
  
  ### Query String Parsing Section
  
  ## reactive conductor to set the file pathway for 
  passedUrlData <- reactive({
    query <- session$clientData$url_search
    query <- returnDataUrl(query)
    paste(query) # Return a string of data
  })
  
  queryLen <- reactive(({parseQueryString(session$clientData$url_search)}))
  
  # -------------------------------------------------------------------
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
  output$downloadData <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste(ifelse(is.null(input$file1),"Crowdsource",tools::file_path_sans_ext(input$file1)),".csv")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file' if the dataSet is not null
    content = function(file) {
      # Write to a file specified by the 'file' argument
      write.table(current_dataSet()$pep, file, sep = ",",
                  row.names = FALSE)
    }
  )
  
  # -------------------------------------------------------------------
  
  ### Final return of the server DF for other modules to access
  return(current_dataSet_server_side)
}



# -------------------------------------------------------------------

### Scatters Page Display Modular Page

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
    validate(need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded "))
    progress <- shiny::Progress$new()# Create a Progress object
    on.exit(progress$close())# Make sure it closes when we exit this reactive, even if there's an error
    progress$set(message = "Making plot", value = 0.99)
    plotting_func(current_dataSet_server_side(),slider_range())

  })
}


# -------------------------------------------------------------------

### Module for single scatter with options to filter and choose cats

#UI
singleScatPageInput <- function(id){
  ns <- NS(id)
  
  tagList(
    fluidRow(column(width = 9,
                    box(title = "Query Plotter Output",width = 12, 
                        box(width = 12,
                            innerVariableGgplotInput(ns("variableGgPlot"))
                        )
                    )
    ),
    column(width = 3,
           box(title = "Query Builder", 
               width = 12,
               status = "primary",
               solidHeader = FALSE,
               collapsible = TRUE,
               background = "navy",
               box(width = 12,
                   status = "primary",
                   solidHeader = FALSE,
                   collapsible = FALSE,
                   background = "navy",
                   uiOutput(ns("choose_columns_X"))),
               box(width = 12,
                   status = "primary",
                   solidHeader = FALSE,
                   collapsible = FALSE,
                   background = "navy",
                   uiOutput(ns("choose_columns_Y"))),
               box(width = 12,
                   status = "primary",
                   solidHeader = FALSE,
                   collapsible = FALSE,
                   background = "navy",
                   uiOutput(ns("choose_range"))),
               box(width = 12,
                   status = "primary",
                   solidHeader = FALSE,
                   collapsible = FALSE,
                   background = "navy",
                   radioButtons(ns("zAxis"), "z axis:",
                                c("Score" = FALSE,
                                  "Score and Decoy" = TRUE))),
               box(width = 12,
                   status = "primary",
                   solidHeader = FALSE,
                   collapsible = FALSE,
                   background = "navy",
                   checkboxGroupInput(ns("axisScale"), "Scale:",c("X axis Log" = "x","Y axis Log" = "y")))
           ),
           box(width = 12,
               solidHeader = TRUE,
               collapsible = TRUE,
               h5("This plot is a Scatter where selection of the X and Y axis can be changed."),
               h5("Also available is the option to filter the score range of points displayed and separate point colour and shape by whether it is a Decoy"))
    ))
  )
}

#Server
singleScatPage <- function(input, output, session, current_dataSet_server_side){
  ns <- session$ns
  
  # -------------------------------------------------------------------
  
  ## Outputs
  
  #X axis selection choice
  output$choose_columns_X <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(current_dataSet_server_side()$pep)){
      validate(need(FALSE, "No data set uploaded "))
      return()
    }
    selectInput(ns("choose_columns_X"), "Choose X axis", as.list(getColNames(current_dataSet_server_side()$pep)),
                selected = "RT") #Create the drop down list
  })
  
  
  #Y axis selection choice
  output$choose_columns_Y <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(current_dataSet_server_side()$pep)){
      validate(need(FALSE, "No data set uploaded "))
      return()
    }
    selectInput(ns("choose_columns_Y"), "Choose Y axis", as.list(getColNames(current_dataSet_server_side()$pep)),selected = "Mass") #Create the drop down list
  })
  
  
  
  # scoring range choice
  output$choose_range <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(current_dataSet_server_side())){
      validate(need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded "))
      return()
    }

    #Get Range
    range_to_use <- get_score_range(current_dataSet_server_side())
    #Create the Score slider
    sliderInput(ns("score_slider"), "Score Range:", range_to_use[1], range_to_use[2], range_to_use)
    
  })
  
  
  # -------------------------------------------------------------------
  
  #Reactive X axis selection
  x_axis_selection <- reactive({
    validate(need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded "))
    input$choose_columns_X
  })
  
  #Reactive Y axis selection
  y_axis_selection <- reactive({
    validate(need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded "))
    input$choose_columns_Y
  })
  
  #Reactive slider input
  slider_range <- reactive({
    validate(need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded "))
    input$score_slider
  })
  
  #Reactive decoy input
  zChoice <- reactive({
    validate(need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded "))
    input$zAxis
  })
  
  axisScale <- reactive({input$axisScale})
  
  # Call the module for plotting
  callModule(innerVariableGgplot,"variableGgPlot",current_dataSet_server_side,slider_range,zChoice,x_axis_selection,y_axis_selection,axisScale)
  
}


## Inner module for the variable ggplot
innerVariableGgplotInput <- function(id){
  ns <- NS(id)
  tagList(
    plotOutput(ns("plot"),height = 600)
  )
}

innerVariableGgplot <- function(input, output, session, current_dataSet_server_side,slider_range,zChoice,x_axis_selection,y_axis_selection,axisScale){
  output$plot <- renderPlot({
    validate(need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded "))
    progress <- shiny::Progress$new() # Create a Progress object
    on.exit(progress$close()) # Make sure it closes when we exit this reactive, even if there's an error
    progress$set(message = "Making plot", value = 0.99)
    plotVariableScatter(current_dataSet_server_side(),slider_range(),zChoice(),x_axis_selection(),y_axis_selection(),axisScale())
    
  })
}

# -------------------------------------------------------------------

###Modules for histogram page

## Inner Module (per hist comb with slider box)
innerBarInput <- function(id){
  ns <- NS(id)
  tagList(
    box(width=12,
        plotlyOutput(ns("plot"))
    )
  )
}

## Inner bar Server Module
innerBar <- function(input, output, session, current_dataSet_server_side, title, yAxisLab,col_id){
  output$plot <- renderPlotly({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    plot_bar(current_dataSet_server_side(),title,yAxisLab,col_id)
  })
}

## Inner Hist UI
innerHistInput <- function(id){
  ns <- NS(id)
  tagList(
    box(width=12,
        plotlyOutput(ns("plot")),
        sliderInput(ns("hist_slider"), "Number of Bins:", 20, 40, 2)
    )
  )
}

##Inner Hist Server
innerHist <- function(input, output, session, current_dataSet_server_side, title, yAxisLab,col_id){
  output$plot <- renderPlotly({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    plot_hist(current_dataSet_server_side(),input$hist_slider,title,yAxisLab,col_id)
  })
}

## Hist Outer Page Module
histInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  tagList(
    fluidRow(column(width = 6,
                    innerHistInput(ns("ppm")),
                    innerBarInput(ns("charge")),
                    box(width = 12,
                        status = "info", solidHeader = TRUE,
                        collapsible = TRUE,
                        h5("This page contains multiple plots for the distributions, with range sliders to adjust the bin width of the plots"))
    ),
    column(width = 6,innerHistInput(ns("mz")),
           innerHistInput(ns("RetT"))
    ))
  )
}

#hist srver module
hist <- function(input, output, session, current_dataSet_server_side) {
  callModule(innerHist,"ppm",current_dataSet_server_side,"Histogram of Mass Error","Mass Error","ppm")
  callModule(innerHist,"mz",current_dataSet_server_side,"Histogram of Mass over Charge","M/z","m.z")
  callModule(innerBar,"charge",current_dataSet_server_side,"Bar graph of Charge","Charge","z")
  callModule(innerHist,"RetT",current_dataSet_server_side,"Histogram of Retention Time","RT","RT")
}



# -------------------------------------------------------------------

### Cleavages Module Page

#UI
cleavPageInput <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 12,
             box(title = "Bar Graph of Cleavages",
                 plotlyPlotInput(ns("cleavBar"))
             ),
             box(title = "Box Plot for Score Distrubution amongst cleavages",
                 plotlyPlotInput(ns("cleavBox"))
             )
             ),
      column(width = 9,
             box(width = 12,
                 status = "info", solidHeader = TRUE,
                 collapsible = TRUE,
                 h5("Shown in the first figure is the relative count of peptide with x cleavages. Shown in the second figure is a box plot of peptides, per number of cleavages, against score.")
                 )
             ),
      column(width = 3,
             box(title = "Query Builder", 
                 width = 12,
                 status = "primary",
                 solidHeader = FALSE,
                 collapsible = TRUE,
                 background = "navy",
                 box(width = 12,
                     status = "primary",
                     solidHeader = FALSE,
                     collapsible = FALSE,
                     background = "navy",
                     checkboxInput(ns("checkbox"), label = "Show Decoys", value = TRUE))
             )
        )
    )
  )
}

#Server
cleavPage <- function(input,output,session,current_dataSet_server_side){
  toggleDecoy <- reactive({input$checkbox})
  callModule(plotlyPlot2,"cleavBar",plot_bar_cleav,current_dataSet_server_side,toggleDecoy)# Bar of Cleavage distrubution
  callModule(plotlyPlot2,"cleavBox",plot_box_clev_by_score,current_dataSet_server_side,toggleDecoy)# Box plot of Cleavage distrubution amongst Score
}

# -------------------------------------------------------------------

### PTM page module: with a bar graph and a pie chart in plotly

#UI
ptmPageInput <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 12,
             box(width = 6,
                 plotlyPlotInput(ns("ptmBar"))
             ),
             box(width = 6,
                 plotlyPlotInput(ns("ptmModCountBar"))
             ),
             box(width = 9,
                 status = "info", solidHeader = TRUE,
                 collapsible = TRUE,
                 h5("The peptides with modifications bar chart is showing the total count of peptides with that associated modification. (Note: peptides can contain multiple modifications of the same type and/or combinations). 
                    The second plot is a bar chart, displaying how many peptides contain x modifications (Note: x ranged fixed at max 9 modifications for visualisation).")
                 )
             )
    )
  )
}

#Server
ptmPage <- function(input,output,session,current_dataSet_server_side){
  callModule(plotlyPlot,"ptmBar",plot_bar_ptm,current_dataSet_server_side)# Plot: Bar of PTM
  callModule(plotlyPlot,"ptmModCountBar",plot_bar_ptm_mod_count,current_dataSet_server_side)# Plot: Bar of PTM modification count
}

# -------------------------------------------------------------------

### Decoy Display

#UI
decoyDisInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(width = 8,
             box(title = "Scatter plot of Mass error (ppm) and score (-10lgP)",width = 12,
                 plotlyPlotInput(ns("scorePpmScat"))),
             box(title = "Box Plot Marginal Denstity of Mass Error",width = 12,
                 plotlyPlotHeightInput(ns("scoreBox"),200)),
             box(width = 12,
                 status = "info", solidHeader = TRUE,
                 collapsible = TRUE,
                 h5("The above plots show the distribution of Decoy and Target peptides, against score and ppm.
                    Included is a query builder where you can adjust the opacity and filter by Decoy and Targets.
                    "))
      ),
      column(width = 4,
             box(title = "Box Plot Marginal Denstity of Score", width = 12,
                 plotlyPlotInput(ns("ppmBox"))),
             box(title = "Query Builder", 
                 width = 12,
                 status = "primary",
                 solidHeader = FALSE,
                 collapsible = TRUE,
                 background = "navy",
                 box(width = 12,
                     status = "primary",
                     solidHeader = FALSE,
                     collapsible = FALSE,
                     background = "navy",
                     sliderInput(ns("alpha"), "Point Opacity:",
                                 min = 0, max = 1,
                                 value = 0.5)),
                 box(width = 12,
                     status = "primary",
                     solidHeader = FALSE,
                     collapsible = FALSE,
                     background = "navy",
                     radioButtons(ns("pointsDisplayed"), "Points Displayed:",
                                  c("Combined" = 1,
                                    "Target" = 2,
                                    "Decoy" = 3)))
             )
      )
    )
  )
}

#Server
decoyDis <- function(input, output, session, current_dataSet_server_side) {
  opacity <- reactive({input$alpha}) #Reactive cond for alpha
  pointsToDisplay <- reactive({input$pointsDisplayed}) #Reactive cond for points display option
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Making plot", value = 0.99)
  callModule(plotlyPlot3,"scorePpmScat",plot_scat_score_ppm_by_decoy,current_dataSet_server_side,opacity,pointsToDisplay)# Plot: Scatter Plot of score and ppm
  callModule(plotlyPlot,"scoreBox",plot_box_marg_score,current_dataSet_server_side)# Plot: Box plot Score
  callModule(plotlyPlot,"ppmBox",plot_box_marg_ppm,current_dataSet_server_side)# Plot: Box plot Score
}


# -------------------------------------------------------------------

### Module for the score page

#UI
scorePageDisplayInput <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(width = 12,plotlyPlotInput(ns("scorePep"))),
      box(width = 12,
          status = "info", solidHeader = TRUE,
          collapsible = TRUE,
          h5("This page displays 3 plots. The first above is the distribution of peptides, by score, where bins are stacked by whether the peptide is a Decoy or not. The second plot below is the FDR curve and the third is the Q value curve. 
A query builder tab on the right allows for custom setting of the FDR percentage, adjusting this changes the the position of the markers on the first two plots and the PSM count.
             ")),
      
      column(width = 9,
             tabBox(
               title = "Curve Plots",width = 12,
               # The id lets us use input$tabset1 on the server to find the current tab
               id = "tabset1",
               tabPanel("FDR Curve", 
                        plotlyPlotInput(ns("FDRcurvePep"))
               ),
               tabPanel("Q Value Curve", 
                        plotlyPlotInput(ns("QcurvePep"))
               )
             )
      ),
      
      column(width = 3,
             box(title = "Query Builder", 
                 width = 12,
                 status = "primary",
                 solidHeader = FALSE,
                 collapsible = TRUE,
                 background = "navy",
                 box(width = 12,
                     status = "primary",
                     solidHeader = FALSE,
                     collapsible = FALSE,
                     background = "navy",
                     numericInput(ns("FDRpercent"), label = "FDR Percent", value = 1,min = 0.1,max = 100,step = 0.1)),
                 box(width = 12,
                     status = "primary",
                     solidHeader = FALSE,
                     collapsible = FALSE,
                     background = "navy",
                     textOutput(ns("hitsAboveFDR"))),
                 box(width = 12,
                     status = "primary",
                     solidHeader = FALSE,
                     collapsible = FALSE,
                     background = "navy",
                     plotlyPlotHeightInput(ns("fdrPie"),200))
             )
      )
    )
  )
}

#Server
scorePageDisplay <- function(input,output,session,current_dataSet_server_side){
  #Set input for reactive
  fdrPercent <- reactive({
    input$FDRpercent
  })
  
  # Out text rendering
  output$hitsAboveFDR <- renderText({ 
    validate(
      need(!is.null(current_dataSet_server_side()$pep), "No data set uploaded ")
    )
    paste("PSM Count: ",getIntercept(current_dataSet_server_side()$pep,input$FDRpercent))
  })
  
  callModule(FDRPlot,"FDRcurvePep",plot_FDR_curve,current_dataSet_server_side,fdrPercent)#false discovery rate of peptides
  callModule(FDRPlot,"QcurvePep",plot_Q_curve,current_dataSet_server_side,fdrPercent)#q val curve
  callModule(FDRPlot,"scorePep",plot_hist_score,current_dataSet_server_side,fdrPercent)#Score distrubution page
  callModule(FDRPlot,"fdrPie",plotIdentFdr,current_dataSet_server_side,fdrPercent)#FDR pie
}

#Inner Module Server for the FDR curve plot
FDRPlot <- function(input, output, session,plotting_func,current_dataSet_server_side,fdrPercent){
  output$plot <- renderPlotly({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    plotting_func(current_dataSet_server_side(),fdrPercent())
  })
}

# -------------------------------------------------------------------

### Module to produce  plotly plot

#UI
plotlyPlotInput <- function(id){
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


#Server with multiple arguements - for the scatter of decoys
plotlyPlot3 <- function(input, output, session, ploting_function,current_dataSet_server_side,opacity,pointsToDisplay)  {
  output$plot <- renderPlotly({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0.99)
    ploting_function(current_dataSet_server_side(),opacity(),pointsToDisplay())
  })
}


# -------------------------------------------------------------------