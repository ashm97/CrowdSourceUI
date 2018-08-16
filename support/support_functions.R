##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Supporting functions for the app
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

# -------------------------------------------------------------------

### Function to return a subsetted dataframe with new columns calculating
### for statistical anal. Takes only one column of the orignal data frame ('FP')
### which is used to determine between Decoys & hits and takes the string
### used to class the FP rows (default = 'DECOY')

get_df_stat <- function(df_FP_only, decoy_string){
  
  df_FP_only <- data.frame(df_FP_only)
  colnames(df_FP_only) <- c("Accession")
  
  # Add an FP column which represents is a DECOY or hit, 1 for decoys
  df_FP_only <- transform(df_FP_only, FP = ifelse(grepl(decoy_string, Accession),1,0))
  
  # Add another column for Decoy or Target
  df_FP_only <- transform(df_FP_only, Decoy = ifelse(grepl(decoy_string, Accession),"Decoy","Target"))
  
  # Add cumulative column to the data frame
  df_FP_only[,"Cum.FP"] <- cumsum(df_FP_only$FP)
  
  # Add a column for the num of hits above a thresh
  df_FP_only$Hits.Above.Tresh <- seq(1,nrow(df_FP_only),1)
  
  # Add a column for TP which is = Hits.Above.Tresh - 2*Cum.FP
  df_FP_only <- transform(df_FP_only, TP = Hits.Above.Tresh - 2*Cum.FP)
  
  ### Add Column for FDR
  df_FP_only <- transform(df_FP_only, FDR = Cum.FP / (TP  + Cum.FP))
  
  # Add a column for the Q value
  df_FP_only$Q.val <- return_Q_Val(df_FP_only$FDR,nrow(df_FP_only))
  
  return(df_FP_only)
}

### Stats function where there exists an isDecoy Col

getStatsDfExistingDecoy <- function(df_FP_only){
  
  df_FP_only <- data.frame(df_FP_only)
  colnames(df_FP_only) <- c("isDecoy")
  
  df_FP_only <- transform(df_FP_only, FP = ifelse(isDecoy,1,0))
  df_FP_only <- transform(df_FP_only, Decoy = ifelse(isDecoy,"Decoy","Target"))
  df_FP_only$isDecoy <- NULL
  df_FP_only[,"Cum.FP"] <- cumsum(df_FP_only$FP)
  df_FP_only$Hits.Above.Tresh <- seq(1,nrow(df_FP_only),1)
  df_FP_only <- transform(df_FP_only, TP = Hits.Above.Tresh - 2*Cum.FP)
  df_FP_only <- transform(df_FP_only, FDR = Cum.FP / (TP  + Cum.FP))
  df_FP_only$Q.val <- return_Q_Val(df_FP_only$FDR,nrow(df_FP_only))
  
  return(df_FP_only)
}

# -------------------------------------------------------------------

### Function to calculate the column for Q value:
### Cycles through taking the current position to the end selecting the
### min value from FDR.                                                   

return_Q_Val <- function(dataCol,n){
  return_col <- rep(0,n)
  for(i in 1:n){
    return_col[i] <- min(dataCol[i:n])
  }
  return(return_col)
}


# -------------------------------------------------------------------

### Function to count the cleavages in each entry row
# Returns a column of the num of cleavages in the adjacent column

get_pep_cleav <- function(peptides_to_process){
  
  #create the functions internal dataframe to handle 
  clev_df <- data.frame(peptides_to_process)
  colnames(clev_df) <- c("Peptide")
  
  #create a new column cutting the front off of any matching peptides
  clev_df <- transform(clev_df, mod_p = gsub("^R|^K","A", Peptide))
  
  #Create a new column cutting the end off of any matching peptides
  clev_df <- transform(clev_df, mod_p_2 = gsub("R$|K$","A", mod_p))
  
  #Create a col replacing where K followed by a P
  clev_df <- transform(clev_df, mod_p_3 = gsub("KP","A", mod_p_2))
  
  #Create a new col which has a count of the number of strings with a space seperating them
  clev_df <- transform(clev_df, pep_clev_count = str_count(mod_p_3,"R|K"))
  
  return(clev_df$pep_clev_count)
}


# -------------------------------------------------------------------

## Function to get the current data Set

get_current_dataSet <- function(in_File,passedUrlData,queryLen){
  # Selecting from the default or the uploaded
  if (is.null(in_File)&(length(queryLen)<1)){ # if both no file has been uploaded and no file passed by the url
    
    #Return data as a list without mods
    returnList <- list("pep" = NULL, "mod" = NULL)
    return(returnList) # give the default dataSet
    
  }else{
    #if no file has been uploaded and there is a passed URL file, then set in_File to the passedUrl and validate
    if(is.null(in_File)&(length(queryLen)>=1)){
      
      returnList <- list("pep" = get_url_dat(passedUrlData), "mod" = NULL)
      return(returnList)
    }
    #check in File
    if(validate_file(in_File)){
      if(grepl("csv",in_File$datapath)){
        
        returnList <- list("pep" = handleFileCsv(in_File), "mod" = NULL)
        return(returnList)
        
      }else if(grepl("mzid",in_File$datapath)){ #IF .mzid
        
        returnList <- list("pep" = handleFileMzid(in_File), "mod" = NULL)
        return(returnList)
      }else{
        #na (gzip currently which is handled the same as mzid)
        
        returnList <- list("pep" = handleFileMzid(in_File), "mod" = NULL)
        return(returnList)
      }
      
    }else{
      shinyalert(title = "Bad upload file!",text = "Please upload a file of: .csv .mzid .gz", type = "warning")
      
      returnList <- list("pep" = returnDefaultDf(default_dataSet), "mod" = NULL)
      return(returnList) # give the default dataSet
    }
  }
}


# -------------------------------------------------------------------

## Function to create a new dataframe for behinds the scenes which has the score set, has been ordered by score

get_dataSet_withScore <- function(df_to_use,selected_col){
  
  #Create a new column of Charge
  try(df_to_use <- transform(df_to_use, z = round(Mass / m.z)))
  
  #Round the ppm col to 4dp
  try(df_to_use <- transform(df_to_use, ppm = round(ppm,digits = 4)))
  
  #set the score column
  colnames(df_to_use)[colnames(df_to_use)==selected_col] <- "Score"
  
  #order by score
  df_to_use <- df_to_use[rev(order(df_to_use$Score)),]
  
  return(df_to_use)
  
  
}

# -------------------------------------------------------------------

## Function to return min, max and a default 

get_score_range <- function(df_to_use){
  
  rang <- c(min(df_to_use$pep$Score),max(df_to_use$pep$Score))
  
  return(rang)
  
}

# -------------------------------------------------------------------

## mzid and gz file input handler

handleFileMzid <- function(inputFile){
  mzid <- openIDfile(inputFile$datapath)
  
  #Score dataframe
  scr <- score(mzid)
  mzid_df <- data.frame(psms(mzid),scr[,2])
  
  
  #Format the cols for the app
  colnames(mzid_df)[colnames(mzid_df)=="scr...2."] <- "mzid.Scoring"
  colnames(mzid_df)[colnames(mzid_df)=="DatabaseAccess"] <- "Accession"
  colnames(mzid_df)[colnames(mzid_df)=="sequence"] <- "Peptide"
  colnames(mzid_df)[colnames(mzid_df)=="experimentalMassToCharge"] <- "m.z"
  
  # Add a col for ppm
  mzid_df <- transform(mzid_df, ppm = ((m.z - calculatedMassToCharge)*1000000)/m.z)
  
  #remove this col so easier to view - ask andy if this is okay
  #mzid_df$DatabaseSeq <- NULL
  mzid_df$counted_cleavages <- as.factor(get_pep_cleav(mzid_df$Peptide))
  return(mzid_df)
}

# -------------------------------------------------------------------

## Function to return the CSV dataset

handleFileCsv <- function(in_File){
  
  uploaded_df <- read.csv(in_File$datapath, header=TRUE)
  
  ## Count the cleavages and add to the data set as a new column
  uploaded_df$counted_cleavages <- as.factor(get_pep_cleav(uploaded_df$Peptide))
  
  return(uploaded_df)
}

# -------------------------------------------------------------------

## Function to upload data from web

get_url_dat <- function(webAddress){
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Reading in Data", value = 0.99)
  webDat <- read.csv(url(webAddress))
  webDat$counted_cleavages <- as.factor(get_pep_cleav(webDat$Peptide))
  return(webDat)
}

# -------------------------------------------------------------------

## Function to check the file uploaded is valid

validate_file <- function(inputFile){
  
  #First check to make sure of the right file format
  if(grepl(".csv",inputFile$datapath)|grepl(".mzid",inputFile$datapath)|grepl(".gz",inputFile$datapath)){
    #Contains one of the right formats
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# -------------------------------------------------------------------

## function to check if score the score column is a numeric

checkScoreColNum <- function(df_to_use,selected_col){
  
  #set the score column
  colnames(df_to_use)[colnames(df_to_use)==selected_col] <- "Check"
  
  return(!is.numeric(df_to_use$Check))
}

# -------------------------------------------------------------------

## function to return the url correct form

returnDataUrl <- function(string){
  string <- sub("[?]","",string)
  string <- sub("data=","",string)
  return(string)
}

# -------------------------------------------------------------------

## Function to heck the columns required exist and if not return an error message instead of plot

checkDataUploaded <- function(df){
  validate(
    need(is.null(df), erMessage)
  )
}

checkColExist <- function(dfCol,erMessage){
  validate(
    need(dfCol != "", erMessage)
  )
}

# -------------------------------------------------------------------

## Function to return the current server side version of the dataset

returnCurrentServerDF <- function(current_dataSet,col,decoyString){
  
  #Set the score col
  returnPepDf <- get_dataSet_withScore(current_dataSet()$pep,col)
  
  #get the stats cols
  if("isDecoy" %in% colnames(returnPepDf)){
    statsDf <-getStatsDfExistingDecoy(returnPepDf$isDecoy)
  }else{
    statsDf <- get_df_stat(returnPepDf$Accession,decoyString)
  }
  
  #Combine pep and stats df
  returnPepDf <- cbind.data.frame(returnPepDf,statsDf)
  
  #return all
  returnList <- list("pep" = returnPepDf, "mod" = current_dataSet()$mod)
  return(returnList)
  
}


