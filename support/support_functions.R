##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Supporting functions for the app
## Date: 23.08.2018
## Author: Ashleigh Myall
##################################################



# -------------------------------------------------------------------

### Function to return a subsetted dataframe with new columns calculating

# for statistical anal. Takes only one column of the orignal data frame ('FP')
# which is used to determine between Decoys & hits and takes the string
# used to class the FP rows (default = 'DECOY')

get_df_stat <- function(df_FP_only, decoy_string){
  df_FP_only <- data.frame(df_FP_only)
  colnames(df_FP_only) <- c("Accession")
  df_FP_only <- transform(df_FP_only, FP = ifelse(grepl(decoy_string, Accession),1,0))  # Add an FP column which represents is a DECOY or hit, 1 for decoys
  df_FP_only <- transform(df_FP_only, Decoy = ifelse(grepl(decoy_string, Accession),"Decoy","Target")) # Add another column for Decoy or Target
  df_FP_only[,"Cum.FP"] <- cumsum(df_FP_only$FP) # Add cumulative column to the data frame
  df_FP_only$Hits.Above.Tresh <- seq(1,nrow(df_FP_only),1) # Add a column for the num of hits above a thresh
  df_FP_only <- transform(df_FP_only, TP = Hits.Above.Tresh - 2*Cum.FP)  # Add a column for TP which is = Hits.Above.Tresh - 2*Cum.FP
  df_FP_only <- transform(df_FP_only, FDR = Cum.FP / (TP  + Cum.FP)) # Add Column for FDR
  df_FP_only$Q.val <- return_Q_Val(df_FP_only$FDR,nrow(df_FP_only))  # Add a column for the Q value
  return(df_FP_only)
}

# -------------------------------------------------------------------

### Stats function where there exists an isDecoy Col

# Same as previous but uses the isDecoy column instead of regular expressions to identify decoys

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

### Function to calculate the column for Q value

# Cycles through taking the current position to the end selecting the
# min value from FDR.                                                   

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
  clev_df <- data.frame(peptides_to_process) #create the functions internal dataframe to handle 
  colnames(clev_df) <- c("Peptide")
  clev_df <- transform(clev_df, mod_p = gsub("^R|^K","A", Peptide)) #create a new column cutting the front off of any matching peptides
  clev_df <- transform(clev_df, mod_p_2 = gsub("R$|K$","A", mod_p)) #Create a new column cutting the end off of any matching peptides
  clev_df <- transform(clev_df, mod_p_3 = gsub("KP","A", mod_p_2)) #Create a col replacing where K followed by a P so not identified as a cleav
  clev_df <- transform(clev_df, pep_clev_count = str_count(mod_p_3,"R|K")) #Create a new col which has a count of the number of strings with a space seperating them
  return(clev_df$pep_clev_count)
}


# -------------------------------------------------------------------

### Function to get the current data Set

# Takes the file target and returns the first proccesed daataframe for use. It returns a list, first component is 
# the peptide list and second is the modication dataframe. This function autodetects file type and handles acordingly

get_current_dataSet <- function(in_File,passedUrlData,queryLen){
  
  # Selecting from the default or the uploaded
  if (is.null(in_File)&(length(queryLen)<1)){ # if both no file has been uploaded and no file passed by the url
    returnList <- list("pep" = NULL, "mod" = NULL) #Return data as a list without mods
    return(returnList) # give the default dataSet
    
  }else{
    #if no file has been uploaded and there is a passed URL file, then set in_File to the passedUrl and validate
    if(is.null(in_File)&(length(queryLen)>=1)){
      peptideDf <- get_url_dat(passedUrlData)
      returnList <- list("pep" = peptideDf, "mod" = getModCsv(peptideDf))
      return(returnList)
    }
    #check in File
    if(validate_file(in_File)){
      if(grepl("csv",in_File$datapath)){
        peptideDf <- handleFileCsv(in_File)
        returnList <- list("pep" = peptideDf, "mod" = getModCsv(peptideDf))
        return(returnList)
        
      }else if(grepl("mzid",in_File$datapath)){ #IF .mzid
        returnList <- list("pep" = handleFileMzid(in_File), "mod" = getMzidMod(in_File))
        return(returnList)
      }else{
        #na (gzip currently which is handled the same as mzid)
        returnList <- list("pep" = handleFileMzid(in_File), "mod" = getMzidMod(in_File))
        return(returnList)
      }
      
    }else{
      shinyalert(title = "Bad upload file!",text = "Please upload a file of: .csv .mzid .gz", type = "warning")
      
      returnList <- list("pep" = NULL, "mod" = NULL)
      return(returnList) # give the default dataSet
    }
  }
}


# -------------------------------------------------------------------

### Function to create a new dataframe for behinds the scenes which has the score set, has been ordered by score

# This creates the servers version of the data frame with a few exra cols, and also ordered by score as needed by the
# stats function for calcualted FDR correctly. It takes the column to set to score as an arguement. An original version
# of the primary dataSet is passed each time to this function so that mutliple columns aren't renamed score.

get_dataSet_withScore <- function(df_to_use,selected_col){
  try(df_to_use <- transform(df_to_use, z = round(Mass / m.z))) #Create a new column of Charge
  try(df_to_use <- transform(df_to_use, Mass = round(z * m.z))) #Create a new column of Charge
  try(df_to_use <- transform(df_to_use, ppm = round(ppm,digits = 4))) #Round the ppm col to 4dp
  colnames(df_to_use)[colnames(df_to_use)==selected_col] <- "Score" #set the score column
  df_to_use <- df_to_use[rev(order(df_to_use$Score)),] #order by score
  df_to_use$counted_cleavages <- as.factor(get_pep_cleav(df_to_use$Peptide)) # Count the cleavages and add to the data set as a new column
  return(df_to_use)
}


# -------------------------------------------------------------------

### Function to return min, max and a default of the score range

get_score_range <- function(df_to_use){
  rang <- c(min(df_to_use$pep$Score),max(df_to_use$pep$Score))
  return(rang)
  
}


# -------------------------------------------------------------------

### Function for mzid and gz file input handler:

# Takes an mzid target, imports, merges across the score and peptide list, then filters by max rank=1 and subsets
# for only unique peptide sequences - selecting the top scoring entry - where score is by default always the
# second column of the score df. Then renames some columns and calculates ppm and gets the cleav count

handleFileMzid <- function(inputFile){
  mzid <- openIDfile(inputFile$datapath)
  mzid_df <- merge(psms(mzid), score(mzid), by="spectrumID")  #merge by spectrumID
  mzid_df <- filter(mzid_df , rank == 1) #subset only rank 1
  
  #subet for unqiue peptides
  colnames(mzid_df )[colnames(mzid_df )==colnames(score(mzid)[2])] <- "ScoreFilterCol"
  mzid_df  <- mzid_df  %>% group_by(sequence) %>% slice(which.max(ScoreFilterCol))
  
  #Format the cols for the app
  colnames(mzid_df)[colnames(mzid_df)=="ScoreFilterCol"] <- colnames(score(mzid)[2])
  mzid_df <- getFormtedColDf(mzid_df)
  
  # Add extra columns
  mzid_df <- transform(mzid_df, ppm = ((m.z - calculatedMassToCharge)*1000000)/m.z) # Add a col for ppm
  
  return(mzid_df)
}


# -------------------------------------------------------------------

### Function to get the modification df from mzIdenML fies

getMzidMod <- function(inputFile){
  mods <- modifications(openIDfile(inputFile$datapath))
  return(mods)
}


# -------------------------------------------------------------------

### Function to return the CSV dataset

# Firstly the entire peptide column is set to uppercase for comparison (some inpput files were seen to have mixed).
# Tries to detect a recognised colname for filtering the unqiue peptides by max score. If not then tries to use the
# first column containing the string 'Score'. If no match found we do not filter - could select a peptide with the 
# low score

handleFileCsv <- function(in_File){
  
  uploaded_df <- read.csv(in_File$datapath, header=TRUE)
  
  #If the df does not contains more than 5 rows than we assume the wrong format and proceed to check for tab del with skip=0 & skip =1
  if(!NCOL(uploaded_df) >= 5){
    tryCatch({
        uploaded_df <- suppressWarnings(read_delim(in_File$datapath, "\t", #This is especially for MSAmanda
                                                   escape_double = FALSE, trim_ws = TRUE, 
                                                   skip = 0))
      },error = function(e) {stop(safeError(e))})# return a safeError if a parsing error occurs
    
    
    if(!NCOL(uploaded_df) >= 5){
      tryCatch({
        uploaded_df <- suppressWarnings(read_delim(in_File$datapath, "\t", 
                                                   escape_double = FALSE, trim_ws = TRUE, 
                                                   skip = 1))
      },error = function(e) {stop(safeError(e))})# return a safeError if a parsing error occurs
      
      
      if(!NCOL(uploaded_df) >= 5){
        shinyalert(title = "Warning!",text = "Could not read in CSV. Please ensure comma or tab delimited", type = "warning")
        return(NULL) #could read in CSV
      }
    }
  }
  
  uploaded_df <- getFormtedColDf(uploaded_df) #Rename cols
  uploaded_df$Peptide = toupper(uploaded_df$Peptide) #upper case the enitre col
  try(uploaded_df <- filter(uploaded_df , Rank == 1),silent = TRUE) #subset only rank 1
  try(uploaded_df <- filter(uploaded_df , rank == 1),silent = TRUE)
  uploaded_df <- getFilteredUniqueDfPep(uploaded_df) #filter for best scoring peptide
  return(uploaded_df)
}


# -------------------------------------------------------------------

### Fucntion to filter for unique peptides

# Filter for best scoring peptide sequence by score, if cannot get Score col then no filtering is done.

getFilteredUniqueDfPep <- function(df){
  if("X.10lgP" %in% colnames(df)){ #This is the general case for crowdsource data (the score col)
    df <- df %>% group_by(Peptide) %>% slice(which.max(X.10lgP))
  }else if(grep("Score", colnames(df))){
    ## R makes this problamatic as it's difficutl to select a colname by a variable
    # Soltuion: set variable column to fixed col name then return it to original name later
    colnames <- colnames(df)
    col_mat <- grepl("Score",colnames)
    selected_default <- colnames[which(col_mat)]
    colnames(df)[colnames(df)== selected_default[1] ] <- "ScoreFilterCol"
    df <- df %>% group_by(Peptide) %>% slice(which.max(ScoreFilterCol))
    colnames(df)[colnames(df)=="ScoreFilterCol"] <- selected_default[1]
    
  }else{
    #no match found, so rather than filter at random we choose to leave
  }
  return(df)
}


# -------------------------------------------------------------------

### Function to read in data from a web URL

get_url_dat <- function(webAddress){
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Reading in Data", value = 0.99)
  webDat <- read.csv(url(webAddress))
  webDat$counted_cleavages <- as.factor(get_pep_cleav(webDat$Peptide))
  return(webDat)
}


# -------------------------------------------------------------------

### Function to check the file uploaded is valid

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

### function to check if score the score column is a numeric

checkScoreColNum <- function(df_to_use,selected_col){
  colnames(df_to_use)[colnames(df_to_use)==selected_col] <- "Check" #set the score column
  return(!is.numeric(df_to_use$Check))
}


# -------------------------------------------------------------------

### Function to return the url correct form

returnDataUrl <- function(string){
  string <- sub("[?]","",string)
  string <- sub("id=","",string)
  return(paste("http://pgb.liv.ac.uk/~andrew/crowdsource-server/src/public_html/results",string,"psm.csv",sep = "/"))
}


# -------------------------------------------------------------------

### Function to return the url for Server Data csv

returnServerDataCsvUrl <- function(string){
  string <- sub("[?]","",string)
  string <- sub("id=","",string)
  return(paste("http://pgb.liv.ac.uk/~andrew/crowdsource-server/src/public_html/results",string,"locations.csv",sep = "/"))
}


# -------------------------------------------------------------------

### Functions to check the columns required exist and if not return an error message instead of plot

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

### Function to return the current server side version of the dataset

# cbinds the server dataframe with a calculated stats df for plots like FDR and Q curve.
# function of choice depends whether the column isDecoy exists.

returnCurrentServerDF <- function(current_dataSet,col,decoyString){
  returnPepDf <- get_dataSet_withScore(current_dataSet()$pep,col) #Set the score col
  if("isDecoy" %in% colnames(returnPepDf)){ #get the stats cols
    statsDf <-getStatsDfExistingDecoy(returnPepDf$isDecoy)
  }else{
    statsDf <- get_df_stat(returnPepDf$Accession,decoyString)
  }
  returnPepDf <- cbind.data.frame(returnPepDf,statsDf) #Combine pep and stats df
  returnList <- list("pep" = returnPepDf, "mod" = current_dataSet()$mod) #return all
  return(returnList)
  
}


# -------------------------------------------------------------------

### Function to return numeric colnames as a list without the stats DF elements

getColNames <- function(df){
  colNames <- names(Filter(is.numeric,df)) # Get the data set with the appropriate name
  colNames <- colNames[1:(length(colNames)-6)] # remove the last 6 cols which are ffrom the stats calc
  return(colNames)
}


# -------------------------------------------------------------------

### Function for v&h line on the fdr curve

vline <- function(x = 0, color = "grey") {
  list(type = "line", y0 = 0, y1 = 1, yref = "paper",x0 = x, x1 = x, line = list(color = color))
}

hline <- function(y = 0, color = "grey") {
  list(type = "line", x0 = 0, x1 = 1, xref = "paper",y0 = y, y1 = y, line = list(color = color))
}


# -------------------------------------------------------------------

### Function to get intercept

# Approximates the point where FDR percent meets the FDR curve, takes two arguemnts, the df containing the FDR curve
# and the input of FDR percentage

getIntercept <- function(df,fdr){
  intercepDF <- subset(df, Q.val > 0) #no zero interpolation
  interCep <- approx(x = as.numeric(intercepDF$Q.val), y = as.numeric(intercepDF$TP), xout = fdr/100)
  return(round(interCep$y))
}


# -------------------------------------------------------------------

### Function to return a dataframe with counts of mods, only pass it a dataframe of ID and name

# This counts how many mods are present per type of modification and returns a df for plotting a bar chart

getModCount <- function(modificationDf,numOfPep){
  colnames(modificationDf) <- c("spectrumID", "name")
  uniqueMods <- unique(modificationDf[ , 1:2 ] ) #unique 
  pepIdWithAnyMod <- data.frame(unique(modificationDf[ , 1] )) #peptide ids with any modification
  pepWithoutMod <- numOfPep-nrow(pepIdWithAnyMod) # num of peptides with no mods
  returnDf <- rbind.data.frame(data.frame("Var1" = "no mods", "Freq" = pepWithoutMod),as.data.frame(table(unlist(uniqueMods$name))))
  colnames(returnDf) <- c("Modification", "Frequency")
  return(returnDf)
}


# -------------------------------------------------------------------

### Function to return a modification dataframe for the non .mzid files (.csv) enables to treat it the same for ptm plots

# Takes the modification column and spectrum ID then creates a new df where each individual mod in the PTM col
# has its own entry with the relevant spectrum ID to be used as a key. Note: this is for processing csv into the same
# format as mzIdentML mod dataframes

getModCsv <- function(peptideDf){
  ptmList <- as.data.frame(peptideDf$PTM) #First step is to isolate the ptm column from the peptide list
  ptmList$spectrumID <- seq.int(nrow(ptmList)) #add an ID col
  ptmList <- na.omit(ptmList) #remove na rows
  colnames(ptmList) <- c("name", "spectrumID")
  mods <- cSplit(ptmList, "name", sep = ";", direction = "long") #split so each mod is its own row
  return(mods)
}


# -------------------------------------------------------------------

### Function to return a dataframe with the count of modifications per peptide

getModCountPerPep <- function(modificationDf,numOfPep){
  colnames(modificationDf) <- c("spectrumID", "name")
  
  pepIdWithAnyMod <- data.frame(unique(modificationDf$spectrumID))
  pepWithoutMod <- numOfPep-nrow(pepIdWithAnyMod)
  returnDf <- data.frame("Var1" = 0, "Freq" = pepWithoutMod)
  
  countPerPep <- as.data.frame(table(unlist(modificationDf$spectrumID)))
  countPerNum <- as.data.frame(table(unlist(countPerPep$Freq)))
  
  returnDf <- rbind.data.frame(returnDf,countPerNum)
  colnames(returnDf) <- c("Modifications", "Frequency")
  return(returnDf)
}


# -------------------------------------------------------------------

### Function to take a dataframe (input file) and rename the columns in a correct format

# It became problatic to take files with differeing column names so here is a function to sweep through ranges
# of potential column names and format that in an acceptable way

getFormtedColDf <- function(df){
  
  #General .mzid names
  try(colnames(df)[colnames(df)=="DatabaseAccess"] <- "Accession")
  try(colnames(df)[colnames(df)=="sequence"] <- "Peptide")
  try(colnames(df)[colnames(df)=="experimentalMassToCharge"] <- "m.z")
  try(colnames(df)[colnames(df)=="chargeState"] <- "z")
  
  #MS Amanda renames
  try(colnames(df)[colnames(df)== "Sequence" ] <- "Peptide")
  try(colnames(df)[colnames(df)== "Protein Accessions" ] <- "Accession")
  try(colnames(df)[colnames(df)== "Modifications" ] <- "PTM")
  try(colnames(df)[colnames(df)== "Charge" ] <- "z")
  
  return(df)
}


# -------------------------------------------------------------------

### Function to return column

# Checks the columns for which should be chosen as the intial scoring column
# works from known commonly used score columns present in a range of file uploads

getInitScoreCol <- function(df){
  colnames <- names(df) # Get the data set with the appropriate name
  
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
  }else if(grep("Score", colnames)){
    col_mat <- grepl("Score",colnames)
    selected_default_vec <- colnames[which(col_mat)]
    selected_default <- selected_default_vec[1]
  }else{
    selected_default <- NULL
  }
  
  return(selected_default)
}








