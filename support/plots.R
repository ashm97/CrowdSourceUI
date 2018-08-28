##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Plot functions
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

# -------------------------------------------------------------------

## Function to plot scatter with varibles

plotVariableScatter <- function(current_dataSet_server_side,score_range,seperateByDecoy,x_col_id,y_col_id,axisScale){

  #Subset for the score slider range
  df_to_plot <- subset(current_dataSet_server_side$pep, Score > score_range[1] & Score < score_range[2])
  
  if(seperateByDecoy){  #if seperate by decoy is true
    
    # dummy up data
    dat1<- subset(df_to_plot,grepl("Decoy",Decoy)) #Decoy entries
    dat2<- subset(df_to_plot,!grepl("Decoy",Decoy))  #Target entries
    
    gg_to_plot <- ggplot() 
    gg_to_plot <- gg_to_plot + geom_point(data=dat2, aes_string(x=x_col_id, y=y_col_id, color="Score"), shape=22, size=3)
    gg_to_plot <- gg_to_plot + scale_color_gradient(high="plum1", low="darkmagenta",name="Target Score")
    gg_to_plot <- gg_to_plot +  geom_point(data= dat1, aes_string(x=x_col_id, y=y_col_id, shape="shp", fill="Score"), shape=21, size=2,alpha=0.5)
    gg_to_plot <- gg_to_plot + scale_fill_gradient(high="darkgreen", low="aquamarine",name="Decoy Score") + theme_bw()
    
    if("x" %in% axisScale){
      gg_to_plot <- gg_to_plot + scale_x_continuous(trans='log2') 
    }
    
    if("y" %in% axisScale){
      gg_to_plot <- gg_to_plot + scale_y_continuous(trans='log2')
    }
    
    return(gg_to_plot)
    
  }else{  #if seperate by deocy was false
    
    gg_to_plot <- ggplot(df_to_plot, aes_string(x_col_id,y_col_id, colour = "Score"))+ 
      geom_count(show.legend=T,alpha = 0.5)+
      theme_bw()+
      scale_colour_gradient2(name = "Score", low = "yellow", mid = "aquamarine",
                             high = "darkmagenta")+
      labs(size="Point Density")
    
    if("x" %in% axisScale){
      gg_to_plot <- gg_to_plot + scale_x_continuous(trans='log2') 
    }
    
    if("y" %in% axisScale){
      gg_to_plot <- gg_to_plot + scale_y_continuous(trans='log2')
    }
    
    return(gg_to_plot)
    
  }
  
  
  
}




# -------------------------------------------------------------------

## Function to produce generic hist plot with inputs

plot_hist <- function(current_dataSet_server_side,no_bins,title_of_plot,yLabel,col_id){
  df_to_plot <- data.frame(current_dataSet_server_side$pep)
  if(col_id %in% colnames(df_to_plot)){
    gg_to_plot <- ggplot(df_to_plot,aes_string(x=col_id)) +
      geom_histogram(fill="blue", col="blue", alpha = .2, bins = no_bins) +
      #geom_density(aes(y=..count..), colour="red", adjust=4) +
      theme_bw() +
      labs(title=title_of_plot, x=yLabel, y="Count")
    
    ggplotly(gg_to_plot, tooltip = df_to_plot$ppm)
    
  }else{
    validate(
      need(FALSE, paste("Column: ",yLabel," does not exist"))
    )
  }
  
}

# -------------------------------------------------------------------

## Function to produce generic bar plot without inputs

plot_bar <- function(current_dataSet_server_side,title_of_plot,yLabel,col_id){
  df_to_plot <- data.frame(current_dataSet_server_side$pep)
  #Checking the col name exist
  if(col_id %in% names(df_to_plot)){
    gg_to_plot <- ggplot(data=df_to_plot, aes_string(x=col_id)) +
      geom_bar(fill="blue", col="blue", alpha = .2)+
      theme_bw() +
      labs(title=title_of_plot, x=yLabel, y="Count")
    
    return(gg_to_plot)
  }else{
    validate(
      need(FALSE, paste("Column: ",yLabel," does not exist"))
    )
  }
}

# -------------------------------------------------------------------

## Function to plot bar of cleavages

plot_bar_cleav <- function(current_dataSet_server_side,decoyToggle){
  #if decoy decoy = TRUE show all entries 
  if(decoyToggle){
    df_to_plot <- data.frame(current_dataSet_server_side$pep)
  }else{  #subset the DB for only target entries
    df_to_plot<- subset(current_dataSet_server_side$pep,!grepl("Decoy",Decoy))  #Target entries
  }
  checkColExist(df_to_plot$counted_cleavages, "Missing Column: Cleavages Count")
  #rename for the plot hover over
  colnames(df_to_plot)[colnames(df_to_plot)=="counted_cleavages"] <- "Cleavages"
  
  gg_to_plot <- gg_to_plot <- ggplot(data=df_to_plot, aes(Cleavages)) +
    geom_bar(fill="blue", col="blue", alpha = .2)+
    theme_bw() +
    labs(x="Cleavages", y="Count")
  return(gg_to_plot)
  
}


# -------------------------------------------------------------------

## Function to plot the box plot of cleavages by score

plot_box_clev_by_score <- function(current_dataSet_server_side,decoyToggle){
  #if decoy decoy = TRUE show all entries 
  if(decoyToggle){
    df_to_plot <- data.frame(current_dataSet_server_side$pep)
  }else{  #subset the DB for only target entries
    df_to_plot<- subset(current_dataSet_server_side$pep,!grepl("Decoy",Decoy))  #Target entries
  }
  
  checkColExist(df_to_plot$counted_cleavages, "Missing Column: Cleavages Count")
  plot_ly(df_to_plot, y = ~Score, color = ~counted_cleavages, type = "box")%>%
    layout(xaxis = list(title = "Num Cleavages"),
           yaxis = list(title = "Score"))
}


# -------------------------------------------------------------------

## Function to plot a barplot of PTMs

plot_bar_ptm <- function(current_dataSet_server_side){
  
  checkColExist(current_dataSet_server_side$mod$spectrumID, "Missing Mod Data")
  checkColExist(current_dataSet_server_side$mod$name, "Missing Mod Data")
  
  barDf <- getModCount(data.frame(current_dataSet_server_side$mod$spectrumID,current_dataSet_server_side$mod$name),nrow(current_dataSet_server_side$pep))
  
  
  
  gg_to_plot <- ggplot(data=barDf, aes(Modification,Frequency)) +
    geom_col(fill="blue", col="blue", alpha = .2,width = 0.5)+
    theme_bw() +
    labs(title = "Peptides with Modification ",x=" ", y="Count")
  return(gg_to_plot)
  
  
}


# -------------------------------------------------------------------

## Funcition to plot bar of ptm modification count

plot_bar_ptm_mod_count <- function(current_dataSet_server_side){
  checkColExist(current_dataSet_server_side$mod$spectrumID, "Missing Mod Data")
  checkColExist(current_dataSet_server_side$mod$name, "Missing Mod Data")
  
  plotDf <- getModCountPerPep(data.frame(current_dataSet_server_side$mod$spectrumID,current_dataSet_server_side$mod$name),nrow(current_dataSet_server_side$pep))
  
  gg_to_plot <- ggplot(data=plotDf[1:10,], aes(Modifications,Frequency)) +
    geom_col(fill="blue", col="blue", alpha = .2)+
    theme_bw() +
    labs(title = "Number of Mods per Peptide ",x=" ", y="Modifications")
  return(gg_to_plot)
  
}


# -------------------------------------------------------------------

## Funcition to plot scatter of score and ppm, coloured by decoy

plot_scat_score_ppm_by_decoy <- function(current_dataSet_server_side,pointOpacity,pointsToDisplay){
  checkColExist(current_dataSet_server_side$pep$ppm, "Missing Column: ppm")
  
  df_to_plot <- data.frame(current_dataSet_server_side$pep$ppm,current_dataSet_server_side$pep$Score,
                           current_dataSet_server_side$pep$Decoy)
  
  colnames(df_to_plot) <- c("ppm","Score","Decoy")
  
  if(pointsToDisplay == 2){ #Display only Targets
    df_to_plot <- subset(df_to_plot,grepl("Target",Decoy))
    
    plot_ly(
      df_to_plot, x = ~ppm, y = ~Score, type = 'scatter',
      color = ~Decoy, colors = "mediumturquoise", symbol = ~Decoy, symbols = 'o',
      marker=list( size=10 , opacity=pointOpacity)
      
    )
    
    
  }else if (pointsToDisplay == 3){  #Display only Decoys
    df_to_plot <- subset(df_to_plot,grepl("Decoy",Decoy))
    
    plot_ly(
      df_to_plot, x = ~ppm, y = ~Score, type = 'scatter',
      color = ~Decoy, colors = 'indianred', symbol = ~Decoy, symbols ='x',
      marker=list( size=10 , opacity=pointOpacity)
      
    )
    
  }else{  #Display all points
    plot_ly(
      df_to_plot, x = ~ppm, y = ~Score, type = 'scatter',
      color = ~Decoy,colors = c('indianred',"mediumturquoise"), symbol = ~Decoy, symbols = c('x','o'),
      marker=list( size=10 , opacity=pointOpacity)
      
    )
  }
  
}


# -------------------------------------------------------------------

## Function to plot Box plot for Marginal Density of ppm

plot_box_marg_ppm <- function(current_dataSet_server_side){
  checkColExist(current_dataSet_server_side$pep$ppm, "Missing Column: ppm")
  df_to_plot <- data.frame(current_dataSet_server_side$pep$ppm,current_dataSet_server_side$pep$Score,
                           current_dataSet_server_side$pep$Decoy)
  colnames(df_to_plot) <- c("ppm","Score","Decoy")
  
  plot_ly(type = 'box') %>%
    layout(xaxis = list(showgrid = F),
           yaxis = list(showgrid = F))%>%
    add_boxplot(y = df_to_plot$Score[df_to_plot$Decoy == "Decoy"], name = "Decoy", boxpoints = FALSE,
                marker = list(color = 'indianred'),
                line = list(color = 'indianred')) %>%
    add_boxplot(y = df_to_plot$Score[df_to_plot$Decoy == "Target"], name = "Target", boxpoints = FALSE,
                marker = list(color = 'mediumturquoise'),
                line = list(color = 'mediumturquoise')) %>%
    layout(yaxis = list(title = "Score"))
  
  
}


# -------------------------------------------------------------------

## Function to plot Box plot for Marginal Density of Score (-10lgP)

plot_box_marg_score <- function(current_dataSet_server_side){
  checkColExist(current_dataSet_server_side$pep$ppm, "Missing Column: ppm")
  df_to_plot <- data.frame(current_dataSet_server_side$pep$ppm,current_dataSet_server_side$pep$Score,
                           current_dataSet_server_side$pep$Decoy)
  colnames(df_to_plot) <- c("ppm","Score","Decoy")
  
  plot_ly(type = 'box') %>%
    layout(yaxis = list(showgrid = F),
           xaxis = list(showgrid = F))%>%
    add_boxplot(x = df_to_plot$ppm[df_to_plot$Decoy == "Decoy"], name = "Decoy", boxpoints = FALSE,
                marker = list(color = 'indianred'),
                line = list(color = 'indianred')) %>%
    add_boxplot(x = df_to_plot$ppm[df_to_plot$Decoy == "Target"], name = "Target", boxpoints = FALSE,
                marker = list(color = 'mediumturquoise'),
                line = list(color = 'mediumturquoise')) %>%
    layout(xaxis = list(title = "Mass Error"))
  
}


# -------------------------------------------------------------------

## Function plot histogram of score, with deocy vs target

plot_hist_score <- function(current_dataSet_server_side,fdrPercent){
  
  checkColExist(current_dataSet_server_side$pep$Decoy, "Missing Column: Decoy")
  df_to_plot <- data.frame(current_dataSet_server_side$pep$Decoy,current_dataSet_server_side$pep$Score)
  names(df_to_plot)[names(df_to_plot)=="current_dataSet_server_side.pep.Decoy"] <- "Type"
  names(df_to_plot)[names(df_to_plot)=="current_dataSet_server_side.pep.Score"] <- "Score"
  
  xInterceptLine <- current_dataSet_server_side$pep$Score[getIntercept(current_dataSet_server_side$pep,fdrPercent)]
  
  plot <- ggplot(df_to_plot, aes(Score)) +
    geom_histogram(aes(fill = Type),binwidth=2,col="white") +
    theme_bw() +
    labs(title="Distribution of peptides", y="Count", x="Score") +
    scale_fill_manual(name="Peptide", 
                      values = c("Target"='mediumturquoise', 
                                 "Decoy"='indianred'))+ 
    geom_vline(xintercept = xInterceptLine, linetype="dotted")
  
  return(plot)
}


# -------------------------------------------------------------------

## Function to plot FDR curve

plot_FDR_curve <- function(current_dataSet_server_side,fdrPercent){
  validate(
    need(!is.null(current_dataSet_server_side$pep), "No data set uploaded ")
  )
  
  
  df <- current_dataSet_server_side$pep
  
  plot_ly(df, x = ~TP, y = ~FDR, mode = 'lines',line = list(color = 'rgb(205, 12, 24)', width = 4)) %>%
    layout(shapes = list(vline(getIntercept(df,fdrPercent)),hline(fdrPercent/100))) %>% 
    layout(xaxis = list(title = "Num of peptide spectrum matchings"),
           yaxis = list(range = c(0, .05)))
  
}


# -------------------------------------------------------------------

## Function to plot Q curve

plot_Q_curve <- function(current_dataSet_server_side){
  validate(
    need(!is.null(current_dataSet_server_side$pep), "No data set uploaded ")
  )
  df <- current_dataSet_server_side$pep
  plot_ly(df, x = ~Score, y = ~Q.val, type = 'scatter', mode = 'lines',line = list(color = 'rgb(205, 12, 24)', width = 4)) %>% 
    layout(xaxis = list(title = "Score"),
           yaxis = list(range = c(0, .05)))
  
}


# -------------------------------------------------------------------

## Plotly Pie chart of Identified at Set FDR

plotIdentFdr <- function(current_dataSet_server_side,fdrPercent){
  
  validate(
    need(!is.null(current_dataSet_server_side$pep), "No data set uploaded ")
  )
  
  count <- getIntercept(current_dataSet_server_side$pep,fdrPercent)
  
  df <- data.frame(hits=c("below", "above"), 
                   columnofvalues=c(count,nrow(current_dataSet_server_side$pep)-count))
  
  colors <- c('mediumturquoise',  'indianred')
  
  plot_ly(df, labels = ~hits, values = ~columnofvalues, type = 'pie',
          textposition = 'inside',
          textinfo = 'label+percent',
          insidetextfont = list(color = '#FFFFFF'),
          hoverinfo = 'text',
          text = ~paste('Hits: ', columnofvalues),
          marker = list(colors = colors,
                        line = list(color = '#FFFFFF', width = 1)),
          showlegend = FALSE)
  
  
}

# -------------------------------------------------------------------{}