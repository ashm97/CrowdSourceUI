##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Script managing each display tab
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################


# Second tab 'Data' content
data_tab <- tabItem(tabName = "data",
                    dataPageInput("datP")
                    
)


# Scatter tab
scatters_tab <- tabItem(tabName = "scatters",
                        scatDisInput("scats")
)

# Hist tab
histograms_tab <- tabItem(tabName = "histograms",
                          histInput("hist")
                          
)

cleavages_tab <- tabItem(tabName = "cleavages",
                         cleavPageInput("cleav")
)

ptm_tab <- tabItem(tabName = "ptm",
                   ptmPageInput("ptmPage")
)

decoy_tab <- tabItem(tabName = "decoy",
                     decoyDisInput("dec")
)

score_tab <- tabItem(tabName = "score",
                     scorePageDisplayInput("scoreDis")
)

info_tab <- tabItem(tabName = "info",
                    h2("Information Page"),
                    box(width = 12, h4("Placeholder text is the label for possible content in a text box. It can normally be found when there are prompts to fill out a form. It's the hint that tells you 'Last name' or the format with which to enter your birthdate or phone number. Placeholder text typically exists as a hint to fill in actual text."))
                    )