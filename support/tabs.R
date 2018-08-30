##################################################
## Project: Omics Shiny Search Results Application
## Script purpose: Script managing each display tab
## Date: 06.08.2018
## Author: Ashleigh Myall
##################################################

# Main tab
main_tab <- tabItem(tabNam = "main",
                    homePageInput("mainPage")
                    )

# Second tab 'Data' content
data_tab <- tabItem(tabName = "data",
                    dataPageInput("datP")
                    
)

# Scatter tab
scatter_tab <- tabItem(tabName = "scatter",
                       singleScatPageInput("scat")
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
                    box(width = 8, title = "Background",
                        h4("The crowdsourcing search engine project has been to develop a distributed computing network that performs peptide identification by scoring MS/MS spectra against peptides derived from a protein sequence database."),
                        h4("Developed by Andy Jones and Andrew Collins the project has hoped to share the load and compute time of large searches amongst many nodes, while maintaining similar performance to other search engines (for example: MS Amanda and MS-GF+)."),
                        h4("Created to aid the search engine, this shiny web interface allows for easy visualising of the results, with a high degree of interactivity for data exploration. An additional feature beyond the crowdsourcing project requirements, has been the ability to upload results from a wide range of search engines. Compatibility has been a focus; both mzIdentML as well as CSV peptide search results can be uploaded and interpreted.")
                        ),
                    box(width=4,title = "Profiles")
                    )