####################################################################################################
##
## Project: Omics Shiny Search Results Application
## Script purpose: Description of Project Code
## Date: 03.09.2018
## Author: Ashleigh Myall
##
####################################################################################################

#---- Description ----#

The application is for a web based, interactive visualisation of the results for peptide search engines 
in the shotgun proteomics pipeline, focusing compatibility with the crowdsourcing distributed engine output 
and mzIdentML file formats.


#---- Installation ----#

For offline usage unpack all the files into a directory. Open R studio, set the working directory and source 
from Global.R after installing the packages called inside the Global.R script.


#---- Usage ----#

For general use navigate to the Data page first. Upload a file which meets the format requirements (.csv .mzid .gz) 
containing peptide search results (Note: Include a column labeled sequence or peptide). You may view the 
uploaded peptide dataframe in the table on the same page. From this page also the decoy term can be set, 
which is taken from the accession column if no isDecoy column is already present. A drop down menu for score 
column selection also allows the user to specify what to be used as a scoring metric (Note: for mzIdentML files 
the whole scoring dataframe is binded to the peptide list).

The application also allows for passing the URL of a mzid or csv. This has been developed to specifically interface 
with the results from the crowdsource. And such is takes /?id=num. Where num is the number in the URL of the 
directory result from our crowdsource developers file URL. This URL also directs the application to a location 
for server data. Which is used for plotting points on the leaflet map and user statistics calculations on the 
homepage. For zn mzid file to be passed, we save the file as temp in th dat file while the application reads it in as mzR
doesnt support URL as a source dest

Visualizations are done in a mixture of shiny and ggplot. Many offer interactivity such as zooming and mouse 
over information. Also query creators are placed next to some plots which can be used to change displays.


#---- Structure  ----#

Split across multiple scripts the code accesses functions which has been located by category. The application makes
use of module in a large part. All modules are placed in the module script (both UI and server Modules). Each 
page has a UI module and Server module. Which within can call other inner modules for other displays like plots.

The UI is split across multiple scripts. With UI.R at the highest level holding the tabs to display. Then body.R 
which UI.R access, then the tab.R script, which calls the UI module for each tab.

Much of the script revolves around the data page module which returns the data in the global environment for the 
server to access and display plots for.



#---- Credits ----#

Ashleigh Myall - R coding 
Andy Jones - Project Supervisor
Andrew Collins - crowdsource search engine software developer


####################################################################################################