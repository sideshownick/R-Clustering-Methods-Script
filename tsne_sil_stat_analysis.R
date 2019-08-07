rm(list=ls()) # to clean the environment

#install.packages("foreach")
#install.packages("doParallel")


#LIBRARIES USED IN THIS SCRIPT
library(ggplot2)      #used for ggplot
library(ggpubr)       #used for ggarrange
library(factoextra)   #used for eclust function
require(gridExtra)    #used for grid.arrange
library(Rtsne)        #used for t-SNE algorithm
library(RColorBrewer) #used for color palettes for ggplot
library(ggthemr)      #theme for plotting in ggplot
library(ggthemes)     #theme for plotting in ggplot
library(plyr)         #used to calc mean for df
library(rlist)        #for appending lists
library(parallel) #detects number of cores in the system

library(doParallel)
library(foreach)


#TESTING GITHUB FUNCTIONALITY
#----@@@@@@@@@@@----RUNNING SCRIPT ----@@@@@@@@@@@@@@@@@@@@

#set current folder as your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

source("functions.R")

region_analysed = "Bath and North East Somerset"

# create a direction for saving RDS files
direction <- sprintf("./tsne_sil_stat_analysis_%s", region_analysed)
dir.create(file.path(direction), showWarnings = FALSE) 

df_selection <- region_selector(region_analysed)


#defining range for number of clusters
range_def <- seq(6, 8, by=1)

#defining the number of iterations to be performed
no_iterations = 4

numCores <- detectCores() # detects how many cores the system has
registerDoParallel(numCores)  # use multicore, set to the number of our cores

# Return a vector
system.time({
foreach (i=range_def) %dopar% {
  df_mc <- silhouettes1(df_selection,no_iterations,i)
  saveRDS(df_mc,file = sprintf("./tsne_sil_stat_analysis_%s/%.0f_cl - tsne sil_coeff.rds", region_analysed,i))
}
})

