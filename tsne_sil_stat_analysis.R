#rm(list=ls()) # to clean the environment

#LIBRARIES TO BE INSTALLED IN CASE ERRORS APPEAR FOR PARALLELISED RUNNING
# install.packages("doSNOW")
# install.packages("doParallel") 
# install.packages("doMPI")


libdir='../Rpackages'

.libPaths( c( libdir , .libPaths() ) )


#NECESSARY LIBRARIES FOR THIS SCRIPT
library(Rtsne)        #used for t-SNE algorithm
library(factoextra)   # used for clustering algorithms
library(doParallel)   # parallelised scripting
library(foreach)      # for each parallelised loops



#---------------------RUNNING SCRIPT --------------------------------------

#set current folder as your working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#getwd()

source("functions.R")

region_analysed = "Bath and North East Somerset"

# create a direction for saving RDS files
direction <- sprintf("./tsne_sil_stat_analysis_%s", region_analysed)
dir.create(file.path(direction), showWarnings = FALSE) 


df_selection <- region_selector(region_analysed)


#defining range for number of clusters
range_def <- seq(6, 10, by=1)

#define range for numbe of iterations by specifying different seeds
seed_range <- seq(1,1000,by=1)

#clustering methods definition
clustm <- c("kmeans", "pam","agnes", "clara", "diana")


# Setting up parallelised script running
numCores <- detectCores() # detects how many cores the system has

cl <- makeCluster(numCores[1]-1) #not to overload the computer
registerDoParallel(cl)

system.time({
for (clust_no in range_def){
df_new <- data.frame()
for (method in clustm){  
df_tsne_sil_stat <-   foreach (i = seed_range, .combine = rbind) %dopar% { 
  
                      .libPaths( c( libdir , .libPaths() ) )
                      library(Rtsne)   
                      library(factoextra) 
  
                      tsne_df <- tsne_transform(df_selection, region_analysed, perplexity = 10, max_iter = 10000, seed = i,
                                                full_name = FALSE)
                      clust.res <- eclust(tsne_df, method, k=clust_no, graph = FALSE)
                      clust.sil <- fviz_silhouette(clust.res, ggtheme = theme_minimal())
                      sil_title <- clust.sil$labels$title
                      sil_coeff <- as.numeric(unlist(regmatches(sil_title,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",sil_title))))
                      df <- data.frame(method, clust_no, sil_coeff)
                      return(df)
                    }

df_new <- rbind(df_new,df_tsne_sil_stat)
}
  df_new$method <- as.character(df_new$method)
  df_new$clust_no <-as.numeric(as.character(df_new$clust_no))
  df_new$sil_coeff <-as.double(as.character(df_new$sil_coeff))
saveRDS(df_new,file = sprintf("./tsne_sil_stat_analysis_%s/tsne sil_coeff - %.0f clusters.rds", region_analysed, clust_no))
}
})

stopImplicitCluster() # clean up the cluster


# CHECK WHAT HAS BEEN WRITTEN TO THE RDS FILES
# library(rlist)
# list1 <- list()
# for (clust_no in range_def){
#   list1 <- list.append(list1,readRDS(file = sprintf("./tsne_sil_stat_analysis_%s/tsne sil_coeff - %.0f clusters.rds", region_analysed, clust_no)))
# }
# list1



