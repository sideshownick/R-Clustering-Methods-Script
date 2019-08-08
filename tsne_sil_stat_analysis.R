#rm(list=ls()) # to clean the environment

# install.packages("doSNOW")
# install.packages("doParallel") 
# install.packages("doMPI")

library(Rtsne)        #used for t-SNE algorithm
library(factoextra)   # used for clustering algorithms
library(doParallel)   # parallelised scripting
library(foreach)      # for each parallelised loops


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
range_def <- seq(6, 10, by=1)

#define range for numbe of iterations by specifying different seeds
seed_range <- seq(1,10,by=1)

#clustering methods definition
clustm <- c("kmeans", "pam","agnes", "clara", "diana")


# Setting up parallelised script running
numCores <- detectCores() # detects how many cores the system has

cl <- makeCluster(numCores[1]-1) #not to overload the computer
registerDoParallel(cl)



system.time({
df_tsne_sil_stat <- foreach (clust_no = range_def, .combine = rbind)%:% #no of clusters
                    foreach(method = clustm, .combine = rbind) %:% #5 methods
                    foreach (i = seed_range, .combine = rbind, .packages=c('Rtsne','factoextra')) %dopar% { #a lot of iterations
                      tsne_df <- tsne_transform(df_selection, region_analysed, perplexity = 10, max_iter = 10000, seed = i,
                                                full_name = FALSE)
                      clust.res <- eclust(tsne_df, method, k=clust_no, graph = FALSE)
                      clust.sil <- fviz_silhouette(clust.res, ggtheme = theme_minimal())
                      sil_title <- clust.sil$labels$title
                      sil_coeff <- as.numeric(unlist(regmatches(sil_title,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",sil_title))))
                      df_pca <- data.frame(method, clust_no, sil_coeff)
                      return(df_pca)
                    }
})
stopImplicitCluster() # clean up the cluster

df_tsne_sil_stat$method <- as.character(df_tsne_sil_stat$method)
df_tsne_sil_stat$clust_no <-as.numeric(as.character(df_tsne_sil_stat$clust_no))
df_tsne_sil_stat$sil_coeff <-as.double(as.character(df_tsne_sil_stat$sil_coeff))

saveRDS(df_tsne_sil_stat,file = sprintf("./tsne_sil_stat_analysis_%s/tsne sil_coeff.rds", region_analysed))






