#rm(list=ls()) # to clean the environment
libdir='./Rpackages'
.libPaths( c( libdir , .libPaths() ) )

#r <- getOption("repos")
#r["CRAN"] <- "http://www.stats.bris.ac.uk/R/"
#options(repos = r)
#LIBRARIES TO BE INSTALLED IN CASE ERRORS APPEAR FOR PARALLELISED RUNNIN# 
#install.packages("doSNOW")
#install.packages("doParallel")
#install.packages("doMPI")


#NECESSARY LIBRARIES FOR THIS SCRIPT
library(Rtsne)        #used for t-SNE algorithm
library(factoextra)   # used for clustering algorithms
library(doParallel)   # parallelised scripting
library(foreach)      # for each parallelised loops


#---------------------RUNNING SCRIPT --------------------------------------

#set current folder as your working directory * only for RStudio
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#getwd()

source("functions.R")

region_analysed = "Bath and North East Somerset"

# create a direction for saving RDS files
#direction <- sprintf("./tsne_sil_stat_analysis_%s", region_analysed)
#dir.create(file.path(direction), showWarnings = FALSE) 


# Setting up parallelised script running
numCores <- 16 #set manually (10) #detectCores() # detects how many cores the system has

#cl <- makeCluster(cores=numCores, outfile="") #numCores[1]/2) #not to overload the computer
#registerDoParallel(cl) #(cl)
registerDoParallel(cores=numCores)


df_selection <- region_selector(region_analysed)

#defining range for number of clusters
range_def <- seq(6, 10, by=1)

#define range for numbe of iterations by specifying different seeds
N <- 1024 #total number of seeds (power of two to divide neatly)
num <- numCores #number by which to subdivide jobs [2,4,8,16]
seed_range <- seq(1,N,by=1)
outer_seed_range <- seq(1,N/num,by=1) #number of batches to send to cores
inner_seed_range <- seq(1,num,by=1) #jobs in cores' memory at one time

#clustering methods definition
clustm <- c("kmeans", "pam","agnes", "clara", "diana")


dtime = system.time({
for (clust_no in range_def){
    df_new <- data.frame()
    for (method in clustm){
        for (j in outer_seed_range) { #batches of parallel jobs
            df_tsne_sil_stat <- foreach (i = inner_seed_range, .combine = rbind) %dopar% 
            {     #parallel jobs
                  tsne_df <- tsne_transform(df_selection, 
                    region_analysed, perplexity = 10, max_iter = 10000, 
                    seed = j*num + i, full_name = FALSE)
                  clust.res <- eclust(tsne_df, method, k=clust_no, 
                    graph = FALSE)
                  clust.sil <- fviz_silhouette(clust.res, 
                    ggtheme = theme_minimal())
                  sil_title <- clust.sil$labels$title
                  sil_coeff <- as.numeric(unlist(regmatches(sil_title, 
                    gregexpr("[[:digit:]]+\\.*[[:digit:]]*",sil_title))))
                  df <- data.frame(method, clust_no, sil_coeff)
                  return(df)
            }
            df_new <- rbind(df_new,df_tsne_sil_stat)
        }
      df_new$method <- as.character(df_new$method)
      df_new$clust_no <-as.numeric(as.character(df_new$clust_no))
      df_new$sil_coeff <-as.double(as.character(df_new$sil_coeff))
    saveRDS(df_new,file = sprintf("./tsne_sil_stat_analysis_%s/tsne sil_coeff - %.0f clusters.rds", region_analysed, clust_no))
    write.csv(df_new, file = sprintf("./tsne_sil_stat_analysis_%s/sample_data_%.0f_cl.csv",region_analysed,clust_no))
    }
}
})

stopImplicitCluster() # clean up the cluster

#write.csv(df_new, file = "sample_data.csv")

# CHECK WHAT HAS BEEN WRITTEN TO THE RDS FILES
# CHECK=FALSE
# if(CHECK){
#  library(rlist)
#  list1 <- list()
#  for (clust_no in range_def){
#    list1 <- list.append(list1,readRDS(file = sprintf("./tsne_sil_stat_analysis_%s/tsne sil_coeff - %.0f clusters.rds", region_analysed, clust_no)))
#  }
#  print(list1)
#  print(dtime)
# }
# list1

# list1 <- list()
# for (clust_no in range_def){
#   list1 <- list.append(list1,readRDS(file = sprintf("./tsne_sil_stat_analysis_%s/tsne sil_coeff - %.0f clusters.rds", region_analysed, clust_no)))
# }
# print(list1)
