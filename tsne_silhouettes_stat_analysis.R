#rm(list=ls()) # to clean the environment

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


#----@@@@@@@@@@@----RUNNING SCRIPT ----@@@@@@@@@@@@@@@@@@@@

#set current folder as your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

source("functions.R")

region_analysed = "Bath and North East Somerset"

# create output folder
direction <- sprintf("./%s_CLUST_OUTPUTS", region_analysed)
dir.create(file.path(direction), showWarnings = FALSE) 

df_selection <- region_selector(region_analysed)


# create a direction for saving RDS files
direction <- sprintf("./%s_CLUST_OUTPUTS/RDS_tsne_sil_coeff", region_analysed)
dir.create(file.path(direction), showWarnings = FALSE) 

# create a direction for saving CSV files to store mean median and mean
direction <- sprintf("./%s_CLUST_OUTPUTS/MMM_tsne_sil_coeff", region_analysed)
dir.create(file.path(direction), showWarnings = FALSE) 



#defining range for number of clusters
range_def <- seq(6, 10, by=1)

#defining the number of iterations to be performed
no_iterations = 1000

# ---- FUNCTION TO EXTRACT DATA POINTS FOR SILHOUETTE COEFFICIENTS ------
# also saves mean median and mode into a csv file, for reference

for (i in range_def){
df_mc <- silhouettes1(df_selection,no_iterations,i)

# --- SAVE THE MEAN MEDIAN AND MODE FOR ALL 5 CLUSTERING METHODS ----
clustm <- c("kmeans", "pam", "agnes", "clara", "diana")

mmm_df <- data.frame()
for (method in clustm){
  mask_2 = df_mc$method == method
  df_mc_summary <- df_mc[which(mask_2), ]
  res.mean <- mean(df_mc_summary[,3])
  res.median <- median(df_mc_summary[,3])
  res.mode <- getmode(df_mc_summary[,3])
  x <- cbind(res.mean, res.median, res.mode)
  rownames(x) <- method
  mmm_df <- rbind(mmm_df,x)
}

write.csv(mmm_df, file = sprintf("%s_CLUST_OUTPUTS/MMM_tsne_sil_coeff/Mean - median - mode df %.0f clusters.csv", region_analysed, i))
saveRDS(df_mc,file = sprintf("./%s_CLUST_OUTPUTS/RDS_tsne_sil_coeff/%.0f_cl - tsne sil_coeff.rds", region_analysed,i))
}



# ----- READ THE RDS FILE AND MANAGE PLOTTING OF RESULTS ----- 
list1 <- list()
for (i in range_def){
  list1 <- list.append(list1,readRDS(file = sprintf("./%s_CLUST_OUTPUTS/RDS_tsne_sil_coeff/%.0f_cl - tsne sil_coeff.rds", region_analysed,i)))
}


#plotting for comparison of silhouette coefficients on multiple runs of tsne dim red algorithm
axis_text_size = 8
legend_text_size = 7
title_size = 8
leg_size = 15
leg_width =10


#define the plotting parameters as themes for boxplot and 
theme_1 = theme(legend.position="right", legend.title = element_text(size = axis_text_size),
                axis.title = element_text(size=axis_text_size),
                axis.text.x=element_blank(),
                axis.title.x = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_text(size = axis_text_size),
                legend.text=element_text(size=legend_text_size, margin = margin(t = 0)),
                plot.title = element_text(hjust = 0.5, size=title_size, face="bold"),
                axis.line=element_blank(),
                legend.key.size = unit(leg_size, "pt"),
                legend.key.width = unit(leg_width, "pt"),
                panel.grid.major.x = element_line(size = 0.1),
                panel.grid.major.y = element_line(size = 0.25),
                aspect.ratio = 3/4)

#define the plotting parameters as themes for density plots
theme_2 = theme(legend.position="right", legend.title = element_text(size = axis_text_size),
                axis.title = element_text(size=axis_text_size),
                axis.text.x=element_text(size = axis_text_size),
                axis.text = element_text(size = axis_text_size),
                legend.text=element_text(size=legend_text_size, margin = margin(t = 0)),
                plot.title = element_text(hjust = 0.5, size=title_size, face="bold"),
                axis.line=element_blank(),
                legend.key.size = unit(leg_size, "pt"),
                legend.key.width = unit(leg_width, "pt"),
                panel.grid.major.x = element_line(size = 0.1),
                panel.grid.major.y = element_line(size = 0.25))
                

myplots <- list()
k <- 1
for (i in range_def){
  df <- as.data.frame(list1[i-5])
  res.boxplot_jitter <- ggplot(df, aes(x=method, y=sil_coeff, fill=method)) +
    geom_boxplot()  + ggtitle(sprintf("t-SNE boxplot %.0f clusters",i)) + theme_1 #geom_jitter(shape=16, position=position_jitter(0.1))
    

  mu <- ddply(df, "method", summarise, grp.mean=mean(sil_coeff)) # calculate the mean for each method
  
  
  res.density <- ggplot(df, aes(x=sil_coeff, color = method, fill = method)) + 
    geom_density(alpha=.05)+
    scale_fill_manual(values=brewer.pal(n=5, name="Set1")) + geom_vline(data=mu, aes(xintercept=grp.mean, color=method), linetype="dashed") +
    scale_color_manual(values=brewer.pal(n=5, name="Set1")) + ggtitle(sprintf("t-SNE density plot %.0f clusters",i)) + theme_2
  res.density
  

  myplots[[k]] <- res.boxplot_jitter
  myplots[[k+1]] <- res.density
  k <- k + 2
  myplots
}

plot_width = 150
plot_height = 240

library(grid)
m_plots <- grid.arrange(grobs = myplots, ncol=2, nrow=(length(myplots)/2)) #, top = textGrob(paste(dim_red_method,region_analysed, "-", nclust, "clusters"), gp=gpar(fontsize=20,font=8)))
ggsave(filename = sprintf("%s_CLUST_OUTPUTS/t-sne clustering stat analysis summary.pdf",region_analysed), 
       plot = m_plots, width = plot_width, height = plot_height, units = "mm")

