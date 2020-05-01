setwd("/Users/Vartika_Bisht/Individual_Project")
source("Libraries_Needed.R")
source("Borrowed_Functions.R")
source("Penalty_Function.R")
source("Incorporate_Groups.R")
source("Missing_val.R")

#=========================================================================================================
#Input all data
' Dataset_1 , Dataset_2 , Dataset_3 , Stimulated_Dataset '
# If Stimulated_Dataset :
set.seed(10)
data_analysis = "Stimulated_Dataset"

source("Input_data.R")

#=========================================================================================================
# Module 1
source("Module_1.R")

#=========================================================================================================
# Module 2
source("Module_2.R")

#=========================================================================================================
# Module 3

#cluster sizes
cluster_size = c()

for (i in names(weights)) {
  cs = sum(charToRaw(i) == charToRaw("~")) + 1
  cluster_size = c(cluster_size,cs)
}

imp_variable <- data.frame(names(weights),weights,as.numeric(beta_initial),as.numeric(best_alasso_coef1[-1]),as.factor(cluster_size))
names(imp_variable) <- c("Microbiome",'weights','Initial_Betas', 'Adaptive_Lasso_Results', 'Cluster_Size')
write.csv(imp_variable, "All_results.csv")

'# Intractive Visualisation 
p <- imp_variable %>%
  mutate(text = paste("Weights: ", weights, "\nCluster Size: ", cluster_size , "\nMicrobiome :", rownames(imp_variable) , "\nInitial Betas :",beta_initial,"\nAdaptive Lasso Results: ", best_alasso_coef1[-1], sep="")) %>%
  
# Classic ggplot
ggplot(aes(x=weights, y=Initial_Betas, size = Adaptive_Lasso_Results, text=text), scale="std") +
  geom_point(alpha=0.7, aes( colour = Cluster_Size) , show.legend = TRUE) +
  scale_size(range = c(1.4, 15), name="Cluster Size") +
  scale_color_viridis(discrete=TRUE, guide=FALSE) +
  theme_minimal()+
  theme(legend.position="bottom")

pp <- ggplotly(p, tooltip="text")
pp'


#Only for simulated data
if(data_analysis == "Stimulated_Dataset"){
  seed = seq(15,100,5)
  df_list = c(imp_variable)
  for (s in seed) {
    set.seed(s)
    source("Input_data.R")
    source("Module_1.R")
    source("Module_2.R")
    
    #cluster sizes
    cluster_size = c()
    for (i in names(weights)) {
      cs = sum(charToRaw(i) == charToRaw("~")) + 1
      cluster_size = c(cluster_size,cs)
    }
    
    imp_variable <- data.frame(names(weights),weights,as.numeric(beta_initial),as.numeric(best_alasso_coef1[-1]),as.factor(cluster_size))
    names(imp_variable) <- c("Microbiome",'weights','Initial_Betas', 'Adaptive_Lasso_Results', 'Cluster_Size')
    
    df_list = c(df_list,imp_variable)
  }
}
