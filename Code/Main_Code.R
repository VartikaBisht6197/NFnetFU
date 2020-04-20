setwd("/Users/Vartika_Bisht/Individual_Project")
source("Libraries_Needed.R")
source("Borrowed_Functions.R")
source("Penalty_Function.R")
source("Incorporate_Groups.R")


#Input all data
#Data_Set_1 <- read.csv("Data_Set_1.csv")

stimulated_data <- data[1:20,2:71]

#Choose Microbiome Data
#data1 <- data.matrix(Data_Set_1[3:22])
data1 <- data.matrix(stimulated_data)

#Create labels for prediction 
#label_dat <- label_creator(Data_Set_1)
label_dat <- as.numeric(class_y[1:20])

#Learn the data to give rules
microbiome_ANFIS <- frbs.learn(data1, range.data = NULL, method.type = c("ANFIS"),control = list())

#Save Integer format of rules
rules_int <- microbiome_ANFIS[["rule.data.num"]]
colnames(rules_int) <- microbiome_ANFIS[["colnames.var"]]
write.csv(rules_int, "rules_int.csv")

#Scale and Save the Integer format of rules
scaled_rules_int <- scale(microbiome_ANFIS[["rule.data.num"]])
colnames(scaled_rules_int) <- microbiome_ANFIS[["colnames.var"]]
write.csv(scaled_rules_int, "scaled_rules_int.csv")

#Correlation plot for integer rules
jpeg("Correlation plot for integer rules.jpg", width = 1000, height = 1000)
corrplot(cor(rules_int), type = "upper")
dev.off()
#Correlation plot for original Data
jpeg("Correlation plot for original Data.jpg", width = 1000, height = 1000)
corrplot(cor(data1), type = "upper")
dev.off()

#Check smilarity/Correlation between new "rule-based" data and original data
cormat_TSK <- round(cor(rules_int),2)
cormat_data <- round(cor(data1),2)
print(paste0("Similarity between new rule-based data and original data : ", format( cor(c(cormat_TSK), c(cormat_data)), digits = 3, format = 'f' ) ) )

#Create clusters using hierarchical clustering (dist = manhattan)
kNNdistplot(t(scaled_rules_int), k = 2)
res <- dbscan(t(scaled_rules_int),eps = 2, minPts = 2)

label_clust <- as.matrix(res[["cluster"]])
rownames(label_clust) <- colnames(scaled_rules_int)


#Group together the clusters as list of list
groups_we_need <- list()
for (i in 1:max(label_clust)) {
  group_we_need <- microbiome_ANFIS[["colnames.var"]][which(label_clust %in% i)]
  groups_we_need <- list.append(groups_we_need,group_we_need)
}

#Incorporate cluster to form new data frame
#new_data1 = incorporate_groups(as.data.frame(Data_Set_1[3:22]),groups_we_need)
new_data1 = incorporate_groups(as.data.frame(stimulated_data),groups_we_need)

#Save the new data frame
write.csv(new_data1,"DataFrame_module7.csv")

#Correlation plot for Clusterd New Data
jpeg("Correlation plot for Clusterd New Data.jpg", width = 1000, height = 1000)
corrplot(cor(new_data1), type = "upper")
dev.off()

#calculate priors required for adaptive Lasso
prior_call <- cal_penalty(new_data1,label_dat)

#Load and save the new data frame
data_new <- prior_call[1][[1]]
write.csv(data_new, "data_new_for_LASSO.csv")

#Correlation plot for Clusterd New Data afte NA handle (if-valid)
if( length(colnames(new_data1)) != length(colnames(data_new)) ){
  print("Need for NA handle")
  jpeg("Correlation plot for Clusterd New Data with NA handled.jpg", width = 1000, height = 1000)
  corrplot(cor(data_new), type = "upper")
  dev.off()
}


#Load and save weights
weights <- in_2_dec(prior_call[2][[1]],2)
write.csv(weights, "weights.csv")

#load  penalty
coef = prior_call[3][[1]]
write.csv(coef, "coef.csv")

#Load and save initial betas
beta_initial = in_2_dec(prior_call[3][[1]][-1],5)
write.csv(beta_initial, "beta_initial.csv")

#Adaptive Lasso
alasso1_cv <- cv.glmnet(x = data_new , y = label_dat ,type.measure = "mse",nfold = 5,alpha = 0,penalty.factor = 1/abs(coef),keep = TRUE)
best_alasso_coef1 <- coef(alasso1_cv, s = alasso1_cv$lambda.min)
write.csv(best_alasso_coef1[-1], "Lasso_result.csv")

#Save weight network Information
weight_network = Network_start_end(weights,"euclidean")
write.csv(weight_network, "weight_network.csv", row.names=T)

#Height map using weights 
jpeg("Height map for distance between weights.jpg", width = 1000, height = 1000)
heatmap.2(as.matrix(dist(weights)),
          main = "Weights", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(15,12),     # widens margins around plot
          )
dev.off()


imp_variable <- data.frame(names(weights),weights,beta_initial,best_alasso_coef1[-1])
names(imp_variable) <- c("Microbiome",'weights','Initial_Betas', 'Adaptive_Lasso_Results')
write.csv(imp_variable, "All_results.csv")

p <- imp_variable %>%
  arrange(desc(Adaptive_Lasso_Results)) %>%
  mutate(text = paste("Weights: ", weights, "\nAdaptive Lasso Results: ", best_alasso_coef1[-1], "\nMicrobiome :",names(weights), "\nInitial_Betas :",beta_initial, sep="")) %>%
  
# Classic ggplot
ggplot(aes(x=weights, y=Initial_Betas, size = Adaptive_Lasso_Results, color = Microbiome, text=text), scale="globalminmax") +
  geom_point(alpha=0.7) +
  scale_size(range = c(1.4, 15), name="Adaptive Lasso Results") +
  scale_color_viridis(discrete=TRUE, guide=FALSE) +
  theme_minimal()+
  theme(legend.position="none")

pp <- ggplotly(p, tooltip="text")
pp
