source("Libraries_Needed.R")
source("Borrowed_Functions.R")
source("Penalty_Function.R")
source("Incorporate_Groups.R")

#Input all data
Data_Set_1 <- read.csv("Data_Set_1.csv")

#Choose Microbiome Data
data1 <- data.matrix(Data_Set_1[3:22])

#Create labels for prediction 
label_dat <- label_creator(Data_Set_1)

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
corr_plot(rules_int)
#Correlation plot for scaled integer rules
corr_plot(scaled_rules_int)
#Correlation plot for original Data
corr_plot(data1)

#Check smilarity/Correlation between new "rule-based" data and original data
cormat_TSK <- round(cor(rules_int),2)
cormat_data <- round(cor(data1),2)
print(paste0("Similarity between new rule-based data and original data : ", format( cor(c(cormat_TSK), c(cormat_data)), digits = 3, format = 'f' ) ) )

#Save start-end-edge netwrok format data frame
distance_TSK <- Network_start_end(t(scaled_rules_int),"manhattan")
write.csv(as.data.frame( distance_TSK ), "distance_TSK.csv")

#Create clusters using hierarchical clustering (dist = manhattan)
cluster_dist <- dist(t(scaled_rules_int), method = "manhattan", diag = FALSE, upper = FALSE, p = 2)
he_clust <- hclust(cluster_dist)
clusters <- color_branches( as.dendrogram(he_clust) , h = 3)
par(mfrow=c(1,1), mar=c(10,0.5,3,1))
plot(clusters)
label_clust <- as.matrix(cutree(clusters, k = NULL, h = 3))
clusplot(t(rules_int), label_clust, color=TRUE, shade=TRUE,
         labels=2, lines=0)


#Group together the clusters as list of list
groups_we_need <- list()
for (i in 1:max(label_clust)) {
  group_we_need <- microbiome_ANFIS[["colnames.var"]][which(label_clust %in% i)]
  groups_we_need <- list.append(groups_we_need,group_we_need)
}

#Incorporate cluster to form new data frame
new_data1 = incorporate_groups(as.data.frame(Data_Set_1[3:22]),groups_we_need)

#Save the new data frame
write.csv(new_data1,"DataFrame_module7.csv")

#calculate priors required for adaptive Lasso
prior_call <- cal_penalty(new_data1,label_dat)

#Load and save the new data frame
data_new <- prior_call[1][[1]]
write.csv(data_new, "data_new_for_LASSO.csv")

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
heatmap.2(as.matrix(dist(weights)),
          main = "Weights", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(15,12),     # widens margins around plot
          )



