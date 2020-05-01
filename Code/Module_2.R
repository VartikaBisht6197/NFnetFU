#Create clusters using hierarchical clustering (dist = manhattan)
kNNdistplot(t(scaled_rules_int), k = 2)
res <- dbscan(t(scaled_rules_int),eps = 3.5, minPts = 2)

label_clust <- as.matrix(res[["cluster"]])
rownames(label_clust) <- colnames(scaled_rules_int)


#Group together the clusters as list of list
groups_we_need <- list()
for (i in 1:max(label_clust)) {
  group_we_need <- microbiome_ANFIS[["colnames.var"]][which(label_clust %in% i)]
  groups_we_need <- list.append(groups_we_need,group_we_need)
}

#Incorporate cluster to form new data frame
new_data1 = incorporate_groups(as.data.frame(rules_int),groups_we_need)
#new_data1 = incorporate_groups(as.data.frame(stimulated_data),groups_we_need)

#Missing Value heatplot : Original data
miss_val <- missing_val_plot(data1)
miss_original <- melt(data1)
jpeg("Missing Value heatplot Original Data.jpg", width = 1500, height = 1000)
ggplot(miss_original, aes(Var1, Var2, fill= miss_val)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="black") +
  theme(text = element_text(size=20), legend.title = element_blank()) +
  xlab("Samples") + ylab("Microbiome")
dev.off()

#Missing Value heatplot : Rule based data
miss_val <- missing_val_plot(rules_int)
miss_rule <- melt(rules_int)
jpeg("Missing Value heatplot Rule Based Data.jpg", width = 1500, height = 1000)
ggplot(miss_rule, aes(Var1, Var2, fill= miss_val)) + 
  geom_tile() +
  scale_fill_gradient(high="black") +
  theme(text = element_text(size=20), legend.title = element_blank()) +
  xlab("Samples") + ylab("Microbiome")
dev.off()

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

#Heat map using weights 
jpeg("Height map for distance between weights.jpg", width = 1000, height = 1000)
heatmap.2(as.matrix(dist(weights)),
          main = "Weights", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(15,12),     # widens margins around plot
)
dev.off()
