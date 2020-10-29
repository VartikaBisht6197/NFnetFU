## Input : Scaled Rule Based Matrix (scaled_rules_int)

#DBSCAN
epsilons <- c(0.0001,0.001,0.01,0.1,seq(1,abs(max(scaled_rules_int)),0.5))
for (i in epsilons) {
  res <- dbscan(t(scaled_rules_int),eps = i, minPts = 2)
  if(max(res$cluster)!=min(res$cluster)){
    if(best_DBSCAN_E(res,rules_int,label_dat)){
      break
    }
  }
}
cluster_found <- max(unique(res$cluster))
print(sprintf("Epsilon value used : %s",i))
print(sprintf("%s cluster(s) found!",max(res$cluster)))
print("Clustering Done!")


label_clust <- as.matrix(res[["cluster"]])
rownames(label_clust) <- colnames(scaled_rules_int)
colnames(label_clust) <- c("Cluster Number")
write.csv(label_clust,"Feature's cluster number.csv")
print("Feature's cluster number saved")

groups_we_need <- list()
if(cluster_found > 0){
#Group together the clusters as list of list
print("Grouping Highly Colinear Features Together :-")
for (i in 1:max(label_clust)) {
  group_we_need <- rownames(label_clust)[which(label_clust %in% i)]
  #print(c(sprintf("Group %s :",i),group_we_need))
  groups_we_need <- list.append(groups_we_need,group_we_need)
}


#Incorporate cluster to form new data frame
print("Clubbing features in a group together")
#Input : data frame , list of clusters , print PCA Loadings used to combine the features
incorporate_groups_res <- incorporate_groups(as.data.frame(rules_int),groups_we_need,0)
#Output : data frame as matrix (new_data1) , PCA loadings for the group (PCA_loadings)
new_data1 <- incorporate_groups_res[1][[1]]
PCA_loadings <- as.data.frame(t(incorporate_groups_res[2][[1]]))
colnames(PCA_loadings) <- as.matrix(PCA_loadings[1,])
PCA_loadings <- PCA_loadings[-1,]
PCA_loadings <- as.matrix(PCA_loadings)
rownames(PCA_loadings) <- c("PCA Loadings")
print("Features Clubbed and incorporated in a new Data Frame!")


#Save the new data frame and PCA Loadings
write.csv(new_data1,"Rule Based matrix with Colinearity Handled.csv")
print("Rule Based matrix with Colinearity Handled saved")
write.csv(PCA_loadings,"PCA Loadings used to combine groups.csv" )
print("PCA Loadings used to combine groups saved")

#Correlation plot for Clusterd New Data
tiff("Correlation plot for New Clusterd Data.tiff", width = 10, height = 10, units = 'in', res = 300)
par(cex = 0.7)
corrplot(cor(new_data1), type = "upper")
dev.off()

#Significantly Correlated : No color = (p value cut off 0.05)
tiff("p value Correlation plot for clustered new.tiff", width = 10, height = 10, units = 'in', res = 300)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
par(cex = 0.4)
p.mat <- rcorr(new_data1, type ="spearman")
corrplot(p.mat$r, method = "color", col = col(200),number.cex = .7,
         type = "upper", addCoef.col = "black",tl.col = "black",
         p.mat = p.mat$P, sig.level = 0.05, insig = "blank", tl.srt = 90, diag = TRUE)
dev.off()
} else {print("No Clusters found")
  new_data1 <- as.data.frame(rules_int)}
## Output : Rule Based matrix with Colinearity Handled (new_data1) and PCA Loadings used to combine groups (PCA_loadings) 

