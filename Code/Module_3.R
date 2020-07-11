## Input : Rule Based matrix with Colinearity Handled (new_data1) and PCA Loadings used to combine groups (PCA_loadings)

#calculate priors required for adaptive Lasso
prior_call <- cal_penalty(new_data1,label_dat)

#Load and save the new data frame
data_new <- prior_call[1][[1]]
write.csv(data_new, "Data Frame used for Adaptive Lasso.csv")

#Load and save weights
weights <- in_2_dec(prior_call[2][[1]],2)
write.csv(weights, "Intial Weights for Adaptive LASSO.csv")

#load  penalty
coef <- prior_call[3][[1]]
write.csv(coef, "Coefficient for Adaptive LASSO.csv")

#Load and save initial betas
beta_initial = in_2_dec(prior_call[3][[1]][-1],5)
write.csv(beta_initial, "Intial Betas for Adaptive LASSO.csv")

#Adaptive Lasso
alasso1_cv <- cv.glmnet(x = data_new , y = label_dat ,type.measure = "mse",nfold = 5,alpha = 0,penalty.factor = 1/abs(coef),keep = TRUE)
best_alasso_coef1 <- coef(alasso1_cv, s = alasso1_cv$lambda.min)
ADCoeff <- as.data.frame(as.matrix(scale(best_alasso_coef1[-1])))
rownames(ADCoeff) <- names(weights)
colnames(ADCoeff) <- c("Adaptive Lasso Coefficient")
write.csv(ADCoeff, "Adaptive Lasso Coefficient.csv")

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

#Heat map using Adaptive LASSO 
jpeg("Height map for distance between Adaptive LASSO.jpg", width = 1000, height = 1000)
heatmap.2(as.matrix(dist(ADCoeff)),
          main = "Adaptive LASSO ", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(15,12),     # widens margins around plot
)
dev.off()


#Calculating Overall scores for all features
index <- 1
names_features <- c()
ADres <- c()
for(i in regexpr('~', rownames(ADCoeff))){
  if(i > 0){
    features <- strsplit(rownames(ADCoeff)[index],"~")[[1]]
    val <- ADCoeff[rownames(ADCoeff)[index],]
    for(j in features){
      names_features <- c(names_features, j)
      ADres <- c(ADres,(as.numeric(PCA_loadings[,j])* val))
    }
  } else{
    names_features <- c(names_features,rownames(ADCoeff)[index])
    ADres <- c(ADres,ADCoeff[rownames(ADCoeff)[index],])
  }
  index <- index + 1
}

feature_scores <- as.data.frame(ADres)
rownames(feature_scores) <- names_features
feature_scores <- arrange(feature_scores, ADres)
write.csv(feature_scores,"Feature scores.csv")

#Histogram of Feature Scores
tiff("Feature Scores.tiff", width = 10, height = 10, units = 'in', res = 300)
par(cex = 0.7)
name <- rownames(feature_scores)
feature_scores %>%
  mutate(name=factor(name, levels=name)) %>%
  ggplot( aes(x=name, y=feature_scores$ADres)) +
  geom_segment( aes(xend=name, yend=0)) +
  geom_bar(stat="identity", fill="red", alpha=.6, width=.4)+
  coord_flip() +
  theme_bw() +
  xlab("") + ylab("Feature Score")
dev.off()

print("Feature scores computed and saved")

## Output : Feature scores (feature_scores)


