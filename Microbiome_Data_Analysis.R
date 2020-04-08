library(readxl)
library(reshape2)
library(ggplot2)
library(rlist)
library(RCy3)
library(reshape2)
library(frbs)
library(philentropy)
library(dendextend)
library(cluster)
library(ggfortify)
library(glmnet)

#-------------------------------------------------------------------------------------------

#Function to chage a list of float numbers to n decimal places float number
#Input : list , int (decimal places)
#Output : list
in_2_dec <- function(list_given,n){
  j = 1
  for (i in list_given){
    list_given[j] = formatC(i, digits = n, format = "f")
    j = j + 1
  }
  return(list_given)
}

#Function ti plot correlation plot
#Input : data frame
#Output : plot
corr_plot <- function(col_var_data){
  cormat_TSK <- round(cor(col_var_data),2)
  melted_cormat_TSK <- melt(cormat_TSK)
  ggplot(data = melted_cormat_TSK, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()}

#Function to create start-end-edge node network data frame
#Input : data frame
#Output : data frame
Network_start_end <- function(df_need){
  inDist <- dist(t(scaled_rules_int), method = "manhattan", diag = FALSE, upper = FALSE, p = 2)
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  d = data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
  return(d)
}

#Function to Incorporate clusters in the new data frame
#Input : data frame , list of clusters
#Output : data frame
incorporate_groups <- function(df,group_list){
  for (i in group_list) {
    if(length(i) > 1){
      df_rest = df[, !(names(df) %in% i)]
      df_group = df[, (names(df) %in% i)]
      res.pca <- prcomp(df_group)
      var_explain = summary(res.pca)[["importance"]][2][1]
      if(var_explain <0.80){print("Variance Explained by PC1 for the following group is < 80%")
        print(i) 
        print(var_explain)}
      pca.loadings <- res.pca$rotation
      col_name = c()
      new_col = 0
      k = 1
      for (j in colnames(df_group)) {
        new_col = new_col + (df[,j]*pca.loadings[k][1])
        col_name = c(col_name,substr(j,1,5))
        k = k + 1
      }
      col_name = paste(col_name, sep = "", collapse = "~")
      df_rest$col_name = new_col
      names(df_rest)[names(df_rest) == "col_name"] <- col_name
      df = df_rest}
    
  }
  return(data.matrix(df))
}

#Create labels for data using factors
#Input : list - num[]
#Output : list - num[]
label_creator <- function(Data_Set_1){
  label_dat <- as.vector(Data_Set_1[2])
  factors <- factor(label_dat$Class)
  label_dat <- as.numeric(factors)
  return(label_dat)}

#Initial penalty calculator
#Input : data frame(features) ,data frame(estimate as factor)
#Output : data frame, weights, penalty , initial betas
cal_penalty <- function(data,label_dat){
  
  #Scale the data
  n <- nrow(data)
  label_dat_mean <- mean(label_dat)
  label_dat <- label_dat - label_dat_mean
  data_mean <- colMeans(data)
  data_norm <- sqrt(n-1)*apply(data,2,sd)
  data <- scale(data, center = data_mean, scale = data_norm)
  
  
  #fit a logistic regression model to the data to calculate the weights
  lm.fit <- glm(label_dat ~ data)
  beta.init <- coef(lm.fit)[-1]
  beta.init_new <- beta.init[!is.na(beta.init)]
  beta.init_NA <- beta.init[-match(beta.init_new,beta.init)]
  w  <- abs(beta.init_new) 
  
  #list NA 
  data_NA <- list()
  j = 1
  for (i in names(beta.init_NA)) {
    data_NA[j] <- substr(i,5,nchar(i))
    j = j + 1
  }
  
  #Indicate NA values found
  if(length(data_NA) == 0){print("No NA coefficient encounter while fitting logistic regression model to get weights")}
  else{print("NA values found for the following elements")
    print(data_NA)}

  #Remove features with NA values
  data_new <- data[,setdiff(colnames(data), data_NA)]
  
  #Scale the data by 1/weight[]
  data_scale <- scale(data_new, center=FALSE, scale=1/w)  
  
  #Calculate the initial coefficient
  lasso.fit <- cv.glmnet(data_scale, label_dat, family = "gaussian", alpha = 0, standardize = FALSE, nfolds = 5)
  beta <- predict(lasso.fit, data_scale, type="coefficients", s="lambda.min")[-1]
  beta_initial <- beta
  
  #calculate penalties
  beta <- beta * w
  beta <- matrix(beta, nrow=1)
  data_mean <- matrix(data_mean)
  b0 <- label_dat_mean -  beta%*%data_mean 
  coef <- cbind(b0, beta)
  
  #Name the colums in weight matrix
  names(w) <- colnames(data_new)
  
  return_val <- list(data_new,w,coef,beta_initial)
  
  return(return_val)
}

#-------------------------------------------------------------------------------------------

#Input all data
Data_Set_1 <- read_excel("Data_Set_1.xlsx")

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
print("Similarity between new rule-based data and original data")
print( cor(c(cormat_TSK), c(cormat_data)) )

#Save start-end-edge netwrok format data frame
distance_TSK <- Network_start_end(rules_int)
write.csv(as.data.frame( distance_TSK ), "distance_TSK.csv")

#Create clusters using hierarchical clustering (dist = manhattan)
cluster_dist <- dist(t(scaled_rules_int), method = "manhattan", diag = FALSE, upper = FALSE, p = 2)
he_clust <- hclust(cluster_dist)
clusters <- color_branches( as.dendrogram(he_clust) , h = 3)
par(mfrow=c(1,1), mar=c(10,0.5,3,1))
plot(clusters)
plot(he_clust)
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

#par(mfrow=c(1,1), mar=c(3,1,3,2))
#plot(beta_initial,weights,col="lightblue", pch=19, cex=2)
#text(beta_initial,weights, label = 1:length(weights), cex=0.9, font=2 )

#data_plot <- data.frame(beta_initial,weights)
#ggplot(data_plot, aes(x=beta_initial, y=weights, group = 1)) + geom_point(size=10, shape=23)+  geom_text(label=1:length(weights))

# why abs w?
# penalty.factor = 1/abs(coef) ? abs?
