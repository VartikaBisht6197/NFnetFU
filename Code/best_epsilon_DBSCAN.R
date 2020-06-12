
best_DBSCAN_E <- function(res,rules_int,label_dat){
  label_clust <- as.matrix(res[["cluster"]])
  rownames(label_clust) <- colnames(scaled_rules_int)
  colnames(label_clust) <- c("Cluster Number")
  groups_we_need <- list()
  for (i in 1:max(label_clust)) {
    group_we_need <- rownames(label_clust)[which(label_clust %in% i)]
    groups_we_need <- list.append(groups_we_need,group_we_need)
  }
  new_data1 <- incorporate_groups(as.data.frame(rules_int),groups_we_need,0)[[1]]
  #Scale the data
  n <- nrow(new_data1)
  label_dat_mean <- mean(label_dat)
  label_dat <- label_dat - label_dat_mean
  data_mean <- colMeans(new_data1)
  data_norm <- sqrt(n-1)*apply(new_data1,2,sd)
  new_data1 <- scale(new_data1, center = data_mean, scale = data_norm)
  
  
  #fit a logistic regression model to the data to calculate the weights
  lm.fit <- glm(label_dat ~ new_data1)
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
  if(length(data_NA) == 0){return(TRUE)}else{return(FALSE)}
}