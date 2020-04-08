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
